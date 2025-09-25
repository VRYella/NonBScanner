#!/usr/bin/env python3
"""
Motif detection utilities and common functions
- Shared scoring and normalization functions
- Sequence processing utilities  
- Output formatting and standardization
- Common statistical and analytical functions

This module provides common utilities used by both A-philic DNA and Z-DNA
scanners, ensuring consistent behavior and reducing code duplication.

Author: Dr. Venkata Rajesh Yella
Integration: 2024 - NBDFinder Common Utilities
"""

import re
import math
import numpy as np
from typing import Dict, List, Any, Tuple, Optional
import sys
import os


def parse_fasta(fasta_str: str) -> str:
    """
    Parse FASTA string to clean DNA sequence.
    
    | Parameter  | Type | Description                      | Range      |
    |------------|------|----------------------------------|------------|
    | fasta_str  | str  | FASTA formatted string          | any        |
    | return     | str  | Clean DNA sequence              | any        |
    """
    return "".join([line.strip() for line in fasta_str.split('\n') 
                   if not line.startswith(">")]).upper().replace(" ", "").replace("U", "T")


def wrap_sequence(seq: str, width: int = 60) -> str:
    """
    Format sequence with line breaks for display.
    
    | Parameter | Type | Description                    | Default | Range      |
    |-----------|------|--------------------------------|---------|------------|
    | seq       | str  | DNA sequence to format         | -       | any        |
    | width     | int  | Characters per line            | 60      | 10-200     |
    | return    | str  | Formatted sequence with breaks | -       | any        |
    """
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


def gc_content(seq: str) -> float:
    """
    Calculate GC content percentage.
    
    | Parameter | Type | Description                    | Range        |
    |-----------|------|--------------------------------|--------------|
    | seq       | str  | DNA sequence to analyze        | any          |
    | return    | float| GC content percentage          | 0.0 - 100.0  |
    """
    if not seq:
        return 0.0
    gc = seq.count('G') + seq.count('C')
    return (gc / len(seq)) * 100


def reverse_complement(seq: str) -> str:
    """
    Generate reverse complement of DNA sequence.
    
    | Parameter | Type | Description                    | Range      |
    |-----------|------|--------------------------------|------------|
    | seq       | str  | DNA sequence to complement     | any        |
    | return    | str  | Reverse complement sequence    | any        |
    """
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement_map.get(base, 'N') for base in seq[::-1])


def is_palindrome(seq: str) -> bool:
    """
    Check if sequence is palindromic (reads same on both strands).
    
    | Parameter | Type | Description                    | Range      |
    |-----------|------|--------------------------------|------------|
    | seq       | str  | DNA sequence to check          | any        |
    | return    | bool | True if palindromic            | True/False |
    """
    return seq == reverse_complement(seq)


def normalize_score_linear(raw_score: float, min_val: float, max_val: float) -> float:
    """
    Normalize score using linear scaling to 0-100 range.
    
    | Parameter | Type | Description                    | Range        |
    |-----------|------|--------------------------------|--------------|
    | raw_score | float| Raw score to normalize         | any          |
    | min_val   | float| Minimum expected score         | any          |
    | max_val   | float| Maximum expected score         | any          |
    | return    | float| Normalized score (0-100)       | 0.0 - 100.0  |
    """
    if max_val <= min_val:
        return 50.0  # Default middle value
    
    normalized = ((raw_score - min_val) / (max_val - min_val)) * 100.0
    return max(0.0, min(100.0, normalized))


def normalize_score_sigmoid(raw_score: float, midpoint: float, steepness: float = 1.0) -> float:
    """
    Normalize score using sigmoid function for smooth transitions.
    
    | Parameter | Type | Description                    | Default | Range        |
    |-----------|------|--------------------------------|---------|--------------|
    | raw_score | float| Raw score to normalize         | -       | any          |
    | midpoint  | float| Score value at sigmoid center  | -       | any          |
    | steepness | float| Sigmoid curve steepness        | 1.0     | 0.1 - 10.0   |
    | return    | float| Normalized score (0-100)       | -       | 0.0 - 100.0  |
    """
    try:
        exp_val = math.exp(-steepness * (raw_score - midpoint))
        sigmoid = 1.0 / (1.0 + exp_val)
        return sigmoid * 100.0
    except (OverflowError, ZeroDivisionError):
        return 50.0  # Default middle value


def calculate_sequence_complexity(seq: str) -> float:
    """
    Calculate sequence complexity using Shannon entropy.
    
    | Parameter | Type | Description                    | Range        |
    |-----------|------|--------------------------------|--------------|
    | seq       | str  | DNA sequence to analyze        | any          |
    | return    | float| Complexity score (entropy)     | 0.0 - 2.0    |
    """
    if not seq:
        return 0.0
    
    # Count base frequencies
    base_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    for base in seq:
        if base in base_counts:
            base_counts[base] += 1
    
    # Calculate Shannon entropy
    length = len(seq)
    entropy = 0.0
    
    for count in base_counts.values():
        if count > 0:
            frequency = count / length
            entropy -= frequency * math.log2(frequency)
    
    return entropy


def find_tandem_repeats(seq: str, min_unit_len: int = 2, max_unit_len: int = 10) -> List[Dict[str, Any]]:
    """
    Find tandem repeat regions in sequence.
    
    | Parameter     | Type | Description                    | Default | Range      |
    |---------------|------|--------------------------------|---------|------------|
    | seq           | str  | DNA sequence to search         | -       | any        |
    | min_unit_len  | int  | Minimum repeat unit length     | 2       | 1-20       |
    | max_unit_len  | int  | Maximum repeat unit length     | 10      | min_unit+  |
    | return        | list | List of tandem repeat regions  | []      | any        |
    """
    repeats = []
    seq = seq.upper()
    
    for unit_len in range(min_unit_len, min(max_unit_len + 1, len(seq) // 2 + 1)):
        for start in range(len(seq) - unit_len + 1):
            unit = seq[start:start + unit_len]
            
            # Count consecutive repeats
            count = 1
            pos = start + unit_len
            
            while pos + unit_len <= len(seq) and seq[pos:pos + unit_len] == unit:
                count += 1
                pos += unit_len
            
            # Only record if we have at least 2 repeats
            if count >= 2:
                repeat_len = count * unit_len
                repeats.append({
                    'Start': start + 1,  # 1-based
                    'End': start + repeat_len,
                    'Length': repeat_len,
                    'Unit': unit,
                    'Unit_Length': unit_len,
                    'Copy_Number': count,
                    'Sequence': seq[start:start + repeat_len]
                })
    
    return repeats


def merge_overlapping_motifs(motifs: List[Dict[str, Any]], max_distance: int = 10) -> List[Dict[str, Any]]:
    """
    Merge nearby or overlapping motifs to reduce redundancy.
    
    | Parameter    | Type | Description                    | Default | Range      |
    |--------------|------|--------------------------------|---------|------------|
    | motifs       | list | List of motif dictionaries    | -       | any        |
    | max_distance | int  | Max distance to merge motifs   | 10      | 0-50       |
    | return       | list | Merged motif list              | []      | any        |
    """
    if not motifs:
        return []
    
    # Sort motifs by start position
    sorted_motifs = sorted(motifs, key=lambda x: x.get('Start', 0))
    merged = []
    current = sorted_motifs[0].copy()
    
    for motif in sorted_motifs[1:]:
        current_end = current.get('End', current.get('Start', 0))
        motif_start = motif.get('Start', 0)
        
        # Check if motifs should be merged
        if motif_start - current_end <= max_distance:
            # Merge motifs - keep the best score
            current['End'] = max(current_end, motif.get('End', motif_start))
            current['Length'] = current['End'] - current.get('Start', 0) + 1
            
            # Keep the higher score
            current_score = current.get('Raw_Score', current.get('Score', 0))
            motif_score = motif.get('Raw_Score', motif.get('Score', 0))
            if motif_score > current_score:
                current['Raw_Score'] = motif_score
                current['Score'] = motif_score
                current['Classification'] = motif.get('Classification', current.get('Classification'))
        else:
            # Start new region
            merged.append(current)
            current = motif.copy()
    
    # Add the last motif
    merged.append(current)
    
    return merged


def filter_motifs_by_quality(motifs: List[Dict[str, Any]], 
                            min_score: float = 0.0, 
                            min_length: int = 10,
                            max_overlap: float = 0.5) -> List[Dict[str, Any]]:
    """
    Filter motifs based on quality criteria.
    
    | Parameter   | Type | Description                    | Default | Range      |
    |-------------|------|--------------------------------|---------|------------|
    | motifs      | list | List of motif dictionaries    | -       | any        |
    | min_score   | float| Minimum score threshold        | 0.0     | any        |
    | min_length  | int  | Minimum motif length           | 10      | 1-100      |
    | max_overlap | float| Max allowed overlap fraction   | 0.5     | 0.0 - 1.0  |
    | return      | list | Filtered motif list            | []      | any        |
    """
    filtered = []
    
    for motif in motifs:
        # Check score threshold
        score = motif.get('Raw_Score', motif.get('Score', 0))
        if score < min_score:
            continue
            
        # Check length threshold
        length = motif.get('Length', 0)
        if length < min_length:
            continue
            
        # Check for excessive overlap with existing motifs
        overlaps_too_much = False
        for existing in filtered:
            overlap = calculate_overlap_fraction(motif, existing)
            if overlap > max_overlap:
                overlaps_too_much = True
                break
        
        if not overlaps_too_much:
            filtered.append(motif)
    
    return filtered


def calculate_overlap_fraction(motif1: Dict[str, Any], motif2: Dict[str, Any]) -> float:
    """
    Calculate fraction of overlap between two motifs.
    
    | Parameter | Type | Description                    | Range      |
    |-----------|------|--------------------------------|------------|
    | motif1    | dict | First motif dictionary         | any        |
    | motif2    | dict | Second motif dictionary        | any        |
    | return    | float| Overlap fraction (0-1)         | 0.0 - 1.0  |
    """
    start1, end1 = motif1.get('Start', 0), motif1.get('End', 0)
    start2, end2 = motif2.get('Start', 0), motif2.get('End', 0)
    
    # Calculate overlap
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    
    if overlap_start >= overlap_end:
        return 0.0
    
    overlap_length = overlap_end - overlap_start + 1
    min_length = min(end1 - start1 + 1, end2 - start2 + 1)
    
    return overlap_length / min_length if min_length > 0 else 0.0


def standardize_motif_output(motif: Dict[str, Any], sequence_name: str = "", motif_id: int = 0) -> Dict[str, Any]:
    """
    Standardize motif output format for consistent reporting.
    
    | Parameter     | Type | Description                    | Default | Range      |
    |---------------|------|--------------------------------|---------|------------|
    | motif         | dict | Motif dictionary to standardize| -       | any        |
    | sequence_name | str  | Name of source sequence        | ""      | any        |
    | motif_id      | int  | Unique motif identifier        | 0       | any        |
    | return        | dict | Standardized motif dictionary  | {}      | any        |
    """
    standardized = motif.copy()
    
    # Ensure required fields are present
    if 'Sequence_Name' not in standardized:
        standardized['Sequence_Name'] = sequence_name
    
    if 'Motif_ID' not in standardized:
        motif_class = standardized.get('Class', 'Unknown')
        standardized['Motif_ID'] = f"{motif_class}_{motif_id}"
    
    # Standardize coordinate system (1-based)
    if 'Start' in standardized and standardized['Start'] < 1:
        standardized['Start'] = standardized['Start'] + 1
    
    # Ensure score fields are present
    if 'Raw_Score' not in standardized and 'Score' in standardized:
        standardized['Raw_Score'] = standardized['Score']
    
    if 'Score' not in standardized and 'Raw_Score' in standardized:
        standardized['Score'] = standardized['Raw_Score']
    
    # Add normalized score if not present
    if 'Normalized_Score' not in standardized:
        raw_score = standardized.get('Raw_Score', 0)
        standardized['Normalized_Score'] = normalize_score_linear(raw_score, 0, 100)
    
    # Ensure length is calculated correctly
    if 'Length' not in standardized or standardized['Length'] <= 0:
        start = standardized.get('Start', 1)
        end = standardized.get('End', start)
        standardized['Length'] = max(1, end - start + 1)
    
    return standardized


def validate_motif(motif: Dict[str, Any], seq_length: int) -> bool:
    """
    Validate that a motif has required fields and valid coordinates.
    
    | Parameter  | Type | Description                    | Range      |
    |------------|------|--------------------------------|------------|
    | motif      | dict | Motif dictionary to validate   | any        |
    | seq_length | int  | Length of source sequence      | >0         |
    | return     | bool | True if motif is valid         | True/False |
    """
    # Check required fields
    required_fields = ['Start', 'End', 'Length', 'Class']
    for field in required_fields:
        if field not in motif:
            return False
    
    # Check coordinate validity
    start = motif.get('Start', 0)
    end = motif.get('End', 0)
    length = motif.get('Length', 0)
    
    if start < 1 or end < start or end > seq_length:
        return False
    
    if length != (end - start + 1):
        return False
    
    return True


# === Statistics and analysis functions ===

def calculate_motif_density(motifs: List[Dict[str, Any]], sequence_length: int) -> float:
    """
    Calculate motif density (motifs per kilobase).
    
    | Parameter       | Type | Description                    | Range        |
    |-----------------|------|--------------------------------|--------------|
    | motifs          | list | List of motif dictionaries    | any          |
    | sequence_length | int  | Length of analyzed sequence    | >0           |
    | return          | float| Motifs per kilobase            | 0.0+         |
    """
    if sequence_length <= 0:
        return 0.0
    
    return (len(motifs) * 1000.0) / sequence_length


def get_motif_distribution(motifs: List[Dict[str, Any]]) -> Dict[str, int]:
    """
    Get distribution of motifs by class/subclass.
    
    | Parameter | Type | Description                    | Range      |
    |-----------|------|--------------------------------|------------|
    | motifs    | list | List of motif dictionaries    | any        |
    | return    | dict | Distribution by motif type     | any        |
    """
    distribution = {}
    
    for motif in motifs:
        motif_class = motif.get('Class', 'Unknown')
        subclass = motif.get('Subclass', '')
        
        key = f"{motif_class}" + (f" ({subclass})" if subclass else "")
        distribution[key] = distribution.get(key, 0) + 1
    
    return distribution


# === Example usage and testing ===
if __name__ == "__main__":
    # Test sequence processing functions
    test_seq = "ATCGATCGATCG"
    print(f"Original: {test_seq}")
    print(f"Reverse complement: {reverse_complement(test_seq)}")
    print(f"GC content: {gc_content(test_seq):.1f}%")
    print(f"Complexity: {calculate_sequence_complexity(test_seq):.2f}")
    print(f"Is palindrome: {is_palindrome('ATCGAT')}")
    
    # Test normalization
    raw_scores = [10, 25, 50, 75, 100]
    print("\nScore normalization:")
    for score in raw_scores:
        linear = normalize_score_linear(score, 0, 100)
        sigmoid = normalize_score_sigmoid(score, 50, 0.1)
        print(f"  Raw: {score} -> Linear: {linear:.1f}, Sigmoid: {sigmoid:.1f}")