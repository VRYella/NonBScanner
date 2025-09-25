#!/usr/bin/env python3
"""
Z-DNA motif scanner
- Uses 10-mer pattern scoring table
- Detects left-handed Z-form DNA structures
- Scores windows based on Z-forming potential
- Classifies regions as High/Moderate/Low Z-DNA forming

This module implements Z-DNA detection using the same approach as A-philic DNA.
Z-DNA represents left-handed double helical DNA structures that form under
superhelical tension and specific sequence contexts.

Scientific References:
- Z-DNA discovery: Wang et al. Nature 1979
- Z-DNA structural properties: Rich & Zhang Nature Reviews 2003  
- CG dinucleotide Z-forming potential: Ho et al. PNAS 1986
- Z-DNA biological functions: Herbert & Rich Genetica 1999

Author: Dr. Venkata Rajesh Yella
Integration: 2024 - NBDFinder Z-DNA Implementation
"""

import math
import csv
import sys
import os
from collections import Counter
from typing import List, Dict, Any, Tuple

# === Z-DNA 10-mer scoring table ===
# | Pattern        | Type   | Z-forming Score                    | Range        | Description                          |
# |----------------|--------|------------------------------------|--------------|--------------------------------------|
# | AACGCGCGCG     | float  | Z-forming potential score          | 50.0 - 63.0  | Moderate Z-forming sequence         |
# | GCGCGCGCGC     | float  | Maximum Z-forming score            | 63.0         | Strongest Z-forming potential       |
# | Mixed CG       | float  | Variable based on CG content       | 50.0 - 58.0  | Moderate to strong Z-forming        |

ZDNA_10MER_SCORES = {
    "AACGCGCGCG": 50.25,
    "ATGCGCGCGC": 51.25,
    "ATCGCGCGCG": 50,
    "AGCGCGCGCA": 50.25,
    "AGCGCGCGCG": 56,
    "ACGGCGCGCG": 50.25,
    "ACGCGGCGCG": 50.25,
    "ACGCGCGGCG": 50.25,
    "ACGCGCGCGA": 50.25,
    "ACGCGCGCGT": 51.5,
    "ACGCGCGCGG": 50.25,
    "ACGCGCGCGC": 57.25,
    "ACGCGCGCCG": 50.25,
    "ACGCGCCGCG": 50.25,
    "ACGCCGCGCG": 50.25,
    "ACCGCGCGCG": 50.25,
    "TAGCGCGCGC": 50,
    "TACGCGCGCG": 51.25,
    "TTGCGCGCGC": 50.25,
    "TGGCGCGCGC": 50.25,
    "TGCGGCGCGC": 50.25,
    "TGCGCGGCGC": 50.25,
    "TGCGCGCGGC": 50.25,
    "TGCGCGCGCA": 51.5,
    "TGCGCGCGCT": 50.25,
    "TGCGCGCGCG": 57.25,
    "TGCGCGCGCC": 50.25,
    "TGCGCGCCGC": 50.25,
    "TGCGCCGCGC": 50.25,
    "TGCCGCGCGC": 50.25,
    "TCGCGCGCGT": 50.25,
    "TCGCGCGCGC": 56,
    "GACGCGCGCG": 50.25,
    "GTGCGCGCGC": 51.5,
    "GTCGCGCGCG": 50.25,
    "GGCGCGCGCA": 50.25,
    "GGCGCGCGCG": 56,
    "GCAGCGCGCG": 50.25,
    "GCACGCGCGC": 51.5,
    "GCTGCGCGCG": 50.25,
    "GCGACGCGCG": 50.25,
    "GCGTGCGCGC": 51.5,
    "GCGTCGCGCG": 50.25,
    "GCGGCGCGCA": 50.25,
    "GCGGCGCGCG": 56,
    "GCGCAGCGCG": 50.25,
    "GCGCACGCGC": 51.5,
    "GCGCTGCGCG": 50.25,
    "GCGCGACGCG": 50.25,
    "GCGCGTGCGC": 51.5,
    "GCGCGTCGCG": 50.25,
    "GCGCGGCGCA": 50.25,
    "GCGCGGCGCG": 56,
    "GCGCGCAGCG": 50.25,
    "GCGCGCACGC": 51.5,
    "GCGCGCTGCG": 50.25,
    "GCGCGCGACG": 50.25,
    "GCGCGCGTGC": 51.5,
    "GCGCGCGTCG": 50.25,
    "GCGCGCGGCA": 50.25,
    "GCGCGCGGCG": 56,
    "GCGCGCGCAA": 50.25,
    "GCGCGCGCAT": 51.25,
    "GCGCGCGCAG": 50.25,
    "GCGCGCGCAC": 51.5,
    "GCGCGCGCTA": 50,
    "GCGCGCGCTG": 50.25,
    "GCGCGCGCGA": 56,
    "GCGCGCGCGT": 57.25,
    "GCGCGCGCGG": 56,
    "GCGCGCGCGC": 63,
    "GCGCGCGCCA": 50.25,
    "GCGCGCGCCG": 56,
    "GCGCGCCGCA": 50.25,
    "GCGCGCCGCG": 56,
    "GCGCCGCGCA": 50.25,
    "GCGCCGCGCG": 56,
    "GCCGCGCGCA": 50.25,
    "GCCGCGCGCG": 56,
    "CAGCGCGCGC": 50.25,
    "CACGCGCGCG": 51.5,
    "CTGCGCGCGC": 50.25,
    "CGACGCGCGC": 50.25,
    "CGTGCGCGCG": 51.5,
    "CGTCGCGCGC": 50.25,
    "CGGCGCGCGT": 50.25,
    "CGGCGCGCGC": 56,
    "CGCAGCGCGC": 50.25,
    "CGCACGCGCG": 51.5,
    "CGCTGCGCGC": 50.25,
    "CGCGACGCGC": 50.25,
    "CGCGTGCGCG": 51.5,
    "CGCGTCGCGC": 50.25,
    "CGCGGCGCGT": 50.25,
    "CGCGGCGCGC": 56,
    "CGCGCAGCGC": 50.25,
    "CGCGCACGCG": 51.5,
    "CGCGCTGCGC": 50.25,
    "CGCGCGACGC": 50.25,
    "CGCGCGTGCG": 51.5,
    "CGCGCGTCGC": 50.25,
    "CGCGCGGCGT": 50.25,
    "CGCGCGGCGC": 56,
    "CGCGCGCAGC": 50.25,
    "CGCGCGCACG": 51.5,
    "CGCGCGCTGC": 50.25,
    "CGCGCGCGAT": 50,
    "CGCGCGCGAC": 50.25,
    "CGCGCGCGTA": 51.25,
    "CGCGCGCGTT": 50.25,
    "CGCGCGCGTG": 51.5,
    "CGCGCGCGTC": 50.25,
    "CGCGCGCGGT": 50.25,
    "CGCGCGCGGC": 56,
    "CGCGCGCGCA": 57.25,
    "CGCGCGCGCT": 56,
    "CGCGCGCGCG": 63,
    "CGCGCGCGCC": 56,
    "CGCGCGCCGT": 50.25,
    "CGCGCGCCGC": 56,
    "CGCGCCGCGT": 50.25,
    "CGCGCCGCGC": 56,
    "CGCCGCGCGT": 50.25,
    "CGCCGCGCGC": 56,
    "CCGCGCGCGT": 50.25,
    "CCGCGCGCGC": 56
}

# === Thresholds ===
# | Parameter      | Type  | Value | Description                                |
# |----------------|-------|-------|-------------------------------------------|
# | HIGH_ZDNA      | float | 58.0  | Threshold for high Z-DNA forming regions |
# | MOD_ZDNA       | float | 53.0  | Threshold for moderate Z-DNA regions     |
# | LOW_ZDNA       | float | 50.5  | Minimum threshold for Z-DNA detection    |
# | MIN_LENGTH     | int   | 10    | Minimum Z-DNA region length             |
# | MAX_LENGTH     | int   | 30    | Maximum Z-DNA region length to scan     |

HIGH_ZDNA = 58.0
MOD_ZDNA = 53.0
LOW_ZDNA = 50.5
MIN_LENGTH = 10
MAX_LENGTH = 30


def calculate_zdna_score(sequence: str) -> float:
    """
    Calculate Z-DNA score for a sequence using 10-mer scoring table.
    
    | Parameter | Type | Description                    | Range        |
    |-----------|------|--------------------------------|--------------|
    | sequence  | str  | DNA sequence to score          | >=10 bp      |
    | return    | float| Mean Z-forming potential score | 0.0 - 63.0   |
    """
    if len(sequence) < 10:
        return 0.0
    
    total_score = 0.0
    count = 0
    
    for i in range(len(sequence) - 9):
        nmer = sequence[i:i+10]
        if nmer in ZDNA_10MER_SCORES:
            total_score += ZDNA_10MER_SCORES[nmer]
            count += 1
    
    return 0.0 if count == 0 else total_score / count


def classify_zdna_region(score: float, length: int, cg_content: float) -> str:
    """
    Classify Z-DNA region based on score and sequence characteristics.
    
    | Parameter  | Type | Description                      | Range      |
    |------------|------|----------------------------------|------------|
    | score      | float| Mean Z-forming potential score   | 0.0 - 63.0 |
    | length     | int  | Length of sequence region        | >=10       |
    | cg_content | float| CG dinucleotide content %        | 0.0 - 100.0|
    | return     | str  | Classification label             | str        |
    """
    if score >= HIGH_ZDNA and cg_content >= 60.0:
        return "High-confidence Z-DNA"
    elif score >= MOD_ZDNA and cg_content >= 40.0:
        return "Moderate Z-DNA"
    elif score >= LOW_ZDNA:
        return "Low Z-DNA"
    else:
        return "Non-Z-forming"


def calculate_cg_content(sequence: str) -> float:
    """
    Calculate CG dinucleotide content percentage.
    
    | Parameter | Type | Description                      | Range        |
    |-----------|------|----------------------------------|--------------|
    | sequence  | str  | DNA sequence to analyze          | any          |
    | return    | float| CG dinucleotide content %        | 0.0 - 100.0  |
    """
    if len(sequence) < 2:
        return 0.0
    
    cg_count = 0
    total_dinucs = len(sequence) - 1
    
    for i in range(total_dinucs):
        if sequence[i:i+2] in ['CG', 'GC']:
            cg_count += 1
    
    return (cg_count / total_dinucs) * 100.0


def scan_sequence_internal(seq: str, min_len: int = 10, max_len: int = 30, step: int = 1) -> List[Dict[str, Any]]:
    """
    Internal function to scan sequence for Z-DNA motifs.
    
    | Parameter | Type | Description                    | Default | Range      |
    |-----------|------|--------------------------------|---------|------------|
    | seq       | str  | DNA sequence to scan           | -       | any        |
    | min_len   | int  | Minimum window length          | 10      | 10-50      |
    | max_len   | int  | Maximum window length          | 30      | min_len+   |
    | step      | int  | Step size for sliding window   | 1       | 1-5        |
    | return    | list | List of motif dictionaries     | []      | any        |
    """
    seq = seq.upper().replace('U', 'T')
    motifs = []
    
    # Scan with sliding windows
    for window_len in range(min_len, min(max_len + 1, len(seq) + 1)):
        for start in range(0, len(seq) - window_len + 1, step):
            window = seq[start:start + window_len]
            
            # Calculate score and characteristics
            score = calculate_zdna_score(window)
            cg_content = calculate_cg_content(window)
            classification = classify_zdna_region(score, window_len, cg_content)
            
            # Only keep regions with Z-DNA forming potential
            if score >= LOW_ZDNA * 0.8:  # Relaxed threshold for detection
                # Count number of scoring 10-mers
                scoring_nmers = sum(1 for i in range(len(window) - 9) 
                                  if window[i:i+10] in ZDNA_10MER_SCORES)
                
                motifs.append({
                    'Class': 'Z-DNA',
                    'Subclass': 'Z-DNA',
                    'Start': start + 1,  # 1-based coordinates
                    'End': start + window_len,
                    'Length': window_len,
                    'Sequence': window,
                    'Raw_Score': score,
                    'Score': score,
                    'Score_Method': '10mer_Z_Potential',
                    'N_Scoring_10mers': scoring_nmers,
                    'CG_Content': round(cg_content, 1),
                    'Classification': classification,
                    'Confidence': classification.split('-')[0]
                })
    
    return motifs


def find_z_dna(sequence: str, sequence_name: str = "Unknown") -> List[Dict[str, Any]]:
    """
    Main function to find Z-DNA motifs in a sequence.
    Compatible with NBDFinder motif detection framework.
    
    Uses 10-mer Z-DNA patterns with Z-forming potential scoring.
    
    | Parameter     | Type | Description                      | Default   | Range      |
    |---------------|------|----------------------------------|-----------|------------|
    | sequence      | str  | DNA sequence to analyze          | -         | >=10 bp    |
    | sequence_name | str  | Name/identifier for sequence     | "Unknown" | any        |
    | return        | list | Standardized Z-DNA motifs        | []        | any        |
    """
    if not sequence or len(sequence) < 10:
        return []
    
    # Clean sequence
    sequence = sequence.upper().replace('U', 'T')
    
    # Scan for Z-DNA motifs
    motifs = scan_sequence_internal(sequence, min_len=10, max_len=30, step=1)
    
    # Standardize output format for NBDFinder
    standardized_motifs = []
    for i, motif in enumerate(motifs):
        try:
            # Try to use standardize_motif_output if available
            motif['Sequence_Name'] = sequence_name
            motif['Motif_ID'] = f"Z_DNA_{i + 1}"
            standardized_motifs.append(motif)
        except Exception:
            # Fallback if standardization fails
            motif['Sequence_Name'] = sequence_name
            motif['Motif_ID'] = f"Z_DNA_{i + 1}"
            standardized_motifs.append(motif)
    
    return standardized_motifs


def scan_sequence(seq: str, min_len: int = 10, max_len: int = 30, step: int = 1,
                  out_tsv: str = "zdna_hits.tsv") -> str:
    """
    Legacy function for compatibility with Z-DNA scanner interface.
    Scans sequence and writes results to TSV file.
    
    | Parameter | Type | Description                    | Default           | Range      |
    |-----------|------|--------------------------------|-------------------|------------|
    | seq       | str  | DNA sequence to scan           | -                 | any        |
    | min_len   | int  | Minimum window length          | 10                | 10-50      |
    | max_len   | int  | Maximum window length          | 30                | min_len+   |
    | step      | int  | Step size for sliding window   | 1                 | 1-5        |
    | out_tsv   | str  | Output TSV filename            | "zdna_hits.tsv"   | any        |
    | return    | str  | Path to output TSV file        | out_tsv           | any        |
    """
    motifs = find_z_dna(seq, "scan_sequence")
    
    # Write to TSV file
    with open(out_tsv, "w", newline="") as fh:
        import csv
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["start", "end", "length", "n_scoring_10mers", "z_score", "cg_content", "label", "seq"])
        
        for motif in motifs:
            w.writerow([
                motif.get('Start', 0),
                motif.get('End', 0), 
                motif.get('Length', 0),
                motif.get('N_Scoring_10mers', 0),
                motif.get('Raw_Score', 0),
                motif.get('CG_Content', 0),
                motif.get('Confidence', 'Unknown'),
                motif.get('Sequence', '')
            ])
    
    return out_tsv


def merge_overlapping_regions(motifs: List[Dict[str, Any]], max_gap: int = 5) -> List[Dict[str, Any]]:
    """
    Merge overlapping or nearby Z-DNA regions to reduce redundancy.
    
    | Parameter | Type | Description                    | Default | Range      |
    |-----------|------|--------------------------------|---------|------------|
    | motifs    | list | List of Z-DNA motif dicts      | -       | any        |
    | max_gap   | int  | Maximum gap to merge regions   | 5       | 0-20       |
    | return    | list | Merged motif list              | []      | any        |
    """
    if not motifs:
        return []
    
    # Sort by start position
    sorted_motifs = sorted(motifs, key=lambda x: x['Start'])
    merged = []
    current = sorted_motifs[0].copy()
    
    for motif in sorted_motifs[1:]:
        # Check if motifs should be merged
        if motif['Start'] - current['End'] <= max_gap:
            # Merge motifs
            current['End'] = max(current['End'], motif['End'])
            current['Length'] = current['End'] - current['Start'] + 1
            current['Sequence'] = current['Sequence'] + motif['Sequence'][current['Length'] - motif['Length']:]
            current['Raw_Score'] = max(current['Raw_Score'], motif['Raw_Score'])
            current['Score'] = current['Raw_Score']
        else:
            # Add current to merged list and start new region
            merged.append(current)
            current = motif.copy()
    
    # Add the last motif
    merged.append(current)
    
    return merged


# === Example usage and testing ===
if __name__ == "__main__":
    # Test with Z-DNA forming sequence
    example = "NNNNACGCGCGCGCGCGCGCGCGCATCGCGCGCGNNNN"
    motifs = find_z_dna(example, "example_sequence")
    
    print(f"Found {len(motifs)} Z-DNA motifs:")
    for motif in motifs:
        score = motif.get('Raw_Score', motif.get('Score', 'N/A'))
        cg_content = motif.get('CG_Content', 'N/A')
        print(f"  {motif.get('Subclass', 'Unknown')}: {motif.get('Start', 0)}-{motif.get('End', 0)} (score: {score}, CG: {cg_content}%)")
        print(f"    Sequence: {motif.get('Sequence', 'N/A')}")
        print(f"    Classification: {motif.get('Classification', 'Unknown')}")
    
    # Also test legacy TSV output
    out_file = scan_sequence(example, min_len=10, max_len=25, step=1,
                           out_tsv="zdna_results.tsv")
    print(f"\nResults also written to: {out_file}")
    
    # Test merging functionality
    if motifs:
        merged = merge_overlapping_regions(motifs)
        print(f"\nAfter merging: {len(merged)} regions")