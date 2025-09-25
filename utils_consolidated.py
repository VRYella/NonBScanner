"""
NBDScanner Utilities - Helper Functions & I/O Operations
========================================================

Consolidated utility functions for sequence processing, I/O operations,
statistical calculations, and data formatting for the NBDScanner system.

UTILITY FUNCTIONS TABLE:
========================
Category          | Functions                    | Description
------------------|------------------------------|----------------------------------
Sequence I/O      | parse_fasta, write_fasta    | FASTA file handling
Sequence Utils    | reverse_complement, gc_content | Basic sequence operations  
Statistics        | get_basic_stats, motif_stats | Sequence and motif statistics
Formatting        | wrap, format_results        | Output formatting
Validation        | validate_sequence, quality_check | Input validation
Export            | to_bed, to_csv, to_json     | Data export formats

Author: Dr. Venkata Rajesh Yella
License: MIT
Version: 2024.1
"""

import re
import os
import json
import csv
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Union, Tuple
from collections import Counter, defaultdict
from io import StringIO
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# SEQUENCE I/O OPERATIONS
# =============================================================================

def parse_fasta(fasta_content: str) -> Dict[str, str]:
    """
    Parse FASTA format content into sequences dictionary
    
    Args:
        fasta_content: FASTA format string content
        
    Returns:
        Dictionary of {sequence_name: sequence}
    """
    sequences = {}
    current_name = None
    current_seq = []
    
    lines = fasta_content.strip().split('\n')
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        if line.startswith('>'):
            # Save previous sequence
            if current_name and current_seq:
                sequences[current_name] = ''.join(current_seq)
            
            # Start new sequence
            current_name = line[1:].strip()
            if not current_name:
                current_name = f"sequence_{len(sequences) + 1}"
            current_seq = []
        else:
            # Add to current sequence
            current_seq.append(line.upper())
    
    # Save last sequence
    if current_name and current_seq:
        sequences[current_name] = ''.join(current_seq)
    
    return sequences

def write_fasta(sequences: Dict[str, str], filename: str) -> bool:
    """
    Write sequences to FASTA format file
    
    Args:
        sequences: Dictionary of {name: sequence}
        filename: Output filename
        
    Returns:
        True if successful
    """
    try:
        with open(filename, 'w') as f:
            for name, seq in sequences.items():
                f.write(f">{name}\n")
                # Write sequence in 80-character lines
                for i in range(0, len(seq), 80):
                    f.write(f"{seq[i:i+80]}\n")
        return True
    except Exception as e:
        print(f"Error writing FASTA file {filename}: {e}")
        return False

def read_fasta_file(filename: str) -> Dict[str, str]:
    """
    Read FASTA file and return sequences dictionary
    
    Args:
        filename: Path to FASTA file
        
    Returns:
        Dictionary of {sequence_name: sequence}
    """
    try:
        with open(filename, 'r') as f:
            content = f.read()
        return parse_fasta(content)
    except Exception as e:
        print(f"Error reading FASTA file {filename}: {e}")
        return {}

# =============================================================================
# SEQUENCE MANIPULATION & VALIDATION
# =============================================================================

def validate_sequence(sequence: str) -> Tuple[bool, str]:
    """
    Validate DNA sequence format and content
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        (is_valid, error_message)
    """
    if not sequence:
        return False, "Empty sequence"
    
    if not isinstance(sequence, str):
        return False, "Sequence must be string"
    
    # Check for valid DNA characters
    valid_chars = set('ATGCRYSWKMBDHVN-')  # Include ambiguous bases
    invalid_chars = set(sequence.upper()) - valid_chars
    
    if invalid_chars:
        return False, f"Invalid characters found: {invalid_chars}"
    
    if len(sequence) < 10:
        return False, "Sequence too short (minimum 10 bp)"
    
    return True, "Valid sequence"

def reverse_complement(sequence: str) -> str:
    """
    Generate reverse complement of DNA sequence
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Reverse complement sequence
    """
    complement_map = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'R': 'Y', 'Y': 'R', 'S': 'W', 'W': 'S',
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
        'D': 'H', 'H': 'D', 'N': 'N', '-': '-'
    }
    
    return ''.join(complement_map.get(base.upper(), base) for base in reversed(sequence))

def gc_content(sequence: str) -> float:
    """
    Calculate GC content of sequence
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        GC content as percentage (0-100)
    """
    if not sequence:
        return 0.0
    
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100

def at_content(sequence: str) -> float:
    """
    Calculate AT content of sequence
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        AT content as percentage (0-100)
    """
    return 100.0 - gc_content(sequence)

def is_palindrome(sequence: str) -> bool:
    """
    Check if sequence is palindromic (same as reverse complement)
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        True if palindromic
    """
    return sequence.upper() == reverse_complement(sequence).upper()

def calculate_tm(sequence: str) -> float:
    """
    Calculate melting temperature using nearest neighbor method (simplified)
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Melting temperature in Celsius
    """
    if len(sequence) < 2:
        return 0.0
    
    # Simplified Tm calculation
    gc = gc_content(sequence)
    length = len(sequence)
    
    if length <= 13:
        # For short sequences
        tm = (sequence.upper().count('A') + sequence.upper().count('T')) * 2 + \
             (sequence.upper().count('G') + sequence.upper().count('C')) * 4
    else:
        # For longer sequences
        tm = 64.9 + 41 * (gc / 100) - 650 / length
    
    return tm

def shuffle_sequence(sequence: str, seed: Optional[int] = None) -> str:
    """
    Generate shuffled version of sequence (preserving composition)
    
    Args:
        sequence: DNA sequence string
        seed: Random seed for reproducibility
        
    Returns:
        Shuffled sequence
    """
    import random
    if seed is not None:
        random.seed(seed)
    
    seq_list = list(sequence.upper())
    random.shuffle(seq_list)
    return ''.join(seq_list)

# =============================================================================
# STATISTICS & ANALYSIS FUNCTIONS
# =============================================================================

def get_basic_stats(sequence: str, motifs: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
    """
    Calculate basic sequence statistics
    
    Args:
        sequence: DNA sequence string
        motifs: Optional list of detected motifs
        
    Returns:
        Dictionary of statistics
    """
    if not sequence:
        return {}
    
    seq = sequence.upper()
    length = len(seq)
    
    # Base composition
    base_counts = Counter(seq)
    
    stats = {
        'Length': length,
        'A': base_counts.get('A', 0),
        'T': base_counts.get('T', 0),
        'G': base_counts.get('G', 0),
        'C': base_counts.get('C', 0),
        'N': base_counts.get('N', 0),
        'GC%': round(gc_content(seq), 2),
        'AT%': round(at_content(seq), 2),
        'Tm': round(calculate_tm(seq), 1)
    }
    
    # Motif statistics if provided
    if motifs:
        stats.update(calculate_motif_statistics(motifs, length))
    
    return stats

def calculate_motif_statistics(motifs: List[Dict[str, Any]], sequence_length: int) -> Dict[str, Any]:
    """
    Calculate comprehensive motif statistics
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of analyzed sequence
        
    Returns:
        Dictionary of motif statistics
    """
    if not motifs:
        return {
            'Total_Motifs': 0,
            'Coverage%': 0.0,
            'Density': 0.0,
            'Classes_Detected': 0,
            'Subclasses_Detected': 0
        }
    
    # Count by class and subclass
    class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    subclass_counts = Counter(m.get('Subclass', 'Unknown') for m in motifs)
    
    # Calculate coverage
    covered_positions = set()
    for motif in motifs:
        start = motif.get('Start', 0) - 1  # Convert to 0-based
        end = motif.get('End', 0)
        covered_positions.update(range(start, end))
    
    coverage_percent = (len(covered_positions) / sequence_length * 100) if sequence_length > 0 else 0
    density = len(motifs) / (sequence_length / 1000) if sequence_length > 0 else 0  # Motifs per kb
    
    stats = {
        'Total_Motifs': len(motifs),
        'Coverage%': round(coverage_percent, 2),
        'Density': round(density, 2),
        'Classes_Detected': len(class_counts),
        'Subclasses_Detected': len(subclass_counts),
        'Class_Distribution': dict(class_counts),
        'Subclass_Distribution': dict(subclass_counts)
    }
    
    # Score statistics
    scores = [m.get('Score', 0) for m in motifs if isinstance(m.get('Score'), (int, float))]
    if scores:
        stats.update({
            'Score_Mean': round(np.mean(scores), 3),
            'Score_Std': round(np.std(scores), 3),
            'Score_Min': round(min(scores), 3),
            'Score_Max': round(max(scores), 3)
        })
    
    # Length statistics
    lengths = [m.get('Length', 0) for m in motifs if isinstance(m.get('Length'), int)]
    if lengths:
        stats.update({
            'Length_Mean': round(np.mean(lengths), 1),
            'Length_Std': round(np.std(lengths), 1),
            'Length_Min': min(lengths),
            'Length_Max': max(lengths)
        })
    
    return stats

# =============================================================================
# FORMATTING & OUTPUT UTILITIES
# =============================================================================

def wrap(sequence: str, width: int = 80) -> str:
    """
    Wrap sequence to specified width
    
    Args:
        sequence: DNA sequence string
        width: Line width for wrapping
        
    Returns:
        Wrapped sequence string
    """
    if not sequence or width <= 0:
        return sequence
    
    return '\n'.join(sequence[i:i+width] for i in range(0, len(sequence), width))

def format_motif_rows(motifs: List[Dict[str, Any]]) -> List[List[str]]:
    """
    Format motifs for tabular display
    
    Args:
        motifs: List of motif dictionaries
        
    Returns:
        List of formatted rows
    """
    if not motifs:
        return []
    
    rows = []
    headers = ['Class', 'Subclass', 'Start', 'End', 'Length', 'Sequence', 'Score', 'Strand']
    
    for motif in motifs:
        row = []
        for header in headers:
            value = motif.get(header, 'N/A')
            if header == 'Sequence' and isinstance(value, str) and len(value) > 50:
                value = value[:47] + '...'
            elif header == 'Score' and isinstance(value, (int, float)):
                value = f"{value:.3f}"
            row.append(str(value))
        rows.append(row)
    
    return rows

def create_summary_table(sequences: Dict[str, str], results: Dict[str, List[Dict[str, Any]]]) -> pd.DataFrame:
    """
    Create summary table for multiple sequence analysis
    
    Args:
        sequences: Dictionary of {name: sequence}
        results: Dictionary of {name: motifs_list}
        
    Returns:
        Summary DataFrame
    """
    summary_data = []
    
    for name, sequence in sequences.items():
        motifs = results.get(name, [])
        stats = get_basic_stats(sequence, motifs)
        
        summary_data.append({
            'Sequence_Name': name,
            'Length_bp': stats.get('Length', 0),
            'GC_Content': stats.get('GC%', 0),
            'Total_Motifs': stats.get('Total_Motifs', 0),
            'Coverage_Percent': stats.get('Coverage%', 0),
            'Motif_Density': stats.get('Density', 0),
            'Classes_Detected': stats.get('Classes_Detected', 0),
            'Subclasses_Detected': stats.get('Subclasses_Detected', 0)
        })
    
    return pd.DataFrame(summary_data)

# =============================================================================
# DATA EXPORT FUNCTIONS
# =============================================================================

def export_to_bed(motifs: List[Dict[str, Any]], sequence_name: str = "sequence", 
                  filename: Optional[str] = None) -> str:
    """
    Export motifs to BED format
    
    Args:
        motifs: List of motif dictionaries
        sequence_name: Name of the sequence
        filename: Optional output filename
        
    Returns:
        BED format string
    """
    bed_lines = []
    bed_lines.append("track name=NBDScanner_motifs description=\"Non-B DNA motifs\" itemRgb=On")
    
    # Color mapping for different classes
    class_colors = {
        'Curved_DNA': '255,182,193',      # Light pink
        'Slipped_DNA': '255,218,185',     # Peach
        'Cruciform': '173,216,230',       # Light blue
        'R-Loop': '144,238,144',          # Light green
        'Triplex': '221,160,221',         # Plum
        'G-Quadruplex': '255,215,0',      # Gold
        'i-Motif': '255,165,0',           # Orange
        'Z-DNA': '138,43,226',            # Blue violet
        'A-philic_DNA': '230,230,250',    # Lavender
        'Hybrid': '192,192,192',          # Silver
        'Non-B_DNA_Clusters': '128,128,128'  # Gray
    }
    
    for motif in motifs:
        chrom = sequence_name
        start = max(0, motif.get('Start', 1) - 1)  # Convert to 0-based
        end = motif.get('End', start + 1)
        name = f"{motif.get('Class', 'Unknown')}_{motif.get('Subclass', 'Unknown')}"
        score = int(min(1000, max(0, motif.get('Score', 0) * 1000)))  # Scale to 0-1000
        strand = motif.get('Strand', '+')
        color = class_colors.get(motif.get('Class'), '128,128,128')
        
        bed_line = f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{start}\t{end}\t{color}"
        bed_lines.append(bed_line)
    
    bed_content = '\n'.join(bed_lines)
    
    if filename:
        try:
            with open(filename, 'w') as f:
                f.write(bed_content)
        except Exception as e:
            print(f"Error writing BED file {filename}: {e}")
    
    return bed_content

def export_to_csv(motifs: List[Dict[str, Any]], filename: Optional[str] = None) -> str:
    """
    Export motifs to CSV format
    
    Args:
        motifs: List of motif dictionaries
        filename: Optional output filename
        
    Returns:
        CSV format string
    """
    if not motifs:
        return "No motifs to export"
    
    # Get all unique keys from motifs
    all_keys = set()
    for motif in motifs:
        all_keys.update(motif.keys())
    
    # Standard column order
    standard_columns = ['ID', 'Sequence_Name', 'Class', 'Subclass', 'Start', 'End', 
                       'Length', 'Sequence', 'Score', 'Strand', 'Method']
    
    # Arrange columns with standard ones first
    columns = [col for col in standard_columns if col in all_keys]
    columns.extend([col for col in sorted(all_keys) if col not in standard_columns])
    
    output = StringIO()
    writer = csv.DictWriter(output, fieldnames=columns)
    writer.writeheader()
    
    for motif in motifs:
        # Ensure all columns have values
        row = {col: motif.get(col, '') for col in columns}
        writer.writerow(row)
    
    csv_content = output.getvalue()
    output.close()
    
    if filename:
        try:
            with open(filename, 'w', newline='') as f:
                f.write(csv_content)
        except Exception as e:
            print(f"Error writing CSV file {filename}: {e}")
    
    return csv_content

def export_to_json(motifs: List[Dict[str, Any]], filename: Optional[str] = None, 
                   pretty: bool = True) -> str:
    """
    Export motifs to JSON format
    
    Args:
        motifs: List of motif dictionaries
        filename: Optional output filename
        pretty: Whether to format JSON prettily
        
    Returns:
        JSON format string
    """
    json_data = {
        'version': '2024.1',
        'analysis_type': 'NBDScanner_Non-B_DNA_Analysis',
        'total_motifs': len(motifs),
        'motifs': motifs
    }
    
    if pretty:
        json_content = json.dumps(json_data, indent=2, ensure_ascii=False)
    else:
        json_content = json.dumps(json_data, ensure_ascii=False)
    
    if filename:
        try:
            with open(filename, 'w') as f:
                f.write(json_content)
        except Exception as e:
            print(f"Error writing JSON file {filename}: {e}")
    
    return json_content

def export_to_gff3(motifs: List[Dict[str, Any]], sequence_name: str = "sequence", 
                   filename: Optional[str] = None) -> str:
    """
    Export motifs to GFF3 format
    
    Args:
        motifs: List of motif dictionaries
        sequence_name: Name of the sequence
        filename: Optional output filename
        
    Returns:
        GFF3 format string
    """
    gff_lines = []
    gff_lines.append("##gff-version 3")
    gff_lines.append(f"##sequence-region {sequence_name} 1 {len(sequence_name)}")
    
    for i, motif in enumerate(motifs, 1):
        seqid = sequence_name
        source = "NBDScanner"
        feature_type = "Non_B_DNA_motif"
        start = motif.get('Start', 1)
        end = motif.get('End', start)
        score = motif.get('Score', '.')
        strand = motif.get('Strand', '+')
        phase = "."
        
        # Attributes
        attributes = [
            f"ID=motif_{i}",
            f"Name={motif.get('Class', 'Unknown')}_{motif.get('Subclass', 'Unknown')}",
            f"motif_class={motif.get('Class', 'Unknown')}",
            f"motif_subclass={motif.get('Subclass', 'Unknown')}",
            f"length={motif.get('Length', 0)}",
            f"method={motif.get('Method', 'NBDScanner')}"
        ]
        
        attributes_str = ';'.join(attributes)
        
        gff_line = f"{seqid}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes_str}"
        gff_lines.append(gff_line)
    
    gff_content = '\n'.join(gff_lines)
    
    if filename:
        try:
            with open(filename, 'w') as f:
                f.write(gff_content)
        except Exception as e:
            print(f"Error writing GFF3 file {filename}: {e}")
    
    return gff_content

# =============================================================================
# QUALITY CONTROL & FILTERING
# =============================================================================

def quality_check_motifs(motifs: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Perform quality checks on detected motifs
    
    Args:
        motifs: List of motif dictionaries
        
    Returns:
        Quality check report
    """
    if not motifs:
        return {'status': 'No motifs to check', 'passed': True, 'issues': []}
    
    issues = []
    
    # Check for required fields
    required_fields = ['Class', 'Subclass', 'Start', 'End', 'Sequence']
    for i, motif in enumerate(motifs):
        missing_fields = [field for field in required_fields if field not in motif]
        if missing_fields:
            issues.append(f"Motif {i+1}: Missing fields {missing_fields}")
    
    # Check coordinate consistency
    for i, motif in enumerate(motifs):
        start = motif.get('Start')
        end = motif.get('End')
        length = motif.get('Length')
        sequence = motif.get('Sequence', '')
        
        if start and end and start >= end:
            issues.append(f"Motif {i+1}: Invalid coordinates (start >= end)")
        
        if start and end and length and (end - start + 1) != length:
            issues.append(f"Motif {i+1}: Length inconsistent with coordinates")
        
        if sequence and length and len(sequence) != length:
            issues.append(f"Motif {i+1}: Sequence length doesn't match reported length")
    
    # Check for overlaps within same class
    class_groups = defaultdict(list)
    for motif in motifs:
        class_groups[motif.get('Class')].append(motif)
    
    for class_name, class_motifs in class_groups.items():
        for i, motif1 in enumerate(class_motifs):
            for motif2 in class_motifs[i+1:]:
                if (motif1.get('Start', 0) < motif2.get('End', 0) and 
                    motif2.get('Start', 0) < motif1.get('End', 0)):
                    issues.append(f"Overlapping motifs in class {class_name}")
                    break
    
    report = {
        'total_motifs': len(motifs),
        'issues_found': len(issues),
        'passed': len(issues) == 0,
        'issues': issues[:10],  # Limit to first 10 issues
        'status': 'PASSED' if len(issues) == 0 else f'FAILED ({len(issues)} issues)'
    }
    
    return report

def filter_motifs_by_score(motifs: List[Dict[str, Any]], min_score: float = 0.0) -> List[Dict[str, Any]]:
    """
    Filter motifs by minimum score threshold
    
    Args:
        motifs: List of motif dictionaries
        min_score: Minimum score threshold
        
    Returns:
        Filtered motifs list
    """
    return [m for m in motifs if m.get('Score', 0) >= min_score]

def filter_motifs_by_length(motifs: List[Dict[str, Any]], 
                           min_length: int = 0, max_length: int = float('inf')) -> List[Dict[str, Any]]:
    """
    Filter motifs by length range
    
    Args:
        motifs: List of motif dictionaries
        min_length: Minimum length
        max_length: Maximum length
        
    Returns:
        Filtered motifs list
    """
    return [m for m in motifs if min_length <= m.get('Length', 0) <= max_length]

def filter_motifs_by_class(motifs: List[Dict[str, Any]], 
                          allowed_classes: List[str]) -> List[Dict[str, Any]]:
    """
    Filter motifs by allowed classes
    
    Args:
        motifs: List of motif dictionaries
        allowed_classes: List of allowed class names
        
    Returns:
        Filtered motifs list
    """
    return [m for m in motifs if m.get('Class') in allowed_classes]

# =============================================================================
# TESTING & EXAMPLES
# =============================================================================

def test_utilities():
    """Test utility functions with example data"""
    print("Testing NBDScanner utilities...")
    
    # Test sequence validation
    test_sequences = [
        "ATGCATGCATGC",     # Valid
        "ATGCXYZ",          # Invalid characters
        "ATG",              # Too short
        "",                 # Empty
    ]
    
    print("\nSequence validation tests:")
    for seq in test_sequences:
        valid, msg = validate_sequence(seq)
        print(f"  '{seq[:20]}': {'VALID' if valid else 'INVALID'} - {msg}")
    
    # Test basic statistics
    test_seq = "GGGTTAGGGTTAGGGTTAGGGAAAAATTTTCGCGCGCGCG"
    stats = get_basic_stats(test_seq)
    print(f"\nBasic stats for test sequence:")
    for key, value in stats.items():
        print(f"  {key}: {value}")
    
    # Test FASTA parsing
    fasta_content = """>test_sequence
GGGTTAGGGTTAGGGTTAGGG
>another_sequence
AAAAATTTTCCCCGGGG"""
    
    sequences = parse_fasta(fasta_content)
    print(f"\nFASTA parsing test: {len(sequences)} sequences parsed")
    for name, seq in sequences.items():
        print(f"  {name}: {seq}")
    
    # Test motif formatting
    example_motifs = [
        {
            'Class': 'G-Quadruplex',
            'Subclass': 'Canonical G4',
            'Start': 1,
            'End': 21,
            'Length': 21,
            'Sequence': 'GGGTTAGGGTTAGGGTTAGGG',
            'Score': 0.857,
            'Strand': '+'
        }
    ]
    
    quality_report = quality_check_motifs(example_motifs)
    print(f"\nQuality check: {quality_report['status']}")
    
    print("✓ All utility tests completed")

if __name__ == "__main__":
    test_utilities()