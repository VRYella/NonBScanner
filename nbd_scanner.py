#!/usr/bin/env python3
"""
NonBScanner - Non-B DNA Structure Detection Tool
- Command-line interface for A-philic DNA and Z-DNA detection
- Processes FASTA files and individual sequences
- Outputs results in multiple formats (TSV, JSON, summary)
- Supports batch processing and filtering options

This is the main CLI entry point for the consolidated NonBScanner tool,
providing unified access to A-philic DNA and Z-DNA detection capabilities.

Usage Examples:
  python nbd_scanner.py --sequence "ATCGATCGATCG" --detect aphilic
  python nbd_scanner.py --fasta input.fa --detect zdna --output results.tsv
  python nbd_scanner.py --fasta input.fa --detect both --format json

Author: Dr. Venkata Rajesh Yella
Integration: 2024 - NBDFinder Consolidated CLI
"""

import argparse
import sys
import os
import json
from typing import List, Dict, Any, Optional
from pathlib import Path

# Import our scanners and utilities
from aphilic_scanner import find_a_philic_dna, scan_sequence as scan_aphilic
from zdna_scanner import find_z_dna, scan_sequence as scan_zdna
from motif_utils import (
    parse_fasta, merge_overlapping_motifs, filter_motifs_by_quality,
    calculate_motif_density, get_motif_distribution, standardize_motif_output
)


def parse_fasta_file(filename: str) -> List[Dict[str, str]]:
    """
    Parse FASTA file and return list of sequences.
    
    | Parameter | Type | Description                    | Range      |
    |-----------|------|--------------------------------|------------|
    | filename  | str  | Path to FASTA file             | any        |
    | return    | list | List of {name, sequence} dicts | any        |
    """
    sequences = []
    current_name = ""
    current_seq = ""
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_name:
                        sequences.append({
                            'name': current_name,
                            'sequence': current_seq.upper().replace('U', 'T')
                        })
                    current_name = line[1:].split()[0]  # Take first word after >
                    current_seq = ""
                elif line:
                    current_seq += line
            
            # Add last sequence
            if current_name:
                sequences.append({
                    'name': current_name,
                    'sequence': current_seq.upper().replace('U', 'T')
                })
    
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file '{filename}': {e}")
        sys.exit(1)
    
    return sequences


def detect_motifs(sequence: str, sequence_name: str, detect_type: str, **kwargs) -> List[Dict[str, Any]]:
    """
    Detect motifs in sequence based on specified type.
    
    | Parameter     | Type | Description                    | Default | Range           |
    |---------------|------|--------------------------------|---------|-----------------|
    | sequence      | str  | DNA sequence to analyze        | -       | any             |
    | sequence_name | str  | Name/ID for sequence           | -       | any             |
    | detect_type   | str  | Type of detection to perform   | -       | aphilic/zdna/both|
    | return        | list | List of detected motifs        | []      | any             |
    """
    motifs = []
    
    if detect_type in ['aphilic', 'both']:
        try:
            aphilic_motifs = find_a_philic_dna(sequence, sequence_name)
            motifs.extend(aphilic_motifs)
        except Exception as e:
            print(f"Warning: A-philic detection failed for {sequence_name}: {e}")
    
    if detect_type in ['zdna', 'both']:
        try:
            zdna_motifs = find_z_dna(sequence, sequence_name)
            motifs.extend(zdna_motifs)
        except Exception as e:
            print(f"Warning: Z-DNA detection failed for {sequence_name}: {e}")
    
    return motifs


def format_output_tsv(motifs: List[Dict[str, Any]], output_file: str):
    """
    Write motifs to TSV format file.
    
    | Parameter    | Type | Description                    | Range      |
    |--------------|------|--------------------------------|------------|
    | motifs       | list | List of motif dictionaries    | any        |
    | output_file  | str  | Output file path               | any        |
    """
    import csv
    
    with open(output_file, 'w', newline='') as f:
        if not motifs:
            return
            
        # Get all possible field names
        all_fields = set()
        for motif in motifs:
            all_fields.update(motif.keys())
        
        # Standard field order
        standard_fields = ['Sequence_Name', 'Class', 'Subclass', 'Start', 'End', 'Length', 
                          'Sequence', 'Raw_Score', 'Score', 'Classification', 'Confidence']
        
        # Order fields - standard first, then others
        ordered_fields = []
        for field in standard_fields:
            if field in all_fields:
                ordered_fields.append(field)
                all_fields.remove(field)
        ordered_fields.extend(sorted(all_fields))
        
        writer = csv.DictWriter(f, fieldnames=ordered_fields, delimiter='\t')
        writer.writeheader()
        
        for motif in motifs:
            writer.writerow(motif)


def format_output_json(motifs: List[Dict[str, Any]], output_file: str):
    """
    Write motifs to JSON format file.
    
    | Parameter    | Type | Description                    | Range      |
    |--------------|------|--------------------------------|------------|
    | motifs       | list | List of motif dictionaries    | any        |
    | output_file  | str  | Output file path               | any        |
    """
    with open(output_file, 'w') as f:
        json.dump(motifs, f, indent=2, default=str)


def format_output_summary(motifs: List[Dict[str, Any]], sequence_info: List[Dict], output_file: str):
    """
    Write summary report to text file.
    
    | Parameter     | Type | Description                    | Range      |
    |---------------|------|--------------------------------|------------|
    | motifs        | list | List of motif dictionaries    | any        |
    | sequence_info | list | List of sequence info dicts   | any        |
    | output_file   | str  | Output file path               | any        |
    """
    with open(output_file, 'w') as f:
        f.write("NonBScanner Analysis Summary\n")
        f.write("=" * 50 + "\n\n")
        
        # Overall statistics
        total_sequences = len(sequence_info)
        total_length = sum(len(seq['sequence']) for seq in sequence_info)
        total_motifs = len(motifs)
        
        f.write(f"Sequences analyzed: {total_sequences}\n")
        f.write(f"Total sequence length: {total_length:,} bp\n")
        f.write(f"Total motifs found: {total_motifs}\n")
        
        if total_length > 0:
            density = calculate_motif_density(motifs, total_length)
            f.write(f"Motif density: {density:.2f} motifs/kb\n")
        
        f.write("\n")
        
        # Motif distribution
        if motifs:
            distribution = get_motif_distribution(motifs)
            f.write("Motif Distribution:\n")
            f.write("-" * 20 + "\n")
            for motif_type, count in sorted(distribution.items()):
                percentage = (count / total_motifs) * 100
                f.write(f"{motif_type}: {count} ({percentage:.1f}%)\n")
            f.write("\n")
        
        # Per-sequence statistics
        f.write("Per-Sequence Results:\n")
        f.write("-" * 21 + "\n")
        
        for seq_info in sequence_info:
            seq_name = seq_info['name']
            seq_length = len(seq_info['sequence'])
            seq_motifs = [m for m in motifs if m.get('Sequence_Name', '') == seq_name]
            
            f.write(f"{seq_name}:\n")
            f.write(f"  Length: {seq_length:,} bp\n")
            f.write(f"  Motifs: {len(seq_motifs)}\n")
            
            if seq_motifs and seq_length > 0:
                seq_density = calculate_motif_density(seq_motifs, seq_length)
                f.write(f"  Density: {seq_density:.2f} motifs/kb\n")
                
                seq_distribution = get_motif_distribution(seq_motifs)
                for motif_type, count in seq_distribution.items():
                    f.write(f"    {motif_type}: {count}\n")
            
            f.write("\n")


def main():
    """
    Main CLI function with argument parsing and execution logic.
    
    | Feature         | Type | Description                          | Options                    |
    |-----------------|------|--------------------------------------|----------------------------|
    | Input           | str  | Sequence input method                | --sequence, --fasta        |
    | Detection       | str  | Type of motifs to detect             | aphilic, zdna, both        |
    | Output format   | str  | Output format type                   | tsv, json, summary         |
    | Filtering       | bool | Enable quality filtering             | --filter                   |
    | Merging         | bool | Merge overlapping regions            | --merge                    |
    """
    parser = argparse.ArgumentParser(
        description='NonBScanner - Detect A-philic DNA and Z-DNA motifs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Scan single sequence for A-philic DNA
  python nbd_scanner.py --sequence "ATCGATCGATCG" --detect aphilic
  
  # Scan FASTA file for Z-DNA and output to TSV
  python nbd_scanner.py --fasta input.fa --detect zdna --output results.tsv
  
  # Scan for both motif types with filtering and merging
  python nbd_scanner.py --fasta input.fa --detect both --filter --merge
  
  # Output detailed summary report
  python nbd_scanner.py --fasta input.fa --detect both --format summary
        """
    )
    
    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--sequence', '-s', 
                           help='Single DNA sequence to analyze')
    input_group.add_argument('--fasta', '-f', 
                           help='FASTA file containing sequences to analyze')
    
    # Detection options
    parser.add_argument('--detect', '-d', 
                       choices=['aphilic', 'zdna', 'both'], 
                       default='both',
                       help='Type of motifs to detect (default: both)')
    
    # Output options
    parser.add_argument('--output', '-o', 
                       help='Output file name (default: stdout or auto-generated)')
    parser.add_argument('--format', '-fmt', 
                       choices=['tsv', 'json', 'summary'], 
                       default='tsv',
                       help='Output format (default: tsv)')
    
    # Processing options
    parser.add_argument('--filter', action='store_true',
                       help='Apply quality filtering to results')
    parser.add_argument('--merge', action='store_true',
                       help='Merge overlapping or nearby motifs')
    parser.add_argument('--min-score', type=float, default=0.0,
                       help='Minimum score threshold for filtering')
    parser.add_argument('--min-length', type=int, default=10,
                       help='Minimum motif length for filtering')
    
    # Scanning parameters  
    parser.add_argument('--max-length', type=int, default=30,
                       help='Maximum motif length to scan')
    parser.add_argument('--step-size', type=int, default=1,
                       help='Step size for sliding window')
    
    # Other options
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose output')
    parser.add_argument('--version', action='version', version='NonBScanner 1.0')
    
    args = parser.parse_args()
    
    # Prepare sequences
    sequences = []
    if args.sequence:
        sequences = [{'name': 'input_sequence', 'sequence': args.sequence.upper().replace('U', 'T')}]
    else:
        sequences = parse_fasta_file(args.fasta)
    
    if not sequences:
        print("Error: No sequences to analyze")
        sys.exit(1)
    
    if args.verbose:
        print(f"Analyzing {len(sequences)} sequence(s) for {args.detect} motifs...")
    
    # Process sequences
    all_motifs = []
    for seq_info in sequences:
        if args.verbose:
            print(f"Processing {seq_info['name']} ({len(seq_info['sequence'])} bp)...")
        
        motifs = detect_motifs(
            seq_info['sequence'], 
            seq_info['name'], 
            args.detect,
            max_length=args.max_length,
            step_size=args.step_size
        )
        
        all_motifs.extend(motifs)
    
    # Apply filtering if requested
    if args.filter:
        if args.verbose:
            print(f"Applying quality filters (min_score={args.min_score}, min_length={args.min_length})...")
        
        all_motifs = filter_motifs_by_quality(
            all_motifs,
            min_score=args.min_score,
            min_length=args.min_length
        )
    
    # Apply merging if requested
    if args.merge:
        if args.verbose:
            print("Merging overlapping regions...")
        
        all_motifs = merge_overlapping_motifs(all_motifs)
    
    # Generate output filename if not provided
    output_file = args.output
    if not output_file:
        if args.sequence:
            base_name = "sequence_results"
        else:
            base_name = Path(args.fasta).stem + "_results"
        
        extensions = {'tsv': '.tsv', 'json': '.json', 'summary': '.txt'}
        output_file = base_name + extensions[args.format]
    
    # Output results
    if args.verbose:
        print(f"Found {len(all_motifs)} motifs total")
        print(f"Writing results to {output_file} in {args.format} format...")
    
    if args.format == 'tsv':
        format_output_tsv(all_motifs, output_file)
    elif args.format == 'json':
        format_output_json(all_motifs, output_file)
    elif args.format == 'summary':
        format_output_summary(all_motifs, sequences, output_file)
    
    # Console output for single sequence
    if args.sequence and not args.verbose:
        print(f"Found {len(all_motifs)} motifs in sequence")
        for motif in all_motifs[:10]:  # Show first 10
            print(f"  {motif.get('Class', 'Unknown')}: {motif.get('Start', 0)}-{motif.get('End', 0)} "
                  f"(score: {motif.get('Score', 0):.2f})")
        if len(all_motifs) > 10:
            print(f"  ... and {len(all_motifs) - 10} more (see output file for details)")
    
    if args.verbose:
        print("Analysis complete!")


if __name__ == "__main__":
    main()