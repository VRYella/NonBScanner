"""
FASTA Slicing Utility
====================

Standalone utility for extracting regions from FASTA files with
windowing and overlap capabilities.
"""

import argparse
import sys
from pathlib import Path
from typing import List, Tuple, Dict
import re

from ..io.fasta import FastaReader, write_fasta

def main():
    """Main entry point for FASTA slicing."""
    parser = argparse.ArgumentParser(
        description="Extract regions from FASTA files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract specific regions
  slice-fasta --fasta genome.fa --regions chr1:1000-2000,chr2:5000-6000
  
  # Create sliding windows
  slice-fasta --fasta genome.fa --windows --window-size 10000 --overlap 1000
  
  # Extract around motif coordinates
  slice-fasta --fasta genome.fa --motifs motifs.bed --flank 500
        """
    )
    
    # Input/Output
    parser.add_argument('--fasta', '-f', required=True, type=Path,
                       help='Input FASTA file')
    parser.add_argument('--output', '-o', type=Path, default=Path('sliced.fasta'),
                       help='Output FASTA file')
    
    # Region specification methods (mutually exclusive)
    region_group = parser.add_mutually_exclusive_group(required=True)
    region_group.add_argument('--regions', '-r',
                             help='Comma-separated regions (chr:start-end)')
    region_group.add_argument('--windows', action='store_true',
                             help='Create sliding windows')
    region_group.add_argument('--motifs', type=Path,
                             help='BED file with motif coordinates')
    
    # Window parameters
    parser.add_argument('--window-size', type=int, default=10000,
                       help='Window size for sliding windows')
    parser.add_argument('--overlap', type=int, default=1000,
                       help='Overlap between windows')
    parser.add_argument('--step-size', type=int,
                       help='Step size (alternative to overlap)')
    
    # Motif flanking
    parser.add_argument('--flank', type=int, default=0,
                       help='Flanking sequence around motifs')
    parser.add_argument('--upstream', type=int,
                       help='Upstream flanking (overrides --flank)')
    parser.add_argument('--downstream', type=int,
                       help='Downstream flanking (overrides --flank)')
    
    # Filtering
    parser.add_argument('--min-length', type=int, default=1,
                       help='Minimum sequence length')
    parser.add_argument('--max-length', type=int,
                       help='Maximum sequence length')
    parser.add_argument('--sequences', 
                       help='Comma-separated sequence names to process')
    
    # Output options
    parser.add_argument('--naming', choices=['coordinates', 'sequential', 'original'],
                       default='coordinates',
                       help='Naming scheme for extracted sequences')
    parser.add_argument('--line-width', type=int, default=80,
                       help='Line width for FASTA output')
    
    # Verbosity
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    try:
        return slice_fasta_main(args)
    except Exception as e:
        print(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1

def slice_fasta_main(args) -> int:
    """Main slicing logic."""
    if not args.fasta.exists():
        print(f"Error: FASTA file not found: {args.fasta}")
        return 1
    
    # Determine sequences to process
    target_sequences = None
    if args.sequences:
        target_sequences = [seq.strip() for seq in args.sequences.split(',')]
    
    extracted_sequences = {}
    
    with FastaReader(args.fasta) as reader:
        sequence_names = reader.get_sequence_names()
        
        # Filter sequences if specified
        if target_sequences:
            sequence_names = [name for name in sequence_names if name in target_sequences]
        
        if args.verbose:
            print(f"Processing {len(sequence_names)} sequences")
        
        for seq_name in sequence_names:
            if args.verbose:
                print(f"Processing sequence: {seq_name}")
            
            if args.regions:
                # Extract specific regions
                regions = parse_regions(args.regions, seq_name)
                for region_name, start, end in regions:
                    try:
                        sequence = reader.get_sequence_slice(seq_name, start, end)
                        if _passes_length_filter(sequence, args):
                            extracted_sequences[region_name] = sequence
                    except Exception as e:
                        print(f"Warning: Could not extract {region_name}: {e}")
            
            elif args.windows:
                # Create sliding windows
                step_size = args.step_size if args.step_size else (args.window_size - args.overlap)
                
                for start, end, window_seq in reader.create_windows(
                    seq_name, args.window_size, args.overlap):
                    
                    if _passes_length_filter(window_seq, args):
                        window_name = _generate_window_name(seq_name, start, end, args.naming)
                        extracted_sequences[window_name] = window_seq
            
            elif args.motifs:
                # Extract around motif coordinates
                motif_regions = parse_bed_file(args.motifs, seq_name)
                
                upstream = args.upstream if args.upstream is not None else args.flank
                downstream = args.downstream if args.downstream is not None else args.flank
                
                for motif_name, motif_start, motif_end in motif_regions:
                    extract_start = max(0, motif_start - upstream)
                    extract_end = motif_end + downstream
                    
                    try:
                        sequence = reader.get_sequence_slice(seq_name, extract_start, extract_end)
                        if _passes_length_filter(sequence, args):
                            region_name = f"{motif_name}_flank{upstream}_{downstream}"
                            extracted_sequences[region_name] = sequence
                    except Exception as e:
                        print(f"Warning: Could not extract around {motif_name}: {e}")
    
    if not extracted_sequences:
        print("No sequences extracted")
        return 1
    
    # Write output
    write_fasta(extracted_sequences, args.output, args.line_width)
    
    if args.verbose:
        print(f"Extracted {len(extracted_sequences)} sequences to {args.output}")
        
        # Length statistics
        lengths = [len(seq) for seq in extracted_sequences.values()]
        print(f"Length range: {min(lengths)}-{max(lengths)} bp")
        print(f"Mean length: {sum(lengths)/len(lengths):.1f} bp")
    
    return 0

def parse_regions(regions_str: str, default_seq: str = None) -> List[Tuple[str, int, int]]:
    """Parse region specifications."""
    regions = []
    
    for region_str in regions_str.split(','):
        region_str = region_str.strip()
        
        if ':' in region_str and '-' in region_str:
            # Format: chr:start-end
            seq_name, coords = region_str.split(':', 1)
            start_str, end_str = coords.split('-', 1)
            start = int(start_str) - 1  # Convert to 0-based
            end = int(end_str)
            region_name = region_str
        
        elif '-' in region_str and default_seq:
            # Format: start-end (use default sequence)
            start_str, end_str = region_str.split('-', 1)
            start = int(start_str) - 1  # Convert to 0-based
            end = int(end_str)
            seq_name = default_seq
            region_name = f"{seq_name}:{start+1}-{end}"
        
        else:
            raise ValueError(f"Invalid region format: {region_str}")
        
        regions.append((region_name, start, end))
    
    return regions

def parse_bed_file(bed_path: Path, target_seq: str = None) -> List[Tuple[str, int, int]]:
    """Parse BED file for motif coordinates."""
    regions = []
    
    with open(bed_path) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('track'):
                continue
            
            try:
                fields = line.split('\t')
                if len(fields) < 3:
                    continue
                
                seq_name = fields[0]
                start = int(fields[1])  # BED is 0-based
                end = int(fields[2])
                
                # Use name field if available, otherwise generate
                if len(fields) >= 4:
                    motif_name = fields[3]
                else:
                    motif_name = f"motif_{line_num}"
                
                # Filter by target sequence if specified
                if target_seq is None or seq_name == target_seq:
                    regions.append((motif_name, start, end))
                
            except (ValueError, IndexError) as e:
                print(f"Warning: Skipping invalid BED line {line_num}: {e}")
    
    return regions

def _passes_length_filter(sequence: str, args) -> bool:
    """Check if sequence passes length filters."""
    length = len(sequence)
    
    if length < args.min_length:
        return False
    
    if args.max_length and length > args.max_length:
        return False
    
    return True

def _generate_window_name(seq_name: str, start: int, end: int, naming: str) -> str:
    """Generate name for window sequence."""
    if naming == 'coordinates':
        return f"{seq_name}:{start+1}-{end}"
    elif naming == 'sequential':
        # This would need a counter - simplified here
        return f"{seq_name}_window_{start//1000}k"
    elif naming == 'original':
        return seq_name
    else:
        return f"{seq_name}:{start+1}-{end}"

if __name__ == '__main__':
    sys.exit(main())