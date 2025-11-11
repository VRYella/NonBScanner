#!/usr/bin/env python3
"""
Genome-Scale Analysis Example
==============================
Demonstrates how to use the genome-scale scanner for large sequences.
"""

import time
import argparse
from genome_scale_scanner import analyze_genome_sequence
import pandas as pd


def read_fasta(filename):
    """Read FASTA file and return header and sequence"""
    print(f"Reading FASTA file: {filename}")
    with open(filename, 'r') as f:
        lines = f.readlines()
        header = lines[0].strip()[1:]  # Remove '>'
        sequence = ''.join(line.strip().upper() for line in lines[1:])
    return header, sequence


def export_results(motifs, output_prefix):
    """Export results to multiple formats"""
    print(f"\nExporting results with prefix: {output_prefix}")
    
    # Convert to DataFrame
    df = pd.DataFrame(motifs)
    
    # CSV export
    csv_file = f"{output_prefix}_motifs.csv"
    df.to_csv(csv_file, index=False)
    print(f"  ✓ Saved CSV: {csv_file}")
    
    # BED format export (for genome browsers)
    bed_file = f"{output_prefix}_motifs.bed"
    bed_df = df[['Sequence_Name', 'Start', 'End', 'ID', 'Score', 'Strand']].copy()
    bed_df['Start'] = bed_df['Start'] - 1  # BED is 0-based
    bed_df['Score'] = (bed_df['Score'] * 1000).astype(int)  # BED scores are 0-1000
    bed_df.to_csv(bed_file, sep='\t', header=False, index=False)
    print(f"  ✓ Saved BED: {bed_file}")
    
    # Class-specific exports
    for motif_class in df['Class'].unique():
        class_df = df[df['Class'] == motif_class]
        class_file = f"{output_prefix}_{motif_class}.csv"
        class_df.to_csv(class_file, index=False)
        print(f"  ✓ Saved {motif_class}: {class_file} ({len(class_df)} motifs)")
    
    # Summary statistics
    summary_file = f"{output_prefix}_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("NonBScanner Genome-Scale Analysis Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Sequence: {df['Sequence_Name'].iloc[0]}\n")
        f.write(f"Total motifs detected: {len(df)}\n\n")
        f.write("Class distribution:\n")
        for cls, count in df['Class'].value_counts().items():
            f.write(f"  {cls:25} {count:6} motifs\n")
        f.write("\n")
        f.write(f"Average score: {df['Score'].mean():.3f}\n")
        f.write(f"Average length: {df['Length'].mean():.1f} bp\n")
    print(f"  ✓ Saved summary: {summary_file}")


def analyze_and_export(fasta_file, output_prefix, chunk_size=1_000_000, 
                       enable_hybrid=False):
    """
    Main analysis workflow
    
    Args:
        fasta_file: Input FASTA file path
        output_prefix: Prefix for output files
        chunk_size: Chunk size for processing (bp)
        enable_hybrid: Whether to detect hybrid/cluster motifs
    """
    # Read input
    header, sequence = read_fasta(fasta_file)
    seq_len = len(sequence)
    
    print(f"\nSequence: {header}")
    print(f"Length: {seq_len:,} bp ({seq_len/1_000_000:.2f} MB)")
    
    # Analyze
    print("\nStarting genome-scale analysis...")
    start_time = time.time()
    
    motifs = analyze_genome_sequence(
        sequence, 
        header,
        chunk_size=chunk_size,
        enable_hybrid_cluster=enable_hybrid
    )
    
    analysis_time = time.time() - start_time
    
    # Export results
    if motifs:
        export_results(motifs, output_prefix)
    else:
        print("\n⚠ No motifs detected!")
    
    # Final summary
    print(f"\n{'='*80}")
    print("Analysis complete!")
    print(f"Total time: {analysis_time:.2f} seconds ({analysis_time/60:.2f} minutes)")
    print(f"Throughput: {seq_len/analysis_time:,.0f} bp/s")
    print(f"Motifs detected: {len(motifs)}")
    print('='*80)


def main():
    parser = argparse.ArgumentParser(
        description='Genome-scale Non-B DNA motif analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze a genome
  python genome_scale_example.py genome.fasta -o results/genome
  
  # Use larger chunks for faster processing
  python genome_scale_example.py large_genome.fasta -o results/large -c 5000000
  
  # Enable hybrid/cluster detection (slower)
  python genome_scale_example.py genome.fasta -o results/genome --hybrid
        """
    )
    
    parser.add_argument('input', help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, 
                       help='Output prefix for results files')
    parser.add_argument('-c', '--chunk-size', type=int, default=1_000_000,
                       help='Chunk size for processing (bp, default: 1000000)')
    parser.add_argument('--hybrid', action='store_true',
                       help='Enable hybrid/cluster motif detection (slower)')
    
    args = parser.parse_args()
    
    # Run analysis
    analyze_and_export(
        args.input, 
        args.output,
        chunk_size=args.chunk_size,
        enable_hybrid=args.hybrid
    )


if __name__ == "__main__":
    main()
