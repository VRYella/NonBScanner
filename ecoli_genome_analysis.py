#!/usr/bin/env python3
"""
E. coli K-12 MG1655 Genome Analysis with NBDScanner
=====================================================

This script downloads the complete E. coli K-12 MG1655 genome from NCBI
and performs comprehensive Non-B DNA motif analysis.

Author: Dr. Venkata Rajesh Yella
Date: 2024
"""

import os
import sys
import time
import json
from pathlib import Path
from collections import Counter, defaultdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Import NBDScanner modules
from utils.nbdscanner import analyze_sequence, get_motif_classification_info
from utils.utils import parse_fasta, get_basic_stats, export_to_csv, export_to_bed, export_to_json
from utils.visualization import (
    plot_motif_distribution, plot_coverage_map, plot_score_distribution,
    plot_length_distribution, plot_nested_pie_chart, save_all_plots
)

# Try to import Biopython for NCBI downloads
try:
    from Bio import Entrez, SeqIO
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False
    print("Warning: Biopython not available. Cannot download from NCBI.")

# Configuration
NCBI_EMAIL = "yvrajesh_bt@kluniversity.in"
ECOLI_ACCESSION = "U00096.3"  # E. coli K-12 MG1655 complete genome
OUTPUT_DIR = Path("ecoli_analysis_results")
CHUNK_SIZE = 100000  # Process genome in chunks for memory efficiency

def download_ecoli_genome(accession=ECOLI_ACCESSION, output_file=None):
    """
    Download E. coli genome from NCBI.
    
    Args:
        accession: NCBI accession number
        output_file: Path to save the genome FASTA file
    
    Returns:
        Path to downloaded genome file
    """
    if not BIO_AVAILABLE:
        raise ImportError("Biopython is required to download genome from NCBI")
    
    if output_file is None:
        output_file = OUTPUT_DIR / f"{accession}.fasta"
    
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Check if already downloaded
    if output_file.exists():
        print(f"✓ Genome file already exists: {output_file}")
        return output_file
    
    print(f"Downloading E. coli genome {accession} from NCBI...")
    Entrez.email = NCBI_EMAIL
    
    try:
        # Fetch the genome
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        genome_data = handle.read()
        handle.close()
        
        # Save to file
        with open(output_file, 'w') as f:
            f.write(genome_data)
        
        print(f"✓ Downloaded genome to: {output_file}")
        return output_file
        
    except Exception as e:
        print(f"✗ Error downloading genome: {e}")
        raise

def analyze_ecoli_genome(genome_file, chunk_size=CHUNK_SIZE):
    """
    Analyze E. coli genome for Non-B DNA motifs.
    
    Args:
        genome_file: Path to E. coli genome FASTA file
        chunk_size: Size of chunks to process (for large genomes)
    
    Returns:
        Dictionary containing analysis results and metadata
    """
    print(f"\n{'='*70}")
    print(f"E. coli K-12 MG1655 Genome Analysis with NBDScanner")
    print(f"{'='*70}\n")
    
    # Load genome
    print("Loading genome sequence...")
    
    # Read file content
    with open(genome_file, 'r') as f:
        fasta_content = f.read()
    
    genome_data = parse_fasta(fasta_content)
    
    if not genome_data:
        raise ValueError("No sequences found in genome file")
    
    # Get the first (and usually only) sequence
    seq_name = list(genome_data.keys())[0]
    sequence = list(genome_data.values())[0]
    
    seq_length = len(sequence)
    print(f"Sequence: {seq_name}")
    print(f"Length: {seq_length:,} bp")
    
    # Basic statistics
    basic_stats = get_basic_stats(sequence)
    print(f"GC Content: {basic_stats['GC%']:.2f}%")
    print(f"AT Content: {basic_stats['AT%']:.2f}%")
    
    # Start timing
    start_time = time.time()
    
    # Analyze genome
    print(f"\nStarting motif detection...")
    print(f"This may take several minutes for a {seq_length:,} bp genome...")
    
    results = analyze_sequence(sequence, seq_name)
    
    # Calculate elapsed time
    elapsed_time = time.time() - start_time
    processing_speed = seq_length / elapsed_time if elapsed_time > 0 else 0
    
    print(f"\n{'='*70}")
    print(f"Analysis Complete!")
    print(f"{'='*70}")
    print(f"Time elapsed: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
    print(f"Processing speed: {processing_speed:,.0f} bp/second")
    print(f"Total motifs detected: {len(results)}")
    
    # Separate regular motifs from hybrid/cluster
    regular_motifs = [m for m in results if m.get('Class') not in ['Hybrid', 'Non-B_DNA_Clusters']]
    hybrid_cluster_motifs = [m for m in results if m.get('Class') in ['Hybrid', 'Non-B_DNA_Clusters']]
    
    print(f"  - Regular motifs: {len(regular_motifs)}")
    print(f"  - Hybrid/Cluster motifs: {len(hybrid_cluster_motifs)}")
    
    # Motif class breakdown
    print(f"\nMotif Classes Detected:")
    class_counts = Counter(m.get('Class', 'Unknown') for m in regular_motifs)
    for motif_class, count in sorted(class_counts.items(), key=lambda x: -x[1]):
        percentage = (count / len(regular_motifs) * 100) if regular_motifs else 0
        print(f"  {motif_class:30s}: {count:6d} ({percentage:5.2f}%)")
    
    # Subclass breakdown
    print(f"\nTop 10 Motif Subclasses:")
    subclass_counts = Counter(m.get('Subclass', 'Unknown') for m in regular_motifs)
    for i, (subclass, count) in enumerate(subclass_counts.most_common(10), 1):
        percentage = (count / len(regular_motifs) * 100) if regular_motifs else 0
        print(f"  {i:2d}. {subclass:30s}: {count:6d} ({percentage:5.2f}%)")
    
    # Coverage statistics
    total_coverage = sum(m['End'] - m['Start'] for m in regular_motifs)
    coverage_pct = (total_coverage / seq_length * 100) if seq_length > 0 else 0
    density = (len(regular_motifs) / seq_length * 1000) if seq_length > 0 else 0
    
    print(f"\nCoverage Statistics:")
    print(f"  Total coverage: {total_coverage:,} bp ({coverage_pct:.2f}%)")
    print(f"  Motif density: {density:.2f} motifs/kb")
    
    # Prepare results dictionary
    analysis_results = {
        'metadata': {
            'sequence_name': seq_name,
            'sequence_length': seq_length,
            'gc_content': basic_stats['GC%'],
            'at_content': basic_stats['AT%'],
            'analysis_time': elapsed_time,
            'processing_speed': processing_speed,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        },
        'motifs': results,
        'regular_motifs': regular_motifs,
        'hybrid_cluster_motifs': hybrid_cluster_motifs,
        'statistics': {
            'total_motifs': len(results),
            'regular_motifs': len(regular_motifs),
            'hybrid_cluster_motifs': len(hybrid_cluster_motifs),
            'class_counts': dict(class_counts),
            'subclass_counts': dict(subclass_counts),
            'coverage_bp': total_coverage,
            'coverage_percent': coverage_pct,
            'density_per_kb': density
        }
    }
    
    return analysis_results

def generate_visualizations(results, output_dir=OUTPUT_DIR):
    """
    Generate comprehensive visualizations for E. coli analysis.
    
    Args:
        results: Analysis results dictionary
        output_dir: Directory to save plots
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*70}")
    print(f"Generating Visualizations")
    print(f"{'='*70}\n")
    
    motifs = results['regular_motifs']
    seq_length = results['metadata']['sequence_length']
    seq_name = results['metadata']['sequence_name']
    
    if not motifs:
        print("No motifs to visualize")
        return
    
    # Set publication style
    plt.style.use('seaborn-v0_8-darkgrid')
    sns.set_palette("husl")
    
    plots = []
    
    # 1. Motif class distribution
    print("1. Creating motif class distribution plot...")
    fig1 = plot_motif_distribution(motifs, by='Class', title=f"Motif Classes - {seq_name}")
    fig1.savefig(output_dir / "motif_class_distribution.png", dpi=300, bbox_inches='tight')
    plots.append(fig1)
    plt.close(fig1)
    
    # 2. Motif subclass distribution
    print("2. Creating motif subclass distribution plot...")
    fig2 = plot_motif_distribution(motifs, by='Subclass', title=f"Motif Subclasses - {seq_name}")
    fig2.savefig(output_dir / "motif_subclass_distribution.png", dpi=300, bbox_inches='tight')
    plots.append(fig2)
    plt.close(fig2)
    
    # 3. Coverage map
    print("3. Creating coverage map...")
    fig3 = plot_coverage_map(motifs, seq_length, title=f"Motif Coverage - {seq_name}")
    fig3.savefig(output_dir / "motif_coverage_map.png", dpi=300, bbox_inches='tight')
    plots.append(fig3)
    plt.close(fig3)
    
    # 4. Score distribution
    print("4. Creating score distribution plot...")
    fig4 = plot_score_distribution(motifs, by_class=True, title="Score Distribution by Class")
    fig4.savefig(output_dir / "score_distribution.png", dpi=300, bbox_inches='tight')
    plots.append(fig4)
    plt.close(fig4)
    
    # 5. Length distribution
    print("5. Creating length distribution plot...")
    fig5 = plot_length_distribution(motifs, by_class=True, title="Length Distribution by Class")
    fig5.savefig(output_dir / "length_distribution.png", dpi=300, bbox_inches='tight')
    plots.append(fig5)
    plt.close(fig5)
    
    # 6. Nested pie chart
    print("6. Creating nested pie chart...")
    fig6 = plot_nested_pie_chart(motifs, title=f"Class-Subclass Distribution - {seq_name}")
    fig6.savefig(output_dir / "nested_pie_chart.png", dpi=300, bbox_inches='tight')
    plots.append(fig6)
    plt.close(fig6)
    
    print(f"\n✓ All visualizations saved to: {output_dir}")
    
    return plots

def export_results(results, output_dir=OUTPUT_DIR):
    """
    Export analysis results in multiple formats.
    
    Args:
        results: Analysis results dictionary
        output_dir: Directory to save exports
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*70}")
    print(f"Exporting Results")
    print(f"{'='*70}\n")
    
    motifs = results['regular_motifs']
    seq_name = results['metadata']['sequence_name']
    
    # 1. CSV export
    print("1. Exporting to CSV...")
    csv_data = export_to_csv(motifs)
    csv_file = output_dir / "ecoli_motifs.csv"
    with open(csv_file, 'w') as f:
        f.write(csv_data)
    print(f"   ✓ CSV saved to: {csv_file}")
    
    # 2. BED export
    print("2. Exporting to BED...")
    bed_data = export_to_bed(motifs, seq_name)
    bed_file = output_dir / "ecoli_motifs.bed"
    with open(bed_file, 'w') as f:
        f.write(bed_data)
    print(f"   ✓ BED saved to: {bed_file}")
    
    # 3. JSON export
    print("3. Exporting to JSON...")
    json_data = export_to_json(motifs, pretty=True)
    json_file = output_dir / "ecoli_motifs.json"
    with open(json_file, 'w') as f:
        f.write(json_data)
    print(f"   ✓ JSON saved to: {json_file}")
    
    # 4. Summary statistics
    print("4. Exporting summary statistics...")
    summary_file = output_dir / "analysis_summary.json"
    with open(summary_file, 'w') as f:
        json.dump({
            'metadata': results['metadata'],
            'statistics': results['statistics']
        }, f, indent=2)
    print(f"   ✓ Summary saved to: {summary_file}")
    
    # 5. Detailed report
    print("5. Generating detailed report...")
    report_file = output_dir / "analysis_report.txt"
    with open(report_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("E. coli K-12 MG1655 Non-B DNA Motif Analysis Report\n")
        f.write("="*70 + "\n\n")
        
        f.write("METADATA\n")
        f.write("-"*70 + "\n")
        for key, value in results['metadata'].items():
            f.write(f"{key:25s}: {value}\n")
        
        f.write("\n" + "="*70 + "\n")
        f.write("STATISTICS\n")
        f.write("-"*70 + "\n")
        for key, value in results['statistics'].items():
            if isinstance(value, dict):
                f.write(f"\n{key}:\n")
                for k, v in sorted(value.items(), key=lambda x: -x[1] if isinstance(x[1], (int, float)) else 0):
                    f.write(f"  {k:30s}: {v}\n")
            else:
                f.write(f"{key:25s}: {value}\n")
        
        f.write("\n" + "="*70 + "\n")
    
    print(f"   ✓ Report saved to: {report_file}")
    print(f"\n✓ All exports saved to: {output_dir}")

def main():
    """Main execution function."""
    print("""
    ╔═══════════════════════════════════════════════════════════════════╗
    ║  E. coli K-12 MG1655 Genome Analysis with NBDScanner             ║
    ║  Non-B DNA Motif Detection and Characterization                  ║
    ╚═══════════════════════════════════════════════════════════════════╝
    """)
    
    try:
        # Step 1: Download genome or use local file
        print("\nStep 1: Downloading E. coli genome from NCBI...")
        
        # Check if sample file exists
        sample_file = Path("ecoli_sample.fasta")
        if sample_file.exists():
            print(f"✓ Using local sample file: {sample_file}")
            genome_file = sample_file
        else:
            try:
                genome_file = download_ecoli_genome()
            except Exception as e:
                print(f"✗ Could not download from NCBI: {e}")
                print("ℹ️  Using sample E. coli sequence instead...")
                genome_file = sample_file
                if not genome_file.exists():
                    raise ValueError("No genome file available")
        
        # Step 2: Analyze genome
        print("\nStep 2: Analyzing genome for Non-B DNA motifs...")
        results = analyze_ecoli_genome(genome_file)
        
        # Step 3: Generate visualizations
        print("\nStep 3: Generating visualizations...")
        generate_visualizations(results)
        
        # Step 4: Export results
        print("\nStep 4: Exporting results...")
        export_results(results)
        
        print(f"\n{'='*70}")
        print(f"Analysis Complete!")
        print(f"{'='*70}")
        print(f"\nResults saved to: {OUTPUT_DIR}")
        print(f"\nFiles generated:")
        print(f"  - Genome: {genome_file}")
        print(f"  - Visualizations: {OUTPUT_DIR}/*.png")
        print(f"  - Data exports: {OUTPUT_DIR}/*.csv, *.bed, *.json")
        print(f"  - Report: {OUTPUT_DIR}/analysis_report.txt")
        print(f"\n✓ All done!")
        
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
