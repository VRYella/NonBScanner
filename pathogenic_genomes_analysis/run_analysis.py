#!/usr/bin/env python3
"""
Comprehensive NonBScanner analysis of pathogenic genomes.
This script performs:
1. Non-B DNA motif detection on all genomes
2. Statistical analysis and comparison
3. Publication-quality visualizations
4. Excel export with detailed results
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Add parent directory to path to import nonbscanner
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import nonbscanner as nbs
from utilities import (
    read_fasta_file,
    export_to_excel,
    analyze_class_subclass_detection,
    print_detection_report,
    calculate_motif_statistics
)
from visualizations import (
    plot_class_analysis_comprehensive,
    plot_subclass_analysis_comprehensive,
    plot_score_statistics_by_class,
    plot_length_statistics_by_class
)

# Set matplotlib style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

# Directories
DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
FIGURES_DIR = Path("figures")
TABLES_DIR = Path("tables")

# Create directories
for dir_path in [RESULTS_DIR, FIGURES_DIR, TABLES_DIR]:
    dir_path.mkdir(exist_ok=True)

def analyze_all_genomes():
    """Analyze all pathogenic genomes with NonBScanner."""
    print("=" * 80)
    print("NONBSCANNER ANALYSIS OF PATHOGENIC GENOMES")
    print("=" * 80)
    print()
    
    # Find all FASTA files
    fasta_files = sorted(DATA_DIR.glob("*.fasta"))
    
    if not fasta_files:
        print("ERROR: No FASTA files found in data directory!")
        return None
    
    print(f"Found {len(fasta_files)} genome(s) to analyze:")
    for f in fasta_files:
        print(f"  - {f.name}")
    print()
    
    # Store all results
    all_results = {}
    all_motifs_combined = []
    
    # Analyze each genome
    for fasta_file in fasta_files:
        genome_name = fasta_file.stem
        print(f"\n{'='*80}")
        print(f"Analyzing: {genome_name}")
        print(f"{'='*80}")
        
        # Read sequences
        sequences = read_fasta_file(str(fasta_file))
        
        if not sequences:
            print(f"  WARNING: No sequences found in {fasta_file.name}")
            continue
        
        # Analyze each sequence in the file
        genome_motifs = []
        for seq_name, sequence in sequences.items():
            print(f"\nSequence: {seq_name}")
            print(f"Length: {len(sequence):,} bp")
            
            # Run NonBScanner (with fast mode if available)
            try:
                motifs = nbs.analyze_sequence(sequence, seq_name, use_fast_mode=True)
            except:
                # Fallback to standard mode
                motifs = nbs.analyze_sequence(sequence, seq_name, use_fast_mode=False)
            
            print(f"Detected: {len(motifs)} Non-B DNA motifs")
            
            # Add genome identifier
            for motif in motifs:
                motif['Genome'] = genome_name
            
            genome_motifs.extend(motifs)
        
        all_results[genome_name] = genome_motifs
        all_motifs_combined.extend(genome_motifs)
        
        # Export individual genome results
        if genome_motifs:
            excel_file = RESULTS_DIR / f"{genome_name}_results.xlsx"
            export_to_excel(genome_motifs, str(excel_file))
            print(f"\n✓ Exported results to {excel_file}")
    
    print(f"\n{'='*80}")
    print(f"ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"Total genomes analyzed: {len(all_results)}")
    print(f"Total motifs detected: {len(all_motifs_combined)}")
    print()
    
    return all_results, all_motifs_combined

def generate_comparative_analysis(all_results, all_motifs_combined):
    """Generate comparative analysis across genomes."""
    print("\n" + "="*80)
    print("GENERATING COMPARATIVE ANALYSIS")
    print("="*80)
    
    # Create comprehensive Excel export
    print("\nExporting combined results...")
    combined_excel = RESULTS_DIR / "all_genomes_combined.xlsx"
    export_to_excel(all_motifs_combined, str(combined_excel))
    print(f"✓ Saved to {combined_excel}")
    
    # Summary statistics by genome
    print("\nGenerating summary statistics...")
    summary_data = []
    
    for genome_name, motifs in all_results.items():
        if not motifs:
            continue
        
        # Get sequence length (from first motif)
        seq_length = max([m.get('Sequence_Length', 0) for m in motifs])
        
        # Count by class
        class_counts = {}
        for motif in motifs:
            cls = motif.get('Class', 'Unknown')
            class_counts[cls] = class_counts.get(cls, 0) + 1
        
        summary_data.append({
            'Genome': genome_name,
            'Length (bp)': seq_length,
            'Total Motifs': len(motifs),
            'Motifs per kb': (len(motifs) / seq_length * 1000) if seq_length > 0 else 0,
            'Unique Classes': len(class_counts),
            **class_counts
        })
    
    # Create summary DataFrame
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.fillna(0)
    
    # Save summary table
    summary_csv = TABLES_DIR / "genome_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"✓ Summary saved to {summary_csv}")
    
    # Print summary
    print("\n" + "="*80)
    print("SUMMARY BY GENOME")
    print("="*80)
    print(summary_df.to_string(index=False))
    print()
    
    return summary_df

def generate_visualizations(all_motifs_combined, summary_df):
    """Generate publication-quality visualizations."""
    print("\n" + "="*80)
    print("GENERATING PUBLICATION-QUALITY VISUALIZATIONS")
    print("="*80)
    print()
    
    if not all_motifs_combined:
        print("No motifs to visualize!")
        return
    
    # 1. Comprehensive class analysis
    print("1. Comprehensive class analysis...")
    fig = plot_class_analysis_comprehensive(all_motifs_combined)
    fig_path = FIGURES_DIR / "fig1_class_analysis.png"
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"   ✓ Saved to {fig_path}")
    
    # 2. Comprehensive subclass analysis
    print("2. Comprehensive subclass analysis...")
    fig = plot_subclass_analysis_comprehensive(all_motifs_combined)
    fig_path = FIGURES_DIR / "fig2_subclass_analysis.png"
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"   ✓ Saved to {fig_path}")
    
    # 3. Score statistics by class
    print("3. Score statistics by class...")
    fig = plot_score_statistics_by_class(all_motifs_combined)
    fig_path = FIGURES_DIR / "fig3_score_statistics.png"
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"   ✓ Saved to {fig_path}")
    
    # 4. Length statistics by class
    print("4. Length statistics by class...")
    fig = plot_length_statistics_by_class(all_motifs_combined)
    fig_path = FIGURES_DIR / "fig4_length_statistics.png"
    fig.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"   ✓ Saved to {fig_path}")
    
    # 5. Genome comparison bar chart
    print("5. Genome comparison chart...")
    create_genome_comparison_chart(summary_df)
    
    # 6. Motif density heatmap
    print("6. Motif density heatmap...")
    create_motif_density_heatmap(all_motifs_combined)
    
    print("\n✓ All visualizations generated!")

def create_genome_comparison_chart(summary_df):
    """Create comparison chart across genomes."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Total motifs per genome
    ax1 = axes[0]
    ax1.bar(summary_df['Genome'], summary_df['Total Motifs'], 
            color=sns.color_palette("husl", len(summary_df)))
    ax1.set_xlabel('Genome', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Total Non-B DNA Motifs', fontsize=12, fontweight='bold')
    ax1.set_title('Total Non-B DNA Motifs by Genome', fontsize=14, fontweight='bold')
    ax1.tick_params(axis='x', rotation=45)
    ax1.grid(axis='y', alpha=0.3)
    
    # Motifs per kb (normalized by genome size)
    ax2 = axes[1]
    ax2.bar(summary_df['Genome'], summary_df['Motifs per kb'],
            color=sns.color_palette("husl", len(summary_df)))
    ax2.set_xlabel('Genome', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Motifs per kb', fontsize=12, fontweight='bold')
    ax2.set_title('Non-B DNA Density (Motifs/kb)', fontsize=14, fontweight='bold')
    ax2.tick_params(axis='x', rotation=45)
    ax2.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    fig_path = FIGURES_DIR / "fig5_genome_comparison.png"
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   ✓ Saved to {fig_path}")

def create_motif_density_heatmap(all_motifs_combined):
    """Create heatmap of motif class distribution across genomes."""
    # Create pivot table
    df = pd.DataFrame(all_motifs_combined)
    
    if df.empty or 'Genome' not in df.columns or 'Class' not in df.columns:
        print("   ⚠ Insufficient data for heatmap")
        return
    
    pivot = pd.crosstab(df['Genome'], df['Class'])
    
    # Create heatmap
    plt.figure(figsize=(14, 8))
    sns.heatmap(pivot, annot=True, fmt='d', cmap='YlOrRd', 
                cbar_kws={'label': 'Number of Motifs'},
                linewidths=0.5, linecolor='gray')
    plt.xlabel('Non-B DNA Class', fontsize=12, fontweight='bold')
    plt.ylabel('Genome', fontsize=12, fontweight='bold')
    plt.title('Non-B DNA Motif Distribution Across Pathogenic Genomes', 
              fontsize=14, fontweight='bold', pad=20)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    fig_path = FIGURES_DIR / "fig6_motif_heatmap.png"
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   ✓ Saved to {fig_path}")

def main():
    """Main analysis pipeline."""
    print("\n" + "▓"*80)
    print("▓" + " "*78 + "▓")
    print("▓" + "  NONBSCANNER: PATHOGENIC GENOME COMPARATIVE ANALYSIS  ".center(78) + "▓")
    print("▓" + " "*78 + "▓")
    print("▓"*80 + "\n")
    
    # Step 1: Analyze all genomes
    results = analyze_all_genomes()
    if results is None:
        return
    
    all_results, all_motifs_combined = results
    
    # Step 2: Comparative analysis
    summary_df = generate_comparative_analysis(all_results, all_motifs_combined)
    
    # Step 3: Generate visualizations
    generate_visualizations(all_motifs_combined, summary_df)
    
    # Step 4: Print class/subclass detection report
    print("\n" + "="*80)
    print("NON-B DNA CLASS/SUBCLASS DETECTION REPORT")
    print("="*80)
    detection_report = analyze_class_subclass_detection(all_motifs_combined)
    report_text = print_detection_report(detection_report)
    print(report_text)
    
    # Save report
    report_file = TABLES_DIR / "detection_report.txt"
    with open(report_file, 'w') as f:
        f.write(report_text)
    print(f"\n✓ Detection report saved to {report_file}")
    
    print("\n" + "▓"*80)
    print("▓" + " "*78 + "▓")
    print("▓" + "  ANALYSIS COMPLETE - CHECK results/, figures/, tables/  ".center(78) + "▓")
    print("▓" + " "*78 + "▓")
    print("▓"*80 + "\n")

if __name__ == "__main__":
    main()
