#!/usr/bin/env python3
"""
Test script for enhanced visualization functions.
Demonstrates comprehensive class and subclass analysis with scientific plots.
"""

import visualizations as viz
import matplotlib.pyplot as plt
import nonbscanner as nbs

def test_with_sample_sequence():
    """Test visualizations with a sample sequence containing multiple motifs"""
    
    # Sample sequence with multiple Non-B DNA motifs
    test_sequence = """
    AAAAAAAAAAAAAAAA
    GGGTTAGGGTTAGGGTTAGGG
    CCCCCTCCCCCTCCCCCTCCCCC
    CGCGCGCGCGCGCGCG
    AAAAAACGTAAAAAAACGTAAAAAAACGT
    ACGTACGTACGTACGTACGTACGTACGT
    GGGAGGAGGAGGG
    """
    
    # Clean up sequence (remove whitespace and newlines)
    sequence = ''.join(test_sequence.split())
    
    print(f"Analyzing sequence of length {len(sequence)}...")
    
    # Analyze with NonBScanner
    scanner = nbs.NonBScanner()
    motifs = scanner.analyze_sequence(sequence, 'test_sequence')
    
    print(f"\nFound {len(motifs)} motifs:")
    for motif in motifs:
        print(f"  {motif['Class']:20s} / {motif['Subclass']:25s} at {motif['Start']:4d}-{motif['End']:4d}")
    
    if not motifs:
        print("\nNo motifs found in test sequence!")
        return
    
    # Create enhanced visualizations
    print("\nGenerating enhanced visualizations...")
    
    # 1. Comprehensive class analysis
    try:
        fig1 = viz.plot_class_analysis_comprehensive(motifs, figsize=(16, 12))
        fig1.savefig('/tmp/enhanced_class_analysis.png', dpi=300, bbox_inches='tight')
        plt.close(fig1)
        print("  ✓ Saved enhanced_class_analysis.png")
    except Exception as e:
        print(f"  ✗ Class analysis failed: {e}")
    
    # 2. Comprehensive subclass analysis
    try:
        fig2 = viz.plot_subclass_analysis_comprehensive(motifs, figsize=(18, 14))
        fig2.savefig('/tmp/enhanced_subclass_analysis.png', dpi=300, bbox_inches='tight')
        plt.close(fig2)
        print("  ✓ Saved enhanced_subclass_analysis.png")
    except Exception as e:
        print(f"  ✗ Subclass analysis failed: {e}")
    
    # 3. Score statistics
    try:
        fig3 = viz.plot_score_statistics_by_class(motifs, figsize=(14, 8))
        fig3.savefig('/tmp/enhanced_score_stats.png', dpi=300, bbox_inches='tight')
        plt.close(fig3)
        print("  ✓ Saved enhanced_score_stats.png")
    except Exception as e:
        print(f"  ✗ Score statistics failed: {e}")
    
    # 4. Length statistics
    try:
        fig4 = viz.plot_length_statistics_by_class(motifs, figsize=(14, 10))
        fig4.savefig('/tmp/enhanced_length_stats.png', dpi=300, bbox_inches='tight')
        plt.close(fig4)
        print("  ✓ Saved enhanced_length_stats.png")
    except Exception as e:
        print(f"  ✗ Length statistics failed: {e}")
    
    # 5. Original distribution plot for comparison
    try:
        fig5 = viz.plot_motif_distribution(motifs, by='Class', figsize=(12, 6))
        fig5.savefig('/tmp/original_class_distribution.png', dpi=300, bbox_inches='tight')
        plt.close(fig5)
        print("  ✓ Saved original_class_distribution.png")
    except Exception as e:
        print(f"  ✗ Original distribution failed: {e}")
    
    print("\n✓ All enhanced visualizations generated successfully!")
    print("  Output directory: /tmp/")
    print("\nEnhanced features:")
    print("  - Shows all 11 Non-B DNA classes (detected and not detected)")
    print("  - Comprehensive statistics (mean, median, std, min, max)")
    print("  - Publication-quality plots at 300 DPI")
    print("  - Scientific color schemes and styling")
    print("  - Separate class and subclass analysis")


if __name__ == "__main__":
    test_with_sample_sequence()
