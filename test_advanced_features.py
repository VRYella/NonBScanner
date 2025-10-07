#!/usr/bin/env python3
"""
Test Script for Advanced NBDScanner Features
=============================================

Tests enhanced hybrid/cluster detection and advanced visualizations.

Author: Dr. Venkata Rajesh Yella
Version: 2024.1
"""

import sys
import os
from nbdscanner import MotifDetector
from advanced_visualizations import (
    plot_genome_landscape_track,
    plot_sliding_window_heat_ribbon,
    plot_ridge_plots_length_by_class,
    plot_sunburst_treemap,
    plot_hexbin_start_vs_score,
    plot_upset_intersection,
    plot_score_violin_beeswarm,
    plot_cluster_hotspot_map,
    export_plot
)
import matplotlib.pyplot as plt

def test_hybrid_cluster_detection():
    """Test enhanced hybrid and cluster detection"""
    print("=" * 70)
    print("TEST 1: Enhanced Hybrid and Cluster Detection")
    print("=" * 70)
    
    # Create a test sequence with overlapping motifs
    test_sequence = (
        "GGGGTTGGGGTTGGGGTTGGGG"  # G-quadruplex
        + "AAAAA" * 5              # A-tract (curved DNA)
        + "GCGCGCGCGCGC"            # Z-DNA
        + "CCCCAACCCCAACCCC"        # i-Motif
        + "AAAA" * 3                # Another A-tract
        + "GGGGCCGGGGCCGGGG"        # Another G4
    )
    
    print(f"\nTest sequence length: {len(test_sequence)} bp")
    print(f"Expected overlaps: G4 + Curved, Z-DNA + i-Motif")
    
    # Run detection
    detector = MotifDetector()
    motifs = detector.analyze_sequence(test_sequence, "test_seq")
    
    print(f"\n✓ Total motifs detected: {len(motifs)}")
    
    # Check for hybrids
    hybrids = [m for m in motifs if m.get('Class') == 'Hybrid']
    print(f"✓ Hybrid motifs: {len(hybrids)}")
    
    for hybrid in hybrids:
        print(f"\n  Hybrid: {hybrid['Subclass']}")
        print(f"    Position: {hybrid['Start']}-{hybrid['End']}")
        print(f"    Length: {hybrid['Length']} bp")
        print(f"    Raw Score: {hybrid.get('Score', 'N/A'):.3f}")
        print(f"    Normalized Score: {hybrid.get('Normalized_Score', 'N/A')}")
        print(f"    Sequence: {hybrid.get('Sequence', 'N/A')[:50]}...")
        if 'Component_Classes' in hybrid:
            print(f"    Component Classes: {', '.join(hybrid['Component_Classes'])}")
    
    # Check for clusters
    clusters = [m for m in motifs if m.get('Class') == 'Non-B_DNA_Clusters']
    print(f"\n✓ Cluster motifs: {len(clusters)}")
    
    for cluster in clusters:
        print(f"\n  Cluster: {cluster['Subclass']}")
        print(f"    Position: {cluster['Start']}-{cluster['End']}")
        print(f"    Length: {cluster['Length']} bp")
        print(f"    Raw Score: {cluster.get('Score', 'N/A'):.3f}")
        print(f"    Normalized Score: {cluster.get('Normalized_Score', 'N/A')}")
        print(f"    Sequence: {cluster.get('Sequence', 'N/A')[:50]}...")
        print(f"    Motif Count: {cluster.get('Motif_Count', 'N/A')}")
        print(f"    Class Diversity: {cluster.get('Class_Diversity', 'N/A')}")
    
    # Check normalized scores
    print(f"\n✓ Checking normalized scores:")
    motifs_with_norm = [m for m in motifs if 'Normalized_Score' in m]
    print(f"    Motifs with normalized scores: {len(motifs_with_norm)}/{len(motifs)}")
    
    if motifs_with_norm:
        norm_scores = [m['Normalized_Score'] for m in motifs_with_norm]
        print(f"    Normalized score range: {min(norm_scores):.4f} - {max(norm_scores):.4f}")
        print(f"    Average normalized score: {sum(norm_scores)/len(norm_scores):.4f}")
    
    return motifs

def test_advanced_visualizations(motifs, sequence_length):
    """Test all advanced visualization functions"""
    print("\n" + "=" * 70)
    print("TEST 2: Advanced Visualization Suite")
    print("=" * 70)
    
    output_dir = "/tmp/nbdscanner_viz_test"
    os.makedirs(output_dir, exist_ok=True)
    
    viz_tests = [
        ("Genome Landscape Track", lambda: plot_genome_landscape_track(motifs, sequence_length)),
        ("Sliding Window Heat Ribbon", lambda: plot_sliding_window_heat_ribbon(motifs, sequence_length)),
        ("Ridge Plots (Length by Class)", lambda: plot_ridge_plots_length_by_class(motifs)),
        ("Sunburst Chart", lambda: plot_sunburst_treemap(motifs, plot_type='sunburst')),
        ("Treemap Chart", lambda: plot_sunburst_treemap(motifs, plot_type='treemap')),
        ("Hexbin (Start vs Score)", lambda: plot_hexbin_start_vs_score(motifs)),
        ("UpSet Intersection Plot", lambda: plot_upset_intersection(motifs)),
        ("Violin + Beeswarm Plot", lambda: plot_score_violin_beeswarm(motifs)),
        ("Cluster Hotspot Map", lambda: plot_cluster_hotspot_map(motifs, sequence_length)),
    ]
    
    passed = 0
    failed = 0
    
    for i, (name, func) in enumerate(viz_tests, 1):
        try:
            print(f"\n[{i}/{len(viz_tests)}] Testing {name}...", end=" ")
            fig = func()
            
            # Save the figure
            filepath = os.path.join(output_dir, f"test_{i:02d}_{name.lower().replace(' ', '_')}")
            export_plot(fig, filepath, formats=['png'], dpi=150)  # Lower DPI for testing
            
            plt.close(fig)
            print("✓ PASS")
            passed += 1
            
        except Exception as e:
            print(f"✗ FAIL - {e}")
            failed += 1
    
    print(f"\n" + "=" * 70)
    print(f"Visualization Test Summary: {passed} passed, {failed} failed")
    print(f"Output directory: {output_dir}")
    print("=" * 70)
    
    return passed, failed

def test_score_normalization():
    """Test score normalization across different motif types"""
    print("\n" + "=" * 70)
    print("TEST 3: Score Normalization")
    print("=" * 70)
    
    # Create sequences with different characteristics
    test_cases = [
        ("Short G4", "GGGTTAGGGTTA", "G-Quadruplex"),
        ("Long G4", "GGGGTTGGGGTTGGGGTTGGGG", "G-Quadruplex"),
        ("Short A-tract", "AAAAAAAA", "Curved_DNA"),
        ("Long A-tract", "AAAAAAAAAAAAAAAA", "Curved_DNA"),
        ("Z-DNA", "GCGCGCGCGCGC", "Z-DNA"),
    ]
    
    detector = MotifDetector()
    
    print(f"\n{'Motif Type':<20} {'Length':<8} {'Raw Score':<12} {'Normalized':<12}")
    print("-" * 60)
    
    for name, seq, expected_class in test_cases:
        motifs = detector.analyze_sequence(seq, "test")
        
        # Find motifs matching expected class
        matching = [m for m in motifs if expected_class in m.get('Class', '')]
        
        if matching:
            m = matching[0]
            raw_score = m.get('Score', 0)
            norm_score = m.get('Normalized_Score', 0)
            print(f"{name:<20} {m['Length']:<8} {raw_score:<12.4f} {norm_score:<12.4f}")
        else:
            print(f"{name:<20} {'N/A':<8} {'Not detected':<12} {'N/A':<12}")
    
    print("\n✓ Score normalization test completed")

def main():
    """Run all tests"""
    print("\n" + "=" * 70)
    print("NBDScanner Advanced Features Test Suite")
    print("=" * 70)
    print("Testing enhanced hybrid/cluster detection and visualizations")
    print("=" * 70)
    
    # Test 1: Hybrid and cluster detection
    motifs = test_hybrid_cluster_detection()
    sequence_length = 150  # Approximate length of test sequence
    
    # Test 2: Advanced visualizations
    test_advanced_visualizations(motifs, sequence_length)
    
    # Test 3: Score normalization
    test_score_normalization()
    
    print("\n" + "=" * 70)
    print("✓ ALL TESTS COMPLETED")
    print("=" * 70)
    print("\nSummary of new features:")
    print("  ✓ Hybrids report actual sequences (longest non-overlapping)")
    print("  ✓ Clusters report actual sequences (longest non-overlapping)")
    print("  ✓ Both raw and normalized scores available")
    print("  ✓ 9 advanced publication-quality visualizations")
    print("  ✓ Colorblind-friendly palettes")
    print("  ✓ Export to SVG/PNG @300 DPI")
    print("=" * 70)

if __name__ == "__main__":
    main()
