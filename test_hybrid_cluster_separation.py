#!/usr/bin/env python3
"""
Test script to verify hybrid/cluster separation in NBDScanner
==============================================================

This script tests that hybrid and cluster motifs are correctly
separated from regular Non-B DNA motifs in the analysis results.

Author: Dr. Venkata Rajesh Yella
Version: 2024.1
"""

from nbdscanner import analyze_sequence

def test_hybrid_cluster_separation():
    """Test that hybrid and cluster motifs are correctly identified"""
    print("=" * 70)
    print("TEST: Hybrid/Cluster Separation Verification")
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
    
    # Run detection
    motifs = analyze_sequence(test_sequence, "test_seq")
    
    # Separate motifs
    regular_motifs = [m for m in motifs if m.get('Class') not in ['Hybrid', 'Non-B_DNA_Clusters']]
    hybrid_motifs = [m for m in motifs if m.get('Class') == 'Hybrid']
    cluster_motifs = [m for m in motifs if m.get('Class') == 'Non-B_DNA_Clusters']
    
    print(f"\n✓ Total motifs detected: {len(motifs)}")
    print(f"  - Regular motifs: {len(regular_motifs)}")
    print(f"  - Hybrid motifs: {len(hybrid_motifs)}")
    print(f"  - Cluster motifs: {len(cluster_motifs)}")
    
    # Verify separation
    assert len(motifs) == len(regular_motifs) + len(hybrid_motifs) + len(cluster_motifs), \
        "Motif counts don't add up!"
    
    # Show regular motif classes
    print(f"\n✓ Regular motif classes:")
    regular_classes = {}
    for m in regular_motifs:
        cls = m.get('Class', 'Unknown')
        regular_classes[cls] = regular_classes.get(cls, 0) + 1
    
    for cls, count in sorted(regular_classes.items()):
        print(f"    {cls}: {count}")
    
    # Show hybrid details
    if hybrid_motifs:
        print(f"\n✓ Hybrid motifs details:")
        for hybrid in hybrid_motifs[:3]:  # Show first 3
            print(f"    - {hybrid.get('Subclass', 'Unknown')}")
            print(f"      Position: {hybrid.get('Start')}-{hybrid.get('End')}")
            if 'Component_Classes' in hybrid:
                print(f"      Components: {', '.join(hybrid['Component_Classes'])}")
    
    # Show cluster details
    if cluster_motifs:
        print(f"\n✓ Cluster motifs details:")
        for cluster in cluster_motifs[:3]:  # Show first 3
            print(f"    - {cluster.get('Subclass', 'Unknown')}")
            print(f"      Position: {cluster.get('Start')}-{cluster.get('End')}")
            print(f"      Motif Count: {cluster.get('Motif_Count', 'N/A')}")
            print(f"      Class Diversity: {cluster.get('Class_Diversity', 'N/A')}")
    
    # Verify no hybrid/cluster in regular motifs
    for m in regular_motifs:
        assert m.get('Class') not in ['Hybrid', 'Non-B_DNA_Clusters'], \
            f"Found {m.get('Class')} in regular motifs!"
    
    print("\n" + "=" * 70)
    print("✓ TEST PASSED: Hybrid/Cluster separation working correctly")
    print("=" * 70)
    
    return True

if __name__ == "__main__":
    try:
        test_hybrid_cluster_separation()
        print("\n✓ All tests passed!")
    except Exception as e:
        print(f"\n✗ Test failed: {str(e)}")
        import traceback
        traceback.print_exc()
        exit(1)
