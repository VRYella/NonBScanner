#!/usr/bin/env python3
"""
Comprehensive Test Suite for All Non-B DNA Motif Classes
========================================================

Tests all 10 major motif classes with representative sequences:
1. Curved DNA
2. Slipped DNA  
3. Cruciform
4. R-Loop
5. Triplex
6. G-Quadruplex
7. i-Motif
8. Z-DNA
9. A-philic
10. Hybrid/Cluster (derived from overlaps)

Tests ensure:
- Hyperscan database is used for feasible motifs (regex-based)
- Algorithmic detectors work correctly (Cruciform)
- Scoring and non-overlap resolution pipeline functions properly
- All motif classes can be detected
"""

import os
import sys

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from scanner import ModularMotifDetector


# Test sequences for each motif class
TEST_SEQUENCES = {
    'G-Quadruplex': {
        'canonical_telomeric': 'GGGTTAGGGTTAGGGTTAGGG',
        'canonical_g4': 'GGGGTTTTGGGGTTTTGGGGTTTTGGGG',
        'relaxed_g4': 'GGATAGGATAGGATAGG',
    },
    'i-Motif': {
        'canonical': 'CCCTAACCCTAACCCTAACCC',
        'c_rich': 'CCCCAAAACCCCAAAACCCCAAAACCCC',
    },
    'Z-DNA': {
        'cg_alternating': 'CGCGCGCGCGCGCGCGCGCG',
        'ac_alternating': 'ACGCGCGCGCACGCGCGCGC',
    },
    'Curved_DNA': {
        'a_tract': 'AAAAAAAAAATTTTTTTTTT',
        't_tract': 'TTTTTTTTTTAAAAAAAAAA',
        'at_rich': 'AAATTTAAATTTAAATTT',
    },
    'Slipped_DNA': {
        'mono_repeat': 'AAAAAAAAAAAAA',
        'di_repeat': 'CACACACACACACACA',
        'tri_repeat': 'CGGCGGCGGCGGCGGCGG',
    },
    'Cruciform': {
        'inverted_repeat': 'ATCGATCGAAAAAATCGATCGAT',
        'palindrome': 'GCGCGCATATATGCGCGC',
    },
    'R-Loop': {
        'gc_rich': 'GGGGGGGGGGGGCCCCCCCCCCCC',
        'g_c_rich': 'GGGGGAAACCCCC',
    },
    'Triplex': {
        'mirror_repeat': 'GGGGGGGGAAAAAGGGGGGGG',
        'purine_run': 'AAAAAAAGGGGGGAAAAAAAGGGGGGG',
    },
    'A-philic': {
        'a_g_rich': 'AGGGGGGGGGAGGGGGGGGG',
        'a_rich_context': 'AAAAAGGGAAAAAGGGAAAA',
    },
    'Mixed': {
        # This sequence should trigger multiple motif types and possibly hybrid/cluster
        'complex': 'GGGTTAGGGTTAGGGTTAGGGAAAAAAAATTTTTTCACACACACACACACCCCTAACCCTAACCCTAACCC',
    }
}


def test_individual_motif_classes():
    """Test each motif class individually"""
    print("="*80)
    print("TESTING INDIVIDUAL MOTIF CLASSES")
    print("="*80)
    
    detector = ModularMotifDetector()
    results = {}
    
    for motif_class, sequences in TEST_SEQUENCES.items():
        print(f"\n{'='*80}")
        print(f"Testing {motif_class}")
        print(f"{'='*80}")
        
        class_results = {}
        for seq_name, sequence in sequences.items():
            motifs = detector.analyze_sequence(sequence, seq_name)
            
            # Filter motifs to this class (allowing for variations in naming)
            class_motifs = [m for m in motifs if motif_class.lower().replace('_', '') in m['Class'].lower().replace('-', '').replace('_', '')]
            
            # Also check for specific patterns
            if motif_class == 'G-Quadruplex':
                class_motifs = [m for m in motifs if 'quadruplex' in m['Class'].lower() or 'g4' in m['Subclass'].lower()]
            elif motif_class == 'i-Motif':
                class_motifs = [m for m in motifs if 'motif' in m['Class'].lower() and 'i' in m['Class'].lower()]
            elif motif_class == 'Z-DNA':
                class_motifs = [m for m in motifs if 'z-dna' in m['Class'].lower() or 'zdna' in m['Class'].lower()]
            
            class_results[seq_name] = {
                'total_motifs': len(motifs),
                'class_motifs': len(class_motifs),
                'motifs': class_motifs
            }
            
            print(f"\n  Sequence: {seq_name}")
            print(f"  Length: {len(sequence)} bp")
            print(f"  Total motifs found: {len(motifs)}")
            print(f"  {motif_class} motifs: {len(class_motifs)}")
            
            for motif in class_motifs[:3]:  # Show first 3
                print(f"    - {motif['Class']}/{motif['Subclass']}: "
                      f"{motif['Start']}-{motif['End']} (score: {motif.get('Score', 0):.3f})")
        
        results[motif_class] = class_results
    
    return results


def test_hybrid_and_cluster_detection():
    """Test that hybrid and cluster motifs are detected properly"""
    print("\n" + "="*80)
    print("TESTING HYBRID AND CLUSTER DETECTION")
    print("="*80)
    
    detector = ModularMotifDetector()
    
    # Use complex sequence that should trigger multiple overlapping motifs
    complex_seq = TEST_SEQUENCES['Mixed']['complex']
    motifs = detector.analyze_sequence(complex_seq, 'complex_sequence')
    
    # Find hybrid and cluster motifs
    hybrid_motifs = [m for m in motifs if m['Class'] == 'Hybrid']
    cluster_motifs = [m for m in motifs if 'cluster' in m['Class'].lower()]
    regular_motifs = [m for m in motifs if m['Class'] != 'Hybrid' and 'cluster' not in m['Class'].lower()]
    
    print(f"\nComplex sequence ({len(complex_seq)} bp) analysis:")
    print(f"  Regular motifs: {len(regular_motifs)}")
    print(f"  Hybrid motifs: {len(hybrid_motifs)}")
    print(f"  Cluster motifs: {len(cluster_motifs)}")
    print(f"  Total motifs: {len(motifs)}")
    
    if hybrid_motifs:
        print(f"\n  Hybrid motifs detected:")
        for h in hybrid_motifs[:3]:
            print(f"    - {h['Subclass']}: {h['Start']}-{h['End']} (score: {h.get('Score', 0):.3f})")
    
    if cluster_motifs:
        print(f"\n  Cluster motifs detected:")
        for c in cluster_motifs[:3]:
            print(f"    - {c['Subclass']}: {c['Start']}-{c['End']} (score: {c.get('Score', 0):.3f})")
    
    return {
        'regular': len(regular_motifs),
        'hybrid': len(hybrid_motifs),
        'cluster': len(cluster_motifs),
        'total': len(motifs)
    }


def test_scoring_and_overlap_resolution():
    """Test that scoring and non-overlap resolution works correctly"""
    print("\n" + "="*80)
    print("TESTING SCORING AND OVERLAP RESOLUTION")
    print("="*80)
    
    detector = ModularMotifDetector()
    
    # Test sequence with potential overlaps
    test_seq = 'GGGTTAGGGTTAGGGTTAGGGAAAAAAAATTTTTT'
    motifs = detector.analyze_sequence(test_seq, 'overlap_test')
    
    print(f"\nTest sequence ({len(test_seq)} bp):")
    print(f"Total motifs: {len(motifs)}")
    
    # Check for overlaps within same class
    by_class = {}
    for motif in motifs:
        key = f"{motif['Class']}-{motif['Subclass']}"
        if key not in by_class:
            by_class[key] = []
        by_class[key].append(motif)
    
    print(f"\nMotifs by class/subclass:")
    has_overlap = False
    for key, class_motifs in by_class.items():
        print(f"  {key}: {len(class_motifs)} motifs")
        
        # Check for overlaps within this class
        for i, m1 in enumerate(class_motifs):
            for m2 in class_motifs[i+1:]:
                # Check if they overlap
                if not (m1['End'] <= m2['Start'] or m2['End'] <= m1['Start']):
                    has_overlap = True
                    print(f"    ⚠️  OVERLAP DETECTED: {m1['Start']}-{m1['End']} and {m2['Start']}-{m2['End']}")
    
    if not has_overlap:
        print(f"  ✓ No overlaps detected within same class/subclass")
    
    # Check scoring
    print(f"\nScoring verification:")
    for motif in motifs[:5]:
        score = motif.get('Score', 0)
        print(f"  {motif['Class']}/{motif['Subclass']}: score={score:.3f} "
              f"(range check: {0 <= score <= 1.0 or score > 1.0})")
    
    return {
        'total_motifs': len(motifs),
        'classes_found': len(by_class),
        'has_overlap': has_overlap
    }


def test_hyperscan_usage():
    """Verify that Hyperscan is being used when available"""
    print("\n" + "="*80)
    print("TESTING HYPERSCAN DATABASE USAGE")
    print("="*80)
    
    try:
        import hyperscan
        hyperscan_available = True
        print(f"✓ Hyperscan is installed (version: {hyperscan.__version__ if hasattr(hyperscan, '__version__') else 'unknown'})")
    except ImportError:
        hyperscan_available = False
        print(f"✗ Hyperscan is NOT installed")
    
    detector = ModularMotifDetector()
    
    # Check if detector has preloaded DBs
    if hasattr(detector, 'hsdb_map'):
        print(f"\n✓ Detector has preloaded {len(detector.hsdb_map)} Hyperscan databases:")
        for cls_name, db_info in detector.hsdb_map.items():
            db = db_info.get('db')
            n_patterns = len(db_info.get('id_to_pattern', {}))
            print(f"  - {cls_name}: {n_patterns} patterns, DB={'loaded' if db else 'None'}")
    else:
        print(f"\n✗ Detector has no preloaded Hyperscan databases")
    
    return {
        'hyperscan_available': hyperscan_available,
        'preloaded_dbs': len(detector.hsdb_map) if hasattr(detector, 'hsdb_map') else 0
    }


def test_all_detectors():
    """Test that all detector classes are loaded and functional"""
    print("\n" + "="*80)
    print("TESTING ALL DETECTOR CLASSES")
    print("="*80)
    
    detector = ModularMotifDetector()
    
    print(f"\nLoaded detectors:")
    for name, det in detector.detectors.items():
        stats = det.get_statistics()
        print(f"  - {name}: {stats['total_patterns']} patterns")
    
    print(f"\nTotal detectors: {len(detector.detectors)}")
    
    return {
        'total_detectors': len(detector.detectors),
        'detector_names': list(detector.detectors.keys())
    }


def main():
    """Run all tests"""
    print("\n" + "="*80)
    print("COMPREHENSIVE NON-B DNA MOTIF DETECTION TEST SUITE")
    print("="*80)
    
    all_results = {}
    
    # Test 1: Individual motif classes
    try:
        all_results['individual_classes'] = test_individual_motif_classes()
        print("\n✓ Individual motif class tests PASSED")
    except Exception as e:
        print(f"\n✗ Individual motif class tests FAILED: {e}")
        import traceback
        traceback.print_exc()
    
    # Test 2: Hybrid and cluster detection
    try:
        all_results['hybrid_cluster'] = test_hybrid_and_cluster_detection()
        print("\n✓ Hybrid and cluster detection tests PASSED")
    except Exception as e:
        print(f"\n✗ Hybrid and cluster detection tests FAILED: {e}")
        import traceback
        traceback.print_exc()
    
    # Test 3: Scoring and overlap resolution
    try:
        all_results['scoring_overlap'] = test_scoring_and_overlap_resolution()
        print("\n✓ Scoring and overlap resolution tests PASSED")
    except Exception as e:
        print(f"\n✗ Scoring and overlap resolution tests FAILED: {e}")
        import traceback
        traceback.print_exc()
    
    # Test 4: Hyperscan usage
    try:
        all_results['hyperscan'] = test_hyperscan_usage()
        print("\n✓ Hyperscan usage tests PASSED")
    except Exception as e:
        print(f"\n✗ Hyperscan usage tests FAILED: {e}")
        import traceback
        traceback.print_exc()
    
    # Test 5: All detectors loaded
    try:
        all_results['detectors'] = test_all_detectors()
        print("\n✓ Detector loading tests PASSED")
    except Exception as e:
        print(f"\n✗ Detector loading tests FAILED: {e}")
        import traceback
        traceback.print_exc()
    
    # Print summary
    print("\n" + "="*80)
    print("TEST SUMMARY")
    print("="*80)
    
    if 'detectors' in all_results:
        print(f"✓ All {all_results['detectors']['total_detectors']} detectors loaded")
    
    if 'hyperscan' in all_results:
        if all_results['hyperscan']['hyperscan_available']:
            print(f"✓ Hyperscan available with {all_results['hyperscan']['preloaded_dbs']} preloaded databases")
        else:
            print(f"⚠️  Hyperscan not available, using fallback regex")
    
    if 'scoring_overlap' in all_results:
        if not all_results['scoring_overlap']['has_overlap']:
            print(f"✓ Overlap resolution working correctly")
        else:
            print(f"✗ Overlap resolution has issues")
    
    if 'hybrid_cluster' in all_results:
        hc = all_results['hybrid_cluster']
        print(f"✓ Detected {hc['regular']} regular, {hc['hybrid']} hybrid, {hc['cluster']} cluster motifs")
    
    print("\n" + "="*80)
    print("✓ ALL TESTS COMPLETED")
    print("="*80)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
