#!/usr/bin/env python3
"""
Test script to verify all motif detectors report non-overlapping motifs.

This script validates that:
1. No detector reports overlapping motifs within the same class/subclass
2. All detectors properly handle overlap removal
3. The output is consistent and reliable

Usage:
    python test_non_overlapping_motifs.py
"""

import sys
from typing import List, Dict, Any


def check_overlaps(motifs: List[Dict[str, Any]], detector_name: str) -> bool:
    """
    Check if any motifs overlap within the same class/subclass.
    Returns True if no overlaps found, False otherwise.
    """
    if len(motifs) <= 1:
        return True
    
    # Group by class and subclass
    from collections import defaultdict
    groups = defaultdict(list)
    for motif in motifs:
        # Handle both dict formats (with Start/End or start/end)
        start = motif.get('Start', motif.get('start', 0))
        end = motif.get('End', motif.get('end', 0))
        subclass = motif.get('Subclass', motif.get('class_name', 'unknown'))
        
        groups[subclass].append({
            'start': start,
            'end': end,
            'motif': motif
        })
    
    overlaps_found = False
    for subclass, group_motifs in groups.items():
        # Check for overlaps within this subclass
        for i in range(len(group_motifs) - 1):
            for j in range(i + 1, len(group_motifs)):
                m1, m2 = group_motifs[i], group_motifs[j]
                # Two motifs overlap if their regions overlap
                if not (m1['end'] <= m2['start'] or m2['end'] <= m1['start']):
                    print(f"  ✗ OVERLAP in {detector_name} ({subclass}): "
                          f"motif {i} [{m1['start']}, {m1['end']}) overlaps with "
                          f"motif {j} [{m2['start']}, {m2['end']})")
                    overlaps_found = True
    
    return not overlaps_found


def test_triplex_detector():
    """Test Triplex detector for non-overlapping motifs"""
    print("\n" + "=" * 70)
    print("TEST: Triplex Detector")
    print("=" * 70)
    
    from motif_detection.triplex_detector import TriplexDetector
    
    # Test sequence with GAA repeat (should trigger sticky DNA pattern)
    seq = 'GAAGAAGAAGAAGAAGAAGAAGAAGAA'
    detector = TriplexDetector()
    motifs = detector.detect_motifs(seq, 'test_triplex')
    
    print(f"Sequence: {seq}")
    print(f"Found {len(motifs)} motif(s)")
    
    for i, m in enumerate(motifs):
        print(f"  {i+1}. {m['Subclass']}: [{m['Start']}, {m['End']}) len={m['Length']} score={m['Score']}")
    
    result = check_overlaps(motifs, 'Triplex')
    if result:
        print("  ✓ No overlaps detected")
    
    return result


def test_cruciform_detector():
    """Test Cruciform detector for non-overlapping motifs"""
    print("\n" + "=" * 70)
    print("TEST: Cruciform Detector")
    print("=" * 70)
    
    from motif_detection.cruciform_detector import CruciformDetector
    
    # Test sequence with inverted repeats
    seq = 'AAAAAAAATTTTTTTTAGCTAGCT' * 2
    detector = CruciformDetector()
    motifs = detector.detect_motifs(seq, 'test_cruciform')
    
    print(f"Sequence length: {len(seq)} bp")
    print(f"Found {len(motifs)} motif(s)")
    
    for i, m in enumerate(motifs[:5]):  # Show first 5
        print(f"  {i+1}. {m['Subclass']}: [{m['Start']}, {m['End']}) len={m['Length']} score={m['Score']}")
    
    result = check_overlaps(motifs, 'Cruciform')
    if result:
        print("  ✓ No overlaps detected")
    
    return result


def test_rloop_detector():
    """Test R-loop detector for non-overlapping motifs"""
    print("\n" + "=" * 70)
    print("TEST: R-Loop Detector")
    print("=" * 70)
    
    from motif_detection.r_loop_detector import RLoopDetector
    
    # Test sequence with GC-rich regions
    seq = 'GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC'
    detector = RLoopDetector()
    motifs = detector.detect_motifs(seq, 'test_rloop')
    
    print(f"Sequence: {seq}")
    print(f"Found {len(motifs)} motif(s)")
    
    for i, m in enumerate(motifs):
        print(f"  {i+1}. {m['Subclass']}: [{m['Start']}, {m['End']}) len={m['Length']} score={m['Score']}")
    
    result = check_overlaps(motifs, 'R-loop')
    if result:
        print("  ✓ No overlaps detected")
    
    return result


def test_g4_detector():
    """Test G-Quadruplex detector for non-overlapping motifs"""
    print("\n" + "=" * 70)
    print("TEST: G-Quadruplex Detector")
    print("=" * 70)
    
    from motif_detection.g_quadruplex_detector import GQuadruplexDetector
    
    # Test sequence with G4 patterns
    seq = 'GGGAGGGAGGGAGGGATCGGGATCGGG'
    detector = GQuadruplexDetector()
    results = detector.annotate_sequence(seq)
    
    print(f"Sequence: {seq}")
    print(f"Found {len(results)} region(s)")
    
    for i, r in enumerate(results):
        print(f"  {i+1}. {r['class_name']}: [{r['start']}, {r['end']}) len={r['length']} score={r['score']:.3f}")
    
    result = check_overlaps(results, 'G-Quadruplex')
    if result:
        print("  ✓ No overlaps detected")
    
    return result


def test_slipped_dna_detector():
    """Test Slipped DNA detector for non-overlapping motifs"""
    print("\n" + "=" * 70)
    print("TEST: Slipped DNA Detector")
    print("=" * 70)
    
    from motif_detection.slipped_dna_detector import SlippedDNADetector
    
    # Test sequence with STRs
    seq = 'AAAAAAAAAAATTTTTTTTTTTGCGCGCGCGC'
    detector = SlippedDNADetector()
    results = detector.annotate_sequence(seq)
    
    print(f"Sequence: {seq}")
    print(f"Found {len(results)} region(s)")
    
    for i, r in enumerate(results):
        print(f"  {i+1}. {r['class_name']}: [{r['start']}, {r['end']}) len={r['length']} score={r['score']:.3f}")
    
    result = check_overlaps(results, 'Slipped DNA')
    if result:
        print("  ✓ No overlaps detected")
    
    return result


def test_curved_dna_detector():
    """Test Curved DNA detector for non-overlapping motifs"""
    print("\n" + "=" * 70)
    print("TEST: Curved DNA Detector")
    print("=" * 70)
    
    from motif_detection.curved_dna_detector import CurvedDNADetector
    from collections import defaultdict
    
    # Test sequence with A-tracts
    seq = 'AAAAAAAATGCAAAAAAAATGCAAAAAAAATGCAAAAAAA'
    detector = CurvedDNADetector()
    motifs = detector.detect_motifs(seq, 'test_curved')
    
    print(f"Sequence length: {len(seq)} bp")
    print(f"Found {len(motifs)} motif(s)")
    
    for i, m in enumerate(motifs[:5]):  # Show first 5
        print(f"  {i+1}. {m['Subclass']}: [{m['Start']}, {m['End']}) len={m['Length']} score={m['Score']}")
    
    # Check overlaps within each subclass (not across subclasses)
    groups = defaultdict(list)
    for m in motifs:
        groups[m['Subclass']].append(m)
    
    all_good = True
    for subclass, group in groups.items():
        overlaps = False
        for i in range(len(group) - 1):
            for j in range(i + 1, len(group)):
                m1, m2 = group[i], group[j]
                if not (m1['End'] <= m2['Start'] or m2['End'] <= m1['Start']):
                    print(f"  ✗ OVERLAP in {subclass}: motif {i} overlaps with motif {j}")
                    overlaps = False
                    all_good = False
        
    if all_good:
        print("  ✓ No overlaps detected within subclasses")
    
    return all_good


def main():
    """Run all tests and report results."""
    print("\n" + "=" * 70)
    print("NON-OVERLAPPING MOTIF DETECTION - TEST SUITE")
    print("=" * 70)
    print("\nThis test suite validates that all motif detectors:")
    print("  1. Report non-overlapping motifs within the same class/subclass")
    print("  2. Properly handle overlap removal")
    print("  3. Produce consistent and reliable output")
    
    tests = [
        ("Triplex Detector", test_triplex_detector),
        ("Cruciform Detector", test_cruciform_detector),
        ("R-Loop Detector", test_rloop_detector),
        ("G-Quadruplex Detector", test_g4_detector),
        ("Slipped DNA Detector", test_slipped_dna_detector),
        ("Curved DNA Detector", test_curved_dna_detector),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result, None))
        except Exception as e:
            print(f"\n  ✗ TEST FAILED with exception: {e}")
            import traceback
            traceback.print_exc()
            results.append((test_name, False, str(e)))
    
    # Print summary
    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)
    
    passed = sum(1 for _, result, _ in results if result)
    total = len(results)
    
    for test_name, result, error in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {test_name}")
        if error:
            print(f"      Error: {error}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n✓ All tests passed! All detectors report non-overlapping motifs.")
        return 0
    else:
        print(f"\n✗ {total - passed} test(s) failed.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
