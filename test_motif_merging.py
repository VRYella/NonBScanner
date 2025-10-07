#!/usr/bin/env python3
"""
Test and demonstration script for A-philic and Z-DNA detector merging behavior.

This script validates that:
1. Overlapping/adjacent 10-mer matches are ALWAYS merged into single regions
2. No duplicate or split reporting occurs
3. The output includes all required fields (start, end, length, sequence, n_10mers, score)
4. Both detect_motifs() and annotate_sequence() produce merged output

Usage:
    python test_motif_merging.py
"""

import sys
from typing import List, Dict, Any

from motif_detection.a_philic_detector import APhilicDetector
from motif_detection.z_dna_detector import ZDNADetector


def validate_merged_output(regions: List[Dict[str, Any]], detector_name: str, 
                          expected_min_regions: int = 1, expected_max_regions: int = 1) -> bool:
    """
    Validate that the output represents properly merged regions.
    
    Args:
        regions: List of region dictionaries from annotate_sequence()
        detector_name: Name of the detector being tested
        expected_min_regions: Minimum expected number of regions
        expected_max_regions: Maximum expected number of regions
    
    Returns:
        True if validation passes, False otherwise
    """
    print(f"\n  Validating {detector_name} output...")
    
    # Check we have the expected number of regions
    num_regions = len(regions)
    if not (expected_min_regions <= num_regions <= expected_max_regions):
        print(f"    ✗ FAIL: Expected {expected_min_regions}-{expected_max_regions} regions, got {num_regions}")
        return False
    print(f"    ✓ Found {num_regions} region(s) as expected")
    
    # Validate each region has required fields
    required_fields = ['start', 'end', 'length', 'n_10mers']
    for i, region in enumerate(regions):
        for field in required_fields:
            if field not in region:
                print(f"    ✗ FAIL: Region {i} missing required field '{field}'")
                return False
        
        # Validate that length matches coordinates
        expected_len = region['end'] - region['start']
        if region['length'] != expected_len:
            print(f"    ✗ FAIL: Region {i} length mismatch: {region['length']} != {expected_len}")
            return False
        
        # Validate that we have contributing 10-mers
        if region['n_10mers'] < 1:
            print(f"    ✗ FAIL: Region {i} has no contributing 10-mers")
            return False
        
        print(f"    ✓ Region {i}: [{region['start']}, {region['end']}), "
              f"length={region['length']}, n_10mers={region['n_10mers']}")
    
    # Check for overlapping regions (there should be none)
    for i in range(len(regions) - 1):
        if regions[i]['end'] > regions[i+1]['start']:
            print(f"    ✗ FAIL: Regions {i} and {i+1} overlap!")
            return False
    
    print(f"    ✓ All validation checks passed")
    return True


def test_aphilic_overlapping_10mers():
    """
    Test A-philic detector with a sequence containing overlapping 10-mer matches.
    Expected: Single merged region containing all overlapping matches.
    """
    print("=" * 70)
    print("TEST 1: A-philic Detector - Overlapping 10-mers")
    print("=" * 70)
    
    detector = APhilicDetector()
    
    # Sequence with multiple overlapping A-philic 10-mers
    # AGGGGGGGGG at position 0 (score: 2.702)
    # GGGGGGGGGC at position 1 (score: 2.517)
    # etc. - all overlapping by 9bp
    test_seq = "AGGGGGGGGGCCCCCCCCCTAGGGGGGGGC"
    print(f"Sequence: {test_seq}")
    print(f"Length: {len(test_seq)} bp")
    
    # Test annotate_sequence
    print("\n1. Testing annotate_sequence()...")
    annotations = detector.annotate_sequence(test_seq)
    
    if not validate_merged_output(annotations, "A-philic annotate_sequence", 
                                  expected_min_regions=1, expected_max_regions=1):
        return False
    
    # Test detect_motifs
    print("\n2. Testing detect_motifs()...")
    motifs = detector.detect_motifs(test_seq, "test_seq")
    
    if len(motifs) != len(annotations):
        print(f"  ✗ FAIL: detect_motifs returned {len(motifs)} motifs, "
              f"but annotate_sequence returned {len(annotations)} regions")
        return False
    
    for motif in motifs:
        print(f"  ✓ Motif: {motif['ID']}, start={motif['Start']}, "
              f"end={motif['End']}, length={motif['Length']}, score={motif['Score']}")
    
    print("\n✓ TEST 1 PASSED: A-philic detector correctly merges overlapping 10-mers")
    return True


def test_zdna_overlapping_10mers():
    """
    Test Z-DNA detector with a sequence containing overlapping 10-mer matches.
    Expected: Single merged region containing all overlapping matches.
    """
    print("\n" + "=" * 70)
    print("TEST 2: Z-DNA Detector - Overlapping 10-mers")
    print("=" * 70)
    
    detector = ZDNADetector()
    
    # Sequence with multiple overlapping Z-DNA 10-mers
    # GCGCGCGCGC at positions 0, 2, 4, 6, 8, 10 (score: 63.0 each)
    # CGCGCGCGCG at positions 1, 3, 5, 7, 9 (score: 63.0 each)
    test_seq = "GCGCGCGCGCGCGCGCGCGC"
    print(f"Sequence: {test_seq}")
    print(f"Length: {len(test_seq)} bp")
    
    # Test annotate_sequence
    print("\n1. Testing annotate_sequence()...")
    annotations = detector.annotate_sequence(test_seq)
    
    if not validate_merged_output(annotations, "Z-DNA annotate_sequence", 
                                  expected_min_regions=1, expected_max_regions=1):
        return False
    
    # Test detect_motifs
    print("\n2. Testing detect_motifs()...")
    motifs = detector.detect_motifs(test_seq, "test_seq")
    
    if len(motifs) != len(annotations):
        print(f"  ✗ FAIL: detect_motifs returned {len(motifs)} motifs, "
              f"but annotate_sequence returned {len(annotations)} regions")
        return False
    
    for motif in motifs:
        print(f"  ✓ Motif: {motif['ID']}, start={motif['Start']}, "
              f"end={motif['End']}, length={motif['Length']}, score={motif['Score']}")
    
    print("\n✓ TEST 2 PASSED: Z-DNA detector correctly merges overlapping 10-mers")
    return True


def test_aphilic_separated_regions():
    """
    Test A-philic detector with separated 10-mer clusters.
    Expected: Multiple merged regions (one per cluster), not merged together.
    """
    print("\n" + "=" * 70)
    print("TEST 3: A-philic Detector - Separated Regions")
    print("=" * 70)
    
    detector = APhilicDetector()
    
    # Two clusters of A-philic 10-mers separated by non-matching sequence
    test_seq = "AGGGGGGGGGCCCCCAAAAAAAAAAAAAAGGGGGGGGGCCCCCC"
    print(f"Sequence: {test_seq}")
    print(f"Length: {len(test_seq)} bp")
    
    # Test annotate_sequence - should get 2 separate regions
    print("\n1. Testing annotate_sequence()...")
    annotations = detector.annotate_sequence(test_seq)
    
    if not validate_merged_output(annotations, "A-philic annotate_sequence", 
                                  expected_min_regions=2, expected_max_regions=2):
        return False
    
    print("\n✓ TEST 3 PASSED: A-philic detector correctly handles separated regions")
    return True


def test_empty_sequence():
    """
    Test both detectors with sequences containing no matches.
    Expected: Empty output (no regions, no motifs).
    """
    print("\n" + "=" * 70)
    print("TEST 4: Empty Sequences (No Matches)")
    print("=" * 70)
    
    aphilic_detector = APhilicDetector()
    zdna_detector = ZDNADetector()
    
    # Sequence with no A-philic or Z-DNA motifs
    test_seq = "AAAAAAAAAAAAAAAAAAAAAA"
    print(f"Sequence: {test_seq}")
    
    # Test A-philic
    print("\n1. Testing A-philic detector...")
    annotations = aphilic_detector.annotate_sequence(test_seq)
    motifs = aphilic_detector.detect_motifs(test_seq, "test")
    
    if len(annotations) != 0 or len(motifs) != 0:
        print(f"  ✗ FAIL: Expected empty output, got {len(annotations)} annotations, {len(motifs)} motifs")
        return False
    print(f"  ✓ A-philic: {len(annotations)} annotations, {len(motifs)} motifs (as expected)")
    
    # Test Z-DNA
    print("\n2. Testing Z-DNA detector...")
    annotations = zdna_detector.annotate_sequence(test_seq)
    motifs = zdna_detector.detect_motifs(test_seq, "test")
    
    if len(annotations) != 0 or len(motifs) != 0:
        print(f"  ✗ FAIL: Expected empty output, got {len(annotations)} annotations, {len(motifs)} motifs")
        return False
    print(f"  ✓ Z-DNA: {len(annotations)} annotations, {len(motifs)} motifs (as expected)")
    
    print("\n✓ TEST 4 PASSED: Both detectors correctly handle sequences with no matches")
    return True


def main():
    """Run all tests and report results."""
    print("\n" + "=" * 70)
    print("MOTIF DETECTOR MERGING BEHAVIOR - TEST SUITE")
    print("=" * 70)
    print("\nThis test suite validates that A-philic and Z-DNA detectors:")
    print("  1. Always merge overlapping/adjacent 10-mer matches")
    print("  2. Never output duplicate or split 10-mers")
    print("  3. Include all required fields in output")
    print("  4. Produce consistent results from detect_motifs() and annotate_sequence()")
    
    tests = [
        ("A-philic overlapping 10-mers", test_aphilic_overlapping_10mers),
        ("Z-DNA overlapping 10-mers", test_zdna_overlapping_10mers),
        ("A-philic separated regions", test_aphilic_separated_regions),
        ("Empty sequences", test_empty_sequence),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            results.append((test_name, test_func()))
        except Exception as e:
            print(f"\n✗ TEST FAILED with exception: {e}")
            import traceback
            traceback.print_exc()
            results.append((test_name, False))
    
    # Print summary
    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✓ PASSED" if result else "✗ FAILED"
        print(f"{status}: {test_name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 ALL TESTS PASSED! Detectors correctly merge overlapping 10-mers.")
        return 0
    else:
        print(f"\n❌ {total - passed} test(s) failed.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
