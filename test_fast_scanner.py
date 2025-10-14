#!/usr/bin/env python3
"""
Test Suite for Fast Scanner and Canonicalize Motif Utilities
=============================================================

This test suite validates:
1. Fast scanner G4 detection functionality
2. Canonicalize motif key normalization
3. Integration with existing motif detection pipeline

Author: Dr. Venkata Rajesh Yella
"""

import sys
import os

# Add project root to path for imports
project_root = os.path.dirname(os.path.abspath(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from utils.fast_scanner import (
    find_g4_candidates, extract_g4_features, 
    score_g4_candidates, fast_scan_and_score_g4
)
from utils.canonicalize_motif import canonicalize_motif


def test_fast_scanner():
    """Test fast_scanner G4 detection"""
    print("\n" + "=" * 80)
    print("Testing Fast Scanner G4 Detection")
    print("=" * 80)
    
    # Test 1: Human telomeric G4
    print("\n[Test 1] Human telomeric G-quadruplex")
    telomeric = "TTAGGGTTAGGGTTAGGGTTAGGG" * 2
    results = fast_scan_and_score_g4(telomeric, "telomeric_test")
    print(f"  Found {len(results)} G4 candidates")
    assert len(results) > 0, "Failed to detect telomeric G4"
    for r in results[:2]:
        print(f"  - Start: {r['Start']}, End: {r['End']}, Score: {r['Score']:.3f}")
        assert r['Class'] == 'G-Quadruplex', "Wrong class"
        assert r['Subclass'] == 'G4_Canonical', "Wrong subclass"
    print("  ✓ PASS")
    
    # Test 2: Canonical G4 pattern
    print("\n[Test 2] Canonical G4 with 4 tracts")
    canonical = "GGGGTTTTGGGGTTTTGGGGTTTTGGGG"
    results = fast_scan_and_score_g4(canonical, "canonical_test")
    print(f"  Found {len(results)} G4 candidates")
    assert len(results) > 0, "Failed to detect canonical G4"
    print(f"  - Score: {results[0]['Score']:.3f}")
    print("  ✓ PASS")
    
    # Test 3: No G4 sequence
    print("\n[Test 3] Sequence without G4 motifs")
    no_g4 = "ATATATATATAT"
    results = fast_scan_and_score_g4(no_g4, "no_g4_test")
    print(f"  Found {len(results)} G4 candidates")
    assert len(results) == 0, "False positive: detected G4 in non-G4 sequence"
    print("  ✓ PASS")
    
    # Test 4: Feature extraction
    print("\n[Test 4] G4 feature extraction")
    g4_seq = "GGGGTTGGGGTTGGGGTTGGGG"
    features = extract_g4_features(g4_seq)
    print(f"  Length: {features['length']}")
    print(f"  G-tracts: {features['g_tracts']}")
    print(f"  Loop lengths: {features['loop_lengths']}")
    print(f"  Mean tract: {features['mean_tract']:.2f}")
    print(f"  Mean loop: {features['mean_loop']:.2f}")
    assert len(features['g_tracts']) >= 4, "Failed to extract G-tracts"
    print("  ✓ PASS")
    
    print("\n" + "=" * 80)
    print("Fast Scanner Tests: ALL PASSED")
    print("=" * 80)
    return True


def test_canonicalize_motif():
    """Test canonicalize_motif normalization"""
    print("\n" + "=" * 80)
    print("Testing Canonicalize Motif Utility")
    print("=" * 80)
    
    # Test 1: Type/Subtype -> Class/Subclass
    print("\n[Test 1] Type/Subtype normalization")
    motif1 = {
        'Type': 'G-Quadruplex',
        'Subtype': 'Canonical',
        'Start': 10,
        'End': 30,
        'Score': 0.95
    }
    result1 = canonicalize_motif(motif1)
    print(f"  Input: Type={motif1['Type']}, Subtype={motif1['Subtype']}")
    print(f"  Output: Class={result1['Class']}, Subclass={result1['Subclass']}")
    assert result1['Class'] == 'G-Quadruplex', "Failed to map Type to Class"
    assert result1['Subclass'] == 'Canonical', "Failed to map Subtype to Subclass"
    assert result1['Length'] == 20, "Failed to calculate length"
    print("  ✓ PASS")
    
    # Test 2: Actual Score with space
    print("\n[Test 2] Actual Score normalization")
    motif2 = {
        'Class': 'i-Motif',
        'Subclass': 'Classic',
        'Start': 100,
        'End': 120,
        'Actual Score': 0.88
    }
    result2 = canonicalize_motif(motif2)
    print(f"  Input: Actual Score={motif2['Actual Score']}")
    print(f"  Output: Actual_Score={result2['Actual_Score']}")
    assert result2['Actual_Score'] == 0.88, "Failed to normalize Actual Score"
    assert result2['Score'] == 0.88, "Failed to propagate score"
    print("  ✓ PASS")
    
    # Test 3: Missing fields handling
    print("\n[Test 3] Missing fields with defaults")
    motif3 = {
        'Start': 50,
        'End': 70
    }
    result3 = canonicalize_motif(motif3)
    print(f"  Class: {result3['Class']} (should be 'Unknown')")
    print(f"  Subclass: {result3['Subclass']} (should be 'Other')")
    print(f"  Length: {result3['Length']}")
    assert result3['Class'] == 'Unknown', "Failed default Class"
    assert result3['Subclass'] == 'Other', "Failed default Subclass"
    assert result3['Length'] == 20, "Failed to calculate length from Start/End"
    print("  ✓ PASS")
    
    # Test 4: matched_seq handling
    print("\n[Test 4] matched_seq -> Motif mapping")
    motif4 = {
        'Class': 'G-Quadruplex',
        'Start': 0,
        'End': 20,
        'matched_seq': 'GGGGTTGGGGTTGGGGTTGG'
    }
    result4 = canonicalize_motif(motif4)
    print(f"  Input: matched_seq={motif4['matched_seq'][:20]}")
    print(f"  Output: Motif={result4['Motif'][:20]}")
    assert result4['Motif'] == motif4['matched_seq'], "Failed to map matched_seq to Motif"
    print("  ✓ PASS")
    
    # Test 5: All fields present
    print("\n[Test 5] Complete motif with all fields")
    motif5 = {
        'Class': 'Z-DNA',
        'Subclass': 'ZB_DNA',
        'Start': 200,
        'End': 250,
        'Length': 50,
        'Score': 0.92,
        'Actual_Score': 0.92,
        'Normalized_Score': 0.85,
        'Sequence_Name': 'test_sequence',
        'Motif': 'CGCGCGCG'
    }
    result5 = canonicalize_motif(motif5)
    print(f"  All fields preserved correctly")
    assert result5['Class'] == 'Z-DNA', "Failed to preserve Class"
    assert result5['Subclass'] == 'ZB_DNA', "Failed to preserve Subclass"
    assert result5['Score'] == 0.92, "Failed to preserve Score"
    assert result5['Normalized_Score'] == 0.85, "Failed to preserve Normalized_Score"
    print("  ✓ PASS")
    
    print("\n" + "=" * 80)
    print("Canonicalize Motif Tests: ALL PASSED")
    print("=" * 80)
    return True


def test_integration():
    """Test integration with existing pipeline"""
    print("\n" + "=" * 80)
    print("Testing Integration with Existing Pipeline")
    print("=" * 80)
    
    # Test 1: Fast scanner output format compatibility
    print("\n[Test 1] Fast scanner output format")
    seq = "GGGGTTTTGGGGTTTTGGGGTTTTGGGG"
    results = fast_scan_and_score_g4(seq, "integration_test")
    
    if results:
        motif = results[0]
        print(f"  Checking required fields...")
        required_fields = ['Class', 'Subclass', 'Start', 'End', 'Length', 'Score', 'Sequence_Name']
        for field in required_fields:
            assert field in motif, f"Missing required field: {field}"
            print(f"    - {field}: {motif[field]}")
        print("  ✓ PASS")
    else:
        print("  ⚠ WARNING: No G4 detected, skipping format check")
    
    # Test 2: Canonicalization of fast scanner output
    print("\n[Test 2] Canonicalize fast scanner output")
    if results:
        normalized = canonicalize_motif(results[0])
        print(f"  Original keys: {list(results[0].keys())}")
        print(f"  Normalized keys: {list(normalized.keys())}")
        assert 'Class' in normalized, "Missing Class after normalization"
        assert 'Subclass' in normalized, "Missing Subclass after normalization"
        print("  ✓ PASS")
    
    print("\n" + "=" * 80)
    print("Integration Tests: ALL PASSED")
    print("=" * 80)
    return True


def main():
    """Run all tests"""
    print("\n" + "=" * 80)
    print("FAST SCANNER AND CANONICALIZE MOTIF TEST SUITE")
    print("=" * 80)
    
    all_passed = True
    
    try:
        all_passed = test_fast_scanner() and all_passed
    except Exception as e:
        print(f"\n✗ Fast Scanner tests FAILED: {e}")
        all_passed = False
    
    try:
        all_passed = test_canonicalize_motif() and all_passed
    except Exception as e:
        print(f"\n✗ Canonicalize Motif tests FAILED: {e}")
        all_passed = False
    
    try:
        all_passed = test_integration() and all_passed
    except Exception as e:
        print(f"\n✗ Integration tests FAILED: {e}")
        all_passed = False
    
    # Summary
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    if all_passed:
        print("✓ ALL TESTS PASSED")
        print("=" * 80 + "\n")
        return 0
    else:
        print("✗ SOME TESTS FAILED")
        print("=" * 80 + "\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
