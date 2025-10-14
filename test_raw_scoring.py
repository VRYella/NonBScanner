#!/usr/bin/env python3
"""
Test script for raw scoring system (without normalized scores)
Tests:
1. Raw scores are present and scientifically accurate
2. No Normalized_Score field in new detections
3. No overlaps within same subclass
4. Performance is maintained
5. Different motif classes have appropriate score ranges
"""

import time
from collections import defaultdict, Counter
from utils.modular_scanner import ModularMotifDetector

def test_raw_scores_present():
    """Test that raw scores are present and Normalized_Score is not added"""
    print("\n" + "=" * 80)
    print("TEST 1: Raw Scores Present (No Normalized_Score)")
    print("=" * 80)
    
    scanner = ModularMotifDetector()
    seq = 'GGGGTTGGGGTTGGGGTTGGGGAAAAAAAAAATTTTTTTTTTCGCGCGCGCGATATATATATAT'
    motifs = scanner.analyze_sequence(seq, 'test_raw_scores')
    
    print(f"Found {len(motifs)} motifs")
    
    has_score = 0
    has_normalized = 0
    
    for motif in motifs:
        if 'Score' in motif:
            has_score += 1
        if 'Normalized_Score' in motif and motif.get('Normalized_Score', 0) != 0:
            has_normalized += 1
            print(f"  WARNING: Found Normalized_Score={motif['Normalized_Score']} in {motif.get('Class')}")
    
    print(f"\nMotifs with Score field: {has_score}/{len(motifs)}")
    print(f"Motifs with non-zero Normalized_Score: {has_normalized}/{len(motifs)}")
    
    assert has_score == len(motifs), "All motifs should have Score field"
    assert has_normalized == 0, "No motifs should have non-zero Normalized_Score"
    
    print("✓ PASS: All motifs have raw scores only")
    return True


def test_score_ranges():
    """Test that different motif classes have appropriate score ranges"""
    print("\n" + "=" * 80)
    print("TEST 2: Class-Specific Score Ranges")
    print("=" * 80)
    
    scanner = ModularMotifDetector()
    
    # Create sequence with multiple motif types
    seq = ('GGGGTTGGGGTTGGGGTTGGGG' +  # G4
           'AAAAAAAAAAAAAAAAAA' +        # A-philic
           'CGCGCGCGCGCGCGCGCG' +        # Z-DNA
           'CCCCTTCCCCTTCCCCTTCCCC' +    # i-Motif
           'ATATATATATATATAT')           # Various
    
    motifs = scanner.analyze_sequence(seq, 'test_ranges')
    
    # Collect scores by class
    class_scores = defaultdict(list)
    for m in motifs:
        if m.get('Class') not in ['Hybrid', 'Non-B_DNA_Clusters']:  # Exclude derived classes
            class_scores[m.get('Class')].append(m.get('Score'))
    
    print("\nScore ranges by motif class:")
    for cls in sorted(class_scores.keys()):
        scores = class_scores[cls]
        if scores:
            print(f"  {cls:20s}: min={min(scores):8.3f}, max={max(scores):8.3f}, n={len(scores)}")
    
    # Verify expected ranges
    if 'G-Quadruplex' in class_scores:
        g4_scores = class_scores['G-Quadruplex']
        assert all(0 <= s <= 5 for s in g4_scores), "G4 scores should be in reasonable range"
    
    if 'Z-DNA' in class_scores:
        zdna_scores = class_scores['Z-DNA']
        assert any(s > 10 for s in zdna_scores), "Z-DNA should have high raw scores"
    
    if 'A-philic_DNA' in class_scores:
        aphi_scores = class_scores['A-philic_DNA']
        assert any(s > 5 for s in aphi_scores), "A-philic should have raw scores > 5"
    
    print("✓ PASS: Class-specific score ranges are scientifically accurate")
    return True


def test_no_overlaps_in_subclass():
    """Test that there are no overlaps within same subclass"""
    print("\n" + "=" * 80)
    print("TEST 3: No Overlaps Within Same Subclass")
    print("=" * 80)
    
    scanner = ModularMotifDetector()
    
    # Create a sequence with potential overlapping patterns
    seq = 'GGGGTTGGGGTTGGGGTTGGGG' * 5  # Repeated G4 patterns
    motifs = scanner.analyze_sequence(seq, 'test_overlaps')
    
    print(f"Found {len(motifs)} motifs")
    
    # Group by class/subclass
    groups = defaultdict(list)
    for m in motifs:
        if m.get('Class') not in ['Hybrid', 'Non-B_DNA_Clusters']:  # Exclude derived classes
            key = f"{m.get('Class')}/{m.get('Subclass')}"
            groups[key].append((m.get('Start'), m.get('End')))
    
    overlaps_found = 0
    for key, positions in groups.items():
        for i in range(len(positions)):
            for j in range(i+1, len(positions)):
                start1, end1 = positions[i]
                start2, end2 = positions[j]
                # Check if they overlap
                if not (end1 <= start2 or end2 <= start1):
                    overlap_start = max(start1, start2)
                    overlap_end = min(end1, end2)
                    overlap_len = overlap_end - overlap_start
                    print(f"  OVERLAP in {key}: [{start1}-{end1}] and [{start2}-{end2}] = {overlap_len}bp")
                    overlaps_found += 1
    
    print(f"\nTotal overlaps found: {overlaps_found}")
    assert overlaps_found == 0, "No overlaps should exist within same subclass"
    
    print("✓ PASS: No overlaps within same subclass")
    return True


def test_performance():
    """Test that performance is maintained"""
    print("\n" + "=" * 80)
    print("TEST 4: Performance Test")
    print("=" * 80)
    
    scanner = ModularMotifDetector()
    
    # Generate a 10KB test sequence
    seq = ('GGGGTTGGGGTTGGGGTTGGGG' + 
           'AAAAAAAAAAAA' + 
           'TTTTTTTTTTTT' +
           'CGCGCGCGCGCG' +
           'ATATATATATAT' +
           'CCCCTTCCCCTTCCCCTTCCCC') * 100
    
    seq_len = len(seq)
    print(f"Test sequence length: {seq_len} bp")
    
    start = time.time()
    motifs = scanner.analyze_sequence(seq, 'perf_test')
    elapsed = time.time() - start
    
    speed = seq_len / elapsed
    print(f"Time taken: {elapsed:.3f} seconds")
    print(f"Speed: {speed:.0f} bp/second")
    print(f"Motifs found: {len(motifs)}")
    
    # Performance should be at least 100 bp/second
    assert speed > 100, f"Performance too slow: {speed:.0f} bp/s (expected >100 bp/s)"
    
    print("✓ PASS: Performance is maintained")
    return True


def test_detailed_analysis():
    """Test that analysis is detailed with multiple classes and subclasses"""
    print("\n" + "=" * 80)
    print("TEST 5: Detailed Analysis")
    print("=" * 80)
    
    scanner = ModularMotifDetector()
    
    # Create diverse sequence
    seq = ('GGGGTTGGGGTTGGGGTTGGGG' +  # G4
           'AAAAAAAAAAAAAA' +            # Curved/A-philic
           'CGCGCGCGCGCGCG' +            # Z-DNA
           'CCCCTTCCCCTTCCCCTTCCCC' +    # i-Motif
           'ATATATATATATATAT' +          # Various
           'TTTTTTTTTTTTTT')             # Curved
    
    motifs = scanner.analyze_sequence(seq, 'test_detailed')
    
    # Count classes and subclasses
    classes = Counter(m.get('Class') for m in motifs)
    subclasses = Counter(f"{m.get('Class')}/{m.get('Subclass')}" for m in motifs)
    
    print(f"\nTotal motifs: {len(motifs)}")
    print(f"Unique classes: {len(classes)}")
    print(f"Unique subclasses: {len(subclasses)}")
    
    print("\nClasses detected:")
    for cls, count in classes.most_common():
        print(f"  {cls}: {count}")
    
    print(f"\nSubclasses detected: {len(subclasses)}")
    for subcls, count in sorted(subclasses.items())[:10]:
        print(f"  {subcls}: {count}")
    
    # Should detect multiple classes
    assert len(classes) >= 3, "Should detect at least 3 different motif classes"
    
    print("✓ PASS: Detailed analysis with multiple classes")
    return True


def main():
    """Run all tests"""
    print("=" * 80)
    print("RAW SCORING SYSTEM TESTS")
    print("Testing improvements: Raw scores, No normalization, No overlaps")
    print("=" * 80)
    
    tests = [
        test_raw_scores_present,
        test_score_ranges,
        test_no_overlaps_in_subclass,
        test_performance,
        test_detailed_analysis
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
        except AssertionError as e:
            print(f"✗ FAIL: {e}")
            failed += 1
        except Exception as e:
            print(f"✗ ERROR: {e}")
            failed += 1
    
    print("\n" + "=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"Tests passed: {passed}/{len(tests)}")
    print(f"Tests failed: {failed}/{len(tests)}")
    
    if failed == 0:
        print("✓ ALL TESTS PASSED")
        return 0
    else:
        print("✗ SOME TESTS FAILED")
        return 1


if __name__ == '__main__':
    exit(main())
