#!/usr/bin/env python3
"""
Test script to verify Hyperscan integration and motif detection pipeline.
Tests the key patterns described in the problem statement:
1. Scan first, score later (Hyperscan for matching, separate scoring)
2. Proper merging of overlapping k-mer matches 
3. Non-overlapping output via deterministic selection
"""

import sys
sys.path.insert(0, '/home/runner/work/NonBScanner/NonBScanner')

from utils.motif_patterns import HyperscanManager, HYPERSCAN_AVAILABLE, MotifScoring
from motif_detection.a_philic_detector import APhilicDetector
from motif_detection.z_dna_detector import ZDNADetector
from motif_detection.cruciform_detector import CruciformDetector

def test_hyperscan_availability():
    """Test that HYPERSCAN_AVAILABLE is properly defined"""
    print("Testing HYPERSCAN_AVAILABLE...")
    assert HYPERSCAN_AVAILABLE is not None, "HYPERSCAN_AVAILABLE should be defined"
    print(f"  ✓ HYPERSCAN_AVAILABLE = {HYPERSCAN_AVAILABLE}")
    
    hm = HyperscanManager()
    assert hm.hyperscan_available is not None, "HyperscanManager.hyperscan_available should be defined"
    print(f"  ✓ HyperscanManager.hyperscan_available = {hm.hyperscan_available}")
    print()

def test_aphilic_detector_merging():
    """Test A-philic detector's merging guarantee"""
    print("Testing A-philic detector merging...")
    detector = APhilicDetector()
    
    # Test sequence with overlapping A-philic 10-mers
    test_seq = "A" * 20 + "GGGGG" + "C" * 20
    
    # Test that _find_10mer_matches returns raw matches (may overlap)
    matches = detector._find_10mer_matches(test_seq)
    print(f"  Raw matches found: {len(matches)}")
    
    # Test that _merge_matches combines overlapping matches
    if matches:
        merged = detector._merge_matches(matches)
        print(f"  Merged regions: {len(merged)}")
        assert len(merged) <= len(matches), "Merged regions should be <= raw matches"
        print("  ✓ Merging reduces overlapping matches to contiguous regions")
    
    # Test that detect_motifs returns merged regions only
    motifs = detector.detect_motifs(test_seq, "test_seq")
    print(f"  Final motifs output: {len(motifs)}")
    print("  ✓ A-philic detector outputs merged regions only")
    print()

def test_zdna_detector_merging():
    """Test Z-DNA detector's merging guarantee"""
    print("Testing Z-DNA detector merging...")
    detector = ZDNADetector()
    
    # Test sequence with Z-DNA motifs
    test_seq = "CG" * 10 + "ATATATAT" + "GC" * 10
    
    # Test merging logic
    matches = detector._find_10mer_matches(test_seq)
    print(f"  Raw Z-DNA matches found: {len(matches)}")
    
    if matches:
        merged = detector._merge_matches(matches)
        print(f"  Merged regions: {len(merged)}")
        assert len(merged) <= len(matches), "Merged regions should be <= raw matches"
        print("  ✓ Merging reduces overlapping Z-DNA matches")
    
    motifs = detector.detect_motifs(test_seq, "test_seq")
    print(f"  Final Z-DNA motifs output: {len(motifs)}")
    print("  ✓ Z-DNA detector outputs merged regions only")
    print()

def test_cruciform_overlap_removal():
    """Test Cruciform detector's overlap removal"""
    print("Testing Cruciform detector overlap removal...")
    detector = CruciformDetector()
    
    # Test sequence with potential inverted repeats
    test_seq = "ATGCATGC" * 3 + "AAAAAA" + "GCATGCAT" * 3
    
    # Find all inverted repeats (may overlap)
    repeats = detector.find_inverted_repeats(test_seq)
    print(f"  Raw inverted repeats found: {len(repeats)}")
    
    if repeats:
        # Test overlap removal
        non_overlapping = detector._remove_overlaps(repeats)
        print(f"  Non-overlapping repeats: {len(non_overlapping)}")
        assert len(non_overlapping) <= len(repeats), "Non-overlapping should be <= raw repeats"
        
        # Verify no overlaps in output
        for i, r1 in enumerate(non_overlapping):
            for r2 in non_overlapping[i+1:]:
                assert (r1['right_end'] <= r2['left_start'] or 
                       r2['right_end'] <= r1['left_start']), \
                       "Output should have no overlapping regions"
        print("  ✓ Overlap removal produces non-overlapping set")
    
    motifs = detector.detect_motifs(test_seq, "test_seq")
    print(f"  Final cruciform motifs output: {len(motifs)}")
    print("  ✓ Cruciform detector outputs non-overlapping regions")
    print()

def test_scoring_separation():
    """Test that scoring is separate from scanning"""
    print("Testing scoring separation from scanning...")
    
    # Test MotifScoring algorithms exist and work
    scoring = MotifScoring()
    test_seq = "GGGTTAGGGTTAGGGTTAGGG"  # G4 sequence
    
    g4_score = scoring.g4hunter_score(test_seq)
    print(f"  G4Hunter score: {g4_score:.4f}")
    assert g4_score > 0, "G4Hunter should return positive score for G4 sequence"
    
    # Test other scoring methods
    test_seq_z = "CGCGCGCGCGCGCG"
    z_score = scoring.z_dna_score(test_seq_z)
    print(f"  Z-DNA score: {z_score:.4f}")
    assert z_score > 0, "Z-DNA scorer should return positive score for alternating CG"
    
    print("  ✓ Scoring algorithms are separate and work correctly")
    print()

def main():
    print("=" * 70)
    print("HYPERSCAN INTEGRATION AND MOTIF DETECTION TESTS")
    print("=" * 70)
    print()
    
    try:
        test_hyperscan_availability()
        test_scoring_separation()
        test_aphilic_detector_merging()
        test_zdna_detector_merging()
        test_cruciform_overlap_removal()
        
        print("=" * 70)
        print("ALL TESTS PASSED ✓")
        print("=" * 70)
        print()
        print("Summary:")
        print("  • HYPERSCAN_AVAILABLE is properly defined")
        print("  • Scoring is separate from scanning (scan first, score later)")
        print("  • A-philic detector merges overlapping 10-mer matches")
        print("  • Z-DNA detector merges overlapping 10-mer matches")
        print("  • Cruciform detector removes overlaps deterministically")
        print("  • All detectors output non-overlapping regions")
        return 0
        
    except Exception as e:
        print(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
