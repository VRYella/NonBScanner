#!/usr/bin/env python3
"""
Test suite for RLoopDetector (QmRLFS-finder implementation)
"""

from detectors import RLoopDetector


def test_basic_initialization():
    """Test that RLoopDetector initializes correctly"""
    detector = RLoopDetector()
    assert detector is not None
    assert detector.get_motif_class_name() == "R-Loop"
    print("✓ Basic initialization test passed")


def test_model_1_detection():
    """Test Model 1 (G{3,}) pattern detection"""
    detector = RLoopDetector()
    
    # Sequence with clear Model 1 pattern (G{3,})
    # Multiple G-tracts of 3+ separated by short spacers
    test_seq = 'GGGAAAAGGGATCGGGAAATGGGCCC'
    motifs = detector.detect_motifs(test_seq, 'test_m1')
    
    # Should find at least one R-loop region
    assert len(motifs) >= 1, f"Expected at least 1 motif, found {len(motifs)}"
    
    # Check that it's classified as QmRLFS-m1
    for motif in motifs:
        assert motif['Subclass'] in ['QmRLFS-m1', 'QmRLFS-m2']
        assert motif['RIZ_G_Percent'] >= 50.0  # Should meet minimum G% threshold
        assert motif['RIZ_3G_Tracts'] >= 2  # Should have multiple G-tracts
    
    print(f"✓ Model 1 detection test passed ({len(motifs)} motifs found)")


def test_model_2_detection():
    """Test Model 2 (G{4,}) pattern detection"""
    detector = RLoopDetector()
    
    # Sequence with clear Model 2 pattern (G{4,})
    test_seq = 'GGGGATCGGGGAATCGGGGAAATGGGGATC'
    motifs = detector.detect_motifs(test_seq, 'test_m2')
    
    # Should find at least one R-loop region
    assert len(motifs) >= 1, f"Expected at least 1 motif, found {len(motifs)}"
    
    # Check that it's classified correctly
    for motif in motifs:
        assert motif['Subclass'] in ['QmRLFS-m1', 'QmRLFS-m2']
        assert motif['RIZ_G_Percent'] >= 50.0
        if motif['Subclass'] == 'QmRLFS-m2':
            assert motif['RIZ_4G_Tracts'] >= 2  # Should have 4+ G-tracts
    
    print(f"✓ Model 2 detection test passed ({len(motifs)} motifs found)")


def test_riz_rez_structure():
    """Test that RIZ and REZ are properly detected"""
    detector = RLoopDetector()
    
    # Create sequence with RIZ and potential REZ
    # RIZ: G-rich pattern, followed by linker, then REZ: more G-rich region
    riz_seq = 'GGGAATCGGGAATCGGGAAA'
    linker = 'A' * 60
    rez_seq = 'G' * 100  # Very G-rich REZ
    test_seq = riz_seq + linker + rez_seq
    
    motifs = detector.detect_motifs(test_seq, 'test_riz_rez')
    
    if len(motifs) > 0:
        # At least one motif should have REZ information
        has_rez = any(m.get('REZ_Length', 0) > 0 for m in motifs)
        
        for motif in motifs:
            # Check RIZ properties
            assert motif['RIZ_Start'] is not None
            assert motif['RIZ_End'] is not None
            assert motif['RIZ_Length'] > 0
            assert motif['RIZ_G_Percent'] >= 50.0
            
            # If REZ exists, check its properties
            if motif.get('REZ_Length', 0) > 0:
                assert motif['REZ_G_Percent'] >= 40.0  # REZ threshold is lower
                assert motif['Linker_Length'] >= 0
                print(f"  Found motif with REZ: RIZ_len={motif['RIZ_Length']}, "
                      f"REZ_len={motif['REZ_Length']}, Linker={motif['Linker_Length']}")
    
    print(f"✓ RIZ/REZ structure test passed ({len(motifs)} motifs found)")


def test_g_percent_threshold():
    """Test that G% thresholds are enforced"""
    detector = RLoopDetector()
    
    # Low G% sequence - should not be detected
    low_g_seq = 'ATATATATATAT' * 10
    motifs_low = detector.detect_motifs(low_g_seq, 'test_low_g')
    assert len(motifs_low) == 0, "Low G% sequence should not produce motifs"
    
    # High G% sequence - should be detected
    high_g_seq = 'GGGAATGGGAATGGGAATGGGAAT' * 2
    motifs_high = detector.detect_motifs(high_g_seq, 'test_high_g')
    
    for motif in motifs_high:
        assert motif['RIZ_G_Percent'] >= 50.0
    
    print(f"✓ G% threshold test passed (low_g={len(motifs_low)}, high_g={len(motifs_high)} motifs)")


def test_overlap_removal():
    """Test that overlapping motifs are properly handled"""
    detector = RLoopDetector()
    
    # Create sequence that could produce overlapping matches
    test_seq = 'GGGAAAGGGAAAGGGAAAGGGAAA' * 3
    motifs = detector.detect_motifs(test_seq, 'test_overlap')
    
    # Check for overlaps
    for i, motif1 in enumerate(motifs):
        for motif2 in motifs[i+1:]:
            # Motifs should not overlap
            overlap = not (motif1['End'] <= motif2['Start'] or 
                          motif1['Start'] >= motif2['End'])
            assert not overlap, f"Found overlapping motifs: {motif1['Start']}-{motif1['End']} and {motif2['Start']}-{motif2['End']}"
    
    print(f"✓ Overlap removal test passed ({len(motifs)} non-overlapping motifs)")


def test_hyperscan_fallback():
    """Test that detector works with and without hyperscan"""
    detector = RLoopDetector()
    
    test_seq = 'GGGAATCGGGAATCGGGAAA' * 2
    
    # Get results (will use hyperscan if available, regex otherwise)
    motifs = detector.detect_motifs(test_seq, 'test_fallback')
    
    # Should get consistent results regardless of hyperscan availability
    assert len(motifs) >= 1
    
    hs_status = "enabled" if detector.hs_db is not None else "disabled"
    print(f"✓ Hyperscan fallback test passed (hyperscan {hs_status}, {len(motifs)} motifs)")


def test_annotate_sequence():
    """Test the annotate_sequence method"""
    detector = RLoopDetector()
    
    test_seq = 'GGGAAAGGGAAAGGGAAA' * 2
    annotations = detector.annotate_sequence(test_seq, models=['qmrlfs_model_1'])
    
    # Should get annotation results
    assert isinstance(annotations, list)
    
    for ann in annotations:
        # Check that annotation has required fields
        assert 'riz_start' in ann
        assert 'riz_end' in ann
        assert 'riz_sequence' in ann
        assert 'riz_perc_g' in ann
        assert 'model' in ann
    
    print(f"✓ Annotate sequence test passed ({len(annotations)} annotations)")


def test_calculate_score():
    """Test score calculation"""
    detector = RLoopDetector()
    
    # High G% sequence
    high_g_seq = 'GGGAAAGGGAAAGGGAAA'
    score_high = detector.calculate_score(high_g_seq, None)
    
    # Low G% sequence
    low_g_seq = 'ATATATATATAT'
    score_low = detector.calculate_score(low_g_seq, None)
    
    # High G% should score higher
    assert score_high > score_low, f"High G% ({score_high}) should score higher than low G% ({score_low})"
    
    print(f"✓ Calculate score test passed (high={score_high:.3f}, low={score_low:.3f})")


def run_all_tests():
    """Run all test functions"""
    print("\n" + "="*60)
    print("Running RLoopDetector Test Suite")
    print("="*60 + "\n")
    
    tests = [
        test_basic_initialization,
        test_model_1_detection,
        test_model_2_detection,
        test_riz_rez_structure,
        test_g_percent_threshold,
        test_overlap_removal,
        test_hyperscan_fallback,
        test_annotate_sequence,
        test_calculate_score,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            print(f"✗ {test.__name__} FAILED: {e}")
            failed += 1
        except Exception as e:
            print(f"✗ {test.__name__} ERROR: {e}")
            failed += 1
    
    print("\n" + "="*60)
    print(f"Test Results: {passed} passed, {failed} failed")
    print("="*60 + "\n")
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    exit(0 if success else 1)
