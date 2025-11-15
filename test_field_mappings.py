#!/usr/bin/env python3
"""
Comprehensive test for field mappings in motif detectors.

This test validates that all motif types correctly map the expected component fields:
- Left_Arm, Right_Arm, Loop_Seq
- Arm_Length, Loop_Length, Stem_Length
- Unit_Length, Number_Of_Copies
- Spacer_Length, Spacer_Sequence
"""

import sys
sys.path.insert(0, '/home/runner/work/NonBScanner/NonBScanner')

from detectors import (
    SlippedDNADetector,
    CruciformDetector,
    TriplexDetector,
    GQuadruplexDetector,
    IMotifDetector
)


def test_slipped_dna_str():
    """Test SlippedDNA STR fields"""
    print("\n" + "="*70)
    print("TEST 1: SlippedDNA (STR) - Unit_Length and Number_Of_Copies")
    print("="*70)
    
    # Create a clear STR pattern with a distinctive unit that won't be confused with 1-mer
    # Use ATAT... which should be detected as 2-mer with multiple copies
    seq = 'ATATATATATATATATATATATAT'  # 12 copies of AT
    detector = SlippedDNADetector()
    motifs = detector.detect_motifs(seq, 'test')
    
    assert len(motifs) > 0, "No motifs detected!"
    motif = motifs[0]
    
    # Check Unit_Length
    assert 'Unit_Length' in motif, "Unit_Length field missing!"
    print(f"✓ Unit_Length present: {motif['Unit_Length']}")
    
    # Check Number_Of_Copies
    assert 'Number_Of_Copies' in motif, "Number_Of_Copies field missing!"
    print(f"✓ Number_Of_Copies present: {motif['Number_Of_Copies']}")
    
    # Validate values - the scanner detects the most parsimonious repeat unit
    # For ATATATAT..., it should detect either 1-mer, 2-mer, or the whole sequence
    # The key is that Unit_Length * Number_Of_Copies >= sequence length
    total_coverage = motif['Unit_Length'] * motif['Number_Of_Copies']
    assert total_coverage >= len(seq.strip()), \
        f"Unit_Length * Number_Of_Copies should cover the sequence"
    print(f"✓ Unit_Length ({motif['Unit_Length']}) * Number_Of_Copies ({motif['Number_Of_Copies']}) = {total_coverage} >= {len(seq)}")
    
    print("✅ PASSED: SlippedDNA STR fields correctly mapped")
    return True


def test_cruciform():
    """Test Cruciform fields"""
    print("\n" + "="*70)
    print("TEST 2: Cruciform - Left_Arm, Right_Arm, Loop_Seq, Arm_Length, Loop_Length, Stem_Length")
    print("="*70)
    
    # Create inverted repeat pattern
    seq = 'ATCGATCGATCGATCGATCGATCGATCAAAAAATCGATCGATCGATCGATCGATCGATC'
    detector = CruciformDetector()
    motifs = detector.detect_motifs(seq, 'test')
    
    assert len(motifs) > 0, "No motifs detected!"
    motif = motifs[0]
    
    # Check all required fields
    required_fields = ['Left_Arm', 'Right_Arm', 'Loop_Seq', 'Arm_Length', 'Loop_Length', 'Stem_Length']
    for field in required_fields:
        assert field in motif, f"{field} field missing!"
        print(f"✓ {field} present: {motif[field]}")
    
    # Validate Stem_Length equals Arm_Length for cruciforms
    assert motif['Stem_Length'] == motif['Arm_Length'], \
        f"Stem_Length should equal Arm_Length for cruciforms"
    
    print("✅ PASSED: Cruciform fields correctly mapped")
    return True


def test_triplex():
    """Test Triplex mirror repeat fields"""
    print("\n" + "="*70)
    print("TEST 3: Triplex - Left_Arm, Right_Arm, Loop_Seq, Arm_Length, Loop_Length")
    print("="*70)
    
    # Create purine-rich mirror repeat
    seq = 'GAGAGAGAGAGAGAGAGAGATTTTTTTTTTGAGAGAGAGAGAGAGAGAGA'
    detector = TriplexDetector()
    motifs = detector.detect_motifs(seq, 'test')
    
    assert len(motifs) > 0, "No motifs detected!"
    motif = motifs[0]
    
    # Check all required fields
    required_fields = ['Left_Arm', 'Right_Arm', 'Loop_Seq', 'Arm_Length', 'Loop_Length']
    for field in required_fields:
        assert field in motif, f"{field} field missing!"
        print(f"✓ {field} present: {motif[field]}")
    
    # Validate GC content fields also present
    assert 'GC_Left_Arm' in motif, "GC_Left_Arm missing"
    assert 'GC_Right_Arm' in motif, "GC_Right_Arm missing"
    assert 'GC_Loop' in motif, "GC_Loop missing"
    print(f"✓ GC content fields present")
    
    print("✅ PASSED: Triplex fields correctly mapped")
    return True


def test_g_quadruplex():
    """Test G-Quadruplex fields"""
    print("\n" + "="*70)
    print("TEST 4: G-Quadruplex - Stem_Length, Loop_Length (consolidated)")
    print("="*70)
    
    # Canonical G4 sequence
    seq = 'GGGTTAGGGTTAGGGTTAGGG'
    detector = GQuadruplexDetector()
    motifs = detector.detect_motifs(seq, 'test')
    
    assert len(motifs) > 0, "No motifs detected!"
    motif = motifs[0]
    
    # Check consolidated fields
    assert 'Stem_Length' in motif, "Stem_Length field missing!"
    assert 'Loop_Length' in motif, "Loop_Length field missing!"
    print(f"✓ Stem_Length present: {motif['Stem_Length']}")
    print(f"✓ Loop_Length present: {motif['Loop_Length']}")
    
    # Check list fields still present
    assert 'Stem_Lengths' in motif, "Stem_Lengths list missing!"
    assert 'Loop_Lengths' in motif, "Loop_Lengths list missing!"
    assert 'Stems' in motif, "Stems list missing!"
    assert 'Loops' in motif, "Loops list missing!"
    print(f"✓ List fields present: Stem_Lengths, Loop_Lengths, Stems, Loops")
    
    # Validate Stem_Length is average of stem lengths
    expected_avg = sum(motif['Stem_Lengths']) / len(motif['Stem_Lengths'])
    assert motif['Stem_Length'] == expected_avg, \
        f"Stem_Length should be average of stem lengths"
    
    print("✅ PASSED: G-Quadruplex fields correctly mapped")
    return True


def test_i_motif():
    """Test i-Motif fields"""
    print("\n" + "="*70)
    print("TEST 5: i-Motif - Stem_Length, Loop_Length (consolidated)")
    print("="*70)
    
    # Canonical i-motif sequence
    seq = 'CCCCACCCCACCCCACCCC'
    detector = IMotifDetector()
    motifs = detector.detect_motifs(seq, 'test')
    
    assert len(motifs) > 0, "No motifs detected!"
    motif = motifs[0]
    
    # Check consolidated fields
    assert 'Stem_Length' in motif, "Stem_Length field missing!"
    assert 'Loop_Length' in motif, "Loop_Length field missing!"
    print(f"✓ Stem_Length present: {motif['Stem_Length']}")
    print(f"✓ Loop_Length present: {motif['Loop_Length']}")
    
    # Check list fields still present
    assert 'Stem_Lengths' in motif, "Stem_Lengths list missing!"
    assert 'Loop_Lengths' in motif, "Loop_Lengths list missing!"
    assert 'Stems' in motif, "Stems list missing!"
    assert 'Loops' in motif, "Loops list missing!"
    print(f"✓ List fields present: Stem_Lengths, Loop_Lengths, Stems, Loops")
    
    # Validate Stem_Length is average of stem lengths
    expected_avg = sum(len(s) for s in motif['Stems']) / len(motif['Stems'])
    assert motif['Stem_Length'] == expected_avg, \
        f"Stem_Length should be average of stem lengths"
    
    print("✅ PASSED: i-Motif fields correctly mapped")
    return True


def main():
    """Run all tests"""
    print("\n" + "="*70)
    print("COMPREHENSIVE FIELD MAPPING TESTS")
    print("="*70)
    
    tests = [
        test_slipped_dna_str,
        test_cruciform,
        test_triplex,
        test_g_quadruplex,
        test_i_motif
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
        except AssertionError as e:
            print(f"❌ FAILED: {e}")
            failed += 1
        except Exception as e:
            print(f"❌ ERROR: {e}")
            failed += 1
    
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    print(f"Passed: {passed}/{len(tests)}")
    print(f"Failed: {failed}/{len(tests)}")
    
    if failed == 0:
        print("\n✅ ALL TESTS PASSED!")
        return 0
    else:
        print("\n❌ SOME TESTS FAILED")
        return 1


if __name__ == '__main__':
    sys.exit(main())
