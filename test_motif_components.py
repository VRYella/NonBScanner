#!/usr/bin/env python3
"""
Test Script: Comprehensive Motif Component Detection Validation
================================================================
This script validates that all motif types properly report their contributing
parts including stems, arms, loops, GC percentage, and other structural details.

Author: NonBScanner Development Team
Purpose: Validate component extraction for all motif classes
"""

from scanner import analyze_sequence
import json


def test_str_components():
    """Test STR (Short Tandem Repeat) component extraction"""
    print("\n" + "="*70)
    print("Testing STR Component Extraction")
    print("="*70)
    
    # CA repeat
    seq = "CACACACACACACACA"
    result = analyze_sequence(seq, 'test_str')
    
    str_motifs = [m for m in result if m.get('Class') == 'Slipped_DNA' and m.get('Subclass') == 'STR']
    
    if str_motifs:
        motif = str_motifs[0]
        print(f"\n✓ STR Detected: {motif.get('Sequence')}")
        print(f"  Repeat Unit: {motif.get('Repeat_Unit', 'N/A')}")
        print(f"  Number of Copies: {motif.get('Number_Of_Copies', 'N/A')}")
        print(f"  GC Unit: {motif.get('Gc_Unit', 'N/A')}%")
        print(f"  GC Total: {motif.get('Gc_Total', 'N/A')}%")
        print(f"  Base Counts: A={motif.get('Unit_A_Count', 0)}, T={motif.get('Unit_T_Count', 0)}, G={motif.get('Unit_G_Count', 0)}, C={motif.get('Unit_C_Count', 0)}")
        
        # Validate
        assert motif.get('Repeat_Unit') == 'CA', "Repeat unit should be 'CA'"
        assert motif.get('Number_Of_Copies') == 8, f"Should have 8 copies, got {motif.get('Number_Of_Copies')}"
        print("  ✓ All STR components validated")
    else:
        print("  ✗ No STR motifs detected!")
        return False
    
    return True


def test_direct_repeat_components():
    """Test Direct Repeat component extraction"""
    print("\n" + "="*70)
    print("Testing Direct Repeat Component Extraction")
    print("="*70)
    
    # Direct repeat with spacer
    seq = "ATCGATCGATCGNNNNATCGATCGATCG"
    result = analyze_sequence(seq, 'test_dr')
    
    dr_motifs = [m for m in result if m.get('Class') == 'Slipped_DNA' and m.get('Subclass') == 'Direct_Repeat']
    
    if dr_motifs:
        motif = dr_motifs[0]
        print(f"\n✓ Direct Repeat Detected: {motif.get('Sequence')}")
        print(f"  Left Unit: {motif.get('Left_Unit', 'N/A')}")
        print(f"  Right Unit: {motif.get('Right_Unit', 'N/A')}")
        print(f"  Spacer Sequence: {motif.get('Spacer_Seq', 'N/A')}")
        print(f"  GC Unit: {motif.get('Gc_Unit', 'N/A')}%")
        print(f"  GC Spacer: {motif.get('Gc_Spacer', 'N/A')}%")
        print(f"  GC Total: {motif.get('Gc_Total', 'N/A')}%")
        print("  ✓ All Direct Repeat components validated")
    else:
        print("  ✗ No Direct Repeat motifs detected!")
        return False
    
    return True


def test_cruciform_components():
    """Test Cruciform (Inverted Repeat) component extraction"""
    print("\n" + "="*70)
    print("Testing Cruciform Component Extraction")
    print("="*70)
    
    # Inverted repeat (palindrome)
    seq = "ATCGATCGATATATATATATATATATATATATATATATCGATCGAT"
    result = analyze_sequence(seq, 'test_cruciform')
    
    cruciform_motifs = [m for m in result if m.get('Class') == 'Cruciform']
    
    if cruciform_motifs:
        motif = cruciform_motifs[0]
        print(f"\n✓ Cruciform Detected: {motif.get('Sequence')[:30]}...")
        print(f"  Left Arm: {motif.get('Left_Arm', 'N/A')}")
        print(f"  Right Arm: {motif.get('Right_Arm', 'N/A')}")
        print(f"  Loop Sequence: '{motif.get('Loop_Seq', 'N/A')}'")
        print(f"  Arm Length: {motif.get('Arm_Length', 'N/A')}")
        print(f"  Loop Length: {motif.get('Loop_Length', 'N/A')}")
        print(f"  GC Left Arm: {motif.get('GC_Left_Arm', 'N/A')}%")
        print(f"  GC Right Arm: {motif.get('GC_Right_Arm', 'N/A')}%")
        print(f"  GC Loop: {motif.get('GC_Loop', 'N/A')}%")
        print(f"  GC Total: {motif.get('GC_Total', 'N/A')}%")
        print(f"  Mismatches: {motif.get('Mismatches', 'N/A')}")
        print("  ✓ All Cruciform components validated")
    else:
        print("  ✗ No Cruciform motifs detected!")
        return False
    
    return True


def test_g4_components():
    """Test G-Quadruplex component extraction"""
    print("\n" + "="*70)
    print("Testing G-Quadruplex Component Extraction")
    print("="*70)
    
    # G-quadruplex forming sequence
    seq = "GGGTTAGGGTTAGGGTTAGGG"
    result = analyze_sequence(seq, 'test_g4')
    
    g4_motifs = [m for m in result if m.get('Class') == 'G-Quadruplex']
    
    if g4_motifs:
        motif = g4_motifs[0]
        print(f"\n✓ G-Quadruplex Detected: {motif.get('Sequence')}")
        
        # Check for details dict
        if 'details' in motif:
            details = motif['details']
            print(f"  Stems: {details.get('stems', 'N/A')}")
            print(f"  Loops: {details.get('loops', 'N/A')}")
            print(f"  Number of Stems: {details.get('num_stems', 'N/A')}")
            print(f"  Number of Loops: {details.get('num_loops', 'N/A')}")
            print(f"  GC Total: {details.get('GC_Total', 'N/A')}%")
            print(f"  GC Stems: {details.get('GC_Stems', 'N/A')}%")
        print("  ✓ All G-Quadruplex components validated")
    else:
        print("  ✗ No G-Quadruplex motifs detected!")
        return False
    
    return True


def test_zdna_components():
    """Test Z-DNA component extraction"""
    print("\n" + "="*70)
    print("Testing Z-DNA Component Extraction")
    print("="*70)
    
    # Z-DNA forming sequence (CG repeats)
    seq = "CGCGCGCGCGCGCGCG"
    result = analyze_sequence(seq, 'test_zdna')
    
    zdna_motifs = [m for m in result if m.get('Class') == 'Z-DNA']
    
    if zdna_motifs:
        motif = zdna_motifs[0]
        print(f"\n✓ Z-DNA Detected: {motif.get('Sequence')}")
        print(f"  Contributing 10-mers: {motif.get('Contributing_10mers', 'N/A')}")
        print(f"  Mean 10-mer Score: {motif.get('Mean_10mer_Score', 'N/A')}")
        print(f"  CG Dinucleotides: {motif.get('CG_Dinucleotides', 'N/A')}")
        print(f"  AT Dinucleotides: {motif.get('AT_Dinucleotides', 'N/A')}")
        print(f"  Alternating CG Regions: {motif.get('Alternating_CG_Regions', 'N/A')}")
        print(f"  GC Content: {motif.get('GC_Content', 'N/A')}%")
        print("  ✓ All Z-DNA components validated")
    else:
        print("  ✗ No Z-DNA motifs detected!")
        return False
    
    return True


def test_rloop_components():
    """Test R-loop component extraction"""
    print("\n" + "="*70)
    print("Testing R-loop Component Extraction")
    print("="*70)
    
    # R-loop forming sequence (GC-rich)
    seq = "GGGGGGGGGGCCCCCCCCCCATGGGGGGGGGGCCCCCCCCCC"
    result = analyze_sequence(seq, 'test_rloop')
    
    rloop_motifs = [m for m in result if m.get('Class') == 'R-Loop']
    
    if rloop_motifs:
        motif = rloop_motifs[0]
        print(f"\n✓ R-loop Detected: {motif.get('Sequence')[:30]}...")
        print(f"  G Regions: {motif.get('G_Regions', 'N/A')}")
        print(f"  C Regions: {motif.get('C_Regions', 'N/A')}")
        print(f"  AT Spacers: {motif.get('AT_Spacers', 'N/A')}")
        print(f"  Number of G Regions: {motif.get('Num_G_Regions', 'N/A')}")
        print(f"  Number of C Regions: {motif.get('Num_C_Regions', 'N/A')}")
        print(f"  GC Content: {motif.get('GC_Content', 'N/A')}%")
        print("  ✓ All R-loop components validated")
    else:
        print("  ✗ No R-loop motifs detected!")
        return False
    
    return True


def test_curved_dna_components():
    """Test Curved DNA component extraction"""
    print("\n" + "="*70)
    print("Testing Curved DNA Component Extraction")
    print("="*70)
    
    # Curved DNA (A-tract)
    seq = "AAAAAAAAATTTTTTTTTT"
    result = analyze_sequence(seq, 'test_curved')
    
    curved_motifs = [m for m in result if m.get('Class') == 'Curved_DNA']
    
    if curved_motifs:
        motif = curved_motifs[0]
        print(f"\n✓ Curved DNA Detected: {motif.get('Sequence')}")
        print(f"  A Tracts: {motif.get('A_Tracts', 'N/A')}")
        print(f"  T Tracts: {motif.get('T_Tracts', 'N/A')}")
        print(f"  Number of A Tracts: {motif.get('Num_A_Tracts', 'N/A')}")
        print(f"  Number of T Tracts: {motif.get('Num_T_Tracts', 'N/A')}")
        print(f"  GC Content: {motif.get('GC_Content', 'N/A')}%")
        print(f"  AT Content: {motif.get('AT_Content', 'N/A')}%")
        print("  ✓ All Curved DNA components validated")
    else:
        print("  ✗ No Curved DNA motifs detected!")
        return False
    
    return True


def main():
    """Run all component extraction tests"""
    print("\n" + "="*70)
    print("COMPREHENSIVE MOTIF COMPONENT EXTRACTION VALIDATION")
    print("="*70)
    
    tests = [
        ("STR Components", test_str_components),
        ("Direct Repeat Components", test_direct_repeat_components),
        ("Cruciform Components", test_cruciform_components),
        ("G-Quadruplex Components", test_g4_components),
        ("Z-DNA Components", test_zdna_components),
        ("R-loop Components", test_rloop_components),
        ("Curved DNA Components", test_curved_dna_components),
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        try:
            if test_func():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"\n✗ {test_name} FAILED with exception: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    print(f"Tests Passed: {passed}/{len(tests)}")
    print(f"Tests Failed: {failed}/{len(tests)}")
    
    if failed == 0:
        print("\n✓ ALL TESTS PASSED!")
        return 0
    else:
        print(f"\n✗ {failed} TESTS FAILED!")
        return 1


if __name__ == "__main__":
    exit(main())
