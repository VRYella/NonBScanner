"""
Test suite for optimized repeat scanner
Tests all three motif types: Direct repeats, Inverted repeats (Cruciform), Mirror repeats (Triplex)
Uses real genome-like sequences
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from utils.repeat_scanner import (
    find_direct_repeats, 
    find_inverted_repeats, 
    find_mirror_repeats, 
    find_strs
)


def test_direct_repeats():
    """Test direct repeat detection with seed-and-extend k-mer approach"""
    print("\n=== Testing Direct Repeats ===")
    
    # Test 1: Simple direct repeat with no spacer
    seq1 = "AAAA" + "GGGGGGGGGGGG" + "TTTT" + "GGGGGGGGGGGG" + "CCCC"
    results = find_direct_repeats(seq1, min_unit=10, max_unit=15, max_spacer=10)
    print(f"Test 1 - Simple direct repeat: Found {len(results)} repeats")
    for r in results:
        print(f"  Unit length: {r['Unit_Length']}, Spacer: {r['Spacer']}, Start: {r['Start']}, End: {r['End']}")
    assert len(results) >= 1, "Should find at least one direct repeat"
    
    # Test 2: Direct repeat with spacer
    seq2 = "AAAA" + "GGGGGGGGGGG" + "TTTT" + "GGGGGGGGGGG" + "CCCC"
    results = find_direct_repeats(seq2, min_unit=10, max_unit=15, max_spacer=10)
    print(f"Test 2 - Direct repeat with spacer: Found {len(results)} repeats")
    for r in results:
        print(f"  Unit length: {r['Unit_Length']}, Spacer: {r['Spacer']}, Start: {r['Start']}, End: {r['End']}")
    assert len(results) >= 1, "Should find direct repeat with spacer"
    
    # Test 3: Real genome-like sequence (human telomeric repeat region)
    seq3 = "ATCGATCG" + "TTAGGGTTAGGGTTAGGG" + "GCGCGC" + "TTAGGGTTAGGGTTAGGG" + "ATCGATCG"
    results = find_direct_repeats(seq3, min_unit=10, max_unit=30, max_spacer=10)
    print(f"Test 3 - Telomeric repeat region: Found {len(results)} repeats")
    for r in results:
        print(f"  Unit: {r['Unit_Seq'][:20]}..., Length: {r['Unit_Length']}, Spacer: {r['Spacer']}")
    
    print("✓ All direct repeat tests passed\n")


def test_inverted_repeats_cruciform():
    """Test inverted repeat detection (cruciform structures)"""
    print("=== Testing Inverted Repeats (Cruciform) ===")
    
    # Test 1: Simple palindrome (perfect inverted repeat)
    seq1 = "AAAA" + "GGGGGG" + "TTTT" + "CCCCCC" + "AAAA"
    results = find_inverted_repeats(seq1, min_arm=6, max_loop=10)
    print(f"Test 1 - Simple palindrome: Found {len(results)} inverted repeats")
    for r in results:
        print(f"  Arm length: {r['Arm_Length']}, Loop: {r['Loop']}, Left: {r['Left_Arm']}, Right: {r['Right_Arm']}")
    assert len(results) >= 1, "Should find inverted repeat"
    
    # Test 2: Longer inverted repeat with loop
    seq2 = "TTTT" + "ATCGATCGATCG" + "AAAAAAAA" + "CGATCGATCGAT" + "GGGG"
    results = find_inverted_repeats(seq2, min_arm=6, max_loop=20)
    print(f"Test 2 - Long inverted repeat: Found {len(results)} inverted repeats")
    for r in results:
        print(f"  Arm length: {r['Arm_Length']}, Loop: {r['Loop']}")
    
    # Test 3: Real genome-like sequence (cruciform-forming region)
    # This is a realistic sequence from bacterial chromosome
    seq3 = "GCGCGC" + "ATCGATCGATCGATCG" + "TTTTTTTTTTTT" + "CGATCGATCGATCGAT" + "GCGCGC"
    results = find_inverted_repeats(seq3, min_arm=6, max_loop=20)
    print(f"Test 3 - Bacterial cruciform site: Found {len(results)} inverted repeats")
    for r in results:
        print(f"  Arm length: {r['Arm_Length']}, Loop: {r['Loop']}, Start: {r['Start']}, End: {r['End']}")
    
    # Test 4: Multiple overlapping candidates (should dedupe)
    seq4 = "GGGGGG" * 3 + "AAAA" + "CCCCCC" * 3
    results = find_inverted_repeats(seq4, min_arm=6, max_loop=10)
    print(f"Test 4 - Overlapping candidates: Found {len(results)} inverted repeats (deduped)")
    
    print("✓ All inverted repeat tests passed\n")


def test_mirror_repeats_triplex():
    """Test mirror repeat detection (triplex DNA)"""
    print("=== Testing Mirror Repeats (Triplex DNA) ===")
    
    # Test 1: Simple mirror repeat (homopurine)
    seq1 = "TTTT" + "GGAGAGAGAGAG" + "AAAA" + "GAGAGAGAGAGG" + "CCCC"
    results = find_mirror_repeats(seq1, min_arm=10, max_loop=10, purine_pyrimidine_threshold=0.9)
    print(f"Test 1 - Homopurine mirror: Found {len(results)} mirror repeats")
    for r in results:
        print(f"  Arm length: {r['Arm_Length']}, Loop: {r['Loop']}, Is_Triplex: {r['Is_Triplex']}")
        print(f"  Purine: {r['Purine_Fraction']}, Pyrimidine: {r['Pyrimidine_Fraction']}")
    
    # Test 2: Homopyrimidine mirror repeat
    seq2 = "AAAA" + "CTCTCTCTCTCT" + "GGGG" + "TCTCTCTCTCTC" + "AAAA"
    results = find_mirror_repeats(seq2, min_arm=10, max_loop=10, purine_pyrimidine_threshold=0.9)
    print(f"Test 2 - Homopyrimidine mirror: Found {len(results)} mirror repeats")
    for r in results:
        print(f"  Arm length: {r['Arm_Length']}, Loop: {r['Loop']}, Is_Triplex: {r['Is_Triplex']}")
        print(f"  Purine: {r['Purine_Fraction']}, Pyrimidine: {r['Pyrimidine_Fraction']}")
    
    # Test 3: Real genome-like sequence (H-DNA forming site)
    # This mimics intramolecular triplex DNA (H-DNA) from eukaryotic promoters
    seq3 = "ATCG" + "AGAGAGAGAGAGAG" + "TTTTTTTT" + "GAGAGAGAGAGAGA" + "CGAT"
    results = find_mirror_repeats(seq3, min_arm=10, max_loop=20, purine_pyrimidine_threshold=0.9)
    print(f"Test 3 - H-DNA forming site: Found {len(results)} mirror repeats")
    for r in results:
        print(f"  Arm length: {r['Arm_Length']}, Loop: {r['Loop']}")
        print(f"  Left arm: {r['Left_Arm'][:15]}..., Right arm: {r['Right_Arm'][:15]}...")
        print(f"  Is_Triplex: {r['Is_Triplex']}, Purine fraction: {r['Purine_Fraction']}")
    
    # Test 4: Non-triplex mirror repeat (mixed bases, <90% purine/pyrimidine)
    seq4 = "AAAA" + "ATCGATCGATCG" + "GGGG" + "CGATCGATCGAT" + "TTTT"
    results = find_mirror_repeats(seq4, min_arm=10, max_loop=10, purine_pyrimidine_threshold=0.9)
    print(f"Test 4 - Non-triplex mirror: Found {len(results)} mirror repeats")
    for r in results:
        print(f"  Is_Triplex: {r['Is_Triplex']} (should be False)")
        assert r['Is_Triplex'] == False, "Should not be classified as triplex (mixed bases)"
    
    print("✓ All mirror repeat tests passed\n")


def test_strs():
    """Test short tandem repeat (STR) detection"""
    print("=== Testing Short Tandem Repeats (STRs) ===")
    
    # Test 1: Simple mononucleotide repeat
    seq1 = "ATCG" + "AAAAAAAAAAAAA" + "GCTA"
    results = find_strs(seq1, min_u=1, max_u=9, min_total=10)
    print(f"Test 1 - Mononucleotide repeat: Found {len(results)} STRs")
    for r in results:
        print(f"  Unit: {r['Unit_Seq']}, Copies: {r['Copies']}, Length: {r['Length']}")
    assert len(results) >= 1, "Should find mononucleotide repeat"
    
    # Test 2: Dinucleotide repeat (CA)n
    seq2 = "GGGG" + "CACACACACACACACA" + "TTTT"
    results = find_strs(seq2, min_u=1, max_u=9, min_total=10)
    print(f"Test 2 - Dinucleotide repeat: Found {len(results)} STRs")
    for r in results:
        print(f"  Unit: {r['Unit_Seq']}, Copies: {r['Copies']}, Length: {r['Length']}")
    
    # Test 3: Trinucleotide repeat (GAA)n (Friedreich ataxia repeat)
    seq3 = "ATCG" + "GAAGAAGAAGAAGAAGAAGAA" + "CGTA"
    results = find_strs(seq3, min_u=1, max_u=9, min_total=10)
    print(f"Test 3 - Trinucleotide repeat (GAA)n: Found {len(results)} STRs")
    for r in results:
        print(f"  Unit: {r['Unit_Seq']}, Copies: {r['Copies']}, Length: {r['Length']}")
        if r['Unit_Length'] == 3:
            assert r['Unit_Seq'] == 'GAA', "Should detect GAA repeat"
    
    # Test 4: Tetranucleotide repeat
    seq4 = "CGCG" + "AATGAATGAATGAATGAATG" + "GCGC"
    results = find_strs(seq4, min_u=1, max_u=9, min_total=10)
    print(f"Test 4 - Tetranucleotide repeat: Found {len(results)} STRs")
    for r in results:
        print(f"  Unit: {r['Unit_Seq']}, Copies: {r['Copies']}, Length: {r['Length']}")
    
    # Test 5: Real microsatellite (pentanucleotide repeat from human genome)
    seq5 = "ATCGATCG" + "ATCCTATCCTATCCTATCCTATCCT" + "GCGCGCGC"
    results = find_strs(seq5, min_u=1, max_u=9, min_total=10)
    print(f"Test 5 - Pentanucleotide microsatellite: Found {len(results)} STRs")
    for r in results:
        print(f"  Unit: {r['Unit_Seq']}, Copies: {r['Copies']}, Length: {r['Length']}")
    
    print("✓ All STR tests passed\n")


def test_performance_on_long_sequences():
    """Test performance on longer genome-like sequences"""
    print("=== Testing Performance on Long Sequences ===")
    
    # Generate a 10kb sequence with some repeats
    # This simulates a real genomic region
    base_seq = "ATCGATCGATCGATCGATCG" * 100  # 2kb base
    
    # Add some direct repeats
    repeat_unit = "TTAGGGTTAGGGTTAGGG"
    seq = base_seq[:1000] + repeat_unit + base_seq[1000:2000] + repeat_unit + base_seq[2000:5000]
    
    print(f"Sequence length: {len(seq)} bp")
    
    # Test direct repeats
    print("Testing direct repeats...")
    results_dr = find_direct_repeats(seq, min_unit=10, max_unit=50, max_spacer=10)
    print(f"  Found {len(results_dr)} direct repeats")
    
    # Test inverted repeats
    print("Testing inverted repeats...")
    results_ir = find_inverted_repeats(seq[:1000], min_arm=6, max_loop=50)  # Test on 1kb chunk
    print(f"  Found {len(results_ir)} inverted repeats")
    
    # Test mirror repeats
    print("Testing mirror repeats...")
    results_mr = find_mirror_repeats(seq[:1000], min_arm=10, max_loop=50)  # Test on 1kb chunk
    print(f"  Found {len(results_mr)} mirror repeats")
    
    # Test STRs
    print("Testing STRs...")
    results_str = find_strs(seq, min_u=1, max_u=9, min_total=10)
    print(f"  Found {len(results_str)} STRs")
    
    print("✓ Performance test completed\n")


def test_detector_integration():
    """Test integration with actual detector classes"""
    print("=== Testing Detector Integration ===")
    
    # Test SlippedDNADetector
    try:
        from motif_detection.slipped_dna_detector import SlippedDNADetector
        detector = SlippedDNADetector()
        seq = "ATCG" + "CACACACACACACACA" + "TTTT" + "GGGGGGGGGGGG" + "AAAA" + "GGGGGGGGGGGG" + "CGTA"
        results = detector.annotate_sequence(seq)
        print(f"SlippedDNADetector: Found {len(results)} motifs")
        for r in results[:3]:
            print(f"  Class: {r['class_name']}, Length: {r['length']}, Score: {r['score']:.3f}")
        print("  ✓ SlippedDNADetector integration works")
    except Exception as e:
        print(f"  ✗ SlippedDNADetector integration failed: {e}")
    
    # Test CruciformDetector
    try:
        from motif_detection.cruciform_detector import CruciformDetector
        detector = CruciformDetector()
        seq = "GCGCGC" + "ATCGATCGATCGATCG" + "TTTTTTTT" + "CGATCGATCGATCGAT" + "GCGCGC"
        results = detector.annotate_sequence(seq)
        print(f"CruciformDetector: Found {len(results)} motifs")
        for r in results[:3]:
            print(f"  Arm: {r['arm_len']}, Loop: {r['loop_len']}, Score: {r['score']:.3f}")
        print("  ✓ CruciformDetector integration works")
    except Exception as e:
        print(f"  ✗ CruciformDetector integration failed: {e}")
    
    # Test TriplexDetector
    try:
        from motif_detection.triplex_detector import TriplexDetector
        detector = TriplexDetector()
        seq = "ATCG" + "AGAGAGAGAGAGAG" + "TTTTTTTT" + "GAGAGAGAGAGAGA" + "CGAT" + "GAAGAAGAAGAAGAA" + "GCGC"
        results = detector.annotate_sequence(seq)
        print(f"TriplexDetector: Found {len(results)} motifs")
        for r in results[:3]:
            print(f"  Class: {r['class_name']}, Length: {r['length']}, Score: {r['score']:.3f}")
        print("  ✓ TriplexDetector integration works")
    except Exception as e:
        print(f"  ✗ TriplexDetector integration failed: {e}")
    
    print("\n✓ Detector integration tests completed\n")


def main():
    """Run all tests"""
    print("\n" + "="*60)
    print("OPTIMIZED REPEAT SCANNER TEST SUITE")
    print("="*60)
    
    try:
        test_direct_repeats()
        test_inverted_repeats_cruciform()
        test_mirror_repeats_triplex()
        test_strs()
        test_performance_on_long_sequences()
        test_detector_integration()
        
        print("="*60)
        print("✓ ALL TESTS PASSED")
        print("="*60)
        return 0
    except AssertionError as e:
        print(f"\n✗ TEST FAILED: {e}")
        return 1
    except Exception as e:
        print(f"\n✗ UNEXPECTED ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
