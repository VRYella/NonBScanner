#!/usr/bin/env python3
"""
Test suite for scanner_agent.py and coverage/density/enrichment calculations.

COORDINATE SYSTEM:
    NonBScanner uses 1-based INCLUSIVE coordinates:
    - Start: 1-based position (first base = 1)
    - End: 1-based position (INCLUSIVE - last base in motif)
    - Length: End - Start + 1
    
    Example: Start=1, End=20, Length=20 represents 20 bases
             In Python: range(0, 20) covers these positions

This test validates:
1. ParallelScanner chunking and deduplication
2. Coverage calculations with overlaps
3. Density calculations with overlaps
4. Enrichment calculations
5. Integration with existing scanner
"""

import sys
import time
import numpy as np
from typing import List, Dict, Any

# Test imports
try:
    from scanner_agent import ParallelScanner, scan_genome_parallel
    from utilities import (
        calculate_genomic_density, 
        calculate_positional_density,
        calculate_motif_statistics
    )
    from nonbscanner import analyze_sequence
    print("✅ All imports successful")
except ImportError as e:
    print(f"❌ Import failed: {e}")
    sys.exit(1)


def test_parallel_scanner():
    """Test ParallelScanner chunking and deduplication"""
    print("\n" + "="*70)
    print("TEST 1: ParallelScanner Chunking and Deduplication")
    print("="*70)
    
    # Create test genome
    test_genome = "GGGTTAGGGTTAGGGTTAGGG" * 50 + "AAAAATTTTAAAAATTTT" * 25
    genome_length = len(test_genome)
    
    print(f"Genome length: {genome_length:,} bp")
    
    # Create scanner
    scanner = ParallelScanner(test_genome, hs_db=None)
    stats = scanner.get_statistics()
    
    print(f"Configuration:")
    print(f"  Chunk size: {stats['chunk_size']:,} bp")
    print(f"  Overlap size: {stats['overlap_size']:,} bp")
    print(f"  Number of chunks: {stats['num_chunks']}")
    print(f"  Workers: {stats['num_workers']}")
    print(f"  Using Hyperscan: {stats['using_hyperscan']}")
    
    # Run scan
    start_time = time.time()
    motifs = scanner.run_scan()
    scan_time = time.time() - start_time
    
    print(f"\nResults:")
    print(f"  Motifs found: {len(motifs)}")
    print(f"  Scan time: {scan_time:.3f}s")
    print(f"  Speed: {genome_length/scan_time:,.0f} bp/s")
    
    # Verify motifs are tuples
    if motifs:
        sample = motifs[0]
        print(f"  Sample motif format: {type(sample)} with {len(sample)} elements")
        assert len(sample) == 3, "Motif should be (start, end, pattern_id)"
        print("  ✅ Motif format correct")
    
    print("\n✅ TEST 1 PASSED")
    return True


def test_coverage_with_overlaps():
    """Test coverage calculations handle overlaps correctly"""
    print("\n" + "="*70)
    print("TEST 2: Coverage Calculation with Overlaps")
    print("="*70)
    
    # Create test cases
    test_cases = [
        {
            'name': 'No overlaps',
            'motifs': [
                {'Class': 'G4', 'Start': 1, 'End': 20, 'Length': 20},
                {'Class': 'G4', 'Start': 31, 'End': 50, 'Length': 20},  # Non-overlapping (31-50 = 20 positions in 1-based inclusive)
            ],
            'seq_len': 100,
            'expected_coverage': 40.0  # 20 + 20 = 40
        },
        {
            'name': 'Full overlap (identical)',
            'motifs': [
                {'Class': 'G4', 'Start': 1, 'End': 20, 'Length': 20},
                {'Class': 'G4', 'Start': 1, 'End': 20, 'Length': 20},
            ],
            'seq_len': 100,
            'expected_coverage': 20.0
        },
        {
            'name': 'Partial overlap',
            'motifs': [
                {'Class': 'G4', 'Start': 1, 'End': 30, 'Length': 30},
                {'Class': 'G4', 'Start': 20, 'End': 50, 'Length': 30},  # 10bp overlap
            ],
            'seq_len': 100,
            'expected_coverage': 50.0  # 1-50 = 50bp
        },
        {
            'name': 'Multiple overlapping classes',
            'motifs': [
                {'Class': 'G4', 'Start': 1, 'End': 40, 'Length': 40},
                {'Class': 'Curved', 'Start': 30, 'End': 60, 'Length': 30},
                {'Class': 'Z_DNA', 'Start': 50, 'End': 80, 'Length': 30},
            ],
            'seq_len': 100,
            'expected_coverage': 80.0  # 1-80 = 80bp
        },
    ]
    
    all_passed = True
    for test_case in test_cases:
        print(f"\nTest: {test_case['name']}")
        print(f"  Motifs: {len(test_case['motifs'])}")
        print(f"  Sum of lengths: {sum(m['Length'] for m in test_case['motifs'])} bp")
        
        density = calculate_genomic_density(
            test_case['motifs'], 
            test_case['seq_len'], 
            by_class=False
        )
        
        coverage = density['Overall']
        expected = test_case['expected_coverage']
        
        print(f"  Expected coverage: {expected}%")
        print(f"  Calculated coverage: {coverage}%")
        
        # Check if within tolerance
        tolerance = 0.01
        passed = abs(coverage - expected) < tolerance
        
        if passed:
            print(f"  ✅ PASSED")
        else:
            print(f"  ❌ FAILED (difference: {abs(coverage - expected):.4f}%)")
            all_passed = False
        
        # Ensure coverage never exceeds 100%
        assert coverage <= 100.0, f"Coverage {coverage}% exceeds 100%!"
    
    if all_passed:
        print("\n✅ TEST 2 PASSED")
    else:
        print("\n❌ TEST 2 FAILED")
    
    return all_passed


def test_density_calculations():
    """Test positional density calculations"""
    print("\n" + "="*70)
    print("TEST 3: Positional Density Calculations")
    print("="*70)
    
    # Test motifs
    motifs = [
        {'Class': 'G4', 'Subclass': 'Canonical', 'Start': 1, 'End': 20, 'Length': 20},
        {'Class': 'G4', 'Subclass': 'Relaxed', 'Start': 30, 'End': 50, 'Length': 20},
        {'Class': 'Curved', 'Subclass': 'A-tract', 'Start': 60, 'End': 80, 'Length': 20},
    ]
    
    sequence_length = 1000  # 1kb
    
    print(f"Sequence length: {sequence_length} bp (1 kb)")
    print(f"Number of motifs: {len(motifs)}")
    
    # Test kbp density
    density_kbp = calculate_positional_density(motifs, sequence_length, unit='kbp', by_class=False)
    print(f"\nDensity (motifs/kb): {density_kbp['Overall']}")
    assert density_kbp['Overall'] == 3.0, "Should be 3 motifs per kb"
    print("  ✅ kbp density correct")
    
    # Test Mbp density
    sequence_length_mb = 1000000  # 1Mb
    density_mbp = calculate_positional_density(motifs, sequence_length_mb, unit='Mbp', by_class=False)
    print(f"Density (motifs/Mb): {density_mbp['Overall']}")
    assert density_mbp['Overall'] == 3.0, "Should be 3 motifs per Mb"
    print("  ✅ Mbp density correct")
    
    # Test by-class density
    density_by_class = calculate_positional_density(motifs, sequence_length, unit='kbp', by_class=True)
    print(f"\nDensity by class:")
    for cls, dens in sorted(density_by_class.items()):
        print(f"  {cls}: {dens} motifs/kb")
    
    assert density_by_class['G4'] == 2.0, "Should be 2 G4 motifs per kb"
    assert density_by_class['Curved'] == 1.0, "Should be 1 Curved motif per kb"
    print("  ✅ By-class density correct")
    
    print("\n✅ TEST 3 PASSED")
    return True


def test_motif_statistics():
    """Test comprehensive motif statistics"""
    print("\n" + "="*70)
    print("TEST 4: Motif Statistics Integration")
    print("="*70)
    
    # Test with overlapping motifs
    # Note: Start and End are 1-based inclusive coordinates
    motifs = [
        {'Class': 'G4', 'Subclass': 'Canonical', 'Start': 1, 'End': 50, 'Length': 50, 'Score': 0.8},
        {'Class': 'G4', 'Subclass': 'Relaxed', 'Start': 40, 'End': 80, 'Length': 41, 'Score': 0.7},  # Overlaps (positions 40-80 inclusive = 41 bases)
        {'Class': 'Curved', 'Subclass': 'A-tract', 'Start': 100, 'End': 130, 'Length': 31, 'Score': 0.9},  # 100-130 inclusive = 31 bases
    ]
    
    sequence_length = 200
    
    print(f"Sequence length: {sequence_length} bp")
    print(f"Motifs: {len(motifs)}")
    print(f"Sum of lengths: {sum(m['Length'] for m in motifs)} bp = {sum(m['Length'] for m in motifs)/ sequence_length * 100:.1f}%")
    
    stats = calculate_motif_statistics(motifs, sequence_length)
    
    print(f"\nStatistics:")
    print(f"  Total motifs: {stats['Total_Motifs']}")
    print(f"  Coverage: {stats['Coverage%']}%")
    print(f"  Density: {stats['Density']} motifs/kb")
    print(f"  Classes detected: {stats['Classes_Detected']}")
    print(f"  Subclasses detected: {stats['Subclasses_Detected']}")
    
    # Verify coverage is correct (should be less than sum of lengths due to overlap)
    # In 0-based half-open: [0,50) + [39,80) + [99,130)
    # Unique positions: 0-49 (50), 39-79 (overlap 39-49 = 10), 99-129 (31)
    # Total unique: 50 + (80-39-10) + 31 = 50 + 31 + 31 = 112? No wait...
    # Simpler: range(0,50) ∪ range(39,80) ∪ range(99,130)
    # = [0,79] ∪ [99,129] = 80 + 31 = 111 positions
    expected_coverage = 111 / 200 * 100  # 55.5%
    assert abs(stats['Coverage%'] - expected_coverage) < 0.1, f"Coverage should be ~{expected_coverage}%"
    print(f"  ✅ Coverage correct: {stats['Coverage%']}% (expected ~{expected_coverage}%)")
    
    # Ensure coverage never exceeds 100%
    assert stats['Coverage%'] <= 100.0, "Coverage exceeds 100%!"
    print("  ✅ Coverage capped at 100%")
    
    print("\n✅ TEST 4 PASSED")
    return True


def test_integration_with_scanner():
    """Test integration with existing analyze_sequence"""
    print("\n" + "="*70)
    print("TEST 5: Integration with Existing Scanner")
    print("="*70)
    
    # Test sequence with known motifs
    test_seq = "GGGTTAGGGTTAGGGTTAGGG" * 5 + "AAAAATTTTAAAAATTTT" * 3
    
    print(f"Test sequence length: {len(test_seq)} bp")
    
    # Run analysis
    start_time = time.time()
    motifs = analyze_sequence(test_seq, "integration_test")
    analysis_time = time.time() - start_time
    
    print(f"Analysis time: {analysis_time:.3f}s")
    print(f"Motifs found: {len(motifs)}")
    
    if motifs:
        # Calculate statistics
        stats = calculate_motif_statistics(motifs, len(test_seq))
        
        print(f"\nStatistics:")
        print(f"  Coverage: {stats['Coverage%']}%")
        print(f"  Density: {stats['Density']} motifs/kb")
        print(f"  Classes: {stats['Classes_Detected']}")
        
        # Verify coverage is reasonable
        assert 0 <= stats['Coverage%'] <= 100, "Coverage out of bounds!"
        print("  ✅ Coverage within bounds")
        
        # Check motif structure
        sample = motifs[0]
        required_fields = ['Class', 'Subclass', 'Start', 'End', 'Length']
        for field in required_fields:
            assert field in sample, f"Missing field: {field}"
        print("  ✅ Motif structure correct")
    
    print("\n✅ TEST 5 PASSED")
    return True


def run_all_tests():
    """Run all tests"""
    print("\n" + "="*70)
    print("SCANNER AGENT & COVERAGE CALCULATION TEST SUITE")
    print("="*70)
    
    tests = [
        ("ParallelScanner", test_parallel_scanner),
        ("Coverage with Overlaps", test_coverage_with_overlaps),
        ("Density Calculations", test_density_calculations),
        ("Motif Statistics", test_motif_statistics),
        ("Scanner Integration", test_integration_with_scanner),
    ]
    
    results = []
    for name, test_func in tests:
        try:
            result = test_func()
            results.append((name, result))
        except Exception as e:
            print(f"\n❌ TEST FAILED WITH EXCEPTION: {e}")
            import traceback
            traceback.print_exc()
            results.append((name, False))
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for name, result in results:
        status = "✅ PASSED" if result else "❌ FAILED"
        print(f"{name:.<50} {status}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 ALL TESTS PASSED! 🎉")
        return 0
    else:
        print("\n⚠️  SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    exit_code = run_all_tests()
    sys.exit(exit_code)
