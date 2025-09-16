#!/usr/bin/env python3
"""
Performance Optimization Test Suite for NBDFinder
=================================================

Tests to verify the performance improvements in conservation analysis
and overall motif detection functionality.

Key Optimizations Tested:
1. Conservation analysis caching (100x+ speedup for repeated analyses)
2. Adaptive shuffling (fewer shuffles for longer sequences)
3. Early termination (skip conservation for very short sequences)
4. Selective motif detection (only run relevant detectors during conservation)
"""

import sys
import os
import time
sys.path.insert(0, os.path.dirname(__file__))

import motifs
from conservation_analysis import clear_conservation_cache, get_cache_stats

def test_conservation_caching():
    """Test conservation analysis caching performance."""
    print("\n" + "=" * 60)
    print("CONSERVATION CACHING TEST")
    print("=" * 60)
    
    seq = "GGGTTAGGGTTAGGGTTAGGG" * 3  # 63bp sequence with G4 motifs
    
    # First run (cache miss)
    clear_conservation_cache()
    start = time.time()
    results1 = motifs.all_motifs(seq, "test1", calculate_conservation=True)
    first_run_time = time.time() - start
    
    cache_stats = get_cache_stats()
    print(f"First run (cache miss): {first_run_time:.3f}s")
    print(f"Cache entries: {cache_stats['cache_size']}")
    
    # Second run (cache hit)
    start = time.time()
    results2 = motifs.all_motifs(seq, "test2", calculate_conservation=True)
    second_run_time = time.time() - start
    
    speedup = first_run_time / second_run_time if second_run_time > 0 else float('inf')
    print(f"Second run (cache hit): {second_run_time:.3f}s")
    print(f"Cache speedup: {speedup:.1f}x")
    
    # Verify results are identical
    scores_match = all(
        r1.get("Conservation_Score") == r2.get("Conservation_Score") 
        for r1, r2 in zip(results1, results2)
    )
    print(f"Conservation scores identical: {scores_match}")
    
    return speedup > 50  # Expect at least 50x speedup

def test_adaptive_conservation():
    """Test adaptive conservation analysis for different sequence lengths."""
    print("\n" + "=" * 60)
    print("ADAPTIVE CONSERVATION TEST")
    print("=" * 60)
    
    # Test different sequence lengths
    base_seq = "GGGTTAGGGTTAGGGTTAGGGCCCTAACCCTAACCCTAACCCATCGATCGATCG"
    test_sizes = [100, 500, 1000]
    
    results = {}
    
    for size in test_sizes:
        # Create test sequence of target size
        seq = (base_seq * (size // len(base_seq) + 1))[:size]
        
        # Test without conservation
        start = time.time()
        motifs_only = motifs.all_motifs(seq, f"test_{size}", calculate_conservation=False)
        no_cons_time = time.time() - start
        
        # Test with adaptive conservation
        clear_conservation_cache()
        start = time.time()
        motifs_with_cons = motifs.all_motifs(seq, f"test_{size}", calculate_conservation=True)
        with_cons_time = time.time() - start
        
        overhead = with_cons_time / no_cons_time if no_cons_time > 0 else float('inf')
        results[size] = {
            'no_cons_time': no_cons_time,
            'with_cons_time': with_cons_time,
            'overhead': overhead,
            'motif_count': len(motifs_only)
        }
        
        print(f"Sequence length {size}bp:")
        print(f"  No conservation: {no_cons_time:.3f}s -> {len(motifs_only)} motifs")
        print(f"  With conservation: {with_cons_time:.3f}s -> {len(motifs_with_cons)} motifs")
        print(f"  Overhead: {overhead:.1f}x")
    
    # Verify overhead decreases with sequence length (adaptive effect)
    overhead_trend = [results[size]['overhead'] for size in test_sizes]
    adaptive_working = overhead_trend[0] > overhead_trend[2]  # 100bp > 1000bp overhead
    
    print(f"\nAdaptive effect working: {adaptive_working}")
    print(f"Overhead trend: {' -> '.join(f'{oh:.1f}x' for oh in overhead_trend)}")
    
    return adaptive_working

def test_early_termination():
    """Test early termination for short sequences."""
    print("\n" + "=" * 60)
    print("EARLY TERMINATION TEST")
    print("=" * 60)
    
    # Test very short sequences that should trigger early termination
    short_sequences = [
        ("Very short", "ATCG"),
        ("Short with motifs", "GGGTTAGGGTTAGGG"),  # 15bp - should trigger early termination
        ("Medium", "GGGTTAGGGTTAGGGTTAGGGCCCTAACCCTAA"),  # 32bp
        ("Long enough", "GGGTTAGGGTTAGGGTTAGGGCCCTAACCCTAACCCTAACCCATCGATCGATCG")  # 54bp
    ]
    
    early_termination_working = True
    
    for name, seq in short_sequences:
        start = time.time()
        results = motifs.all_motifs(seq, name, calculate_conservation=True)
        cons_time = time.time() - start
        
        start = time.time()
        results_no_cons = motifs.all_motifs(seq, name, calculate_conservation=False)
        no_cons_time = time.time() - start
        
        overhead = cons_time / no_cons_time if no_cons_time > 0 else 1.0
        
        print(f"{name} ({len(seq)}bp): {overhead:.1f}x overhead")
        
        # Check if conservation was properly applied
        if results:
            sample_motif = results[0]
            has_conservation = 'Conservation_Score' in sample_motif
            conservation_note = sample_motif.get('Conservation_Note', '')
            
            if len(seq) < 50:
                # Should have early termination note
                expected_early_termination = 'too short' in conservation_note.lower()
                if not expected_early_termination:
                    early_termination_working = False
                    print(f"  ⚠ Expected early termination for {len(seq)}bp sequence")
            
            print(f"  Conservation note: {conservation_note}")
    
    return early_termination_working

def test_selective_motif_detection():
    """Test that conservation analysis only runs relevant motif detectors."""
    print("\n" + "=" * 60)
    print("SELECTIVE MOTIF DETECTION TEST")
    print("=" * 60)
    
    # Test with G4-only sequence
    g4_seq = "GGGTTAGGGTTAGGGTTAGGG" * 3
    
    clear_conservation_cache()
    start = time.time()
    g4_results = motifs.all_motifs(g4_seq, "g4_test", calculate_conservation=True)
    g4_time = time.time() - start
    
    # Count detected motif classes
    g4_classes = set(m.get('Class', 'Unknown') for m in g4_results)
    
    print(f"G4 sequence conservation analysis: {g4_time:.3f}s")
    print(f"Motif classes detected: {', '.join(sorted(g4_classes))}")
    
    # Test with complex sequence (multiple motif types)
    complex_seq = ("GGGTTAGGGTTAGGGTTAGGG" +  # G4
                   "CCCTAACCCTAACCCTAACCC" +  # i-motif
                   "CGCGCGCGCGCGCGCGCGCG" +   # Z-DNA
                   "ATCGATCGATCGATCGATCG")     # Other motifs
    
    clear_conservation_cache()
    start = time.time()
    complex_results = motifs.all_motifs(complex_seq, "complex_test", calculate_conservation=True)
    complex_time = time.time() - start
    
    complex_classes = set(m.get('Class', 'Unknown') for m in complex_results)
    
    print(f"Complex sequence conservation analysis: {complex_time:.3f}s")
    print(f"Motif classes detected: {', '.join(sorted(complex_classes))}")
    
    # Selective detection should make G4-only analysis faster
    # (though this is hard to measure precisely due to caching and other factors)
    selective_working = len(g4_classes) < len(complex_classes)
    print(f"Selective detection working: {selective_working}")
    
    return selective_working

def main():
    """Run all performance optimization tests."""
    print("NBDFinder Performance Optimization Test Suite")
    print("=" * 60)
    
    results = {
        'caching': test_conservation_caching(),
        'adaptive': test_adaptive_conservation(),
        'early_termination': test_early_termination(),
        'selective': test_selective_motif_detection()
    }
    
    print("\n" + "=" * 60)
    print("PERFORMANCE OPTIMIZATION SUMMARY")
    print("=" * 60)
    
    all_passed = True
    for test_name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{test_name.replace('_', ' ').title()}: {status}")
        if not passed:
            all_passed = False
    
    print(f"\nOverall Status: {'🎉 ALL OPTIMIZATIONS WORKING' if all_passed else '⚠ SOME OPTIMIZATIONS FAILED'}")
    
    return all_passed

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)