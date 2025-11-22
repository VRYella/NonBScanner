#!/usr/bin/env python3
"""
Benchmark script for NonBScanner two-layer architecture.

Tests performance improvements of the new architecture:
- Single-threaded vs multi-threaded
- Standard mode vs fast mode
- Different sequence sizes
- Realistic genomic sequences
"""

import time
import sys
from typing import List, Dict, Any

try:
    import nonbscanner as nbs
    from two_layer_scanner import TwoLayerScanner
except ImportError:
    print("Error: Cannot import NonBScanner modules")
    sys.exit(1)


def generate_test_sequences() -> Dict[str, str]:
    """Generate test sequences of various sizes"""
    
    # Motif-rich unit (contains G4, i-motif, Z-DNA, A-tracts, etc.)
    unit = "GGGTTAGGGTTAGGGTTAGGGAAAAATTTTCGCGCGCGCGATATATATATCCCCTAACCCTAACCCTAACCC"
    
    sequences = {
        "Small (720bp)": unit * 10,
        "Medium (7.2kb)": unit * 100,
        "Large (72kb)": unit * 1000,
        "XLarge (720kb)": unit * 10000,
    }
    
    return sequences


def benchmark_mode(sequence: str, seq_name: str, mode: str) -> tuple:
    """
    Benchmark a specific mode.
    
    Returns:
        (num_motifs, elapsed_time, speed_bps)
    """
    start = time.time()
    
    if mode == "standard":
        motifs = nbs.analyze_sequence(sequence, seq_name, use_fast_mode=False)
    elif mode == "fast_single":
        scanner = TwoLayerScanner()
        motifs = scanner.analyze_sequence(sequence, seq_name, use_parallel=False)
    elif mode == "fast_multi":
        scanner = TwoLayerScanner()
        motifs = scanner.analyze_sequence(sequence, seq_name, use_parallel=True)
    else:
        raise ValueError(f"Unknown mode: {mode}")
    
    elapsed = time.time() - start
    speed = len(sequence) / elapsed if elapsed > 0 else 0
    
    return len(motifs), elapsed, speed


def run_benchmarks():
    """Run comprehensive benchmarks"""
    
    print("="*80)
    print("NonBScanner Two-Layer Architecture Benchmark")
    print("="*80)
    
    sequences = generate_test_sequences()
    modes = ["standard", "fast_single", "fast_multi"]
    mode_names = {
        "standard": "Standard Mode",
        "fast_single": "Fast Mode (Single-threaded)",
        "fast_multi": "Fast Mode (Multi-threaded)"
    }
    
    results = {}
    
    for seq_name, sequence in sequences.items():
        print(f"\n{seq_name}")
        print("-" * 80)
        print(f"Sequence length: {len(sequence):,} bp")
        
        results[seq_name] = {}
        
        for mode in modes:
            try:
                num_motifs, elapsed, speed = benchmark_mode(sequence, seq_name, mode)
                results[seq_name][mode] = {
                    'motifs': num_motifs,
                    'time': elapsed,
                    'speed': speed
                }
                
                print(f"\n{mode_names[mode]}:")
                print(f"  Motifs found: {num_motifs:,}")
                print(f"  Time: {elapsed:.4f}s")
                print(f"  Speed: {speed:,.0f} bp/s")
                
            except Exception as e:
                print(f"\n{mode_names[mode]}: ERROR - {e}")
                results[seq_name][mode] = None
    
    # Print summary
    print("\n" + "="*80)
    print("SUMMARY - Speedup Factors")
    print("="*80)
    
    for seq_name in sequences.keys():
        print(f"\n{seq_name}:")
        
        if results[seq_name].get("standard") and results[seq_name].get("fast_single"):
            std_time = results[seq_name]["standard"]["time"]
            fast_single_time = results[seq_name]["fast_single"]["time"]
            if fast_single_time > 0:
                speedup_single = std_time / fast_single_time
                print(f"  Fast (single) vs Standard: {speedup_single:.1f}x faster")
        
        if results[seq_name].get("standard") and results[seq_name].get("fast_multi"):
            std_time = results[seq_name]["standard"]["time"]
            fast_multi_time = results[seq_name]["fast_multi"]["time"]
            if fast_multi_time > 0:
                speedup_multi = std_time / fast_multi_time
                print(f"  Fast (multi) vs Standard: {speedup_multi:.1f}x faster")
        
        if results[seq_name].get("fast_single") and results[seq_name].get("fast_multi"):
            fast_single_time = results[seq_name]["fast_single"]["time"]
            fast_multi_time = results[seq_name]["fast_multi"]["time"]
            if fast_multi_time > 0:
                speedup_parallel = fast_single_time / fast_multi_time
                print(f"  Parallel speedup: {speedup_parallel:.1f}x")
    
    print("\n" + "="*80)
    print("Benchmark complete!")
    print("="*80)


if __name__ == "__main__":
    run_benchmarks()
