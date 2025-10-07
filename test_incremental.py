#!/usr/bin/env python3
"""Test with progressively larger sequences to find bottleneck"""

import time
from modular_scanner import analyze_sequence

def test_size(size, name):
    print(f"\nTesting {name} ({size} bp)...")
    seq = "ATCG" * (size // 4)
    start = time.time()
    motifs = analyze_sequence(seq, name)
    elapsed = time.time() - start
    print(f"  Time: {elapsed:.3f}s, Motifs: {len(motifs)}, Speed: {size/elapsed:.0f} bp/s")
    return elapsed

print("Testing incrementally larger sequences...")
test_size(100, "100bp")
test_size(500, "500bp")
test_size(1000, "1Kbp")
test_size(5000, "5Kbp")
test_size(10000, "10Kbp")
