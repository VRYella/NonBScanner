#!/usr/bin/env python3
"""Profile individual detectors to find bottleneck"""

import time
from motif_detection import (
    CurvedDNADetector,
    SlippedDNADetector,
    CruciformDetector,
    RLoopDetector,
    TriplexDetector,
    GQuadruplexDetector,
    IMotifDetector,
    ZDNADetector,
    APhilicDetector
)

# Test sequence
seq = "ATCG" * 1250  # 5000 bp

detectors = {
    'CurvedDNA': CurvedDNADetector(),
    'SlippedDNA': SlippedDNADetector(),
    'Cruciform': CruciformDetector(),
    'RLoop': RLoopDetector(),
    'Triplex': TriplexDetector(),
    'GQuadruplex': GQuadruplexDetector(),
    'IMotif': IMotifDetector(),
    'ZDNA': ZDNADetector(),
    'APhilic': APhilicDetector()
}

print("Testing 5000 bp sequence on individual detectors...\n")

for name, detector in detectors.items():
    print(f"{name}...")
    start = time.time()
    try:
        motifs = detector.detect_motifs(seq, "test")
        elapsed = time.time() - start
        print(f"  Time: {elapsed:.3f}s, Motifs: {len(motifs)}")
    except Exception as e:
        elapsed = time.time() - start
        print(f"  ERROR after {elapsed:.3f}s: {str(e)[:100]}")
