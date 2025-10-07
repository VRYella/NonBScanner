#!/usr/bin/env python3
"""Test 100K sequence with fast detectors only"""

import time
import random
from modular_scanner import ModularMotifDetector

def generate_test_sequence(length: int) -> str:
    """Generate realistic test sequence with motifs."""
    motif_templates = [
        'GGGTTAGGGTTAGGGTTAGGG',  # G4
        'CGCGCGCGCGCG',  # Z-DNA
        'AAAAAAAA',  # Curved
        'GGGCCCGGG',  # R-loop
    ]
    
    sequence_parts = []
    remaining = length
    motif_frequency = 0.05
    
    while remaining > 0:
        if random.random() < motif_frequency and remaining > 50:
            motif = random.choice(motif_templates)
            if len(motif) <= remaining:
                sequence_parts.append(motif)
                remaining -= len(motif)
        else:
            chunk_size = min(random.randint(20, 100), remaining)
            random_chunk = ''.join(random.choices('ATGC', k=chunk_size))
            sequence_parts.append(random_chunk)
            remaining -= chunk_size
    
    return ''.join(sequence_parts)

print("\n" + "="*80)
print("NBDScanner Performance Test - 100,000 bp (Fast Detectors Only)")
print("="*80)

# Generate 100K test sequence
print("\nGenerating 100,000 bp test sequence...")
test_seq = generate_test_sequence(100000)
print(f"Generated sequence: {len(test_seq)} bp")

# Create detector and disable slow ones
detector = ModularMotifDetector()
print("\nDisabling slow detectors (Cruciform, SlippedDNA)...")
del detector.detectors['cruciform']
del detector.detectors['slipped_dna']

print(f"Active detectors: {list(detector.detectors.keys())}")

# Run analysis
print("\nRunning motif detection analysis...")
start_time = time.time()
motifs = detector.analyze_sequence(test_seq, "test_100k", use_pure_python=False)
elapsed = time.time() - start_time

# Results
print(f"\n{'='*80}")
print("RESULTS")
print(f"{'='*80}")
print(f"Analysis time: {elapsed:.3f} seconds")
print(f"Total motifs found: {len(motifs)}")
print(f"Performance: {len(test_seq)/elapsed:.0f} bp/second")

# Count by class
class_counts = {}
for m in motifs:
    cls = m.get('Class', 'Unknown')
    class_counts[cls] = class_counts.get(cls, 0) + 1

print(f"\nMotif distribution by class:")
for cls, count in sorted(class_counts.items(), key=lambda x: -x[1]):
    print(f"  {cls}: {count}")

print(f"\n{'='*80}\n")
