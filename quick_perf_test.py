#!/usr/bin/env python3
"""
Quick Performance Test for NBDScanner
Tests on 100,000 nucleotide sequence to identify best approach
"""

import time
import random

# Import modular scanner (appears to be the main implementation)
from modular_scanner import analyze_sequence

def generate_test_sequence(length: int) -> str:
    """Generate realistic test sequence with motifs."""
    motif_templates = {
        'g4': ['GGGTTAGGGTTAGGGTTAGGG', 'GGGAGGGTGGGAGGGT'],
        'curved': ['AAAAAAAA', 'AAAAAATTTT'],
        'z_dna': ['CGCGCGCGCGCG', 'GCGCGCGCGCGC'],
        'str': ['CACACACACACA', 'CGCGCGCGCG'],
    }
    
    sequence_parts = []
    remaining = length
    motif_frequency = 0.10  # 10% motifs
    
    while remaining > 0:
        if random.random() < motif_frequency and remaining > 50:
            motif_type = random.choice(list(motif_templates.keys()))
            motif = random.choice(motif_templates[motif_type])
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
print("NBDScanner Performance Test - 100,000 Nucleotide Sequence")
print("="*80)

# Generate 100K test sequence
print("\nGenerating 100,000 bp test sequence...")
test_seq = generate_test_sequence(100000)
print(f"Generated sequence: {len(test_seq)} bp")

# Run analysis
print("\nRunning motif detection analysis...")
start_time = time.time()
motifs = analyze_sequence(test_seq, "test_100k")
elapsed = time.time() - start_time

# Results
print(f"\n{'='*80}")
print("RESULTS")
print(f"{'='*80}")
print(f"Analysis time: {elapsed:.3f} seconds")
print(f"Total motifs found: {len(motifs)}")

# Count by class
class_counts = {}
for m in motifs:
    cls = m.get('Class', 'Unknown')
    class_counts[cls] = class_counts.get(cls, 0) + 1

print(f"\nMotif distribution by class:")
for cls, count in sorted(class_counts.items(), key=lambda x: -x[1]):
    print(f"  {cls}: {count}")

print(f"\nPerformance: {len(test_seq)/elapsed:.0f} bp/second")
print(f"{'='*80}\n")
