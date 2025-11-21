#!/usr/bin/env python3
"""
Quick test to measure actual performance on 1M bp sequence
Uses the fastest approach (parallel chunked processing)
"""

import time
import random
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
import nonbscanner as nbs

def generate_sequence(length=1_000_000):
    random.seed(42)
    bases = ['A', 'T', 'G', 'C']
    weights = [0.30, 0.30, 0.20, 0.20]
    sequence = ''.join(random.choices(bases, weights=weights, k=length))
    
    # Insert motifs
    motifs = ['GGGTTAGGGTTAGGGTTAGGG', 'CGCGCGCGCGCGCGCGCGCG', 'AAAAAAAAAATTTTTTTTT',
              'ATGATGATGATGATGATG', 'CCCTACCCTACCCTACCC', 'GAAAGAAAGAAAGAAA']
    for _ in range(100):
        pattern = random.choice(motifs)
        pos = random.randint(0, length - len(pattern))
        sequence = sequence[:pos] + pattern + sequence[pos + len(pattern):]
    
    return sequence[:length]

def process_chunk(args):
    chunk, name, offset = args
    motifs = nbs.analyze_sequence(chunk, name)
    for m in motifs:
        m['Start'] += offset
        m['End'] += offset
    return motifs

print("Testing NonBScanner on 1M bp sequence...")
print("=" * 60)

sequence = generate_sequence(1_000_000)
print(f"Generated {len(sequence):,} bp sequence")

# Parallel chunked approach
chunk_size = 100_000
overlap = 2000
num_workers = 4

chunks = []
i = 0
while i < len(sequence):
    end = min(i + chunk_size, len(sequence))
    chunks.append((sequence[i:end], f"chunk_{len(chunks)}", i))
    i += chunk_size - overlap
    if end >= len(sequence):
        break

print(f"Processing with {len(chunks)} chunks, {num_workers} workers...")

start = time.time()
with ProcessPoolExecutor(max_workers=num_workers) as executor:
    results = list(executor.map(process_chunk, chunks))

all_motifs = []
for chunk_motifs in results:
    all_motifs.extend(chunk_motifs)

# Deduplicate
all_motifs.sort(key=lambda x: (x['Start'], x['End']))
deduplicated = []
for motif in all_motifs:
    if not deduplicated or not (abs(motif['Start'] - deduplicated[-1]['Start']) < 100 and 
                                 motif['Class'] == deduplicated[-1]['Class']):
        deduplicated.append(motif)

elapsed = time.time() - start

print("=" * 60)
print(f"✓ COMPLETED in {elapsed:.2f} seconds")
print(f"  Motifs detected: {len(deduplicated)}")
print(f"  Processing speed: {len(sequence)/elapsed:,.0f} bp/second")
print("=" * 60)
