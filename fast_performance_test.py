"""
Fast Performance Testing for NonBScanner on 1M bp sequence
Uses only optimized approaches that can handle large sequences efficiently
"""

import time
import random
import json
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
import sys

# Import NonBScanner modules
import nonbscanner as nbs
from utilities import get_basic_stats


def generate_test_sequence(length: int = 1_000_000, seed: int = 42) -> str:
    """Generate a realistic DNA sequence for testing."""
    random.seed(seed)
    
    # Base composition with realistic GC content (~40%)
    bases = ['A', 'T', 'G', 'C']
    weights = [0.30, 0.30, 0.20, 0.20]
    
    # Generate base sequence
    sequence = ''.join(random.choices(bases, weights=weights, k=length))
    
    # Insert known motif patterns
    motif_patterns = [
        'GGGTTAGGGTTAGGGTTAGGG',  # G-quadruplex
        'CGCGCGCGCGCGCGCGCGCG',  # Z-DNA
        'AAAAAAAAAATTTTTTTTT',   # A-tract (Curved DNA)
        'ATGATGATGATGATGATG',    # Direct repeat
        'ATGCATGCATGCATGCAT',    # Inverted repeat
        'GAAAGAAAGAAAGAAA',      # Mirror repeat
        'CCCTACCCTACCCTACCC',    # i-Motif
        'CAGCAGCAGCAGCAGCAG',    # STR
        'GGGGAAAGGGGGAAAGGGGGAAA', # R-loop
        'AAAAAAAAAAAAAAAAAA'     # A-philic
    ]
    
    # Insert motifs at random positions
    num_insertions = min(100, length // 10000)
    for _ in range(num_insertions):
        pattern = random.choices(motif_patterns)[0]
        pos = random.randint(0, length - len(pattern))
        sequence = sequence[:pos] + pattern + sequence[pos + len(pattern):]
    
    return sequence[:length]


def parallel_chunk_analysis(sequence: str, name: str, chunk_size: int = 100_000, num_workers: int = None) -> dict:
    """
    Parallel processing with overlapping chunks - optimized version.
    """
    if num_workers is None:
        num_workers = min(mp.cpu_count(), 8)  # Limit to 8 workers
    
    start_time = time.time()
    
    # Overlap to catch motifs at boundaries
    overlap = 2000
    
    # Create chunks
    chunks = []
    chunk_positions = []
    i = 0
    chunk_id = 0
    
    while i < len(sequence):
        end = min(i + chunk_size, len(sequence))
        chunk = sequence[i:end]
        chunks.append((chunk, f"{name}_chunk_{chunk_id}", i))
        chunk_positions.append(i)
        chunk_id += 1
        i += chunk_size - overlap
        
        if end >= len(sequence):
            break
    
    print(f"  Created {len(chunks)} chunks with {num_workers} workers")
    
    # Process chunks in parallel
    def process_chunk(args):
        chunk, chunk_name, offset = args
        motifs = nbs.analyze_sequence(chunk, chunk_name)
        # Adjust positions
        for motif in motifs:
            motif['Start'] += offset
            motif['End'] += offset
            motif['Sequence_Name'] = name
        return motifs
    
    all_motifs = []
    
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        chunk_results = list(executor.map(process_chunk, chunks))
        for chunk_motifs in chunk_results:
            all_motifs.extend(chunk_motifs)
    
    # Simple deduplication
    all_motifs.sort(key=lambda x: (x['Start'], x['End']))
    deduplicated = []
    for motif in all_motifs:
        if not deduplicated:
            deduplicated.append(motif)
        else:
            last = deduplicated[-1]
            if not (abs(motif['Start'] - last['Start']) < 100 and 
                    motif['Class'] == last['Class']):
                deduplicated.append(motif)
    
    elapsed = time.time() - start_time
    
    return {
        'name': f'Parallel Chunks ({num_workers} workers)',
        'time': elapsed,
        'motif_count': len(deduplicated),
        'speed_bp_per_sec': len(sequence) / elapsed if elapsed > 0 else 0,
        'motifs': deduplicated,
        'chunks': len(chunks),
        'workers': num_workers
    }


def sequential_chunk_analysis(sequence: str, name: str, chunk_size: int = 150_000) -> dict:
    """
    Sequential chunked processing - memory efficient.
    """
    start_time = time.time()
    
    overlap = 2000
    all_motifs = []
    i = 0
    chunk_id = 0
    
    total_chunks = (len(sequence) + chunk_size - overlap - 1) // (chunk_size - overlap)
    
    while i < len(sequence):
        end = min(i + chunk_size, len(sequence))
        chunk = sequence[i:end]
        
        chunk_motifs = nbs.analyze_sequence(chunk, f"{name}_chunk_{chunk_id}")
        
        # Adjust positions
        for motif in chunk_motifs:
            motif['Start'] += i
            motif['End'] += i
            motif['Sequence_Name'] = name
        
        all_motifs.extend(chunk_motifs)
        
        if chunk_id % 5 == 0:
            sys.stdout.write(f"\r  Processing chunk {chunk_id + 1}/{total_chunks}")
            sys.stdout.flush()
        
        chunk_id += 1
        i += chunk_size - overlap
        
        if end >= len(sequence):
            break
    
    print(f"\r  Processed {chunk_id} chunks" + " " * 20)
    
    # Deduplicate
    all_motifs.sort(key=lambda x: (x['Start'], x['End']))
    deduplicated = []
    for motif in all_motifs:
        if not deduplicated:
            deduplicated.append(motif)
        else:
            last = deduplicated[-1]
            if not (abs(motif['Start'] - last['Start']) < 100 and 
                    motif['Class'] == last['Class']):
                deduplicated.append(motif)
    
    elapsed = time.time() - start_time
    
    return {
        'name': 'Sequential Chunks',
        'time': elapsed,
        'motif_count': len(deduplicated),
        'speed_bp_per_sec': len(sequence) / elapsed if elapsed > 0 else 0,
        'motifs': deduplicated,
        'chunks': chunk_id
    }


def main():
    """Run optimized performance tests on 1M bp sequence."""
    print("=" * 80)
    print("NonBScanner - Optimized Performance Test for 1M bp Sequence")
    print("=" * 80)
    print()
    
    # Generate test sequence
    print("Generating 1,000,000 bp test sequence...")
    test_sequence = generate_test_sequence(1_000_000)
    test_name = "test_1M_bp"
    
    stats = get_basic_stats(test_sequence)
    print(f"✓ Sequence generated: {len(test_sequence):,} bp")
    print(f"  GC Content: {stats['GC%']:.2f}%")
    print()
    
    results = []
    
    # Test 1: Sequential Chunks (most memory efficient)
    print("\n" + "=" * 80)
    print("Test 1: Sequential Chunked Processing")
    print("=" * 80)
    try:
        result = sequential_chunk_analysis(test_sequence, test_name, chunk_size=150_000)
        results.append(result)
        print(f"✓ Completed in {result['time']:.2f} seconds")
        print(f"  Motifs detected: {result['motif_count']}")
        print(f"  Processing speed: {result['speed_bp_per_sec']:,.0f} bp/second")
    except Exception as e:
        print(f"✗ Failed: {str(e)}")
    
    # Test 2: Parallel Chunks with different worker counts
    for num_workers in [2, 4, mp.cpu_count()]:
        print("\n" + "=" * 80)
        print(f"Test: Parallel Chunked Processing ({num_workers} workers)")
        print("=" * 80)
        try:
            result = parallel_chunk_analysis(test_sequence, test_name, 
                                            chunk_size=100_000, 
                                            num_workers=num_workers)
            results.append(result)
            print(f"✓ Completed in {result['time']:.2f} seconds")
            print(f"  Motifs detected: {result['motif_count']}")
            print(f"  Processing speed: {result['speed_bp_per_sec']:,.0f} bp/second")
        except Exception as e:
            print(f"✗ Failed: {str(e)}")
    
    # Summary
    print("\n\n" + "=" * 80)
    print("PERFORMANCE SUMMARY - 1 MILLION BASE PAIR SEQUENCE")
    print("=" * 80)
    print()
    
    if results:
        # Sort by time (fastest first)
        results.sort(key=lambda x: x['time'])
        
        print(f"{'Approach':<40} {'Time (s)':<12} {'Motifs':<10} {'Speed (bp/s)':<15}")
        print("-" * 80)
        
        for result in results:
            print(f"{result['name']:<40} {result['time']:<12.2f} {result['motif_count']:<10} {result['speed_bp_per_sec']:>14,.0f}")
        
        print()
        print("=" * 80)
        print("🏆 BEST APPROACH:")
        print(f"   {results[0]['name']}")
        print(f"   Time: {results[0]['time']:.2f} seconds")
        print(f"   Speed: {results[0]['speed_bp_per_sec']:,.0f} bp/second")
        print(f"   Motifs: {results[0]['motif_count']}")
        print("=" * 80)
        
        # Check consistency
        print("\nMotif Detection Consistency:")
        motif_counts = [r['motif_count'] for r in results]
        print(f"  Min motifs: {min(motif_counts)}")
        print(f"  Max motifs: {max(motif_counts)}")
        print(f"  Variance: {max(motif_counts) - min(motif_counts)}")
        
        # Save results
        output_file = '/tmp/performance_results_1M.json'
        with open(output_file, 'w') as f:
            summary = []
            for r in results:
                summary.append({
                    'name': r['name'],
                    'time': r['time'],
                    'motif_count': r['motif_count'],
                    'speed_bp_per_sec': r['speed_bp_per_sec'],
                    'chunks': r.get('chunks', 0),
                    'workers': r.get('workers', 1)
                })
            json.dump({
                'sequence_length': 1_000_000,
                'results': summary,
                'best_approach': results[0]['name'],
                'best_time': results[0]['time']
            }, f, indent=2)
        
        print(f"\n✓ Results saved to: {output_file}")
    else:
        print("No successful results to report")
    
    return results


if __name__ == "__main__":
    results = main()
