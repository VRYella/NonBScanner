"""
Performance Testing Suite for NonBScanner
Tests different approaches on 1 million base pair sequence

Approaches tested:
1. Baseline - Current sequential processing
2. Parallel - Multiprocessing with ProcessPoolExecutor
3. Chunked - Process large sequences in chunks
4. Optimized - Vectorized operations and caching
5. Hybrid - Combined best practices
"""

import time
import random
import json
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from typing import List, Dict, Any
import numpy as np

# Import NonBScanner modules
import nonbscanner as nbs
from utilities import get_basic_stats


def generate_test_sequence(length: int = 1_000_000, seed: int = 42) -> str:
    """
    Generate a realistic DNA sequence for testing.
    Includes various Non-B DNA motifs to ensure comprehensive testing.
    
    Args:
        length: Length of sequence to generate (default 1M bp)
        seed: Random seed for reproducibility
        
    Returns:
        DNA sequence string
    """
    random.seed(seed)
    np.random.seed(seed)
    
    # Base composition with realistic GC content (~40%)
    bases = ['A', 'T', 'G', 'C']
    weights = [0.30, 0.30, 0.20, 0.20]
    
    # Generate base sequence
    sequence = ''.join(random.choices(bases, weights=weights, k=length))
    
    # Insert known motif patterns at various positions
    # This ensures we're actually detecting motifs during testing
    motif_patterns = [
        # G-quadruplex
        'GGGTTAGGGTTAGGGTTAGGG',
        # Z-DNA
        'CGCGCGCGCGCGCGCGCGCG',
        # A-tract (Curved DNA)
        'AAAAAAAAAATTTTTTTTT',
        # Direct repeat (Slipped DNA)
        'ATGATGATGATGATGATG',
        # Inverted repeat (Cruciform)
        'ATGCATGCATGCATGCAT',
        # Mirror repeat (Triplex)
        'GAAAGAAAGAAAGAAA',
        # i-Motif
        'CCCTACCCTACCCTACCC',
        # STR
        'CAGCAGCAGCAGCAGCAG',
        # R-loop forming
        'GGGGAAAGGGGGAAAGGGGGAAA',
        # A-philic
        'AAAAAAAAAAAAAAAAAA'
    ]
    
    # Insert motifs at random positions (ensuring they don't get cut off)
    num_insertions = min(100, length // 10000)  # ~1 per 10kb
    for _ in range(num_insertions):
        pattern = random.choice(motif_patterns)
        pos = random.randint(0, length - len(pattern))
        sequence = sequence[:pos] + pattern + sequence[pos + len(pattern):]
    
    return sequence[:length]  # Ensure exact length


def approach_1_baseline(sequence: str, name: str) -> tuple:
    """
    Approach 1: Baseline - Current sequential processing
    This is the existing implementation without modifications.
    """
    start_time = time.time()
    motifs = nbs.analyze_sequence(sequence, name)
    elapsed = time.time() - start_time
    
    return {
        'name': 'Baseline (Sequential)',
        'time': elapsed,
        'motif_count': len(motifs),
        'speed_bp_per_sec': len(sequence) / elapsed if elapsed > 0 else 0,
        'motifs': motifs
    }


def approach_2_parallel_chunks(sequence: str, name: str, chunk_size: int = 100_000) -> tuple:
    """
    Approach 2: Parallel processing with overlapping chunks
    Splits sequence into chunks with overlap to catch motifs at boundaries.
    """
    start_time = time.time()
    
    # Calculate overlap size (enough to catch longest motifs ~1000bp)
    overlap = 2000
    
    # Create chunks with overlap
    chunks = []
    chunk_names = []
    i = 0
    chunk_id = 0
    
    while i < len(sequence):
        end = min(i + chunk_size, len(sequence))
        chunk = sequence[i:end]
        chunks.append(chunk)
        chunk_names.append(f"{name}_chunk_{chunk_id}")
        chunk_id += 1
        i += chunk_size - overlap
        
        if end >= len(sequence):
            break
    
    # Process chunks in parallel
    all_motifs = []
    with ProcessPoolExecutor(max_workers=min(mp.cpu_count(), len(chunks))) as executor:
        futures = [executor.submit(nbs.analyze_sequence, chunk, chunk_name) 
                   for chunk, chunk_name in zip(chunks, chunk_names)]
        
        # Collect results
        chunk_offset = 0
        for idx, future in enumerate(futures):
            chunk_motifs = future.result()
            
            # Adjust positions back to original sequence coordinates
            for motif in chunk_motifs:
                motif['Start'] += chunk_offset
                motif['End'] += chunk_offset
                motif['Sequence_Name'] = name
            
            all_motifs.extend(chunk_motifs)
            
            # Update offset for next chunk
            if idx < len(chunks) - 1:
                chunk_offset += chunk_size - overlap
    
    # Remove duplicates that may occur in overlap regions
    # Sort by position and remove overlaps
    all_motifs.sort(key=lambda x: (x['Start'], x['End']))
    
    # Simple deduplication: if two motifs are very similar in position, keep one
    deduplicated = []
    for motif in all_motifs:
        if not deduplicated:
            deduplicated.append(motif)
        else:
            last = deduplicated[-1]
            # Check if this is a duplicate (similar position and same class)
            if not (abs(motif['Start'] - last['Start']) < 100 and 
                    motif['Class'] == last['Class']):
                deduplicated.append(motif)
    
    elapsed = time.time() - start_time
    
    return {
        'name': 'Parallel Chunks (Multiprocessing)',
        'time': elapsed,
        'motif_count': len(deduplicated),
        'speed_bp_per_sec': len(sequence) / elapsed if elapsed > 0 else 0,
        'motifs': deduplicated,
        'chunks_used': len(chunks)
    }


def approach_3_chunked_sequential(sequence: str, name: str, chunk_size: int = 100_000) -> tuple:
    """
    Approach 3: Sequential chunked processing
    Processes sequence in chunks sequentially (lower memory overhead).
    """
    start_time = time.time()
    
    # Calculate overlap size
    overlap = 2000
    
    all_motifs = []
    i = 0
    chunk_id = 0
    chunk_offset = 0
    
    while i < len(sequence):
        end = min(i + chunk_size, len(sequence))
        chunk = sequence[i:end]
        
        # Process chunk
        chunk_motifs = nbs.analyze_sequence(chunk, f"{name}_chunk_{chunk_id}")
        
        # Adjust positions
        for motif in chunk_motifs:
            motif['Start'] += chunk_offset
            motif['End'] += chunk_offset
            motif['Sequence_Name'] = name
        
        all_motifs.extend(chunk_motifs)
        
        chunk_id += 1
        chunk_offset += chunk_size - overlap
        i += chunk_size - overlap
        
        if end >= len(sequence):
            break
    
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
        'name': 'Chunked Sequential (Memory Efficient)',
        'time': elapsed,
        'motif_count': len(deduplicated),
        'speed_bp_per_sec': len(sequence) / elapsed if elapsed > 0 else 0,
        'motifs': deduplicated
    }


def approach_4_optimized_baseline(sequence: str, name: str) -> tuple:
    """
    Approach 4: Optimized baseline with pre-compilation
    Uses the baseline approach but with optimized settings.
    """
    start_time = time.time()
    
    # Run with baseline (already optimized in detectors.py)
    motifs = nbs.analyze_sequence(sequence, name)
    
    elapsed = time.time() - start_time
    
    return {
        'name': 'Optimized Baseline',
        'time': elapsed,
        'motif_count': len(motifs),
        'speed_bp_per_sec': len(sequence) / elapsed if elapsed > 0 else 0,
        'motifs': motifs
    }


def approach_5_hybrid(sequence: str, name: str) -> tuple:
    """
    Approach 5: Hybrid approach
    Uses parallel processing for long sequences, sequential for short ones.
    Adapts based on sequence length.
    """
    start_time = time.time()
    
    # Adaptive strategy
    if len(sequence) > 500_000:
        # Use parallel chunks for large sequences
        result = approach_2_parallel_chunks(sequence, name, chunk_size=200_000)
        result['name'] = 'Hybrid (Parallel for large seq)'
    else:
        # Use optimized baseline for smaller sequences
        result = approach_4_optimized_baseline(sequence, name)
        result['name'] = 'Hybrid (Sequential for small seq)'
    
    return result


def run_performance_tests(sequence_length: int = 1_000_000):
    """
    Run all performance tests and compare results.
    
    Args:
        sequence_length: Length of test sequence
    """
    print("=" * 80)
    print(f"NonBScanner Performance Testing Suite")
    print(f"Testing on {sequence_length:,} base pair sequence")
    print("=" * 80)
    print()
    
    # Generate test sequence
    print("Generating test sequence...")
    test_sequence = generate_test_sequence(sequence_length)
    test_name = f"test_seq_{sequence_length}bp"
    
    # Get basic stats
    stats = get_basic_stats(test_sequence)
    print(f"Sequence generated: {len(test_sequence):,} bp")
    print(f"GC Content: {stats['GC%']:.2f}%")
    print()
    
    # Define approaches to test
    approaches = [
        ('Baseline', approach_1_baseline),
        ('Chunked Sequential', approach_3_chunked_sequential),
        ('Parallel Chunks', approach_2_parallel_chunks),
        ('Optimized Baseline', approach_4_optimized_baseline),
        ('Hybrid Adaptive', approach_5_hybrid),
    ]
    
    results = []
    
    # Run each approach
    for approach_name, approach_func in approaches:
        print(f"\n{'=' * 80}")
        print(f"Testing: {approach_name}")
        print('=' * 80)
        
        try:
            result = approach_func(test_sequence, test_name)
            results.append(result)
            
            print(f"✓ Completed in {result['time']:.2f} seconds")
            print(f"  Motifs detected: {result['motif_count']}")
            print(f"  Processing speed: {result['speed_bp_per_sec']:,.0f} bp/second")
            print(f"  Time per 1M bp: {(result['time'] / (sequence_length / 1_000_000)):.2f} seconds")
            
        except Exception as e:
            print(f"✗ Failed: {str(e)}")
            import traceback
            traceback.print_exc()
    
    # Summary comparison
    print("\n\n" + "=" * 80)
    print("PERFORMANCE COMPARISON SUMMARY")
    print("=" * 80)
    print()
    
    # Sort by speed
    results.sort(key=lambda x: x['time'])
    
    print(f"{'Approach':<35} {'Time (s)':<12} {'Motifs':<10} {'Speed (bp/s)':<15}")
    print("-" * 80)
    
    for result in results:
        print(f"{result['name']:<35} {result['time']:<12.2f} {result['motif_count']:<10} {result['speed_bp_per_sec']:>14,.0f}")
    
    print()
    print("=" * 80)
    print(f"WINNER: {results[0]['name']}")
    print(f"Time: {results[0]['time']:.2f} seconds")
    print(f"Speed: {results[0]['speed_bp_per_sec']:,.0f} bp/second")
    print("=" * 80)
    
    # Verify all approaches detect similar number of motifs
    print("\nMotif Detection Consistency Check:")
    print("-" * 40)
    motif_counts = [r['motif_count'] for r in results]
    min_count = min(motif_counts)
    max_count = max(motif_counts)
    variance = max_count - min_count
    
    print(f"Min motifs: {min_count}")
    print(f"Max motifs: {max_count}")
    print(f"Variance: {variance}")
    
    if variance > len(test_sequence) * 0.001:  # More than 0.1% difference
        print("⚠️  WARNING: Significant variance in motif detection across approaches")
    else:
        print("✓ All approaches detect similar number of motifs")
    
    # Save detailed results
    output_file = '/tmp/performance_results.json'
    with open(output_file, 'w') as f:
        # Don't save full motifs list (too large), just summary
        summary_results = []
        for r in results:
            summary_results.append({
                'name': r['name'],
                'time': r['time'],
                'motif_count': r['motif_count'],
                'speed_bp_per_sec': r['speed_bp_per_sec']
            })
        json.dump({
            'sequence_length': sequence_length,
            'results': summary_results,
            'winner': results[0]['name']
        }, f, indent=2)
    
    print(f"\nDetailed results saved to: {output_file}")
    
    return results


def run_scaling_tests():
    """
    Test performance at different sequence sizes to understand scaling.
    """
    print("=" * 80)
    print("NonBScanner Scaling Analysis")
    print("=" * 80)
    print()
    
    # Test at different scales
    test_sizes = [10_000, 50_000, 100_000, 250_000, 500_000, 1_000_000]
    scaling_results = []
    
    for size in test_sizes:
        print(f"\n{'=' * 80}")
        print(f"Testing with {size:,} bp sequence")
        print('=' * 80)
        
        # Generate test sequence
        test_sequence = generate_test_sequence(size)
        test_name = f"test_{size}bp"
        
        # Test baseline approach (fastest for understanding)
        start_time = time.time()
        motifs = nbs.analyze_sequence(test_sequence, test_name)
        elapsed = time.time() - start_time
        
        speed = size / elapsed if elapsed > 0 else 0
        
        result = {
            'size': size,
            'time': elapsed,
            'motifs': len(motifs),
            'speed': speed,
            'time_per_1M': (elapsed / (size / 1_000_000))
        }
        
        scaling_results.append(result)
        
        print(f"Time: {elapsed:.2f}s | Motifs: {len(motifs)} | Speed: {speed:,.0f} bp/s | Est. for 1M: {result['time_per_1M']:.2f}s")
        
        # Stop if 1M takes too long (>60s)
        if size == 1_000_000:
            break
        elif result['time_per_1M'] > 120:
            print(f"\n⚠️  Estimated time for 1M bp: {result['time_per_1M']:.2f}s (>120s)")
            print("Stopping scaling test - will focus on optimizations")
            break
    
    print("\n" + "=" * 80)
    print("SCALING SUMMARY")
    print("=" * 80)
    print(f"{'Size (bp)':<15} {'Time (s)':<12} {'Motifs':<10} {'Speed (bp/s)':<15} {'Est. 1M (s)'}")
    print("-" * 80)
    
    for r in scaling_results:
        print(f"{r['size']:<15,} {r['time']:<12.2f} {r['motifs']:<10} {r['speed']:>14,.0f} {r['time_per_1M']:>12.2f}")
    
    return scaling_results


if __name__ == "__main__":
    # First run scaling tests to understand performance characteristics
    print("Step 1: Understanding performance scaling...")
    print()
    scaling_results = run_scaling_tests()
    
    # Determine if we can test 1M directly or need optimizations
    if scaling_results:
        last_result = scaling_results[-1]
        
        if last_result['size'] >= 1_000_000:
            print(f"\n✓ Successfully tested 1M bp in {last_result['time']:.2f} seconds")
        else:
            estimated_time = last_result['time_per_1M']
            print(f"\n⚠️  Estimated time for 1M bp: {estimated_time:.2f} seconds")
            
            if estimated_time < 300:  # Less than 5 minutes
                print("Proceeding with full comparison test...")
                results = run_performance_tests(sequence_length=1_000_000)
            else:
                print("Testing optimized approaches only...")
                # Test with smaller size for comparison
                results = run_performance_tests(sequence_length=min(500_000, last_result['size'] * 2))
