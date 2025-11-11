#!/usr/bin/env python3
"""
Performance Test for NonBScanner
=================================
Test performance on sequences of various sizes up to 100MB
"""

import time
import random
from scanner import analyze_sequence
import sys

def generate_random_dna_sequence(length):
    """Generate a random DNA sequence with motif patterns"""
    bases = ['A', 'C', 'G', 'T']
    
    # Generate realistic sequence with some motifs
    sequence = []
    i = 0
    while i < length:
        # Random chance to add a motif
        rand = random.random()
        if rand < 0.01:  # 1% chance to add G4 pattern
            sequence.append('GGG' + 'ATCG' * random.randint(1, 3) + 'GGG' + 'ATCG' * random.randint(1, 3) + 'GGG' + 'ATCG' * random.randint(1, 3) + 'GGG')
            i += 40
        elif rand < 0.02:  # 1% chance to add A-tract
            sequence.append('A' * random.randint(5, 10))
            i += 8
        elif rand < 0.03:  # 1% chance to add STR
            unit = ''.join(random.choices(bases, k=random.randint(2, 4)))
            sequence.append(unit * random.randint(3, 6))
            i += 15
        else:
            sequence.append(random.choice(bases))
            i += 1
    
    return ''.join(sequence[:length])

def test_performance(size_mb):
    """Test performance on a sequence of given size in MB"""
    # 1 MB = ~1,000,000 bp
    size_bp = int(size_mb * 1_000_000)
    
    print(f"\n{'='*70}")
    print(f"Testing sequence size: {size_mb} MB ({size_bp:,} bp)")
    print('='*70)
    
    # Generate sequence
    print(f"Generating {size_bp:,} bp sequence...")
    seq_start = time.time()
    sequence = generate_random_dna_sequence(size_bp)
    seq_time = time.time() - seq_start
    print(f"  Generated in {seq_time:.2f} seconds")
    
    # Analyze sequence
    print(f"Analyzing sequence...")
    analysis_start = time.time()
    try:
        motifs = analyze_sequence(sequence, f"test_{size_mb}MB")
        analysis_time = time.time() - analysis_start
        
        # Calculate throughput
        throughput = size_bp / analysis_time
        
        print(f"\n  ✓ Analysis completed successfully!")
        print(f"  Time taken: {analysis_time:.2f} seconds")
        print(f"  Throughput: {throughput:,.0f} bp/s ({throughput/1_000_000:.2f} MB/s)")
        print(f"  Motifs detected: {len(motifs)}")
        
        # Show memory usage estimate
        import sys
        motif_memory = sys.getsizeof(motifs) / 1_000_000
        print(f"  Motif data size: {motif_memory:.2f} MB")
        
        return {
            'size_mb': size_mb,
            'size_bp': size_bp,
            'time_seconds': analysis_time,
            'throughput_bp_s': throughput,
            'motifs_found': len(motifs),
            'success': True
        }
        
    except Exception as e:
        analysis_time = time.time() - analysis_start
        print(f"\n  ✗ Analysis failed!")
        print(f"  Error: {str(e)}")
        print(f"  Time before failure: {analysis_time:.2f} seconds")
        
        return {
            'size_mb': size_mb,
            'size_bp': size_bp,
            'time_seconds': analysis_time,
            'error': str(e),
            'success': False
        }

def main():
    # Test sizes
    test_sizes = [0.1, 0.5, 1, 5, 10, 50, 100]  # MB
    
    if len(sys.argv) > 1:
        # Custom size provided
        test_sizes = [float(arg) for arg in sys.argv[1:]]
    
    print("\n" + "="*70)
    print("NonBScanner Performance Test")
    print("="*70)
    print(f"Testing sequence sizes: {', '.join(str(s) + ' MB' for s in test_sizes)}")
    
    results = []
    for size in test_sizes:
        result = test_performance(size)
        results.append(result)
        
        # Stop if we encounter an error
        if not result['success']:
            print(f"\n⚠ Stopping tests due to error at {size} MB")
            break
    
    # Summary
    print("\n" + "="*70)
    print("PERFORMANCE SUMMARY")
    print("="*70)
    print(f"{'Size (MB)':>10} {'Time (s)':>12} {'Throughput (bp/s)':>20} {'Motifs':>10} {'Status':>10}")
    print("-"*70)
    
    for result in results:
        if result['success']:
            print(f"{result['size_mb']:>10.1f} {result['time_seconds']:>12.2f} "
                  f"{result['throughput_bp_s']:>20,.0f} {result['motifs_found']:>10} "
                  f"{'✓':>10}")
        else:
            print(f"{result['size_mb']:>10.1f} {result['time_seconds']:>12.2f} "
                  f"{'N/A':>20} {'N/A':>10} {'✗':>10}")
    
    print("="*70)
    
    # Performance recommendations
    successful_results = [r for r in results if r['success']]
    if successful_results:
        avg_throughput = sum(r['throughput_bp_s'] for r in successful_results) / len(successful_results)
        print(f"\nAverage throughput: {avg_throughput:,.0f} bp/s ({avg_throughput/1_000_000:.2f} MB/s)")
        
        # Estimate time for 100MB
        if results[-1]['size_mb'] < 100:
            estimated_time = 100_000_000 / avg_throughput
            print(f"Estimated time for 100 MB: {estimated_time:.2f} seconds ({estimated_time/60:.2f} minutes)")

if __name__ == "__main__":
    main()
