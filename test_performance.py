#!/usr/bin/env python3
"""
Performance Testing Suite for NBDScanner
=========================================

This script performs rigorous performance testing on 100,000 nucleotide sequences
to identify the best performing implementation and optimization opportunities.

Test Coverage:
- Pattern matching efficiency
- Memory usage
- Execution time
- Detection accuracy
- Scalability analysis

Author: Dr. Venkata Rajesh Yella
"""

import time
import random
import tracemalloc
from typing import Dict, Any, List, Tuple
import numpy as np

# Import both implementations
from nbdscanner import analyze_sequence as nbdscanner_analyze
from modular_scanner import analyze_sequence as modular_analyze


def generate_test_sequence(length: int, motif_density: str = 'medium') -> str:
    """
    Generate realistic test sequence with various Non-B DNA motifs.
    
    Args:
        length: Sequence length in base pairs
        motif_density: 'low', 'medium', or 'high' motif density
    
    Returns:
        DNA sequence string
    """
    sequence_parts = []
    remaining = length
    
    # Define motif templates with realistic structures
    motif_templates = {
        'g4': ['GGGTTAGGGTTAGGGTTAGGG', 'GGGAGGGTGGGAGGGT'],
        'curved': ['AAAAAAAA', 'AAAAAATTTT', 'AAAAAAAAAAATTTTTTTTTT'],
        'z_dna': ['CGCGCGCGCGCG', 'GCGCGCGCGCGC'],
        'str': ['CACACACACACA', 'CGCGCGCGCG', 'GAGAGAGAGA'],
        'triplex': ['GAAGAAGAAGAAGAA', 'TCCTCCTCCTCC'],
        'r_loop': ['GGGCCCGGGCCC', 'GGGGGCCCCC'],
        'a_philic': ['AAAAAAAAAA'] * 3,
        'i_motif': ['CCCTTCCCTTCCCTTCCC', 'CCCTCCCTCCCTCCC']
    }
    
    # Set motif frequency based on density
    if motif_density == 'low':
        motif_frequency = 0.05  # 5% motifs
    elif motif_density == 'medium':
        motif_frequency = 0.15  # 15% motifs
    else:  # high
        motif_frequency = 0.30  # 30% motifs
    
    # Generate sequence with interspersed motifs
    while remaining > 0:
        if random.random() < motif_frequency and remaining > 50:
            # Insert a motif
            motif_type = random.choice(list(motif_templates.keys()))
            motif = random.choice(motif_templates[motif_type])
            if len(motif) <= remaining:
                sequence_parts.append(motif)
                remaining -= len(motif)
        else:
            # Insert random DNA
            chunk_size = min(random.randint(20, 100), remaining)
            random_chunk = ''.join(random.choices('ATGC', k=chunk_size))
            sequence_parts.append(random_chunk)
            remaining -= chunk_size
    
    return ''.join(sequence_parts)


def measure_performance(
    func: callable,
    sequence: str,
    sequence_name: str,
    runs: int = 3
) -> Dict[str, Any]:
    """
    Measure performance metrics for a detection function.
    
    Args:
        func: Function to test
        sequence: Test sequence
        sequence_name: Sequence identifier
        runs: Number of test runs
    
    Returns:
        Performance metrics dictionary
    """
    times = []
    memory_peaks = []
    results_list = []
    
    for i in range(runs):
        # Start memory tracking
        tracemalloc.start()
        
        # Measure execution time
        start_time = time.time()
        try:
            result = func(sequence, sequence_name)
            elapsed = time.time() - start_time
            times.append(elapsed)
            results_list.append(result)
            
            # Get peak memory usage
            current, peak = tracemalloc.get_traced_memory()
            memory_peaks.append(peak / 1024 / 1024)  # Convert to MB
        except Exception as e:
            print(f"Error in run {i+1}: {str(e)}")
            times.append(float('inf'))
            memory_peaks.append(float('inf'))
            results_list.append([])
        finally:
            tracemalloc.stop()
    
    # Calculate statistics
    avg_time = np.mean(times)
    std_time = np.std(times)
    avg_memory = np.mean(memory_peaks)
    
    # Get motif counts from first successful run
    motif_count = len(results_list[0]) if results_list[0] else 0
    
    return {
        'avg_time': avg_time,
        'std_time': std_time,
        'min_time': min(times),
        'max_time': max(times),
        'avg_memory_mb': avg_memory,
        'motif_count': motif_count,
        'all_times': times,
        'all_memory': memory_peaks
    }


def run_comprehensive_test():
    """
    Run comprehensive performance testing suite.
    """
    print("=" * 80)
    print("NBDScanner Performance Testing Suite")
    print("=" * 80)
    print()
    
    # Test configurations - start with smaller tests
    test_configs = [
        {'length': 1000, 'density': 'medium', 'name': '1K_medium'},
        {'length': 5000, 'density': 'medium', 'name': '5K_medium'},
        {'length': 10000, 'density': 'low', 'name': '10K_low'},
        {'length': 100000, 'density': 'low', 'name': '100K_low'},
    ]
    
    results_summary = []
    
    for config in test_configs:
        print(f"\nTest: {config['name']} ({config['length']} bp, {config['density']} density)")
        print("-" * 80)
        
        # Generate test sequence
        print(f"Generating test sequence...")
        test_seq = generate_test_sequence(config['length'], config['density'])
        print(f"Generated {len(test_seq)} bp sequence")
        
        # Test nbdscanner implementation
        print("\nTesting nbdscanner.py implementation...")
        nbdscanner_metrics = measure_performance(
            nbdscanner_analyze,
            test_seq,
            config['name'],
            runs=1
        )
        
        print(f"  Average time: {nbdscanner_metrics['avg_time']:.3f}s ± {nbdscanner_metrics['std_time']:.3f}s")
        print(f"  Memory usage: {nbdscanner_metrics['avg_memory_mb']:.2f} MB")
        print(f"  Motifs found: {nbdscanner_metrics['motif_count']}")
        
        # Test modular_scanner implementation
        print("\nTesting modular_scanner.py implementation...")
        modular_metrics = measure_performance(
            modular_analyze,
            test_seq,
            config['name'],
            runs=1
        )
        
        print(f"  Average time: {modular_metrics['avg_time']:.3f}s ± {modular_metrics['std_time']:.3f}s")
        print(f"  Memory usage: {modular_metrics['avg_memory_mb']:.2f} MB")
        print(f"  Motifs found: {modular_metrics['motif_count']}")
        
        # Calculate performance comparison
        speedup = nbdscanner_metrics['avg_time'] / modular_metrics['avg_time']
        memory_ratio = nbdscanner_metrics['avg_memory_mb'] / modular_metrics['avg_memory_mb']
        
        print(f"\nComparison:")
        print(f"  Speed ratio (nbdscanner/modular): {speedup:.2f}x")
        print(f"  Memory ratio (nbdscanner/modular): {memory_ratio:.2f}x")
        
        if speedup > 1.1:
            print(f"  ✓ modular_scanner is {speedup:.2f}x FASTER")
        elif speedup < 0.9:
            print(f"  ✓ nbdscanner is {1/speedup:.2f}x FASTER")
        else:
            print(f"  ≈ Similar performance")
        
        results_summary.append({
            'config': config['name'],
            'length': config['length'],
            'nbdscanner_time': nbdscanner_metrics['avg_time'],
            'modular_time': modular_metrics['avg_time'],
            'nbdscanner_memory': nbdscanner_metrics['avg_memory_mb'],
            'modular_memory': modular_metrics['avg_memory_mb'],
            'nbdscanner_motifs': nbdscanner_metrics['motif_count'],
            'modular_motifs': modular_metrics['motif_count'],
            'speedup': speedup
        })
    
    # Print summary
    print("\n" + "=" * 80)
    print("PERFORMANCE SUMMARY")
    print("=" * 80)
    print()
    print(f"{'Test':<20} {'Length':<10} {'NBD Time':<12} {'Mod Time':<12} {'Speedup':<10}")
    print("-" * 80)
    
    for result in results_summary:
        print(f"{result['config']:<20} {result['length']:<10} "
              f"{result['nbdscanner_time']:>10.3f}s {result['modular_time']:>10.3f}s "
              f"{result['speedup']:>8.2f}x")
    
    # Determine winner
    avg_speedup = np.mean([r['speedup'] for r in results_summary])
    
    print("\n" + "=" * 80)
    print("RECOMMENDATION")
    print("=" * 80)
    
    if avg_speedup > 1.1:
        print(f"\n✓ RECOMMENDED: Use modular_scanner.py")
        print(f"  Average speedup: {avg_speedup:.2f}x faster than nbdscanner.py")
    elif avg_speedup < 0.9:
        print(f"\n✓ RECOMMENDED: Use nbdscanner.py")
        print(f"  Average speedup: {1/avg_speedup:.2f}x faster than modular_scanner.py")
    else:
        print(f"\n≈ Both implementations have similar performance")
        print(f"  Consider using modular_scanner.py for better code organization")
    
    print("\n" + "=" * 80)
    
    return results_summary


def test_specific_motif_classes():
    """
    Test performance of specific motif class detectors.
    """
    print("\n" + "=" * 80)
    print("MOTIF CLASS SPECIFIC PERFORMANCE TESTS")
    print("=" * 80)
    
    # Generate sequences enriched for specific motifs
    test_cases = {
        'G4-rich': 'ATCG' * 10 + 'GGGTTAGGGTTAGGGTTAGGG' * 50 + 'ATCG' * 10,
        'Z-DNA-rich': 'ATCG' * 10 + 'CGCGCGCGCGCGCG' * 50 + 'ATCG' * 10,
        'Curved-rich': 'ATCG' * 10 + 'AAAAAAAATTTTTTTT' * 50 + 'ATCG' * 10,
        'Mixed': generate_test_sequence(10000, 'high')
    }
    
    for name, seq in test_cases.items():
        print(f"\nTest case: {name} ({len(seq)} bp)")
        print("-" * 40)
        
        start = time.time()
        motifs = modular_analyze(seq, name)
        elapsed = time.time() - start
        
        print(f"  Time: {elapsed:.3f}s")
        print(f"  Motifs found: {len(motifs)}")
        
        # Count by class
        class_counts = {}
        for m in motifs:
            cls = m.get('Class', 'Unknown')
            class_counts[cls] = class_counts.get(cls, 0) + 1
        
        print(f"  Class distribution:")
        for cls, count in sorted(class_counts.items(), key=lambda x: -x[1]):
            print(f"    {cls}: {count}")


if __name__ == "__main__":
    print("\n🧬 Starting NBDScanner Performance Testing Suite\n")
    
    try:
        # Run comprehensive tests
        results = run_comprehensive_test()
        
        # Run motif-class specific tests
        test_specific_motif_classes()
        
        print("\n✓ Performance testing completed successfully!\n")
        
    except Exception as e:
        print(f"\n✗ Error during testing: {str(e)}")
        import traceback
        traceback.print_exc()
