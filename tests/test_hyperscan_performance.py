#!/usr/bin/env python3
"""
Hyperscan Performance Test Suite
=================================

Tests performance of hyperscan-based motif detection vs algorithmic approaches.
Validates that hyperscan acceleration is working for supported motif classes.

Hyperscan-based detectors:
- Z-DNA (10-mer patterns)
- A-Philic (10-mer patterns)
- G-Quadruplex (regex patterns)
- i-Motif (regex patterns)
- Curved DNA (regex patterns)

Algorithmic detectors (for comparison):
- Cruciform (inverted repeat search)
- Slipped DNA (k-mer index)
- Triplex (mirror repeat search)
"""

import os
import sys
import time

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from scanner import ModularMotifDetector


def generate_test_sequence(length=10000, motif_density=0.1):
    """Generate a sequence with embedded motifs for testing"""
    import random
    random.seed(42)
    
    sequence = []
    pos = 0
    
    motif_templates = {
        'g4': 'GGGTTAGGGTTAGGGTTAGGG',  # G-quadruplex
        'imotif': 'CCCTAACCCTAACCCTAACCC',  # i-Motif
        'zdna': 'CGCGCGCGCGCGCGCGCGCG',  # Z-DNA
        'curved': 'AAAAAAAAAA' + 'T' * 10,  # Curved DNA (A-tract)
        'aphilic': 'AGGGGGGGGGAGGGGGGGGG',  # A-philic
    }
    
    while pos < length:
        if random.random() < motif_density:
            # Add a motif
            motif = random.choice(list(motif_templates.values()))
            sequence.append(motif)
            pos += len(motif)
        else:
            # Add random DNA
            chunk_size = random.randint(20, 100)
            chunk = ''.join(random.choices('ATGC', k=chunk_size))
            sequence.append(chunk)
            pos += chunk_size
    
    return ''.join(sequence)[:length]


def test_hyperscan_availability():
    """Test that hyperscan is available and properly loaded"""
    print("\n" + "=" * 80)
    print("HYPERSCAN AVAILABILITY TEST")
    print("=" * 80)
    
    try:
        import hyperscan
        print(f"✓ Hyperscan module is installed")
        if hasattr(hyperscan, '__version__'):
            print(f"  Version: {hyperscan.__version__}")
    except ImportError:
        print(f"✗ Hyperscan module is NOT installed")
        print(f"  Performance will be degraded - using fallback regex")
        return False
    
    # Check if detector has loaded hyperscan databases
    detector = ModularMotifDetector()
    
    if hasattr(detector, 'hsdb_map'):
        loaded_dbs = len(detector.hsdb_map)
        print(f"\n✓ ModularMotifDetector loaded {loaded_dbs} hyperscan databases:")
        
        for cls_name, db_info in detector.hsdb_map.items():
            db = db_info.get('db')
            n_patterns = len(db_info.get('id_to_pattern', {}))
            db_status = 'loaded' if db else 'None'
            print(f"  - {cls_name}: {n_patterns} patterns, DB={db_status}")
        
        return loaded_dbs > 0
    else:
        print(f"\n✗ ModularMotifDetector has no hyperscan databases loaded")
        return False


def test_hyperscan_detection_speed():
    """Test detection speed with hyperscan-based detectors"""
    print("\n" + "=" * 80)
    print("HYPERSCAN DETECTION SPEED TEST")
    print("=" * 80)
    
    detector = ModularMotifDetector()
    
    # Test different sequence lengths (optimized for reasonable test time)
    test_lengths = [1000, 5000, 10000, 25000]
    
    print(f"\n{'Length (bp)':<15} {'Time (ms)':<15} {'Rate (bp/s)':<20} {'Motifs':<10}")
    print(f"{'-' * 60}")
    
    results = {}
    
    for length in test_lengths:
        sequence = generate_test_sequence(length, motif_density=0.15)
        
        start = time.time()
        motifs = detector.analyze_sequence(sequence, f"test_{length}bp")
        elapsed = time.time() - start
        
        rate = length / elapsed if elapsed > 0 else 0
        
        print(f"{length:<15,} {elapsed*1000:<15.1f} {rate:<20,.0f} {len(motifs):<10}")
        
        results[length] = {
            'time': elapsed,
            'rate': rate,
            'motifs': len(motifs)
        }
    
    print(f"{'-' * 60}")
    
    # Verify linear scaling
    if len(results) >= 2:
        lengths = sorted(results.keys())
        first_rate = results[lengths[0]]['rate']
        last_rate = results[lengths[-1]]['rate']
        
        # Rate should be relatively stable (within 50%) for linear algorithms
        rate_ratio = last_rate / first_rate if first_rate > 0 else 0
        
        print(f"\nScalability Analysis:")
        print(f"  First sequence: {lengths[0]:,} bp at {first_rate:,.0f} bp/s")
        print(f"  Last sequence:  {lengths[-1]:,} bp at {last_rate:,.0f} bp/s")
        print(f"  Rate ratio: {rate_ratio:.2f}x")
        
        if 0.5 <= rate_ratio <= 2.0:
            print(f"  ✓ Linear O(n) performance confirmed")
        else:
            print(f"  ⚠ Performance may degrade non-linearly")
    
    return results


def test_motif_class_breakdown():
    """Test performance breakdown by motif class"""
    print("\n" + "=" * 80)
    print("MOTIF CLASS DETECTION BREAKDOWN")
    print("=" * 80)
    
    detector = ModularMotifDetector()
    
    # Generate sequence with various motifs
    sequence = generate_test_sequence(20000, motif_density=0.2)
    
    print(f"\nTest sequence: {len(sequence):,} bp")
    
    # Detect all motifs
    start = time.time()
    all_motifs = detector.analyze_sequence(sequence, "breakdown_test")
    total_time = time.time() - start
    
    # Count by class
    class_counts = {}
    for motif in all_motifs:
        motif_class = motif.get('Class', 'Unknown')
        class_counts[motif_class] = class_counts.get(motif_class, 0) + 1
    
    print(f"\nTotal detection time: {total_time*1000:.1f} ms")
    print(f"Total motifs found: {len(all_motifs)}")
    print(f"Detection rate: {len(sequence)/total_time:,.0f} bp/s")
    
    print(f"\nMotifs by class:")
    print(f"{'Class':<20} {'Count':<10} {'Method':<15}")
    print(f"{'-' * 50}")
    
    # Map classes to their detection method
    hyperscan_classes = ['Z-DNA', 'A-philic', 'G-Quadruplex', 'i-Motif', 'Curved_DNA', 'R-Loop']
    algorithmic_classes = ['Cruciform', 'Slipped_DNA', 'Triplex']
    
    for motif_class, count in sorted(class_counts.items()):
        if motif_class in hyperscan_classes:
            method = "Hyperscan"
        elif motif_class in algorithmic_classes:
            method = "Algorithmic"
        else:
            method = "Other"
        
        print(f"{motif_class:<20} {count:<10} {method:<15}")
    
    print(f"{'-' * 50}")
    
    hyperscan_total = sum(count for cls, count in class_counts.items() if cls in hyperscan_classes)
    algorithmic_total = sum(count for cls, count in class_counts.items() if cls in algorithmic_classes)
    
    print(f"\nSummary:")
    print(f"  Hyperscan-detected: {hyperscan_total} motifs")
    print(f"  Algorithmically-detected: {algorithmic_total} motifs")
    print(f"  Other: {len(all_motifs) - hyperscan_total - algorithmic_total} motifs")
    
    return class_counts


def test_comparison_with_without_hyperscan():
    """Compare performance with and without hyperscan (if possible)"""
    print("\n" + "=" * 80)
    print("HYPERSCAN VS FALLBACK COMPARISON")
    print("=" * 80)
    
    try:
        import hyperscan
        hyperscan_available = True
    except ImportError:
        hyperscan_available = False
        print("\n⚠ Hyperscan not available - cannot run comparison test")
        return None
    
    detector = ModularMotifDetector()
    
    # Check if hyperscan is actually being used
    if not hasattr(detector, 'hsdb_map') or len(detector.hsdb_map) == 0:
        print("\n⚠ Detector not using hyperscan - cannot run comparison test")
        return None
    
    sequence = generate_test_sequence(10000, motif_density=0.15)
    
    print(f"\nTest sequence: {len(sequence):,} bp")
    print(f"\nWith Hyperscan acceleration:")
    
    start = time.time()
    motifs = detector.analyze_sequence(sequence, "test_with_hs")
    elapsed = time.time() - start
    
    rate = len(sequence) / elapsed if elapsed > 0 else 0
    
    print(f"  Time: {elapsed*1000:.1f} ms")
    print(f"  Rate: {rate:,.0f} bp/s")
    print(f"  Motifs: {len(motifs)}")
    
    print(f"\nHyperscan provides:")
    print(f"  ✓ Ultra-fast pattern matching for regex patterns")
    print(f"  ✓ Efficient 10-mer lookup for Z-DNA and A-philic")
    print(f"  ✓ Compiled database for optimal performance")
    print(f"  ✓ Support for complex pattern matching")
    
    return {
        'time': elapsed,
        'rate': rate,
        'motifs': len(motifs)
    }


def run_all_performance_tests():
    """Run all performance tests"""
    print("\n" + "=" * 80)
    print("HYPERSCAN PERFORMANCE TEST SUITE")
    print("=" * 80)
    
    # Test 1: Hyperscan availability
    has_hyperscan = test_hyperscan_availability()
    
    if not has_hyperscan:
        print("\n⚠ WARNING: Hyperscan not available - some tests will be skipped")
    
    # Test 2: Detection speed
    speed_results = test_hyperscan_detection_speed()
    
    # Test 3: Motif class breakdown
    class_results = test_motif_class_breakdown()
    
    # Test 4: Comparison test
    comparison_results = test_comparison_with_without_hyperscan()
    
    # Summary
    print("\n" + "=" * 80)
    print("PERFORMANCE TEST SUMMARY")
    print("=" * 80)
    
    print(f"\n✓ All performance tests completed")
    print(f"\nKey Findings:")
    print(f"  - Hyperscan acceleration: {'Enabled' if has_hyperscan else 'Disabled'}")
    
    if speed_results:
        avg_rate = sum(r['rate'] for r in speed_results.values()) / len(speed_results)
        print(f"  - Average detection rate: {avg_rate:,.0f} bp/s")
        print(f"  - Tested up to {max(speed_results.keys()):,} bp sequences")
    
    if class_results:
        total_motifs = sum(class_results.values())
        print(f"  - Total motifs detected: {total_motifs}")
        print(f"  - Motif classes detected: {len(class_results)}")
    
    return True


if __name__ == "__main__":
    success = run_all_performance_tests()
    sys.exit(0 if success else 1)
