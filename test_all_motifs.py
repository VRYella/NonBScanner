#!/usr/bin/env python3
"""
Comprehensive test script for all Non-B DNA motif classes.

This script tests each motif detector with appropriate test sequences
to ensure they are working correctly.
"""

from motif_detection.a_philic_detector import APhilicDetector
from motif_detection.z_dna_detector import ZDNADetector
from motif_detection.i_motif_detector import IMotifDetector
from motif_detection.g_quadruplex_detector import GQuadruplexDetector
from motif_detection.triplex_detector import TriplexDetector
from motif_detection.cruciform_detector import CruciformDetector
from motif_detection.slipped_dna_detector import SlippedDNADetector
from motif_detection.curved_dna_detector import CurvedDNADetector
from motif_detection.r_loop_detector import RLoopDetector


def test_motif_class(name, detector, test_sequences):
    """
    Test a motif detector with multiple sequences.
    
    Args:
        name: Name of the motif class
        detector: Instance of the detector
        test_sequences: List of (description, sequence) tuples
    
    Returns:
        Dict with test results
    """
    print(f"\n{'=' * 80}")
    print(f"Testing {name}")
    print(f"{'=' * 80}")
    
    results = {
        'class': name,
        'total_tests': len(test_sequences),
        'passed': 0,
        'failed': 0,
        'details': []
    }
    
    for desc, seq in test_sequences:
        try:
            motifs = detector.detect_motifs(seq, 'test')
            status = '✓ PASS' if motifs else '✗ FAIL'
            
            if motifs:
                results['passed'] += 1
            else:
                results['failed'] += 1
            
            print(f"{status:10s} {desc:40s} Found: {len(motifs):2d} motifs")
            
            # Show details for detected motifs
            for motif in motifs[:3]:  # Show first 3
                subclass = motif.get('Subclass', 'N/A')
                start = motif.get('Start', 0)
                end = motif.get('End', 0)
                score = motif.get('Score', 0)
                seq_fragment = motif.get('Sequence', '')[:40]
                print(f"           -> {subclass:30s} [{start:3d}-{end:3d}] Score: {score:.3f}")
                print(f"              Sequence: {seq_fragment}...")
            
            results['details'].append({
                'description': desc,
                'sequence': seq,
                'found': len(motifs),
                'passed': len(motifs) > 0
            })
            
        except Exception as e:
            results['failed'] += 1
            print(f"✗ ERROR   {desc:40s} {str(e)[:40]}")
            results['details'].append({
                'description': desc,
                'sequence': seq,
                'error': str(e),
                'passed': False
            })
    
    return results


def main():
    """Run comprehensive tests for all motif classes."""
    
    print("\n" + "=" * 80)
    print("COMPREHENSIVE NON-B DNA MOTIF DETECTOR TEST SUITE")
    print("=" * 80)
    
    all_results = []
    
    # 1. A-philic DNA Tests
    a_philic_tests = [
        ("Problem statement sequence", "AGGGGGGGGGAGGGGGGGGC"),
        ("Multiple G-runs", "AGGGGGGGGGTGGGGGGGGC"),
        ("C-rich (complement)", "CCCCCCCCCCCCCCCCCCCC"),
        ("Mixed GC-rich", "GGGGGCCCCCGGGGGCCCCC"),
    ]
    results = test_motif_class("A-philic DNA", APhilicDetector(), a_philic_tests)
    all_results.append(results)
    
    # 2. Z-DNA Tests
    z_dna_tests = [
        ("CG alternating repeats", "CGCGCGCGCGCGCG"),
        ("Longer CG repeat", "CGCGCGCGCGCGCGCGCGCG"),
        ("Embedded CG repeat", "ATCGCGCGCGCGCGAT"),
    ]
    results = test_motif_class("Z-DNA", ZDNADetector(), z_dna_tests)
    all_results.append(results)
    
    # 3. i-Motif Tests
    i_motif_tests = [
        ("C-rich telomeric", "CCCTAACCCTAACCCTAACCCT"),
        ("C-tract region", "CCCCTTCCCCTTCCCC"),
        ("Multiple C-runs", "CCCAACCCAACCCAACCC"),
    ]
    results = test_motif_class("i-Motif", IMotifDetector(), i_motif_tests)
    all_results.append(results)
    
    # 4. G-Quadruplex Tests
    g4_tests = [
        ("G4 telomeric", "GGGTTAGGGTTAGGGTTAGGG"),
        ("Canonical G4", "GGGGTTTTGGGGTTTTGGGG"),
        ("Long loop G4", "GGGGAAAAAAAGGGGAAAAAAAGGGG"),
    ]
    results = test_motif_class("G-Quadruplex", GQuadruplexDetector(), g4_tests)
    all_results.append(results)
    
    # 5. Triplex Tests
    triplex_tests = [
        ("Polypurine-Polypyrimidine", "AGAGAGAGAGAGAGAGAGAGAGAGAGAG"),
        ("GA repeat", "GAGAGAGAGAGAGAGAGAGAGA"),
        ("Purine tract", "AAAAGGGGAAAAGGGGAAAA"),
    ]
    results = test_motif_class("Triplex", TriplexDetector(), triplex_tests)
    all_results.append(results)
    
    # 6. Cruciform Tests
    cruciform_tests = [
        ("Inverted repeat", "AAAATTTTAAAATTTT"),
        ("Palindrome", "ATCGATCGATCGATCG"),
        ("Long inverted repeat", "AAAAATTTTGGGGGCCCCCCCCCCCCCCGGGGG"),
    ]
    results = test_motif_class("Cruciform", CruciformDetector(), cruciform_tests)
    all_results.append(results)
    
    # 7. Slipped DNA Tests
    slipped_tests = [
        ("CAG repeat (Huntington's)", "CAGCAGCAGCAGCAGCAG"),
        ("CGG repeat (Fragile X)", "CGGCGGCGGCGGCGGCGG"),
        ("GAA repeat", "GAAGAAGAAGAAGAAGAA"),
    ]
    results = test_motif_class("Slipped DNA", SlippedDNADetector(), slipped_tests)
    all_results.append(results)
    
    # 8. Curved DNA Tests
    curved_tests = [
        ("A-tract with phasing", "AAAAAAAATAAAAAAA"),
        ("Multiple A-tracts", "AAAAATTTTTAAAAATTTTTAAAAA"),
        ("AT-rich region", "AAAAAATTTTTTAAAAAATTTTTT"),
    ]
    results = test_motif_class("Curved DNA", CurvedDNADetector(), curved_tests)
    all_results.append(results)
    
    # 9. R-Loop Tests
    r_loop_tests = [
        ("G-rich skew", "GGGGGAAAAAGGGGGAAAAAGGGGGAAAAA"),
        ("Long G-skew", "GGGGGGGGAAAAAAAAAAGGGGGGGGAAAAAAAAAA"),
        ("Multiple G clusters", "GGGGGTTGGGGGTTGGGGG"),
    ]
    results = test_motif_class("R-Loop", RLoopDetector(), r_loop_tests)
    all_results.append(results)
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    total_classes = len(all_results)
    total_tests = sum(r['total_tests'] for r in all_results)
    total_passed = sum(r['passed'] for r in all_results)
    total_failed = sum(r['failed'] for r in all_results)
    
    for result in all_results:
        status = "✓" if result['failed'] == 0 else "✗"
        print(f"{status} {result['class']:20s} "
              f"Passed: {result['passed']:2d}/{result['total_tests']:2d} "
              f"Failed: {result['failed']:2d}")
    
    print(f"\n{'=' * 80}")
    print(f"Total Classes: {total_classes}")
    print(f"Total Tests:   {total_tests}")
    print(f"Total Passed:  {total_passed}")
    print(f"Total Failed:  {total_failed}")
    print(f"Success Rate:  {100 * total_passed / total_tests:.1f}%")
    print(f"{'=' * 80}\n")
    
    return all_results


if __name__ == "__main__":
    results = main()
