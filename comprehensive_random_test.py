#!/usr/bin/env python3
"""
Comprehensive random test suite for all Non-B DNA motif classes.
Generates random sequences with specific motif patterns and verifies detection.
"""

import sys
import os
import random

# Add project root to path
project_root = os.path.dirname(os.path.abspath(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from utils.detector_registry import get_all_detectors_with_test_sequences


def generate_random_dna(length, seed=None):
    """Generate random DNA sequence of given length.
    
    Args:
        length: Length of sequence to generate
        seed: Optional seed for reproducibility
    """
    if seed is not None:
        random.seed(seed)
    return ''.join(random.choices('ATGC', k=length))


def generate_test_sequences(use_seed=True):
    """Generate comprehensive test sequences for each motif class.
    
    Args:
        use_seed: If True, use fixed seed for reproducibility
    """
    # Set seed for reproducible test sequences
    if use_seed:
        random.seed(42)
    
    test_cases = {
        'G-Quadruplex': [
            # Canonical G4 - tight loops
            ('Canonical G4 (short loops)', 'GGG' + generate_random_dna(2) + 'GGG' + generate_random_dna(2) + 'GGG' + generate_random_dna(2) + 'GGG'),
            # Canonical G4 - medium loops
            ('Canonical G4 (med loops)', 'GGGG' + generate_random_dna(5) + 'GGGG' + generate_random_dna(4) + 'GGGG' + generate_random_dna(3) + 'GGGG'),
            # Relaxed G4 - longer loops
            ('Relaxed G4 (long loops)', 'GGG' + generate_random_dna(10) + 'GGG' + generate_random_dna(9) + 'GGG' + generate_random_dna(8) + 'GGG'),
            # Telomeric pattern
            ('Telomeric G4', 'GGGTTA' * 4),
            # Multiple G runs
            ('Multi-run G4', 'GGGGG' + generate_random_dna(3) + 'GGG' + generate_random_dna(4) + 'GGGG' + generate_random_dna(2) + 'GGGGG'),
        ],
        
        'i-Motif': [
            # Canonical i-motif
            ('Canonical i-motif (short)', 'CCC' + generate_random_dna(2) + 'CCC' + generate_random_dna(2) + 'CCC' + generate_random_dna(2) + 'CCC'),
            # Telomeric C-rich
            ('Telomeric i-motif', 'CCCTAA' * 4),
            # Relaxed i-motif
            ('Relaxed i-motif', 'CCCC' + generate_random_dna(7) + 'CCC' + generate_random_dna(6) + 'CCC' + generate_random_dna(5) + 'CCCC'),
            # AC-motif
            ('AC-motif pattern', 'AAACCCCAAACCCAAACCCC'),
        ],
        
        'Curved_DNA': [
            # Phased A-tracts (global curvature)
            ('Phased A-tracts 1', 'AAAA' + generate_random_dna(11) + 'AAAA' + generate_random_dna(11) + 'AAAA'),
            ('Phased A-tracts 2', 'AAAAA' + generate_random_dna(10) + 'AAAA' + generate_random_dna(11) + 'AAAAA' + generate_random_dna(10) + 'AAAA'),
            # Long A-tract (local curvature)
            ('Long A-tract', generate_random_dna(10) + 'AAAAAAAA' + generate_random_dna(10)),
            # Long T-tract
            ('Long T-tract', generate_random_dna(10) + 'TTTTTTTT' + generate_random_dna(10)),
            # Mixed AT-rich
            ('Mixed AT pattern', 'AAAA' + 'TT' + 'AAAA' + generate_random_dna(11) + 'AAAA' + 'TT' + 'AAAA'),
        ],
        
        'Z-DNA': [
            # Perfect CG alternating
            ('CG alternating (perfect)', 'CG' * 10),
            # Longer CG repeat
            ('CG repeat (long)', 'CG' * 15),
            # Embedded CG repeat
            ('Embedded CG', generate_random_dna(5) + 'CGCGCGCGCGCG' + generate_random_dna(5)),
            # CA alternating
            ('CA alternating', 'CA' * 8),
        ],
        
        'Slipped_DNA': [
            # CAG repeat (Huntington's)
            ('CAG repeat (HD)', 'CAG' * 10),
            # CGG repeat (Fragile X)
            ('CGG repeat (FX)', 'CGG' * 10),
            # GAA repeat (Friedreich's)
            ('GAA repeat (FA)', 'GAA' * 10),
            # CTG repeat
            ('CTG repeat', 'CTG' * 8),
            # Direct repeat with spacer
            ('Direct repeat', 'ATCGATCG' + generate_random_dna(20) + 'ATCGATCG'),
        ],
        
        'Cruciform': [
            # Perfect palindrome
            ('Perfect palindrome', 'ATCGATCG' + generate_random_dna(10) + 'CGATCGAT'),
            # Longer inverted repeat
            ('Long inverted repeat', 'ATCGATCGATCG' + generate_random_dna(15) + 'CGATCGATCGAT'),
            # Short palindrome
            ('Short palindrome', 'ATCGCGAT'),
            # AT-rich palindrome  
            ('AT palindrome', 'AAATTTT' + generate_random_dna(8) + 'AAAATTTT'),
        ],
        
        'R-Loop': [
            # G-rich skew
            ('G-rich skew 1', 'GGGGGG' + generate_random_dna(5) + 'GGGGGG' + generate_random_dna(5) + 'GGGGGG'),
            # Long G-cluster
            ('Long G-cluster', 'GGGGGGGGG' + generate_random_dna(10) + 'GGGGGGGG'),
            # Multiple G-runs
            ('Multi G-runs', 'GGGGG' + 'AA' + 'GGGGG' + 'TT' + 'GGGGG' + 'CC' + 'GGGGG'),
            # GC-rich region
            ('GC-rich', 'GCGCGCGCGC' + generate_random_dna(5) + 'GCGCGCGC'),
        ],
        
        'Triplex': [
            # Polypurine-polypyrimidine
            ('Pu-Py mirror', 'AGAGAGAGAGAGAGAG' + generate_random_dna(10) + 'CTCTCTCTCTCTCTCT'),
            # GA repeat
            ('GA repeat', 'GAGAGAGAGAGAGA'),
            # AG repeat  
            ('AG repeat', 'AGAGAGAGAGAGAG'),
            # Purine tract
            ('Purine tract', 'AAAGGGAAAGGGAAA'),
            # Mirror with spacer
            ('Mirror spacer', 'GGAGGAGGAG' + generate_random_dna(15) + 'CTCCTCCTCC'),
        ],
        
        'A-philic_DNA': [
            # A-rich region
            ('A-rich region', 'AAAAAAAAAAAAAAAAAAAA'),
            # AT-rich
            ('AT-rich', 'ATATATATATATATATAT'),
            # Multiple A-tracts
            ('Multi A-tracts', 'AAAAA' + 'TT' + 'AAAAA' + 'GG' + 'AAAAA' + 'CC' + 'AAAAA'),
            # A-rich with spacing
            ('A-rich spaced', 'AAAA' + generate_random_dna(3) + 'AAAA' + generate_random_dna(3) + 'AAAA'),
        ],
    }
    
    return test_cases


def run_comprehensive_test():
    """Run comprehensive test on all detectors with random sequences."""
    
    print("\n" + "=" * 80)
    print("COMPREHENSIVE RANDOM SEQUENCE TEST SUITE")
    print("=" * 80)
    
    # Get all detectors
    detectors_with_tests = get_all_detectors_with_test_sequences()
    detector_map = {name: (cls, tests) for name, cls, tests in detectors_with_tests}
    
    # Generate test sequences
    test_sequences = generate_test_sequences()
    
    all_results = []
    
    # Create explicit mapping between motif classes and detector names
    motif_to_detector = {
        'G-Quadruplex': 'G-Quadruplex',
        'i-Motif': 'i-Motif',
        'Curved_DNA': 'Curved_DNA',
        'Z-DNA': 'Z-DNA',
        'Slipped_DNA': 'Slipped_DNA',
        'Cruciform': 'Cruciform',
        'R-Loop': 'R-Loop',
        'Triplex': 'Triplex',
        'A-philic_DNA': 'A-philic_DNA',
    }
    
    for motif_class, test_cases in test_sequences.items():
        # Find matching detector using explicit mapping
        detector_name = motif_to_detector.get(motif_class)
        
        if not detector_name or detector_name not in detector_map:
            print(f"\n⚠️  No detector found for {motif_class}")
            continue
        
        detector_cls, _ = detector_map[detector_name]
        detector = detector_cls()
        
        print(f"\n{'=' * 80}")
        print(f"Testing {motif_class}")
        print(f"{'=' * 80}")
        
        results = {
            'class': motif_class,
            'total_tests': len(test_cases),
            'passed': 0,
            'failed': 0,
            'details': []
        }
        
        for desc, seq in test_cases:
            try:
                motifs = detector.detect_motifs(seq, 'test')
                status = '✓ PASS' if motifs else '✗ FAIL'
                
                if motifs:
                    results['passed'] += 1
                else:
                    results['failed'] += 1
                
                print(f"{status:10s} {desc:40s} Found: {len(motifs):2d} motifs")
                
                # Show details for first motif
                if motifs:
                    m = motifs[0]
                    subclass = m.get('Subclass', 'N/A')
                    start = m.get('Start', 0)
                    end = m.get('End', 0)
                    score = m.get('Score', 0)
                    print(f"           -> {subclass:30s} [{start:3d}-{end:3d}] Score: {score:.3f}")
                    print(f"              Sequence: {seq[:40]}...")
                
                results['details'].append({
                    'description': desc,
                    'found': len(motifs),
                    'passed': len(motifs) > 0
                })
                
            except Exception as e:
                results['failed'] += 1
                print(f"✗ ERROR   {desc:40s} {str(e)[:40]}")
                results['details'].append({
                    'description': desc,
                    'error': str(e),
                    'passed': False
                })
        
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
        success_rate = (result['passed'] / result['total_tests'] * 100) if result['total_tests'] > 0 else 0
        print(f"{status} {result['class']:25s} "
              f"Passed: {result['passed']:2d}/{result['total_tests']:2d} "
              f"({success_rate:5.1f}%)")
    
    print(f"\n{'=' * 80}")
    print(f"Total Classes Tested: {total_classes}")
    print(f"Total Tests:          {total_tests}")
    print(f"Total Passed:         {total_passed}")
    print(f"Total Failed:         {total_failed}")
    print(f"Overall Success Rate: {100 * total_passed / total_tests:.1f}%")
    print(f"{'=' * 80}\n")
    
    return all_results, (total_passed == total_tests)


if __name__ == "__main__":
    results, all_passed = run_comprehensive_test()
    sys.exit(0 if all_passed else 1)
