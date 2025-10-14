#!/usr/bin/env python3
"""
Comprehensive test script for all Non-B DNA motif classes.

This script tests each motif detector with appropriate test sequences
to ensure they are working correctly.

Uses automatic detector discovery to ensure all detectors are tested.
"""

import sys
import os
import importlib.util

# Add project root to path for imports
project_root = os.path.dirname(os.path.abspath(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# Import detector_registry directly without loading the whole utils package
spec = importlib.util.spec_from_file_location(
    "detector_registry",
    os.path.join(project_root, "utils", "detector_registry.py")
)
detector_registry = importlib.util.module_from_spec(spec)
spec.loader.exec_module(detector_registry)


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
    
    # Automatically discover all detectors and their test sequences
    detectors_with_tests = detector_registry.get_all_detectors_with_test_sequences()
    
    print(f"\nAutomatically discovered {len(detectors_with_tests)} detector classes")
    print("=" * 80)
    
    all_results = []
    
    # Test each detector with its test sequences
    for display_name, detector_class, test_sequences in detectors_with_tests:
        if not test_sequences:
            print(f"\nWarning: No test sequences defined for {display_name}")
            continue
        
        # Instantiate the detector
        detector = detector_class()
        
        # Run tests
        results = test_motif_class(display_name, detector, test_sequences)
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
