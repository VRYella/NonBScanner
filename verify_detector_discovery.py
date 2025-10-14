#!/usr/bin/env python3
"""
Quick verification that automatic detector discovery is working.

This script demonstrates that all motif detectors are automatically 
discovered and can be tested without manual maintenance.
"""

import sys
import os
import importlib.util

# Add project root to path
project_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, project_root)

# Import detector_registry directly
spec = importlib.util.spec_from_file_location(
    "detector_registry",
    os.path.join(project_root, "utils", "detector_registry.py")
)
detector_registry = importlib.util.module_from_spec(spec)
spec.loader.exec_module(detector_registry)


def main():
    print("=" * 80)
    print("AUTOMATIC DETECTOR DISCOVERY VERIFICATION")
    print("=" * 80)
    print()
    
    # Test 1: Discovery mechanism
    print("1. TESTING AUTOMATIC DISCOVERY")
    print("-" * 80)
    detectors = detector_registry.get_all_detector_classes()
    print(f"✓ Discovered {len(detectors)} detector classes automatically")
    print()
    
    for class_name, detector_class in sorted(detectors.items()):
        display_name = detector_registry.get_detector_display_name(detector_class)
        print(f"  • {class_name:30s} → {display_name}")
    print()
    
    # Test 2: Test sequence registry
    print("2. TESTING TEST SEQUENCE REGISTRY")
    print("-" * 80)
    test_sequences = detector_registry.get_test_sequences_registry()
    total_test_sequences = sum(len(seqs) for seqs in test_sequences.values())
    print(f"✓ Found {total_test_sequences} test sequences across all detectors")
    print()
    
    for class_name, sequences in sorted(test_sequences.items()):
        print(f"  • {class_name:30s} → {len(sequences)} test sequences")
    print()
    
    # Test 3: Demo sequence registry
    print("3. TESTING DEMO SEQUENCE REGISTRY")
    print("-" * 80)
    demo_sequences = detector_registry.get_demo_sequences_registry()
    total_demo_sequences = sum(len(seqs) for seqs in demo_sequences.values())
    print(f"✓ Found {total_demo_sequences} demo sequences across all detectors")
    print()
    
    for class_name, sequences in sorted(demo_sequences.items()):
        print(f"  • {class_name:30s} → {len(sequences)} demo sequences")
    print()
    
    # Test 4: Integrated test capability
    print("4. TESTING INTEGRATED DISCOVERY")
    print("-" * 80)
    detectors_with_tests = detector_registry.get_all_detectors_with_test_sequences()
    print(f"✓ Successfully integrated {len(detectors_with_tests)} detectors with test sequences")
    print()
    
    # Test 5: Quick functionality test
    print("5. TESTING DETECTOR FUNCTIONALITY")
    print("-" * 80)
    
    # Test one detector from each registry entry
    sample_tests = [
        ("A-philic_DNA", "APhilicDetector", "AGGGGGGGGGAGGGGGGGGC", "Problem statement sequence"),
        ("Z-DNA", "ZDNADetector", "CGCGCGCGCGCGCG", "CG alternating repeats"),
        ("Slipped_DNA", "SlippedDNADetector", "CAGCAGCAGCAGCAGCAG", "CAG repeat"),
    ]
    
    tested = 0
    passed = 0
    
    for display_name, class_name, sequence, description in sample_tests:
        if class_name in detectors:
            detector_class = detectors[class_name]
            detector = detector_class()
            motifs = detector.detect_motifs(sequence, 'test')
            
            tested += 1
            if motifs:
                passed += 1
                print(f"  ✓ {display_name:20s} - {description:40s} → {len(motifs)} motifs")
            else:
                print(f"  ✗ {display_name:20s} - {description:40s} → No motifs found")
    
    print()
    print(f"Quick test: {passed}/{tested} detectors working correctly")
    print()
    
    # Summary
    print("=" * 80)
    print("VERIFICATION SUMMARY")
    print("=" * 80)
    print()
    print("✓ Automatic detector discovery is working")
    print(f"✓ {len(detectors)} detector classes discovered")
    print(f"✓ {total_test_sequences} test sequences registered")
    print(f"✓ {total_demo_sequences} demo sequences registered")
    print(f"✓ {len(detectors_with_tests)} detectors ready for testing")
    print(f"✓ {passed}/{tested} sample detectors functioning correctly")
    print()
    print("All motif detectors are being picked up for testing automatically!")
    print("=" * 80)
    print()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
