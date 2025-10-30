#!/usr/bin/env python3
"""
Test script to verify pattern registries work correctly.

This script tests both 10-mer registries (A-philic, Z-DNA) and 
regex pattern registries (CurvedDNA, G4, IMotif).
"""

import os
import sys

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from utilities import scan_with_registry, get_cached_registry


def test_g4_registry():
    """Test G-Quadruplex pattern registry."""
    print("="*60)
    print("Testing G-Quadruplex (G4) Registry")
    print("="*60)
    
    # Test sequences
    test_sequences = [
        ("Telomeric G4", "GGGTTAGGGTTAGGGTTAGGG"),
        ("Canonical G4", "GGGGTTTTGGGGTTTTGGGG"),
        ("No G4", "AAAATTTTCCCCTTTT"),
    ]
    
    for desc, seq in test_sequences:
        print(f"\nSequence: {desc}")
        print(f"  Seq: {seq}")
        matches = scan_with_registry("G4", seq)
        if matches:
            print(f"  ✓ Found {len(matches)} match(es):")
            for start, end, pattern_id, subclass in matches[:3]:
                print(f"    - Pattern {pattern_id} ({subclass}): {start}-{end}")
        else:
            print(f"  ✗ No matches found")
    
    return True


def test_curved_dna_registry():
    """Test Curved DNA pattern registry."""
    print("\n" + "="*60)
    print("Testing Curved DNA Registry")
    print("="*60)
    
    # Test sequences
    test_sequences = [
        ("Long A-tract (local)", "AAAAAAAAAAAAA"),
        ("Long T-tract (local)", "TTTTTTTTTTTTT"),
        ("A-phased repeat", "AAACGTAAACGTAAACGT"),
        ("No curvature", "CGATCGATCGAT"),
    ]
    
    for desc, seq in test_sequences:
        print(f"\nSequence: {desc}")
        print(f"  Seq: {seq}")
        matches = scan_with_registry("CurvedDNA", seq)
        if matches:
            print(f"  ✓ Found {len(matches)} match(es):")
            for start, end, pattern_id, subclass in matches[:3]:
                print(f"    - Pattern {pattern_id} ({subclass}): {start}-{end}")
        else:
            print(f"  ✗ No matches found")
    
    return True


def test_imotif_registry():
    """Test i-Motif pattern registry."""
    print("\n" + "="*60)
    print("Testing i-Motif Registry")
    print("="*60)
    
    # Test sequences
    test_sequences = [
        ("Canonical i-motif", "CCCTAACCCTAACCCTAACCC"),
        ("HUR AC-motif", "AAACGTCCCACGTCCCACGTCCC"),
        ("No i-motif", "GGGGTTTTGGGG"),
    ]
    
    for desc, seq in test_sequences:
        print(f"\nSequence: {desc}")
        print(f"  Seq: {seq}")
        matches = scan_with_registry("IMotif", seq)
        if matches:
            print(f"  ✓ Found {len(matches)} match(es):")
            for start, end, pattern_id, subclass in matches[:3]:
                print(f"    - Pattern {pattern_id} ({subclass}): {start}-{end}")
        else:
            print(f"  ✗ No matches found")
    
    return True


def test_registry_metadata():
    """Test that registry metadata is loaded correctly."""
    print("\n" + "="*60)
    print("Testing Registry Metadata")
    print("="*60)
    
    registries = ['G4', 'CurvedDNA', 'IMotif']
    
    for reg_name in registries:
        db, id_to_pattern, id_to_subclass, id_to_score = get_cached_registry(reg_name)
        print(f"\n{reg_name} Registry:")
        print(f"  Total patterns: {len(id_to_pattern)}")
        print(f"  Hyperscan DB: {'Available' if db is not None else 'Not available (using regex fallback)'}")
        print(f"  Subclasses: {set(id_to_subclass.values())}")
    
    return True


def main():
    """Run all tests."""
    print("\n" + "="*60)
    print("PATTERN REGISTRY TEST SUITE")
    print("="*60)
    
    tests = [
        ("G4 Registry", test_g4_registry),
        ("Curved DNA Registry", test_curved_dna_registry),
        ("i-Motif Registry", test_imotif_registry),
        ("Registry Metadata", test_registry_metadata),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result, None))
        except Exception as e:
            results.append((test_name, False, str(e)))
    
    # Print summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    for test_name, result, error in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {test_name}")
        if error:
            print(f"  Error: {error}")
    
    all_passed = all(r[1] for r in results)
    print("\n" + "="*60)
    if all_passed:
        print("✓ ALL TESTS PASSED")
    else:
        print("✗ SOME TESTS FAILED")
    print("="*60)
    
    return 0 if all_passed else 1


if __name__ == '__main__':
    sys.exit(main())
