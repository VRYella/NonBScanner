#!/usr/bin/env python3
"""
Comprehensive test suite for the updated regex registry.
Tests pattern validity, scientific accuracy, and Hyperscan compatibility.
"""

import re
import sys
import time
from collections import defaultdict

# Add the project root to path
sys.path.append('.')

from core.regex_registry import (
    ALL_PATTERNS, 
    get_patterns_for_motif, 
    get_all_hyperscan_patterns,
    validate_patterns,
    get_pattern_info,
    G_QUADRUPLEX_PATTERNS,
    I_MOTIF_PATTERNS,
    Z_DNA_PATTERNS,
    CURVED_DNA_PATTERNS,
    TRIPLEX_PATTERNS,
    CRUCIFORM_PATTERNS,
    R_LOOP_PATTERNS,
    SLIPPED_DNA_PATTERNS
)

def test_pattern_validity():
    """Test that all regex patterns are syntactically valid."""
    print("Testing regex pattern validity...")
    invalid_patterns = []
    
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        for pattern_type, patterns in pattern_dict.items():
            for i, pattern_tuple in enumerate(patterns):
                if pattern_tuple:
                    try:
                        regex_pattern = pattern_tuple[0]
                        re.compile(regex_pattern)
                    except re.error as e:
                        invalid_patterns.append((motif_class, pattern_type, i, regex_pattern, str(e)))
    
    if invalid_patterns:
        print(f"❌ Found {len(invalid_patterns)} invalid patterns:")
        for motif, ptype, idx, pattern, error in invalid_patterns:
            print(f"   {motif}.{ptype}[{idx}]: {pattern[:50]}... - {error}")
        return False
    else:
        print("✅ All patterns are syntactically valid")
        return True

def test_pattern_completeness():
    """Test that all pattern tuples have the required 9 elements."""
    print("\nTesting pattern tuple completeness...")
    incomplete_patterns = []
    
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        for pattern_type, patterns in pattern_dict.items():
            for i, pattern_tuple in enumerate(patterns):
                if pattern_tuple and len(pattern_tuple) != 9:
                    incomplete_patterns.append((motif_class, pattern_type, i, len(pattern_tuple)))
    
    if incomplete_patterns:
        print(f"❌ Found {len(incomplete_patterns)} incomplete pattern tuples:")
        for motif, ptype, idx, length in incomplete_patterns:
            print(f"   {motif}.{ptype}[{idx}]: has {length} elements (expected 9)")
        return False
    else:
        print("✅ All pattern tuples are complete (9 elements each)")
        return True

def test_pattern_ids_unique():
    """Test that all pattern IDs are unique across the registry."""
    print("\nTesting pattern ID uniqueness...")
    id_map = defaultdict(list)
    
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        for pattern_type, patterns in pattern_dict.items():
            for i, pattern_tuple in enumerate(patterns):
                if pattern_tuple:
                    pattern_id = pattern_tuple[1]
                    id_map[pattern_id].append((motif_class, pattern_type, i))
    
    duplicates = {pid: locations for pid, locations in id_map.items() if len(locations) > 1}
    
    if duplicates:
        print(f"❌ Found {len(duplicates)} duplicate pattern IDs:")
        for pid, locations in duplicates.items():
            print(f"   ID {pid}: {locations}")
        return False
    else:
        print(f"✅ All pattern IDs are unique ({len(id_map)} total)")
        return True

def test_scientific_patterns():
    """Test some key scientific patterns against known sequences."""
    print("\nTesting scientific pattern accuracy...")
    
    test_cases = [
        # G-quadruplex tests
        {
            'motif': 'g_quadruplex',
            'sequence': 'GGGTTAGGGTTAGGGTTAGGG',  # Canonical G4
            'expected_matches': True,
            'description': 'Human telomeric G4 sequence'
        },
        {
            'motif': 'g_quadruplex', 
            'sequence': 'GGGAAAGGGAAAGGGAAAGGG',  # Relaxed G4
            'expected_matches': True,
            'description': 'Relaxed G4 with A-rich loops'
        },
        # I-motif tests
        {
            'motif': 'i_motif',
            'sequence': 'CCCTAACCCTAACCCTAACCC',  # Canonical i-motif
            'expected_matches': True,
            'description': 'Human telomeric i-motif sequence'
        },
        # Z-DNA tests
        {
            'motif': 'z_dna',
            'sequence': 'CGCGCGCGCGCGCGCGCGCG',  # CG alternating
            'expected_matches': True,
            'description': 'Classical Z-DNA sequence'
        },
        # Curved DNA tests
        {
            'motif': 'curved_dna',
            'sequence': 'AAAAAAAATTTTTTT',  # A-tract + T-tract
            'expected_matches': True,
            'description': 'A-tract and T-tract'
        }
    ]
    
    failed_tests = []
    
    for test in test_cases:
        motif_patterns = get_patterns_for_motif(test['motif'])
        found_match = False
        
        for pattern_type, patterns in motif_patterns.items():
            for pattern_tuple in patterns:
                if pattern_tuple:
                    regex_pattern = pattern_tuple[0]
                    if re.search(regex_pattern, test['sequence'], re.IGNORECASE):
                        found_match = True
                        break
            if found_match:
                break
        
        if found_match != test['expected_matches']:
            failed_tests.append(test)
    
    if failed_tests:
        print(f"❌ {len(failed_tests)} scientific pattern tests failed:")
        for test in failed_tests:
            print(f"   {test['description']}: Expected {test['expected_matches']}, got {not test['expected_matches']}")
        return False
    else:
        print(f"✅ All {len(test_cases)} scientific pattern tests passed")
        return True

def test_hyperscan_compatibility():
    """Test basic Hyperscan compatibility requirements."""
    print("\nTesting Hyperscan compatibility...")
    
    try:
        import hyperscan
        hyperscan_available = True
    except ImportError:
        print("⚠️  Hyperscan not available, testing basic compatibility only")
        hyperscan_available = False
    
    incompatible_patterns = []
    
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        for pattern_type, patterns in pattern_dict.items():
            for i, pattern_tuple in enumerate(patterns):
                if pattern_tuple:
                    regex_pattern = pattern_tuple[0]
                    
                    # Check for common Hyperscan incompatibilities
                    if '\\b' in regex_pattern or '\\B' in regex_pattern:
                        incompatible_patterns.append((motif_class, pattern_type, i, "word boundaries"))
                    elif '(?=' in regex_pattern or '(?!' in regex_pattern:
                        incompatible_patterns.append((motif_class, pattern_type, i, "lookahead assertions"))
                    elif '(?<=' in regex_pattern or '(?<!' in regex_pattern:
                        incompatible_patterns.append((motif_class, pattern_type, i, "lookbehind assertions"))
                    elif '\\1' in regex_pattern or '\\2' in regex_pattern:
                        incompatible_patterns.append((motif_class, pattern_type, i, "back-references"))
    
    if incompatible_patterns:
        print(f"⚠️  Found {len(incompatible_patterns)} potentially incompatible patterns:")
        for motif, ptype, idx, issue in incompatible_patterns:
            print(f"   {motif}.{ptype}[{idx}]: {issue}")
        return False
    else:
        print("✅ All patterns appear Hyperscan-compatible")
        return True

def test_performance_benchmark():
    """Basic performance benchmark for pattern compilation."""
    print("\nRunning performance benchmark...")
    
    all_patterns = get_all_hyperscan_patterns()
    
    # Test pattern compilation time
    start_time = time.time()
    for pattern, _ in all_patterns[:50]:  # Test first 50 patterns
        re.compile(pattern)
    compile_time = time.time() - start_time
    
    print(f"✅ Compiled {min(50, len(all_patterns))} patterns in {compile_time:.3f} seconds")
    print(f"   Average: {compile_time/min(50, len(all_patterns))*1000:.2f} ms per pattern")
    
    return True

def print_registry_summary():
    """Print a summary of the regex registry contents."""
    print("\n" + "="*60)
    print("REGEX REGISTRY SUMMARY")
    print("="*60)
    
    total_patterns = 0
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        motif_total = sum(len(patterns) for patterns in pattern_dict.values())
        total_patterns += motif_total
        print(f"{motif_class.upper():<20} {motif_total:>3} patterns")
        
        for pattern_type, patterns in pattern_dict.items():
            print(f"  ├─ {pattern_type:<25} {len(patterns):>2} patterns")
    
    print("-" * 60)
    print(f"{'TOTAL':<20} {total_patterns:>3} patterns")
    
    # Scientific references count
    scoring_methods = set()
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        for pattern_type, patterns in pattern_dict.items():
            for pattern_tuple in patterns:
                if pattern_tuple and len(pattern_tuple) >= 9:
                    scoring_methods.add(pattern_tuple[8])
    
    print(f"{'SCORING METHODS':<20} {len(scoring_methods):>3} different")
    print("="*60)

def main():
    """Run all tests and print results."""
    print("REGEX REGISTRY VALIDATION SUITE")
    print("="*50)
    
    tests = [
        test_pattern_validity,
        test_pattern_completeness, 
        test_pattern_ids_unique,
        test_scientific_patterns,
        test_hyperscan_compatibility,
        test_performance_benchmark
    ]
    
    results = []
    for test in tests:
        try:
            result = test()
            results.append(result)
        except Exception as e:
            print(f"❌ Test {test.__name__} failed with exception: {e}")
            results.append(False)
    
    print("\n" + "="*50)
    print("TEST RESULTS SUMMARY")
    print("="*50)
    
    passed = sum(results)
    total = len(results)
    
    for i, (test, result) in enumerate(zip(tests, results)):
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{test.__name__:<30} {status}")
    
    print("-" * 50)
    print(f"TOTAL: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 ALL TESTS PASSED! Registry is ready for production use.")
    else:
        print("⚠️  Some tests failed. Please review the issues above.")
    
    print_registry_summary()
    
    return passed == total

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)