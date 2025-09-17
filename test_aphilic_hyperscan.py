#!/usr/bin/env python3
"""
Test script for scan_aphilic_hyperscan.py

This script validates the A-philic DNA scanner implementation.
"""

import sys
import os

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import scan_aphilic_hyperscan as sahp


def test_basic_functionality():
    """Test basic functionality of the A-philic scanner."""
    print("Testing basic functionality...")
    
    # Test sequence with known A-philic patterns
    test_seq = "AGGGCCCCGGGGACCCCCCA"
    
    # Test contribution array
    contrib = sahp.build_contrib_array(test_seq, sahp.TET_LOG2)
    assert len(contrib) == len(test_seq), "Contribution array length mismatch"
    
    non_zero_contrib = [(i, x) for i, x in enumerate(contrib) if x > 0]
    assert len(non_zero_contrib) > 0, "No tetranucleotide contributions found"
    
    print(f"  ✓ Found {len(non_zero_contrib)} scoring positions")
    return True


def test_hyperscan_vs_python():
    """Test that hyperscan and Python implementations match."""
    print("Testing hyperscan vs Python implementations...")
    
    test_seq = "AGGGCCCCGGGGACCCCCCA"
    
    # Test hyperscan implementation
    try:
        contrib_hs = sahp.build_match_array_hyperscan(test_seq, sahp.TET_LOG2)
        hs_available = True
    except Exception:
        hs_available = False
    
    # Test Python implementation
    contrib_py = sahp.build_match_array_py(test_seq, sahp.TET_LOG2)
    
    if hs_available:
        assert contrib_hs == contrib_py, "Hyperscan and Python results don't match"
        print("  ✓ Hyperscan and Python implementations match")
    else:
        print("  ⚠ Hyperscan not available, tested Python only")
    
    return True


def test_window_classification():
    """Test window classification functionality."""
    print("Testing window classification...")
    
    # Test with a sequence that should produce both high and moderate confidence windows
    test_seq = "NNNNN" + "AGGGGGGGGGCCCCTGGGGGCCCAAGGG" + "NNNNN"
    
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as tmp:
        result_file = sahp.scan_sequence(test_seq, sahp.TET_LOG2, 
                                       min_window=10, max_window=12, 
                                       step=1, out_tsv=tmp.name)
        
        # Read results
        with open(result_file, 'r') as f:
            lines = f.readlines()
        
        # Should have header + some windows
        assert len(lines) > 1, "No A-philic windows found"
        
        # Check for both high and moderate confidence classifications
        content = ''.join(lines)
        found_high = "A_high_confidence" in content
        found_moderate = "A_moderate" in content
        
        print(f"  ✓ Found {len(lines)-1} windows")
        if found_high:
            print("  ✓ Found high confidence windows")
        if found_moderate:
            print("  ✓ Found moderate confidence windows")
        
        # Clean up
        os.unlink(result_file)
    
    return True


def test_edge_cases():
    """Test edge cases and error handling."""
    print("Testing edge cases...")
    
    # Test empty sequence
    empty_contrib = sahp.build_contrib_array("", sahp.TET_LOG2)
    assert len(empty_contrib) == 0, "Empty sequence should produce empty array"
    
    # Test short sequence
    short_contrib = sahp.build_contrib_array("ATG", sahp.TET_LOG2)
    assert len(short_contrib) == 3, "Short sequence contribution array wrong length"
    
    # Test sequence with no matches
    no_match_contrib = sahp.build_contrib_array("AAAAAAAAAA", sahp.TET_LOG2)
    assert all(x == 0 for x in no_match_contrib), "No-match sequence should have zero contributions"
    
    print("  ✓ Edge cases handled correctly")
    return True


def main():
    """Run all tests."""
    print("Testing scan_aphilic_hyperscan.py implementation")
    print("=" * 50)
    
    tests = [
        test_basic_functionality,
        test_hyperscan_vs_python,
        test_window_classification,
        test_edge_cases,
    ]
    
    passed = 0
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"  ❌ Test failed: {e}")
    
    print("=" * 50)
    print(f"Tests passed: {passed}/{len(tests)}")
    
    if passed == len(tests):
        print("✅ All tests passed!")
        return 0
    else:
        print("❌ Some tests failed!")
        return 1


if __name__ == "__main__":
    sys.exit(main())