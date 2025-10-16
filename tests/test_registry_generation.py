#!/usr/bin/env python3
"""
Test registry generation functionality.

Tests:
- Registry generation from detector classes
- Pattern validation
- File format correctness (txt, pkl, json)
- Deterministic ID assignment
"""
import os
import sys
import tempfile
import shutil
import pickle
import json

# Ensure we can import from the repo
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tools.generate_all_registries import (
    discover_detector_modules,
    normalize_pattern_table,
    validate_pattern,
    save_pattern_text,
    save_registry,
    generate_class_registry
)


def test_discovery():
    """Test that detector discovery works"""
    print("TEST: Detector discovery")
    print("-" * 60)
    
    discovered = discover_detector_modules("motif_detection")
    
    assert len(discovered) >= 2, "Should discover at least ZDNA and APhilic"
    assert "ZDNA" in discovered or "Zdna" in discovered, "Should discover ZDNA"
    assert "APhilic" in discovered or "Aphilic" in discovered, "Should discover APhilic"
    
    for class_name, (module, detector_cls, table, attr_name) in discovered.items():
        assert isinstance(table, dict), f"{class_name} table should be dict"
        assert len(table) > 0, f"{class_name} table should not be empty"
        print(f"  ✓ {class_name}: {len(table)} patterns from {attr_name}")
    
    print("  ✓ Discovery test passed")
    print()
    return discovered


def test_pattern_validation():
    """Test pattern validation logic"""
    print("TEST: Pattern validation")
    print("-" * 60)
    
    # Valid patterns
    valid_patterns = ["ACGT", "AAAA", "CGCG", "ATCGATCG", "NNNACGT"]
    for p in valid_patterns:
        is_valid, warning = validate_pattern(p)
        assert is_valid, f"Pattern '{p}' should be valid"
        print(f"  ✓ '{p}' is valid")
    
    # Invalid patterns
    invalid_patterns = ["ACGTZ", "12345", "ACGT!", "ACGT-"]
    for p in invalid_patterns:
        is_valid, warning = validate_pattern(p)
        assert not is_valid, f"Pattern '{p}' should be invalid"
        print(f"  ✓ '{p}' is invalid (as expected)")
    
    print("  ✓ Validation test passed")
    print()


def test_normalization():
    """Test pattern table normalization"""
    print("TEST: Pattern normalization")
    print("-" * 60)
    
    # Test table with various formats
    test_table = {
        "acgt": 1.5,
        "CGTA": 2.0,
        "gggg": 3.5,
        "TTTT": 4.0,
    }
    
    items = normalize_pattern_table(test_table)
    
    # Check that patterns are sorted
    patterns = [p for (i, p, s) in items]
    assert patterns == sorted(patterns), "Patterns should be sorted"
    print(f"  ✓ Patterns are sorted: {patterns}")
    
    # Check that IDs are sequential
    ids = [i for (i, p, s) in items]
    assert ids == list(range(len(ids))), "IDs should be sequential"
    print(f"  ✓ IDs are sequential: {ids}")
    
    # Check that patterns are uppercased
    for i, pattern, score in items:
        assert pattern == pattern.upper(), "Patterns should be uppercase"
    print(f"  ✓ Patterns are uppercased")
    
    # Check that scores are float
    for i, pattern, score in items:
        assert isinstance(score, float), "Scores should be float"
    print(f"  ✓ Scores are float")
    
    print("  ✓ Normalization test passed")
    print()


def test_registry_files():
    """Test that registry files are created correctly"""
    print("TEST: Registry file generation")
    print("-" * 60)
    
    # Create temp directory
    tmpdir = tempfile.mkdtemp()
    try:
        # Test data
        test_table = {
            "ACGT": 1.0,
            "CGTA": 2.0,
            "GGGG": 3.0,
            "TTTT": 4.0,
        }
        
        items = normalize_pattern_table(test_table)
        class_name = "TestClass"
        
        # Save pattern text
        save_pattern_text(tmpdir, class_name, items)
        txt_path = os.path.join(tmpdir, f"{class_name}_patterns.txt")
        assert os.path.exists(txt_path), "Pattern text file should exist"
        
        with open(txt_path, 'r') as f:
            lines = [line.strip() for line in f]
        assert len(lines) == len(items), "Should have one line per pattern"
        print(f"  ✓ Pattern text file created: {len(lines)} lines")
        
        # Save registry
        meta = {"test": "metadata"}
        save_registry(tmpdir, class_name, items, meta)
        
        # Check pickle
        pkl_path = os.path.join(tmpdir, f"{class_name}_registry.pkl")
        assert os.path.exists(pkl_path), "Pickle file should exist"
        with open(pkl_path, 'rb') as f:
            pkl_data = pickle.load(f)
        assert pkl_data['class'] == class_name
        assert pkl_data['n_patterns'] == len(items)
        assert len(pkl_data['patterns']) == len(items)
        print(f"  ✓ Pickle file created with {len(pkl_data['patterns'])} patterns")
        
        # Check JSON
        json_path = os.path.join(tmpdir, f"{class_name}_registry.json")
        assert os.path.exists(json_path), "JSON file should exist"
        with open(json_path, 'r') as f:
            json_data = json.load(f)
        assert json_data['class'] == class_name
        assert json_data['n_patterns'] == len(items)
        assert len(json_data['patterns']) == len(items)
        print(f"  ✓ JSON file created with {len(json_data['patterns'])} patterns")
        
        # Check pattern format (should use 'tenmer' key)
        first_pattern = pkl_data['patterns'][0]
        assert 'id' in first_pattern
        assert 'tenmer' in first_pattern, "Pattern should use 'tenmer' key"
        assert 'score' in first_pattern
        print(f"  ✓ Pattern format correct: {first_pattern}")
        
        print("  ✓ File generation test passed")
        print()
        
    finally:
        shutil.rmtree(tmpdir)


def test_end_to_end():
    """Test end-to-end registry generation"""
    print("TEST: End-to-end generation")
    print("-" * 60)
    
    tmpdir = tempfile.mkdtemp()
    try:
        # Generate a small test registry
        test_table = {
            "AAAAAAAAAA": 10.0,
            "CCCCCCCCCC": 20.0,
            "GGGGGGGGGG": 30.0,
            "TTTTTTTTTT": 40.0,
        }
        
        meta = {"source": "test", "note": "test data"}
        generate_class_registry(tmpdir, "SmallTest", test_table, meta, shard_size=50000)
        
        # Verify files exist
        txt_path = os.path.join(tmpdir, "SmallTest_patterns.txt")
        pkl_path = os.path.join(tmpdir, "SmallTest_registry.pkl")
        json_path = os.path.join(tmpdir, "SmallTest_registry.json")
        
        assert os.path.exists(txt_path), "Text file should exist"
        assert os.path.exists(pkl_path), "Pickle file should exist"
        assert os.path.exists(json_path), "JSON file should exist"
        
        # Load and verify registry can be loaded
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'utils'))
        from load_hsdb import load_db_for_class
        
        db, id_to_ten, id_to_score = load_db_for_class("SmallTest", tmpdir)
        assert len(id_to_ten) == 4, "Should have 4 patterns"
        assert len(id_to_score) == 4, "Should have 4 scores"
        print(f"  ✓ Loaded {len(id_to_ten)} patterns via load_db_for_class")
        
        # Verify pattern order is deterministic
        patterns_list = [id_to_ten[i] for i in sorted(id_to_ten.keys())]
        assert patterns_list == sorted(patterns_list), "Patterns should be sorted"
        print(f"  ✓ Patterns are deterministically ordered: {patterns_list}")
        
        print("  ✓ End-to-end test passed")
        print()
        
    finally:
        shutil.rmtree(tmpdir)


def test_loader_compatibility():
    """Test that generated registries are compatible with load_hsdb"""
    print("TEST: Loader compatibility")
    print("-" * 60)
    
    tmpdir = tempfile.mkdtemp()
    try:
        # Use real detector data
        from motif_detection.z_dna_detector import ZDNADetector
        table = ZDNADetector.TENMER_SCORE
        
        meta = {"source": "ZDNADetector.TENMER_SCORE"}
        generate_class_registry(tmpdir, "ZDNA", table, meta, shard_size=50000)
        
        # Load via load_hsdb
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'utils'))
        from load_hsdb import load_db_for_class
        
        db, id_to_ten, id_to_score = load_db_for_class("ZDNA", tmpdir)
        
        assert len(id_to_ten) == len(table), f"Should load all {len(table)} patterns"
        print(f"  ✓ Loaded {len(id_to_ten)} patterns")
        
        # Verify some patterns
        for pattern_id in [0, 10, 50]:
            if pattern_id < len(id_to_ten):
                pattern = id_to_ten[pattern_id]
                score = id_to_score[pattern_id]
                assert pattern in table, f"Pattern {pattern} should be in original table"
                assert score == table[pattern], f"Score should match original"
        print(f"  ✓ Pattern scores match original table")
        
        print("  ✓ Loader compatibility test passed")
        print()
        
    finally:
        shutil.rmtree(tmpdir)


def main():
    """Run all tests"""
    print("="*70)
    print("REGISTRY GENERATION TESTS")
    print("="*70)
    print()
    
    try:
        test_discovery()
        test_pattern_validation()
        test_normalization()
        test_registry_files()
        test_end_to_end()
        test_loader_compatibility()
        
        print("="*70)
        print("ALL TESTS PASSED ✓")
        print("="*70)
        return 0
        
    except AssertionError as e:
        print()
        print("="*70)
        print(f"TEST FAILED: {e}")
        print("="*70)
        return 1
    except Exception as e:
        print()
        print("="*70)
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        print("="*70)
        return 1


if __name__ == "__main__":
    sys.exit(main())
