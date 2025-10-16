#!/usr/bin/env python3
"""
Test Hyperscan registry loader integration.

This test validates:
1. Registry generation (tools/generate_class_hsdb.py)
2. Registry loading (utils/load_hsdb.py)
3. Integration with motif_patterns module
4. Integration with modular_scanner
"""

import os
import sys
import tempfile
import shutil

# Ensure we can import from the repo
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_registry_generation():
    """Test that registries can be generated"""
    print("=" * 70)
    print("TEST 1: Registry Generation")
    print("=" * 70)
    
    # Create temporary directory for test registry
    test_dir = tempfile.mkdtemp()
    try:
        # Generate registries
        from tools.generate_class_hsdb import main as gen_main
        
        # Temporarily override sys.argv to pass arguments
        old_argv = sys.argv
        sys.argv = ['generate_class_hsdb.py', '--out', test_dir]
        
        try:
            gen_main()
        except SystemExit:
            pass  # main() may call sys.exit()
        finally:
            sys.argv = old_argv
        
        # Check that files were created
        expected_files = [
            'ZDNA_patterns.txt',
            'ZDNA_registry.pkl',
            'ZDNA_registry.json',
            'APhilic_patterns.txt',
            'APhilic_registry.pkl',
            'APhilic_registry.json',
        ]
        
        for fname in expected_files:
            fpath = os.path.join(test_dir, fname)
            assert os.path.isfile(fpath), f"Missing file: {fname}"
            print(f"  ✓ {fname} generated")
        
        print("✓ Registry generation PASSED\n")
        return test_dir
    except Exception as e:
        print(f"✗ Registry generation FAILED: {e}\n")
        shutil.rmtree(test_dir, ignore_errors=True)
        raise


def test_registry_loading(test_dir):
    """Test that registries can be loaded"""
    print("=" * 70)
    print("TEST 2: Registry Loading")
    print("=" * 70)
    
    from utils.load_hsdb import load_db_for_class
    
    # Test Z-DNA loading
    db, id_to_ten, id_to_score = load_db_for_class('ZDNA', test_dir)
    print(f"  ✓ Loaded ZDNA: {len(id_to_ten)} patterns")
    assert len(id_to_ten) > 0, "ZDNA registry should have patterns"
    assert len(id_to_ten) == len(id_to_score), "Pattern and score counts should match"
    
    # Test A-philic loading
    db, id_to_ten, id_to_score = load_db_for_class('APhilic', test_dir)
    print(f"  ✓ Loaded APhilic: {len(id_to_ten)} patterns")
    assert len(id_to_ten) > 0, "APhilic registry should have patterns"
    
    print("✓ Registry loading PASSED\n")


def test_motif_patterns_integration(test_dir):
    """Test integration with motif_patterns module"""
    print("=" * 70)
    print("TEST 3: motif_patterns Integration")
    print("=" * 70)
    
    from utils.motif_patterns import get_hs_db_for_class, get_pattern_registry
    
    # Test get_hs_db_for_class
    db, id_to_ten, id_to_score = get_hs_db_for_class('ZDNA', test_dir)
    print(f"  ✓ get_hs_db_for_class('ZDNA'): {len(id_to_ten)} patterns")
    
    # Test get_pattern_registry
    registry = get_pattern_registry('ZDNA', test_dir)
    assert 'patterns' in registry, "Registry should have 'patterns' key"
    assert 'class' in registry, "Registry should have 'class' key"
    assert registry['class'] == 'ZDNA', "Class should be ZDNA"
    print(f"  ✓ get_pattern_registry('ZDNA'): {registry['n_patterns']} patterns")
    
    # Test caching
    db2, id_to_ten2, id_to_score2 = get_hs_db_for_class('ZDNA', test_dir)
    assert id_to_ten is id_to_ten2, "Should return cached result"
    print(f"  ✓ Caching works correctly")
    
    print("✓ motif_patterns integration PASSED\n")


def test_modular_scanner_integration(test_dir):
    """Test integration with modular_scanner"""
    print("=" * 70)
    print("TEST 4: modular_scanner Integration")
    print("=" * 70)
    
    # Set environment variable to use test registry
    old_env = os.environ.get('NBD_REGISTRY_DIR')
    os.environ['NBD_REGISTRY_DIR'] = test_dir
    
    try:
        from utils.modular_scanner import ModularMotifDetector
        
        # Initialize scanner (should preload DBs)
        scanner = ModularMotifDetector()
        print(f"  ✓ Scanner initialized")
        
        # Check that hsdb_map was created
        assert hasattr(scanner, 'hsdb_map'), "Scanner should have hsdb_map"
        print(f"  ✓ hsdb_map attribute created")
        
        # Check that expected classes are in hsdb_map
        assert 'ZDNA' in scanner.hsdb_map, "ZDNA should be in hsdb_map"
        assert 'APhilic' in scanner.hsdb_map, "APhilic should be in hsdb_map"
        print(f"  ✓ ZDNA and APhilic loaded into hsdb_map")
        
        # Check that DBs have the expected structure
        zdna_info = scanner.hsdb_map['ZDNA']
        assert 'db' in zdna_info, "ZDNA info should have 'db' key"
        assert 'id_to_ten' in zdna_info, "ZDNA info should have 'id_to_ten' key"
        assert 'id_to_score' in zdna_info, "ZDNA info should have 'id_to_score' key"
        print(f"  ✓ DB info structure is correct")
        
        # Test that scanner can still analyze sequences
        test_seq = "CGCGCGCGCGCGCGCG"
        motifs = scanner.analyze_sequence(test_seq, "test_seq")
        print(f"  ✓ Scanner can analyze sequences: found {len(motifs)} motifs")
        
        print("✓ modular_scanner integration PASSED\n")
    finally:
        # Restore environment
        if old_env is None:
            os.environ.pop('NBD_REGISTRY_DIR', None)
        else:
            os.environ['NBD_REGISTRY_DIR'] = old_env


def test_detector_class_attribute():
    """Test that detector classes have HS_DB_INFO set"""
    print("=" * 70)
    print("TEST 5: Detector Class Attribute")
    print("=" * 70)
    
    from motif_detection.z_dna_detector import ZDNADetector
    from motif_detection.a_philic_detector import APhilicDetector
    
    # Check if HS_DB_INFO was set by scanner initialization
    if hasattr(ZDNADetector, 'HS_DB_INFO'):
        print(f"  ✓ ZDNADetector.HS_DB_INFO is set")
        info = ZDNADetector.HS_DB_INFO
        assert 'id_to_ten' in info, "HS_DB_INFO should have id_to_ten"
        print(f"    - Contains {len(info['id_to_ten'])} patterns")
    else:
        print(f"  ℹ ZDNADetector.HS_DB_INFO not set (optional)")
    
    if hasattr(APhilicDetector, 'HS_DB_INFO'):
        print(f"  ✓ APhilicDetector.HS_DB_INFO is set")
        info = APhilicDetector.HS_DB_INFO
        assert 'id_to_ten' in info, "HS_DB_INFO should have id_to_ten"
        print(f"    - Contains {len(info['id_to_ten'])} patterns")
    else:
        print(f"  ℹ APhilicDetector.HS_DB_INFO not set (optional)")
    
    print("✓ Detector class attribute test PASSED\n")


def main():
    """Run all tests"""
    print("\n" + "=" * 70)
    print("HYPERSCAN REGISTRY LOADER INTEGRATION TEST SUITE")
    print("=" * 70 + "\n")
    
    test_dir = None
    try:
        # Run tests
        test_dir = test_registry_generation()
        test_registry_loading(test_dir)
        test_motif_patterns_integration(test_dir)
        test_modular_scanner_integration(test_dir)
        test_detector_class_attribute()
        
        print("=" * 70)
        print("ALL TESTS PASSED ✓")
        print("=" * 70)
        return 0
    except Exception as e:
        print("\n" + "=" * 70)
        print(f"TEST SUITE FAILED: {e}")
        print("=" * 70)
        import traceback
        traceback.print_exc()
        return 1
    finally:
        # Cleanup
        if test_dir and os.path.isdir(test_dir):
            shutil.rmtree(test_dir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
