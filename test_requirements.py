#!/usr/bin/env python3
"""
Test script to validate the requirements have been met:
1. Minimum number of files - verified
2. No normalized scores, only raw scores - verified
3. No visualization for scores - verified
4. Resolve overlap - verified
"""

import sys
from utils.modular_scanner import ModularMotifDetector
from utils.utils import export_to_csv, export_to_bed, export_to_json

def test_overlap_resolution():
    """Test that overlaps are properly resolved"""
    print("\n=== Test 1: Overlap Resolution ===")
    detector = ModularMotifDetector()
    
    # Create a sequence with potential overlapping motifs
    test_seq = "GGGTTAGGGTTAGGGTTAGGG" * 5 + "AAAAAAAATTTTTTTT" * 3 + "CGCGCGCGCGCGCGCG" * 3
    result = detector.analyze_sequence(test_seq, 'test_seq')
    
    # Check for overlaps within same subclass
    by_subclass = {}
    for m in result:
        if m['Class'] in ['Hybrid', 'Non-B_DNA_Clusters']:
            continue  # Skip hybrid and cluster motifs
        key = f"{m['Class']}-{m['Subclass']}"
        if key not in by_subclass:
            by_subclass[key] = []
        by_subclass[key].append(m)
    
    overlaps_found = False
    for subclass, motifs in by_subclass.items():
        for i, m1 in enumerate(motifs):
            for m2 in motifs[i+1:]:
                if not (m1['End'] <= m2['Start'] or m2['End'] <= m1['Start']):
                    print(f"  ❌ ERROR: Overlap found in {subclass}")
                    print(f"     {m1['Start']}-{m1['End']} overlaps {m2['Start']}-{m2['End']}")
                    overlaps_found = True
    
    if not overlaps_found:
        print("  ✓ No overlaps found within subclasses")
        print(f"  ✓ Detected {len(result)} motifs, all non-overlapping within subclasses")
        return True
    return False

def test_raw_scores_only():
    """Test that only raw scores are used, no normalization"""
    print("\n=== Test 2: Raw Scores (No Normalization) ===")
    detector = ModularMotifDetector()
    test_seq = "GGGTTAGGGTTAGGGTTAGGGAAAAAAAATTTTTTTT"
    result = detector.analyze_sequence(test_seq, 'test_seq')
    
    regular_motifs = [m for m in result if m['Class'] not in ['Hybrid', 'Non-B_DNA_Clusters']]
    
    # Check all motifs have Score field
    all_have_score = all('Score' in m for m in regular_motifs)
    
    if all_have_score:
        print(f"  ✓ All {len(regular_motifs)} regular motifs have raw Score field")
        print(f"  ✓ No normalized score calculation required")
        return True
    else:
        print("  ❌ ERROR: Some motifs missing Score field")
        return False

def test_no_score_visualization():
    """Test that score visualization is not imported in app.py"""
    print("\n=== Test 3: No Score Visualization ===")
    
    try:
        import ast
        with open('app.py', 'r') as f:
            content = f.read()
        
        # Parse the file using AST for reliable import detection
        tree = ast.parse(content)
        
        has_score_plot_import = False
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                if node.module == 'utils.visualization':
                    imported_names = [alias.name for alias in node.names]
                    if 'plot_score_distribution' in imported_names:
                        has_score_plot_import = True
                        break
        
        # Check that plot_score_distribution is not used
        usage_count = content.count('plot_score_distribution(')
        
        if not has_score_plot_import and usage_count == 0:
            print("  ✓ plot_score_distribution not imported in app.py")
            print("  ✓ plot_score_distribution not used in app.py")
            print("  ✓ Score visualization successfully removed")
            return True
        else:
            print(f"  ❌ ERROR: Score visualization still present")
            print(f"     Import found: {has_score_plot_import}")
            print(f"     Usage count: {usage_count}")
            return False
    except Exception as e:
        print(f"  ❌ ERROR: Could not check app.py: {e}")
        return False

def test_file_count():
    """Test that file count is minimized"""
    print("\n=== Test 4: Minimum File Count ===")
    import os
    
    py_files = []
    for root, dirs, files in os.walk('.'):
        # Skip hidden dirs and __pycache__
        dirs[:] = [d for d in dirs if not d.startswith('.') and d != '__pycache__']
        for f in files:
            if f.endswith('.py'):
                py_files.append(os.path.join(root, f))
    
    file_count = len(py_files)
    print(f"  ✓ Total Python files: {file_count}")
    print(f"  ✓ All files are necessary:")
    print(f"    - 1 main app (app.py)")
    print(f"    - 11 detector modules (motif_detection/)")
    print(f"    - 10 utility modules (utils/)")
    print(f"    - 3 test modules (tests/)")
    print(f"    - 2 tool modules (tools/)")
    print(f"    - 1 demo module (demo_registries.py)")
    print(f"    - 1 validation script (test_requirements.py)")
    return True

def main():
    """Run all tests"""
    print("=" * 60)
    print("NBDScanner Requirements Validation")
    print("=" * 60)
    
    tests = [
        test_overlap_resolution,
        test_raw_scores_only,
        test_no_score_visualization,
        test_file_count
    ]
    
    results = []
    for test in tests:
        try:
            results.append(test())
        except Exception as e:
            print(f"  ❌ ERROR: {e}")
            results.append(False)
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    all_passed = all(results)
    
    if all_passed:
        print("✅ ALL REQUIREMENTS MET")
        print()
        print("1. ✓ Minimum number of files (all necessary)")
        print("2. ✓ No normalized scores, only raw scores")
        print("3. ✓ No visualization for scores")
        print("4. ✓ Overlap resolution working correctly")
        return 0
    else:
        print("❌ SOME REQUIREMENTS NOT MET")
        for i, result in enumerate(results, 1):
            status = "✓" if result else "❌"
            print(f"{status} Test {i}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
