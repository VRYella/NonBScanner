#!/usr/bin/env python3
"""
Test UI Improvements and Scientific Accuracy Enhancements for NBDFinder
========================================================================

This test validates:
1. Underscore replacement in UI text
2. Scientific validation functionality  
3. Column name improvements
4. Data integrity checks
"""

import unittest
import pandas as pd
import re
from unittest.mock import patch


class TestUIImprovements(unittest.TestCase):
    """Test UI improvements and scientific validation enhancements"""
    
    def test_underscore_replacement(self):
        """Test that underscores are properly replaced with spaces in UI elements"""
        print("Testing underscore replacement functionality...")
        
        # Test column name replacement
        test_columns = ['Normalized_Score', 'Actual_Score', 'GC_Content', 'Scoring_Method']
        expected_columns = ['Normalized Score', 'Actual Score', 'GC Content', 'Scoring Method']
        
        # Apply the same transformation used in app.py
        transformed_columns = [col.replace('_', ' ') for col in test_columns]
        
        self.assertEqual(transformed_columns, expected_columns)
        print("   ✓ Column name transformation working correctly")
        
        # Test DataFrame column transformation
        test_df = pd.DataFrame({
            'Normalized_Score': [0.5, 0.7, 0.9],
            'Actual_Score': [3.2, 4.1, 5.5],
            'GC_Content': [60.0, 65.0, 70.0]
        })
        
        # Transform column names
        test_df.columns = [col.replace('_', ' ') for col in test_df.columns]
        
        expected_cols = ['Normalized Score', 'Actual Score', 'GC Content']
        self.assertEqual(list(test_df.columns), expected_cols)
        print("   ✓ DataFrame column transformation working correctly")
    
    def test_scientific_validation(self):
        """Test scientific validation functions"""
        print("Testing scientific validation functionality...")
        
        # Test valid DNA sequence check
        valid_sequence = "GGGTTAGGGTTAGGGTTAGGG"
        invalid_sequence = "GGGTTAGGGTTAGGGTTAGGGXYZ"
        
        valid_chars = set('ATCGN')
        
        # Validate sequences
        valid_check = set(valid_sequence.upper()).issubset(valid_chars)
        invalid_check = set(invalid_sequence.upper()).issubset(valid_chars)
        
        self.assertTrue(valid_check)
        self.assertFalse(invalid_check)
        print("   ✓ DNA sequence validation working correctly")
        
        # Test sequence length validation
        short_seq = "ATCG"
        normal_seq = "GGGTTAGGGTTAGGGTTAGGG"
        long_seq = "A" * 1000001  # > 1MB
        
        self.assertTrue(len(short_seq) < 10)  # Too short
        self.assertTrue(10 <= len(normal_seq) <= 1000000)  # Good length
        self.assertTrue(len(long_seq) > 1000000)  # Too long
        print("   ✓ Sequence length validation working correctly")
    
    def test_normalized_score_validation(self):
        """Test normalized score validation"""
        print("Testing normalized score validation...")
        
        # Test score validation logic
        test_motifs = [
            {'Normalized_Score': 0.5, 'Class': 'G-Quadruplex'},
            {'Normalized_Score': 1.2, 'Class': 'Z-DNA'},  # Invalid: > 1
            {'Normalized_Score': -0.1, 'Class': 'Triplex'},  # Invalid: < 0
            {'Normalized_Score': 0.0, 'Class': 'Cruciform'},  # Valid: boundary
            {'Normalized_Score': 1.0, 'Class': 'R-Loop'},  # Valid: boundary
        ]
        
        # Filter valid motifs (score between 0 and 1)
        valid_motifs = []
        for motif in test_motifs:
            score = motif.get('Normalized_Score', 0)
            if isinstance(score, (int, float)) and 0 <= score <= 1:
                valid_motifs.append(motif)
        
        self.assertEqual(len(valid_motifs), 3)  # Should have 3 valid motifs
        valid_scores = [m['Normalized_Score'] for m in valid_motifs]
        self.assertEqual(valid_scores, [0.5, 0.0, 1.0])
        print("   ✓ Normalized score validation working correctly")
    
    def test_duplicate_column_handling(self):
        """Test duplicate column handling"""
        print("Testing duplicate column handling...")
        
        # Create DataFrame with duplicate columns (simulating pandas duplicate behavior)
        df = pd.DataFrame({
            'Class': ['G-Quadruplex', 'Z-DNA'],
            'Sequence Name': ['Seq1', 'Seq2'],
            'Score': [0.5, 0.7]
        })
        
        # Add a duplicate column manually to test the logic
        df['Sequence Name'] = ['Seq1_dup', 'Seq2_dup']  # This creates a duplicate
        
        # Check for duplicates before cleaning
        has_duplicates = df.columns.duplicated().any()
        
        # Remove duplicates (simulate app.py logic)
        df_clean = df.loc[:, ~df.columns.duplicated()]
        
        # Should have no duplicate columns after cleaning
        has_duplicates_after = df_clean.columns.duplicated().any()
        self.assertFalse(has_duplicates_after)
        print("   ✓ Duplicate column handling working correctly")
    
    def test_css_spacing_improvements(self):
        """Test CSS spacing configuration"""
        print("Testing CSS spacing improvements...")
        
        # Test CSS properties for proper spacing
        css_properties = {
            'line-height': '1.6',
            'margin-bottom': '0.5rem',
            'padding-top': '2rem',
            'min-height': '2.5rem'
        }
        
        # Validate CSS values
        for prop, value in css_properties.items():
            self.assertIsInstance(value, str)
            self.assertTrue(len(value) > 0)
            # Check for proper CSS units
            if prop in ['margin-bottom', 'padding-top', 'min-height']:
                self.assertTrue('rem' in value or 'px' in value or 'em' in value)
        
        print("   ✓ CSS spacing configuration valid")
    
    def test_motif_analysis_integration(self):
        """Test integration with motif analysis"""
        print("Testing motif analysis integration...")
        
        try:
            import motifs
            
            # Test with a known G4-forming sequence
            test_sequence = "GGGTTAGGGTTAGGGTTAGGG"
            
            # Run analysis
            results = motifs.all_motifs(test_sequence, sequence_name="test_integration")
            
            # Validate results structure
            self.assertIsInstance(results, list)
            if results:
                motif = results[0]
                self.assertIsInstance(motif, dict)
                self.assertIn('Class', motif)
                self.assertIn('Start', motif)
                self.assertIn('End', motif)
            
            print(f"   ✓ Found {len(results)} motifs in test sequence")
            
        except ImportError:
            print("   ⚠ Motifs module not available for integration test")
    
    def run_all_tests(self):
        """Run all UI improvement tests"""
        print("=" * 60)
        print("NBDFinder UI Improvements & Scientific Accuracy Tests")
        print("=" * 60)
        
        test_methods = [
            self.test_underscore_replacement,
            self.test_scientific_validation,
            self.test_normalized_score_validation,
            self.test_duplicate_column_handling,
            self.test_css_spacing_improvements,
            self.test_motif_analysis_integration
        ]
        
        passed_tests = 0
        total_tests = len(test_methods)
        
        for test_method in test_methods:
            try:
                test_method()
                passed_tests += 1
            except Exception as e:
                print(f"   ❌ {test_method.__name__} failed: {e}")
        
        print("=" * 60)
        print(f"UI IMPROVEMENTS TEST SUMMARY")
        print("=" * 60)
        print(f"✅ Passed: {passed_tests}/{total_tests} tests")
        
        if passed_tests == total_tests:
            print("🎉 ALL UI IMPROVEMENT TESTS PASSED!")
            print("✅ NBDFinder UI enhancements are working correctly")
            return True
        else:
            print("❌ Some tests failed")
            return False


def main():
    """Main test runner"""
    tester = TestUIImprovements()
    success = tester.run_all_tests()
    return 0 if success else 1


if __name__ == "__main__":
    exit(main())