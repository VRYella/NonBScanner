#!/usr/bin/env python3
"""
Comprehensive tests for QmRLFS R-loop detection integration
"""

import unittest
import sys
import os

# Add current directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from qmrlfs_finder import QmRLFSDetector, detect_r_loop_qmrlfs, qmrlfs_to_nbdscanner_format
from motif_patterns import MotifScoring, PatternRegistry


class TestQmRLFSDetector(unittest.TestCase):
    """Test cases for QmRLFSDetector class"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.detector = QmRLFSDetector(models="m1,m2", quick_mode=False)
        self.quick_detector = QmRLFSDetector(models="m1,m2", quick_mode=True)
        
        # Test sequences
        self.test_sequences = {
            'simple_rlfs': 'GGGGATTTTGGGGCCCCGGGGAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCAAAAAAAAATTTTTTTTT',
            'weak_signal': 'ATATATATATGCGCGCATATATATGCGCGCATATATGCGCGC',
            'g4_like': 'GGGTTTAGGGTTAGGGTTAGGGAAAACCCAAACCCAAACCCAAACCC',
            'high_gc': 'GGGCCCGGGATGCCCGGGAAAGGGCCCTTTGGGCCCAAAGGGCCCTTTGGGCCCAAA',
            'no_rlfs': 'ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT'
        }
    
    def test_detector_initialization(self):
        """Test QmRLFSDetector initialization"""
        # Test default parameters
        detector = QmRLFSDetector()
        self.assertIn('m1', detector.models)
        self.assertIn('m2', detector.models)
        self.assertEqual(detector.min_perc_g_riz, 50.0)
        self.assertEqual(detector.min_perc_g_rez, 40.0)
        
        # Test custom parameters
        custom_detector = QmRLFSDetector(
            models="m1",
            min_perc_g_riz=60.0,
            min_perc_g_rez=35.0,
            quick_mode=True
        )
        self.assertIn('m1', custom_detector.models)
        self.assertNotIn('m2', custom_detector.models)
        self.assertEqual(custom_detector.min_perc_g_riz, 60.0)
        self.assertEqual(custom_detector.min_perc_g_rez, 35.0)
        self.assertTrue(custom_detector.quick_mode)
    
    def test_percent_g_calculation(self):
        """Test G percentage calculation"""
        test_cases = [
            ('GGGG', 100.0),
            ('GGGGAAAA', 50.0),
            ('AAAA', 0.0),
            ('ATCG', 25.0),
            ('', 0.0)
        ]
        
        for seq, expected in test_cases:
            with self.subTest(sequence=seq):
                result = self.detector.percent_g(seq)
                self.assertAlmostEqual(result, expected, places=1)
    
    def test_riz_search_model_1(self):
        """Test RIZ search with model 1"""
        seq = self.test_sequences['simple_rlfs']
        results = self.detector.riz_search(seq, 'm1')
        
        self.assertIsInstance(results, list)
        self.assertGreater(len(results), 0, "Should find at least one RIZ region")
        
        for riz in results:
            self.assertIn('start', riz)
            self.assertIn('end', riz)
            self.assertIn('perc_g', riz)
            self.assertIn('seq', riz)
            self.assertGreaterEqual(riz['perc_g'], self.detector.min_perc_g_riz)
            self.assertLess(riz['start'], riz['end'])
    
    def test_riz_search_model_2(self):
        """Test RIZ search with model 2"""
        seq = self.test_sequences['simple_rlfs']
        results = self.detector.riz_search(seq, 'm2')
        
        self.assertIsInstance(results, list)
        
        for riz in results:
            self.assertIn('model', riz)
            self.assertEqual(riz['model'], 'm2')
            self.assertGreaterEqual(riz['perc_g'], self.detector.min_perc_g_riz)
    
    def test_rez_search(self):
        """Test REZ search functionality"""
        seq = self.test_sequences['simple_rlfs']
        # Find a RIZ first
        riz_results = self.detector.riz_search(seq, 'm1')
        self.assertGreater(len(riz_results), 0)
        
        riz = riz_results[0]
        remaining_seq = seq[riz['end']:]
        
        if len(remaining_seq) >= 20:  # Only test if there's enough sequence
            rez_result = self.detector.rez_search(riz['end'], remaining_seq, len(seq))
            
            if rez_result:  # If REZ found
                self.assertIn('start', rez_result)
                self.assertIn('end', rez_result)
                self.assertIn('perc_g', rez_result)
                self.assertIn('seq', rez_result)
                self.assertLess(rez_result['start'], rez_result['end'])
    
    def test_search_rlfs_single_strand(self):
        """Test RLFS search on single strand"""
        seq = self.test_sequences['simple_rlfs']
        results = self.detector.search_rlfs(seq, "+", "test_seq")
        
        self.assertIsInstance(results, list)
        
        for result in results:
            self.assertIn('sequence_name', result)
            self.assertIn('model', result)
            self.assertIn('strand', result)
            self.assertIn('qmrlfs_score', result)
            self.assertIn('start_RIZ', result)
            self.assertIn('end_RIZ', result)
            self.assertIn('start_REZ', result)
            self.assertIn('end_REZ', result)
            
            self.assertEqual(result['strand'], '+')
            self.assertEqual(result['sequence_name'], 'test_seq')
            self.assertGreaterEqual(result['qmrlfs_score'], 0.0)
            self.assertLessEqual(result['qmrlfs_score'], 1.0)
    
    def test_analyze_sequence_both_strands(self):
        """Test sequence analysis on both strands"""
        seq = self.test_sequences['g4_like']
        results = self.detector.analyze_sequence(seq, "test_seq", analyze_both_strands=True)
        
        self.assertIsInstance(results, list)
        
        # Check for both strand representations
        strands = set(result['strand'] for result in results)
        if len(results) > 1:
            # Should have results from both strands or multiple from one strand
            self.assertTrue('+' in strands or '-' in strands)
    
    def test_qmrlfs_score_calculation(self):
        """Test QmRLFS scoring calculation"""
        # Mock RIZ and REZ data
        riz = {
            'perc_g': 60.0,
            'length': 20,
            'G3s': 3,
            'G4s': 2
        }
        rez = {
            'perc_g': 45.0,
            'length': 100
        }
        
        score = self.detector.calculate_qmrlfs_score(riz, rez, 10)
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)
        self.assertIsInstance(score, float)
    
    def test_quick_mode_performance(self):
        """Test that quick mode runs faster (basic test)"""
        seq = self.test_sequences['simple_rlfs']
        
        import time
        
        # Test regular mode
        start_time = time.time()
        results_regular = self.detector.analyze_sequence(seq, "test", analyze_both_strands=False)
        regular_time = time.time() - start_time
        
        # Test quick mode
        start_time = time.time()
        results_quick = self.quick_detector.analyze_sequence(seq, "test", analyze_both_strands=False)
        quick_time = time.time() - start_time
        
        # Quick mode should be faster or at least not significantly slower
        self.assertLessEqual(quick_time, regular_time * 2)  # Allow some tolerance
        
        # Both should find some results
        if results_regular:
            self.assertGreater(len(results_quick), 0, "Quick mode should also find results")
    
    def test_empty_sequence(self):
        """Test handling of empty sequences"""
        results = self.detector.analyze_sequence("", "empty")
        self.assertEqual(len(results), 0)
        
        results = self.detector.analyze_sequence("ATCG", "short")  # Too short
        self.assertEqual(len(results), 0)
    
    def test_no_rlfs_sequence(self):
        """Test sequence with no RLFS"""
        seq = self.test_sequences['no_rlfs']
        results = self.detector.analyze_sequence(seq, "no_rlfs")
        
        # Should return empty list or very low scores
        if results:
            for result in results:
                self.assertLess(result['qmrlfs_score'], 0.3)  # Very low score


class TestConvenienceFunctions(unittest.TestCase):
    """Test convenience functions for NBDScanner integration"""
    
    def test_detect_r_loop_qmrlfs(self):
        """Test detect_r_loop_qmrlfs convenience function"""
        seq = 'GGGGATTTTGGGGCCCCGGGGAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCAAAAAAAAATTTTTTTTT'
        results = detect_r_loop_qmrlfs(seq, "test_seq", models="m1,m2")
        
        self.assertIsInstance(results, list)
        if results:
            self.assertIn('sequence_name', results[0])
            self.assertIn('qmrlfs_score', results[0])
    
    def test_qmrlfs_to_nbdscanner_format(self):
        """Test conversion to NBDScanner format"""
        seq = 'GGGGATTTTGGGGCCCCGGGGAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCAAAAAAAAATTTTTTTTT'
        qmrlfs_results = detect_r_loop_qmrlfs(seq, "test_seq")
        
        if qmrlfs_results:
            motifs = qmrlfs_to_nbdscanner_format(qmrlfs_results)
            
            self.assertIsInstance(motifs, list)
            self.assertGreater(len(motifs), 0)
            
            for motif in motifs:
                # Check required NBDScanner format fields
                required_fields = ['ID', 'Class', 'Subclass', 'Start', 'End', 'Length', 
                                 'Strand', 'Raw_Score', 'Normalized_Score', 'Method']
                for field in required_fields:
                    self.assertIn(field, motif)
                
                self.assertEqual(motif['Class'], 'R-loop')
                self.assertIn('QmRLFS', motif['Subclass'])
                self.assertEqual(motif['Method'], 'QmRLFS_Pure_Python')


class TestMotifPatternsIntegration(unittest.TestCase):
    """Test integration with motif_patterns.py"""
    
    def test_qmrlfs_patterns_in_registry(self):
        """Test that QmRLFS patterns are in the pattern registry"""
        all_patterns = PatternRegistry.get_all_patterns()
        
        self.assertIn('r_loop', all_patterns)
        r_loop_patterns = all_patterns['r_loop']
        
        # Check for QmRLFS pattern groups
        self.assertIn('qmrlfs_model_1', r_loop_patterns)
        self.assertIn('qmrlfs_model_2', r_loop_patterns)
        
        # Check pattern content
        m1_patterns = r_loop_patterns['qmrlfs_model_1']
        m2_patterns = r_loop_patterns['qmrlfs_model_2']
        
        self.assertGreater(len(m1_patterns), 0)
        self.assertGreater(len(m2_patterns), 0)
    
    def test_qmrlfs_subclass_mapping(self):
        """Test QmRLFS subclasses in mapping"""
        subclass_mapping = PatternRegistry.get_subclass_mapping()
        
        self.assertIn('r_loop', subclass_mapping)
        r_loop_subclasses = subclass_mapping['r_loop']
        
        self.assertIn('QmRLFS-m1', r_loop_subclasses)
        self.assertIn('QmRLFS-m2', r_loop_subclasses)
    
    def test_qmrlfs_scoring_function(self):
        """Test QmRLFS scoring function in MotifScoring"""
        test_sequences = [
            'GGGGATTTTGGGGCCCCGGGGAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCAAAAAAAAATTTTTTTTT',
            'ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT',
            'GGGTTTAGGGTTAGGGTTAGGGAAAACCCAAACCCAAACCCAAACCC'
        ]
        
        for seq in test_sequences:
            with self.subTest(sequence=seq[:30]):
                score = MotifScoring.qmrlfs_score(seq)
                self.assertIsInstance(score, float)
                self.assertGreaterEqual(score, 0.0)
                self.assertLessEqual(score, 1.0)
    
    def test_scoring_methods_list(self):
        """Test that qmrlfs_score is in the scoring methods list"""
        from motif_patterns import get_pattern_statistics
        stats = get_pattern_statistics()
        
        self.assertIn('qmrlfs_score', stats['scoring_methods'])


class TestQmRLFSStaticMethods(unittest.TestCase):
    """Test static methods and utilities"""
    
    def test_qmrlfs_score_function_static(self):
        """Test static QmRLFS scoring function"""
        test_cases = [
            ('GGGGATTTTGGGGCCCCGGGGAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCAAAAAAAAATTTTTTTTT', True),
            ('ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT', False),
            ('GGGTTTAGGGTTAGGGTTAGGGAAAACCCAAACCCAAACCCAAACCC', True),
            ('', False),
            ('ATCG', False)  # Too short
        ]
        
        for seq, should_have_score in test_cases:
            with self.subTest(sequence=seq[:30]):
                score = QmRLFSDetector.qmrlfs_score_function(seq)
                self.assertIsInstance(score, float)
                self.assertGreaterEqual(score, 0.0)
                self.assertLessEqual(score, 1.0)
                
                if should_have_score:
                    self.assertGreater(score, 0.0)
                # Note: We don't test for exactly 0.0 for negative cases
                # because the algorithm might still find weak signals


def run_tests():
    """Run all tests"""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    test_classes = [
        TestQmRLFSDetector,
        TestConvenienceFunctions,
        TestMotifPatternsIntegration,
        TestQmRLFSStaticMethods
    ]
    
    for test_class in test_classes:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result.wasSuccessful()


if __name__ == '__main__':
    print("Running QmRLFS Integration Tests")
    print("=" * 50)
    
    success = run_tests()
    
    if success:
        print("\n" + "=" * 50)
        print("✅ All QmRLFS tests passed successfully!")
        sys.exit(0)
    else:
        print("\n" + "=" * 50)
        print("❌ Some tests failed!")
        sys.exit(1)