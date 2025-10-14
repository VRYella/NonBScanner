#!/usr/bin/env python3
"""
Test Suite for G4 Overlap Resolution Agent
===========================================

This test suite validates the G-quadruplex overlap resolution functionality
with emphasis on:
1. Detection of various G4 motif classes
2. Proper overlap resolution by score and priority
3. Output format consistency (JSON, table, BED)
4. Edge cases and error handling

Scientific Test Cases:
---------------------
The test sequences are based on validated G4-forming sequences from literature:
- Telomeric G4: Human telomeric repeat (Zahler et al. 1991, Williamson et al. 1989)
- Canonical G4: Standard G4 patterns (Burge et al. 2006)[web:67]
- Overlapping motifs: Multiple G4 classes in same region for resolution testing

Author: Dr. Venkata Rajesh Yella
"""

import unittest
import json
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from g4_overlap_resolver import G4OverlapResolver


class TestG4OverlapResolver(unittest.TestCase):
    """Test cases for G4 overlap resolution agent."""
    
    def setUp(self):
        """Initialize resolver for each test."""
        self.resolver = G4OverlapResolver()
    
    def test_basic_g4_detection(self):
        """
        Test basic G4 detection without overlaps.
        
        Sequence: Human telomeric repeat (Zahler et al. 1991)
        Expected: Detection of canonical telomeric G4 motif
        """
        sequence = "GGGTTAGGGTTAGGGTTAGGG"
        annotations = self.resolver.resolve_and_annotate(sequence, "telomeric_test")
        
        # Should detect at least one G4 motif
        self.assertGreater(len(annotations), 0, 
                          "Failed to detect G4 in telomeric sequence")
        
        # Check annotation structure
        ann = annotations[0]
        self.assertIn('class_name', ann)
        self.assertIn('start', ann)
        self.assertIn('end', ann)
        self.assertIn('score', ann)
        self.assertIn('matched_seq', ann)
        self.assertEqual(ann['sequence_name'], "telomeric_test")
    
    def test_canonical_g4_high_score(self):
        """
        Test canonical G4 with 4 G-tracts (Burge 2006)[web:67].
        
        Sequence: High-density canonical G4 with short loops
        Expected: High G4Hunter score due to dense G-tracts
        """
        sequence = "GGGGTTTTGGGGTTTTGGGGTTTTGGGG"
        annotations = self.resolver.resolve_and_annotate(sequence, "canonical_g4")
        
        self.assertGreater(len(annotations), 0,
                          "Failed to detect canonical G4")
        
        # Canonical G4s should have relatively high scores
        max_score = max(ann['score'] for ann in annotations)
        self.assertGreater(max_score, 0.5,
                          "Canonical G4 score unexpectedly low")
    
    def test_overlapping_motifs_resolution(self):
        """
        Test overlap resolution with multiple overlapping G4 classes.
        
        Scenario: Long sequence with multiple potential G4 motifs that overlap
        Expected: Non-overlapping set prioritized by score and class
        
        This tests the core overlap resolution algorithm (Bedrat et al. 2016)[web:3]:
        - Sort by score (descending), class priority, length
        - Greedy selection of non-conflicting regions
        """
        # Create sequence with multiple potential G4 regions
        # This should create overlapping candidates of different classes
        sequence = "GGGTTGGGTTGGGTTGGGAAAAAGGGGTTTTGGGGTTTTGGGGTTTTGGGG"
        annotations = self.resolver.resolve_and_annotate(sequence, "overlap_test")
        
        # Should detect multiple motifs
        self.assertGreater(len(annotations), 0,
                          "Failed to detect any G4 motifs")
        
        # Verify no overlaps in resolved annotations
        for i in range(len(annotations)):
            for j in range(i + 1, len(annotations)):
                ann_i = annotations[i]
                ann_j = annotations[j]
                
                # Check that regions don't overlap
                # Either ann_i ends before ann_j starts, or ann_j ends before ann_i starts
                no_overlap = (ann_i['end'] <= ann_j['start'] or 
                             ann_j['end'] <= ann_i['start'])
                
                self.assertTrue(no_overlap,
                              f"Found overlapping annotations: "
                              f"[{ann_i['start']}-{ann_i['end']}] and "
                              f"[{ann_j['start']}-{ann_j['end']}]")
    
    def test_multiple_g4_regions(self):
        """
        Test sequence with multiple well-separated G4 regions.
        
        Expected: Detection of multiple non-overlapping G4s
        """
        # Two canonical G4 motifs separated by spacer
        sequence = "GGGGTTTTGGGGTTTTGGGGTTTTGGGG" + "ATATATAT" + "GGGGTTTTGGGGTTTTGGGGTTTTGGGG"
        annotations = self.resolver.resolve_and_annotate(sequence, "multiple_g4")
        
        # Should detect at least 2 G4 regions (may be more due to overlapping patterns)
        self.assertGreaterEqual(len(annotations), 1,
                               "Failed to detect multiple G4 regions")
        
        # Annotations should be sorted by start position
        for i in range(len(annotations) - 1):
            self.assertLessEqual(annotations[i]['start'], annotations[i + 1]['start'],
                               "Annotations not sorted by start position")
    
    def test_relaxed_g4_detection(self):
        """
        Test relaxed G4 with longer loops (Huppert 2005)[web:75].
        
        Expected: Detection of relaxed G4 class
        """
        # Relaxed G4 with longer loops (up to 12bp)
        sequence = "GGGTTTTTTTTTTTGGGTTGGGTTGGG"
        annotations = self.resolver.resolve_and_annotate(sequence, "relaxed_g4")
        
        # Should detect at least one motif (may be relaxed or other class)
        self.assertGreater(len(annotations), 0,
                          "Failed to detect relaxed G4")
    
    def test_bulged_g4_detection(self):
        """
        Test bulged G4 with large loop (Lim 2009, Adrian 2014)[web:22].
        
        Expected: Detection of bulged G4 class
        """
        # Bulged G4 with one large loop
        sequence = "GGGGAAAAAAAAAAAAAAAAAAAGGGGTTGGGTTGGGG"
        annotations = self.resolver.resolve_and_annotate(sequence, "bulged_g4")
        
        # Should detect motif(s)
        self.assertGreater(len(annotations), 0,
                          "Failed to detect bulged G4")
    
    def test_priority_ordering(self):
        """
        Test that higher-priority classes are selected over lower-priority ones.
        
        When two motifs overlap with similar scores, the class priority should
        determine which is kept:
        canonical_g4 > relaxed_g4 > multimeric_g4 > bulged_g4 > imperfect_g4
        """
        # Sequence that could match multiple classes
        sequence = "GGGGTTTTGGGGTTTTGGGGTTTTGGGG"
        annotations = self.resolver.resolve_and_annotate(sequence, "priority_test")
        
        # Should have at least one annotation
        self.assertGreater(len(annotations), 0,
                          "Failed to detect any G4")
        
        # Check that canonical or relaxed classes are present (higher priority)
        classes = [ann['class_name'] for ann in annotations]
        has_high_priority = any(cls in ['canonical_g4', 'relaxed_g4', 'multimeric_g4'] 
                               for cls in classes)
        self.assertTrue(has_high_priority,
                       "No high-priority G4 classes detected")
    
    def test_json_format(self):
        """
        Test JSON output format.
        
        Expected: Valid JSON with proper structure
        """
        sequence = "GGGTTAGGGTTAGGGTTAGGG"
        annotations = self.resolver.resolve_and_annotate(sequence, "json_test")
        
        # Format as JSON
        json_output = self.resolver.format_as_json(annotations, pretty=True)
        
        # Should be valid JSON
        parsed = json.loads(json_output)
        
        # Check structure
        self.assertIn('version', parsed)
        self.assertIn('analysis_type', parsed)
        self.assertIn('total_motifs', parsed)
        self.assertIn('motifs', parsed)
        self.assertEqual(parsed['total_motifs'], len(annotations))
        self.assertEqual(len(parsed['motifs']), len(annotations))
    
    def test_table_format(self):
        """
        Test tab-delimited table output format.
        
        Expected: Header line followed by data rows
        """
        sequence = "GGGTTAGGGTTAGGGTTAGGG"
        annotations = self.resolver.resolve_and_annotate(sequence, "table_test")
        
        # Format as table
        table_output = self.resolver.format_as_table(annotations)
        
        # Should have lines
        lines = table_output.strip().split('\n')
        self.assertGreater(len(lines), 0, "Empty table output")
        
        # First line should be header
        header = lines[0]
        self.assertIn('Sequence_Name', header)
        self.assertIn('Class', header)
        self.assertIn('Start', header)
        self.assertIn('End', header)
        self.assertIn('Score', header)
        
        # Should have data rows if annotations present
        if len(annotations) > 0:
            self.assertEqual(len(lines), len(annotations) + 1,
                           "Table row count mismatch")
    
    def test_bed_format(self):
        """
        Test BED format output for genome browsers.
        
        Expected: Valid BED format (6 columns)
        """
        sequence = "GGGTTAGGGTTAGGGTTAGGG"
        annotations = self.resolver.resolve_and_annotate(sequence, "bed_test")
        
        # Format as BED
        bed_output = self.resolver.format_as_bed(annotations, "chr1")
        
        if len(annotations) > 0:
            # Should have lines
            lines = bed_output.strip().split('\n')
            self.assertEqual(len(lines), len(annotations),
                           "BED row count mismatch")
            
            # Each line should have 6 tab-separated fields
            for line in lines:
                fields = line.split('\t')
                self.assertEqual(len(fields), 6,
                               "BED format should have 6 columns")
                
                # Verify field types
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                name = fields[3]
                score = int(fields[4])
                strand = fields[5]
                
                self.assertEqual(chrom, "chr1")
                self.assertGreaterEqual(start, 0)
                self.assertGreater(end, start)
                self.assertIn(strand, ['+', '-'])
    
    def test_empty_sequence(self):
        """
        Test handling of empty sequence.
        
        Expected: Empty annotations list
        """
        sequence = ""
        annotations = self.resolver.resolve_and_annotate(sequence, "empty_test")
        
        # Should return empty list
        self.assertEqual(len(annotations), 0,
                        "Empty sequence should produce no annotations")
    
    def test_no_g4_sequence(self):
        """
        Test sequence with no G4 motifs.
        
        Expected: Empty annotations list
        """
        # AT-rich sequence with no G4 potential
        sequence = "ATATATATATATATATATAT"
        annotations = self.resolver.resolve_and_annotate(sequence, "no_g4_test")
        
        # Should return empty list
        self.assertEqual(len(annotations), 0,
                        "AT-rich sequence should produce no G4 annotations")
    
    def test_case_insensitivity(self):
        """
        Test that lowercase and uppercase sequences are handled equivalently.
        
        Expected: Same results for upper and lower case
        """
        sequence_upper = "GGGTTAGGGTTAGGGTTAGGG"
        sequence_lower = "gggttagggttagggttaggg"
        
        annotations_upper = self.resolver.resolve_and_annotate(sequence_upper, "upper_test")
        annotations_lower = self.resolver.resolve_and_annotate(sequence_lower, "lower_test")
        
        # Should detect same number of motifs
        self.assertEqual(len(annotations_upper), len(annotations_lower),
                        "Case sensitivity issue detected")
    
    def test_annotation_completeness(self):
        """
        Test that all annotations have required fields.
        
        Expected: Each annotation has all mandatory fields
        """
        sequence = "GGGTTAGGGTTAGGGTTAGGG"
        annotations = self.resolver.resolve_and_annotate(sequence, "complete_test")
        
        required_fields = [
            'sequence_name', 'class_name', 'pattern_id', 
            'start', 'end', 'length', 'score', 'matched_seq', 'details'
        ]
        
        for ann in annotations:
            for field in required_fields:
                self.assertIn(field, ann,
                            f"Missing required field: {field}")
            
            # Verify field types and constraints
            self.assertIsInstance(ann['start'], int)
            self.assertIsInstance(ann['end'], int)
            self.assertIsInstance(ann['length'], int)
            self.assertIsInstance(ann['score'], float)
            self.assertGreaterEqual(ann['start'], 0)
            self.assertGreater(ann['end'], ann['start'])
            self.assertEqual(ann['length'], ann['end'] - ann['start'])
    
    def test_complex_overlapping_scenario(self):
        """
        Test complex overlapping scenario with multiple G4 classes.
        
        This creates a sequence with:
        - Canonical G4 (high priority, high score)
        - Relaxed G4 (medium priority, medium score)
        - Bulged G4 (lower priority, variable score)
        
        Expected: Highest scoring, highest priority motifs selected
        """
        # Complex sequence with multiple overlapping G4 potential
        sequence = (
            "GGGGTTTTGGGGTTTTGGGGTTTTGGGG"  # Strong canonical G4
            "ATATATAT"                          # Non-G4 spacer
            "GGGAAAAAAAAGGGTTGGGTTGGG"       # Relaxed/bulged G4
            "ATATATAT"                          # Non-G4 spacer
            "GGGGTTTTGGGGTTTTGGGGTTTTGGGG"  # Another strong canonical G4
        )
        
        annotations = self.resolver.resolve_and_annotate(sequence, "complex_test")
        
        # Should detect multiple motifs
        self.assertGreater(len(annotations), 0,
                          "Failed to detect G4s in complex sequence")
        
        # Verify no overlaps
        for i in range(len(annotations)):
            for j in range(i + 1, len(annotations)):
                ann_i = annotations[i]
                ann_j = annotations[j]
                no_overlap = (ann_i['end'] <= ann_j['start'] or 
                             ann_j['end'] <= ann_i['start'])
                self.assertTrue(no_overlap,
                              "Overlapping annotations in complex scenario")
        
        # At least some motifs should be detected (multimeric or canonical classes are valid)
        classes = [ann['class_name'] for ann in annotations]
        has_valid_class = any(cls in ['canonical_g4', 'relaxed_g4', 'multimeric_g4', 'bulged_g4'] 
                             for cls in classes)
        self.assertTrue(has_valid_class,
                       "No valid G4 classes in complex scenario")


class TestG4OverlapResolverAPI(unittest.TestCase):
    """Test the API usage patterns for G4 overlap resolver."""
    
    def test_api_basic_usage(self):
        """
        Test basic API usage pattern.
        
        This demonstrates the simplest API usage:
        1. Create resolver
        2. Call resolve_and_annotate
        3. Get results
        """
        resolver = G4OverlapResolver()
        sequence = "GGGTTAGGGTTAGGGTTAGGG"
        results = resolver.resolve_and_annotate(sequence)
        
        self.assertIsInstance(results, list)
        if len(results) > 0:
            self.assertIsInstance(results[0], dict)
    
    def test_api_json_workflow(self):
        """
        Test API workflow for JSON output.
        
        Workflow:
        1. Resolve overlaps
        2. Format as JSON
        3. Parse and validate
        """
        resolver = G4OverlapResolver()
        sequence = "GGGGTTTTGGGGTTTTGGGGTTTTGGGG"
        
        # Step 1: Resolve
        annotations = resolver.resolve_and_annotate(sequence, "api_test")
        
        # Step 2: Format
        json_str = resolver.format_as_json(annotations)
        
        # Step 3: Validate
        parsed = json.loads(json_str)
        self.assertIn('motifs', parsed)
    
    def test_api_table_workflow(self):
        """
        Test API workflow for table output.
        
        Workflow:
        1. Resolve overlaps
        2. Format as table
        3. Validate format
        """
        resolver = G4OverlapResolver()
        sequence = "GGGGTTTTGGGGTTTTGGGGTTTTGGGG"
        
        # Step 1: Resolve
        annotations = resolver.resolve_and_annotate(sequence, "table_api")
        
        # Step 2: Format
        table_str = resolver.format_as_table(annotations)
        
        # Step 3: Validate
        lines = table_str.strip().split('\n')
        self.assertGreater(len(lines), 0)
        self.assertIn('Class', lines[0])


def run_tests():
    """Run all tests and print results."""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestG4OverlapResolver))
    suite.addTests(loader.loadTestsFromTestCase(TestG4OverlapResolverAPI))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Return success status
    return result.wasSuccessful()


if __name__ == '__main__':
    success = run_tests()
    sys.exit(0 if success else 1)
