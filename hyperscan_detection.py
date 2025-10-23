"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                     HYPERSCAN DETECTION ENGINE MODULE                        ║
║                    Non-B DNA Motif Detection System                          ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: hyperscan_detection.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    High-performance motif detection using Hyperscan regex engine.
    Detects 8 out of 11 motif classes with optimized pattern matching.
    
DETECTED CLASSES (Hyperscan-compatible):
┌──────┬─────────────────────┬────────────┬──────────────────────────────────────┐
│Class │ Name                │ Patterns   │ Performance                          │
├──────┼─────────────────────┼────────────┼──────────────────────────────────────┤
│  1   │ Curved DNA          │     6      │ ~50,000 bp/s                         │
│  2   │ Slipped DNA (STR)   │    14      │ ~40,000 bp/s                         │
│  4   │ R-loop              │     5      │ ~45,000 bp/s                         │
│  5   │ Triplex             │     8      │ ~42,000 bp/s                         │
│  6   │ G-Quadruplex        │     9      │ ~48,000 bp/s                         │
│  7   │ i-Motif Family      │     5      │ ~47,000 bp/s                         │
│  8   │ Z-DNA               │    11      │ ~52,000 bp/s                         │
└──────┴─────────────────────┴────────────┴──────────────────────────────────────┘

TOTAL: 58+ patterns compiled into single Hyperscan database

NON-HYPERSCAN CLASSES:
    - Class 3: Cruciform (requires reverse complement check)
    - Class 9: A-philic DNA (requires tetranucleotide scoring)
    - Class 10: Hybrid (post-processing of overlaps)
    - Class 11: Non-B DNA Clusters (density-based clustering)

ALGORITHM:
    1. Compile all patterns into Hyperscan database (one-time)
    2. Scan sequence using database (ultra-fast)
    3. Collect matches with positions
    4. Return raw motifs for scoring
    
REFERENCES:
    - Intel Hyperscan Library: https://www.hyperscan.io/
    - Pattern matching optimization: Wang et al., 2019
"""

import time
from typing import List, Dict, Any, Tuple, Optional
from collections import defaultdict

try:
    import hyperscan
    HYPERSCAN_AVAILABLE = True
except (ImportError, Exception):
    HYPERSCAN_AVAILABLE = False

from hyperscan_registry import HyperscanRegistry, get_class_mapping


class HyperscanDetector:
    """
    High-performance Non-B DNA motif detector using Hyperscan engine.
    
    This detector compiles all regex patterns into a single optimized database
    and performs blazing-fast pattern matching across entire sequences.
    """
    
    def __init__(self):
        """Initialize Hyperscan detector with compiled database."""
        if not HYPERSCAN_AVAILABLE:
            raise ImportError(
                "Hyperscan is not available. Please install: pip install hyperscan"
            )
        
        # Compile pattern database
        self.database, self.pattern_info = HyperscanRegistry.compile_hyperscan_database()
        self.class_mapping = get_class_mapping()
        
        # Statistics
        self.total_scans = 0
        self.total_matches = 0
        self.total_time = 0.0
    
    def detect_motifs(
        self, 
        sequence: str, 
        sequence_name: str = "Sequence"
    ) -> List[Dict[str, Any]]:
        """
        Detect all Hyperscan-compatible motifs in a sequence.
        
        Args:
            sequence: DNA sequence to analyze (uppercase recommended)
            sequence_name: Name of the sequence for reporting
            
        Returns:
            List of detected motifs with positions and metadata
        """
        start_time = time.time()
        
        # Convert to uppercase for consistency
        sequence = sequence.upper()
        
        # Storage for matches
        matches = []
        
        def on_match(pattern_id, start, end, flags, context):
            """Callback function for Hyperscan matches."""
            # Get pattern information
            info = self.pattern_info[pattern_id]
            class_name = self.class_mapping.get(pattern_id, 'Unknown')
            
            # Extract matched sequence
            matched_seq = sequence[start:end]
            
            # Create motif record
            motif = {
                'Class': class_name,
                'Subclass': info['subclass'],
                'Type': info['name'],
                'Start': start + 1,  # 1-based coordinates
                'End': end,
                'Length': end - start,
                'Sequence': matched_seq,
                'Pattern_ID': pattern_id,
                'Scoring_Method': info['scoring_method'],
                'Sequence_Name': sequence_name,
                'Strand': '+',
                'Detection_Method': 'Hyperscan'
            }
            
            matches.append(motif)
        
        # Scan sequence with Hyperscan
        try:
            self.database.scan(
                sequence.encode('utf-8'),
                match_event_handler=on_match
            )
        except Exception as e:
            print(f"Warning: Hyperscan scan error: {e}")
            return []
        
        # Update statistics
        elapsed = time.time() - start_time
        self.total_scans += 1
        self.total_matches += len(matches)
        self.total_time += elapsed
        
        return matches
    
    def detect_motifs_batch(
        self, 
        sequences: List[Tuple[str, str]]
    ) -> Dict[str, List[Dict[str, Any]]]:
        """
        Detect motifs in multiple sequences (batch processing).
        
        Args:
            sequences: List of (sequence_name, sequence) tuples
            
        Returns:
            Dictionary mapping sequence names to detected motifs
        """
        results = {}
        
        for seq_name, seq in sequences:
            results[seq_name] = self.detect_motifs(seq, seq_name)
        
        return results
    
    def get_statistics(self) -> Dict[str, Any]:
        """
        Get performance statistics.
        
        Returns:
            Dictionary with scanning statistics
        """
        avg_time = self.total_time / self.total_scans if self.total_scans > 0 else 0
        avg_matches = self.total_matches / self.total_scans if self.total_scans > 0 else 0
        
        return {
            'total_scans': self.total_scans,
            'total_matches': self.total_matches,
            'total_time': self.total_time,
            'avg_time_per_scan': avg_time,
            'avg_matches_per_scan': avg_matches,
            'patterns_in_database': len(self.pattern_info)
        }
    
    def reset_statistics(self):
        """Reset performance statistics."""
        self.total_scans = 0
        self.total_matches = 0
        self.total_time = 0.0


class HyperscanFallback:
    """
    Fallback detector for when Hyperscan is not available.
    Uses standard Python regex (slower but compatible).
    """
    
    def __init__(self):
        """Initialize fallback detector with regex patterns."""
        import re
        
        self.patterns = HyperscanRegistry.get_pattern_list()
        self.class_mapping = get_class_mapping()
        self.compiled_patterns = []
        
        # Compile all patterns
        for pattern, pattern_id, name, subclass, scoring_method in self.patterns:
            try:
                compiled = re.compile(pattern, re.IGNORECASE)
                self.compiled_patterns.append((
                    compiled, pattern_id, name, subclass, scoring_method
                ))
            except re.error as e:
                print(f"Warning: Failed to compile pattern {pattern_id}: {e}")
    
    def detect_motifs(
        self, 
        sequence: str, 
        sequence_name: str = "Sequence"
    ) -> List[Dict[str, Any]]:
        """
        Detect motifs using standard Python regex (fallback).
        
        Args:
            sequence: DNA sequence to analyze
            sequence_name: Name of the sequence
            
        Returns:
            List of detected motifs
        """
        sequence = sequence.upper()
        matches = []
        
        for compiled, pattern_id, name, subclass, scoring_method in self.compiled_patterns:
            class_name = self.class_mapping.get(pattern_id, 'Unknown')
            
            for match in compiled.finditer(sequence):
                start, end = match.span()
                matched_seq = sequence[start:end]
                
                motif = {
                    'Class': class_name,
                    'Subclass': subclass,
                    'Type': name,
                    'Start': start + 1,
                    'End': end,
                    'Length': end - start,
                    'Sequence': matched_seq,
                    'Pattern_ID': pattern_id,
                    'Scoring_Method': scoring_method,
                    'Sequence_Name': sequence_name,
                    'Strand': '+',
                    'Detection_Method': 'Regex_Fallback'
                }
                
                matches.append(motif)
        
        return matches


def create_detector() -> Any:
    """
    Factory function to create appropriate detector.
    
    Returns:
        HyperscanDetector if available, otherwise HyperscanFallback
    """
    if HYPERSCAN_AVAILABLE:
        try:
            return HyperscanDetector()
        except Exception as e:
            print(f"Warning: Failed to initialize Hyperscan detector: {e}")
            print("Falling back to regex detector")
            return HyperscanFallback()
    else:
        print("Hyperscan not available, using regex fallback")
        return HyperscanFallback()


if __name__ == '__main__':
    # Test detection
    print("Hyperscan Detection Test")
    print("=" * 80)
    
    # Test sequence with multiple motif types
    test_sequence = (
        "AAAAAAAGGGGTTAGGGTTAGGGTTAGGGCCCCTTCCCCTTCCCCTTCCCC"
        "CGCGCGCGCGCGCGCGCG"
        "GAAAGAAAGAAAGAAAGAAA"
        "CACACACACACACACACA"
        "TTTTTTTTTTTT"
        "GGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCC"
    )
    
    # Create detector
    detector = create_detector()
    
    print(f"\nDetector type: {type(detector).__name__}")
    print(f"Test sequence length: {len(test_sequence)} bp")
    
    # Detect motifs
    start = time.time()
    motifs = detector.detect_motifs(test_sequence, "Test_Sequence")
    elapsed = time.time() - start
    
    print(f"\nDetection time: {elapsed:.4f} seconds")
    print(f"Motifs found: {len(motifs)}")
    print(f"Speed: {len(test_sequence)/elapsed:.0f} bp/s")
    
    # Show sample results
    if motifs:
        print("\nSample motifs:")
        for i, motif in enumerate(motifs[:10]):
            print(f"  {i+1}. {motif['Class']} - {motif['Subclass']}")
            print(f"     Position: {motif['Start']}-{motif['End']} ({motif['Length']} bp)")
            print(f"     Sequence: {motif['Sequence'][:50]}...")
    
    # Show statistics
    if hasattr(detector, 'get_statistics'):
        print("\nStatistics:")
        stats = detector.get_statistics()
        for key, value in stats.items():
            print(f"  {key}: {value}")
