"""
Slipped DNA Motif Detector
=========================

Detects tandem repeat-induced slippage patterns including:
- Short tandem repeats (STRs)
- Direct repeats

Based on Wells 2005 and microsatellite instability research.
"""

import re
from typing import List, Dict, Any, Tuple
from .base_detector import BaseMotifDetector


class SlippedDNADetector(BaseMotifDetector):
    """Detector for slipped DNA motifs"""
    
    def get_motif_class_name(self) -> str:
        return "Slipped_DNA"
    
    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return slipped DNA patterns"""
        return {
            'short_tandem_repeats': [
                (r'([ATGC])\1{9,}', 'SLP_2_1', 'Mononucleotide repeat', 'STR', 10, 'instability_score', 0.98, 'Replication slippage', 'Schlötterer 2000'),
                (r'([ATGC]{2})\1{4,}', 'SLP_2_2', 'Dinucleotide repeat', 'STR', 10, 'instability_score', 0.95, 'Microsatellite instability', 'Weber 1989'),
                (r'([ATGC]{3})\1{3,}', 'SLP_2_3', 'Trinucleotide repeat', 'STR', 12, 'instability_score', 0.92, 'Expansion diseases', 'Ashley 1993'),
                (r'([ATGC]{4})\1{2,}', 'SLP_2_4', 'Tetranucleotide repeat', 'STR', 12, 'instability_score', 0.85, 'Genetic polymorphisms', 'Edwards 1991'),
                (r'(CA)\1{4,}', 'SLP_2_5', 'CA repeat', 'STR', 10, 'instability_score', 0.95, 'Common microsatellite', 'Weber 1989'),
                (r'(CGG)\1{3,}', 'SLP_2_6', 'CGG repeat', 'STR', 12, 'instability_score', 0.90, 'Fragile X syndrome', 'Verkerk 1991'),
            ],
            'direct_repeats': [
                (r'([ATGC]{5,20})(?:[ATGC]{0,100})\1', 'SLP_2_7', 'Direct repeat', 'Direct Repeat', 10, 'repeat_score', 0.80, 'Recombination hotspots', 'Jeffreys 1985'),
                (r'([ATGC]{10,50})(?:[ATGC]{0,200})\1', 'SLP_2_8', 'Long direct repeat', 'Direct Repeat', 20, 'repeat_score', 0.75, 'Genomic instability', 'Lupski 1998'),
            ]
        }
    
    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        """Calculate instability score for the sequence"""
        scoring_method = pattern_info[5] if len(pattern_info) > 5 else 'instability_score'
        
        if scoring_method == 'instability_score':
            return self._instability_score(sequence)
        elif scoring_method == 'repeat_score':
            return self._repeat_score(sequence)
        else:
            return self._instability_score(sequence)
    
    def _instability_score(self, sequence: str) -> float:
        """
        Repeat instability scoring for slipped DNA structures
        """
        if len(sequence) < 6:
            return 0.0
        
        # Find repeating units
        max_instability = 0
        
        for unit_length in range(1, min(7, len(sequence) // 3)):
            for i in range(len(sequence) - unit_length * 2):
                unit = sequence[i:i + unit_length]
                count = 1
                
                # Count consecutive repeats
                pos = i + unit_length
                while pos + unit_length <= len(sequence) and sequence[pos:pos + unit_length] == unit:
                    count += 1
                    pos += unit_length
                
                if count >= 3:  # At least 3 repeats
                    instability = count * (unit_length ** 0.5)
                    max_instability = max(max_instability, instability)
        
        return min(max_instability / 10, 1.0)
    
    def _repeat_score(self, sequence: str) -> float:
        """
        Direct repeat scoring based on repeat characteristics
        """
        if len(sequence) < 10:
            return 0.0
        
        # Look for direct repeats with spacing
        max_score = 0
        
        # Try different repeat unit sizes
        for unit_size in range(5, min(51, len(sequence) // 2)):
            for i in range(len(sequence) - unit_size):
                unit = sequence[i:i + unit_size]
                
                # Look for the same unit later in the sequence
                for j in range(i + unit_size, len(sequence) - unit_size + 1):
                    if sequence[j:j + unit_size] == unit:
                        # Calculate score based on unit size and spacing
                        spacing = j - i - unit_size
                        score = unit_size / (1 + spacing / 100)  # Penalty for large spacing
                        max_score = max(max_score, score)
        
        return min(max_score / 50, 1.0)  # Normalize to 0-1