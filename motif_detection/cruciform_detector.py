"""
Cruciform DNA Motif Detector
===========================

Detects inverted repeat-induced four-way junctions including:
- Potential palindromes
- Inverted repeat candidates

Based on Lilley 2000 and Sinden 1994.
"""

import re
from typing import List, Dict, Any, Tuple
from .base_detector import BaseMotifDetector


class CruciformDetector(BaseMotifDetector):
    """Detector for cruciform DNA motifs"""
    
    def get_motif_class_name(self) -> str:
        return "Cruciform"
    
    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return cruciform DNA patterns"""
        return {
            'inverted_repeats': [
                (r'([ATGC]{6,20})[ATGC]{0,50}', 'CRU_3_1', 'Potential palindrome', 'Inverted Repeats', 12, 'cruciform_stability', 0.95, 'DNA secondary structure', 'Lilley 2000'),
                (r'([ATGC]{8,15})[ATGC]{2,20}([ATGC]{8,15})', 'CRU_3_2', 'Inverted repeat candidate', 'Inverted Repeats', 16, 'cruciform_stability', 0.80, 'Secondary structure prone', 'Pearson 1996'),
                (r'([ATGC]{4,10})[ATGC]{0,10}([ATGC]{4,10})', 'CRU_3_3', 'Short inverted repeat', 'Inverted Repeats', 8, 'cruciform_stability', 0.70, 'Local secondary structure', 'Sinden 1994'),
            ]
        }
    
    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        """Calculate cruciform stability score for the sequence"""
        scoring_method = pattern_info[5] if len(pattern_info) > 5 else 'cruciform_stability'
        
        if scoring_method == 'cruciform_stability':
            return self._cruciform_stability(sequence)
        else:
            return self._cruciform_stability(sequence)
    
    def _cruciform_stability(self, sequence: str) -> float:
        """
        Cruciform stability based on palindrome characteristics
        """
        if len(sequence) < 8:
            return 0.0
        
        def reverse_complement(seq):
            complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            return ''.join(complement.get(base, base) for base in reversed(seq))
        
        max_palindrome = 0
        
        # Look for palindromic regions
        for i in range(len(sequence)):
            for j in range(i + 8, len(sequence) + 1):
                subseq = sequence[i:j]
                if subseq == reverse_complement(subseq):
                    palindrome_length = len(subseq)
                    stability = palindrome_length ** 1.5 / len(sequence)
                    max_palindrome = max(max_palindrome, stability)
        
        return min(max_palindrome, 1.0)
    
    def passes_quality_threshold(self, sequence: str, score: float, pattern_info: Tuple) -> bool:
        """Enhanced quality check for cruciform structures"""
        # Basic threshold check
        if not super().passes_quality_threshold(sequence, score, pattern_info):
            return False
        
        # Additional check: ensure minimum length
        min_length = pattern_info[4] if len(pattern_info) > 4 else 8
        if len(sequence) < min_length:
            return False
        
        # Check for some degree of palindromicity
        return self._has_palindromic_potential(sequence)
    
    def _has_palindromic_potential(self, sequence: str) -> bool:
        """Check if sequence has potential for palindrome formation"""
        def reverse_complement(seq):
            complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            return ''.join(complement.get(base, base) for base in reversed(seq))
        
        # Check for partial palindromes (at least 50% match)
        rev_comp = reverse_complement(sequence)
        matches = sum(1 for a, b in zip(sequence, rev_comp) if a == b)
        
        return matches / len(sequence) >= 0.5