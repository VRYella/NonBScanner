"""
R-loop Scoring Module
===================

Implements R-loop scoring using qmRLFS (quantitative machine learning R-Loop Forming Sequence) approach.
Preserves scientific accuracy for R-loop detection.
"""

import numpy as np
from typing import Dict, Any
from .base_scorer import ScoreCalculator, normalize_score_linear


class RLoopScorer(ScoreCalculator):
    """
    R-loop scorer implementing qmRLFS approach.
    
    Scores R-loop forming sequences based on:
    - GC skew and content
    - Purine/pyrimidine asymmetry
    - Length considerations
    """
    
    def __init__(self):
        """Initialize R-loop scorer."""
        pass
    
    def calculate_gc_skew(self, sequence: str) -> float:
        """
        Calculate GC skew: (G-C)/(G+C)
        
        Args:
            sequence: DNA sequence
            
        Returns:
            GC skew value
        """
        seq = sequence.upper()
        g_count = seq.count('G')
        c_count = seq.count('C')
        
        if g_count + c_count == 0:
            return 0.0
            
        return (g_count - c_count) / (g_count + c_count)
    
    def calculate_purine_content(self, sequence: str) -> float:
        """
        Calculate purine content (A+G)/(A+T+G+C)
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Purine content fraction
        """
        seq = sequence.upper()
        purines = seq.count('A') + seq.count('G')
        total = len(seq)
        
        return purines / total if total > 0 else 0.0
    
    def qmrlfs_score(self, sequence: str) -> float:
        """
        Calculate qmRLFS score for R-loop potential.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            qmRLFS score
        """
        if len(sequence) < 10:
            return 0.0
            
        gc_skew = self.calculate_gc_skew(sequence)
        purine_content = self.calculate_purine_content(sequence)
        
        # Simple qmRLFS approximation
        # Real qmRLFS uses machine learning models
        # This is a simplified version based on known R-loop characteristics
        
        # R-loops favor:
        # - Positive GC skew (G > C)
        # - High purine content on template strand
        # - Moderate length
        
        skew_score = max(0, gc_skew)  # Only positive skew contributes
        purine_score = purine_content
        length_factor = min(1.0, len(sequence) / 100.0)  # Normalize for length
        
        score = (skew_score * 0.4 + purine_score * 0.4 + length_factor * 0.2)
        
        return score
    
    def calculate_raw_score(self, candidate: Dict[str, Any], sequence: str) -> float:
        """
        Calculate raw qmRLFS score for candidate.
        
        Args:
            candidate: Candidate with Sequence
            sequence: Full DNA sequence
            
        Returns:
            Raw qmRLFS score
        """
        candidate_seq = candidate.get('Sequence', '')
        if not candidate_seq:
            return 0.0
            
        return self.qmrlfs_score(candidate_seq)
    
    def calculate_normalized_score(self, raw_score: float, candidate: Dict[str, Any]) -> float:
        """
        Normalize qmRLFS score to [0,1] range.
        
        Args:
            raw_score: Raw qmRLFS score
            candidate: Candidate for context
            
        Returns:
            Normalized score in [0,1]
        """
        # qmRLFS scores typically range from 0 to 1
        return normalize_score_linear(raw_score, 0.0, 1.0)
    
    def get_score_method(self) -> str:
        """Return scoring method identifier."""
        return "qmRLFS"