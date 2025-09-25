"""
G-Quadruplex Scoring Module
=========================

Implements G4Hunter scoring algorithm for G-quadruplex motifs.
Preserves scientific accuracy while providing normalized scores.
"""

import numpy as np
from typing import Dict, Any
from .base_scorer import ScoreCalculator, normalize_score_linear


class G4Scorer(ScoreCalculator):
    """
    G-Quadruplex scorer implementing G4Hunter algorithm.
    
    Scores based on G-richness and C-richness balance with
    run-length squared scoring for consecutive G/C runs.
    """
    
    def __init__(self, window_size: int = 25):
        """
        Initialize G4 scorer.
        
        Args:
            window_size: Window size for G4Hunter scoring
        """
        self.window_size = window_size
    
    def g4hunter_score(self, sequence: str) -> float:
        """
        Calculate G4Hunter score for sequence.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            G4Hunter score
        """
        seq = sequence.upper()
        if len(seq) < 2:
            return 0.0
            
        # Count G and C runs
        g_score = 0.0
        c_score = 0.0
        
        # G-run scoring
        current_g_run = 0
        for char in seq:
            if char == 'G':
                current_g_run += 1
            else:
                if current_g_run >= 2:
                    g_score += current_g_run ** 2
                current_g_run = 0
        # Final G run
        if current_g_run >= 2:
            g_score += current_g_run ** 2
            
        # C-run scoring  
        current_c_run = 0
        for char in seq:
            if char == 'C':
                current_c_run += 1
            else:
                if current_c_run >= 2:
                    c_score += current_c_run ** 2
                current_c_run = 0
        # Final C run
        if current_c_run >= 2:
            c_score += current_c_run ** 2
            
        # Calculate final score
        if g_score > 0 and c_score > 0:
            return 0.0  # Balanced G/C cancels out
        elif g_score > 0:
            return g_score / len(seq)
        elif c_score > 0:
            return -c_score / len(seq)
        else:
            return 0.0
    
    def calculate_raw_score(self, candidate: Dict[str, Any], sequence: str) -> float:
        """
        Calculate raw G4Hunter score for candidate.
        
        Args:
            candidate: Candidate with Sequence
            sequence: Full DNA sequence
            
        Returns:
            Raw G4Hunter score
        """
        candidate_seq = candidate.get('Sequence', '')
        if not candidate_seq:
            return 0.0
            
        return self.g4hunter_score(candidate_seq)
    
    def calculate_normalized_score(self, raw_score: float, candidate: Dict[str, Any]) -> float:
        """
        Normalize G4Hunter score to [0,1] range.
        
        G4Hunter scores typically range from -2.0 to +2.0:
        - Positive scores indicate G4 potential
        - Negative scores indicate i-motif potential
        - We normalize positive scores to [0,1]
        
        Args:
            raw_score: Raw G4Hunter score
            candidate: Candidate for context
            
        Returns:
            Normalized score in [0,1]
        """
        # Only consider positive scores for G4 potential
        if raw_score <= 0:
            return 0.0
            
        # Typical G4Hunter range is 0 to 2.0 for G4s
        return normalize_score_linear(raw_score, 0.0, 2.0)
    
    def get_score_method(self) -> str:
        """Return scoring method identifier."""
        return "G4Hunter"


class iMotifScorer(ScoreCalculator):
    """
    i-Motif scorer adapted from G4Hunter for C-richness.
    """
    
    def __init__(self):
        """Initialize i-Motif scorer."""
        pass
    
    def imotif_score(self, sequence: str) -> float:
        """
        Calculate i-Motif score based on C-richness.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            i-Motif score (positive values indicate i-motif potential)
        """
        seq = sequence.upper()
        if len(seq) < 2:
            return 0.0
            
        # C-run scoring (adapted from G4Hunter)
        c_score = 0.0
        current_c_run = 0
        
        for char in seq:
            if char == 'C':
                current_c_run += 1
            else:
                if current_c_run >= 2:
                    c_score += current_c_run ** 2
                current_c_run = 0
        
        # Final C run
        if current_c_run >= 2:
            c_score += current_c_run ** 2
            
        return c_score / len(seq) if len(seq) > 0 else 0.0
    
    def calculate_raw_score(self, candidate: Dict[str, Any], sequence: str) -> float:
        """
        Calculate raw i-Motif score.
        
        Args:
            candidate: Candidate with Sequence
            sequence: Full DNA sequence
            
        Returns:
            Raw i-Motif score
        """
        candidate_seq = candidate.get('Sequence', '')
        if not candidate_seq:
            return 0.0
            
        return self.imotif_score(candidate_seq)
    
    def calculate_normalized_score(self, raw_score: float, candidate: Dict[str, Any]) -> float:
        """
        Normalize i-Motif score to [0,1] range.
        
        Args:
            raw_score: Raw i-Motif score
            candidate: Candidate for context
            
        Returns:
            Normalized score in [0,1]
        """
        # i-Motif scores typically range from 0 to 2.0
        return normalize_score_linear(raw_score, 0.0, 2.0)
    
    def get_score_method(self) -> str:
        """Return scoring method identifier."""
        return "iMotif-Hunter"