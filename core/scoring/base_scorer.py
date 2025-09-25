"""
Base Scoring Framework for NBDFinder
===================================

Provides base classes and interfaces for all scoring algorithms.
Ensures consistent scoring interface and standardized normalization.
"""

from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional, Union
import numpy as np


class ScoreCalculator(ABC):
    """
    Abstract base class for all scoring algorithms.
    
    Each scoring algorithm must implement:
    - calculate_raw_score: Returns the raw scientific score
    - calculate_normalized_score: Returns score normalized to [0,1]
    - get_score_method: Returns string identifier of scoring method
    """
    
    @abstractmethod
    def calculate_raw_score(self, candidate: Dict[str, Any], sequence: str) -> float:
        """
        Calculate raw scientific score for a candidate motif.
        
        Args:
            candidate: Candidate motif dict with start, end, sequence, class, subclass
            sequence: Full DNA sequence for context if needed
            
        Returns:
            Raw score (algorithm-specific scale)
        """
        pass
    
    @abstractmethod 
    def calculate_normalized_score(self, raw_score: float, candidate: Dict[str, Any]) -> float:
        """
        Normalize raw score to [0,1] range based on scientific consensus.
        
        Args:
            raw_score: Raw score from calculate_raw_score
            candidate: Candidate motif for context-dependent normalization
            
        Returns:
            Normalized score in [0,1] range
        """
        pass
    
    @abstractmethod
    def get_score_method(self) -> str:
        """Return string identifier for the scoring method."""
        pass
    
    def score_candidate(self, candidate: Dict[str, Any], sequence: str) -> Dict[str, Any]:
        """
        Score a single candidate and return enhanced candidate with scores.
        
        Args:
            candidate: Candidate motif dict
            sequence: Full DNA sequence
            
        Returns:
            Enhanced candidate with raw_score, normalized_score, score_method
        """
        raw_score = self.calculate_raw_score(candidate, sequence)
        normalized_score = self.calculate_normalized_score(raw_score, candidate)
        
        enhanced_candidate = candidate.copy()
        enhanced_candidate.update({
            'Raw_Score': raw_score,
            'Normalized_Score': normalized_score,
            'Score_Method': self.get_score_method()
        })
        
        return enhanced_candidate


class MotifScorer:
    """
    Main interface for scoring candidate motifs.
    
    Manages multiple scoring algorithms and applies appropriate scorer
    based on motif class/subclass.
    """
    
    def __init__(self):
        self.scorers = {}
        self._register_default_scorers()
    
    def _register_default_scorers(self):
        """Register default scoring algorithms for each motif class."""
        # Will be populated as we create specific scorers
        pass
    
    def register_scorer(self, motif_class: str, scorer: ScoreCalculator):
        """
        Register a scoring algorithm for a motif class.
        
        Args:
            motif_class: Motif class identifier (e.g., "Z-DNA", "G-Quadruplex")
            scorer: ScoreCalculator instance
        """
        self.scorers[motif_class] = scorer
    
    def score_candidates(self, candidates: List[Dict[str, Any]], sequence: str) -> List[Dict[str, Any]]:
        """
        Score a list of candidate motifs.
        
        Args:
            candidates: List of candidate motif dicts
            sequence: Full DNA sequence
            
        Returns:
            List of enhanced candidates with scores
        """
        scored_candidates = []
        
        for candidate in candidates:
            motif_class = candidate.get('Class', 'Unknown')
            
            if motif_class in self.scorers:
                scored = self.scorers[motif_class].score_candidate(candidate, sequence)
            else:
                # Default scoring - preserve original if it has scores
                scored = candidate.copy()
                if 'Raw_Score' not in scored:
                    scored['Raw_Score'] = 0.0
                if 'Normalized_Score' not in scored:
                    scored['Normalized_Score'] = 0.0
                if 'Score_Method' not in scored:
                    scored['Score_Method'] = 'Default'
                    
            scored_candidates.append(scored)
        
        return scored_candidates
    
    def get_available_scorers(self) -> List[str]:
        """Return list of available scoring algorithm names."""
        return list(self.scorers.keys())


def normalize_score_linear(raw_score: float, min_val: float, max_val: float) -> float:
    """
    Linear normalization to [0,1] range.
    
    Args:
        raw_score: Raw score value
        min_val: Minimum expected score (maps to 0)
        max_val: Maximum expected score (maps to 1)
        
    Returns:
        Normalized score in [0,1], clamped to range
    """
    if max_val <= min_val:
        return 0.0
    
    normalized = (raw_score - min_val) / (max_val - min_val)
    return max(0.0, min(1.0, normalized))


def normalize_score_sigmoid(raw_score: float, midpoint: float, steepness: float = 1.0) -> float:
    """
    Sigmoid normalization to [0,1] range.
    
    Args:
        raw_score: Raw score value
        midpoint: Score value that maps to 0.5
        steepness: Controls sigmoid steepness (higher = steeper)
        
    Returns:
        Normalized score in [0,1]
    """
    return 1.0 / (1.0 + np.exp(-steepness * (raw_score - midpoint)))