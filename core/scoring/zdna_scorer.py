"""
Z-DNA Scoring Module
==================

Implements Z-DNA scoring using the scientifically validated Z-seeker algorithm
(Ho et al. 1986, Wang et al. 2007) with dinucleotide-based scoring.

Preserves the original z-seeker approach while providing normalized scores
for cross-motif comparison.
"""

import numpy as np
from typing import Dict, Any
from .base_scorer import ScoreCalculator, normalize_score_linear


class ZDNAScorer(ScoreCalculator):
    """
    Z-DNA scorer implementing Z-seeker algorithm.
    
    Uses experimentally validated dinucleotide weights:
    - CG/GC: +7.0 (strong Z-forming potential)
    - AT/TA: +0.5 (weak Z-forming potential)  
    - GT/TG, AC/CA: +1.25 (moderate Z-forming potential)
    - Consecutive AT penalty to avoid false positives
    """
    
    def __init__(self, 
                 GC_weight: float = 7.0, 
                 AT_weight: float = 0.5, 
                 GT_weight: float = 1.25, 
                 AC_weight: float = 1.25,
                 consecutive_AT_scoring: tuple = (0.5, 0.5, 0.5, 0.5, 0.0, 0.0, -5.0, -100.0),
                 mismatch_penalty_type: str = "linear",
                 mismatch_penalty_starting_value: int = 3,
                 mismatch_penalty_linear_delta: int = 3,
                 cadence_reward: float = 0.0):
        """
        Initialize Z-DNA scorer with Z-seeker parameters.
        
        Args:
            GC_weight: Weight for CG/GC dinucleotides
            AT_weight: Weight for AT/TA dinucleotides  
            GT_weight: Weight for GT/TG dinucleotides
            AC_weight: Weight for AC/CA dinucleotides
            consecutive_AT_scoring: Penalty sequence for consecutive AT
            mismatch_penalty_type: Type of mismatch penalty ("linear" or "exponential")
            mismatch_penalty_starting_value: Starting penalty value
            mismatch_penalty_linear_delta: Linear penalty increment
            cadence_reward: Bonus for valid dinucleotides
        """
        self.GC_weight = GC_weight
        self.AT_weight = AT_weight
        self.GT_weight = GT_weight
        self.AC_weight = AC_weight
        self.consecutive_AT_scoring = consecutive_AT_scoring
        self.mismatch_penalty_type = mismatch_penalty_type
        self.mismatch_penalty_starting_value = mismatch_penalty_starting_value
        self.mismatch_penalty_linear_delta = mismatch_penalty_linear_delta
        self.cadence_reward = cadence_reward
    
    def zdna_seeker_scoring_array(self, seq: str) -> np.ndarray:
        """
        Calculate Z-seeker scoring array for sequence.
        
        Args:
            seq: DNA sequence
            
        Returns:
            Array of Z-seeker scores for each position
        """
        seq = seq.upper()
        n = len(seq)
        if n < 2:
            return np.array([0.0])
            
        scoring_array = np.zeros(n-1)
        consecutive_AT_counter = 0
        mismatches_counter = 0
        
        for i in range(n-1):
            t = seq[i:i+2]
            
            if t in ("GC", "CG"):
                scoring_array[i] = self.GC_weight
                mismatches_counter = 0
                consecutive_AT_counter = 0
            elif t in ("GT", "TG"):
                scoring_array[i] = self.GT_weight
                mismatches_counter = 0
                consecutive_AT_counter = 0
            elif t in ("AC", "CA"):
                scoring_array[i] = self.AC_weight
                mismatches_counter = 0
                consecutive_AT_counter = 0
            elif t in ("AT", "TA"):
                adjusted_weight = self.AT_weight
                adjusted_weight += self.consecutive_AT_scoring[consecutive_AT_counter] if consecutive_AT_counter < len(self.consecutive_AT_scoring) else self.consecutive_AT_scoring[-1]
                scoring_array[i] = adjusted_weight
                consecutive_AT_counter += 1
                mismatches_counter = 0
            else:
                mismatches_counter += 1
                consecutive_AT_counter = 0
                if self.mismatch_penalty_type == "exponential":
                    scoring_array[i] = -(self.mismatch_penalty_starting_value**mismatches_counter if mismatches_counter < 15 else 32000.0)
                elif self.mismatch_penalty_type == "linear":
                    scoring_array[i] = -self.mismatch_penalty_starting_value - self.mismatch_penalty_linear_delta * (mismatches_counter - 1)
                else:
                    scoring_array[i] = -10.0
                    
            if t in ("GC", "CG", "GT", "TG", "AC", "CA", "AT", "TA"):
                scoring_array[i] += self.cadence_reward
                
        return scoring_array
    
    def calculate_raw_score(self, candidate: Dict[str, Any], sequence: str) -> float:
        """
        Calculate raw Z-seeker score for candidate region.
        
        Args:
            candidate: Candidate with Start, End, Sequence
            sequence: Full DNA sequence
            
        Returns:
            Raw Z-seeker score (sum of dinucleotide scores)
        """
        # Extract candidate sequence
        start = candidate.get('Start', 1) - 1  # Convert to 0-based
        end = candidate.get('End', 1)
        
        if 'Sequence' in candidate:
            candidate_seq = candidate['Sequence']
        else:
            candidate_seq = sequence[start:end]
        
        # Calculate Z-seeker scores
        scores = self.zdna_seeker_scoring_array(candidate_seq)
        
        # Return sum as raw score
        return float(np.sum(scores))
    
    def calculate_normalized_score(self, raw_score: float, candidate: Dict[str, Any]) -> float:
        """
        Normalize Z-seeker score to [0,1] range.
        
        Based on empirical Z-seeker score ranges:
        - Minimum: ~0 (no Z-forming potential)
        - Maximum: ~length * 7.0 (perfect CG dinucleotides)
        
        Args:
            raw_score: Raw Z-seeker score
            candidate: Candidate for length context
            
        Returns:
            Normalized score in [0,1]
        """
        length = candidate.get('Length', 1)
        if length <= 0:
            return 0.0
            
        # Theoretical maximum: all CG dinucleotides
        max_possible = length * self.GC_weight
        
        # Normalize using linear scaling
        return normalize_score_linear(raw_score, 0.0, max_possible)
    
    def get_score_method(self) -> str:
        """Return scoring method identifier."""
        return "Z-Seeker"


class eGZScorer(ScoreCalculator):
    """
    Scorer for eGZ (Extruded-G) DNA motifs.
    
    Scores CGG repeat expansions based on repeat count and G-fraction.
    """
    
    def calculate_raw_score(self, candidate: Dict[str, Any], sequence: str) -> float:
        """
        Calculate eGZ score based on CGG repeats and G-content.
        
        Args:
            candidate: Candidate with Sequence
            sequence: Full DNA sequence (unused)
            
        Returns:
            Raw eGZ score
        """
        candidate_seq = candidate.get('Sequence', '')
        if not candidate_seq:
            return 0.0
            
        # Count CGG repeats
        n_repeats = len(candidate_seq) // 3
        g_frac = candidate_seq.count('G') / len(candidate_seq) if len(candidate_seq) > 0 else 0
        
        # Score formula: repeat_count * 3 * (1 + 2*G_fraction)
        score = n_repeats * 3 * (1.0 + 2.0 * g_frac)
        
        return float(score)
    
    def calculate_normalized_score(self, raw_score: float, candidate: Dict[str, Any]) -> float:
        """
        Normalize eGZ score to [0,1] range.
        
        Args:
            raw_score: Raw eGZ score
            candidate: Candidate for context
            
        Returns:
            Normalized score in [0,1]
        """
        length = candidate.get('Length', 1)
        if length <= 0:
            return 0.0
            
        # Theoretical maximum for perfect CGG repeats: length * 3
        max_possible = length * 3
        
        return normalize_score_linear(raw_score, 0.0, max_possible)
    
    def get_score_method(self) -> str:
        """Return scoring method identifier."""
        return "eGZ-Repeat"