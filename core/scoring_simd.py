"""
Vectorized SIMD Scoring Functions for NBDFinder
==============================================

High-performance vectorized implementations of scientific scoring algorithms:
- G4Hunter: G-richness and C-richness balance scoring
- i-Motif: C-richness scoring adapted from G4Hunter  
- Z-DNA Seeker: Dinucleotide scoring with penalties

Uses NumPy vectorization and optional Numba JIT compilation for performance.
"""

import numpy as np
from typing import List, Tuple, Optional, Union
from numba import jit, prange
import warnings

try:
    from numba import jit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    # Fallback decorator that does nothing
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return decorator

@jit(nopython=True) if NUMBA_AVAILABLE else jit()
def g4hunter_score_vectorized(sequence: str, window_size: int = 25) -> np.ndarray:
    """
    Vectorized G4Hunter scoring algorithm.
    
    Args:
        sequence: DNA sequence string
        window_size: Sliding window size for scoring
        
    Returns:
        Array of G4Hunter scores for each position
    """
    seq_len = len(sequence)
    if seq_len < window_size:
        return np.zeros(1)
    
    scores = np.zeros(seq_len - window_size + 1)
    
    for i in range(seq_len - window_size + 1):
        window = sequence[i:i + window_size]
        
        # Count G and C runs
        g_score = 0.0
        c_score = 0.0
        
        # G-run scoring
        current_g_run = 0
        for j in range(window_size):
            if window[j] == 'G':
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
        for j in range(window_size):
            if window[j] == 'C':
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
            scores[i] = 0.0  # Balanced G/C cancels out
        elif g_score > 0:
            scores[i] = g_score / window_size
        elif c_score > 0:
            scores[i] = -c_score / window_size
        else:
            scores[i] = 0.0
            
    return scores

@jit(nopython=True) if NUMBA_AVAILABLE else jit()
def imotif_score_vectorized(sequence: str, window_size: int = 25) -> np.ndarray:
    """
    Vectorized i-Motif scoring (C-richness adaptation of G4Hunter).
    
    Args:
        sequence: DNA sequence string
        window_size: Sliding window size for scoring
        
    Returns:
        Array of i-Motif scores for each position
    """
    seq_len = len(sequence)
    if seq_len < window_size:
        return np.zeros(1)
    
    scores = np.zeros(seq_len - window_size + 1)
    
    for i in range(seq_len - window_size + 1):
        window = sequence[i:i + window_size]
        
        # Count C runs (primary) and G runs (interference)
        c_score = 0.0
        g_score = 0.0
        
        # C-run scoring (positive for i-motifs)
        current_c_run = 0
        for j in range(window_size):
            if window[j] == 'C':
                current_c_run += 1
            else:
                if current_c_run >= 2:
                    c_score += current_c_run ** 2
                current_c_run = 0
        # Final C run
        if current_c_run >= 2:
            c_score += current_c_run ** 2
            
        # G-run scoring (interference)
        current_g_run = 0
        for j in range(window_size):
            if window[j] == 'G':
                current_g_run += 1
            else:
                if current_g_run >= 2:
                    g_score += current_g_run ** 2
                current_g_run = 0
        # Final G run
        if current_g_run >= 2:
            g_score += current_g_run ** 2
            
        # Calculate final score (C-runs positive, G-runs reduce score)
        if c_score > 0:
            interference_factor = 1.0 - (g_score / (2 * window_size))
            scores[i] = max(0.0, (c_score / window_size) * interference_factor)
        else:
            scores[i] = 0.0
            
    return scores

# Z-DNA dinucleotide scoring weights
ZDNA_WEIGHTS = {
    'CG': 1.0, 'GC': 1.0,
    'CA': 0.5, 'AC': 0.5, 'TG': 0.5, 'GT': 0.5,
    'CC': 0.3, 'GG': 0.3,
    'CT': 0.2, 'TC': 0.2, 'AG': 0.2, 'GA': 0.2,
    'AT': -0.1, 'TA': -0.1,
    'AA': -0.2, 'TT': -0.2,
}

@jit(nopython=True) if NUMBA_AVAILABLE else jit()
def zdna_score_vectorized(sequence: str, window_size: int = 50) -> np.ndarray:
    """
    Vectorized Z-DNA Seeker scoring algorithm.
    
    Args:
        sequence: DNA sequence string
        window_size: Sliding window size for scoring
        
    Returns:
        Array of Z-DNA scores for each position
    """
    seq_len = len(sequence)
    if seq_len < window_size:
        return np.zeros(1)
    
    scores = np.zeros(seq_len - window_size + 1)
    
    # Create lookup for fast dinucleotide scoring
    dinuc_scores = np.zeros(256)  # ASCII lookup table
    dinuc_scores[ord('C')*16 + ord('G')] = 1.0
    dinuc_scores[ord('G')*16 + ord('C')] = 1.0
    dinuc_scores[ord('C')*16 + ord('A')] = 0.5
    dinuc_scores[ord('A')*16 + ord('C')] = 0.5
    dinuc_scores[ord('T')*16 + ord('G')] = 0.5
    dinuc_scores[ord('G')*16 + ord('T')] = 0.5
    dinuc_scores[ord('C')*16 + ord('C')] = 0.3
    dinuc_scores[ord('G')*16 + ord('G')] = 0.3
    dinuc_scores[ord('C')*16 + ord('T')] = 0.2
    dinuc_scores[ord('T')*16 + ord('C')] = 0.2
    dinuc_scores[ord('A')*16 + ord('G')] = 0.2
    dinuc_scores[ord('G')*16 + ord('A')] = 0.2
    dinuc_scores[ord('A')*16 + ord('T')] = -0.1
    dinuc_scores[ord('T')*16 + ord('A')] = -0.1
    dinuc_scores[ord('A')*16 + ord('A')] = -0.2
    dinuc_scores[ord('T')*16 + ord('T')] = -0.2
    
    for i in range(seq_len - window_size + 1):
        window = sequence[i:i + window_size]
        score = 0.0
        
        # Score all dinucleotides in window
        for j in range(window_size - 1):
            char1 = ord(window[j])
            char2 = ord(window[j + 1])
            score += dinuc_scores[char1 * 16 + char2]
            
        scores[i] = score / (window_size - 1)  # Normalize by dinucleotide count
        
    return scores

def score_sequence_region(sequence: str, start: int, end: int, 
                         scoring_method: str = "g4hunter",
                         **kwargs) -> float:
    """
    Score a specific region of sequence using vectorized algorithms.
    
    Args:
        sequence: Full DNA sequence
        start: Start position (0-based)
        end: End position (0-based, exclusive)
        scoring_method: Scoring algorithm to use
        **kwargs: Additional parameters for scoring functions
        
    Returns:
        Mean score for the region
    """
    region = sequence[start:end].upper()
    
    if scoring_method == "g4hunter":
        scores = g4hunter_score_vectorized(region, kwargs.get('window_size', 25))
    elif scoring_method == "imotif":
        scores = imotif_score_vectorized(region, kwargs.get('window_size', 25))
    elif scoring_method == "zdna":
        scores = zdna_score_vectorized(region, kwargs.get('window_size', 50))
    else:
        raise ValueError(f"Unknown scoring method: {scoring_method}")
    
    return float(np.mean(scores)) if len(scores) > 0 else 0.0

def batch_score_regions(sequence: str, regions: List[Tuple[int, int]], 
                       scoring_method: str = "g4hunter",
                       **kwargs) -> List[float]:
    """
    Score multiple sequence regions efficiently.
    
    Args:
        sequence: Full DNA sequence
        regions: List of (start, end) tuples
        scoring_method: Scoring algorithm to use
        **kwargs: Additional parameters for scoring functions
        
    Returns:
        List of mean scores for each region
    """
    return [score_sequence_region(sequence, start, end, scoring_method, **kwargs)
            for start, end in regions]

__all__ = [
    'g4hunter_score_vectorized',
    'imotif_score_vectorized', 
    'zdna_score_vectorized',
    'score_sequence_region',
    'batch_score_regions',
    'ZDNA_WEIGHTS'
]