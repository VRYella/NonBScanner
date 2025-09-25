"""
NBDFinder Centralized Scoring Module
===================================

Centralized scoring system that accepts candidate motif lists and assigns
both raw and normalized scores based on scientific consensus/range.

This module separates scoring logic from detection logic to allow for:
- Independent scoring of candidate motifs from any detection method
- Standardized score normalization across all motif types
- Easy addition of new scoring algorithms
- Scientific reproducibility through documented scoring methods

Main Classes:
- ScoreCalculator: Base class for all scoring algorithms
- MotifScorer: Main interface for scoring candidate motifs
"""

from .base_scorer import ScoreCalculator, MotifScorer
from .zdna_scorer import ZDNAScorer, eGZScorer
from .g4_scorer import G4Scorer, iMotifScorer
from .imotif_scorer import IMotifScorer
from .rloop_scorer import RLoopScorer

# Create default scorer instance
def create_default_scorer():
    """Create a MotifScorer with default scoring algorithms registered."""
    scorer = MotifScorer()
    
    # Register default scorers for each motif class
    scorer.register_scorer('Z-DNA', ZDNAScorer())
    scorer.register_scorer('eGZ', eGZScorer())
    scorer.register_scorer('G-Quadruplex', G4Scorer())
    scorer.register_scorer('i-Motif', iMotifScorer())
    scorer.register_scorer('R-Loop', RLoopScorer())
    
    return scorer

__all__ = [
    'ScoreCalculator',
    'MotifScorer', 
    'ZDNAScorer',
    'eGZScorer',
    'G4Scorer',
    'iMotifScorer',
    'IMotifScorer',
    'RLoopScorer',
    'create_default_scorer'
]