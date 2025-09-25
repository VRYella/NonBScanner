"""
i-Motif Scoring Module
====================

Implements i-Motif scoring using G4Hunter-style algorithm adapted for C-tracts.
Provides both canonical and AC-motif scoring approaches.
"""

from .g4_scorer import iMotifScorer

# Re-export for convenience
IMotifScorer = iMotifScorer