"""
A-philic DNA Motif Detection (Class 9) - NBDFinder Detector Module
================================================================

This module provides A-philic DNA detection as part of the NBDFinder detector system.
It wraps the existing A-philic DNA implementation from motifs.a_philic_dna to provide
a consistent interface with other detector modules.

A-philic DNA represents DNA sequences with high affinity for specific protein binding
and structural features that favor A-tract formation and protein-DNA interactions.

Scientific References:
- A-tract structural properties: Bolshoy et al. PNAS 1991
- Protein-DNA interactions: Rohs et al. Nature 2009
- Tetranucleotide analysis: Vinogradov Bioinformatics 2003

Author: Dr. Venkata Rajesh Yella
Integration: 2024 - NBDFinder Class 9 Implementation
"""

import sys
import os

# Add current directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from motifs.a_philic_dna import find_a_philic_dna
except ImportError:
    def find_a_philic_dna(sequence: str, sequence_name: str = ""):
        """Fallback function if main A-philic DNA module is not available."""
        return []

def find_a_philic(sequence: str, sequence_name: str = "") -> list:
    """
    Detect A-philic DNA motifs in a sequence.
    
    This function provides a consistent interface for A-philic DNA detection
    within the NBDFinder detector system.
    
    Args:
        sequence: DNA sequence to analyze
        sequence_name: Name/identifier for the sequence
    
    Returns:
        List of A-philic DNA motif dictionaries
    """
    return find_a_philic_dna(sequence, sequence_name)