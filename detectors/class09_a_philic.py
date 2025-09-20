"""
A-philic DNA Detector (Class 9)
==============================

A-philic DNA motif detector for the NBDFinder system.
This detector identifies DNA sequences with high affinity for A-tract formation
and specific protein-DNA interactions.

This module implements Class 9 detection in the 11-class NBDFinder system:
Classes 1-8: Original motif types
Class 9: A-philic DNA (NEW)
Class 10: Hybrid motifs 
Class 11: Clusters

Author: Dr. Venkata Rajesh Yella
Integration: 2024 - NBDFinder Detector Module
"""

import sys
import os
from typing import List, Dict, Any

# Add paths for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from motifs.a_philic_dna import find_a_philic_dna
except ImportError:
    def find_a_philic_dna(sequence: str, sequence_name: str = "Unknown") -> List[Dict[str, Any]]:
        """Fallback function if motifs.a_philic_dna cannot be imported"""
        return []


def find_a_philic_DNA(sequence: str, sequence_name: str = "Unknown") -> List[Dict[str, Any]]:
    """
    Main detector function for A-philic DNA motifs (Class 9).
    
    This function serves as the standardized interface for the detector registry
    and calls the underlying A-philic DNA detection algorithm.
    
    Args:
        sequence: DNA sequence to analyze
        sequence_name: Name/identifier for the sequence
        
    Returns:
        List of A-philic DNA motif dictionaries with standardized format
    """
    return find_a_philic_dna(sequence, sequence_name)


# Export the main detection function
__all__ = ['find_a_philic_DNA']


if __name__ == "__main__":
    # Test the detector
    test_seq = "NNNNAGGGGGGGGGCCCCTGGGGGCCCAAGGGNNNN"
    motifs = find_a_philic_DNA(test_seq, "test_sequence")
    
    print(f"A-philic DNA Detector Test:")
    print(f"Found {len(motifs)} motifs in test sequence")
    for motif in motifs:
        print(f"  {motif.get('Subclass', 'Unknown')}: {motif.get('Start', 0)}-{motif.get('End', 0)} (score: {motif.get('Score', 0)})")