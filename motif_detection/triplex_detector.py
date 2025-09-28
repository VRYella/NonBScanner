"""
Triplex DNA Motif Detector
=========================

Detects three-stranded DNA structures including:
- Triplex-forming sequences
- Sticky DNA repeats

Based on Frank-Kamenetskii 1995 and Sakamoto 1999.
"""

import re
from typing import List, Dict, Any, Tuple
from .base_detector import BaseMotifDetector


class TriplexDetector(BaseMotifDetector):
    """Detector for triplex DNA motifs"""
    
    def get_motif_class_name(self) -> str:
        return "Triplex"
    
    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return triplex DNA patterns"""
        return {
            'triplex_forming_sequences': [
                (r'[GA]{15,}', 'TRX_5_1', 'Homopurine tract', 'Triplex', 15, 'triplex_potential', 0.90, 'H-DNA formation', 'Frank-Kamenetskii 1995'),
                (r'[CT]{15,}', 'TRX_5_2', 'Homopyrimidine tract', 'Triplex', 15, 'triplex_potential', 0.90, 'H-DNA formation', 'Frank-Kamenetskii 1995'),
                (r'(?:GA){6,}[GA]*(?:TC){6,}', 'TRX_5_3', 'Mirror repeat', 'Triplex', 24, 'triplex_potential', 0.85, 'Intermolecular triplex', 'Beal 1996'),
                (r'(?:GAA){4,}', 'TRX_5_4', 'GAA repeat', 'Sticky DNA', 12, 'sticky_dna_score', 0.95, 'Disease-associated repeats', 'Sakamoto 1999'),
                (r'(?:TTC){4,}', 'TRX_5_5', 'TTC repeat', 'Sticky DNA', 12, 'sticky_dna_score', 0.95, 'Disease-associated repeats', 'Sakamoto 1999'),
            ]
        }
    
    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        """Calculate triplex formation score for the sequence"""
        scoring_method = pattern_info[5] if len(pattern_info) > 5 else 'triplex_potential'
        
        if scoring_method == 'triplex_potential':
            return self._triplex_potential(sequence)
        elif scoring_method == 'sticky_dna_score':
            return self._sticky_dna_score(sequence)
        else:
            return self._triplex_potential(sequence)
    
    def _triplex_potential(self, sequence: str) -> float:
        """
        Triplex formation potential (Frank-Kamenetskii & Mirkin, 1995)
        """
        if len(sequence) < 15:
            return 0.0
        
        # Calculate purine content (G+A)
        purine_content = len(re.findall(r'[GA]', sequence)) / len(sequence)
        
        # Calculate pyrimidine content (C+T)  
        pyrimidine_content = len(re.findall(r'[CT]', sequence)) / len(sequence)
        
        # Find homopurine/homopyrimidine tracts
        purine_tracts = re.findall(r'[GA]{10,}', sequence)
        pyrimidine_tracts = re.findall(r'[CT]{10,}', sequence)
        
        # Score based on homogeneity and tract length
        purine_score = sum(len(tract) ** 1.2 for tract in purine_tracts) / len(sequence)
        pyrimidine_score = sum(len(tract) ** 1.2 for tract in pyrimidine_tracts) / len(sequence)
        
        max_score = max(purine_score, pyrimidine_score)
        return min(max_score / (len(sequence) ** 1.2), 1.0)
    
    def _sticky_dna_score(self, sequence: str) -> float:
        """
        Sticky DNA scoring for GAA/TTC repeats
        """
        if len(sequence) < 12:
            return 0.0
        
        # Count GAA and TTC repeats
        gaa_repeats = len(re.findall(r'GAA', sequence))
        ttc_repeats = len(re.findall(r'TTC', sequence))
        
        total_repeats = gaa_repeats + ttc_repeats
        repeat_density = total_repeats * 3 / len(sequence)  # Each repeat is 3 bases
        
        # Bonus for consecutive repeats
        gaa_consecutive = re.findall(r'(?:GAA){2,}', sequence)
        ttc_consecutive = re.findall(r'(?:TTC){2,}', sequence)
        
        consecutive_bonus = sum(len(match) for match in gaa_consecutive + ttc_consecutive) / len(sequence)
        
        return min(repeat_density + consecutive_bonus * 0.3, 1.0)