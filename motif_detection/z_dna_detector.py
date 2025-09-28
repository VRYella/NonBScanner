"""
Z-DNA Motif Detector
===================

Detects left-handed helix structures including:
- Z-DNA canonical
- eGZ DNA (Extruded-G Z-DNA)

Based on Ho 1986 and Z-DNA transition research.
"""

import re
from typing import List, Dict, Any, Tuple
from .base_detector import BaseMotifDetector


class ZDNADetector(BaseMotifDetector):
    """Detector for Z-DNA motifs"""
    
    def get_motif_class_name(self) -> str:
        return "Z-DNA"
    
    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return Z-DNA patterns"""
        return {
            'z_dna_canonical': [
                (r'(?:CG){4,}', 'Z_8_1', 'CG alternating', 'Z-DNA canonical', 8, 'z_dna_score', 0.95, 'Classical Z-DNA', 'Wang 1979'),
                (r'(?:GC){4,}', 'Z_8_2', 'GC alternating', 'Z-DNA canonical', 8, 'z_dna_score', 0.95, 'Classical Z-DNA', 'Rich 1984'),
                (r'(?:[CG][ATGC]){6,}', 'Z_8_3', 'Mixed purine-pyrimidine', 'Z-DNA canonical', 12, 'z_dna_score', 0.85, 'Z-forming potential', 'Nordheim 1981'),
            ],
            'egz_dna': [
                (r'G(?:CG){3,}', 'Z_8_4', 'eGZ motif', 'eGZ DNA', 7, 'z_dna_score', 0.90, 'Extruded-G Z-DNA', 'Schroth 1992'),
                (r'(?:CG){2,}[ATGC](?:CG){2,}', 'Z_8_5', 'Interrupted Z-DNA', 'eGZ DNA', 10, 'z_dna_score', 0.80, 'Z-DNA with interruption', 'Herbert 1996'),
            ]
        }
    
    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        """Calculate Z-DNA formation score for the sequence"""
        scoring_method = pattern_info[5] if len(pattern_info) > 5 else 'z_dna_score'
        
        if scoring_method == 'z_dna_score':
            return self._z_dna_score(sequence)
        else:
            return self._z_dna_score(sequence)
    
    def _z_dna_score(self, sequence: str) -> float:
        """
        Z-DNA scoring based on alternating purine-pyrimidine content (Ho et al., 1986)
        """
        if len(sequence) < 6:
            return 0.0
        
        alternating_score = 0
        total_pairs = 0
        
        for i in range(len(sequence) - 1):
            curr_base = sequence[i]
            next_base = sequence[i + 1]
            total_pairs += 1
            
            # Score alternating purine-pyrimidine pattern
            if ((curr_base in 'AG' and next_base in 'CT') or 
                (curr_base in 'CT' and next_base in 'AG')):
                alternating_score += 1
                
                # Bonus for CG steps (classical Z-DNA)
                if (curr_base == 'C' and next_base == 'G') or (curr_base == 'G' and next_base == 'C'):
                    alternating_score += 0.5
        
        return alternating_score / total_pairs if total_pairs > 0 else 0.0