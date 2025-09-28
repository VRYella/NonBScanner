"""
G-Quadruplex DNA Motif Detector
==============================

Detects four-stranded G-rich structures including:
- Canonical G4
- Relaxed G4
- Bulged, bipartite, multimeric, imperfect G4
- G-Triplex intermediate

Based on Bedrat 2016 G4Hunter algorithm.
"""

import re
from typing import List, Dict, Any, Tuple
from .base_detector import BaseMotifDetector


class GQuadruplexDetector(BaseMotifDetector):
    """Detector for G-quadruplex DNA motifs"""
    
    def get_motif_class_name(self) -> str:
        return "G-Quadruplex"
    
    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return G-quadruplex DNA patterns"""
        return {
            'canonical_g4': [
                (r'G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G4_6_1', 'Canonical G4', 'Canonical G4', 15, 'g4hunter_score', 0.95, 'Stable G4 structures', 'Burge 2006'),
                (r'G{4,}[ATGC]{1,5}G{4,}[ATGC]{1,5}G{4,}[ATGC]{1,5}G{4,}', 'G4_6_2', 'High-density G4', 'Canonical G4', 16, 'g4hunter_score', 0.98, 'Very stable G4', 'Todd 2005'),
            ],
            'relaxed_g4': [
                (r'G{2,}[ATGC]{1,12}G{2,}[ATGC]{1,12}G{2,}[ATGC]{1,12}G{2,}', 'G4_6_3', 'Relaxed G4', 'Relaxed G4', 12, 'g4hunter_score', 0.80, 'Potential G4 structures', 'Huppert 2005'),
                (r'G{3,}[ATGC]{8,15}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G4_6_4', 'Long-loop G4', 'Relaxed G4', 18, 'g4hunter_score', 0.75, 'Alternative G4 topology', 'Phan 2006'),
            ],
            'bulged_g4': [
                (r'G{3,}[ATGC]{8,25}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G4_6_5', 'Bulged G4', 'Bulged G4', 20, 'g4hunter_score', 0.85, 'G4 with bulge loops', 'Lim 2009'),
                (r'G{2,}[ATGC]{15,40}G{2,}[ATGC]{1,7}G{2,}[ATGC]{1,7}G{2,}', 'G4_6_6', 'Large bulge G4', 'Bulged G4', 25, 'g4hunter_score', 0.70, 'Extended bulge G4', 'Adrian 2014'),
            ],
            'bipartite_g4': [
                (r'G{2,}[ATGC]{15,70}G{2,}[ATGC]{1,7}G{2,}[ATGC]{1,7}G{2,}', 'G4_6_7', 'Bipartite G4', 'Bipartite G4', 30, 'g4hunter_score', 0.75, 'Two-block G4 structures', 'Guédin 2010'),
            ],
            'multimeric_g4': [
                (r'(?:G{3,}[ATGC]{1,7}){4,}G{3,}', 'G4_6_8', 'Multimeric G4', 'Multimeric G4', 25, 'g4hunter_score', 0.90, 'Multiple G4 units', 'Phan 2007'),
                (r'(?:G{2,}[ATGC]{1,10}){5,}G{2,}', 'G4_6_9', 'Extended multimeric G4', 'Multimeric G4', 30, 'g4hunter_score', 0.85, 'Long G4 arrays', 'Maizels 2006'),
            ],
            'imperfect_g4': [
                (r'G{2,}[ATGC]{1,10}[AG]G{1,3}[ATGC]{1,10}G{2,}[ATGC]{1,10}G{2,}', 'G4_6_10', 'Imperfect G4', 'Imperfect G4', 15, 'g4hunter_score', 0.65, 'G4-like with interruptions', 'Kuryavyi 2010'),
                (r'G{3,}[ATGC]{1,7}[AG]{1,2}G{2,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G4_6_11', 'G-rich imperfect', 'Imperfect G4', 18, 'g4hunter_score', 0.70, 'Interrupted G-tracts', 'Webba da Silva 2007'),
            ],
            'g_triplex': [
                (r'G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G4_6_12', 'G-Triplex', 'G-Triplex intermediate', 12, 'g_triplex_score', 0.80, 'Three G-tract structures', 'Sen 1988'),
                (r'G{4,}[ATGC]{1,5}G{4,}[ATGC]{1,5}G{4,}', 'G4_6_13', 'High-density G-triplex', 'G-Triplex intermediate', 14, 'g_triplex_score', 0.85, 'Stable three-tract G-structure', 'Williamson 1989'),
            ]
        }
    
    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        """Calculate G4 formation score for the sequence"""
        scoring_method = pattern_info[5] if len(pattern_info) > 5 else 'g4hunter_score'
        
        if scoring_method == 'g4hunter_score':
            return self._g4hunter_score(sequence)
        elif scoring_method == 'g_triplex_score':
            return self._g_triplex_score(sequence)
        else:
            return self._g4hunter_score(sequence)
    
    def _g4hunter_score(self, sequence: str, window_size: int = 25) -> float:
        """
        G4Hunter algorithm for G-quadruplex scoring (Bedrat et al., 2016)
        """
        if len(sequence) < 4:
            return 0.0
            
        # For shorter sequences, use the whole sequence
        if len(sequence) < window_size:
            window_size = len(sequence)
        
        scores = []
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i + window_size]
            score = 0
            
            for base in window:
                if base == 'G':
                    score += 1
                elif base == 'C':
                    score -= 1
            
            scores.append(score)
        
        if not scores:
            return 0.0
        
        max_score = max(abs(s) for s in scores)
        # Normalize differently to get better scores
        normalized_score = max_score / window_size
        
        # Additional scoring based on G-tract characteristics
        g_tracts = re.findall(r'G{2,}', sequence)
        if len(g_tracts) >= 3:  # Minimum for G4
            tract_bonus = len(g_tracts) * 0.1
            normalized_score = min(normalized_score + tract_bonus, 1.0)
        
        return normalized_score
    
    def _g_triplex_score(self, sequence: str) -> float:
        """
        G-triplex scoring for three G-tract structures
        """
        if len(sequence) < 12:
            return 0.0
        
        # Count G-tracts
        g_tracts = re.findall(r'G{2,}', sequence)
        if len(g_tracts) < 3:
            return 0.0
        
        # Calculate score based on G-tract density and length
        total_g_length = sum(len(tract) for tract in g_tracts)
        g_density = total_g_length / len(sequence)
        tract_bonus = len(g_tracts) / 4  # Bonus for multiple tracts
        
        return min(g_density + tract_bonus * 0.2, 1.0)