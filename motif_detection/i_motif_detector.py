"""
i-Motif DNA Motif Detector
=========================

Detects C-rich structures including:
- Canonical i-motif
- Relaxed i-motif
- AC-motif

Based on Zeraati 2018 and pH stability research.
"""

import re
from typing import List, Dict, Any, Tuple
from .base_detector import BaseMotifDetector


class IMotifDetector(BaseMotifDetector):
    """Detector for i-motif DNA motifs"""
    
    def get_motif_class_name(self) -> str:
        return "i-Motif"
    
    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return i-motif DNA patterns"""
        return {
            'canonical_imotif': [
                (r'C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}', 'IM_7_1', 'Canonical i-motif', 'Canonical i-motif', 15, 'imotif_score', 0.95, 'pH-dependent C-rich structure', 'Gehring 1993'),
                (r'C{4,}[ATGC]{1,5}C{4,}[ATGC]{1,5}C{4,}[ATGC]{1,5}C{4,}', 'IM_7_2', 'High-density i-motif', 'Canonical i-motif', 16, 'imotif_score', 0.98, 'Stable i-motif', 'Leroy 1993'),
            ],
            'relaxed_imotif': [
                (r'C{2,}[ATGC]{1,12}C{2,}[ATGC]{1,12}C{2,}[ATGC]{1,12}C{2,}', 'IM_7_3', 'Relaxed i-motif', 'Relaxed i-motif', 12, 'imotif_score', 0.80, 'Potential i-motif structures', 'Mergny 1995'),
                (r'C{3,}[ATGC]{8,15}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}', 'IM_7_4', 'Long-loop i-motif', 'Relaxed i-motif', 18, 'imotif_score', 0.75, 'Alternative i-motif topology', 'Phan 2002'),
            ],
            'ac_motif': [
                (r'(?:AC){6,}', 'IM_7_5', 'AC-motif', 'AC-motif', 12, 'imotif_score', 0.85, 'AC repeat i-motif', 'Catasti 1999'),
                (r'A{2,}(?:CA){3,}(?:A{2,})?', 'IM_7_6', 'A-rich AC-motif', 'AC-motif', 10, 'imotif_score', 0.80, 'A-rich AC repeat', 'Kaushik 2006'),
            ]
        }
    
    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        """Calculate i-motif formation score for the sequence"""
        scoring_method = pattern_info[5] if len(pattern_info) > 5 else 'imotif_score'
        
        if scoring_method == 'imotif_score':
            return self._imotif_score(sequence)
        else:
            return self._imotif_score(sequence)
    
    def _imotif_score(self, sequence: str) -> float:
        """
        i-motif scoring based on C-tract analysis (Zeraati et al., 2018)
        """
        if len(sequence) < 12:
            return 0.0
        
        # Count C-tracts
        c_tracts = re.findall(r'C{2,}', sequence)
        if len(c_tracts) < 3:
            return 0.0
        
        # Calculate score based on C-tract density and length
        total_c_length = sum(len(tract) for tract in c_tracts)
        c_density = total_c_length / len(sequence)
        tract_bonus = len(c_tracts) / 4  # Bonus for multiple tracts
        
        return min(c_density + tract_bonus * 0.2, 1.0)