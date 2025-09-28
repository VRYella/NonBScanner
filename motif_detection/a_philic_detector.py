"""
A-philic DNA Motif Detector
===========================

Detects A-rich structures including:
- A-philic tracts

Based on Gorin 1995 and A-philic propensity research.
"""

import re
from typing import List, Dict, Any, Tuple
from .base_detector import BaseMotifDetector


class APhilicDetector(BaseMotifDetector):
    """Detector for A-philic DNA motifs"""
    
    def get_motif_class_name(self) -> str:
        return "A-philic_DNA"
    
    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return A-philic DNA patterns"""
        return {
            'a_philic_tracts': [
                (r'A{6,}', 'APH_9_1', 'Pure A-tract', 'A-philic DNA', 6, 'a_philic_score', 0.95, 'A-philic tract', 'Gorin 1995'),
                (r'(?:A{3,}[TA]{1,3}){2,}', 'APH_9_2', 'AT-rich A-philic', 'A-philic DNA', 8, 'a_philic_score', 0.85, 'A-rich with T insertions', 'Crothers 1990'),
                (r'A{4,}[ATGC]{1,5}A{4,}', 'APH_9_3', 'Spaced A-tracts', 'A-philic DNA', 10, 'a_philic_score', 0.80, 'Spaced A-rich regions', 'Travers 1989'),
                (r'(?:AAA[TA]){3,}', 'APH_9_4', 'AAA-repeat motif', 'A-philic DNA', 12, 'a_philic_score', 0.90, 'Regular A-repeat pattern', 'Koo 1990'),
            ]
        }
    
    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        """Calculate A-philic formation score for the sequence"""
        scoring_method = pattern_info[5] if len(pattern_info) > 5 else 'a_philic_score'
        
        if scoring_method == 'a_philic_score':
            return self._a_philic_score(sequence)
        else:
            return self._a_philic_score(sequence)
    
    def _a_philic_score(self, sequence: str) -> float:
        """
        A-philic DNA scoring based on A-tract characteristics
        """
        if len(sequence) < 6:
            return 0.0
        
        a_content = len(re.findall(r'A', sequence)) / len(sequence)
        a_tracts = re.findall(r'A{3,}', sequence)
        
        # Score based on A content and tract formation
        tract_score = sum(len(tract) ** 1.2 for tract in a_tracts) / len(sequence) if a_tracts else 0
        
        return min(a_content * 0.4 + tract_score * 0.6, 1.0)