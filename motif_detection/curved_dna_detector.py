"""
Curved DNA Motif Detector
========================

Detects A-tract mediated DNA bending patterns including:
- A-tracts and T-tracts (local curvature)
- Phased A-tracts (global curvature)

Based on Olson et al., 1998 and Crothers 1992.
"""

import re
from typing import List, Dict, Any, Tuple
from .base_detector import BaseMotifDetector


class CurvedDNADetector(BaseMotifDetector):
    """Detector for curved DNA motifs"""
    
    def get_motif_class_name(self) -> str:
        return "Curved_DNA"
    
    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return curved DNA patterns"""
        return {
            'a_tracts': [
                # (pattern, pattern_id, name, subclass, min_len, scoring_method, confidence, biological_significance, reference)
                (r'A{4,15}', 'CRV_1_1', 'A-tract', 'Local Curvature', 4, 'curvature_score', 0.95, 'DNA bending', 'Crothers 1992'),
                (r'T{4,15}', 'CRV_1_2', 'T-tract', 'Local Curvature', 4, 'curvature_score', 0.95, 'DNA bending', 'Crothers 1992'),
                (r'(?:A{3,8}T{1,5}){2,}', 'CRV_1_3', 'AT-rich tract', 'Local Curvature', 8, 'curvature_score', 0.85, 'Sequence-dependent bending', 'Hagerman 1986'),
            ],
            'phased_a_tracts': [
                (r'(?:A{3,8}.{8,12}){3,}A{3,8}', 'CRV_1_4', 'Phased A-tracts', 'Global Curvature', 20, 'phasing_score', 0.90, 'Macroscopic curvature', 'Koo 1986'),
                (r'(?:T{3,8}.{8,12}){3,}T{3,8}', 'CRV_1_5', 'Phased T-tracts', 'Global Curvature', 20, 'phasing_score', 0.90, 'Macroscopic curvature', 'Koo 1986'),
            ]
        }
    
    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        """Calculate curvature score for the sequence"""
        scoring_method = pattern_info[5] if len(pattern_info) > 5 else 'curvature_score'
        
        if scoring_method == 'curvature_score':
            return self._curvature_score(sequence)
        elif scoring_method == 'phasing_score':
            return self._phasing_score(sequence)
        else:
            return self._curvature_score(sequence)
    
    def _curvature_score(self, sequence: str) -> float:
        """
        DNA curvature scoring based on A-tract analysis (Olson et al., 1998)
        """
        if len(sequence) < 4:
            return 0.0
        
        # Find A/T tracts
        a_tracts = re.findall(r'A{3,}', sequence)
        t_tracts = re.findall(r'T{3,}', sequence)
        
        # Calculate curvature based on tract length and frequency
        curvature = 0
        for tract in a_tracts + t_tracts:
            curvature += len(tract) ** 1.5  # Non-linear length dependency
        
        return min(curvature / len(sequence), 1.0)
    
    def _phasing_score(self, sequence: str) -> float:
        """
        Phasing score for periodic A-tract arrangements
        """
        if len(sequence) < 20:
            return 0.0
        
        # Look for periodic A/T tracts with ~10bp spacing
        a_positions = []
        t_positions = []
        
        # Find A-tract positions
        for match in re.finditer(r'A{3,8}', sequence):
            a_positions.append(match.start())
        
        # Find T-tract positions  
        for match in re.finditer(r'T{3,8}', sequence):
            t_positions.append(match.start())
        
        all_positions = sorted(a_positions + t_positions)
        
        if len(all_positions) < 3:
            return 0.0
        
        # Check for periodic spacing (8-12 bp intervals)
        periodic_score = 0
        for i in range(1, len(all_positions)):
            spacing = all_positions[i] - all_positions[i-1]
            if 8 <= spacing <= 12:
                periodic_score += 1
        
        return min(periodic_score / (len(all_positions) - 1), 1.0)