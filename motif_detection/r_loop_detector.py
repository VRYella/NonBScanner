"""
R-loop DNA Motif Detector
========================

Detects RNA-DNA hybrid structures including:
- R-loop formation sites
- QmRLFS models

Based on Aguilera 2012 and Jenjaroenpun 2016.
"""

import re
from typing import List, Dict, Any, Tuple
from .base_detector import BaseMotifDetector


class RLoopDetector(BaseMotifDetector):
    """Detector for R-loop DNA motifs"""
    
    def get_motif_class_name(self) -> str:
        return "R-Loop"
    
    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return R-loop DNA patterns - optimized with non-capturing groups"""
        return {
            'r_loop_formation_sites': [
                (r'[GC]{10,}[AT]{2,10}[GC]{10,}', 'RLP_4_1', 'GC-rich R-loop site', 'R-loop formation sites', 20, 'r_loop_potential', 0.85, 'Transcription-replication conflicts', 'Aguilera 2012'),
                (r'G{5,}[ATGC]{10,100}C{5,}', 'RLP_4_2', 'G-C rich region', 'R-loop formation sites', 20, 'r_loop_potential', 0.80, 'R-loop prone regions', 'Ginno 2012'),
                (r'[GC]{6,}[AT]{1,5}[GC]{6,}', 'RLP_4_3', 'GC-AT pattern', 'R-loop formation sites', 15, 'r_loop_potential', 0.75, 'Transcriptional pausing', 'Skourti-Stathaki 2011'),
            ],
            'qmrlfs_model_1': [
                (r'G{3,}[ATCGU]{1,10}?G{3,}(?:[ATCGU]{1,10}?G{3,}){1,}?', 'QmRLFS_4_1', 'QmRLFS Model 1', 'QmRLFS-m1', 25, 'qmrlfs_score', 0.90, 'RIZ detection with 3+ G tracts', 'Jenjaroenpun 2016'),
            ],
            'qmrlfs_model_2': [
                (r'G{4,}(?:[ATCGU]{1,10}?G{4,}){1,}?', 'QmRLFS_4_2', 'QmRLFS Model 2', 'QmRLFS-m2', 30, 'qmrlfs_score', 0.95, 'RIZ detection with 4+ G tracts', 'Jenjaroenpun 2016'),
            ]
        }
    
    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        """Calculate R-loop formation score for the sequence"""
        scoring_method = pattern_info[5] if len(pattern_info) > 5 else 'r_loop_potential'
        
        if scoring_method == 'r_loop_potential':
            return self._r_loop_potential(sequence)
        elif scoring_method == 'qmrlfs_score':
            return self._qmrlfs_score(sequence)
        else:
            return self._r_loop_potential(sequence)
    
    def _r_loop_potential(self, sequence: str) -> float:
        """
        R-loop formation potential (Aguilera & García-Muse, 2012)
        """
        if len(sequence) < 20:
            return 0.0
        
        # GC skew calculation
        gc_skew = 0
        for base in sequence:
            if base == 'G':
                gc_skew += 1
            elif base == 'C':
                gc_skew -= 1
        
        # Normalize by sequence length
        gc_skew = abs(gc_skew) / len(sequence)
        
        # GC content (R-loops prefer GC-rich regions)
        gc_content = len(re.findall(r'[GC]', sequence)) / len(sequence)
        
        return min(gc_skew * 0.6 + gc_content * 0.4, 1.0)
    
    def _qmrlfs_score(self, sequence: str) -> float:
        """
        Simplified QmRLFS-based R-loop formation scoring
        """
        if len(sequence) < 25:
            return 0.0
        
        # Count G-tracts
        g_tracts = re.findall(r'G{3,}', sequence)
        if len(g_tracts) < 2:
            return 0.0
        
        # Calculate score based on G-tract density and spacing
        total_g_length = sum(len(tract) for tract in g_tracts)
        g_density = total_g_length / len(sequence)
        tract_bonus = len(g_tracts) / 5  # Bonus for multiple tracts
        
        return min(g_density + tract_bonus * 0.3, 1.0)

    def passes_quality_threshold(self, sequence: str, score: float, pattern_info: Tuple) -> bool:
        """Lower threshold for R-loop detection"""
        return score >= 0.3  # Lower threshold for better sensitivity

    def _remove_overlaps(self, motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Remove overlapping motifs, keeping highest scoring non-overlapping set"""
        if not motifs:
            return []
        
        # Sort by score (descending), then by length (descending)
        sorted_motifs = sorted(motifs, 
                              key=lambda x: (-x.get('Score', 0), -x.get('Length', 0)))
        
        non_overlapping = []
        for motif in sorted_motifs:
            # Check if this motif overlaps with any already selected
            overlaps = False
            for selected in non_overlapping:
                # Two motifs overlap if their regions overlap
                if not (motif['End'] <= selected['Start'] or 
                       motif['Start'] >= selected['End']):
                    overlaps = True
                    break
            
            if not overlaps:
                non_overlapping.append(motif)
        
        # Sort by start position for output
        non_overlapping.sort(key=lambda x: x['Start'])
        return non_overlapping

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Override base method to use custom R-loop detection logic"""
        sequence = sequence.upper().strip()
        motifs = []
        
        # Use base method but with lowered thresholds
        base_motifs = super().detect_motifs(sequence, sequence_name)
        
        # Also add simple GC-rich region detection
        # Look for GC-rich regions that might form R-loops
        gc_pattern = re.compile(r'[GC]{8,}', re.IGNORECASE | re.ASCII)
        for match in gc_pattern.finditer(sequence):
            start, end = match.span()
            motif_seq = sequence[start:end]
            
            # Calculate GC content
            gc_content = len(re.findall(r'[GC]', motif_seq)) / len(motif_seq)
            if gc_content >= 0.75:  # At least 75% GC
                score = gc_content * 0.8  # Simple GC-based score
                
                motifs.append({
                    'ID': f"{sequence_name}_RLP_GC_{start+1}",
                    'Sequence_Name': sequence_name,
                    'Class': self.get_motif_class_name(),
                    'Subclass': 'R-loop formation sites',
                    'Start': start + 1,  # 1-based coordinates
                    'End': end,
                    'Length': len(motif_seq),
                    'Sequence': motif_seq,
                    'Score': round(score, 3),
                    'Strand': '+',
                    'Method': 'R-Loop_detection',
                    'Pattern_ID': f'RLP_GC_{start+1}'
                })
        
        motifs.extend(base_motifs)
        
        # Remove overlapping motifs
        motifs = self._remove_overlaps(motifs)
        
        return motifs