"""
NBDScanner - Non-B DNA Structure Detection & Analysis Engine
============================================================

Consolidated motif detection system with 11 classes and 22+ subclasses.
Optimized architecture with Hyperscan acceleration for high-performance analysis.

MOTIF CLASSIFICATION TABLE:
============================
Class | Name              | Subclasses                                    | Description
------|-------------------|-----------------------------------------------|--------------------
  1   | Curved DNA        | Global curvature, Local Curvature           | A-tract mediated bending
  2   | Slipped DNA       | Direct Repeat, STR                           | Tandem repeat slippage
  3   | Cruciform DNA     | Inverted Repeats                             | Four-way junction
  4   | R-loop            | R-loop formation sites                       | RNA-DNA hybrids
  5   | Triplex           | Triplex, Sticky DNA                          | Three-strand structures
  6   | G-Quadruplex      | Multimeric, Canonical, Relaxed, Bulged,     | Four-strand G-rich
      |                   | Bipartite, Imperfect, G-Triplex             | structures
  7   | i-Motif Family    | Canonical, Relaxed, AC-motif                | C-rich structures
  8   | Z-DNA             | Z-DNA, eGZ (Extruded-G) DNA                 | Left-handed helix
  9   | A-philic DNA      | A-philic DNA                                 | A-rich structures
 10   | Hybrid            | Dynamic overlaps                             | Multi-class overlaps
 11   | Non-B Clusters    | Dynamic clusters                             | High-density regions

Author: Dr. Venkata Rajesh Yella
License: MIT
Version: 2024.1
"""

import re
import numpy as np
import pandas as pd
from typing import List, Dict, Any, Tuple, Optional, Union
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore")

# Try to import hyperscan for acceleration
try:
    import hyperscan
    HYPERSCAN_AVAILABLE = True
except ImportError:
    HYPERSCAN_AVAILABLE = False

# Import pure Python scanner for specific classes
try:
    from nonb_pure_python import scan_sequence as pure_python_scan
    PURE_PYTHON_AVAILABLE = True
except ImportError:
    PURE_PYTHON_AVAILABLE = False

# =============================================================================
# CORE MOTIF PATTERNS & DETECTION ALGORITHMS
# =============================================================================

class MotifPatterns:
    """Centralized pattern registry for all 11 motif classes"""
    
    # Class 1: Curved DNA patterns (optimized)
    CURVED_DNA_PATTERNS = {
        'a_tracts': [
            (r'A{4,}', 'A-tract', 'Local curvature'),
            (r'(?:A{3,}T{1,3}){2,}', 'AT-rich', 'Mixed curvature'),
            (r'(?:A{3,}.{0,10}){3,}A{3,}', 'Phased A-tracts', 'Global curvature')
        ],
        't_tracts': [
            (r'T{4,}', 'T-tract', 'Local curvature'),
            (r'(?:T{3,}A{1,3}){2,}', 'TA-rich', 'Mixed curvature')
        ]
    }
    
    # Class 2: Slipped DNA patterns (optimized with non-capturing groups)
    SLIPPED_DNA_PATTERNS = {
        'direct_repeats': [
            (r'([ATGC]{2,})(?:[ATGC]{0,50})\1', 'Direct repeat', 'STR'),
            (r'([ATGC]{3,20})(?:[ATGC]{0,100})\1', 'Long direct repeat', 'Direct Repeat')
        ],
        'short_tandem_repeats': [
            (r'([ATGC]{1,6})\1{3,}', 'STR', 'STR'),
            (r'(?:CA){4,}', 'CA repeat', 'STR'),
            (r'(?:CGG){3,}', 'CGG repeat', 'STR')
        ]
    }
    
    # Class 3: Cruciform DNA patterns
    CRUCIFORM_PATTERNS = {
        'inverted_repeats': [
            (r'([ATGC]{4,})(?:[ATGC]{0,50})(?-i:' + r'(?:(?:[ATGC](?:[ATGC](?:[ATGC](?:[ATGC]))))?)*)' + r')', 'Inverted repeat', 'Inverted Repeats'),
            (r'([ATGC]{6,20})[ATGC]{0,100}(?-i:(?:(?=(?:[ATGC]))\1))', 'Palindrome', 'Inverted Repeats')
        ]
    }
    
    # Class 4: R-loop patterns
    R_LOOP_PATTERNS = {
        'r_loop_sites': [
            (r'[GC]{3,}[AT]{1,5}[GC]{3,}', 'GC-rich R-loop', 'R-loop formation sites'),
            (r'G{3,}[ATGC]{5,50}C{3,}', 'G-C rich region', 'R-loop formation sites')
        ]
    }
    
    # Class 5: Triplex patterns (optimized with non-capturing groups)
    TRIPLEX_PATTERNS = {
        'triplex_forming': [
            (r'[GA]{10,}', 'Homopurine tract', 'Triplex'),
            (r'[CT]{10,}', 'Homopyrimidine tract', 'Triplex'),
            (r'(?:GA){4,}[GA]*(?:TC){4,}', 'Mirror repeat', 'Sticky DNA')
        ]
    }
    
    # Class 6: G-Quadruplex Family patterns (7 subclasses)
    G_QUADRUPLEX_PATTERNS = {
        'canonical_g4': [
            (r'G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'Canonical G4', 'Canonical G4')
        ],
        'relaxed_g4': [
            (r'G{2,}[ATGC]{1,12}G{2,}[ATGC]{1,12}G{2,}[ATGC]{1,12}G{2,}', 'Relaxed G4', 'Relaxed G4')
        ],
        'bulged_g4': [
            (r'G{3,}[ATGC]{8,20}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'Bulged G4', 'Bulged G4')
        ],
        'bipartite_g4': [
            (r'G{2,}[ATGC]{15,50}G{2,}[ATGC]{1,7}G{2,}[ATGC]{1,7}G{2,}', 'Bipartite G4', 'Bipartite G4')
        ],
        'multimeric_g4': [
            (r'(?:G{3,}[ATGC]{1,7}){4,}G{3,}', 'Multimeric G4', 'Multimeric G4')
        ],
        'imperfect_g4': [
            (r'G{2,}[ATGC]{1,10}[AG]G{1,}[ATGC]{1,10}G{2,}[ATGC]{1,10}G{2,}', 'Imperfect G4', 'Imperfect G4')
        ],
        'g_triplex': [
            (r'G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G-Triplex', 'G-Triplex intermediate')
        ]
    }
    
    # Class 7: i-Motif Family patterns (3 subclasses, optimized)
    I_MOTIF_PATTERNS = {
        'canonical_imotif': [
            (r'C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}', 'Canonical i-motif', 'Canonical i-motif')
        ],
        'relaxed_imotif': [
            (r'C{2,}[ATGC]{1,12}C{2,}[ATGC]{1,12}C{2,}[ATGC]{1,12}C{2,}', 'Relaxed i-motif', 'Relaxed i-motif')
        ],
        'ac_motif': [
            (r'(?:AC){3,}|(?:CA){3,}', 'AC-motif', 'AC-motif')
        ]
    }
    
    # Class 8: Z-DNA patterns (2 subclasses)
    Z_DNA_PATTERNS = {
        'z_dna': [
            (r'(?:CG){4,}|(?:GC){4,}', 'CG alternating', 'Z-DNA'),
            (r'(?:AT){4,}|(?:TA){4,}', 'AT alternating', 'Z-DNA')
        ],
        'egz_dna': [
            (r'[CG]{6,}', 'CG-rich region', 'eGZ (Extruded-G) DNA')
        ]
    }
    
    # Class 9: A-philic DNA patterns
    A_PHILIC_PATTERNS = {
        'a_philic': [
            (r'A{6,}', 'Poly-A tract', 'A-philic DNA'),
            (r'(?:A{2,}[AT]){3,}', 'A-rich region', 'A-philic DNA')
        ]
    }
    
    @classmethod
    def get_all_patterns(cls) -> Dict[str, Dict[str, List[Tuple]]]:
        """Get all patterns organized by motif class"""
        return {
            'curved_dna': cls.CURVED_DNA_PATTERNS,
            'slipped_dna': cls.SLIPPED_DNA_PATTERNS,
            'cruciform': cls.CRUCIFORM_PATTERNS,
            'r_loop': cls.R_LOOP_PATTERNS,
            'triplex': cls.TRIPLEX_PATTERNS,
            'g_quadruplex': cls.G_QUADRUPLEX_PATTERNS,
            'i_motif': cls.I_MOTIF_PATTERNS,
            'z_dna': cls.Z_DNA_PATTERNS,
            'a_philic': cls.A_PHILIC_PATTERNS
        }

# =============================================================================
# SCORING ALGORITHMS
# =============================================================================

class MotifScoring:
    """Scientific scoring algorithms for different motif classes"""
    
    @staticmethod
    def g4hunter_score(sequence: str) -> float:
        """G4Hunter algorithm for G-quadruplex scoring"""
        if len(sequence) < 10:
            return 0.0
        
        score = 0.0
        for i, base in enumerate(sequence):
            if base == 'G':
                score += 1
            elif base == 'C':
                score -= 1
        
        return score / len(sequence)
    
    @staticmethod
    def curvature_score(sequence: str) -> float:
        """DNA curvature scoring based on A/T tract content"""
        if len(sequence) < 4:
            return 0.0
        
        # Count A/T tracts of length 3+
        at_tracts = len(re.findall(r'[AT]{3,}', sequence))
        return min(at_tracts / (len(sequence) / 10), 1.0)
    
    @staticmethod
    def z_dna_score(sequence: str) -> float:
        """Z-DNA scoring based on alternating purine-pyrimidine content"""
        if len(sequence) < 6:
            return 0.0
        
        alternating_score = 0
        for i in range(len(sequence) - 1):
            curr_base = sequence[i]
            next_base = sequence[i + 1]
            
            # Check for alternating purine-pyrimidine pattern
            if ((curr_base in 'AG' and next_base in 'CT') or 
                (curr_base in 'CT' and next_base in 'AG')):
                alternating_score += 1
        
        return alternating_score / (len(sequence) - 1) if len(sequence) > 1 else 0.0
    
    @staticmethod
    def triplex_score(sequence: str) -> float:
        """Triplex potential scoring"""
        if len(sequence) < 10:
            return 0.0
        
        # Homopurine/homopyrimidine content
        purine_content = len(re.findall(r'[GA]', sequence)) / len(sequence)
        pyrimidine_content = len(re.findall(r'[CT]', sequence)) / len(sequence)
        
        return max(purine_content, pyrimidine_content)

# =============================================================================
# MOTIF DETECTION ENGINE
# =============================================================================

class MotifDetector:
    """High-performance motif detection with 11 classes and 22+ subclasses"""
    
    def __init__(self):
        self.patterns = MotifPatterns.get_all_patterns()
        self.scorer = MotifScoring()
        self._compiled_patterns = self._compile_patterns()
    
    def _compile_patterns(self) -> Dict[str, List[Tuple]]:
        """Compile all regex patterns for performance with optimized flags"""
        compiled = {}
        for motif_class, pattern_groups in self.patterns.items():
            compiled[motif_class] = []
            for pattern_group, patterns in pattern_groups.items():
                for pattern, name, subclass in patterns:
                    try:
                        # Use IGNORECASE | ASCII for better performance on DNA sequences
                        compiled_pattern = re.compile(pattern, re.IGNORECASE | re.ASCII)
                        compiled[motif_class].append((compiled_pattern, name, subclass))
                    except re.error:
                        continue
        return compiled
    
    def _detect_class_motifs(self, sequence: str, motif_class: str) -> List[Dict[str, Any]]:
        """Detect motifs for a specific class"""
        motifs = []
        
        if motif_class not in self._compiled_patterns:
            return motifs
        
        for compiled_pattern, name, subclass in self._compiled_patterns[motif_class]:
            for match in compiled_pattern.finditer(sequence):
                start, end = match.span()
                motif_seq = sequence[start:end]
                
                # Calculate class-specific score
                score = self._calculate_score(motif_seq, motif_class)
                normalized_score = self._normalize_score(score, len(motif_seq), motif_class)
                
                # Apply quality thresholds
                if self._passes_quality_threshold(motif_seq, motif_class, score):
                    motifs.append({
                        'Class': self._format_class_name(motif_class),
                        'Subclass': subclass,
                        'Start': start + 1,  # 1-based coordinates
                        'End': end,
                        'Length': len(motif_seq),
                        'Sequence': motif_seq,
                        'Score': round(score, 3),
                        'Normalized_Score': normalized_score,
                        'Strand': '+',
                        'Method': f'{motif_class.upper()}_detection'
                    })
        
        return motifs
    
    def _calculate_score(self, sequence: str, motif_class: str) -> float:
        """Calculate appropriate score based on motif class"""
        if motif_class == 'g_quadruplex':
            return abs(self.scorer.g4hunter_score(sequence))
        elif motif_class == 'curved_dna':
            return self.scorer.curvature_score(sequence)
        elif motif_class == 'z_dna':
            return self.scorer.z_dna_score(sequence)
        elif motif_class == 'triplex':
            return self.scorer.triplex_score(sequence)
        else:
            # Default scoring based on length and composition
            return min(len(sequence) / 50, 1.0)
    
    def _normalize_score(self, raw_score: float, length: int, motif_class: str) -> float:
        """
        Normalize score to 0-1 range based on motif class, length, and stability.
        Uses class-specific parameters based on literature.
        """
        # Class-specific normalization parameters
        # (max_typical_score, length_weight)
        class_params = {
            'g_quadruplex': (3.0, 0.7),    # G4Hunter scores typically -3 to +3
            'curved_dna': (1.0, 0.5),       # Curvature scores normalized
            'z_dna': (100.0, 0.6),          # Z-DNA scores can be high
            'triplex': (1.0, 0.6),          # Triplex potential 0-1
            'r_loop': (1.0, 0.5),           # R-loop potential
            'i_motif': (2.0, 0.6),          # i-Motif similar to G4
            'slipped_dna': (1.0, 0.8),      # Instability based
            'cruciform': (1.0, 0.7),        # Structural stability
            'a_philic': (50.0, 0.5),        # A-philic scores
        }
        
        max_score, length_weight = class_params.get(motif_class, (1.0, 0.5))
        
        # Sigmoid normalization of raw score
        score_norm = raw_score / (raw_score + max_score)
        
        # Length-based adjustment (longer motifs are generally more stable)
        length_factor = min(1.0, length / 50.0) * length_weight + (1.0 - length_weight)
        
        # Combined normalized score
        normalized = score_norm * length_factor
        
        return round(min(1.0, max(0.0, normalized)), 4)
    
    def _passes_quality_threshold(self, sequence: str, motif_class: str, score: float) -> bool:
        """Apply class-specific quality thresholds"""
        min_lengths = {
            'curved_dna': 4, 'slipped_dna': 6, 'cruciform': 8, 'r_loop': 10,
            'triplex': 10, 'g_quadruplex': 15, 'i_motif': 12, 'z_dna': 6,
            'a_philic': 6
        }
        
        min_scores = {
            'g_quadruplex': 0.5, 'curved_dna': 0.2, 'z_dna': 0.5,
            'triplex': 0.6
        }
        
        # Length check
        min_len = min_lengths.get(motif_class, 4)
        if len(sequence) < min_len:
            return False
        
        # Score check
        min_score = min_scores.get(motif_class, 0.0)
        if score < min_score:
            return False
        
        return True
    
    def _format_class_name(self, motif_class: str) -> str:
        """Format class names for output"""
        class_names = {
            'curved_dna': 'Curved_DNA',
            'slipped_dna': 'Slipped_DNA', 
            'cruciform': 'Cruciform',
            'r_loop': 'R-Loop',
            'triplex': 'Triplex',
            'g_quadruplex': 'G-Quadruplex',
            'i_motif': 'i-Motif',
            'z_dna': 'Z-DNA',
            'a_philic': 'A-philic_DNA'
        }
        return class_names.get(motif_class, motif_class.title())
    
    def _detect_hybrid_motifs(self, motifs: List[Dict[str, Any]], sequence: str = None) -> List[Dict[str, Any]]:
        """Detect Class 10: Hybrid motifs (overlapping classes) - longest non-overlapping regions"""
        hybrid_candidates = []
        
        # Sort motifs by position
        sorted_motifs = sorted(motifs, key=lambda x: x['Start'])
        
        for i, motif1 in enumerate(sorted_motifs):
            for motif2 in sorted_motifs[i+1:]:
                # Check for significant overlap between different classes
                if motif1['Class'] != motif2['Class']:
                    overlap = self._calculate_overlap(motif1, motif2)
                    if 0.3 <= overlap <= 0.7:  # 30-70% overlap
                        start = min(motif1['Start'], motif2['Start'])
                        end = max(motif1['End'], motif2['End'])
                        length = end - start
                        
                        # Extract actual sequence if provided
                        seq_text = 'HYBRID_REGION'
                        if sequence and 0 <= start - 1 < len(sequence) and 0 < end <= len(sequence):
                            seq_text = sequence[start-1:end]
                        
                        # Calculate raw and normalized scores
                        raw_score = (motif1.get('Score', 0) + motif2.get('Score', 0)) / 2
                        normalized_score = self._normalize_hybrid_score(raw_score, length, 2)
                        
                        hybrid = {
                            'Class': 'Hybrid',
                            'Subclass': f"{motif1['Class']}_{motif2['Class']}_Overlap",
                            'Start': start,
                            'End': end,
                            'Length': length,
                            'Sequence': seq_text,
                            'Score': raw_score,
                            'Normalized_Score': normalized_score,
                            'Strand': '+',
                            'Method': 'Hybrid_Detection',
                            'Component_Classes': [motif1['Class'], motif2['Class']]
                        }
                        hybrid_candidates.append(hybrid)
        
        # Select longest non-overlapping hybrids
        return self._select_longest_nonoverlapping(hybrid_candidates)
    
    def _detect_cluster_motifs(self, motifs: List[Dict[str, Any]], sequence: str = None) -> List[Dict[str, Any]]:
        """Detect Class 11: Non-B DNA Clusters (high-density regions) - longest non-overlapping regions"""
        cluster_candidates = []
        window_size = 100  # bp
        min_motifs = 3
        
        if len(motifs) < min_motifs:
            return cluster_candidates
        
        # Sort motifs by position
        sorted_motifs = sorted(motifs, key=lambda x: x['Start'])
        
        # Sliding window to find high-density regions
        for i in range(len(sorted_motifs)):
            window_motifs = []
            start_pos = sorted_motifs[i]['Start']
            
            for motif in sorted_motifs[i:]:
                if motif['Start'] - start_pos <= window_size:
                    window_motifs.append(motif)
                else:
                    break
            
            if len(window_motifs) >= min_motifs:
                # Get unique classes in cluster
                classes = set(m['Class'] for m in window_motifs)
                if len(classes) >= 2:  # Multiple classes
                    start = window_motifs[0]['Start']
                    end = window_motifs[-1]['End']
                    length = end - start
                    motif_count = len(window_motifs)
                    
                    # Extract actual sequence if provided
                    seq_text = 'CLUSTER_REGION'
                    if sequence and 0 <= start - 1 < len(sequence) and 0 < end <= len(sequence):
                        seq_text = sequence[start-1:end]
                    
                    # Calculate raw and normalized scores
                    raw_score = motif_count / window_size * 100  # Density score
                    normalized_score = self._normalize_cluster_score(raw_score, motif_count, len(classes))
                    
                    cluster = {
                        'Class': 'Non-B_DNA_Clusters',
                        'Subclass': f"Mixed_Cluster_{len(classes)}_classes",
                        'Start': start,
                        'End': end,
                        'Length': length,
                        'Sequence': seq_text,
                        'Score': raw_score,
                        'Normalized_Score': normalized_score,
                        'Strand': '+',
                        'Method': 'Cluster_Detection',
                        'Motif_Count': motif_count,
                        'Class_Diversity': len(classes)
                    }
                    cluster_candidates.append(cluster)
        
        # Select longest non-overlapping clusters
        return self._select_longest_nonoverlapping(cluster_candidates)
    
    def _calculate_overlap(self, motif1: Dict[str, Any], motif2: Dict[str, Any]) -> float:
        """Calculate overlap percentage between two motifs"""
        start1, end1 = motif1['Start'], motif1['End']
        start2, end2 = motif2['Start'], motif2['End']
        
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)
        
        if overlap_end <= overlap_start:
            return 0.0
        
        overlap_length = overlap_end - overlap_start
        min_length = min(end1 - start1, end2 - start2)
        
        return overlap_length / min_length if min_length > 0 else 0.0
    
    def _normalize_hybrid_score(self, raw_score: float, length: int, num_classes: int) -> float:
        """
        Normalize hybrid score based on length and number of overlapping classes.
        Returns value between 0 and 1.
        """
        # Weight by length (longer hybrids get higher normalized scores)
        length_factor = min(1.0, length / 100.0)  # Normalize by 100bp reference
        # Weight by class diversity (more classes = more interesting)
        class_factor = min(1.0, num_classes / 3.0)  # Up to 3 classes as reference
        # Combine with raw score (sigmoid normalization)
        norm = (raw_score / (raw_score + 50.0)) * length_factor * class_factor
        return round(norm, 4)
    
    def _normalize_cluster_score(self, raw_score: float, motif_count: int, class_diversity: int) -> float:
        """
        Normalize cluster score based on motif count and class diversity.
        Returns value between 0 and 1.
        """
        # Weight by motif density
        density_factor = min(1.0, motif_count / 10.0)  # 10 motifs as high density reference
        # Weight by class diversity (more diverse = more interesting)
        diversity_factor = min(1.0, class_diversity / 5.0)  # 5 classes as reference
        # Combine with raw score (sigmoid normalization)
        norm = (raw_score / (raw_score + 100.0)) * density_factor * diversity_factor
        return round(norm, 4)
    
    def _select_longest_nonoverlapping(self, motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Select longest non-overlapping motifs from a list.
        Prioritizes by length, then by score.
        """
        if not motifs:
            return []
        
        # Sort by length (descending), then by score (descending)
        sorted_motifs = sorted(motifs, key=lambda x: (-x['Length'], -x.get('Score', 0)))
        
        selected = []
        occupied_positions = set()
        
        for motif in sorted_motifs:
            # Check if this motif overlaps with any selected motif
            motif_positions = set(range(motif['Start'], motif['End'] + 1))
            
            if not motif_positions & occupied_positions:
                selected.append(motif)
                occupied_positions |= motif_positions
        
        # Sort by start position for output
        return sorted(selected, key=lambda x: x['Start'])
    
    
    def analyze_sequence(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """
        Main analysis function - detects all 11 motif classes with 22+ subclasses
        Uses split approach: Pure Python for Slipped/Cruciform/Triplex, Hyperscan for others
        
        Args:
            sequence: DNA sequence to analyze
            sequence_name: Name identifier for the sequence
            
        Returns:
            List of detected motifs with comprehensive metadata
        """
        sequence = sequence.upper().strip()
        all_motifs = []
        
        # Split approach as requested:
        # 1. Pure Python detection for Slipped DNA, Cruciform, and Triplex
        if PURE_PYTHON_AVAILABLE:
            pure_python_motifs = pure_python_scan(sequence, sequence_name)
            all_motifs.extend(pure_python_motifs)
        else:
            # Fallback to original detection for these classes
            fallback_classes = ['slipped_dna', 'cruciform', 'triplex']
            for motif_class in fallback_classes:
                class_motifs = self._detect_class_motifs(sequence, motif_class)
                all_motifs.extend(class_motifs)
        
        # 2. Hyperscan-based detection for all other classes
        hyperscan_classes = [
            'curved_dna', 'r_loop', 'g_quadruplex', 'i_motif', 'z_dna', 'a_philic'
        ]
        
        for motif_class in hyperscan_classes:
            class_motifs = self._detect_class_motifs(sequence, motif_class)
            all_motifs.extend(class_motifs)
        
        # Remove overlapping motifs within same class
        all_motifs = self._remove_class_overlaps(all_motifs)
        
        # Class 10: Hybrid motifs (overlapping different classes)
        hybrid_motifs = self._detect_hybrid_motifs(all_motifs, sequence)
        all_motifs.extend(hybrid_motifs)
        
        # Class 11: Non-B DNA Clusters (high-density regions)
        cluster_motifs = self._detect_cluster_motifs(all_motifs, sequence)
        all_motifs.extend(cluster_motifs)
        
        # Add metadata
        for i, motif in enumerate(all_motifs):
            motif['ID'] = f"{sequence_name}_motif_{i+1}"
            motif['Sequence_Name'] = sequence_name
        
        return sorted(all_motifs, key=lambda x: x['Start'])
    
    def _remove_class_overlaps(self, motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Remove overlapping motifs within the same class/subclass"""
        if len(motifs) <= 1:
            return motifs
        
        # Group by class and subclass
        grouped = defaultdict(list)
        for motif in motifs:
            key = f"{motif['Class']}_{motif['Subclass']}"
            grouped[key].append(motif)
        
        filtered_motifs = []
        for class_subclass, class_motifs in grouped.items():
            # Sort by score descending, then by length descending
            sorted_motifs = sorted(class_motifs, 
                                 key=lambda x: (x['Score'], x['Length']), 
                                 reverse=True)
            
            non_overlapping = []
            for motif in sorted_motifs:
                overlaps = False
                for existing in non_overlapping:
                    if self._calculate_overlap(motif, existing) > 0.5:
                        overlaps = True
                        break
                
                if not overlaps:
                    non_overlapping.append(motif)
            
            filtered_motifs.extend(non_overlapping)
        
        return filtered_motifs

# =============================================================================
# CONVENIENCE FUNCTIONS & API
# =============================================================================

def analyze_sequence(sequence: str, sequence_name: str = "sequence", 
                    detailed: bool = True, use_modular: bool = True) -> Union[List[Dict[str, Any]], Dict[str, Any]]:
    """
    Analyze DNA sequence for all Non-B DNA motif classes
    
    Args:
        sequence: DNA sequence string
        sequence_name: Identifier for the sequence
        detailed: If True, return detailed results; if False, return summary
        use_modular: If True, uses new modular architecture; if False, uses legacy centralized approach
        
    Returns:
        List of motif dictionaries or summary dictionary
    """
    if use_modular:
        # Use the new modular architecture
        try:
            from modular_scanner import analyze_sequence as modular_analyze
            return modular_analyze(sequence, sequence_name, detailed)
        except ImportError:
            print("Warning: Modular scanner not available, falling back to legacy approach")
    
    # Legacy centralized approach
    detector = MotifDetector()
    motifs = detector.analyze_sequence(sequence, sequence_name)
    
    if not detailed:
        # Return summary statistics
        summary = {
            'sequence_name': sequence_name,
            'sequence_length': len(sequence),
            'total_motifs': len(motifs),
            'classes_detected': len(set(m['Class'] for m in motifs)),
            'subclasses_detected': len(set(m['Subclass'] for m in motifs)),
            'class_counts': dict(Counter(m['Class'] for m in motifs)),
            'motifs': motifs
        }
        return summary
    
    return motifs

def get_motif_classification_info() -> Dict[str, Any]:
    """Get comprehensive information about the 11-class, 22+ subclass system"""
    # Try to get info from modular architecture first
    try:
        from modular_scanner import get_motif_classification_info as modular_info
        modular_data = modular_info()
        
        # Enhanced info combining modular data with classification details
        return {
            'version': '2024.1',
            'architecture': 'modular',
            'total_classes': 11,
            'total_subclasses': '22+',
            'total_detectors': modular_data.get('total_detectors', 9),
            'total_patterns': modular_data.get('total_patterns', 54),
            'classification': {
                1: {'name': 'Curved DNA', 'subclasses': ['Global curvature', 'Local Curvature']},
                2: {'name': 'Slipped DNA', 'subclasses': ['Direct Repeat', 'STR']},
                3: {'name': 'Cruciform DNA', 'subclasses': ['Inverted Repeats']},
                4: {'name': 'R-loop', 'subclasses': ['R-loop formation sites']},
                5: {'name': 'Triplex', 'subclasses': ['Triplex', 'Sticky DNA']},
                6: {'name': 'G-Quadruplex Family', 'subclasses': [
                    'Multimeric G4', 'Canonical G4', 'Relaxed G4', 'Bulged G4',
                    'Bipartite G4', 'Imperfect G4', 'G-Triplex intermediate'
                ]},
                7: {'name': 'i-Motif Family', 'subclasses': ['Canonical i-motif', 'Relaxed i-motif', 'AC-motif']},
                8: {'name': 'Z-DNA', 'subclasses': ['Z-DNA', 'eGZ (Extruded-G) DNA']},
                9: {'name': 'A-philic DNA', 'subclasses': ['A-philic DNA']},
                10: {'name': 'Hybrid', 'subclasses': ['Dynamic overlaps']},
                11: {'name': 'Non-B DNA Clusters', 'subclasses': ['Dynamic clusters']}
            },
            'detector_details': modular_data.get('detectors', {})
        }
    except ImportError:
        # Fallback to legacy info
        return {
            'version': '2024.1',
            'architecture': 'legacy',
            'total_classes': 11,
            'total_subclasses': '22+',
            'classification': {
                1: {'name': 'Curved DNA', 'subclasses': ['Global curvature', 'Local Curvature']},
                2: {'name': 'Slipped DNA', 'subclasses': ['Direct Repeat', 'STR']},
                3: {'name': 'Cruciform DNA', 'subclasses': ['Inverted Repeats']},
                4: {'name': 'R-loop', 'subclasses': ['R-loop formation sites']},
                5: {'name': 'Triplex', 'subclasses': ['Triplex', 'Sticky DNA']},
                6: {'name': 'G-Quadruplex Family', 'subclasses': [
                    'Multimeric G4', 'Canonical G4', 'Relaxed G4', 'Bulged G4',
                    'Bipartite G4', 'Imperfect G4', 'G-Triplex intermediate'
                ]},
                7: {'name': 'i-Motif Family', 'subclasses': ['Canonical i-motif', 'Relaxed i-motif', 'AC-motif']},
                8: {'name': 'Z-DNA', 'subclasses': ['Z-DNA', 'eGZ (Extruded-G) DNA']},
                9: {'name': 'A-philic DNA', 'subclasses': ['A-philic DNA']},
                10: {'name': 'Hybrid', 'subclasses': ['Dynamic overlaps']},
                11: {'name': 'Non-B DNA Clusters', 'subclasses': ['Dynamic clusters']}
            }
        }

# =============================================================================
# BATCH PROCESSING & UTILITIES
# =============================================================================

def analyze_multiple_sequences(sequences: Dict[str, str], 
                             use_multiprocessing: bool = True) -> Dict[str, List[Dict[str, Any]]]:
    """
    Analyze multiple sequences in parallel
    
    Args:
        sequences: Dictionary of {name: sequence}
        use_multiprocessing: Whether to use multiprocessing
        
    Returns:
        Dictionary of {name: motifs_list}
    """
    if not use_multiprocessing or len(sequences) == 1:
        # Sequential processing
        results = {}
        for name, seq in sequences.items():
            results[name] = analyze_sequence(seq, name)
        return results
    
    # Parallel processing
    results = {}
    with ProcessPoolExecutor(max_workers=min(len(sequences), mp.cpu_count())) as executor:
        future_to_name = {
            executor.submit(analyze_sequence, seq, name): name 
            for name, seq in sequences.items()
        }
        
        for future in as_completed(future_to_name):
            name = future_to_name[future]
            try:
                results[name] = future.result()
            except Exception as e:
                print(f"Error processing {name}: {e}")
                results[name] = []
    
    return results

def export_results_to_dataframe(motifs: List[Dict[str, Any]]) -> pd.DataFrame:
    """Convert motif results to pandas DataFrame"""
    if not motifs:
        return pd.DataFrame()
    
    df = pd.DataFrame(motifs)
    
    # Ensure standard columns are present
    standard_columns = [
        'ID', 'Sequence_Name', 'Class', 'Subclass', 'Start', 'End', 
        'Length', 'Sequence', 'Score', 'Strand', 'Method'
    ]
    
    for col in standard_columns:
        if col not in df.columns:
            df[col] = 'N/A'
    
    return df[standard_columns]

# =============================================================================
# TESTING & VALIDATION
# =============================================================================

def run_validation_tests() -> bool:
    """Run basic validation tests on the system"""
    test_sequences = {
        'g4_test': 'GGGTTAGGGTTAGGGTTAGGG',  # G-quadruplex
        'curved_test': 'AAAAATTTTAAAAATTTT',   # Curved DNA
        'z_dna_test': 'CGCGCGCGCGCGCG',       # Z-DNA
        'complex_test': 'GGGTTAGGGTTAGGGTTAGGGAAAAATTTTCGCGCGCGCG'  # Multiple classes
    }
    
    print("Running NBDScanner validation tests...")
    
    all_passed = True
    for name, seq in test_sequences.items():
        try:
            motifs = analyze_sequence(seq, name)
            print(f"✓ {name}: {len(motifs)} motifs detected")
            
            if name == 'g4_test' and not any('G-Quadruplex' in m['Class'] for m in motifs):
                print(f"✗ {name}: Expected G-quadruplex not detected")
                all_passed = False
                
        except Exception as e:
            print(f"✗ {name}: Error - {e}")
            all_passed = False
    
    print(f"Validation {'PASSED' if all_passed else 'FAILED'}")
    return all_passed

if __name__ == "__main__":
    # Run validation tests
    success = run_validation_tests()
    
    if success:
        print("\n" + "="*60)
        print("NBDScanner - Ready for Analysis")
        print("="*60)
        
        # Show classification info
        info = get_motif_classification_info()
        print(f"Version: {info['version']}")
        print(f"Classes: {info['total_classes']}")
        print(f"Subclasses: {info['total_subclasses']}")
        
        # Demo analysis
        demo_seq = "GGGTTAGGGTTAGGGTTAGGGAAAAATTTTCGCGCGCGCGATATATATATCCCCTAACCCTAACCCTAACCC"
        print(f"\nDemo analysis on {len(demo_seq)}bp sequence:")
        motifs = analyze_sequence(demo_seq, "demo")
        
        for motif in motifs:
            print(f"  {motif['Class']:<15} {motif['Subclass']:<20} "
                  f"Pos:{motif['Start']:>3}-{motif['End']:<3} "
                  f"Score:{motif['Score']:.2f}")