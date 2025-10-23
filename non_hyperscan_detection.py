"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                  NON-HYPERSCAN DETECTION ENGINE MODULE                       ║
║                    Non-B DNA Motif Detection System                          ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: non_hyperscan_detection.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    Algorithmic detection for motif classes not feasible with simple regex.
    Implements specialized algorithms for complex structural patterns.

DETECTED CLASSES (Algorithmic):
┌──────┬─────────────────────┬─────────────────────────────────────────────────┐
│Class │ Name                │ Why Not Hyperscan Compatible                    │
├──────┼─────────────────────┼─────────────────────────────────────────────────┤
│  2*  │ Slipped DNA (DR)    │ Requires backtracking for direct repeats       │
│  3   │ Cruciform           │ Requires reverse complement computation         │
│  9   │ A-philic DNA        │ Requires tetranucleotide log2-odds scoring      │
│ 10   │ Hybrid              │ Post-processing of motif overlaps               │
│ 11   │ Non-B Clusters      │ Sliding window density analysis                 │
└──────┴─────────────────────┴─────────────────────────────────────────────────┘

* Class 2 Direct Repeats only; STRs handled by Hyperscan

ALGORITHMS:
    1. Cruciform: Inverted repeat finder with reverse complement
    2. A-philic DNA: Tetranucleotide scoring (Vinogradov, 2003)
    3. Direct Repeats: Suffix tree or dynamic programming
    4. Hybrid: Interval overlap analysis (30-70% threshold)
    5. Clusters: Sliding window with density threshold

REFERENCES:
    - Vinogradov, 2003 (A-philic tetranucleotide propensities)
    - Lilley, 2000 (Cruciform structures)
    - Sinden, 1994 (DNA structure and function)
"""

import re
import numpy as np
from typing import List, Dict, Any, Tuple, Optional
from collections import defaultdict, Counter


class CruciformDetector:
    """
    Detect cruciform-forming inverted repeats.
    
    Requires reverse complement checking, not feasible with Hyperscan.
    """
    
    @staticmethod
    def reverse_complement(sequence: str) -> str:
        """Get reverse complement of DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement.get(b, 'N') for b in reversed(sequence))
    
    @staticmethod
    def find_inverted_repeats(
        sequence: str,
        min_arm_length: int = 6,
        max_arm_length: int = 50,
        max_spacer: int = 50
    ) -> List[Dict[str, Any]]:
        """
        Find inverted repeats (palindromes) in sequence.
        
        Args:
            sequence: DNA sequence
            min_arm_length: Minimum length of palindrome arm
            max_arm_length: Maximum length of palindrome arm
            max_spacer: Maximum spacer between arms
            
        Returns:
            List of cruciform motifs
        """
        sequence = sequence.upper()
        motifs = []
        seq_len = len(sequence)
        
        # Limit search for performance
        if seq_len > 10000:
            print("Warning: Cruciform detection limited to first 10000 bp")
            sequence = sequence[:10000]
            seq_len = 10000
        
        # Search for inverted repeats
        for i in range(seq_len - min_arm_length * 2):
            for arm_len in range(min_arm_length, min(max_arm_length + 1, (seq_len - i) // 2)):
                left_arm = sequence[i:i + arm_len]
                
                # Search for reverse complement in downstream region
                for spacer_len in range(0, min(max_spacer + 1, seq_len - i - 2 * arm_len)):
                    right_start = i + arm_len + spacer_len
                    right_end = right_start + arm_len
                    
                    if right_end > seq_len:
                        break
                    
                    right_arm = sequence[right_start:right_end]
                    
                    # Check if reverse complement
                    if right_arm == CruciformDetector.reverse_complement(left_arm):
                        motif_start = i
                        motif_end = right_end
                        motif_seq = sequence[motif_start:motif_end]
                        
                        motifs.append({
                            'Class': 'Cruciform',
                            'Subclass': 'Inverted Repeats',
                            'Type': 'Palindrome',
                            'Start': motif_start + 1,
                            'End': motif_end,
                            'Length': len(motif_seq),
                            'Sequence': motif_seq,
                            'Arm_Length': arm_len,
                            'Spacer_Length': spacer_len,
                            'Scoring_Method': 'cruciform_stability',
                            'Detection_Method': 'Algorithmic'
                        })
        
        return motifs


class APhilicDetector:
    """
    Detect A-philic DNA using tetranucleotide log2-odds scoring.
    
    Based on Vinogradov (2003) tetranucleotide analysis.
    """
    
    # Tetranucleotide log2-odds scores for A-philic DNA
    # Higher scores indicate A-philic character
    TETRANUCLEOTIDE_SCORES = {
        'AAAA': 2.5, 'AAAT': 1.8, 'AAAG': 1.2, 'AAAC': 0.8,
        'AATA': 1.9, 'AATT': 2.0, 'AATG': 0.9, 'AATC': 0.7,
        'AAGA': 1.3, 'AAGT': 1.1, 'AAGG': 0.5, 'AAGC': 0.4,
        'AACA': 1.0, 'AACT': 0.9, 'AACG': 0.3, 'AACC': 0.2,
        'ATAA': 2.1, 'ATAT': 1.7, 'ATAG': 1.0, 'ATAC': 0.8,
        'ATTA': 2.2, 'ATTT': 1.9, 'ATTG': 1.0, 'ATTC': 0.8,
        'ATGA': 1.1, 'ATGT': 1.0, 'ATGG': 0.4, 'ATGC': 0.3,
        'ATCA': 0.9, 'ATCT': 0.8, 'ATCG': 0.2, 'ATCC': 0.1,
        # Add more tetranucleotides with negative or low scores for non-A-philic
        'GGGG': -1.5, 'CCCC': -1.5, 'CGCG': -1.0, 'GCGC': -1.0,
    }
    
    @staticmethod
    def calculate_aphilic_score(sequence: str, window_size: int = 50) -> float:
        """
        Calculate A-philic propensity score.
        
        Args:
            sequence: DNA sequence
            window_size: Window size for scoring
            
        Returns:
            A-philic score
        """
        if len(sequence) < 4:
            return 0.0
        
        scores = []
        
        # Score all tetranucleotides
        for i in range(len(sequence) - 3):
            tetra = sequence[i:i+4]
            score = APhilicDetector.TETRANUCLEOTIDE_SCORES.get(tetra, 0.0)
            scores.append(score)
        
        return np.mean(scores) if scores else 0.0
    
    @staticmethod
    def detect_aphilic_regions(
        sequence: str,
        min_length: int = 20,
        threshold: float = 0.5
    ) -> List[Dict[str, Any]]:
        """
        Detect A-philic DNA regions.
        
        Args:
            sequence: DNA sequence
            min_length: Minimum region length
            threshold: Minimum score threshold
            
        Returns:
            List of A-philic motifs
        """
        sequence = sequence.upper()
        motifs = []
        
        # Sliding window approach
        window_size = min_length
        
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i + window_size]
            score = APhilicDetector.calculate_aphilic_score(window)
            
            if score >= threshold:
                # Extend window if possible
                extended_end = i + window_size
                while extended_end < len(sequence):
                    next_window = sequence[i:extended_end + 1]
                    next_score = APhilicDetector.calculate_aphilic_score(next_window)
                    if next_score >= threshold:
                        extended_end += 1
                    else:
                        break
                
                motif_seq = sequence[i:extended_end]
                
                motifs.append({
                    'Class': 'A-philic DNA',
                    'Subclass': 'A-philic DNA',
                    'Type': 'A-philic tract',
                    'Start': i + 1,
                    'End': extended_end,
                    'Length': len(motif_seq),
                    'Sequence': motif_seq,
                    'Scoring_Method': 'tetranucleotide_score',
                    'Detection_Method': 'Algorithmic',
                    'Aphilic_Score': round(score, 4)
                })
        
        return motifs


class DirectRepeatDetector:
    """
    Detect direct repeats with variable spacers.
    
    Cannot use simple regex due to backreference limitations.
    """
    
    @staticmethod
    def find_direct_repeats(
        sequence: str,
        min_unit_length: int = 5,
        max_unit_length: int = 50,
        max_spacer: int = 200,
        min_copies: int = 2
    ) -> List[Dict[str, Any]]:
        """
        Find direct repeats in sequence.
        
        Args:
            sequence: DNA sequence
            min_unit_length: Minimum repeat unit length
            max_unit_length: Maximum repeat unit length
            max_spacer: Maximum spacer between repeats
            min_copies: Minimum number of copies
            
        Returns:
            List of direct repeat motifs
        """
        sequence = sequence.upper()
        motifs = []
        seq_len = len(sequence)
        
        # Limit for performance
        if seq_len > 50000:
            print("Warning: Direct repeat detection skipped for sequences >50kb")
            return motifs
        
        # Search for direct repeats
        for unit_len in range(min_unit_length, min(max_unit_length + 1, seq_len // 2)):
            for i in range(seq_len - unit_len):
                unit = sequence[i:i + unit_len]
                
                # Look for second copy downstream
                for j in range(i + unit_len, min(i + unit_len + max_spacer, seq_len - unit_len + 1)):
                    if sequence[j:j + unit_len] == unit:
                        spacer_len = j - (i + unit_len)
                        motif_start = i
                        motif_end = j + unit_len
                        motif_seq = sequence[motif_start:motif_end]
                        
                        motifs.append({
                            'Class': 'Slipped DNA',
                            'Subclass': 'Direct Repeat',
                            'Type': 'Direct repeat',
                            'Start': motif_start + 1,
                            'End': motif_end,
                            'Length': len(motif_seq),
                            'Sequence': motif_seq,
                            'Unit_Length': unit_len,
                            'Spacer_Length': spacer_len,
                            'Scoring_Method': 'repeat_score',
                            'Detection_Method': 'Algorithmic'
                        })
        
        return motifs


class HybridMotifDetector:
    """
    Detect hybrid motifs from overlapping different classes.
    """
    
    @staticmethod
    def calculate_overlap_percentage(motif1: Dict, motif2: Dict) -> float:
        """Calculate overlap percentage between two motifs."""
        start1, end1 = motif1['Start'], motif1['End']
        start2, end2 = motif2['Start'], motif2['End']
        
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)
        overlap_len = max(0, overlap_end - overlap_start)
        
        union_len = max(end1, end2) - min(start1, start2)
        
        return overlap_len / union_len if union_len > 0 else 0.0
    
    @staticmethod
    def detect_hybrid_motifs(
        motifs: List[Dict[str, Any]],
        overlap_threshold: Tuple[float, float] = (0.3, 0.7)
    ) -> List[Dict[str, Any]]:
        """
        Detect hybrid motifs from overlapping classes.
        
        Args:
            motifs: List of all detected motifs
            overlap_threshold: (min, max) overlap percentage for hybrid
            
        Returns:
            List of hybrid motifs
        """
        hybrids = []
        
        # Group motifs by position for efficient overlap checking
        sorted_motifs = sorted(motifs, key=lambda m: m['Start'])
        
        for i, motif1 in enumerate(sorted_motifs):
            for motif2 in sorted_motifs[i+1:]:
                # Stop if no more possible overlaps
                if motif2['Start'] > motif1['End']:
                    break
                
                # Check if different classes
                if motif1['Class'] != motif2['Class']:
                    overlap_pct = HybridMotifDetector.calculate_overlap_percentage(motif1, motif2)
                    
                    if overlap_threshold[0] <= overlap_pct <= overlap_threshold[1]:
                        # Create hybrid motif
                        hybrid_start = min(motif1['Start'], motif2['Start'])
                        hybrid_end = max(motif1['End'], motif2['End'])
                        
                        hybrids.append({
                            'Class': 'Hybrid',
                            'Subclass': 'Multi-class Overlap',
                            'Type': f"{motif1['Class']}-{motif2['Class']} hybrid",
                            'Start': hybrid_start,
                            'End': hybrid_end,
                            'Length': hybrid_end - hybrid_start,
                            'Component_Classes': [motif1['Class'], motif2['Class']],
                            'Overlap_Percentage': round(overlap_pct, 3),
                            'Scoring_Method': 'hybrid_score',
                            'Detection_Method': 'Overlap_Analysis'
                        })
        
        return hybrids


class ClusterDetector:
    """
    Detect Non-B DNA clusters (high-density regions).
    """
    
    @staticmethod
    def detect_clusters(
        motifs: List[Dict[str, Any]],
        window_size: int = 1000,
        min_motifs: int = 3,
        min_diversity: int = 2
    ) -> List[Dict[str, Any]]:
        """
        Detect high-density motif clusters.
        
        Args:
            motifs: List of detected motifs
            window_size: Sliding window size
            min_motifs: Minimum motifs in window
            min_diversity: Minimum number of different classes
            
        Returns:
            List of cluster motifs
        """
        if not motifs:
            return []
        
        clusters = []
        
        # Sort by position
        sorted_motifs = sorted(motifs, key=lambda m: m['Start'])
        
        # Sliding window
        for i, motif in enumerate(sorted_motifs):
            window_start = motif['Start']
            window_end = window_start + window_size
            
            # Collect motifs in window
            window_motifs = []
            for m in sorted_motifs[i:]:
                if m['Start'] <= window_end:
                    window_motifs.append(m)
                else:
                    break
            
            # Check criteria
            if len(window_motifs) >= min_motifs:
                classes = set(m['Class'] for m in window_motifs)
                
                if len(classes) >= min_diversity:
                    actual_start = min(m['Start'] for m in window_motifs)
                    actual_end = max(m['End'] for m in window_motifs)
                    
                    clusters.append({
                        'Class': 'Non-B DNA Clusters',
                        'Subclass': 'Motif Hotspot',
                        'Type': 'Mixed cluster',
                        'Start': actual_start,
                        'End': actual_end,
                        'Length': actual_end - actual_start,
                        'Motif_Count': len(window_motifs),
                        'Class_Diversity': len(classes),
                        'Component_Classes': list(classes),
                        'Scoring_Method': 'cluster_score',
                        'Detection_Method': 'Density_Analysis'
                    })
        
        return clusters


def detect_non_hyperscan_motifs(
    sequence: str,
    sequence_name: str = "Sequence",
    existing_motifs: Optional[List[Dict[str, Any]]] = None
) -> List[Dict[str, Any]]:
    """
    Detect all non-Hyperscan motifs.
    
    Args:
        sequence: DNA sequence
        sequence_name: Sequence name
        existing_motifs: Motifs already detected by Hyperscan (for hybrid/cluster)
        
    Returns:
        List of detected motifs
    """
    all_motifs = []
    
    # 1. Cruciform detection
    print(f"  Detecting Cruciform motifs...")
    cruciforms = CruciformDetector.find_inverted_repeats(sequence)
    for motif in cruciforms:
        motif['Sequence_Name'] = sequence_name
    all_motifs.extend(cruciforms)
    
    # 2. A-philic DNA detection
    print(f"  Detecting A-philic DNA...")
    aphilic = APhilicDetector.detect_aphilic_regions(sequence)
    for motif in aphilic:
        motif['Sequence_Name'] = sequence_name
    all_motifs.extend(aphilic)
    
    # 3. Direct Repeats detection (if sequence is not too large)
    if len(sequence) <= 50000:
        print(f"  Detecting Direct Repeats...")
        direct_repeats = DirectRepeatDetector.find_direct_repeats(sequence)
        for motif in direct_repeats:
            motif['Sequence_Name'] = sequence_name
        all_motifs.extend(direct_repeats)
    
    # 4. Hybrid motifs (if existing motifs provided)
    if existing_motifs:
        print(f"  Detecting Hybrid motifs...")
        combined = existing_motifs + all_motifs
        hybrids = HybridMotifDetector.detect_hybrid_motifs(combined)
        for motif in hybrids:
            motif['Sequence_Name'] = sequence_name
        all_motifs.extend(hybrids)
        
        # 5. Cluster detection
        print(f"  Detecting Non-B DNA Clusters...")
        clusters = ClusterDetector.detect_clusters(combined)
        for motif in clusters:
            motif['Sequence_Name'] = sequence_name
        all_motifs.extend(clusters)
    
    return all_motifs


if __name__ == '__main__':
    # Test non-Hyperscan detection
    print("Non-Hyperscan Detection Test")
    print("=" * 80)
    
    # Test sequence with various motifs
    test_sequence = (
        "ATGCATGC" + "AAAAAAAAAA" * 3 +  # A-philic
        "GCATGCAT" +  # Palindrome
        "ATCGATCG" * 2 +  # Direct repeat
        "GGGGTTTTGGGG" +  # Some G4
        "CCCCAACCCC"  # Some i-motif
    )
    
    print(f"\nTest sequence length: {len(test_sequence)} bp")
    
    # Detect motifs
    motifs = detect_non_hyperscan_motifs(test_sequence, "Test")
    
    print(f"\nMotifs detected: {len(motifs)}")
    for motif in motifs:
        print(f"  - {motif['Class']} ({motif['Subclass']})")
        print(f"    Position: {motif['Start']}-{motif['End']}")
        print(f"    Sequence: {motif.get('Sequence', 'N/A')[:30]}...")
