"""
Base motif detection utilities and common functions.
Shared across all motif detection modules.
"""

import re
import numpy as np
from typing import Dict, List, Any
import sys
import os

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from classification_config import normalize_score, get_motif_limits
except ImportError:
    # Fallback for when classification_config is not available
    def normalize_score(actual_score, motif_length, motif_class, subclass=None):
        """Fallback normalization"""
        return min(100.0, max(0.0, actual_score))
    
    def get_motif_limits(motif_class, subclass=None):
        """Fallback limits"""
        return 10, 200


def parse_fasta(fasta_str: str) -> str:
    """Parse FASTA string to clean DNA sequence"""
    return "".join([line.strip() for line in fasta_str.split('\n') if not line.startswith(">")]).upper().replace(" ", "").replace("U", "T")


def wrap(seq: str, width: int = 60) -> str:
    """Format sequence with line breaks"""
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


def gc_content(seq: str) -> float:
    """Calculate GC content percentage"""
    if not seq:
        return 0.0
    gc = seq.count('G') + seq.count('C')
    return (gc / len(seq)) * 100


def reverse_complement(seq: str) -> str:
    """Generate reverse complement of DNA sequence"""
    return seq.translate(str.maketrans("ATGC", "TACG"))[::-1]


def is_palindrome(seq: str) -> bool:
    """Check if sequence is palindromic"""
    return seq == reverse_complement(seq)


def overlapping_finditer(pattern, seq):
    """Find all overlapping matches of pattern in sequence"""
    regex = re.compile(pattern, re.IGNORECASE)
    pos = 0
    while pos < len(seq):
        m = regex.search(seq, pos)
        if not m:
            break
        yield m
        pos = m.start() + 1


def validate_motif(motif: Dict[str, Any], seq_length: int) -> bool:
    """Validate that a motif has required fields and valid coordinates"""
    required_keys = ["Class", "Subclass", "Start", "End", "Length", "Sequence"]
    if not all(key in motif for key in required_keys):
        return False
    if not (1 <= motif["Start"] <= motif["End"] <= seq_length):
        return False
    if len(motif["Sequence"].replace('\n', '')) == 0:
        return False
    return True


def standardize_motif_output(motif: Dict[str, Any], sequence_name: str = "", motif_id: int = 0) -> Dict[str, Any]:
    """Standardize motif output to required format with normalization and official classification"""
    seq = motif.get("Sequence", "").replace('\n', '')
    actual_score = motif.get("Score", 0)
    motif_class = motif.get("Class", "")
    subclass = motif.get("Subclass", motif.get("Subtype", ""))
    motif_length = motif.get("Length", len(seq))
    
    # Calculate normalized score
    try:
        normalized_score = normalize_score(
            float(actual_score) if actual_score else 0.0,
            motif_length,
            motif_class,
            subclass
        )
    except (ValueError, TypeError):
        normalized_score = 0.0
    
    # Import classification system
    try:
        from motif_classification import get_official_subclass_name, get_motif_id
        official_subclass = get_official_subclass_name(subclass)
        class_motif_id = get_motif_id(subclass)
    except ImportError:
        # Fallback if classification module not available
        official_subclass = subclass
        class_motif_id = "0.0"
    
    return {
        "S.No": motif_id,
        "Sequence_Name": sequence_name,
        "Chromosome/Contig": "",  # To be filled by caller if available
        "Class": motif_class,
        "Subclass": official_subclass,
        "Motif_ID": f"{motif_class}_{class_motif_id}_{motif.get('Start', '')}-{motif.get('End', '')}",
        "Start": motif.get("Start", ""),
        "End": motif.get("End", ""),
        "Length": motif_length,
        "Normalized_Score": normalized_score,
        "Actual_Score": actual_score,
        "Scoring_Method": motif.get("ScoreMethod", "unknown"),
        "GC_Content": round(gc_content(seq), 2) if seq else 0.0,
        "Sequence": wrap(seq),
        "Overlap_Classes": "",  # To be filled by overlap analysis
    }


def select_best_nonoverlapping_motifs(motifs: List[Dict], motif_priority: List[str] = None) -> List[Dict]:
    """
    Enhanced overlap filtering for maximum specificity and minimal redundancy.
    
    Uses official classification system for priority ordering and implements
    stringent overlap filtering with clustering of nearby motifs.
    
    Args:
        motifs: List of motif dictionaries
        motif_priority: Optional priority list (uses official G4 priority if None)
    
    Returns:
        List of selected high-specificity motifs with minimal overlaps and redundancy
    """
    if not motifs:
        return []
    
    # Apply clustering to merge nearby redundant motifs first
    clustered_motifs = cluster_nearby_motifs(motifs)
    
    # Get official priority order from classification system
    if motif_priority is None:
        try:
            from motif_classification import MOTIF_CLASSES
            # Find G-Quadruplex class with priority order
            for class_info in MOTIF_CLASSES.values():
                if class_info.get("class_name") == "G-Quadruplex Family" and "priority_order" in class_info:
                    motif_priority = class_info["priority_order"]
                    break
        except ImportError:
            pass
        
        # Fallback to official names (without underscores)
        if motif_priority is None:
            motif_priority = [
                'Multimeric G4', 'Canonical G4', 'Relaxed G4', 'Bulged G4',
                'Bipartite G4', 'Imperfect G4', 'G-Triplex intermediate'
            ]
    
    # Create priority ranking
    subtype_rank = {subtype: i for i, subtype in enumerate(motif_priority)}
    
    def normalize_subclass_name(subclass):
        """Convert current implementation names to official names"""
        try:
            from motif_classification import CURRENT_TO_OFFICIAL
            return CURRENT_TO_OFFICIAL.get(subclass, subclass)
        except ImportError:
            # Manual mapping as fallback for common G4 types
            mapping = {
                'Multimeric_G4': 'Multimeric G4',
                'Canonical_G4': 'Canonical G4', 
                'Relaxed_G4': 'Relaxed G4',
                'Bulged_G4': 'Bulged G4',
                'Bipartite_G4': 'Bipartite G4',
                'Imperfect_G4': 'Imperfect G4',
                'G-Triplex_intermediate': 'G-Triplex intermediate'
            }
            return mapping.get(subclass, subclass)
    
    def motif_key(m):
        # Get subclass, handling both Subclass and Subtype fields
        raw_subclass = m.get('Subclass', m.get('Subtype', ''))
        normalized_subclass = normalize_subclass_name(raw_subclass)
        
        # Get priority rank
        rank = subtype_rank.get(normalized_subclass, len(subtype_rank))
        
        # Get score with proper priority: Normalized_Score > Score > Actual_Score
        try:
            score = float(m.get('Normalized_Score', m.get('Score', m.get('Actual_Score', 0))))
        except (ValueError, TypeError):
            score = 0.0
        
        length = m.get('Length', 0)
        
        # Return sort key: (Class, -Score, Priority_Rank, -Length)
        # Score is prioritized over subclass rank for within-class selection
        return (m.get('Class', ''), -score, rank, -length)
    
    # Sort motifs by priority (class, then score, then rank, then length)
    sorted_motifs = sorted(clustered_motifs, key=motif_key)
    
    # Enhanced overlap filtering with inter-class conflict resolution
    selected = []
    occupied_per_class = dict()
    global_occupied = set()  # Track all occupied positions globally
    
    for m in sorted_motifs:
        motif_class = m.get('Class', 'Other')
        
        # Validate coordinates
        start = m.get('Start', 0)
        end = m.get('End', 0)
        if start <= 0 or end <= 0 or start > end:
            continue
            
        # Create position range (inclusive)
        region = set(range(start, end + 1))
        
        # Initialize class tracking if needed
        if motif_class not in occupied_per_class:
            occupied_per_class[motif_class] = set()
        
        # Check for overlap within the same class (strict) and globally (lenient)
        intra_class_overlap = not occupied_per_class[motif_class].isdisjoint(region)
        
        # Allow some inter-class overlap for hybrid detection, but limit excessive overlap
        inter_class_overlap_ratio = len(region.intersection(global_occupied)) / len(region)
        excessive_overlap = inter_class_overlap_ratio > 0.5  # Allow up to 50% inter-class overlap
        
        if not intra_class_overlap and not excessive_overlap:
            selected.append(m)
            occupied_per_class[motif_class].update(region)
            global_occupied.update(region)
    
    return selected


def cluster_nearby_motifs(motifs: List[Dict], max_distance: int = 10) -> List[Dict]:
    """
    Cluster nearby motifs of the same class/subclass to reduce redundancy.
    Merges motifs that are within max_distance of each other.
    """
    if not motifs:
        return motifs
    
    # Group motifs by class and subclass
    grouped = {}
    for motif in motifs:
        key = (motif.get('Class', ''), motif.get('Subclass', ''))
        if key not in grouped:
            grouped[key] = []
        grouped[key].append(motif)
    
    clustered_result = []
    
    for (cls, subcls), group in grouped.items():
        if len(group) <= 1:
            clustered_result.extend(group)
            continue
            
        # Sort by position
        group.sort(key=lambda x: (x.get('Start', 0), x.get('End', 0)))
        
        clusters = []
        current_cluster = [group[0]]
        
        for motif in group[1:]:
            # Check if this motif is close to the last one in current cluster
            last_end = current_cluster[-1].get('End', 0)
            current_start = motif.get('Start', 0)
            
            if current_start - last_end <= max_distance:
                current_cluster.append(motif)
            else:
                # Start new cluster
                clusters.append(current_cluster)
                current_cluster = [motif]
        
        clusters.append(current_cluster)
        
        # For each cluster, keep only the best motif (highest score)
        for cluster in clusters:
            if len(cluster) == 1:
                clustered_result.append(cluster[0])
            else:
                # Choose best motif based on score
                best_motif = max(cluster, key=lambda x: float(x.get('Normalized_Score', x.get('Score', x.get('Actual_Score', 0)))))
                clustered_result.append(best_motif)
    
    return clustered_result


def select_best_nonoverlapping_motifs(motifs: List[Dict], motif_priority: List[str] = None) -> List[Dict]:
    """
    Enhanced overlap filtering for maximum specificity and minimal redundancy.
    
    Uses official classification system for priority ordering and implements
    stringent overlap filtering with clustering of nearby motifs.
    
    Args:
        motifs: List of motif dictionaries
        motif_priority: Optional priority list (uses official G4 priority if None)
    
    Returns:
        List of selected high-specificity motifs with minimal overlaps and redundancy
    """
    if not motifs:
        return []
    
    # Apply clustering to merge nearby redundant motifs first
    clustered_motifs = cluster_nearby_motifs(motifs)
    
    # Get official priority order from classification system
    if motif_priority is None:
        try:
            from motif_classification import MOTIF_CLASSES
            # Find G-Quadruplex class with priority order
            for class_info in MOTIF_CLASSES.values():
                if class_info.get("class_name") == "G-Quadruplex Family" and "priority_order" in class_info:
                    motif_priority = class_info["priority_order"]
                    break
        except ImportError:
            pass
        
        # Fallback to official names (without underscores)
        if motif_priority is None:
            motif_priority = [
                'Multimeric G4', 'Canonical G4', 'Relaxed G4', 'Bulged G4',
                'Bipartite G4', 'Imperfect G4', 'G-Triplex intermediate'
            ]
    
    # Create priority ranking
    subtype_rank = {subtype: i for i, subtype in enumerate(motif_priority)}
    
    def normalize_subclass_name(subclass):
        """Convert current implementation names to official names"""
        try:
            from motif_classification import CURRENT_TO_OFFICIAL
            return CURRENT_TO_OFFICIAL.get(subclass, subclass)
        except ImportError:
            # Manual mapping as fallback for common G4 types
            mapping = {
                'Multimeric_G4': 'Multimeric G4',
                'Canonical_G4': 'Canonical G4', 
                'Relaxed_G4': 'Relaxed G4',
                'Bulged_G4': 'Bulged G4',
                'Bipartite_G4': 'Bipartite G4',
                'Imperfect_G4': 'Imperfect G4',
                'G-Triplex_intermediate': 'G-Triplex intermediate'
            }
            return mapping.get(subclass, subclass)
    
    def motif_key(m):
        # Get subclass, handling both Subclass and Subtype fields
        raw_subclass = m.get('Subclass', m.get('Subtype', ''))
        normalized_subclass = normalize_subclass_name(raw_subclass)
        
        # Get priority rank
        rank = subtype_rank.get(normalized_subclass, len(subtype_rank))
        
        # Get score with proper priority: Normalized_Score > Score > Actual_Score
        try:
            score = float(m.get('Normalized_Score', m.get('Score', m.get('Actual_Score', 0))))
        except (ValueError, TypeError):
            score = 0.0
        
        length = m.get('Length', 0)
        
        # Return sort key: (Class, -Score, Priority_Rank, -Length)
        # Score is prioritized over subclass rank for within-class selection
        return (m.get('Class', ''), -score, rank, -length)
    
    # Sort motifs by priority (class, then score, then rank, then length)
    sorted_motifs = sorted(clustered_motifs, key=motif_key)
    
    # Enhanced overlap filtering with inter-class conflict resolution
    selected = []
    occupied_per_class = dict()
    global_occupied = set()  # Track all occupied positions globally
    
    for m in sorted_motifs:
        motif_class = m.get('Class', 'Other')
        
        # Validate coordinates
        start = m.get('Start', 0)
        end = m.get('End', 0)
        if start <= 0 or end <= 0 or start > end:
            continue
            
        # Create position range (inclusive)
        region = set(range(start, end + 1))
        
        # Initialize class tracking if needed
        if motif_class not in occupied_per_class:
            occupied_per_class[motif_class] = set()
        
        # Check for overlap within the same class (strict) and globally (lenient)
        intra_class_overlap = not occupied_per_class[motif_class].isdisjoint(region)
        
        # Allow some inter-class overlap for hybrid detection, but limit excessive overlap
        inter_class_overlap_ratio = len(region.intersection(global_occupied)) / len(region)
        excessive_overlap = inter_class_overlap_ratio > 0.5  # Allow up to 50% inter-class overlap
        
        if not intra_class_overlap and not excessive_overlap:
            selected.append(m)
            occupied_per_class[motif_class].update(region)
            global_occupied.update(region)
    
    return selected