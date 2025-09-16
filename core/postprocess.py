"""
Post-processing for Motif Detection Results
==========================================

Handles non-overlapping filtering, priority-based merging, deduplication,
and final result consolidation for NBDFinder motif detection pipeline.
"""

import numpy as np
from typing import List, Dict, Any, Tuple, Set, Optional
from collections import defaultdict
import itertools

def remove_overlapping_motifs(motifs: List[Dict[str, Any]], 
                             priority_order: List[str] = None) -> List[Dict[str, Any]]:
    """
    Remove overlapping motifs based on priority and score.
    
    Args:
        motifs: List of motif dictionaries with Start, End, Class, Score
        priority_order: Order of motif classes by priority (higher priority first)
        
    Returns:
        Filtered list with non-overlapping motifs
    """
    if not motifs:
        return []
    
    # Default priority order
    if priority_order is None:
        priority_order = [
            "Curved_DNA", "Slipped_DNA", "Cruciform", "R-Loop", "Triplex",
            "G-Quadruplex", "i-Motif", "Z-DNA", "Hybrid", "Cluster"
        ]
    
    # Create priority mapping
    priority_map = {cls: i for i, cls in enumerate(priority_order)}
    
    # Sort motifs by priority, then by score (descending)
    def sort_key(motif):
        cls = motif.get('Class', 'Unknown')
        priority = priority_map.get(cls, len(priority_order))
        score = motif.get('Score', 0)
        return (priority, -score)  # Negative score for descending order
    
    sorted_motifs = sorted(motifs, key=sort_key)
    
    # Greedy non-overlapping selection
    selected = []
    for motif in sorted_motifs:
        if not _overlaps_with_any(motif, selected):
            selected.append(motif)
    
    # Sort final result by position
    return sorted(selected, key=lambda x: x['Start'])

def _overlaps_with_any(motif: Dict[str, Any], motif_list: List[Dict[str, Any]]) -> bool:
    """Check if motif overlaps with any motif in the list."""
    for other in motif_list:
        if _motifs_overlap(motif, other):
            return True
    return False

def _motifs_overlap(motif1: Dict[str, Any], motif2: Dict[str, Any]) -> bool:
    """Check if two motifs overlap in their genomic positions."""
    start1, end1 = motif1['Start'], motif1['End']
    start2, end2 = motif2['Start'], motif2['End']
    
    # No overlap if one ends before the other starts
    return not (end1 <= start2 or end2 <= start1)

def merge_nearby_motifs(motifs: List[Dict[str, Any]], 
                       max_distance: int = 10,
                       same_class_only: bool = True) -> List[Dict[str, Any]]:
    """
    Merge motifs that are close together.
    
    Args:
        motifs: List of motif dictionaries
        max_distance: Maximum distance between motifs to merge
        same_class_only: Only merge motifs of the same class
        
    Returns:
        List with nearby motifs merged
    """
    if not motifs:
        return []
    
    # Sort by position
    sorted_motifs = sorted(motifs, key=lambda x: x['Start'])
    merged = []
    
    i = 0
    while i < len(sorted_motifs):
        current = sorted_motifs[i]
        merge_group = [current]
        
        # Look for nearby motifs to merge
        j = i + 1
        while j < len(sorted_motifs):
            next_motif = sorted_motifs[j]
            
            # Check distance
            distance = next_motif['Start'] - current['End']
            if distance > max_distance:
                break
            
            # Check class compatibility
            if same_class_only and current['Class'] != next_motif['Class']:
                j += 1
                continue
            
            merge_group.append(next_motif)
            current = next_motif  # Update current for distance calculation
            j += 1
        
        # Merge the group or add single motif
        if len(merge_group) > 1:
            merged_motif = _merge_motif_group(merge_group)
            merged.append(merged_motif)
        else:
            merged.append(current)
        
        i = j if len(merge_group) > 1 else i + 1
    
    return merged

def _merge_motif_group(motif_group: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Merge a group of nearby motifs into a single motif."""
    if not motif_group:
        return {}
    
    if len(motif_group) == 1:
        return motif_group[0]
    
    # Calculate merged properties
    merged = motif_group[0].copy()
    
    # Span of all motifs
    merged['Start'] = min(m['Start'] for m in motif_group)
    merged['End'] = max(m['End'] for m in motif_group)
    merged['Length'] = merged['End'] - merged['Start']
    
    # Average/max scores
    scores = [m.get('Score', 0) for m in motif_group]
    merged['Score'] = max(scores)  # Use highest score
    merged['Merged_Count'] = len(motif_group)
    
    # Combine sequences if available
    if all('Sequence' in m for m in motif_group):
        # This is approximate - actual sequence would need parent sequence
        merged['Sequence'] = f"MERGED_{len(motif_group)}_MOTIFS"
    
    return merged

def deduplicate_motifs(motifs: List[Dict[str, Any]], 
                      tolerance: int = 5) -> List[Dict[str, Any]]:
    """
    Remove duplicate motifs (same class/subclass at similar positions).
    
    Args:
        motifs: List of motif dictionaries
        tolerance: Position tolerance for considering motifs duplicates
        
    Returns:
        Deduplicated motif list
    """
    if not motifs:
        return []
    
    # Group by class and subclass
    groups = defaultdict(list)
    for motif in motifs:
        key = (motif.get('Class', ''), motif.get('Subclass', ''))
        groups[key].append(motif)
    
    deduplicated = []
    
    for group in groups.values():
        if len(group) == 1:
            deduplicated.extend(group)
            continue
        
        # Sort by position
        sorted_group = sorted(group, key=lambda x: x['Start'])
        unique_motifs = [sorted_group[0]]
        
        for current in sorted_group[1:]:
            is_duplicate = False
            
            for existing in unique_motifs:
                if _are_duplicate_motifs(current, existing, tolerance):
                    # Keep the one with higher score
                    if current.get('Score', 0) > existing.get('Score', 0):
                        unique_motifs.remove(existing)
                        unique_motifs.append(current)
                    is_duplicate = True
                    break
            
            if not is_duplicate:
                unique_motifs.append(current)
        
        deduplicated.extend(unique_motifs)
    
    return sorted(deduplicated, key=lambda x: x['Start'])

def _are_duplicate_motifs(motif1: Dict[str, Any], motif2: Dict[str, Any], 
                         tolerance: int) -> bool:
    """Check if two motifs are duplicates within tolerance."""
    start_diff = abs(motif1['Start'] - motif2['Start'])
    end_diff = abs(motif1['End'] - motif2['End'])
    
    return start_diff <= tolerance and end_diff <= tolerance

def filter_by_score_threshold(motifs: List[Dict[str, Any]], 
                             min_score: float = 0.1) -> List[Dict[str, Any]]:
    """
    Filter motifs by minimum score threshold.
    
    Args:
        motifs: List of motif dictionaries
        min_score: Minimum score threshold
        
    Returns:
        Filtered motif list
    """
    return [m for m in motifs if m.get('Score', 0) >= min_score]

def filter_by_length_constraints(motifs: List[Dict[str, Any]], 
                                length_limits: Dict[str, Tuple[int, int]] = None) -> List[Dict[str, Any]]:
    """
    Filter motifs by class-specific length constraints.
    
    Args:
        motifs: List of motif dictionaries
        length_limits: Dictionary mapping class names to (min_len, max_len) tuples
        
    Returns:
        Filtered motif list
    """
    if not length_limits:
        return motifs
    
    filtered = []
    for motif in motifs:
        cls = motif.get('Class', '')
        length = motif.get('Length', motif.get('End', 0) - motif.get('Start', 0))
        
        if cls in length_limits:
            min_len, max_len = length_limits[cls]
            if min_len <= length <= max_len:
                filtered.append(motif)
        else:
            # Include motifs with unknown classes
            filtered.append(motif)
    
    return filtered

def calculate_motif_statistics(motifs: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Calculate summary statistics for a motif list.
    
    Args:
        motifs: List of motif dictionaries
        
    Returns:
        Dictionary of statistics
    """
    if not motifs:
        return {
            'total_motifs': 0,
            'unique_classes': 0,
            'total_length': 0,
            'average_score': 0.0,
            'score_range': (0.0, 0.0)
        }
    
    # Basic counts
    total_motifs = len(motifs)
    classes = set(m.get('Class', 'Unknown') for m in motifs)
    unique_classes = len(classes)
    
    # Length statistics
    lengths = [m.get('Length', m.get('End', 0) - m.get('Start', 0)) for m in motifs]
    total_length = sum(lengths)
    
    # Score statistics
    scores = [m.get('Score', 0) for m in motifs]
    average_score = np.mean(scores) if scores else 0.0
    score_range = (min(scores), max(scores)) if scores else (0.0, 0.0)
    
    # Class distribution
    class_counts = defaultdict(int)
    for motif in motifs:
        class_counts[motif.get('Class', 'Unknown')] += 1
    
    return {
        'total_motifs': total_motifs,
        'unique_classes': unique_classes,
        'class_distribution': dict(class_counts),
        'total_length': total_length,
        'average_length': np.mean(lengths) if lengths else 0.0,
        'average_score': float(average_score),
        'score_range': score_range,
        'total_coverage_bp': total_length
    }

def apply_all_postprocessing(motifs: List[Dict[str, Any]], 
                           config: Dict[str, Any] = None) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """
    Apply complete post-processing pipeline to motifs.
    
    Args:
        motifs: Raw motif list
        config: Configuration dictionary with processing parameters
        
    Returns:
        Tuple of (processed_motifs, statistics)
    """
    if config is None:
        config = {}
    
    original_count = len(motifs)
    
    # Step 1: Filter by score threshold
    min_score = config.get('min_score_threshold', 0.1)
    motifs = filter_by_score_threshold(motifs, min_score)
    after_score_filter = len(motifs)
    
    # Step 2: Filter by length constraints
    length_limits = config.get('length_limits', {})
    motifs = filter_by_length_constraints(motifs, length_limits)
    after_length_filter = len(motifs)
    
    # Step 3: Deduplicate
    tolerance = config.get('duplicate_tolerance', 5)
    motifs = deduplicate_motifs(motifs, tolerance)
    after_dedup = len(motifs)
    
    # Step 4: Merge nearby motifs (optional)
    if config.get('merge_nearby', False):
        max_distance = config.get('merge_distance', 10)
        motifs = merge_nearby_motifs(motifs, max_distance)
        after_merge = len(motifs)
    else:
        after_merge = len(motifs)
    
    # Step 5: Remove overlaps based on priority
    priority_order = config.get('priority_order', None)
    motifs = remove_overlapping_motifs(motifs, priority_order)
    final_count = len(motifs)
    
    # Calculate final statistics
    final_stats = calculate_motif_statistics(motifs)
    
    # Add processing statistics
    processing_stats = {
        'original_count': original_count,
        'after_score_filter': after_score_filter,
        'after_length_filter': after_length_filter,
        'after_dedup': after_dedup,
        'after_merge': after_merge,
        'final_count': final_count,
        'filters_applied': {
            'score_threshold': min_score,
            'length_constraints': bool(length_limits),
            'deduplication': True,
            'merge_nearby': config.get('merge_nearby', False),
            'overlap_removal': True
        }
    }
    
    final_stats.update(processing_stats)
    
    return motifs, final_stats

__all__ = [
    'remove_overlapping_motifs',
    'merge_nearby_motifs',
    'deduplicate_motifs',
    'filter_by_score_threshold',
    'filter_by_length_constraints',
    'calculate_motif_statistics',
    'apply_all_postprocessing'
]