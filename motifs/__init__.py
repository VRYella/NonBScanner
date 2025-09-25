"""
NBDFinder Motifs Package
========================
Modular Non-B DNA motif detection with 10 classes

Classes:
1. Curved DNA (Global Array, Local Tract)
2. Slipped DNA (Direct Repeat, STR)
3. Cruciform DNA (Inverted Repeat/HairPin)
4. R-loop
5. Triplex (Triplex, Sticky DNA)
6. G-Quadruplex Family (Multiple variants)
7. i-motif Family (Canonical, Relaxed, AC-motif)
8. Z-DNA (Z-DNA, eGZ)
9. Hybrid (Dynamic overlaps)
10. Non-B DNA Cluster Regions (Dynamic hotspots)
"""

# Import all motif detection functions
from .curved_dna import find_curved_DNA
from .slipped_dna import find_slipped_dna
from .cruciform_dna import find_cruciform
from .r_loop import find_r_loop
from .triplex import find_triplex
from .g_quadruplex import find_g_quadruplex
from .i_motif import find_i_motif
from .z_dna import find_z_dna
from .hybrid import find_hybrid
from .cluster import find_cluster

# Import base utilities
from .base_motif import (
    parse_fasta, wrap, gc_content, reverse_complement, is_palindrome,
    overlapping_finditer, validate_motif, standardize_motif_output,
    select_best_nonoverlapping_motifs
)

# Import visualization
from .visualization import create_all_visualizations

# Import Hyperscan performance utilities
from .hyperscan_manager import clear_hyperscan_cache, get_hyperscan_cache_stats

def all_motifs(seq, nonoverlap=False, report_hotspots=False, sequence_name="Sequence", 
               calculate_conservation=False):
    """
    Find all motifs in a sequence using all detection modules
    
    Args:
        seq: DNA sequence string
        nonoverlap: If True, select best non-overlapping motifs per class
        report_hotspots: If True, also report cluster regions
        sequence_name: Name for the sequence
        calculate_conservation: If True, calculate conservation scores
    
    Returns:
        List of standardized motif dictionaries
    """
    import re
    
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE):
        return []
    
    seq = seq.upper()
    
    # Check cache for existing results
    try:
        from enhanced_cache import get_cache_manager
        cache_manager = get_cache_manager()
        
        # Create cache key based on parameters
        cache_params = {
            'nonoverlap': nonoverlap,
            'report_hotspots': report_hotspots,
            'calculate_conservation': calculate_conservation
        }
        
        cached_result = cache_manager.get_analysis_result(seq, cache_params)
        if cached_result is not None:
            # Update sequence name in cached results
            for motif in cached_result:
                motif['Sequence_Name'] = sequence_name
            return cached_result
    except ImportError:
        cache_manager = None
    
    motif_list = []
    
    # Find all motif types
    motif_list.extend(find_curved_DNA(seq, sequence_name))
    motif_list.extend(find_slipped_dna(seq, sequence_name))
    motif_list.extend(find_cruciform(seq, sequence_name))
    motif_list.extend(find_r_loop(seq, sequence_name))
    motif_list.extend(find_triplex(seq, sequence_name))
    motif_list.extend(find_g_quadruplex(seq, sequence_name))
    motif_list.extend(find_i_motif(seq, sequence_name))
    motif_list.extend(find_z_dna(seq, sequence_name))
    
    # Validate motifs
    motif_list = [m for m in motif_list if validate_motif(m, len(seq))]
    
    # Motifs are already standardized by individual functions, just update sequence names
    for i, motif in enumerate(motif_list):
        motif['Sequence_Name'] = sequence_name
        motif['S.No'] = i + 1
    
    # Add hybrids
    motif_list.extend(find_hybrid(motif_list, seq, sequence_name))
    
    # De-overlap per class if requested
    if nonoverlap:
        motif_list = select_best_nonoverlapping_motifs(motif_list)
    
    # Add hotspots if requested
    if report_hotspots:
        motif_list.extend(find_cluster(motif_list, len(seq), sequence_name))
    
    # Ensure all motifs have sequence name
    for m in motif_list:
        if "Sequence_Name" not in m or not m["Sequence_Name"]:
            m["Sequence_Name"] = sequence_name
    
    # Add conservation analysis if requested
    if calculate_conservation:
        from conservation_analysis import calculate_motif_conservation
        
        # Create a motif finder function for conservation analysis
        def motif_finder_func(seq):
            return all_motifs(seq, sequence_name="shuffled", calculate_conservation=False)
        
        motif_list = calculate_motif_conservation(motif_list, seq, motif_finder_func)
    
    # Store results in cache for future use
    if cache_manager is not None:
        try:
            cache_manager.store_analysis_result(seq, motif_list, cache_params)
        except Exception:
            pass  # Don't fail if caching fails
    
    return motif_list


def get_basic_stats(seq, motifs=None):
    """Calculate basic sequence statistics"""
    seq = seq.upper()
    length = len(seq)
    gc = gc_content(seq)
    at = (seq.count('A') + seq.count('T')) / length * 100 if length else 0
    
    stats = {
        "Length": length,
        "GC%": round(gc, 2),
        "AT%": round(at, 2),
        "A": seq.count('A'),
        "T": seq.count('T'),
        "G": seq.count('G'),
        "C": seq.count('C'),
    }
    
    if motifs is not None:
        # Filter out hybrid and cluster motifs for coverage calculations
        filtered_motifs = [m for m in motifs if m.get('Class') not in ['Hybrid', 'Non-B_DNA_Cluster']]
        
        covered = set()
        for m in filtered_motifs:
            covered.update(range(m['Start'], m['End'] + 1))  # Inclusive range
        coverage_pct = (len(covered) / length * 100) if length else 0
        stats["Motif Coverage %"] = round(coverage_pct, 2)
        
        # Add info about excluded motifs
        excluded_count = len(motifs) - len(filtered_motifs)
        stats["Excluded_Motifs_Count"] = excluded_count
        if excluded_count > 0:
            excluded_types = set(m.get('Class', 'Unknown') for m in motifs if m.get('Class') in ['Hybrid', 'Non-B_DNA_Cluster'])
            stats["Excluded_Motif_Types"] = list(excluded_types)
    
    return stats


def format_motif_rows(motifs):
    """Format motifs for output with standardized column order"""
    ordered = []
    for m in motifs:
        row = {
            "S.No": m.get("S.No", ""),
            "Sequence_Name": m.get("Sequence_Name", ""),
            "Chromosome/Contig": m.get("Chromosome/Contig", ""),
            "Class": m.get("Class", ""),
            "Subclass": m.get("Subclass", ""),
            "Motif_ID": m.get("Motif_ID", ""),
            "Start": m.get("Start", ""),
            "End": m.get("End", ""),
            "Length": m.get("Length", ""),
            "Normalized_Score": m.get("Normalized_Score", ""),
            "Actual_Score": m.get("Actual_Score", ""),
            "Scoring_Method": m.get("Scoring_Method", ""),
            "GC_Content": m.get("GC_Content", ""),
            "Sequence": m.get("Sequence", ""),
            "Overlap_Classes": m.get("Overlap_Classes", "")
        }
        ordered.append(row)
    return ordered


# Version info
__version__ = "1.0.0"
__author__ = "Dr. Venkata Rajesh Yella"


# Export main functions
__all__ = [
    # Main functions
    'all_motifs',
    'get_basic_stats',
    'format_motif_rows',
    
    # Individual motif finders
    'find_curved_DNA',
    'find_slipped_dna', 
    'find_cruciform',
    'find_r_loop',
    'find_triplex',
    'find_g_quadruplex',
    'find_i_motif',
    'find_z_dna',
    'find_hybrid',
    'find_cluster',
    
    # Utilities
    'parse_fasta',
    'wrap',
    'gc_content',
    'reverse_complement',
    'is_palindrome',
    'overlapping_finditer',
    'validate_motif',
    'standardize_motif_output',
    'select_best_nonoverlapping_motifs',
    
    # Visualization
    'create_all_visualizations',
    
    # Performance utilities
    'clear_hyperscan_cache',
    'get_hyperscan_cache_stats'
]