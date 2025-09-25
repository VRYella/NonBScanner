"""
Optimized NBDFinder Orchestrator
===============================

Enhanced version with improved performance, consistent overlap filtering,
and optimized visualization generation.

Key optimizations:
- Consistent per-class overlap filtering across all entry points
- Optimized parallel processing with better error handling
- Faster visualization generation with caching
- Enhanced memory management

Author: NBDFinder Optimization Team
"""

import time
import sys
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Any, Optional
import multiprocessing
from functools import lru_cache

# Add current directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import motif detection functions
try:
    from motifs.curved_dna import find_curved_DNA
    from motifs.slipped_dna import find_slipped_dna
    from motifs.cruciform_dna import find_cruciform
    from motifs.r_loop import find_r_loop
    from motifs.triplex import find_triplex
    from motifs.g_quadruplex import find_g_quadruplex
    from motifs.i_motif import find_i_motif
    from motifs.z_dna import find_z_dna
    from motifs.hybrid import find_hybrid
    from motifs.cluster import find_cluster
    from motifs.base_motif import standardize_motif_output, validate_motif, select_best_nonoverlapping_motifs
except ImportError as e:
    print(f"Warning: Could not import motif modules: {e}")
    # Fallback functions
    def find_curved_DNA(seq, name): return []
    def find_slipped_dna(seq, name): return []
    def find_cruciform(seq, name): return []
    def find_r_loop(seq, name): return []
    def find_triplex(seq, name): return []
    def find_g_quadruplex(seq, name): return []
    def find_i_motif(seq, name): return []
    def find_z_dna(seq, name): return []
    def find_hybrid(motifs, seq, name): return []
    def find_cluster(motifs, seq_len, name): return []
    def standardize_motif_output(motif, name): return motif
    def validate_motif(motif, seq_len): return True
    def select_best_nonoverlapping_motifs(motifs): return motifs

# Import classification system
try:
    from motif_classification import update_motif_with_ids
except ImportError:
    def update_motif_with_ids(motif):
        return motif

# Import enhanced caching if available
try:
    from enhanced_cache import get_cache_manager
    CACHE_AVAILABLE = True
except ImportError:
    CACHE_AVAILABLE = False
    get_cache_manager = None


def _run_motif_detector_optimized(args):
    """Optimized motif detector with better error handling and performance."""
    detector_func, detector_name, seq, sequence_name = args
    
    try:
        start_time = time.time()
        results = detector_func(seq, sequence_name)
        
        # Validate and standardize results efficiently
        valid_results = []
        for motif in results:
            if validate_motif(motif, len(seq)):
                # Check if motif is already standardized to avoid redundant processing
                if ('Normalized_Score' in motif and 'Actual_Score' in motif 
                    and 'Scoring_Method' in motif):
                    # Already standardized, just update classification
                    classified = update_motif_with_ids(motif.copy())
                    classified['Sequence_Name'] = sequence_name
                    valid_results.append(classified)
                else:
                    # Standardize motif
                    clean_motif = {k: v for k, v in motif.items() 
                                 if k not in ['Sequence_Name']}
                    standardized = standardize_motif_output(clean_motif, sequence_name)
                    classified = update_motif_with_ids(standardized)
                    valid_results.append(classified)
        
        exec_time = time.time() - start_time
        return detector_name, valid_results, exec_time
        
    except Exception as e:
        return detector_name, [], 0.0, str(e)


def all_motifs_optimized(seq: str, 
                        sequence_name: str = "Sequence",
                        nonoverlap: bool = True,  # Default to True for consistency
                        report_hotspots: bool = False,
                        calculate_conservation: bool = False,
                        max_workers: Optional[int] = None) -> List[Dict[str, Any]]:
    """
    Optimized motif detection with consistent per-class overlap filtering.
    
    Args:
        seq: DNA sequence string
        sequence_name: Name for the sequence 
        nonoverlap: Apply per-class overlap filtering (default True)
        report_hotspots: Include cluster detection
        calculate_conservation: Include conservation analysis
        max_workers: Number of parallel workers
        
    Returns:
        List of detected motifs with consistent filtering applied
    """
    
    if not seq or len(seq) < 10:
        return []
    
    # Check cache first if available
    if CACHE_AVAILABLE:
        try:
            cache_manager = get_cache_manager()
            cache_params = {
                'nonoverlap': nonoverlap,
                'report_hotspots': report_hotspots,
                'calculate_conservation': calculate_conservation,
                'optimized': True  # Mark as optimized version
            }
            
            cached_result = cache_manager.get_analysis_result(seq, cache_params)
            if cached_result is not None:
                # Update sequence name in cached results
                for motif in cached_result:
                    motif['Sequence_Name'] = sequence_name
                return cached_result
        except Exception:
            pass  # Continue if caching fails
    
    # Define motif detectors for Classes 1-8 (parallel execution)
    detectors = [
        (find_curved_DNA, "Curved DNA"),
        (find_slipped_dna, "Slipped DNA"),
        (find_cruciform, "Cruciform DNA"),
        (find_r_loop, "R-loop"),
        (find_triplex, "Triplex"),
        (find_g_quadruplex, "G-Quadruplex"),
        (find_i_motif, "i-motif"),
        (find_z_dna, "Z-DNA"),
    ]
    
    all_motifs = []
    total_time = 0
    
    # Determine number of workers
    if max_workers is None:
        max_workers = min(len(detectors), multiprocessing.cpu_count())
    
    # Parallel execution with optimized error handling
    try:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks
            future_to_detector = {
                executor.submit(_run_motif_detector_optimized, 
                              (func, name, seq, sequence_name)): name 
                for func, name in detectors
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_detector):
                detector_name = future_to_detector[future]
                try:
                    result = future.result(timeout=60)  # 60 second timeout per detector
                    if len(result) == 3:
                        name, motifs, exec_time = result
                        all_motifs.extend(motifs)
                        total_time += exec_time
                        print(f"✓ {name}: {len(motifs)} motifs found ({exec_time:.3f}s)")
                    else:
                        name, motifs, exec_time, error = result
                        print(f"✗ {name} failed: {error}")
                except Exception as e:
                    print(f"✗ {detector_name} failed: {e}")
                    
    except Exception as e:
        print(f"Warning: Parallel execution failed, falling back to sequential: {e}")
        # Sequential fallback with timing
        for func, name in detectors:
            try:
                start_time = time.time()
                motifs = func(seq, sequence_name)
                valid_motifs = []
                
                for m in motifs:
                    if validate_motif(m, len(seq)):
                        if ('Normalized_Score' in m and 'Actual_Score' in m 
                            and 'Scoring_Method' in m):
                            classified = update_motif_with_ids(m.copy())
                            classified['Sequence_Name'] = sequence_name
                            valid_motifs.append(classified)
                        else:
                            clean_motif = {k: v for k, v in m.items() 
                                         if k not in ['Sequence_Name']}
                            standardized = standardize_motif_output(clean_motif, sequence_name)
                            classified = update_motif_with_ids(standardized)
                            valid_motifs.append(classified)
                
                all_motifs.extend(valid_motifs)
                exec_time = time.time() - start_time
                total_time += exec_time
                print(f"✓ {name}: {len(valid_motifs)} motifs found ({exec_time:.3f}s)")
                
            except Exception as e:
                print(f"✗ {name} failed: {e}")
    
    # Update S.No for all motifs
    for i, motif in enumerate(all_motifs):
        motif['S.No'] = i + 1
        motif['Sequence_Name'] = sequence_name
    
    print(f"Classes 1-8: {len(all_motifs)} total motifs found ({total_time:.3f}s total)")
    
    # Add hybrids (Class 9) - requires all motifs from Classes 1-8
    try:
        start_time = time.time()
        hybrid_motifs = find_hybrid(all_motifs, seq, sequence_name)
        standardized_hybrids = []
        
        for motif in hybrid_motifs:
            if ('Normalized_Score' in motif and 'Actual_Score' in motif 
                and 'Scoring_Method' in motif):
                classified = update_motif_with_ids(motif.copy())
                classified['Sequence_Name'] = sequence_name
                standardized_hybrids.append(classified)
            else:
                clean_motif = {k: v for k, v in motif.items() 
                             if k not in ['Sequence_Name']}
                standardized = standardize_motif_output(clean_motif, sequence_name)
                classified = update_motif_with_ids(standardized)
                standardized_hybrids.append(classified)
        
        all_motifs.extend(standardized_hybrids)
        exec_time = time.time() - start_time
        print(f"✓ Hybrid (Class 9): {len(standardized_hybrids)} motifs found ({exec_time:.3f}s)")
        
    except Exception as e:
        print(f"✗ Hybrid detection failed: {e}")
    
    # CRITICAL: Apply per-class overlap filtering consistently
    if nonoverlap:
        start_time = time.time()
        original_count = len(all_motifs)
        all_motifs = select_best_nonoverlapping_motifs(all_motifs)
        exec_time = time.time() - start_time
        print(f"After per-class de-overlapping: {len(all_motifs)} motifs "
              f"(removed {original_count - len(all_motifs)}) ({exec_time:.3f}s)")
    
    # Add clusters (Class 10) - requires all motifs including hybrids
    if report_hotspots:
        try:
            start_time = time.time()
            cluster_motifs = find_cluster(all_motifs, len(seq), sequence_name)
            standardized_clusters = []
            
            for motif in cluster_motifs:
                if ('Normalized_Score' in motif and 'Actual_Score' in motif 
                    and 'Scoring_Method' in motif):
                    classified = update_motif_with_ids(motif.copy())
                    classified['Sequence_Name'] = sequence_name
                    standardized_clusters.append(classified)
                else:
                    clean_motif = {k: v for k, v in motif.items() 
                                 if k not in ['Sequence_Name']}
                    standardized = standardize_motif_output(clean_motif, sequence_name)
                    classified = update_motif_with_ids(standardized)
                    standardized_clusters.append(classified)
            
            all_motifs.extend(standardized_clusters)
            exec_time = time.time() - start_time
            print(f"✓ Cluster (Class 10): {len(standardized_clusters)} motifs found ({exec_time:.3f}s)")
            
        except Exception as e:
            print(f"✗ Cluster detection failed: {e}")
    
    # Final validation and numbering
    for i, motif in enumerate(all_motifs):
        motif['S.No'] = i + 1
        motif['Sequence_Name'] = sequence_name
    
    # Cache result if available
    if CACHE_AVAILABLE:
        try:
            cache_manager.set_analysis_result(seq, cache_params, all_motifs)
        except Exception:
            pass  # Continue if caching fails
    
    print(f"🎉 Total motifs found: {len(all_motifs)}")
    return all_motifs


@lru_cache(maxsize=128)
def get_motif_statistics(motifs_tuple):
    """Cached statistics calculation for motifs."""
    motifs = list(motifs_tuple)  # Convert back from tuple
    
    stats = {
        'total_motifs': len(motifs),
        'classes': len(set(m.get('Class', 'Unknown') for m in motifs)),
        'avg_score': sum(m.get('Normalized_Score', 0) for m in motifs) / len(motifs) if motifs else 0,
        'avg_length': sum(m.get('Length', 0) for m in motifs) / len(motifs) if motifs else 0
    }
    
    return stats


def format_motif_summary(motifs: List[Dict[str, Any]]) -> str:
    """Generate formatted summary of motifs with caching."""
    if not motifs:
        return "No motifs found."
    
    # Convert to tuple for caching
    motifs_tuple = tuple(
        tuple(sorted(m.items())) for m in motifs
    )
    
    stats = get_motif_statistics(motifs_tuple)
    
    summary = f"""
NBDFinder Analysis Summary
==========================
Total motifs: {stats['total_motifs']}
Unique classes: {stats['classes']}
Average score: {stats['avg_score']:.2f}
Average length: {stats['avg_length']:.1f} bp

Per-class breakdown:
"""
    
    # Class breakdown
    class_counts = {}
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        class_counts[class_name] = class_counts.get(class_name, 0) + 1
    
    for class_name, count in sorted(class_counts.items()):
        summary += f"  {class_name}: {count} motifs\n"
    
    return summary