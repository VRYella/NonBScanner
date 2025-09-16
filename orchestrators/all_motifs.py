"""
NBDFinder Unified Orchestrator - all_motifs_refactored.py
========================================================

Unified detection API that runs all motif detectors (Classes 1–8) in parallel
using ProcessPoolExecutor for maximum performance while preserving scientific
accuracy and independent scoring systems.

Key Features:
- Parallel execution of motif detectors using ProcessPoolExecutor
- Standardized results with standardize_motif_output 
- Official 10-class, 22-subclass taxonomy mapping via classification_config
- Automatic addition of hybrids (Class 9) and clusters (Class 10)
- Hyperscan integration for fast candidate discovery
- Independent scientific scoring systems preserved

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

import re
import sys
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Any, Optional
import multiprocessing

# Add current directory to path for imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

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
    def standardize_motif_output(motif, name, idx): return motif
    def validate_motif(motif, seq_len): return True
    def select_best_nonoverlapping_motifs(motifs): return motifs

# Import classification system
try:
    from motif_classification import update_motif_with_ids
except ImportError:
    def update_motif_with_ids(motif): return motif

# Import enhanced caching if available
try:
    from enhanced_cache import get_cache_manager
    CACHE_AVAILABLE = True
except ImportError:
    CACHE_AVAILABLE = False
    get_cache_manager = None


def _run_motif_detector(args):
    """
    Worker function to run a single motif detector.
    This function is executed in parallel by ProcessPoolExecutor.
    """
    detector_func, seq, sequence_name, detector_name = args
    
    try:
        # Run the detector
        results = detector_func(seq, sequence_name)
        
        # Validate results
        valid_results = []
        for motif in results:
            if validate_motif(motif, len(seq)):
                # Ensure standardization and official classification
                # Remove any existing fields that might conflict to avoid duplicates
                clean_motif = {k: v for k, v in motif.items() if k not in ['Sequence_Name']}
                standardized = standardize_motif_output(clean_motif, sequence_name)
                classified = update_motif_with_ids(standardized)
                valid_results.append(classified)
        
        return detector_name, valid_results
        
    except Exception as e:
        print(f"Warning: {detector_name} detector failed: {e}")
        return detector_name, []


def all_motifs_refactored(seq: str, 
                         sequence_name: str = "Sequence",
                         nonoverlap: bool = False,
                         report_hotspots: bool = False,
                         calculate_conservation: bool = False,
                         max_workers: Optional[int] = None) -> List[Dict[str, Any]]:
    """
    Unified orchestrator for all motif detection using parallel processing.
    
    Args:
        seq: DNA sequence string
        sequence_name: Name for the sequence
        nonoverlap: If True, select best non-overlapping motifs per class
        report_hotspots: If True, also report cluster regions (Class 10)
        calculate_conservation: If True, calculate conservation scores
        max_workers: Maximum number of parallel workers (default: CPU count)
    
    Returns:
        List of standardized motif dictionaries with official classification
    """
    # Input validation
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE):
        return []
    
    seq = seq.upper()
    
    # Check cache for existing results
    if CACHE_AVAILABLE:
        try:
            cache_manager = get_cache_manager()
            cache_params = {
                'nonoverlap': nonoverlap,
                'report_hotspots': report_hotspots,
                'calculate_conservation': calculate_conservation,
                'refactored': True  # Mark as refactored version
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
        (find_z_dna, "Z-DNA")
    ]
    
    # Prepare arguments for parallel execution
    detector_args = [(func, seq, sequence_name, name) for func, name in detectors]
    
    # Set default max_workers
    if max_workers is None:
        max_workers = min(multiprocessing.cpu_count(), len(detectors))
    
    all_motifs = []
    
    # Execute motif detectors in parallel (Classes 1-8)
    try:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit all detector jobs
            future_to_detector = {
                executor.submit(_run_motif_detector, args): args[3] 
                for args in detector_args
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_detector):
                detector_name = future_to_detector[future]
                try:
                    detector_name_result, motifs = future.result()
                    all_motifs.extend(motifs)
                    print(f"✓ {detector_name}: {len(motifs)} motifs found")
                except Exception as e:
                    print(f"✗ {detector_name} failed: {e}")
    
    except Exception as e:
        print(f"Warning: Parallel execution failed, falling back to sequential: {e}")
        # Fallback to sequential execution
        for func, name in detectors:
            try:
                motifs = func(seq, sequence_name)
                valid_motifs = []
                for m in motifs:
                    if validate_motif(m, len(seq)):
                        # Remove conflicting fields to avoid duplicates
                        clean_motif = {k: v for k, v in m.items() if k not in ['Sequence_Name']}
                        standardized = standardize_motif_output(clean_motif, sequence_name)
                        classified = update_motif_with_ids(standardized)
                        valid_motifs.append(classified)
                all_motifs.extend(valid_motifs)
                print(f"✓ {name}: {len(valid_motifs)} motifs found")
            except Exception as e:
                print(f"✗ {name} failed: {e}")
    
    # Update S.No for all motifs
    for i, motif in enumerate(all_motifs):
        motif['S.No'] = i + 1
        motif['Sequence_Name'] = sequence_name
    
    print(f"Classes 1-8: {len(all_motifs)} total motifs found")
    
    # Add hybrids (Class 9) - requires all motifs from Classes 1-8
    try:
        hybrid_motifs = find_hybrid(all_motifs, seq, sequence_name)
        standardized_hybrids = []
        for motif in hybrid_motifs:
            # Remove conflicting fields to avoid duplicates
            clean_motif = {k: v for k, v in motif.items() if k not in ['Sequence_Name']}
            standardized = standardize_motif_output(clean_motif, sequence_name)
            classified = update_motif_with_ids(standardized)
            standardized_hybrids.append(classified)
        
        all_motifs.extend(standardized_hybrids)
        print(f"✓ Hybrid (Class 9): {len(standardized_hybrids)} motifs found")
    except Exception as e:
        print(f"✗ Hybrid detection failed: {e}")
    
    # De-overlap per class if requested
    if nonoverlap:
        original_count = len(all_motifs)
        all_motifs = select_best_nonoverlapping_motifs(all_motifs)
        print(f"After per-class de-overlapping: {len(all_motifs)} motifs (removed {original_count - len(all_motifs)})")
    
    # Add clusters (Class 10) - requires all motifs including hybrids
    if report_hotspots:
        try:
            cluster_motifs = find_cluster(all_motifs, len(seq), sequence_name)
            standardized_clusters = []
            for motif in cluster_motifs:
                # Remove conflicting fields to avoid duplicates
                clean_motif = {k: v for k, v in motif.items() if k not in ['Sequence_Name']}
                standardized = standardize_motif_output(clean_motif, sequence_name)
                classified = update_motif_with_ids(standardized)
                standardized_clusters.append(classified)
            
            all_motifs.extend(standardized_clusters)
            print(f"✓ Cluster (Class 10): {len(standardized_clusters)} motifs found")
        except Exception as e:
            print(f"✗ Cluster detection failed: {e}")
    
    # Final validation and sequence name update
    for motif in all_motifs:
        if "Sequence_Name" not in motif or not motif["Sequence_Name"]:
            motif["Sequence_Name"] = sequence_name
    
    # Add conservation analysis if requested
    if calculate_conservation:
        try:
            from conservation_analysis import calculate_motif_conservation
            
            def motif_finder_func(seq):
                return all_motifs_refactored(seq, "shuffled", calculate_conservation=False)
            
            all_motifs = calculate_motif_conservation(all_motifs, seq, motif_finder_func)
            print(f"✓ Conservation analysis completed")
        except Exception as e:
            print(f"✗ Conservation analysis failed: {e}")
    
    # Store results in cache for future use
    if CACHE_AVAILABLE:
        try:
            cache_manager = get_cache_manager()
            cache_params = {
                'nonoverlap': nonoverlap,
                'report_hotspots': report_hotspots,
                'calculate_conservation': calculate_conservation,
                'refactored': True
            }
            cache_manager.store_analysis_result(seq, all_motifs, cache_params)
        except Exception:
            pass  # Don't fail if caching fails
    
    print(f"🎉 Total motifs found: {len(all_motifs)}")
    return all_motifs


# Convenience functions for compatibility
def get_basic_stats(seq, motifs=None):
    """Calculate basic sequence statistics"""
    from motifs import get_basic_stats as _get_basic_stats
    return _get_basic_stats(seq, motifs)


def format_motif_rows(motifs):
    """Format motifs for output with standardized column order"""
    from motifs import format_motif_rows as _format_motif_rows
    return _format_motif_rows(motifs)


# Export main function
__all__ = ['all_motifs_refactored', 'get_basic_stats', 'format_motif_rows']


if __name__ == "__main__":
    # Test the orchestrator
    test_seq = "GGGAGGGAGGGAGGGATCGATCGATCGAAAAAAAAA" * 3
    print(f"Testing orchestrator with sequence length: {len(test_seq)}")
    
    result = all_motifs_refactored(
        test_seq, 
        sequence_name="Test", 
        report_hotspots=True,
        nonoverlap=False
    )
    
    print(f"\nResults:")
    for i, motif in enumerate(result[:5]):  # Show first 5 motifs
        print(f"{i+1}. {motif.get('Class', 'Unknown')} - {motif.get('Subclass', 'Unknown')} "
              f"at {motif.get('Start', '?')}-{motif.get('End', '?')}")
    
    if len(result) > 5:
        print(f"... and {len(result) - 5} more motifs")