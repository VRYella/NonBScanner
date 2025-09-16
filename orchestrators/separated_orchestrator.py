"""
NBDFinder Separated Orchestrator - Detection then Scoring
========================================================

New orchestrator that implements separated detection and scoring architecture:
1. Run detection modules first to get all candidates (no scoring)
2. Apply centralized scoring to all candidates 
3. Add hybrids and clusters with normalized scoring only
4. Maintain backward compatibility with existing interfaces

This implements the separated architecture while preserving scientific accuracy
and providing both raw and normalized scores.

Author: NBDFinder Team
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

# Import detection functions (candidate-only versions)
try:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from detectors.class08_zdna import find_z_dna_candidates
    # TODO: Add other detection modules as they get separated
    # from detectors.class01_curved import find_curved_dna_candidates
    # from detectors.class02_slipped import find_slipped_dna_candidates
    # ... etc
except ImportError as e:
    print(f"Warning: Could not import detection modules: {e}")
    find_z_dna_candidates = None

# Import scoring system
try:
    from core.scoring import create_default_scorer
except ImportError as e:
    print(f"Warning: Could not import scoring system: {e}")
    create_default_scorer = None

# Import backward compatibility functions
try:
    from motifs.curved_dna import find_curved_DNA
    from motifs.slipped_dna import find_slipped_DNA  # Note: corrected function name
    from motifs.cruciform_dna import find_cruciform
    from motifs.r_loop import find_r_loop
    from motifs.triplex import find_triplex
    from motifs.g_quadruplex import find_g_quadruplex
    from motifs.i_motif import find_i_motif
    from motifs.z_dna import find_z_dna
    from motifs.hybrid import find_hybrid
    from motifs.cluster import find_cluster
    from motifs.base_motif import validate_motif, standardize_motif_output, select_best_nonoverlapping_motifs
except ImportError as e:
    print(f"Warning: Could not import motif functions: {e}")
    # Set fallback None values
    find_curved_DNA = find_slipped_DNA = find_cruciform = None
    find_r_loop = find_triplex = find_g_quadruplex = None 
    find_i_motif = find_z_dna = find_hybrid = find_cluster = None
    validate_motif = standardize_motif_output = select_best_nonoverlapping_motifs = None

# Import classification system
try:
    from motif_classification import update_motif_with_ids
except ImportError:
    def update_motif_with_ids(motif):
        return motif

def _run_detection_module(args):
    """Helper function to run detection module in parallel."""
    detection_func, seq, sequence_name, detector_name = args
    try:
        return detector_name, detection_func(seq, sequence_name)
    except Exception as e:
        print(f"Error in {detector_name}: {e}")
        return detector_name, []

def separated_motif_detection(seq: str, 
                            sequence_name: str = "Sequence",
                            nonoverlap: bool = False,
                            report_hotspots: bool = False,
                            calculate_conservation: bool = False,
                            max_workers: Optional[int] = None,
                            use_separated_architecture: bool = True) -> List[Dict[str, Any]]:
    """
    Unified detection API using separated detection-then-scoring architecture.
    
    Args:
        seq: DNA sequence string
        sequence_name: Name for the sequence
        nonoverlap: If True, select best non-overlapping motifs per class
        report_hotspots: If True, also report cluster regions
        calculate_conservation: If True, calculate conservation scores
        max_workers: Maximum number of parallel workers (None = auto)
        use_separated_architecture: If True, use new separated approach; if False, use legacy
        
    Returns:
        List of standardized motif dictionaries with scores
    """
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE):
        return []
    
    seq = seq.upper()
    
    if not use_separated_architecture:
        # Fallback to legacy combined detection+scoring
        return legacy_all_motifs(seq, sequence_name, nonoverlap, report_hotspots, calculate_conservation, max_workers)
    
    print(f"Running separated detection-then-scoring architecture on {len(seq)} bp sequence")
    
    if max_workers is None:
        max_workers = min(8, multiprocessing.cpu_count())
    
    # Phase 1: Detection (candidates only, no scoring)
    print("Phase 1: Detecting candidate motifs...")
    
    # Available detection modules (candidate-only versions)
    candidate_detectors = []
    
    # Z-DNA detection (separated)
    if 'find_z_dna_candidates' in globals():
        candidate_detectors.append((find_z_dna_candidates, "Z-DNA"))
    
    # TODO: Add other separated detectors as they become available
    # if 'find_curved_dna_candidates' in globals():
    #     candidate_detectors.append((find_curved_dna_candidates, "Curved DNA"))
    
    # For now, use legacy detectors for classes that haven't been separated yet
    legacy_detectors = [
        (find_curved_DNA, "Curved DNA"),
        (find_slipped_dna, "Slipped DNA"), 
        (find_cruciform, "Cruciform"),
        (find_r_loop, "R-Loop"),
        (find_triplex, "Triplex"),
        (find_g_quadruplex, "G-Quadruplex"),
        (find_i_motif, "i-Motif")
    ]
    
    all_candidates = []
    
    # Run separated detection modules
    if candidate_detectors:
        try:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # Submit detection jobs
                future_to_detector = {
                    executor.submit(_run_detection_module, (func, seq, sequence_name, name)): name 
                    for func, name in candidate_detectors
                }
                
                # Collect detection results
                for future in as_completed(future_to_detector):
                    detector_name = future_to_detector[future]
                    try:
                        detector_name_result, candidates = future.result()
                        all_candidates.extend(candidates)
                        print(f"  ✓ {detector_name}: {len(candidates)} candidates detected")
                    except Exception as e:
                        print(f"  ✗ {detector_name} detection failed: {e}")
        
        except Exception as e:
            print(f"Warning: Parallel detection failed, falling back to sequential: {e}")
            for func, name in candidate_detectors:
                try:
                    candidates = func(seq, sequence_name)
                    all_candidates.extend(candidates)
                    print(f"  ✓ {name}: {len(candidates)} candidates detected")
                except Exception as e:
                    print(f"  ✗ {name} detection failed: {e}")
    
    # Run legacy detectors (with built-in scoring)
    legacy_motifs = []
    if legacy_detectors:
        print("  Running legacy detectors (with built-in scoring)...")
        for func, name in legacy_detectors:
            try:
                motifs = func(seq, sequence_name)
                valid_motifs = []
                for m in motifs:
                    if validate_motif(m, len(seq)):
                        clean_motif = {k: v for k, v in m.items() if k not in ['Sequence_Name']}
                        standardized = standardize_motif_output(clean_motif, sequence_name)
                        classified = update_motif_with_ids(standardized)
                        valid_motifs.append(classified)
                legacy_motifs.extend(valid_motifs)
                print(f"  ✓ {name}: {len(valid_motifs)} motifs found")
            except Exception as e:
                print(f"  ✗ {name} failed: {e}")
    
    print(f"Phase 1 complete: {len(all_candidates)} candidates + {len(legacy_motifs)} legacy motifs")
    
    # Phase 2: Centralized Scoring
    print("Phase 2: Applying centralized scoring...")
    
    scored_motifs = []
    
    if all_candidates and create_default_scorer:
        try:
            scorer = create_default_scorer()
            scored_candidates = scorer.score_candidates(all_candidates, seq)
            
            # Standardize and classify scored candidates
            for i, motif in enumerate(scored_candidates):
                clean_motif = {k: v for k, v in motif.items() if k not in ['Sequence_Name']}
                standardized = standardize_motif_output(clean_motif, sequence_name)
                classified = update_motif_with_ids(standardized)
                scored_motifs.append(classified)
            
            print(f"  ✓ Centralized scoring: {len(scored_motifs)} motifs scored")
        except Exception as e:
            print(f"  ✗ Centralized scoring failed: {e}")
            # Fallback: treat candidates as legacy motifs
            scored_motifs.extend(all_candidates)
    else:
        # No centralized scoring available
        scored_motifs.extend(all_candidates)
    
    # Combine all motifs
    all_motifs = scored_motifs + legacy_motifs
    
    # Update S.No for all motifs
    for i, motif in enumerate(all_motifs):
        motif['S.No'] = i + 1
        motif['Sequence_Name'] = sequence_name
    
    print(f"Phase 2 complete: {len(all_motifs)} total scored motifs")
    
    # Phase 3: Add Hybrids and Clusters (normalized scoring only)
    print("Phase 3: Adding hybrids and clusters...")
    
    # Add hybrids (Class 9)
    try:
        hybrid_motifs = find_hybrid(all_motifs, seq, sequence_name)
        standardized_hybrids = []
        for motif in hybrid_motifs:
            clean_motif = {k: v for k, v in motif.items() if k not in ['Sequence_Name']}
            standardized = standardize_motif_output(clean_motif, sequence_name)
            classified = update_motif_with_ids(standardized)
            standardized_hybrids.append(classified)
        
        all_motifs.extend(standardized_hybrids)
        print(f"  ✓ Hybrid (Class 9): {len(standardized_hybrids)} motifs found")
    except Exception as e:
        print(f"  ✗ Hybrid detection failed: {e}")
    
    # Add clusters (Class 10)
    if report_hotspots:
        try:
            cluster_motifs = find_cluster(all_motifs, len(seq), sequence_name)
            standardized_clusters = []
            for motif in cluster_motifs:
                clean_motif = {k: v for k, v in motif.items() if k not in ['Sequence_Name']}
                standardized = standardize_motif_output(clean_motif, sequence_name)
                classified = update_motif_with_ids(standardized)
                standardized_clusters.append(classified)
            
            all_motifs.extend(standardized_clusters)
            print(f"  ✓ Cluster (Class 10): {len(standardized_clusters)} motifs found")
        except Exception as e:
            print(f"  ✗ Cluster detection failed: {e}")
    
    # De-overlap per class if requested
    if nonoverlap:
        print("Applying non-overlap selection...")
        all_motifs = select_best_nonoverlapping_motifs(all_motifs)
        print(f"  ✓ Non-overlap: {len(all_motifs)} motifs after de-overlap")
    
    print(f"🎉 Separated architecture complete: {len(all_motifs)} total motifs")
    return all_motifs


def legacy_all_motifs(seq: str, 
                     sequence_name: str = "Sequence",
                     nonoverlap: bool = False,
                     report_hotspots: bool = False,
                     calculate_conservation: bool = False,
                     max_workers: Optional[int] = None) -> List[Dict[str, Any]]:
    """
    Legacy detection API using combined detection+scoring (backward compatibility).
    """
    # This is essentially the old orchestrator behavior
    # Import the legacy orchestrator function
    try:
        from orchestrators.all_motifs import all_motifs_refactored
        return all_motifs_refactored(seq, sequence_name, nonoverlap, report_hotspots, calculate_conservation, max_workers)
    except ImportError:
        print("Warning: Legacy orchestrator not available")
        return []


# Backward compatibility alias
all_motifs_separated = separated_motif_detection