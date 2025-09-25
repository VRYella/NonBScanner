"""
Z-DNA Motif Detection (Class 8) -- Hyperscan Accelerated

SCIENTIFIC BASIS:
================
Z-DNA is a left-handed double helical form of DNA discovered by Alexander Rich.
It forms under superhelical tension and high salt conditions, characterized by:
- Left-handed helical structure (vs. right-handed B-DNA)
- Alternating purine-pyrimidine sequences favor Z-form transition
- CG dinucleotides are particularly prone to Z-DNA formation
- Role in gene regulation, recombination, and chromatin structure

BIOLOGICAL SIGNIFICANCE:
- Gene expression regulation (Rich & Zhang, 2003)
- Chromatin organization and nucleosome positioning
- Recombination hotspots and genome instability
- Association with active transcription sites

SUBCLASSES DETECTED:
===================
1. Z-DNA (8.1): Classical Z-forming sequences with CG/AT dinucleotides
2. eGZ (Extruded-G) DNA (8.3): CGG repeat expansions forming slipped-out structures

SCORING ALGORITHMS:
==================
Z-seeker algorithm (Ho et al. 1986, Wang et al. 2007):
- Dinucleotide-based scoring with experimentally validated weights
- CG/GC: +7.0 (strong Z-forming potential)
- AT/TA: +0.5 (weak Z-forming potential)  
- GT/TG, AC/CA: +1.25 (moderate Z-forming potential)
- Consecutive AT penalty to avoid false positives
- Sliding window approach with thresholding

HYPERSCAN ACCELERATION:
======================
Uses Hyperscan for CGG repeat detection, Python regex for complex Z-seeker scoring.
Output format: 1-based coordinates for genomic pipeline compatibility.
"""

import hyperscan, numpy as np, re
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from motifs.base_motif import wrap, gc_content, standardize_motif_output
from motifs.hyperscan_manager import optimized_hs_find
from core.regex_registry import get_patterns_for_motif

def standardize_candidate_output(candidate: dict, sequence_name: str = "", motif_id: int = 0) -> dict:
    """
    Standardize candidate output for detection-only mode (no scoring fields).
    
    Args:
        candidate: Candidate motif dict
        sequence_name: Name of the sequence
        motif_id: Motif ID number
        
    Returns:
        Standardized candidate dict without scoring fields
    """
    seq = candidate.get("Sequence", "").replace('\n', '')
    motif_class = candidate.get("Class", "")
    subclass = candidate.get("Subclass", candidate.get("Subtype", ""))
    motif_length = candidate.get("Length", len(seq))
    
    # Import classification system for motif IDs
    try:
        from motif_classification import get_motif_id
        class_motif_id = get_motif_id(subclass)
    except ImportError:
        class_motif_id = "0.0"
    
    return {
        "S.No": motif_id,
        "Sequence_Name": sequence_name,
        "Chromosome/Contig": "",
        "Class": motif_class,
        "Subclass": subclass,
        "Motif_ID": f"{motif_class}_{class_motif_id}_{candidate.get('Start', '')}-{candidate.get('End', '')}",
        "Start": candidate.get("Start", ""),
        "End": candidate.get("End", ""),
        "Length": motif_length,
        "GC_Content": round(gc_content(seq), 2) if seq else 0.0,
        "Sequence": wrap(seq),
        "Overlap_Classes": "",
    }

# === Load patterns from central registry ===
ZDNA_PATTERNS = get_patterns_for_motif('z_dna')

# === Z-DNA Seeker Scoring Algorithm (Ho 1986, Rich 1993, Wang 2007) ===
def zdna_seeker_scoring_array(seq, GC_weight=7.0, AT_weight=0.5, GT_weight=1.25, AC_weight=1.25,
        consecutive_AT_scoring=(0.5, 0.5, 0.5, 0.5, 0.0, 0.0, -5.0, -100.0),
        mismatch_penalty_type="linear",
        mismatch_penalty_starting_value=3,
        mismatch_penalty_linear_delta=3,
        cadence_reward=0.0):
    """
    Generate Z-DNA propensity scores for every dinucleotide in sequence.
    
    SCIENTIFIC BASIS:
    - Experimentally validated dinucleotide weights from crystallographic studies
    - CG/GC dinucleotides have highest Z-forming potential (weight=7.0)
    - AT dinucleotides have weak but positive Z-forming potential (weight=0.5)
    - Consecutive AT sequences penalized to avoid false positives in AT-rich regions
    - Mismatch penalties for non-canonical dinucleotides
    
    ALGORITHM PARAMETERS:
    - GC_weight: Score for CG/GC dinucleotides (experimentally: 7.0)
    - AT_weight: Base score for AT/TA dinucleotides (0.5)
    - consecutive_AT_scoring: Progressive penalty for consecutive AT runs
    - mismatch_penalty: Penalty for non-Z-forming dinucleotides
    
    Returns: NumPy array of per-dinucleotide Z-forming scores
    """
    # BLOCK: Dinucleotide scoring with experimental weights and AT-run penalties
    scoring_array = np.empty(len(seq) - 1, dtype=float); mismatches_counter=0; consecutive_AT_counter=0
    for i in range(len(seq) - 1):
        t = seq[i:i+2].upper()
        if t in ("GC", "CG"): scoring_array[i]=GC_weight; mismatches_counter=0; consecutive_AT_counter=0
        elif t in ("GT", "TG"): scoring_array[i]=GT_weight; mismatches_counter=0; consecutive_AT_counter=0
        elif t in ("AC", "CA"): scoring_array[i]=AC_weight; mismatches_counter=0; consecutive_AT_counter=0
        elif t in ("AT", "TA"):
            adjusted_weight=AT_weight
            adjusted_weight+= consecutive_AT_scoring[consecutive_AT_counter] if consecutive_AT_counter<len(consecutive_AT_scoring) else consecutive_AT_scoring[-1]
            scoring_array[i]=adjusted_weight; consecutive_AT_counter+=1; mismatches_counter=0
        else:
            mismatches_counter+=1; consecutive_AT_counter=0
            if mismatch_penalty_type=="exponential":
                scoring_array[i]=-(mismatch_penalty_starting_value**mismatches_counter if mismatches_counter<15 else 32000.0)
            elif mismatch_penalty_type=="linear":
                scoring_array[i]=-mismatch_penalty_starting_value-mismatch_penalty_linear_delta*(mismatches_counter-1)
            else:
                scoring_array[i]=-10.0
        if t in ("GC","CG","GT","TG","AC","CA","AT","TA"): scoring_array[i]+=cadence_reward
    return scoring_array

#--- Z-DNA motif finder using Z-seeker algorithm (literature-based, sliding window score) ---
def find_zdna_candidates(seq, threshold=50, drop_threshold=50, **kwargs) -> list:
    """
    Find Z-DNA candidate regions using Z-seeker algorithm (detection only, no scoring).
    
    Returns candidate regions without scores for later scoring by centralized scorer.
    """
    seq=seq.upper(); candidates=[]; n=len(seq)
    if n<12: return []
    scoring=zdna_seeker_scoring_array(seq, **kwargs)
    start_idx=0; max_ending_here=scoring[0]; current_max=0; candidate=None; end_idx=1
    for i in range(1, len(scoring)):
        num=scoring[i]
        if num>=max_ending_here+num: start_idx=i; end_idx=i+1; max_ending_here=num
        else: max_ending_here+=num; end_idx=i+1
        if max_ending_here>=threshold and (candidate is None or current_max<max_ending_here):
            candidate=(start_idx,end_idx,max_ending_here); current_max=max_ending_here
        if candidate and (max_ending_here<0 or current_max-max_ending_here>=drop_threshold):
            s,e,score=candidate
            # Return candidate without score - scoring will be done separately
            candidates.append({"Class":"Z-DNA","Subclass":"Z-DNA","Start":s+1,"End":e+1,"Length":e-s+1,
                           "Sequence":wrap(seq[s:e+1])}); candidate=None; max_ending_here=current_max=0
    if candidate:
        s,e,score=candidate
        candidates.append({"Class":"Z-DNA","Subclass":"Z-DNA","Start":s+1,"End":e+1,"Length":e-s+1,
                       "Sequence":wrap(seq[s:e+1])})
    return candidates

# Backward compatibility: keep original function that includes scoring
def find_zdna(seq, threshold=50, drop_threshold=50, **kwargs) -> list:
    # Block: Use Z-seeker score to find Z-DNA regions, per Ho 1986, Wang 2007, with scientific sliding-window thresholding
    seq=seq.upper(); motifs=[]; n=len(seq)
    if n<12: return []
    scoring=zdna_seeker_scoring_array(seq, **kwargs)
    start_idx=0; max_ending_here=scoring[0]; current_max=0; candidate=None; end_idx=1
    for i in range(1, len(scoring)):
        num=scoring[i]
        if num>=max_ending_here+num: start_idx=i; end_idx=i+1; max_ending_here=num
        else: max_ending_here+=num; end_idx=i+1
        if max_ending_here>=threshold and (candidate is None or current_max<max_ending_here):
            candidate=(start_idx,end_idx,max_ending_here); current_max=max_ending_here
        if candidate and (max_ending_here<0 or current_max-max_ending_here>=drop_threshold):
            s,e,score=candidate
            motifs.append({"Class":"Z-DNA","Subclass":"Z-DNA","Start":s+1,"End":e+1,"Length":e-s+1,
                           "Sequence":wrap(seq[s:e+1]),"Score":float(score),"ScoreMethod":"Z-Seeker"}); candidate=None; max_ending_here=current_max=0
    if candidate:
        s,e,score=candidate
        motifs.append({"Class":"Z-DNA","Subclass":"Z-DNA","Start":s+1,"End":e+1,"Length":e-s+1,
                       "Sequence":wrap(seq[s:e+1]),"Score":float(score),"ScoreMethod":"Z-Seeker"})
    return motifs

#--- eGZ motif finder: Extruded-G (CGG)n, using Hyperscan for block-motif scan ---
def find_egz_candidates(seq) -> list:
    """Find eGZ/CGG repeat candidate motifs using optimized Hyperscan scanning (detection only)."""
    if not seq:
        return []
    
    seqU = seq.upper()
    
    # Prepare pattern for optimized Hyperscan
    patterns = [
        (r"(CGG){4,}", 1)
    ]
    
    def egz_callback(id, start, end, flags, ctx, pattern):
        """Optimized callback for eGZ motif detection (candidates only)."""
        motif_seq = seqU[start:end]
        
        return {
            "Class": "Z-DNA", "Subclass": "eGZ", "Start": start+1, "End": end,
            "Length": len(motif_seq), "Sequence": wrap(motif_seq)
        }
    
    # Use optimized Hyperscan manager
    return optimized_hs_find(patterns, seqU, egz_callback)

# Backward compatibility: keep original function that includes scoring
def find_egz_motif(seq) -> list:
    """Find eGZ/CGG repeat motifs using optimized Hyperscan scanning."""
    if not seq:
        return []
    
    seqU = seq.upper()
    
    # Prepare pattern for optimized Hyperscan
    patterns = [
        (r"(CGG){4,}", 1)
    ]
    
    def egz_callback(id, start, end, flags, ctx, pattern):
        """Optimized callback for eGZ motif detection."""
        motif_seq = seqU[start:end]
        n_repeats = len(motif_seq) // 3
        g_frac = motif_seq.count('G') / len(motif_seq)
        score = n_repeats * 3 * (1.0 + 2.0 * g_frac)
        
        return {
            "Class": "Z-DNA", "Subclass": "eGZ", "Start": start+1, "End": end,
            "Length": len(motif_seq), "Sequence": wrap(motif_seq),
            "ScoreMethod": "Repeat_raw", "Score": float(score), "CGG_Repeats": n_repeats
        }
    
    # Use optimized Hyperscan manager
    return optimized_hs_find(patterns, seqU, egz_callback)

#--- Main: Find all Z-DNA and eGZ motifs, output standardized for genomic analysis ---
def find_z_dna(seq: str, sequence_name: str = "") -> list:
    zdna_results=find_zdna(seq); egz_results=find_egz_motif(seq)
    all_results=zdna_results+egz_results
    return [standardize_motif_output(motif, sequence_name, i) for i, motif in enumerate(all_results, 1)]

#--- New: Find candidate Z-DNA regions (detection only, no scoring) ---
def find_z_dna_candidates(seq: str, sequence_name: str = "") -> list:
    """
    Find Z-DNA candidate regions without scoring.
    
    Returns list of candidate regions for later scoring by centralized scorer.
    """
    zdna_candidates = find_zdna_candidates(seq)
    egz_candidates = find_egz_candidates(seq)
    all_candidates = zdna_candidates + egz_candidates
    return [standardize_candidate_output(candidate, sequence_name, i) for i, candidate in enumerate(all_candidates, 1)]

#--- Annotations ---
# - zdna_seeker_scoring_array: core scoring array, weights/penalties from Z-DNA literature.
# - find_zdna_candidates: detection-only Z-seeker, returns candidates without scores for centralized scoring.
# - find_zdna: maximal-scoring subsequence (Z-seeker), with scientific threshold and drop logic (backward compatibility).
# - find_egz_candidates: detection-only Hyperscan scan for eGZ (CGG)n, returns candidates without scores.
# - find_egz_motif: Hyperscan block scan for eGZ (CGG)n, literature length/cutoff, scoring G-bias (backward compatibility).
# - find_z_dna_candidates: combines detection-only functions, output standardized 1-based for centralized scoring.
# - find_z_dna: combines both scoring functions, output standardized 1-based for analysis (backward compatibility).
