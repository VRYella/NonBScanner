"""
Curved DNA Motif Detection (Class 1) - 2024 Literature Aligned, Hyperscan Accelerated
====================================================================================

Detects Poly(A)/Poly(T) tracts and phased global arrays using Hyperscan for speed.
Scores tracts and arrays for DNA curvature, in line with recent literature (see notes below).
- Class 1.1: Global Array (phased/polytract arrays, ~10bp period)
- Class 1.2: Local Tract (isolated long A/T tract, not in array)

References:
- Marini et al, Cell 28:871-879 (1982)
- Crothers et al, Methods Enzymol 212:3-29 (1992)
- Olson et al, PNAS 95:11163 (1998)
- Yella & Bansal, Sci Rep 7:42564 (2017)
- Wang et al, NAR 49:e49 (2021) (for scoring/normalization)
"""

import numpy as np; import re; import hyperscan
from .base_motif import wrap, standardize_motif_output
from .hyperscan_manager import optimized_hs_find
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from core.regex_registry import get_patterns_for_motif

# === Load patterns from central registry ===
CURVED_PATTERNS = get_patterns_for_motif('curved_dna')

# --- Parameter Table ---
"""
| Parameter           | Type   | Description                                                        | Example/Range      |
|---------------------|--------|--------------------------------------------------------------------|--------------------|
| seq                 | str    | DNA sequence to analyze                                            | 'AACCTTAAA...'     |
| min_len             | int    | Minimum tract length for local motif                               | 7 (default)        |
| min_tract_len       | int    | Minimum tract length for global motif                              | 3 (default)        |
| min_repeats         | int    | Minimum number of phased tracts in a global array                  | 3 (default)        |
| min_spacing         | int    | Minimum allowed spacing between phased tracts (bp)                 | 8 (default)        |
| max_spacing         | int    | Maximum allowed spacing between phased tracts (bp)                 | 12 (default)       |
| min_global_score    | float  | Minimum normalized score to report global (array) motif            | 0.2 (default)      |
| min_local_score     | float  | Minimum normalized score to report local (tract) motif             | 0.2 (default)      |
"""

# --- Poly(A)/T tract finder with regex fallback ---
def find_polyA_polyT_tracts(seq: str, min_len: int = 7) -> list:
    """Find poly(A) or poly(T) tracts in sequence. Uses regex fallback for reliability."""
    if not seq:
        return []
    
    seq = seq.upper()
    
    # Try Hyperscan first, fallback to regex if it fails
    try:
        # Prepare patterns for optimized Hyperscan - use registry patterns
        patterns = []
        for pattern_info in CURVED_PATTERNS['poly_a_tracts']:
            patterns.append((pattern_info[0].replace('{7,}', f'{{{min_len},}}'), pattern_info[1]))
        for pattern_info in CURVED_PATTERNS['poly_t_tracts']:
            patterns.append((pattern_info[0].replace('{7,}', f'{{{min_len},}}'), pattern_info[1]))
        
        if not patterns:
            # Fallback to original patterns if registry is empty
            patterns = [
                (f"A{{{min_len},}}", 1),
                (f"T{{{min_len},}}", 2)
            ]
        
        def tract_callback(id, from_, to, flags, ctx, pattern):
            """Optimized callback for tract detection."""
            tract_seq = seq[from_:to]
            tract_type = 'A' if id == 1 else 'T'
            return (from_, to-1, tract_seq, tract_type)
        
        # Use optimized Hyperscan manager
        results = optimized_hs_find(patterns, seq, tract_callback)
        if results and any(r is not None for r in results):
            return [r for r in results if r is not None]
    except Exception:
        pass  # Fall through to regex fallback
    
    # Regex fallback implementation
    results = []
    
    # Find A-tracts
    for match in re.finditer(f'A{{{min_len},}}', seq):
        tract_seq = match.group()
        results.append((match.start(), match.end() - 1, tract_seq, 'A'))
    
    # Find T-tracts  
    for match in re.finditer(f'T{{{min_len},}}', seq):
        tract_seq = match.group()
        results.append((match.start(), match.end() - 1, tract_seq, 'T'))
    
    # Sort by start position
    results.sort(key=lambda x: x[0])
    return results

# --- Updated: Degree-based curvature scoring, AT-run and multitract bonus, normalized ---
def curvature_score(seq: str, tracts=None):
    """
    Returns (raw_score, normalized_score) for the input region.
    - Raw score: sum of dinucleotide bending (Crothers/Olson) + AT-tract and multi-tract bonus.
    - Normalized: 0 (no bend) to 1 (most curved in biologically plausible range).
    """
    seq = seq.upper(); n = len(seq)
    # Dinucleotide bending angles (Crothers/Olson, degrees)
    bend_angles = {'AA':18.9,'AT':14.6,'AG':8.0,'AC':7.2,'TA':16.9,'TT':18.9,'TG':6.1,'TC':8.0,
                   'GA':3.6,'GT':7.2,'GG':5.1,'GC':2.1,'CA':6.1,'CT':3.6,'CG':2.1,'CC':5.1}
    # Sum all bends
    total_bend = sum(bend_angles.get(seq[i:i+2],0.0) for i in range(n-1))
    # AT-tract bonus (Trifonov, Yella&Bansal): reward for each AT run >=4bp, power law for length
    at_runs = re.findall(r"(A{4,}|T{4,})", seq)
    at_bonus = sum(len(run)**1.15 for run in at_runs) * 2.5
    # If known tracts (global array), count for multi-tract bonus
    multi_bonus = 0.0
    if tracts is not None and len(tracts) > 1:
        # Each additional tract adds bonus, as in Yella & Bansal 2017, Wang et al 2021
        multi_bonus = (len(tracts)-1)**1.3 * 12.0
    raw_score = total_bend + at_bonus + multi_bonus
    # Normalization: typical local tract ~100, global arrays up to ~300+
    max_curvature = 350; min_curvature = 0
    norm_score = min(1.0, max(0.0, (raw_score-min_curvature)/(max_curvature-min_curvature)))
    return raw_score, norm_score

# --- Global (Phased Array) finder: phased tracts (~10bp period, 3+ repeats) ---
def find_global_curved_polyA_polyT(seq: str, min_tract_len: int = 3, min_repeats: int = 3,
                                   min_spacing: int = 8, max_spacing: int = 12,
                                   min_global_score: float = 0.2) -> tuple:
    """Find global curved DNA motifs (phased A/T tracts, spaced ~10bp, 3+ in array)"""
    tracts = find_polyA_polyT_tracts(seq, min_tract_len)
    results, apr_regions = [], []
    for i in range(len(tracts)-min_repeats+1):
        group = [tracts[i]]
        for j in range(1, min_repeats):
            prev_center = (tracts[i+j-1][0]+tracts[i+j-1][1])//2
            curr_center = (tracts[i+j][0]+tracts[i+j][1])//2
            spacing = curr_center-prev_center
            if min_spacing <= spacing <= max_spacing: group.append(tracts[i+j])
            else: break
        # Extend group if more phased tracts found
        k = i+len(group)
        while k < len(tracts):
            prev_center = (tracts[k-1][0]+tracts[k-1][1])//2
            curr_center = (tracts[k][0]+tracts[k][1])//2
            spacing = curr_center-prev_center
            if min_spacing <= spacing <= max_spacing: group.append(tracts[k]); k+=1
            else: break
        if len(group) >= min_repeats:
            motif_seq = seq[group[0][0]:group[-1][1]+1]
            raw, norm = curvature_score(motif_seq, group)
            if norm >= min_global_score:
                motif = {
                    "Class":"Curved_DNA","Subclass":"Global_Array",
                    "Start":group[0][0]+1,"End":group[-1][1]+1,
                    "Length":group[-1][1]-group[0][0]+1,
                    "Sequence":wrap(motif_seq),"Score":round(raw,2),"NormScore":round(norm,3),
                    "ScoreMethod":"bend+multiAT","NumTracts":len(group)
                }
                results.append(motif)
                apr_regions.append((motif["Start"], motif["End"]))
    return results, apr_regions

# --- Local tract finder: isolated long tracts (not in global arrays) ---
def find_local_curved_polyA_polyT(seq: str, apr_regions: list, min_len: int = 7,
                                  min_local_score: float = 0.2) -> list:
    """Find local curved DNA motifs (isolated long A/T tracts, not in phased arrays)"""
    tracts = find_polyA_polyT_tracts(seq, min_len)
    results = []
    for start, end, tract_seq, tract_type in tracts:
        s, e = start+1, end+1
        # Exclude if overlaps any global array region
        if not any(r_start <= s <= r_end or r_start <= e <= r_end for r_start, r_end in apr_regions):
            raw, norm = curvature_score(tract_seq)
            if norm >= min_local_score:
                results.append({
                    "Class":"Curved_DNA","Subclass":"Local_Tract",
                    "Start":s,"End":e,"Length":len(tract_seq),
                    "Sequence":wrap(tract_seq),"Score":round(raw,2),
                    "NormScore":round(norm,3),
                    "ScoreMethod":"bend+ATrun","TractType":tract_type
                })
    return results

# --- Main entry: finds all motifs, global and local, normalized scores separate ---
def find_curved_DNA(seq: str, sequence_name: str = "") -> list:
    """Main function to find all curved DNA motifs using Hyperscan for tract finding."""
    global_results, apr_regions = find_global_curved_polyA_polyT(seq)
    local_results = find_local_curved_polyA_polyT(seq, apr_regions)
    all_results = global_results + local_results
    standardized_results = [
        standardize_motif_output(motif, sequence_name, i)
        for i, motif in enumerate(all_results, 1)
    ]
    return standardized_results

# --- Scientific notes ---
# - Uses Hyperscan for high-speed tract finding.
# - Curvature score = sum(dinucleotide bends) + AT-tract (run) bonus + multi-tract bonus (Yella&Bansal, Wang et al).
# - More tracts in phased global array = higher score (reflects greater curvature).
# - Score normalized: 0 (low/no bend) to 1 (high/biologically relevant bend), separate for global/local.
# - Output: both raw and normalized score for each motif, with motif class/subclass and tract count if global.
