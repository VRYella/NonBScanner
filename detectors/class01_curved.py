"""
Curved DNA Motif Detection (Class 1) - Scoring System (2024 Update)
===================================================================

This module detects poly(A)/poly(T) tracts and phased arrays using Hyperscan
for speed, and scores them for intrinsic DNA curvature. The scoring logic is
aligned with experimental and theoretical literature:

1. Global Arrays (Class 1.1: phased A/T tracts, ~10 bp spacing)
   ------------------------------------------------------------
   - Experimental phasing assays (Marini et al., Cell 1982; Crothers et al.,
     Methods Enzymol 1992) show that short A/T runs of 3–6 bp induce bends
     of ~18° when repeated with helical spacing.
   - DNA curvature increases additively with both tract length (up to ~6 bp,
     where saturation occurs) and the number of phased tracts.
   - Arrays of ≥3 phased tracts are the minimal units of curvature; ≥6 phased
     tracts of 6 bp represent a practical experimental ceiling (Olson et al.,
     PNAS 1998).
   - **Scoring rule**: 
        Raw = Σ tract lengths (clamped to 3–6 bp) + number of phased tracts
        Norm = (Raw – 12) / 30, clipped to [0,1]
        where 12 = baseline (3×3 bp + 3 tracts) and 42 = max (6×6 bp + 6 tracts).

2. Local Tracts (Class 1.2: isolated A/T tracts)
   ----------------------------------------------
   - Isolated A/T runs of ≥7 bp bend DNA and exclude nucleosomes
     (Yella & Bansal, Sci Rep 2017; Wang et al., NAR 2021).
   - The curvature effect strengthens with tract length, up to ~20 bp, beyond
     which nucleosome depletion and curvature effects plateau in vivo.
   - **Scoring rule**:
        Raw = tract length (bp)
        Norm = (L – 7) / 13, clipped to [0,1]
        where L = tract length, 7 bp = minimum effective, ≥20 bp = saturation.

Why this design?
----------------
- Reflects **biophysical models** (bend per tract, wedge/roll) and
  **experimental assays** (electrophoretic mobility, nucleosome positioning).
- Keeps **Raw scores** interpretable (tract length and phased count).
- Provides **Norm scores** (0–1) for comparability across genomes and motifs.
- Balances biological realism with computational simplicity.

References
----------
- Marini et al., Cell 28:871–879 (1982) – phased A-tracts bend DNA
- Crothers et al., Methods Enzymol 212:3–29 (1992) – phasing functions
- Olson et al., PNAS 95:11163 (1998) – wedge/roll model confirmation
- Yella & Bansal, Sci Rep 7:42564 (2017) – nucleosome depletion from long tracts
- Wang et al., NAR 49:e49 (2021) – promoter activity scales with tract length
"""

import numpy as np, re, hyperscan
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

# Import from motifs directory where these modules actually exist
from motifs.base_motif import wrap, standardize_motif_output
from motifs.hyperscan_manager import optimized_hs_find
from core.regex_registry import get_patterns_for_motif

# === Load patterns from central registry ===
CURVED_PATTERNS = get_patterns_for_motif('curved_dna')

# --------------------------
# Global scoring constants
# --------------------------
MIN_GLOBAL_TRACT = 3          # bp: minimum tract length to count in phased arrays
CLAMP_GLOBAL_MIN = 3          # bp clamp lower bound
CLAMP_GLOBAL_MAX = 6          # bp clamp upper bound (curvature plateaus)
BASELINE_TRACTS = 3           # theoretical minimum phased tracts
MAX_PHASED_FOR_NORM = 6       # theoretical maximum phased tracts for normalization

# Derived theoretical anchors for normalization
_BASELINE_SUMLEN = BASELINE_TRACTS * CLAMP_GLOBAL_MIN     # 3 tracts × 3 bp = 9
_BASELINE_RAW     = _BASELINE_SUMLEN + BASELINE_TRACTS    # 9 + 3 = 12
_MAX_SUMLEN       = MAX_PHASED_FOR_NORM * CLAMP_GLOBAL_MAX  # 6 × 6 = 36
_MAX_RAW          = _MAX_SUMLEN + MAX_PHASED_FOR_NORM       # 36 + 6 = 42
_NORM_SPAN        = max(1.0, _MAX_RAW - _BASELINE_RAW)      # 30

# -------------------------
# Local scoring constants
# -------------------------
MIN_LOCAL_LEN = 7             # bp; minimum length qualifying as local tract
LOCAL_LEN_CAP = 20            # bp; normalization saturates at ≥20 bp


# --- Parameter Table (unchanged) ---
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


# --- UPDATED GLOBAL SCORER (length + count) ---
def curvature_score(seq: str, tracts=None):
    """
    Returns (raw_score, normalized_score).

    GLOBAL arrays (when `tracts` is provided):
      - For each phased tract, compute Lc = clamp(len, 3..6).
      - Raw score = SUM(Lc) + (# of phased tracts).
      - Norm score = clip((raw - 12) / 30, 0, 1)  where:
          baseline raw = 3 tracts of length 3 → 12
          max raw      = 6 tracts of length 6 → 42

    LOCAL tracts (when `tracts` is None):
      - Raw score = tract length (len(seq)).
      - Norm score = clip((L - 7) / (20 - 7), 0, 1); saturates at L ≥ 20.

    NOTE: We keep this single function signature for API stability.
    """
    if tracts is None or len(tracts) == 0:
        # LOCAL scoring path
        L = len(seq)
        raw = float(L)
        if L <= MIN_LOCAL_LEN:
            norm = 0.0
        else:
            norm = min(1.0, max(0.0, (L - MIN_LOCAL_LEN) / float(max(1, LOCAL_LEN_CAP - MIN_LOCAL_LEN))))
        return raw, norm

    # GLOBAL scoring path
    # Only consider tracts with length >= 3 bp; clamp each to [3,6] for length contribution
    Lsum = 0.0
    count = 0
    for (_s, _e, tseq, _ttype) in tracts:
        L = len(tseq)
        if L >= MIN_GLOBAL_TRACT:
            Lc = max(CLAMP_GLOBAL_MIN, min(CLAMP_GLOBAL_MAX, L))
            Lsum += Lc
            count += 1

    raw = Lsum + count  # length + number of phased tracts
    norm = min(1.0, max(0.0, (raw - _BASELINE_RAW) / _NORM_SPAN))
    return float(raw), float(norm)


# --- Global (Phased Array) finder: phased tracts (~10bp period, 3+ repeats) ---
def find_global_curved_polyA_polyT(seq: str, min_tract_len: int = 3, min_repeats: int = 3,
                                   min_spacing: int = 8, max_spacing: int = 12,
                                   min_global_score: float = 0.2) -> tuple:
    """Find global curved DNA motifs (phased A/T tracts, spaced ~10bp, 3+ in array)."""
    tracts = find_polyA_polyT_tracts(seq, min_tract_len)
    results, apr_regions = [], []
    for i in range(len(tracts)-min_repeats+1):
        group = [tracts[i]]
        # seed with min_repeats phased
        for j in range(1, min_repeats):
            prev_center = (tracts[i+j-1][0]+tracts[i+j-1][1])//2
            curr_center = (tracts[i+j][0]+tracts[i+j][1])//2
            spacing = curr_center-prev_center
            if min_spacing <= spacing <= max_spacing:
                group.append(tracts[i+j])
            else:
                break
        # extend group if more phased tracts found
        k = i+len(group)
        while k < len(tracts):
            prev_center = (tracts[k-1][0]+tracts[k-1][1])//2
            curr_center = (tracts[k][0]+tracts[k][1])//2
            spacing = curr_center-prev_center
            if min_spacing <= spacing <= max_spacing:
                group.append(tracts[k]); k += 1
            else:
                break
        if len(group) >= min_repeats:
            motif_seq = seq[group[0][0]:group[-1][1]+1]
            raw, norm = curvature_score(motif_seq, group)  # length+count global scoring
            if norm >= min_global_score:
                motif = {
                    "Class":"Curved_DNA","Subclass":"Global_Array",
                    "Start":group[0][0]+1,"End":group[-1][1]+1,
                    "Length":group[-1][1]-group[0][0]+1,
                    "Sequence":wrap(motif_seq),
                    "Score":round(raw,3),"NormScore":round(norm,3),
                    "ScoreMethod":"phasedCount+tractLen(3..6)",
                    "NumTracts":len(group)
                }
                results.append(motif)
                apr_regions.append((motif["Start"], motif["End"]))
    return results, apr_regions


# --- Local tract finder: isolated long tracts (not in global arrays) ---
def find_local_curved_polyA_polyT(seq: str, apr_regions: list, min_len: int = MIN_LOCAL_LEN,
                                  min_local_score: float = 0.2) -> list:
    """
    Find local curved DNA motifs (isolated long A/T tracts, not in phased arrays).

    UPDATED local scoring:
      - raw = tract length L
      - norm = (L - 7) / (20 - 7), clipped to [0,1], saturating for L ≥ 20
      - applies to both A- and T-tracts (poly(dA:dT) symmetry)
    """
    tracts = find_polyA_polyT_tracts(seq, min_len)

    # collect local (non-overlapping with global arrays)
    local_candidates = []
    for start, end, tract_seq, tract_type in tracts:
        s, e = start+1, end+1
        if any(r_start <= s <= r_end or r_start <= e <= r_end for r_start, r_end in apr_regions):
            continue
        L = len(tract_seq)
        if L >= MIN_LOCAL_LEN:
            local_candidates.append((s, e, tract_seq, tract_type, L))

    results = []
    for s, e, tract_seq, tract_type, L in local_candidates:
        raw = float(L)
        if L <= MIN_LOCAL_LEN:
            norm = 0.0
        else:
            norm = min(1.0, max(0.0, (L - MIN_LOCAL_LEN) / float(max(1, LOCAL_LEN_CAP - MIN_LOCAL_LEN))))
        if norm >= min_local_score:
            results.append({
                "Class":"Curved_DNA","Subclass":"Local_Tract",
                "Start":s,"End":e,"Length":L,
                "Sequence":wrap(tract_seq),
                "Score":round(raw,3),
                "NormScore":round(norm,3),
                "ScoreMethod":"tractLen(>=7)_min7_max20",
                "TractType":tract_type
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
