#!/usr/bin/env python3
"""
scan_aphilic_hyperscan.py

Fast scanning of sequences for A-philic motifs using Intel Hyperscan (if available).
Produces TSV annotations for windows classified as High / Moderate A-form.

Behavior:
 - Build patterns for all tetranucleotides in table and scan the sequence to get match starts.
 - Sum per-position log2 contributions across windows of requested lengths.
 - Classify windows by sum thresholds (sum-based score reflects motif length).

Author: ChatGPT (Jimmy)
"""

from collections import defaultdict, Counter
import math
import csv
import sys
import shutil

# ----------------------------
# Tetranuc table (user-provided)
# tet -> log2_odds
# ----------------------------
TET_LOG2 = {
    "AGGG": 3.7004,
    "CCCT": 3.7004,
    "CCCC": 2.5361,
    "GGGG": 2.5361,
    "GCCC": 2.4021,
    "GGGC": 2.4021,
    "CCCA": 2.0704,
    "TGGG": 2.0704,
    "CCTA": 1.5850,
    "TAGG": 1.5850,
    "ACCC": 0.8480,
    "CCCG": 0.8480,
    "CGGG": 0.8480,
    "GGGT": 0.8480,
    "CCGG": 0.2591,
    "GCAC": 0.2410,
    "GTGC": 0.2410,
    "GGCC": 0.1343,
}

# Parameters (tweak if desired)
STRONG_LOG2_CUTOFF = 2.0  # tetranuc considered "strong" if log2 >= this
HIGH_MEAN_LOG2 = 2.0      # high-confidence mean threshold
MODERATE_MEAN_LOG2 = 1.0  # moderate mean threshold
# We'll convert mean thresholds to sum thresholds: sum_thresh = mean * n_tets

# Fraction-of-strong-tets required
STRONG_FRACTION_HIGH = 5.0 / 7.0   # for high-confidence
STRONG_FRACTION_MOD = 4.0 / 7.0    # for moderate

# ----------------------------
# Try to import hyperscan; otherwise set fallback
# ----------------------------
USE_HYPERSCAN = False
try:
    import hyperscan
    USE_HYPERSCAN = True
except Exception:
    USE_HYPERSCAN = False

# ----------------------------
# Helper functions
# ----------------------------
def revcomp(s: str) -> str:
    tr = str.maketrans("ACGT", "TGCA")
    return s.translate(tr)[::-1]

def build_match_array_hyperscan(seq: str, tet_log2: dict) -> list:
    """
    Use hyperscan to find all 4-mer matches. Return an array 'contrib' of length len(seq),
    where contrib[i] is the sum of log2 values for tetranucleotides that start at position i.
    """
    try:
        db = hyperscan.Database()
        expressions = []
        ids = []
        id_to_log2 = {}
        
        # Create patterns with proper boundaries to match exactly 4 characters
        for idx, (tet, log2) in enumerate(tet_log2.items()):
            # Use word boundaries or lookahead to match exactly 4 chars
            pattern = f"{tet}(?=[ATGCN]|$)"  # match tet followed by DNA base or end
            expressions.append(pattern.encode())
            ids.append(idx)
            id_to_log2[idx] = log2
        
        # Compile with appropriate flags
        db.compile(expressions=expressions, ids=ids, elements=len(expressions))
        
        # Prepare contributions array
        n = len(seq)
        contrib = [0.0] * n
        
        def on_match(id, start, end, flags, context):
            # Only add contribution if this is exactly a 4-mer match
            if end - start == 4 and start < len(contrib):
                contrib[start] += id_to_log2[id]
            return hyperscan.HS_SUCCESS
        
        # Scan the sequence
        db.scan(seq.encode(), match_event_handler=on_match)
        return contrib
        
    except Exception:
        # Fallback to Python implementation if hyperscan fails
        return build_match_array_py(seq, tet_log2)

def build_match_array_py(seq: str, tet_log2: dict) -> list:
    """Pure-python: for each tet, find all occurrences (overlapping) and add log2 at start."""
    n = len(seq)
    contrib = [0.0] * n
    for tet, log2 in tet_log2.items():
        start = seq.find(tet)
        pos = -1
        if start == -1:
            # try sliding via manual search to capture overlapping
            # naive scan:
            for i in range(0, n - 4 + 1):
                if seq[i:i+4] == tet:
                    contrib[i] += log2
            continue
        # else there are occurrences — still better to do naive loop to catch overlaps
        for i in range(0, n - 4 + 1):
            if seq[i:i+4] == tet:
                contrib[i] += log2
    return contrib

def build_contrib_array(seq: str, tet_log2: dict):
    """
    Build array of contributions per 4-mer start position.
    At position i we put sum of log2 values for 4-mers starting at i (most positions will
    be zero if the 4-mer isn't in table).
    """
    seq = seq.upper()
    if USE_HYPERSCAN:
        try:
            return build_match_array_hyperscan(seq, tet_log2)
        except Exception:
            # fall back to python
            return build_match_array_py(seq, tet_log2)
    else:
        return build_match_array_py(seq, tet_log2)

def classify_window(sum_log2: float, n_tets: int, strong_count: int) -> str:
    """
    Apply sum-based thresholds to classify window.
    """
    high_sum_thresh = HIGH_MEAN_LOG2 * n_tets
    mod_sum_thresh = MODERATE_MEAN_LOG2 * n_tets
    strong_needed_high = math.ceil(n_tets * STRONG_FRACTION_HIGH)
    strong_needed_mod = math.ceil(n_tets * STRONG_FRACTION_MOD)
    if (sum_log2 >= high_sum_thresh) and (strong_count >= strong_needed_high):
        return "A_high_confidence"
    elif (sum_log2 >= mod_sum_thresh) and (strong_count >= strong_needed_mod):
        return "A_moderate"
    else:
        return "not_A"

# ----------------------------
# Main scanning function
# ----------------------------
def scan_sequence(seq: str,
                  tet_log2: dict,
                  min_window: int = 10,
                  max_window: int = 10,
                  step: int = 1,
                  out_tsv: str = "aphilic_windows.tsv"):
    """
    Scan seq for windows of lengths between min_window and max_window (inclusive).
    step: slide start position by this many bases (default 1 for full scan).
    Writes annotated TSV with columns:
      seq_id, start, end, window_len, n_tets, sum_log2, strong_count, label, window_seq
    """
    seq = seq.upper()
    n = len(seq)
    contrib = build_contrib_array(seq, tet_log2)  # contributions at 4-mer starts

    # Build prefix-sum array of contributions at positions 0..n-1
    # But remember only starts 0..n-4 are relevant. We'll compute prefix sums over contrib.
    pref = [0.0] * (n + 1)
    for i in range(n):
        pref[i+1] = pref[i] + contrib[i]

    # Also build an array marking whether a 4-mer at pos i is "strong"
    strong_flag = [0] * n
    for i in range(0, n - 4 + 1):
        tet = seq[i:i+4]
        v = tet_log2.get(tet, 0.0)
        if v >= STRONG_LOG2_CUTOFF:
            strong_flag[i] = 1
    pref_strong = [0] * (n + 1)
    for i in range(n):
        pref_strong[i+1] = pref_strong[i] + strong_flag[i]

    # Open TSV and write header
    with open(out_tsv, "w", newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(["seq_id","start","end","window_len","n_tets","sum_log2","strong_count","label","window_seq"])
        seq_id = "seq1"
        for L in range(min_window, max_window+1):
            n_tets = L - 3
            if n_tets <= 0:
                continue
            for start in range(0, n - L + 1, step):
                end = start + L
                # tet starts are start..end-4 inclusive
                tet_start = start
                tet_end = end - 4  # inclusive index of last 4-mer start
                # sum contributions from contrib[tet_start] .. contrib[tet_end]
                sum_log2 = pref[tet_end+1] - pref[tet_start]
                strong_count = pref_strong[tet_end+1] - pref_strong[tet_start]
                label = classify_window(sum_log2, n_tets, strong_count)
                if label != "not_A":
                    window_seq = seq[start:end]
                    w.writerow([seq_id, start, end, L, n_tets, round(sum_log2,6), strong_count, label, window_seq])
    return out_tsv

# ----------------------------
# Example usage (demo)
# ----------------------------
if __name__ == "__main__":
    # Example sequence — replace with your genome sequence (FASTA parsing optional)
    example_seq = (
        "NNNNN" +
        "AAGTCCAGGGGGGGGGGTTTACG" +
        "TTTAGGGGGGGGGGCCCAAATTT" +
        "ACGTACGTACGT" +
        "TGGGGGGGGGGA" +
        "CCCCGCCCAGGGTGCCC" +
        "GGGATGGGGTGGG" +
        "NNNNN"
    ).upper()

    # provide window range; user asked "10 mer sequence or more" — here we scan 10..20
    OUT = "aphilic_hyperscan_results.tsv"
    out = scan_sequence(example_seq, TET_LOG2, min_window=10, max_window=20, step=1, out_tsv=OUT)
    print("Wrote annotated A-philic windows to:", out)
    print("Hyperscan used?" , USE_HYPERSCAN)