"""
i-Motif Family Motif Detection (Class 7) -- Hyperscan Accelerated
Subclasses: Canonical i-motif (7.1), Relaxed i-motif (7.2), AC-motif (7.3)
"""

import re; import hyperscan  # pip install hyperscan
from motifs.base_motif import wrap, standardize_motif_output
from motifs.hyperscan_manager import optimized_hs_find
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from core.regex_registry import get_patterns_for_motif

# === Load patterns from central registry ===
IMOTIF_PATTERNS = get_patterns_for_motif('i_motif')
from core.regex_registry import get_patterns_for_motif

# --- i-motif scoring using G4Hunter-style algorithm (adapted for C-tracts) ---
def imotif_score(seq):
    """
    G4Hunter-style scoring for i-motifs (Bedrat et al. 2016, adapted for C-richness)
    Algorithm: +1 for C, -1 for G, 0 for A/T; mean score computed.
    For i-motifs, negative scores indicate C-rich regions suitable for i-motif formation.
    Reference: Adapted from Bedrat et al. Nucleic Acids Research 44(4):1746-1759 (2016)
    """
    if not seq or len(seq) == 0:
        return 0.0
    
    # G4Hunter-style scoring: +1 for C, -1 for G, 0 for A/T
    scores = [1 if c == 'C' else -1 if c == 'G' else 0 for c in seq.upper()]
    mean_score = sum(scores) / len(scores) if scores else 0.0
    
    # For i-motifs, we want positive scores (C-rich), so return absolute value
    # Scale by length to give proper weight to longer motifs
    return abs(mean_score) * len(seq)

# --- Hyperscan matcher utility for block-motif finding ---
def hs_find(patterns, seq, subclass_func=None):
    """
    Optimized Hyperscan scanning for i-motif patterns with database caching.
    """
    if not patterns or not seq:
        return []
    
    sequ = seq.upper()
    
    def optimized_callback(id, from_, to, flags, ctx, pattern):
        """Optimized callback for i-motif detection."""
        motif_seq = sequ[from_:to]
        score = pattern[2](motif_seq)
        
        if score > 0:
            c_run_spans = [m.span() for m in re.finditer(r"C{3,}", motif_seq)]
            loops = [c_run_spans[i+1][0]-c_run_spans[i][1] for i in range(len(c_run_spans)-1)] if len(c_run_spans)>1 else []
            
            # Subclass assignment: canonical (all loops 1-7), relaxed (any 8-12), else other
            if subclass_func:
                subclass = subclass_func(loops)
            else:
                subclass = pattern[3]
            
            return {
                "Class": "i-Motif", "Subclass": subclass,
                "Start": from_+1, "End": from_+len(motif_seq),
                "Length": len(motif_seq), "Sequence": wrap(motif_seq),
                "ScoreMethod": "G4Hunter_adapted", "Score": float(score)
            }
        return None
    
    # Use optimized Hyperscan manager
    return optimized_hs_find(patterns, sequ, optimized_callback)

# --- i-motif family finders (all use Hyperscan for primary pattern scan) ---
def find_imotif(seq):
    # Pattern: C3-loop(1-12)-C3-loop(1-12)-C3-loop(1-12)-C3 (per literature, e.g. Zeraati 2018)
    pat = [(pattern[0], pattern[1], pattern[2], pattern[3], imotif_score, None, "G4Hunter_adapted")
           for pattern in IMOTIF_PATTERNS['canonical_imotif']]
    def subclass_func(loops):
        if loops and all(1 <= l <= 7 for l in loops): return "Canonical_iMotif"
        elif loops and any(8 <= l <= 12 for l in loops): return "Relaxed_iMotif"
        else: return "Other_iMotif"
    return hs_find(pat, seq, subclass_func=subclass_func)

def find_ac_motifs(seq):
    # AC-motif: A3-(spacer)-C3-(spacer)-C3-(spacer)-C3 or C3-(spacer)-C3-(spacer)-C3-(spacer)-A3 (per Felsenfeld 1967, Jain 2019)
    # Updated to use G4Hunter-style scoring for consistency
    def ac_score(s):
        # Score based on C-richness and A-tract presence
        scores = [1 if c == 'C' else 0.5 if c == 'A' else -0.5 if c == 'G' else 0 for c in s.upper()]
        return sum(scores)
    
    pat = [(pattern[0], pattern[1], pattern[2], pattern[3], ac_score, pattern[5], pattern[8])
           for pattern in IMOTIF_PATTERNS['ac_motif']]
    return hs_find(pat, seq)

# --- Main entry: all i-motif family (canonical, relaxed, AC) motifs, output standardized ---
def find_i_motif(seq: str, sequence_name: str = "") -> list:
    imotif_results = find_imotif(seq); ac_results = find_ac_motifs(seq)
    all_results = imotif_results + ac_results
    return [standardize_motif_output(motif, sequence_name, i) for i, motif in enumerate(all_results, 1)]

# --- Annotations ---
# - imotif_score: per Zeraati 2018, Jain 2019; C-run compactness, C-fraction, loop bonus.
# - hs_find: utility for block-motif scan, assigns subclass per loop structure.
# - find_imotif: canonical/relaxed/other i-motifs with literature-based loop/scoring, via Hyperscan.
# - find_ac_motifs: AC-motif (A3/C3 alternation), per literature, via Hyperscan.
# - Output: standard, 1-based, for downstream analysis.
