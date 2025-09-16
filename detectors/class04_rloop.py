"""
R-Loop Motif Detection (Class 4)
Subclasses: R-loop (4.1)
"""

import re
from .base_motif import gc_content, wrap, standardize_motif_output
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from regex_registry import get_patterns_for_motif

# === Load patterns from central registry ===
RLOOP_PATTERNS = get_patterns_for_motif('r_loop')


# RLFS models for R-loop prediction - using patterns from registry
RLFS_MODELS = {
    "m1": RLOOP_PATTERNS['rlfs_m1'][0][0] if RLOOP_PATTERNS.get('rlfs_m1') else r"G{3,}[ATGC]{1,10}?G{3,}(?:[ATGC]{1,10}?G{3,}){1,}",
    "m2": RLOOP_PATTERNS['rlfs_m2'][0][0] if RLOOP_PATTERNS.get('rlfs_m2') else r"G{4,}(?:[ATGC]{1,10}?G{4,}){1,}",
}


def find_rez_max(seq, start_pos, max_len=2000, step=100, min_gc=40):
    """Find the maximum GC-rich window for R-loop extension zone"""
    max_window = ""
    for win_start in range(start_pos, min(len(seq), start_pos + max_len), step):
        win_end = min(win_start + step, len(seq))
        window = seq[win_start:win_end]
        if gc_content(window) >= min_gc and len(window) > len(max_window):
            max_window = window
    if max_window:
        return {'seq': max_window, 'end': len(max_window)}
    return None


def find_rlfs(seq: str, models=("m1", "m2"), sequence_name: str = "") -> list:
    """Find R-loop forming sequences"""
    if len(seq) < 100:
        return []
    
    results = []
    for model_name in models:
        pattern = RLFS_MODELS[model_name]
        for m in re.finditer(pattern, seq, re.IGNORECASE):
            riz_seq = m.group(0)
            if gc_content(riz_seq) < 50:
                continue
            rez = find_rez_max(seq, m.end())
            if rez:
                rez_seq = rez['seq']
                concat = riz_seq + rez_seq
                g_runs = len(re.findall(r"G{3,}", concat))
                # Raw stability: GC fraction weight + G-run density scaled by length
                gc_frac = gc_content(concat) / 100.0
                score = (gc_frac * 50.0 + g_runs * 10.0) * (len(concat) ** 0.25)
                results.append({
                    "Class": "R-Loop",
                    "Subclass": f"RLFS_{model_name}",
                    "Start": m.start() + 1,
                    "End": m.start() + len(riz_seq) + rez['end'],
                    "Length": len(riz_seq) + rez['end'],
                    "Sequence": wrap(concat),
                    "ScoreMethod": "QmRLFS_raw",
                    "Score": float(score),
                })
    
    # Standardize output format
    standardized_results = []
    for i, motif in enumerate(results, 1):
        standardized_results.append(standardize_motif_output(motif, sequence_name, i))
    
    return standardized_results


def find_r_loop(seq: str, sequence_name: str = "") -> list:
    """Main function to find R-loop motifs"""
    return find_rlfs(seq, sequence_name=sequence_name)