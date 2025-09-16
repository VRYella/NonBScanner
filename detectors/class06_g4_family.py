"""
G-Quadruplex Family Motif Detection (Class 6) accelerated with Hyperscan.

SCIENTIFIC BASIS:
================
G-quadruplexes are four-stranded DNA structures formed by Hoogsteen base pairing of guanine tetrads.
They are stabilized by monovalent cations (K+, Na+) and play crucial roles in:
- Telomere biology and replication timing (Blackburn & Collins, 2011)
- Gene regulation and transcriptional control (Huppert & Balasubramanian, 2007)
- Genome instability and cancer biology (Maizels, 2015)

SUBCLASSES DETECTED:
===================
1. Canonical G4: G3+N1-7G3+N1-7G3+N1-7G3+ (classic definition, Williamson 2005)
2. Bulged G4: G3+N1-7G2+N1-3G1+N1-7G3+N1-7G3+ (bulges tolerated, Todd et al. 2005)
3. Relaxed G4: G2+N1-12G2+N1-12G2+N1-12G2+ (relaxed criteria, Kikin et al. 2006)
4. Bipartite G4: Two G4-forming regions connected by long spacer (>30bp)
5. Multimeric G4: Four or more G-tracts forming complex structures
6. Imperfect G4: Non-consecutive G-runs with interruptions

HYPERSCAN ACCELERATION:
======================
Uses Intel Hyperscan for primary pattern matching of G-tract arrangements,
followed by G4Hunter scoring (Bedrat et al. 2016) for biological relevance filtering.

OUTPUT FORMAT: 1-based coordinates suitable for genomic analysis pipelines.
"""

import re; import numpy as np; import hyperscan  # pip install hyperscan
from .base_motif import overlapping_finditer, wrap, standardize_motif_output
from .hyperscan_manager import optimized_hs_find
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from regex_registry import get_patterns_for_motif

# === G4Hunter Scoring Algorithm (Bedrat et al. 2016) ===
def g4hunter_score(seq):
    """
    Calculate G4Hunter score for G-quadruplex prediction.
    
    SCIENTIFIC BASIS: G4Hunter algorithm weights G-richness vs C-richness
    to predict G4-forming potential. Score >1.2 indicates high G4 potential.
    
    Algorithm: +1 for G, -1 for C, 0 for A/T; mean score computed.
    Reference: Bedrat et al. Nucleic Acids Research 44(4):1746-1759 (2016)
    """
    scores = [1 if c == 'G' else -1 if c == 'C' else 0 for c in seq.upper()]
    return np.mean(scores) if scores else 0.0

# === Hyperscan Accelerated Pattern Matching for G4 Detection ===
def hs_find(patterns, seq, group=0):
    """
    High-performance pattern matching using Intel Hyperscan for G4 motif detection.
    
    TECHNICAL IMPLEMENTATION:
    - Uses optimized Hyperscan database manager with caching
    - Performs parallel pattern matching on sequence
    - Applies scientific filters (G4Hunter score, G-run count) during callback
    
    PARAMETERS:
    patterns: List of (regex, id, groupno, subclass, score_func, score_scale, min_g_runs, min_g4hunter)
    - regex: Pattern for G-tract arrangement (e.g., G{3,}N{1,7}G{3,}...)
    - score_func: Scientific scoring function (G4Hunter, custom)
    - score_scale: Subclass-specific scaling factor for biological relevance
    - min_g_runs: Minimum number of G3+ tracts required
    - min_g4hunter: Minimum G4Hunter score threshold for inclusion
    
    Returns: List of validated G4 motifs with 1-based coordinates
    """
    if not patterns or not seq:
        return []
    
    sequ = seq.upper()
    
    def optimized_callback(id, from_, to, flags, ctx, pattern):
        """Optimized callback function with pattern access."""
        matched_seq = sequ[from_:to]
        # Use re.search instead of re.match to find pattern anywhere in the matched region
        try:
            m = re.search(pattern[0], matched_seq)
        except:
            return None
        
        if not m: 
            return None
        
        # Get the motif sequence (full match for group 0)
        motif_seq = m.group(0) if pattern[2] == 0 else m.group(pattern[2])
        
        # Calculate actual coordinates within the original sequence
        actual_start = from_ + m.start()
        actual_end = from_ + m.end()
        
        # Calculate score with scaling
        base_score = pattern[4](motif_seq)
        if pattern[5] != 1.0:  # Apply score scaling if specified
            score = base_score * len(motif_seq) * pattern[5]
        else:
            score = base_score * len(motif_seq)
        
        # Apply biological filters
        g_runs = len(re.findall(r"G{3,}", motif_seq))
        if pattern[6] and g_runs < pattern[6]: 
            return None
        if pattern[7] and base_score < pattern[7]: 
            return None
            
        return {
            "Class": "G-Quadruplex", "Subclass": pattern[3], "Start": actual_start+1, "End": actual_end,
            "Length": len(motif_seq), "Sequence": wrap(motif_seq),
            "ScoreMethod": pattern[8], "Score": float(score)
        }
    
    # Use optimized Hyperscan manager
    return optimized_hs_find(patterns, sequ, optimized_callback)

# === Load patterns from central registry ===
G4_PATTERNS = get_patterns_for_motif('g_quadruplex')

# --- All motif finders below use Hyperscan for primary regex matching ---

def find_multimeric_gquadruplex(seq):
    # Multimeric: four or more G3-tracts with up to 12bp loops, G4Hunter >= 0.3
    pat = [(pattern[0], pattern[1], pattern[2], pattern[3], g4hunter_score, pattern[5], pattern[6], pattern[7], pattern[8])
           for pattern in G4_PATTERNS['multimeric_g4']]
    return hs_find(pat, seq)

def find_bipartite_gquadruplex(seq):
    # Bipartite: 8 G3 tracts, special internal loop, at least 8 G3 runs, max score of halves
    pat = [(pattern[0], pattern[1], pattern[2], pattern[3], 
            lambda s: max(g4hunter_score(s[:len(s)//2]), g4hunter_score(s[len(s)//2:])), 
            pattern[5], pattern[6], pattern[7], pattern[8])
           for pattern in G4_PATTERNS['bipartite_g4']]
    return hs_find(pat, seq)

def find_gquadruplex(seq):
    # Canonical: 4 G3 tracts, loops 1–7, G4Hunter >= 0.5 (adjusted for better sensitivity)
    pat = [(pattern[0], pattern[1], pattern[2], pattern[3], g4hunter_score, pattern[5], pattern[6], pattern[7], pattern[8])
           for pattern in G4_PATTERNS['canonical_g4']]
    return hs_find(pat, seq)

def find_relaxed_gquadruplex(seq):
    # Relaxed: as canonical but longer loops 8–12, G4Hunter >=0.3 (more permissive)
    pat = [(pattern[0], pattern[1], pattern[2], pattern[3], g4hunter_score, pattern[5], pattern[6], pattern[7], pattern[8])
           for pattern in G4_PATTERNS['relaxed_g4']]
    return hs_find(pat, seq)

def find_bulged_gquadruplex(seq):
    # Bulged: up to 3nt loops, at least 4 G3 runs, score scale 0.7
    pat = [(pattern[0], pattern[1], pattern[2], pattern[3], g4hunter_score, pattern[5], pattern[6], pattern[7], pattern[8])
           for pattern in G4_PATTERNS['bulged_g4']]
    return hs_find(pat, seq)

def find_imperfect_gquadruplex(seq):
    # Imperfect: one G2 tract, the rest G3, G4Hunter >=0.4 (more permissive)
    pat = [(pattern[0], pattern[1], pattern[2], pattern[3], g4hunter_score, pattern[5], pattern[6], pattern[7], pattern[8])
           for pattern in G4_PATTERNS['imperfect_g4']]
    res = []; [res.extend(hs_find([p], seq)) for p in pat]
    return res

def find_gtriplex(seq):
    # G-triplex: three G3-tracts, loops 1–7, score from G-runs/loop
    def g_triplex_score(s):
        return (sum(len(r) for r in re.findall(r"G{3,}", s))*2.0) + \
               (sum(1/l if l>0 else 0.5 for l in [len(x) for x in re.findall(r"G{3,}(\w{1,7})G{3,}", s)])*5.0)
    
    pat = [(pattern[0], pattern[1], pattern[2], pattern[3], g_triplex_score, pattern[5], pattern[6], pattern[7], pattern[8])
           for pattern in G4_PATTERNS['g_triplex']]
    return hs_find(pat, seq)

# --- Master function: finds all G4-family motifs and standardizes output ---
def find_g_quadruplex(seq: str, sequence_name: str = "") -> list:
    results = []; results.extend(find_multimeric_gquadruplex(seq)); results.extend(find_bipartite_gquadruplex(seq))
    results.extend(find_gquadruplex(seq)); results.extend(find_relaxed_gquadruplex(seq)); results.extend(find_bulged_gquadruplex(seq))
    results.extend(find_imperfect_gquadruplex(seq)); results.extend(find_gtriplex(seq))
    # Standardize output as per NBDFinder conventions; 1-based coordinates
    return [standardize_motif_output(motif, sequence_name, i) for i, motif in enumerate(results, 1)]

# --- Annotations ---
# - Each motif finder is mapped to scientific G4 family definitions (see literature: G4Hunter, Bedrat 2016; Hon 2017 NAR; Kwok 2016).
# - Hyperscan is used for the initial regex scan for maximal performance on large genomes.
# - Each class uses a scoring/thresholding system in line with the literature (G4Hunter, triplex scoring, etc).
# - Output is standardized, with coordinates 1-based, for downstream interoperability.
