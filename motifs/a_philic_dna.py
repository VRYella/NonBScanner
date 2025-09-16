#!/usr/bin/env python3
"""
A-philic DNA motif scanner
- Uses tetranucleotide log2 odds table
- Optionally uses Intel Hyperscan for fast pattern search
- Scores windows >=10bp based on sum of tetranuc log2 odds
- Classifies windows as High-confidence A / Moderate A / not A

This module implements A-philic DNA detection as Class 9 in the NBDFinder system.
A-philic DNA represents DNA sequences with high affinity for specific protein binding
and structural features that favor A-tract formation and protein-DNA interactions.

Scientific References:
- A-tract structural properties: Bolshoy et al. PNAS 1991
- Protein-DNA interactions: Rohs et al. Nature 2009
- Tetranucleotide analysis: Vinogradov Bioinformatics 2003

Author: Dr. Venkata Rajesh Yella
Integration: 2024 - NBDFinder Class 9 Implementation
"""

import math
import csv
import sys
import os
from collections import Counter
from typing import List, Dict, Any, Tuple

# Add path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from motifs.base_motif import standardize_motif_output
except ImportError:
    def standardize_motif_output(motif, sequence_name, motif_id):
        return motif

# --- Tetranuc log2 odds table ---
TET_LOG2 = {
    "AGGG": 3.7004,
    "CCCT": 3.7004,
    "CCCC": 2.5361,
    "GGGG": 2.5361,
    "GCCC": 2.4021,
    "GGGC": 2.4021,
    "CCCA": 2.0704,
    "TGGG": 2.0704,
    "CCTA": 1.585,
    "TAGG": 1.585,
    "ACCC": 0.848,
    "CCCG": 0.848,
    "CGGG": 0.848,
    "GGGT": 0.848,
    "CCGG": 0.2591,
    "GCAC": 0.241,
    "GTGC": 0.241,
    "GGCC": 0.1343,
}

# --- Thresholds ---
STRONG_LOG2 = 2.0
HIGH_MEAN = 2.0
MOD_MEAN = 1.0
FRACTION_HIGH = 5/7
FRACTION_MOD = 4/7

# --- Hyperscan check ---
try:
    import hyperscan
    USE_HS = True
except ImportError:
    USE_HS = False


def build_contrib_array(seq: str, tet_log2: Dict[str, float]) -> List[float]:
    """
    Build contribution array for tetranucleotide log2 odds.
    
    Args:
        seq: DNA sequence (will be converted to uppercase)
        tet_log2: Dictionary of tetranucleotide log2 odds scores
    
    Returns:
        List of contribution scores for each position
    """
    seq = seq.upper()
    n = len(seq)
    contrib = [0.0] * n
    
    for i in range(n - 3):
        tet = seq[i:i+4]
        if tet in tet_log2:
            contrib[i] = tet_log2[tet]
    
    return contrib


def classify_window(sum_log2: float, n_tets: int, strong_count: int) -> str:
    """
    Classify A-philic window based on log2 sum and strong tetranucleotide count.
    
    Args:
        sum_log2: Sum of log2 odds for the window
        n_tets: Number of tetranucleotides in window
        strong_count: Count of strong tetranucleotides (score >= STRONG_LOG2)
    
    Returns:
        Classification string: "A_high_confidence", "A_moderate", or "not_A"
    """
    high_thresh = HIGH_MEAN * n_tets
    mod_thresh = MOD_MEAN * n_tets
    need_high = math.ceil(n_tets * FRACTION_HIGH)
    need_mod = math.ceil(n_tets * FRACTION_MOD)
    
    if sum_log2 >= high_thresh and strong_count >= need_high:
        return "A_high_confidence"
    elif sum_log2 >= mod_thresh and strong_count >= need_mod:
        return "A_moderate"
    return "not_A"


def scan_sequence_internal(seq: str, min_len: int = 10, max_len: int = 15, step: int = 1) -> List[Dict[str, Any]]:
    """
    Internal scanning function that returns A-philic motifs as list of dictionaries.
    
    Args:
        seq: DNA sequence to scan
        min_len: Minimum window length
        max_len: Maximum window length  
        step: Step size for sliding window
    
    Returns:
        List of A-philic motif dictionaries
    """
    seq = seq.upper()
    n = len(seq)
    contrib = build_contrib_array(seq, TET_LOG2)
    
    # Build prefix sums for efficient window scoring
    pref = [0.0] * (n + 1)
    strong_flag = [0] * n
    
    for i in range(n):
        pref[i + 1] = pref[i] + contrib[i]
        if contrib[i] >= STRONG_LOG2:
            strong_flag[i] = 1
    
    pref_strong = [0] * (n + 1)
    for i in range(n):
        pref_strong[i + 1] = pref_strong[i] + strong_flag[i]
    
    # Scan windows and collect A-philic motifs
    motifs = []
    motif_id = 1
    
    for L in range(min_len, max_len + 1):
        n_tets = L - 3
        if n_tets <= 0:
            continue
            
        for s in range(0, n - L + 1, step):
            e = s + L
            sum_log2 = pref[e - 3 + 1] - pref[s]
            strong_c = pref_strong[e - 3 + 1] - pref_strong[s]
            label = classify_window(sum_log2, n_tets, strong_c)
            
            if label != "not_A":
                # Calculate additional metrics
                motif_seq = seq[s:e]
                gc_content = (motif_seq.count('G') + motif_seq.count('C')) / len(motif_seq) if len(motif_seq) > 0 else 0
                
                # Determine subclass based on confidence level
                if label == "A_high_confidence":
                    subclass = "High Confidence A-philic"
                else:
                    subclass = "Moderate A-philic"
                
                motif = {
                    'Class': 'A-philic DNA',
                    'Subclass': subclass,
                    'Start': s + 1,  # 1-based coordinates
                    'End': e,
                    'Length': L,
                    'Sequence': motif_seq,
                    'Score': round(sum_log2, 3),
                    'Actual_Score': round(sum_log2, 3),
                    'Strong_Count': strong_c,
                    'N_Tetranucleotides': n_tets,
                    'GC_Content': round(gc_content, 3),
                    'Classification': label
                }
                
                motifs.append(motif)
                motif_id += 1
    
    return motifs


def find_a_philic_dna(sequence: str, sequence_name: str = "Unknown") -> List[Dict[str, Any]]:
    """
    Main function to find A-philic DNA motifs in a sequence.
    Compatible with NBDFinder motif detection framework.
    
    Args:
        sequence: DNA sequence to analyze
        sequence_name: Name/identifier for the sequence
    
    Returns:
        List of standardized A-philic DNA motif dictionaries
    """
    if not sequence or len(sequence) < 10:
        return []
    
    # Clean sequence
    sequence = sequence.upper().replace('U', 'T')
    
    # Scan for A-philic motifs with extended range for better detection
    motifs = scan_sequence_internal(sequence, min_len=10, max_len=20, step=1)
    
    # Standardize output format for NBDFinder
    standardized_motifs = []
    for i, motif in enumerate(motifs):
        try:
            standardized = standardize_motif_output(motif, sequence_name, i + 1)
            standardized_motifs.append(standardized)
        except Exception:
            # Fallback if standardization fails
            motif['Sequence_Name'] = sequence_name
            motif['Motif_ID'] = f"A_philic_{i + 1}"
            standardized_motifs.append(motif)
    
    return standardized_motifs


def scan_sequence(seq: str, min_len: int = 10, max_len: int = 15, step: int = 1,
                  out_tsv: str = "aphilic_hits.tsv") -> str:
    """
    Legacy function for compatibility with original A-philic scanner interface.
    Scans sequence and writes results to TSV file.
    
    Args:
        seq: DNA sequence to scan
        min_len: Minimum window length
        max_len: Maximum window length
        step: Step size for sliding window
        out_tsv: Output TSV filename
    
    Returns:
        Path to output TSV file
    """
    motifs = scan_sequence_internal(seq, min_len, max_len, step)
    
    # Write to TSV file
    with open(out_tsv, "w", newline="") as fh:
        import csv
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["start", "end", "length", "n_tets", "sum_log2", "strong_count", "label", "seq"])
        
        for motif in motifs:
            w.writerow([
                motif['Start'],
                motif['End'], 
                motif['Length'],
                motif['N_Tetranucleotides'],
                motif['Score'],
                motif['Strong_Count'],
                motif['Classification'],
                motif['Sequence']
            ])
    
    return out_tsv


# Example usage and testing
if __name__ == "__main__":
    example = "NNNNAGGGGGGGGGCCCCTGGGGGCCCAAGGGNNNN"
    motifs = find_a_philic_dna(example, "example_sequence")
    
    print(f"Found {len(motifs)} A-philic DNA motifs:")
    for motif in motifs:
        print(f"  {motif['Subclass']}: {motif['Start']}-{motif['End']} (score: {motif['Score']})")
    
    # Also test legacy TSV output
    out_file = scan_sequence(example, min_len=10, max_len=20, step=1,
                           out_tsv="aphilic_results.tsv")
    print(f"\nResults also written to: {out_file}")