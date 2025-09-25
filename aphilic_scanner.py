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

# Import from the main NBDFinder system
try:
    from detection_scoring import detect_aphilic_motifs, aphilic_score
except ImportError:
    # Fallback implementations
    def detect_aphilic_motifs(sequence, sequence_name):
        return scan_sequence_internal(sequence)
    
    def aphilic_score(sequence, start, end):
        return calculate_aphilic_score(sequence[start:end])

try:
    from motifs.base_motif import standardize_motif_output
except ImportError:
    def standardize_motif_output(motif, sequence_name, motif_id):
        return motif

# === A-philic 10-mer log2 odds table ===
# | Pattern        | Type   | Log2 Odds Score                    | Range        | Description                          |
# |----------------|--------|------------------------------------|--------------|--------------------------------------|
# | ACCCCCCCCA     | float  | Tetranucleotide log2 odds          | 1.0 - 3.0    | High A-philic potential             |
# | CCCCCCCCCC     | float  | Maximum homopolymer score          | 2.5 - 2.7    | Strongest A-tract affinity          |
# | Mixed patterns | float  | Variable based on composition      | 1.0 - 2.5    | Moderate to strong A-philic         |

APHILIC_NMER_LOG2 = {
    "ACCCCCCCCA": 2.2284142857142855,
    "ACCCCCCCCC": 2.294942857142857,
    "ACCCCCCCCG": 2.053785714285714,
    "ACCCCCCCCT": 2.4612714285714286,
    "ACCCCCCCGG": 1.7285,
    "ACCCCCCCTA": 2.3253999999999997,
    "ACCCCCCGGG": 1.4873428571428569,
    "ACCCCCGGGC": 1.4682,
    "ACCCCCGGGG": 1.4873428571428569,
    "ACCCCCGGGT": 1.2461857142857142,
    "ACCCCGGGCC": 1.1250857142857142,
    "ACCCCGGGGC": 1.4682,
    "ACCCCGGGGG": 1.4873428571428569,
    "ACCCCGGGGT": 1.2461857142857142,
    "ACCCGGGCCC": 1.1059428571428571,
    "ACCCGGGGCC": 1.1250857142857142,
    "ACCCGGGGGC": 1.4682,
    "ACCCGGGGGG": 1.4873428571428569,
    "ACCCGGGGGT": 1.2461857142857142,
    "AGGGCCCCCA": 2.2544999999999997,
    "AGGGCCCCCC": 2.321028571428571,
    "AGGGCCCCCG": 2.0798714285714284,
    "AGGGCCCCCT": 2.487357142857143,
    "AGGGCCCCGG": 1.7545857142857142,
    "AGGGCCCCTA": 2.3514857142857144,
    "AGGGCCCGGG": 1.5134285714285713,
    "AGGGGCCCCA": 2.2544999999999997,
    "AGGGGCCCCC": 2.321028571428571,
    "AGGGGCCCCG": 2.0798714285714284,
    "AGGGGCCCCT": 2.487357142857143,
    "AGGGGCCCGG": 1.7545857142857142,
    "AGGGGCCCTA": 2.3514857142857144,
    "AGGGGGCCCA": 2.2544999999999997,
    "AGGGGGCCCC": 2.321028571428571,
    "AGGGGGCCCG": 2.0798714285714284,
    "AGGGGGCCCT": 2.487357142857143,
    "AGGGGGGCCC": 2.321028571428571,
    "AGGGGGGGCC": 2.3401714285714283,
    "AGGGGGGGGC": 2.683285714285714,
    "AGGGGGGGGG": 2.702428571428571,
    "AGGGGGGGGT": 2.4612714285714286,
    "CCCCCCCCCA": 2.4695714285714283,
    "CCCCCCCCCC": 2.5361,
    "CCCCCCCCCG": 2.294942857142857,
    "CCCCCCCCCT": 2.702428571428571,
    "CCCCCCCCGG": 1.9696571428571428,
    "CCCCCCCCTA": 2.5665571428571425,
    "CCCCCCCGGG": 1.7285,
    "CCCCCCGGGC": 1.7093571428571426,
    "CCCCCCGGGG": 1.7285,
    "CCCCCCGGGT": 1.4873428571428569,
    "CCCCCGGGCC": 1.366242857142857,
    "CCCCCGGGGC": 1.7093571428571426,
    "CCCCCGGGGG": 1.7285,
    "CCCCCGGGGT": 1.4873428571428569,
    "CCCCGGGCCC": 1.3471,
    "CCCCGGGGCC": 1.366242857142857,
    "CCCCGGGGGC": 1.7093571428571426,
    "CCCCGGGGGG": 1.7285,
    "CCCCGGGGGT": 1.4873428571428569,
    "CGGGCCCCCA": 1.8470142857142857,
    "CGGGCCCCCC": 1.9135428571428572,
    "CGGGCCCCCG": 1.6723857142857141,
    "CGGGCCCCCT": 2.0798714285714284,
    "CGGGCCCCGG": 1.3471,
    "CGGGCCCCTA": 1.9440000000000002,
    "CGGGCCCGGG": 1.1059428571428571,
    "CGGGGCCCCA": 1.8470142857142857,
    "CGGGGCCCCC": 1.9135428571428572,
    "CGGGGCCCCG": 1.6723857142857141,
    "CGGGGCCCCT": 2.0798714285714284,
    "CGGGGCCCGG": 1.3471,
    "CGGGGCCCTA": 1.9440000000000002,
    "CGGGGGCCCA": 1.8470142857142857,
    "CGGGGGCCCC": 1.9135428571428572,
    "CGGGGGCCCG": 1.6723857142857141,
    "CGGGGGCCCT": 2.0798714285714284,
    "CGGGGGGCCC": 1.9135428571428572,
    "CGGGGGGGCC": 1.9326857142857141,
    "CGGGGGGGGC": 2.2758,
    "CGGGGGGGGG": 2.294942857142857,
    "CGGGGGGGGT": 2.053785714285714,
    "GCCCCCCCCA": 2.450428571428571,
    "GCCCCCCCCC": 2.5169571428571422,
    "GCCCCCCCCG": 2.2758,
    "GCCCCCCCGG": 1.9505142857142856,
    "GCCCCCCCTA": 2.5474142857142854,
    "GCCCCCCGGG": 1.7093571428571426,
    "GCCCCCGGGC": 1.6902142857142857,
    "GCCCCCGGGG": 1.7093571428571426,
    "GCCCCCGGGT": 1.4682,
    "GCCCCGGGCC": 1.366242857142857,
    "GCCCCGGGGC": 1.6902142857142857,
    "GCCCCGGGGG": 1.7093571428571426,
    "GCCCCGGGGT": 1.4682,
    "GCCCGGGCCC": 1.3279571428571428,
    "GCCCGGGGCC": 1.366242857142857,
    "GCCCGGGGGC": 1.6902142857142857,
    "GCCCGGGGGG": 1.7093571428571426,
    "GCCCGGGGGT": 1.4682,
    "GGCCCCCCCA": 2.1073142857142857,
    "GGCCCCCCCC": 2.173842857142857,
    "GGCCCCCCCG": 1.9326857142857141,
    "GGCCCCCCCT": 2.3401714285714283,
    "GGCCCCCCGG": 1.6074,
    "GGCCCCCCTA": 2.2043,
    "GGCCCCCGGG": 1.366242857142857,
    "GGCCCCGGGG": 1.366242857142857,
    "GGCCCCGGGT": 1.1250857142857142,
    "GGCCCGGGCC": 1.0039857142857143,
    "GGCCCGGGGG": 1.366242857142857,
    "GGCCCGGGGT": 1.1250857142857142,
    "GGGCCCCCCA": 2.0881714285714286,
    "GGGCCCCCCC": 2.1546999999999996,
    "GGGCCCCCCG": 1.9135428571428572,
    "GGGCCCCCCT": 2.321028571428571,
    "GGGCCCCCGG": 1.5882571428571428,
    "GGGCCCCCTA": 2.185157142857143,
    "GGGCCCCGGG": 1.3471,
    "GGGCCCGGGC": 1.3279571428571428,
    "GGGCCCGGGG": 1.3471,
    "GGGCCCGGGT": 1.1059428571428571,
    "GGGGCCCCCA": 2.0881714285714286,
    "GGGGCCCCCC": 2.1546999999999996,
    "GGGGCCCCCG": 1.9135428571428572,
    "GGGGCCCCCT": 2.321028571428571,
    "GGGGCCCCGG": 1.5882571428571428,
    "GGGGCCCCTA": 2.185157142857143,
    "GGGGCCCGGG": 1.3471,
    "GGGGGGCCCA": 2.0881714285714286,
    "GGGGGGCCCC": 2.1546999999999996,
    "GGGGGGCCCG": 1.9135428571428572,
    "GGGGGGCCCT": 2.321028571428571,
    "GGGGGGGGCC": 2.173842857142857,
    "GGGGGGGGGC": 2.5169571428571422,
    "GGGGGGGGGG": 2.5361,
    "GGGGGGGGGT": 2.294942857142857,
    "TAGGGCCCCA": 2.1186285714285713,
    "TAGGGCCCCC": 2.185157142857143,
    "TAGGGCCCCG": 1.9440000000000002,
    "TAGGGCCCCT": 2.3514857142857144,
    "TAGGGCCCGG": 1.6187142857142856,
    "TAGGGCCCTA": 2.2156142857142855,
    "TAGGGGCCCA": 2.1186285714285713,
    "TAGGGGCCCC": 2.185157142857143,
    "TAGGGGCCCG": 1.9440000000000002,
    "TAGGGGCCCT": 2.3514857142857144,
    "TAGGGGGCCC": 2.185157142857143,
    "TAGGGGGGCC": 2.2043,
    "TAGGGGGGGC": 2.5474142857142854,
    "TAGGGGGGGG": 2.5665571428571425,
    "TAGGGGGGGT": 2.3253999999999997,
    "TGGGCCCCCA": 2.021642857142857,
    "TGGGCCCCCC": 2.0881714285714286,
    "TGGGCCCCCG": 1.8470142857142857,
    "TGGGCCCCCT": 2.2544999999999997,
    "TGGGCCCCGG": 1.5217285714285713,
    "TGGGCCCCTA": 2.1186285714285713,
    "TGGGCCCGGG": 1.2805714285714287,
    "TGGGGCCCCA": 2.021642857142857,
    "TGGGGCCCCC": 2.0881714285714286,
    "TGGGGCCCCG": 1.8470142857142857,
    "TGGGGCCCCT": 2.2544999999999997,
    "TGGGGCCCGG": 1.5217285714285713,
    "TGGGGCCCTA": 2.1186285714285713,
    "TGGGGGCCCA": 2.021642857142857,
    "TGGGGGCCCC": 2.0881714285714286,
    "TGGGGGCCCG": 1.8470142857142857,
    "TGGGGGCCCT": 2.2544999999999997,
    "TGGGGGGCCC": 2.0881714285714286,
    "TGGGGGGGCC": 2.1073142857142857,
    "TGGGGGGGGC": 2.450428571428571,
    "TGGGGGGGGG": 2.4695714285714283,
    "TGGGGGGGGT": 2.2284142857142855,
}

# === Thresholds ===
# | Parameter      | Type  | Value | Description                                |
# |----------------|-------|-------|-------------------------------------------|
# | STRONG_LOG2    | float | 2.0   | Threshold for strong A-philic motifs     |
# | HIGH_MEAN      | float | 2.0   | Mean score threshold for high confidence |
# | MOD_MEAN       | float | 1.0   | Mean score threshold for moderate conf   |
# | FRACTION_HIGH  | float | 5/7   | Fraction threshold for high class        |
# | FRACTION_MOD   | float | 4/7   | Fraction threshold for moderate class    |

STRONG_LOG2 = 2.0
HIGH_MEAN = 2.0
MOD_MEAN = 1.0
FRACTION_HIGH = 5/7
FRACTION_MOD = 4/7

# === Hyperscan check ===
try:
    import hyperscan
    # Disable hyperscan for now due to API compatibility issues
    USE_HS = False  # Set to False until hyperscan API is fixed
except ImportError:
    USE_HS = False


def calculate_aphilic_score(sequence: str) -> float:
    """
    Calculate A-philic score for a sequence using 10-mer log2 odds table.
    
    | Parameter | Type | Description                    | Range        |
    |-----------|------|--------------------------------|--------------|
    | sequence  | str  | DNA sequence to score          | >=10 bp      |
    | return    | float| Sum of log2 odds scores        | 0.0 - 50.0   |
    """
    if len(sequence) < 10:
        return 0.0
    
    total_score = 0.0
    count = 0
    
    for i in range(len(sequence) - 9):
        nmer = sequence[i:i+10]
        if nmer in APHILIC_NMER_LOG2:
            total_score += APHILIC_NMER_LOG2[nmer]
            count += 1
    
    return total_score if count == 0 else total_score / count


def classify_aphilic_region(score: float, length: int, strong_count: int) -> str:
    """
    Classify A-philic region based on score and characteristics.
    
    | Parameter    | Type | Description                      | Range      |
    |--------------|------|----------------------------------|------------|
    | score        | float| Mean log2 odds score            | 0.0 - 3.0  |
    | length       | int  | Length of sequence region       | >=10       |
    | strong_count | int  | Number of strong scoring 10-mers | 0 - length |
    | return       | str  | Classification label             | str        |
    """
    n_possible = max(1, length - 9)  # Number of possible 10-mers
    strong_fraction = strong_count / n_possible
    
    if score >= HIGH_MEAN and strong_fraction >= FRACTION_HIGH:
        return "High-confidence A-philic"
    elif score >= MOD_MEAN and strong_fraction >= FRACTION_MOD:
        return "Moderate A-philic"
    else:
        return "Low A-philic"


def scan_sequence_internal(seq: str, min_len: int = 10, max_len: int = 15, step: int = 1) -> List[Dict[str, Any]]:
    """
    Internal function to scan sequence for A-philic motifs.
    
    | Parameter | Type | Description                    | Default | Range      |
    |-----------|------|--------------------------------|---------|------------|
    | seq       | str  | DNA sequence to scan           | -       | any        |
    | min_len   | int  | Minimum window length          | 10      | 10-50      |
    | max_len   | int  | Maximum window length          | 15      | min_len+   |
    | step      | int  | Step size for sliding window   | 1       | 1-5        |
    | return    | list | List of motif dictionaries     | []      | any        |
    """
    seq = seq.upper().replace('U', 'T')
    motifs = []
    
    # Scan with sliding windows
    for window_len in range(min_len, min(max_len + 1, len(seq) + 1)):
        for start in range(0, len(seq) - window_len + 1, step):
            window = seq[start:start + window_len]
            
            # Calculate score and characteristics
            score = calculate_aphilic_score(window)
            strong_count = sum(1 for i in range(len(window) - 9) 
                             if window[i:i+10] in APHILIC_NMER_LOG2 
                             and APHILIC_NMER_LOG2[window[i:i+10]] >= STRONG_LOG2)
            
            classification = classify_aphilic_region(score, window_len, strong_count)
            
            # Only keep motifs with some A-philic potential
            if score >= MOD_MEAN * 0.5:  # Relaxed threshold for detection
                motifs.append({
                    'Class': 'A-philic DNA',
                    'Subclass': 'A-philic',
                    'Start': start + 1,  # 1-based coordinates
                    'End': start + window_len,
                    'Length': window_len,
                    'Sequence': window,
                    'Raw_Score': score,
                    'Score': score,
                    'Score_Method': 'Tetranucleotide_Log2_Odds',
                    'N_Tetranucleotides': max(0, window_len - 9),
                    'Strong_Count': strong_count,
                    'Classification': classification,
                    'Confidence': classification.split('-')[0]
                })
    
    return motifs


def find_a_philic_dna(sequence: str, sequence_name: str = "Unknown") -> List[Dict[str, Any]]:
    """
    Main function to find A-philic DNA motifs in a sequence.
    Compatible with NBDFinder motif detection framework.
    
    Uses new 10-mer A-philic patterns with tetranucleotide log2 odds scoring.
    
    | Parameter     | Type | Description                      | Default   | Range      |
    |---------------|------|----------------------------------|-----------|------------|
    | sequence      | str  | DNA sequence to analyze          | -         | >=10 bp    |
    | sequence_name | str  | Name/identifier for sequence     | "Unknown" | any        |
    | return        | list | Standardized A-philic motifs     | []        | any        |
    """
    if not sequence or len(sequence) < 10:  # Changed to 10 for 10-mers
        return []
    
    # Clean sequence
    sequence = sequence.upper().replace('U', 'T')
    
    # Scan for A-philic motifs using NBDFinder detection
    motifs = detect_aphilic_motifs(sequence, sequence_name)
    
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
    
    | Parameter | Type | Description                    | Default           | Range      |
    |-----------|------|--------------------------------|-------------------|------------|
    | seq       | str  | DNA sequence to scan           | -                 | any        |
    | min_len   | int  | Minimum window length          | 10                | 10-50      |
    | max_len   | int  | Maximum window length          | 15                | min_len+   |
    | step      | int  | Step size for sliding window   | 1                 | 1-5        |
    | out_tsv   | str  | Output TSV filename            | "aphilic_hits.tsv"| any        |
    | return    | str  | Path to output TSV file        | out_tsv           | any        |
    """
    motifs = find_a_philic_dna(seq, "scan_sequence")
    
    # Write to TSV file
    with open(out_tsv, "w", newline="") as fh:
        import csv
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["start", "end", "length", "n_tets", "sum_log2", "strong_count", "label", "seq"])
        
        for motif in motifs:
            # Convert to legacy format
            n_tets = motif.get('Length', 10) - 3  # Tetranucleotides in the sequence
            strong_count = 1 if motif.get('Raw_Score', 0) >= STRONG_LOG2 else 0
            
            w.writerow([
                motif.get('Start', 0),
                motif.get('End', 0), 
                motif.get('Length', 0),
                n_tets,
                motif.get('Raw_Score', 0),
                strong_count,
                motif.get('Confidence', 'Unknown'),
                motif.get('Sequence', '')
            ])
    
    return out_tsv


# === Example usage and testing ===
if __name__ == "__main__":
    example = "NNNNAGGGGGGGGGCCCCTGGGGGCCCAAGGGNNNN"
    motifs = find_a_philic_dna(example, "example_sequence")
    
    print(f"Found {len(motifs)} A-philic DNA motifs:")
    for motif in motifs:
        score = motif.get('Raw_Score', motif.get('Score', 'N/A'))
        print(f"  {motif.get('Subclass', 'Unknown')}: {motif.get('Start', 0)}-{motif.get('End', 0)} (score: {score})")
        print(f"    Sequence: {motif.get('Sequence', 'N/A')}")
    
    # Also test legacy TSV output
    out_file = scan_sequence(example, min_len=10, max_len=20, step=1,
                           out_tsv="aphilic_results.tsv")
    print(f"\nResults also written to: {out_file}")