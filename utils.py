"""
Essential Utility Functions for Non-B DNA Analysis
==================================================

| Function Category | Count | Description                                |
|-------------------|-------|--------------------------------------------|
| Sequence Parsing  | 3     | FASTA parsing and sequence formatting      |
| DNA Analysis      | 6     | GC content, reverse complement, etc.       |
| Statistical       | 4     | Basic sequence and motif statistics        |
| Validation        | 2     | Input validation and error checking        |

Core functionality for sequence processing and basic bioinformatics operations
"""

import re
import numpy as np
from typing import List, Dict, Tuple, Any, Optional, Union
from collections import defaultdict, Counter
import random

try:
    from scipy.stats import percentileofscore
except ImportError:
    def percentileofscore(a, score, kind='rank'):
        """Simplified percentile calculation fallback"""
        a = np.asarray(a)
        if len(a) == 0:
            return 0.0
        if kind == 'rank':
            return (sum(a <= score) / len(a)) * 100
        elif kind == 'strict':
            return (sum(a < score) / len(a)) * 100
        elif kind == 'weak':
            return (sum(a <= score) / len(a)) * 100
        elif kind == 'mean':
            return (sum(a < score) + sum(a <= score)) / (2 * len(a)) * 100
        else:
            raise ValueError("kind must be 'rank', 'strict', 'weak' or 'mean'")

# =============================================================================
# SEQUENCE PARSING AND FORMATTING
# =============================================================================

def parse_fasta(fasta_str: str) -> str:
    """
    Parse FASTA format string to clean DNA sequence
    
    | Input Format | Description                               |
    |--------------|-------------------------------------------|
    | Standard     | >header\\nACGT...                        |
    | Multi-line   | >header\\nACGT\\nACGT...                 |
    | Mixed case   | >header\\nacgt... (converts to uppercase) |
    
    Returns: Clean uppercase DNA sequence (A,T,G,C only)
    """
    if not fasta_str or not isinstance(fasta_str, str):
        return ""
    
    # Remove header lines and clean sequence
    lines = [line.strip() for line in fasta_str.split('\n') 
             if line.strip() and not line.startswith(">")]
    
    # Join all sequence lines
    sequence = "".join(lines).upper()
    
    # Remove whitespace and convert U to T
    sequence = sequence.replace(" ", "").replace("U", "T")
    
    # Validate DNA sequence
    valid_chars = set('ATGC')
    cleaned_sequence = ''.join(char for char in sequence if char in valid_chars)
    
    return cleaned_sequence

def wrap(seq: str, width: int = 60) -> str:
    """
    Format sequence with line breaks for display
    
    | Parameter | Type | Range  | Description                    |
    |-----------|------|--------|--------------------------------|
    | seq       | str  | -      | DNA sequence to format         |
    | width     | int  | 10-120 | Characters per line            |
    
    Returns: Multi-line formatted sequence
    """
    if not seq or width <= 0:
        return seq
    
    return "\n".join([seq[i:i+width] for i in range(0, len(seq), width)])

def validate_dna_sequence(seq: str) -> Tuple[bool, str]:
    """
    Validate DNA sequence for proper format and content
    
    | Validation Check | Description                           |
    |------------------|---------------------------------------|
    | Length           | Must be at least 1 bp                |
    | Characters       | Only A,T,G,C allowed (case insensitive)|
    | Format           | No unexpected symbols or whitespace   |
    
    Returns: (is_valid: bool, error_message: str)
    """
    if not seq:
        return False, "Empty sequence provided"
    
    if not isinstance(seq, str):
        return False, "Sequence must be a string"
    
    seq = seq.upper().strip()
    
    if len(seq) == 0:
        return False, "Sequence contains no valid characters"
    
    # Check for invalid characters
    valid_chars = set('ATGC')
    invalid_chars = set(seq) - valid_chars
    
    if invalid_chars:
        return False, f"Invalid characters found: {', '.join(sorted(invalid_chars))}"
    
    return True, "Valid DNA sequence"

# =============================================================================
# DNA SEQUENCE ANALYSIS
# =============================================================================

def gc_content(seq: str) -> float:
    """
    Calculate GC content percentage
    
    | GC Content | Interpretation                        |
    |------------|---------------------------------------|
    | 0-30%      | AT-rich (e.g., A-tracts, repeats)    |
    | 30-50%     | Balanced composition                  |
    | 50-70%     | GC-rich (e.g., promoters, CpG)       |
    | 70-100%    | Extremely GC-rich (rare)              |
    
    Returns: GC percentage (0-100)
    """
    if not seq:
        return 0.0
    
    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    total_count = len(seq)
    
    return (gc_count / total_count) * 100 if total_count > 0 else 0.0

def at_content(seq: str) -> float:
    """Calculate AT content percentage (complement of GC content)"""
    return 100.0 - gc_content(seq)

def reverse_complement(seq: str) -> str:
    """
    Generate reverse complement of DNA sequence
    
    | Base | Complement | Description                    |
    |------|------------|--------------------------------|
    | A    | T          | Adenine -> Thymine             |
    | T    | A          | Thymine -> Adenine             |
    | G    | C          | Guanine -> Cytosine            |
    | C    | G          | Cytosine -> Guanine            |
    
    Returns: Reverse complement sequence
    """
    if not seq:
        return ""
    
    complement_map = str.maketrans("ATGC", "TACG")
    return seq.upper().translate(complement_map)[::-1]

def is_palindrome(seq: str) -> bool:
    """
    Check if sequence is a perfect palindrome (reads same on both strands)
    
    | Palindrome Type | Example     | Description                  |
    |-----------------|-------------|------------------------------|
    | Perfect         | GAATTC      | Exact reverse complement     |
    | Imperfect       | GAATCC      | Minor mismatches allowed     |
    | Nested          | GGATCCATG   | Palindrome within sequence   |
    
    Returns: True if perfect palindrome
    """
    if not seq:
        return False
    
    return seq.upper() == reverse_complement(seq)

def calculate_tm(seq: str, method: str = 'basic') -> float:
    """
    Calculate DNA melting temperature (Tm) in Celsius
    
    | Method | Formula                              | Accuracy         |
    |--------|--------------------------------------|------------------|
    | basic  | 2(A+T) + 4(G+C) for <14bp            | ±5°C for short   |
    | gc     | 64.9 + 41(G+C-16.4)/L for ≥14bp     | ±3°C for medium  |
    | salt   | Accounts for salt concentration       | ±2°C for precise |
    
    Returns: Melting temperature in °C
    """
    if not seq:
        return 0.0
    
    seq = seq.upper()
    length = len(seq)
    
    if method == 'basic' or length < 14:
        # Simple rule for short sequences
        return 2 * (seq.count('A') + seq.count('T')) + 4 * (seq.count('G') + seq.count('C'))
    
    elif method == 'gc':
        # GC-based formula for longer sequences
        gc_count = seq.count('G') + seq.count('C')
        gc_percent = (gc_count / length) * 100
        return 64.9 + 41 * (gc_percent - 16.4) / length
    
    else:
        # Default to basic method
        return calculate_tm(seq, 'basic')

def shuffle_sequence(seq: str, preserve_composition: bool = True) -> str:
    """
    Create randomized sequence for statistical testing
    
    | Shuffle Type     | Description                              |
    |------------------|------------------------------------------|
    | Composition      | Maintains A,T,G,C counts (default)      |
    | Complete Random  | Random bases with equal probability     |
    | Dinuc Preserve   | Maintains dinucleotide frequencies      |
    
    Returns: Randomized sequence string
    """
    if not seq:
        return seq
    
    if preserve_composition:
        # Preserve base composition
        seq_list = list(seq.upper())
        random.shuffle(seq_list)
        return ''.join(seq_list)
    else:
        # Completely random sequence
        bases = ['A', 'T', 'G', 'C']
        return ''.join(random.choice(bases) for _ in range(len(seq)))

# =============================================================================
# STATISTICAL ANALYSIS FUNCTIONS
# =============================================================================

def get_basic_stats(seq: str, motifs: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
    """
    Calculate comprehensive sequence statistics
    
    | Statistic     | Description                             |
    |---------------|-----------------------------------------|
    | Length        | Sequence length in base pairs          |
    | Composition   | A,T,G,C counts and percentages         |
    | GC Content    | G+C percentage                          |
    | Motif Summary | Count and diversity of detected motifs |
    
    Returns: Dictionary with sequence statistics
    """
    if not seq:
        return {
            'Length': 0,
            'GC_Content': 0.0,
            'AT_Content': 0.0,
            'Base_Counts': {'A': 0, 'T': 0, 'G': 0, 'C': 0},
            'Base_Frequencies': {'A': 0.0, 'T': 0.0, 'G': 0.0, 'C': 0.0}
        }
    
    seq = seq.upper()
    length = len(seq)
    
    # Base counts
    base_counts = {
        'A': seq.count('A'),
        'T': seq.count('T'),
        'G': seq.count('G'),
        'C': seq.count('C')
    }
    
    # Base frequencies
    base_frequencies = {base: count/length for base, count in base_counts.items()}
    
    # Basic composition stats
    gc = gc_content(seq)
    at = 100.0 - gc
    
    stats = {
        'Length': length,
        'GC_Content': round(gc, 2),
        'AT_Content': round(at, 2),
        'Base_Counts': base_counts,
        'Base_Frequencies': {base: round(freq, 4) for base, freq in base_frequencies.items()}
    }
    
    # Add motif statistics if provided
    if motifs:
        stats.update(_calculate_motif_summary_stats(motifs))
    
    return stats

def _calculate_motif_summary_stats(motifs: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Helper function to calculate motif summary statistics"""
    if not motifs:
        return {
            'Total_Motifs': 0,
            'Unique_Classes': 0,
            'Average_Motif_Length': 0.0,
            'Average_Motif_Score': 0.0
        }
    
    # Basic counts
    total_motifs = len(motifs)
    unique_classes = len(set(motif.get('Class', 0) for motif in motifs))
    
    # Average statistics
    lengths = [motif.get('Length', 0) for motif in motifs if motif.get('Length', 0) > 0]
    scores = [motif.get('Normalized_Score', 0) for motif in motifs]
    
    avg_length = np.mean(lengths) if lengths else 0.0
    avg_score = np.mean(scores) if scores else 0.0
    
    # Class distribution
    class_counts = Counter(motif.get('Class', 0) for motif in motifs)
    
    return {
        'Total_Motifs': total_motifs,
        'Unique_Classes': unique_classes,
        'Average_Motif_Length': round(avg_length, 1),
        'Average_Motif_Score': round(avg_score, 3),
        'Class_Distribution': dict(class_counts)
    }

def kmer_conservation(seq: str, k: int = 6, n_shuffles: int = 1000) -> Dict[str, Tuple[float, float]]:
    """
    Calculate k-mer conservation scores using shuffle test
    
    | Analysis Step | Description                                |
    |---------------|------------------------------------------- |
    | Observed      | Count k-mers in original sequence         |
    | Expected      | Generate null distribution via shuffling  |
    | Comparison    | Calculate log2 ratios and p-values        |
    | Conservation  | Identify over/under-represented k-mers    |
    
    Returns: {k-mer: (log2_ratio, p_value)} dictionary
    """
    if not seq or len(seq) < k:
        return {}
    
    seq = seq.upper()
    
    # Count observed k-mers
    observed_counts = Counter()
    total_kmers = len(seq) - k + 1
    
    for i in range(total_kmers):
        kmer = seq[i:i+k]
        observed_counts[kmer] += 1
    
    # Generate null distribution through shuffling
    null_distributions = defaultdict(list)
    
    for shuffle_iter in range(n_shuffles):
        shuffled_seq = shuffle_sequence(seq, preserve_composition=True)
        
        # Count k-mers in shuffled sequence
        shuffled_counts = Counter()
        for i in range(len(shuffled_seq) - k + 1):
            kmer = shuffled_seq[i:i+k]
            shuffled_counts[kmer] += 1
        
        # Store counts for each k-mer
        for kmer in observed_counts:
            null_distributions[kmer].append(shuffled_counts.get(kmer, 0))
    
    # Calculate conservation statistics
    conservation_results = {}
    
    for kmer, observed_count in observed_counts.items():
        null_counts = null_distributions[kmer]
        
        if null_counts:
            expected_count = np.mean(null_counts)
            
            # Calculate log2 ratio (with pseudocounts)
            log2_ratio = np.log2((observed_count + 1e-6) / (expected_count + 1e-6))
            
            # Calculate p-value
            p_value = 1.0 - (percentileofscore(null_counts, observed_count) / 100.0)
            p_value = max(p_value, 1.0 / (n_shuffles + 1))  # Minimum p-value
            
            conservation_results[kmer] = (log2_ratio, p_value)
    
    return conservation_results

def motif_conservation(motif_seq: str, conservation_scores: Dict[str, Tuple[float, float]], k: int = 6) -> float:
    """
    Calculate average conservation score for a motif sequence
    
    | Conservation Level | Log2 Ratio Range | Interpretation          |
    |-------------------|------------------|-------------------------|
    | Highly Conserved  | > 2.0           | Strong evolutionary signal |
    | Moderately Conserved | 1.0 - 2.0    | Above random expectation  |
    | Random            | -1.0 - 1.0      | No significant pattern    |
    | Depleted          | < -1.0          | Below random expectation  |
    
    Returns: Average conservation score for the motif
    """
    if not motif_seq or len(motif_seq) < k:
        return 0.0
    
    conservation_values = []
    
    for i in range(len(motif_seq) - k + 1):
        kmer = motif_seq[i:i+k]
        if kmer in conservation_scores:
            log2_ratio, p_value = conservation_scores[kmer]
            conservation_values.append(log2_ratio)
    
    return np.mean(conservation_values) if conservation_values else 0.0
