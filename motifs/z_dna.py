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
from .base_motif import wrap, standardize_motif_output
from .hyperscan_manager import optimized_hs_find
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from core.regex_registry import get_patterns_for_motif

# === Load patterns from central registry ===
ZDNA_PATTERNS = get_patterns_for_motif('z_dna')

# === Z-DNA 10-mer Pattern Database ===
ZDNA_10MER_SCORES = {
    "AACGCGCGCG": 50.25,
    "ACCGCGCGCG": 50.25,
    "ACGCCGCGCG": 50.25,
    "ACGCGCCGCG": 50.25,
    "ACGCGCGCCG": 50.25,
    "ACGCGCGCGA": 50.25,
    "ACGCGCGCGC": 57.25,
    "ACGCGCGCGG": 50.25,
    "ACGCGCGCGT": 51.5,
    "ACGCGGCGCG": 50.25,
    "ACGGCGCGCG": 50.25,
    "AGCGCGCGCA": 50.25,
    "AGCGCGCGCG": 56,
    "ATCGCGCGCG": 50,
    "ATGCGCGCGC": 51.25,
    "CAGCGCGCGC": 50.25,
    "CACGCGCGCG": 51.5,
    "CCGCGCGCGC": 56,
    "CCGCGCGCGT": 50.25,
    "CGACGCGCGC": 50.25,
    "CGCCGCGCGC": 56,
    "CGCCGCGCGT": 50.25,
    "CGCACGCGCG": 51.5,
    "CGCAGCGCGC": 50.25,
    "CGCGACGCGC": 50.25,
    "CGCGCACGCG": 51.5,
    "CGCGCAGCGC": 50.25,
    "CGCGCCGCGC": 56,
    "CGCGCCGCGT": 50.25,
    "CGCGCGACGC": 50.25,
    "CGCGCGCACG": 51.5,
    "CGCGCGCAGC": 50.25,
    "CGCGCGCGAC": 50.25,
    "CGCGCGCGAT": 50,
    "CGCGCGCGCA": 57.25,
    "CGCGCGCGCC": 56,
    "CGCGCGCGCG": 63,
    "CGCGCGCGCT": 56,
    "CGCGCGCGGC": 56,
    "CGCGCGCGGT": 50.25,
    "CGCGCGCGTC": 50.25,
    "CGCGCGCGTG": 51.5,
    "CGCGCGCGTT": 50.25,
    "CGCGCGCGTA": 51.25,
    "CGCGCGCCGC": 56,
    "CGCGCGCCGT": 50.25,
    "CGCGCGGCGC": 56,
    "CGCGCGGCGT": 50.25,
    "CGCGCGTCGC": 50.25,
    "CGCGCGTGCG": 51.5,
    "CGCGCTGCGC": 50.25,
    "CGCGGCGCGC": 56,
    "CGCGGCGCGT": 50.25,
    "CGCGTCGCGC": 50.25,
    "CGCGTGCGCG": 51.5,
    "CGCTGCGCGC": 50.25,
    "CGGCGCGCGC": 56,
    "CGGCGCGCGT": 50.25,
    "CGTCGCGCGC": 50.25,
    "CGTGCGCGCG": 51.5,
    "CTGCGCGCGC": 50.25,
    "GACGCGCGCG": 50.25,
    "GCACGCGCGC": 51.5,
    "GCAGCGCGCG": 50.25,
    "GCCGCGCGCA": 50.25,
    "GCCGCGCGCG": 56,
    "GCGACGCGCG": 50.25,
    "GCGCACGCGC": 51.5,
    "GCGCAGCGCG": 50.25,
    "GCGCCGCGCA": 50.25,
    "GCGCCGCGCG": 56,
    "GCGCGACGCG": 50.25,
    "GCGCGCACGC": 51.5,
    "GCGCGCAGCG": 50.25,
    "GCGCGCGAAA": 50.25,
    "GCGCGCGAAC": 51.5,
    "GCGCGCGAAG": 50.25,
    "GCGCGCGACA": 50.25,
    "GCGCGCGACG": 50.25,
    "GCGCGCGCAG": 50.25,
    "GCGCGCGCAT": 51.25,
    "GCGCGCGCCA": 50.25,
    "GCGCGCGCCG": 56,
    "GCGCGCGCGA": 56,
    "GCGCGCGCGC": 63,
    "GCGCGCGCGG": 56,
    "GCGCGCGCGT": 57.25,
    "GCGCGCGCTA": 50,
    "GCGCGCGCTG": 50.25,
    "GCGCGCGGCA": 50.25,
    "GCGCGCGGCG": 56,
    "GCGCGCGTCG": 50.25,
    "GCGCGCGTGC": 51.5,
    "GCGCGCTGCG": 50.25,
    "GCGCGGCGCA": 50.25,
    "GCGCGGCGCG": 56,
    "GCGCGTCGCG": 50.25,
    "GCGCGTGCGC": 51.5,
    "GCGCTGCGCG": 50.25,
    "GCGGCGCGCA": 50.25,
    "GCGGCGCGCG": 56,
    "GCGTCGCGCG": 50.25,
    "GCGTGCGCGC": 51.5,
    "GCTGCGCGCG": 50.25,
    "GGCGCGCGCA": 50.25,
    "GGCGCGCGCG": 56,
    "GTCGCGCGCG": 50.25,
    "GTGCGCGCGC": 51.5,
    "TACGCGCGCG": 51.25,
    "TAGCGCGCGC": 50,
    "TCGCGCGCGC": 56,
    "TCGCGCGCGT": 50.25,
    "TGCCGCGCGC": 50.25,
    "TGCGCCGCGC": 50.25,
    "TGCGCGCCGC": 50.25,
    "TGCGCGCGCA": 51.5,
    "TGCGCGCGCC": 50.25,
    "TGCGCGCGCG": 57.25,
    "TGCGCGCGCT": 50.25,
    "TGCGCGCGGC": 50.25,
    "TGCGCGGCGC": 50.25,
    "TGCGGCGCGC": 50.25,
    "TGGCGCGCGC": 50.25,
    "TTGCGCGCGC": 50.25,
}

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

#--- Z-DNA 10-mer motif finder using Hyperscan for fast pattern matching ---
def find_zdna_10mer_patterns(seq: str) -> list:
    """
    Find Z-DNA 10-mer patterns using Hyperscan-based fast pattern matching.
    
    Args:
        seq: DNA sequence to scan
        
    Returns:
        List of Z-DNA 10-mer motif matches with positions and scores
    """
    if len(seq) < 10:
        return []
    
    seq = seq.upper()
    matches = []
    
    # Try Hyperscan first for fast pattern matching
    try:
        import hyperscan
        
        # Create patterns for all Z-DNA 10-mers
        patterns = []
        motif_scores = {}
        
        for i, (motif, score) in enumerate(ZDNA_10MER_SCORES.items()):
            patterns.append((motif.encode(), i, hyperscan.HS_FLAG_CASELESS))
            motif_scores[i] = (motif, score)
        
        # Compile the database 
        try:
            db = hyperscan.Database(
                expressions=[p[0] for p in patterns],
                ids=[p[1] for p in patterns],
                flags=[p[2] for p in patterns]
            )
            
            # Create scratch space
            scratch = hyperscan.Scratch(db)
            
        except Exception as compile_error:
            print(f"Hyperscan compilation failed, using regex fallback: {compile_error}")
            return find_zdna_10mer_regex(seq)
        
        # Scan the sequence
        def match_handler(match_id, start, end, flags, context):
            motif, score = motif_scores[match_id]
            matches.append({
                'Class': 'Z-DNA',
                'Subclass': 'Z-DNA-10mer',
                'Start': start + 1,  # 1-based coordinates
                'End': end,
                'Length': len(motif),
                'Sequence': motif,
                'Score': float(score),
                'ScoreMethod': '10mer-Pattern'
            })
        
        # Perform the scan
        db.scan(seq.encode(), match_handler, scratch)
        return matches
        
    except ImportError:
        print("Hyperscan not available, using regex fallback")
        return find_zdna_10mer_regex(seq)
    except Exception as e:
        print(f"Hyperscan failed, falling back to regex: {e}")
        return find_zdna_10mer_regex(seq)

def find_zdna_10mer_regex(seq: str) -> list:
    """
    Fallback regex-based search for Z-DNA 10-mer patterns.
    
    Args:
        seq: DNA sequence to scan
        
    Returns:
        List of Z-DNA 10-mer motif matches with positions and scores  
    """
    import re
    seq = seq.upper()
    matches = []
    
    for motif, score in ZDNA_10MER_SCORES.items():
        pattern = re.compile(motif)
        for match in pattern.finditer(seq):
            matches.append({
                'Class': 'Z-DNA',
                'Subclass': 'Z-DNA-10mer', 
                'Start': match.start() + 1,  # 1-based coordinates
                'End': match.end(),
                'Length': len(motif),
                'Sequence': motif,
                'Score': float(score),
                'ScoreMethod': '10mer-Pattern'
            })
    
    return matches

#--- Main: Find all Z-DNA and eGZ motifs, output standardized for genomic analysis ---
def find_z_dna(seq: str, sequence_name: str = "") -> list:
    zdna_results = find_zdna(seq)
    egz_results = find_egz_motif(seq)
    zdna_10mer_results = find_zdna_10mer_patterns(seq)
    
    all_results = zdna_results + egz_results + zdna_10mer_results
    return [standardize_motif_output(motif, sequence_name, i) for i, motif in enumerate(all_results, 1)]

#--- Annotations ---
# - zdna_seeker_scoring_array: core scoring array, weights/penalties from Z-DNA literature.
# - find_zdna: maximal-scoring subsequence (Z-seeker), with scientific threshold and drop logic.
# - find_egz_motif: Hyperscan block scan for eGZ (CGG)n, literature length/cutoff, scoring G-bias.
# - find_z_dna: combines both, output standardized 1-based for analysis.
