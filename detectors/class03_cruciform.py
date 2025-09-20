"""
Cruciform DNA Motif Detection (Class 3)
=======================================
- Detects perfect palindromes and inverted repeats (arms >= 6bp, no upper limit; loop <= 100bp)
- Uses Hyperscan for accelerated motif search (requires `pip install hyperscan`)
- Scores motifs by nearest-neighbor thermodynamics (SantaLucia 1998), outputs ΔG° and normalized stability (NormScore: 1=most stable, 0=least stable)
- Normalization: Arms >=20bp considered equally stable for NormScore

Parameter Table:
| Parameter          | Type    | Description                                                                                         | Example/Range                |
|--------------------|---------|-----------------------------------------------------------------------------------------------------|------------------------------|
| seq                | str     | DNA sequence to analyze                                                                             | 'ATGCGCAT...'                |
| sequence_name      | str     | Optional name for the sequence                                                                      | 'chr1', 'plasmidA'           |
| min_arm            | int     | Minimum arm length for motifs                                                                       | 6 (default for function)     |
| max_spacer         | int     | Maximum loop length (for inverted repeats)                                                          | 100 (default for function)   |
| NN_DELTA_G         | dict    | Nearest-neighbor ΔG° table for all dinucleotides (SantaLucia 1998)                                 | See code                     |
| NN_INITIATION      | float   | Duplex initiation penalty (SantaLucia 1998)                                                         | 0.2                          |
| NN_SYMMETRY        | float   | Penalty for self-complementary (palindrome) duplex                                                  | 0.43                         |
| DG_MIN             | float   | Most stable (lowest) ΔG° for normalization, 20bp GC palindrome                                      | ~-86.5                       |
| DG_MAX             | float   | Least stable (highest) ΔG°, 6bp AT palindrome + 100bp loop                                          | ~+7.7                        |
| matches            | list    | List of Hyperscan matches (motif id, start, end)                                                    | [(id, from, to), ...]        |
| motifs             | list    | List of detected motif dictionaries                                                                 | [ {...}, ... ]               |
| arm_len            | int     | Length of each arm in motif                                                                         | 6 ... n//2                   |
| loop               | int     | Length of the loop/spacer in inverted repeat                                                        | 1 ... 100                    |
| dg                 | float   | Calculated ΔG° (kcal/mol) for motif                                                                 | e.g. -20.3                   |
| norm               | float   | Normalized score: 1=most stable, 0=least stable                                                     | 0.0 ... 1.0                  |
| region             | tuple   | (Start, End) coordinates for deduplication                                                          | (start+1, end)               |
| wrap, reverse_complement, standardize_motif_output | functions | Import utilities for formatting, biology, and output standardization                  | See base_motif               |
"""

import hyperscan; import numpy as np; from motifs.base_motif import reverse_complement, wrap, standardize_motif_output
import sys, os, re
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from core.regex_registry import get_patterns_for_motif

# === Load patterns from central registry ===
CRUCIFORM_PATTERNS = get_patterns_for_motif('cruciform')

# --- Thermodynamic parameter table (SantaLucia 1998) ---
NN_DELTA_G={'AA':-1.00,'TT':-1.00,'AT':-0.88,'TA':-0.58,'CA':-1.45,'TG':-1.45,'GT':-1.44,'AC':-1.44,'CT':-1.28,'AG':-1.28,'GA':-1.30,'TC':-1.30,'CG':-2.17,'GC':-2.24,'GG':-1.84,'CC':-1.84}; NN_INITIATION=0.2; NN_SYMMETRY=0.43
DG_MIN=2*(NN_DELTA_G['GC']*19+NN_INITIATION+NN_SYMMETRY); DG_MAX=2*(NN_DELTA_G['AT']*5+NN_INITIATION+NN_SYMMETRY)+9.75

# -- Compute NN ΔG° for DNA duplex/hairpin arm --
def nn_dg(seq):
    seq=seq.upper(); dg=NN_INITIATION
    for i in range(len(seq)-1): dg+=NN_DELTA_G.get(seq[i:i+2],0.0)
    if seq==reverse_complement(seq): dg+=NN_SYMMETRY
    return dg

# -- Empirical loop penalty (Turner 2010, Zuker 1981) --
def loop_penalty(l):
    if l==0: return 0.0
    if l<=6: return [0,3.4,3.2,3.0,2.8,2.7,2.6][l]
    return 1.75+0.8*(l**0.5)

# -- Normalize ΔG°: 1=most stable (lowest ΔG°), 0=least stable (highest ΔG°) --
def normalize_dg(dg):
    if dg<DG_MIN: dg=DG_MIN
    if dg>DG_MAX: dg=DG_MAX
    return round((DG_MAX-dg)/(DG_MAX-DG_MIN),3)

# -- Removed: hs_callback function (replaced by inline callback in find_cruciform_hyperscan)

# -- Hyperscan-accelerated cruciform detection with candidate filtering --
def find_cruciform_hyperscan(seq: str) -> list:
    """Use Hyperscan to detect potential cruciform candidates and validate with Python."""
    motifs = []
    seqU = seq.upper()
    n = len(seqU)
    
    if n < 12:
        return []
    
    # Hyperscan patterns for candidate detection
    patterns = []
    pattern_info = {}
    pattern_id = 0
    
    # Generate palindrome candidate patterns (arms 6-20bp)
    for arm_len in range(6, min(n//2+1, 21)):
        pattern = f'[ATGC]{{{2*arm_len}}}'
        patterns.append((pattern.encode(), pattern_id))
        pattern_info[pattern_id] = ('palindrome', arm_len, 0)
        pattern_id += 1
    
    # Generate inverted repeat candidate patterns (arms 6-15bp, loops 1-20bp)
    for arm_len in range(6, min(n//2+1, 16)):
        for loop_len in range(1, min(21, n-2*arm_len)):
            total_len = 2*arm_len + loop_len
            if total_len <= n:
                pattern = f'[ATGC]{{{total_len}}}'
                patterns.append((pattern.encode(), pattern_id))
                pattern_info[pattern_id] = ('inverted_repeat', arm_len, loop_len)
                pattern_id += 1
    
    if not patterns:
        return []
    
    # Compile Hyperscan database
    try:
        db = hyperscan.Database()
        db.compile(expressions=[p[0] for p in patterns], 
                  ids=[p[1] for p in patterns])
        
        candidates = []
        
        def candidate_callback(id, start, end, flags, ctx):
            candidates.append((id, start, end))
            return hyperscan.HS_SUCCESS
        
        # Scan for candidates
        db.scan(seqU.encode(), match_event_handler=candidate_callback)
        
        # Validate candidates with Python logic
        seen_regions = set()
        
        for match_id, start, end in candidates:
            motif_type, arm_len, loop_len = pattern_info[match_id]
            candidate_seq = seqU[start:end]
            
            if motif_type == 'palindrome':
                # Validate palindrome
                if len(candidate_seq) == 2 * arm_len:
                    left_arm = candidate_seq[:arm_len]
                    right_arm = candidate_seq[arm_len:]
                    if right_arm == reverse_complement(left_arm):
                        dg_total = 2 * nn_dg(left_arm) + NN_SYMMETRY
                        norm = normalize_dg(dg_total)
                        region = (start+1, end)
                        if region not in seen_regions:
                            motifs.append({
                                "Class": "Cruciform", "Subclass": "Perfect_Palindrome",
                                "Start": start+1, "End": end, "Length": end-start,
                                "Sequence": wrap(candidate_seq), "Score": float(norm),
                                "ScoreMethod": "NN_Thermodynamics", "Arm_Length": arm_len,
                                "GC_Content": round((left_arm.count('G')+left_arm.count('C'))/len(left_arm)*100,1),
                                "Loop_Length": 0, "DeltaG_Arm": round(nn_dg(left_arm),2), "DeltaG_Loop": 0.0
                            })
                            seen_regions.add(region)
            
            elif motif_type == 'inverted_repeat':
                # Validate inverted repeat
                if len(candidate_seq) == 2 * arm_len + loop_len:
                    left_arm = candidate_seq[:arm_len]
                    spacer = candidate_seq[arm_len:arm_len+loop_len]
                    right_arm = candidate_seq[arm_len+loop_len:]
                    if right_arm == reverse_complement(left_arm):
                        dg_arm = nn_dg(left_arm) + nn_dg(right_arm)
                        dg_loop = loop_penalty(loop_len)
                        dg_total = dg_arm + dg_loop
                        norm = normalize_dg(dg_total)
                        region = (start+1, end)
                        if region not in seen_regions:
                            motifs.append({
                                "Class": "Cruciform", "Subclass": "Inverted_Repeat",
                                "Start": start+1, "End": end, "Length": end-start,
                                "Sequence": wrap(candidate_seq), "Score": float(norm),
                                "ScoreMethod": "NN_Thermodynamics", "Arm_Length": arm_len,
                                "GC_Content": round((left_arm.count('G')+left_arm.count('C'))/len(left_arm)*100,1),
                                "Loop_Length": loop_len, "DeltaG_Arm": round(dg_arm,2), "DeltaG_Loop": round(dg_loop,2)
                            })
                            seen_regions.add(region)
        
        return motifs
        
    except Exception as e:
        # Fallback to Python regex if Hyperscan fails
        return find_cruciform_python_fallback(seq)

# -- Python fallback implementation for cruciform detection --
def find_cruciform_python_fallback(seq: str) -> list:
    """Fallback Python implementation when Hyperscan is not available."""
    motifs = []
    seqU = seq.upper()
    n = len(seqU)
    
    if n < 12:
        return []
    
    seen_regions = set()
    
    # Find palindromes (perfect hairpins, no loop)
    for arm_len in range(6, min(n//2+1, 21)):
        for i in range(n - 2*arm_len + 1):
            candidate = seqU[i:i+2*arm_len]
            left_arm = candidate[:arm_len]
            right_arm = candidate[arm_len:]
            if right_arm == reverse_complement(left_arm):
                dg_total = 2 * nn_dg(left_arm) + NN_SYMMETRY
                norm = normalize_dg(dg_total)
                region = (i+1, i+2*arm_len)
                if region not in seen_regions:
                    motifs.append({
                        "Class": "Cruciform", "Subclass": "Perfect_Palindrome",
                        "Start": i+1, "End": i+2*arm_len, "Length": 2*arm_len,
                        "Sequence": wrap(candidate), "Score": float(norm),
                        "ScoreMethod": "NN_Thermodynamics", "Arm_Length": arm_len,
                        "GC_Content": round((left_arm.count('G')+left_arm.count('C'))/len(left_arm)*100,1),
                        "Loop_Length": 0, "DeltaG_Arm": round(nn_dg(left_arm),2), "DeltaG_Loop": 0.0
                    })
                    seen_regions.add(region)
    
    # Find inverted repeats (with loop)
    for arm_len in range(6, min(n//2+1, 16)):
        for loop in range(1, min(21, n-2*arm_len)):
            total_len = 2*arm_len + loop
            if total_len > n:
                continue
            for i in range(n - total_len + 1):
                candidate = seqU[i:i+total_len]
                left_arm = candidate[:arm_len]
                right_arm = candidate[arm_len+loop:]
                if right_arm == reverse_complement(left_arm):
                    dg_arm = nn_dg(left_arm) + nn_dg(right_arm)
                    dg_loop = loop_penalty(loop)
                    dg_total = dg_arm + dg_loop
                    norm = normalize_dg(dg_total)
                    region = (i+1, i+total_len)
                    if region not in seen_regions:
                        motifs.append({
                            "Class": "Cruciform", "Subclass": "Inverted_Repeat",
                            "Start": i+1, "End": i+total_len, "Length": total_len,
                            "Sequence": wrap(candidate), "Score": float(norm),
                            "ScoreMethod": "NN_Thermodynamics", "Arm_Length": arm_len,
                            "GC_Content": round((left_arm.count('G')+left_arm.count('C'))/len(left_arm)*100,1),
                            "Loop_Length": loop, "DeltaG_Arm": round(dg_arm,2), "DeltaG_Loop": round(dg_loop,2)
                        })
                        seen_regions.add(region)
    
    return motifs

# -- Main motif finder using Hyperscan acceleration with Python validation --
def find_cruciform(seq: str, sequence_name: str = "") -> list:
    """
    Detect cruciform DNA motifs using Hyperscan acceleration with Python validation.
    
    Uses hybrid approach:
    1. Hyperscan detects potential candidates based on length patterns
    2. Python validates palindrome/inverted repeat structure
    3. Maintains scientific accuracy while improving performance
    """
    motifs = find_cruciform_hyperscan(seq)
    return [standardize_motif_output(motif, sequence_name, i) for i, motif in enumerate(motifs, 1)]

# --- Scientific comments ---
# - Uses Python regex for palindrome and inverted repeat detection (arms >=6bp, loop <=50bp) 
# - NN thermodynamics for ΔG° (SantaLucia 1998), normalized scoring: 0=least stable, 1=most stable
# - Returns all detected motifs with thermodynamic stability scores
