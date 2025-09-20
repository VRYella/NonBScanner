"""
Triplex DNA Motif Detection (Class 5) - Hyperscan Accelerated
Subclasses: Triplex (5.1), Sticky DNA (5.2)
"""

import hyperscan, re
from motifs.base_motif import wrap, standardize_motif_output
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from core.regex_registry import get_patterns_for_motif

# === Load patterns from central registry ===
TRIPLEX_PATTERNS = get_patterns_for_motif('triplex')

#--- Purine/pyrimidine fraction calculators (per scientific convention) ---
def purine_fraction(seq): return (seq.count('A')+seq.count('G'))/max(1,len(seq))
def pyrimidine_fraction(seq): return (seq.count('C')+seq.count('T'))/max(1,len(seq))

#--- H-DNA/Triplex finder: hybrid hyperscan + Python regex approach ---
def find_hdna_hyperscan(seq: str) -> list:
    """Enhanced H-DNA detection using Hyperscan for homopurine/homopyrimidine tracts + Python for mirror repeats"""
    motifs = []
    seqU = seq.upper()
    n = len(seqU)
    
    if n < 20:  # Minimum useful size for triplex formation
        return []
    
    try:
        # Hyperscan patterns for homopurine/homopyrimidine tract detection
        patterns = []
        pattern_info = {}
        pattern_id = 0
        
        # Homopurine tracts (A/G only, 15+ bp)
        for tract_len in range(15, min(n+1, 101)):
            pattern = f'[AG]{{{tract_len}}}'
            patterns.append((pattern.encode(), pattern_id))
            pattern_info[pattern_id] = ('homopurine', tract_len)
            pattern_id += 1
        
        # Homopyrimidine tracts (C/T only, 15+ bp)  
        for tract_len in range(15, min(n+1, 101)):
            pattern = f'[CT]{{{tract_len}}}'
            patterns.append((pattern.encode(), pattern_id))
            pattern_info[pattern_id] = ('homopyrimidine', tract_len)
            pattern_id += 1
        
        if patterns:
            # Compile and scan with Hyperscan
            db = hyperscan.Database()
            db.compile(expressions=[p[0] for p in patterns], 
                      ids=[p[1] for p in patterns])
            
            candidates = []
            
            def candidate_callback(id, start, end, flags, ctx):
                candidates.append((id, start, end))
                return hyperscan.HS_SUCCESS
            
            db.scan(seqU.encode(), match_event_handler=candidate_callback)
            
            # Process homopurine/homopyrimidine candidates
            seen_regions = set()
            for match_id, start, end in candidates:
                tract_type, expected_len = pattern_info[match_id]
                tract_seq = seqU[start:end]
                
                if len(tract_seq) >= 15:  # Minimum biological relevance
                    pur_frac = purine_fraction(tract_seq)
                    pyr_frac = pyrimidine_fraction(tract_seq)
                    homogeneity = max(pur_frac, pyr_frac)
                    
                    # High homogeneity threshold for triplex formation
                    if homogeneity >= 0.9:
                        # Enhanced scoring for homogeneous tracts
                        base_stability = len(tract_seq) * homogeneity ** 2
                        length_bonus = len(tract_seq) ** 0.8
                        ph_factor = 1.5 if pyr_frac > 0.8 else 1.0
                        triplex_score = (base_stability + length_bonus) * ph_factor
                        
                        region = (start+1, end)
                        if region not in seen_regions:
                            motifs.append({
                                "Class": "Triplex_DNA",
                                "Subclass": f"Homo{tract_type.split('homo')[1]}_Tract",
                                "Start": start + 1, "End": end, "Length": end - start,
                                "Sequence": wrap(tract_seq), "PurineFrac": round(pur_frac, 2),
                                "PyrimidineFrac": round(pyr_frac, 2), "Homogeneity": round(homogeneity, 3),
                                "Score": float(triplex_score), "ScoreMethod": "Homogeneous_Tract"
                            })
                            seen_regions.add(region)
    
    except Exception:
        pass  # Continue with Python regex approach
    
    # Python regex for mirror repeats (back-references required)
    python_motifs = find_hdna_python_mirror_repeats(seqU)
    motifs.extend(python_motifs)
    
    return motifs

def find_hdna_python_mirror_repeats(seq: str) -> list:
    """Python regex approach for mirror repeats that require back-references"""
    motifs = []
    seqU = seq.upper()
    n = len(seqU)
    
    # Mirror repeat search: scan all rep_len (10–100), all spacers (0–8)
    seen_regions = set()
    for rep_len in range(10, min(101, n//2)):
        for spacer in range(0, 9):
            pattern = fr"([ACGT]{{{rep_len}}})[ACGT]{{{spacer}}}\1"
            for match in re.finditer(pattern, seqU):
                start = match.start()
                end = match.end()
                repeat = match.group(1)
                full_seq = seqU[start:end]
                pur_frac, pyr_frac = purine_fraction(full_seq), pyrimidine_fraction(full_seq)
                is_triplex = (pur_frac >= 0.9 or pyr_frac >= 0.9)
                homogeneity = max(pur_frac, pyr_frac)
                
                # Enhanced triplex scoring based on Frank-Kamenetskii & Mirkin (1995)
                base_stability = len(full_seq) * homogeneity ** 2
                repeat_bonus = rep_len * 0.8
                spacer_penalty = spacer * 2.0
                ph_factor = 1.5 if pyr_frac > 0.8 else 1.0
                triplex_score = (base_stability + repeat_bonus - spacer_penalty) * ph_factor
                
                region = (start+1, end)
                if region not in seen_regions:
                    motifs.append({
                        "Class": "Triplex_DNA" if is_triplex else "Mirror_Repeat",
                        "Subclass": "Mirror_Repeat_Triplex" if is_triplex else "Mirror_Repeat",
                        "Start": start + 1, "End": end, "Length": len(full_seq), "Spacer": spacer,
                        "Sequence": wrap(full_seq), "PurineFrac": round(pur_frac, 2), 
                        "PyrimidineFrac": round(pyr_frac, 2), "Homogeneity": round(homogeneity, 3),
                        "Score": float(triplex_score), "ScoreMethod": "Mirror_Repeat_Enhanced"
                    })
                    seen_regions.add(region)
    return motifs

def find_hdna(seq: str) -> list:
    """Find H-DNA triplex motifs using hybrid Hyperscan + Python approach"""
    return find_hdna_hyperscan(seq)

#--- Sticky DNA finder: GAA/TTC repeats (per Sakamoto 1999), Hyperscan block scan ---
def find_sticky_dna(seq: str) -> list:
    motifs=[]; seqU=seq.replace('\n','').replace(' ','').upper(); db=hyperscan.Database()
    pats=[(r"(?:GAA){59,}",1),(r"(?:TTC){59,}",2)]
    exprs=[p[0].encode() for p in pats]; ids=[p[1] for p in pats]
    def cb(id, start, end, flags, ctx):
        motif_seq=seqU[start:end]; repeat_len=len(motif_seq); repeat_count=repeat_len//3
        at_frac=(motif_seq.count('A')+motif_seq.count('T'))/repeat_len
        
        # Enhanced Sticky DNA scoring based on Sakamoto et al. (1999)
        # GAA/TTC repeats form very stable DNA triplexes
        base_triplex_potential = repeat_count ** 1.2  # Superlinear scaling
        at_stability = (1.0 + at_frac * 0.8)  # AT content affects stability
        length_bonus = repeat_len ** 0.6  # Sublinear length scaling
        sticky_score = base_triplex_potential * at_stability * length_bonus
        
        motifs.append({
            "Class":"Triplex_DNA","Subclass":"Sticky_DNA","Start":start+1,"End":end,
            "Length":repeat_len,"RepeatCount":repeat_count,"Sequence":wrap(motif_seq),
            "AT_Content": round(at_frac, 3), "ScoreMethod":"Sticky_enhanced","Score":float(sticky_score)
        }); return hyperscan.HS_SUCCESS
    db.compile(expressions=exprs, ids=ids)
    db.scan(seqU.encode(), match_event_handler=cb, context=None)
    return motifs

#--- Main entry: all triplex DNA motifs, standardized as per scientific conventions ---
def find_triplex(seq: str, sequence_name: str = "") -> list:
    hdna_results=find_hdna(seq); sticky_results=find_sticky_dna(seq)
    all_results=hdna_results+sticky_results
    return [standardize_motif_output(m,sequence_name,i) for i,m in enumerate(all_results,1)]

#--- Annotations ---
# - find_hdna: mirror repeats, all sizes/spacers, score/purity per literature (Wells 2007), Hyperscan block scan.
# - find_sticky_dna: GAA/TTC long repeats, Sakamoto 1999 threshold, score A/T bias, Hyperscan block scan.
# - find_triplex: combines both, standardized output for downstream genomic analysis.
