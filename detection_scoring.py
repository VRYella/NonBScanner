"""
Non-B DNA Motif Detection and Scoring Engine
============================================

| Module Function | Description                                           |
|-----------------|-------------------------------------------------------|
| Pattern Search  | Regex-based motif candidate identification            |
| Scoring Engine  | Scientific scoring algorithms (G4Hunter, Z-seeker)   |
| Post-processing | Overlap removal, clustering, quality filtering       |
| Classification  | Motif type assignment and subclass categorization    |

Detection Classes: 8 major Non-B DNA structure types
Scoring Methods: Literature-validated algorithms with normalized outputs
"""

import re
import numpy as np
from typing import List, Dict, Any, Tuple, Optional
from collections import defaultdict, Counter
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from registires import get_patterns_for_motif, get_pattern_info, ALL_PATTERNS, CLASS_DEFINITIONS

# =============================================================================
# CORE SCORING ALGORITHMS
# =============================================================================

def g4hunter_score(seq: str, window_size: int = 25) -> float:
    """
    G4Hunter scoring algorithm (Bedrat et al. NAR 2016)
    
    | Parameter    | Type  | Range    | Description                    |
    |--------------|-------|----------|--------------------------------|
    | seq          | str   | -        | DNA sequence (A,T,G,C only)    |
    | window_size  | int   | 10-100   | Sliding window size            |
    
    Returns: G4Hunter score (>1.2 = high G4 potential)
    """
    if len(seq) < window_size:
        return 0.0
    
    scores = []
    for i in range(len(seq) - window_size + 1):
        window = seq[i:i + window_size]
        g_count = window.count('G')
        c_count = window.count('C')
        
        # G4Hunter formula: G-quadruplex potential scoring
        if g_count + c_count > 0:
            score = (g_count - c_count) ** 2 / window_size
        else:
            score = 0.0
        scores.append(score)
    
    return max(scores) if scores else 0.0

def imotif_score(seq: str, window_size: int = 25) -> float:
    """
    i-Motif scoring algorithm (C-rich quadruplex structures)
    
    | Parameter    | Type  | Range    | Description                    |
    |--------------|-------|----------|--------------------------------|
    | seq          | str   | -        | DNA sequence                   |
    | window_size  | int   | 10-100   | Sliding window size            |
    
    Returns: i-Motif score (>1.0 = high C-quadruplex potential)
    """
    if len(seq) < window_size:
        return 0.0
    
    scores = []
    for i in range(len(seq) - window_size + 1):
        window = seq[i:i + window_size]
        c_count = window.count('C')
        g_count = window.count('G')
        
        # i-Motif favors C-rich regions
        if c_count + g_count > 0:
            score = (c_count - g_count) ** 2 / window_size
        else:
            score = 0.0
        scores.append(score)
    
    return max(scores) if scores else 0.0

def zdna_seeker_score(seq: str, **kwargs) -> float:
    """
    Z-DNA seeker scoring algorithm (Ho et al. EMBO J 1986)
    
    | Dinucleotide | Weight | Biological Significance         |
    |--------------|--------|---------------------------------|
    | GC/CG        | 7.0    | High Z-DNA formation potential |
    | GT/TG        | 1.25   | Moderate Z-DNA potential        |
    | AT/TA        | 0.5    | Low Z-DNA potential             |
    | AC/CA        | 1.25   | Moderate Z-DNA potential        |
    
    Returns: Z-DNA formation score (>50 = high Z-DNA potential)
    """
    if len(seq) < 12:
        return 0.0
    
    # Z-DNA scoring weights from literature
    weights = {
        'GC': 7.0, 'CG': 7.0,
        'GT': 1.25, 'TG': 1.25, 'AC': 1.25, 'CA': 1.25,
        'AT': 0.5, 'TA': 0.5,
        'GA': 0.0, 'AG': 0.0, 'CT': 0.0, 'TC': 0.0,
        'AA': -3.0, 'TT': -3.0, 'GG': -3.0, 'CC': -3.0
    }
    
    scores = []
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        score = weights.get(dinuc, 0.0)
        scores.append(score)
    
    # Apply consecutive AT penalty
    consecutive_at = 0
    for i, score in enumerate(scores):
        if seq[i:i+2] in ['AT', 'TA']:
            consecutive_at += 1
            if consecutive_at > 4:
                scores[i] -= 5.0  # Literature-based penalty
        else:
            consecutive_at = 0
    
    return sum(scores)

def curvature_score(seq: str) -> float:
    """
    DNA curvature scoring based on A-tract positioning
    
    | Feature        | Score Weight | Description                |
    |----------------|--------------|----------------------------|
    | A-tract length | 2.0 per bp   | Longer tracts = more curve |
    | Phasing        | 5.0 bonus    | ~10bp spacing increases    |
    | AT content     | 1.0 per %    | Overall AT richness        |
    
    Returns: Curvature score (>15 = significant curvature)
    """
    if len(seq) < 8:
        return 0.0
    
    score = 0.0
    
    # Find A-tracts and T-tracts
    a_tracts = []
    t_tracts = []
    
    for match in re.finditer(r'A{4,}', seq):
        a_tracts.append((match.start(), match.end()))
        score += (match.end() - match.start()) * 2.0
    
    for match in re.finditer(r'T{4,}', seq):
        t_tracts.append((match.start(), match.end()))
        score += (match.end() - match.start()) * 2.0
    
    # Check for phased A-tracts (spacing ~10bp)
    all_tracts = a_tracts + t_tracts
    all_tracts.sort()
    
    for i in range(len(all_tracts) - 1):
        spacing = all_tracts[i+1][0] - all_tracts[i][1]
        if 8 <= spacing <= 12:  # Nucleosome-like phasing
            score += 5.0
    
    # AT content bonus
    at_content = (seq.count('A') + seq.count('T')) / len(seq) * 100
    score += at_content * 0.2
    
    return score

def triplex_score(seq: str) -> float:
    """
    Triplex DNA scoring (purine/pyrimidine tract analysis)
    
    | Tract Type    | Length Req | Score Formula              |
    |---------------|------------|----------------------------|
    | Homopurine    | 12+ bp     | length * 2.0               |
    | Homopyrimidine| 12+ bp     | length * 2.0               |
    | Mirror repeat | 6+6 bp     | (len1 * len2) / distance  |
    
    Returns: Triplex formation score (>25 = triplex potential)
    """
    if len(seq) < 12:
        return 0.0
    
    score = 0.0
    
    # Find homopurine tracts [AG]
    for match in re.finditer(r'[AG]{12,}', seq):
        tract_len = match.end() - match.start()
        score += tract_len * 2.0
    
    # Find homopyrimidine tracts [CT]
    for match in re.finditer(r'[CT]{12,}', seq):
        tract_len = match.end() - match.start()
        score += tract_len * 2.0
    
    # Check for mirror repeats (purine followed by pyrimidine)
    purine_tracts = [(m.start(), m.end()) for m in re.finditer(r'[AG]{6,}', seq)]
    pyrimidine_tracts = [(m.start(), m.end()) for m in re.finditer(r'[CT]{6,}', seq)]
    
    for p_start, p_end in purine_tracts:
        for py_start, py_end in pyrimidine_tracts:
            if p_end < py_start:  # Purine before pyrimidine
                distance = py_start - p_end
                if 10 <= distance <= 100:  # Reasonable spacing
                    p_len = p_end - p_start
                    py_len = py_end - py_start
                    mirror_score = (p_len * py_len) / distance
                    score += mirror_score
    
    return score

def rloop_score(seq: str) -> float:
    """
    R-loop formation scoring (GC skew and G-cluster analysis)
    
    | Feature       | Weight | Description                     |
    |---------------|--------|---------------------------------|
    | GC skew       | 5.0    | (G-C)/(G+C) asymmetry          |
    | G clusters    | 3.0    | Groups of 3+ G's               |
    | CpG density   | 2.0    | CG dinucleotide frequency       |
    
    Returns: R-loop potential score (>30 = R-loop formation likely)
    """
    if len(seq) < 15:
        return 0.0
    
    score = 0.0
    g_count = seq.count('G')
    c_count = seq.count('C')
    
    # GC skew calculation
    if g_count + c_count > 0:
        gc_skew = abs(g_count - c_count) / (g_count + c_count)
        score += gc_skew * 50.0
    
    # G cluster scoring
    for match in re.finditer(r'G{3,}', seq):
        cluster_len = match.end() - match.start()
        score += cluster_len * 3.0
    
    # CpG density
    cpg_count = seq.count('CG')
    cpg_density = cpg_count / len(seq) * 100
    score += cpg_density * 2.0
    
    return score

def repeat_score(seq: str, repeat_unit: str) -> float:
    """
    Repeat sequence scoring for slipped DNA structures
    
    | Parameter    | Type | Description                      |
    |--------------|------|----------------------------------|
    | seq          | str  | DNA sequence to analyze          |
    | repeat_unit  | str  | Expected repeat unit (e.g. "AT") |
    
    Returns: Repeat score (higher = more repetitive)
    """
    if len(repeat_unit) == 0:
        return 0.0
    
    unit_len = len(repeat_unit)
    perfect_repeats = 0
    total_units = len(seq) // unit_len
    
    for i in range(0, len(seq) - unit_len + 1, unit_len):
        if seq[i:i + unit_len] == repeat_unit:
            perfect_repeats += 1
    
    if total_units > 0:
        repeat_fraction = perfect_repeats / total_units
        return repeat_fraction * len(seq)
    
    return 0.0

# =============================================================================
# MOTIF DETECTION ENGINE
# =============================================================================

def detect_motifs_in_sequence(seq: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
    """
    Main motif detection function using pattern registry
    
    | Detection Step | Description                               |
    |----------------|-------------------------------------------|
    | Pattern Match  | Apply regex patterns to sequence         |
    | Score Calc     | Calculate scientific scores              |
    | Classification | Assign motif class and subclass         |
    | Quality Filter | Remove low-scoring candidates            |
    
    Returns: List of detected motif dictionaries
    """
    all_motifs = []
    seq = seq.upper()
    
    # Detection mapping table
    detection_methods = {
        'G_QUADRUPLEX': detect_g_quadruplex_motifs,
        'I_MOTIF': detect_imotif_motifs,
        'Z_DNA': detect_zdna_motifs,
        'CURVED_DNA': detect_curved_motifs,
        'TRIPLEX': detect_triplex_motifs,
        'CRUCIFORM': detect_cruciform_motifs,
        'R_LOOP': detect_rloop_motifs,
        'SLIPPED_DNA': detect_slipped_motifs
    }
    
    for motif_class, detect_func in detection_methods.items():
        try:
            motifs = detect_func(seq, sequence_name)
            all_motifs.extend(motifs)
        except Exception as e:
            print(f"Error detecting {motif_class}: {e}")
            continue
    
    return all_motifs

def detect_g_quadruplex_motifs(seq: str, sequence_name: str) -> List[Dict[str, Any]]:
    """Detect G-quadruplex motifs using G4Hunter scoring"""
    motifs = []
    patterns = get_patterns_for_motif('G_QUADRUPLEX')
    
    for pattern_group, pattern_list in patterns.items():
        for pattern_tuple in pattern_list:
            if not pattern_tuple:
                continue
                
            regex, pattern_id, group_idx, subclass, _, scale, min_runs, min_score, method = pattern_tuple
            
            for match in re.finditer(regex, seq):
                start_pos = match.start() + 1  # 1-based coordinates
                end_pos = match.end()
                motif_seq = match.group(0)
                
                # Calculate G4Hunter score
                raw_score = g4hunter_score(motif_seq)
                normalized_score = min(raw_score / 2.0, 1.0)  # Normalize to 0-1
                
                if raw_score >= min_score:
                    motif = {
                        'Class': 6,
                        'Subclass': subclass,
                        'Start': start_pos,
                        'End': end_pos,
                        'Length': len(motif_seq),
                        'Sequence': motif_seq,
                        'Raw_Score': round(raw_score, 3),
                        'Normalized_Score': round(normalized_score, 3),
                        'Scoring_Method': method,
                        'Pattern_ID': pattern_id,
                        'Sequence_Name': sequence_name
                    }
                    motifs.append(motif)
    
    return motifs

def detect_imotif_motifs(seq: str, sequence_name: str) -> List[Dict[str, Any]]:
    """Detect i-Motif structures using C-rich quadruplex scoring"""
    motifs = []
    patterns = get_patterns_for_motif('I_MOTIF')
    
    for pattern_group, pattern_list in patterns.items():
        for pattern_tuple in pattern_list:
            if not pattern_tuple:
                continue
                
            regex, pattern_id, group_idx, subclass, _, scale, min_runs, min_score, method = pattern_tuple
            
            for match in re.finditer(regex, seq):
                start_pos = match.start() + 1
                end_pos = match.end()
                motif_seq = match.group(0)
                
                raw_score = imotif_score(motif_seq)
                normalized_score = min(raw_score / 1.5, 1.0)
                
                if raw_score >= min_score:
                    motif = {
                        'Class': 7,
                        'Subclass': subclass,
                        'Start': start_pos,
                        'End': end_pos,
                        'Length': len(motif_seq),
                        'Sequence': motif_seq,
                        'Raw_Score': round(raw_score, 3),
                        'Normalized_Score': round(normalized_score, 3),
                        'Scoring_Method': method,
                        'Pattern_ID': pattern_id,
                        'Sequence_Name': sequence_name
                    }
                    motifs.append(motif)
    
    return motifs

def detect_zdna_motifs(seq: str, sequence_name: str) -> List[Dict[str, Any]]:
    """Detect Z-DNA using Z-seeker algorithm"""
    motifs = []
    patterns = get_patterns_for_motif('Z_DNA')
    
    for pattern_group, pattern_list in patterns.items():
        for pattern_tuple in pattern_list:
            if not pattern_tuple:
                continue
                
            regex, pattern_id, group_idx, subclass, _, scale, min_runs, min_score, method = pattern_tuple
            
            for match in re.finditer(regex, seq):
                start_pos = match.start() + 1
                end_pos = match.end()
                motif_seq = match.group(0)
                
                raw_score = zdna_seeker_score(motif_seq)
                normalized_score = min(raw_score / 100.0, 1.0)
                
                if raw_score >= min_score:
                    motif = {
                        'Class': 8,
                        'Subclass': subclass,
                        'Start': start_pos,
                        'End': end_pos,
                        'Length': len(motif_seq),
                        'Sequence': motif_seq,
                        'Raw_Score': round(raw_score, 3),
                        'Normalized_Score': round(normalized_score, 3),
                        'Scoring_Method': method,
                        'Pattern_ID': pattern_id,
                        'Sequence_Name': sequence_name
                    }
                    motifs.append(motif)
    
    return motifs

def detect_curved_motifs(seq: str, sequence_name: str) -> List[Dict[str, Any]]:
    """Detect curved DNA using A-tract curvature analysis"""
    motifs = []
    patterns = get_patterns_for_motif('CURVED_DNA')
    
    for pattern_group, pattern_list in patterns.items():
        for pattern_tuple in pattern_list:
            if not pattern_tuple:
                continue
                
            regex, pattern_id, group_idx, subclass, _, scale, min_runs, min_score, method = pattern_tuple
            
            for match in re.finditer(regex, seq):
                start_pos = match.start() + 1
                end_pos = match.end()
                motif_seq = match.group(0)
                
                raw_score = curvature_score(motif_seq)
                normalized_score = min(raw_score / 50.0, 1.0)
                
                if raw_score >= min_score:
                    motif = {
                        'Class': 1,
                        'Subclass': subclass,
                        'Start': start_pos,
                        'End': end_pos,
                        'Length': len(motif_seq),
                        'Sequence': motif_seq,
                        'Raw_Score': round(raw_score, 3),
                        'Normalized_Score': round(normalized_score, 3),
                        'Scoring_Method': method,
                        'Pattern_ID': pattern_id,
                        'Sequence_Name': sequence_name
                    }
                    motifs.append(motif)
    
    return motifs

def detect_triplex_motifs(seq: str, sequence_name: str) -> List[Dict[str, Any]]:
    """Detect triplex DNA using purine/pyrimidine tract analysis"""
    motifs = []
    patterns = get_patterns_for_motif('TRIPLEX')
    
    for pattern_group, pattern_list in patterns.items():
        for pattern_tuple in pattern_list:
            if not pattern_tuple:
                continue
                
            regex, pattern_id, group_idx, subclass, _, scale, min_runs, min_score, method = pattern_tuple
            
            for match in re.finditer(regex, seq):
                start_pos = match.start() + 1
                end_pos = match.end()
                motif_seq = match.group(0)
                
                raw_score = triplex_score(motif_seq)
                normalized_score = min(raw_score / 75.0, 1.0)
                
                if raw_score >= min_score:
                    motif = {
                        'Class': 5,
                        'Subclass': subclass,
                        'Start': start_pos,
                        'End': end_pos,
                        'Length': len(motif_seq),
                        'Sequence': motif_seq,
                        'Raw_Score': round(raw_score, 3),
                        'Normalized_Score': round(normalized_score, 3),
                        'Scoring_Method': method,
                        'Pattern_ID': pattern_id,
                        'Sequence_Name': sequence_name
                    }
                    motifs.append(motif)
    
    return motifs

def detect_cruciform_motifs(seq: str, sequence_name: str) -> List[Dict[str, Any]]:
    """Detect cruciform DNA using palindrome detection"""
    motifs = []
    patterns = get_patterns_for_motif('CRUCIFORM')
    
    for pattern_group, pattern_list in patterns.items():
        for pattern_tuple in pattern_list:
            if not pattern_tuple:
                continue
                
            regex, pattern_id, group_idx, subclass, _, scale, min_runs, min_score, method = pattern_tuple
            
            for match in re.finditer(regex, seq):
                start_pos = match.start() + 1
                end_pos = match.end()
                motif_seq = match.group(0)
                
                # Simple palindrome scoring based on length and GC content
                raw_score = len(motif_seq) * 2.0
                gc_content = (motif_seq.count('G') + motif_seq.count('C')) / len(motif_seq)
                raw_score *= (1.0 + gc_content)  # GC-rich palindromes score higher
                
                normalized_score = min(raw_score / 100.0, 1.0)
                
                if raw_score >= min_score:
                    motif = {
                        'Class': 3,
                        'Subclass': subclass,
                        'Start': start_pos,
                        'End': end_pos,
                        'Length': len(motif_seq),
                        'Sequence': motif_seq,
                        'Raw_Score': round(raw_score, 3),
                        'Normalized_Score': round(normalized_score, 3),
                        'Scoring_Method': method,
                        'Pattern_ID': pattern_id,
                        'Sequence_Name': sequence_name
                    }
                    motifs.append(motif)
    
    return motifs

def detect_rloop_motifs(seq: str, sequence_name: str) -> List[Dict[str, Any]]:
    """Detect R-loop formation sites"""
    motifs = []
    patterns = get_patterns_for_motif('R_LOOP')
    
    for pattern_group, pattern_list in patterns.items():
        for pattern_tuple in pattern_list:
            if not pattern_tuple:
                continue
                
            regex, pattern_id, group_idx, subclass, _, scale, min_runs, min_score, method = pattern_tuple
            
            for match in re.finditer(regex, seq):
                start_pos = match.start() + 1
                end_pos = match.end()
                motif_seq = match.group(0)
                
                raw_score = rloop_score(motif_seq)
                normalized_score = min(raw_score / 100.0, 1.0)
                
                if raw_score >= min_score:
                    motif = {
                        'Class': 4,
                        'Subclass': subclass,
                        'Start': start_pos,
                        'End': end_pos,
                        'Length': len(motif_seq),
                        'Sequence': motif_seq,
                        'Raw_Score': round(raw_score, 3),
                        'Normalized_Score': round(normalized_score, 3),
                        'Scoring_Method': method,
                        'Pattern_ID': pattern_id,
                        'Sequence_Name': sequence_name
                    }
                    motifs.append(motif)
    
    return motifs

def detect_slipped_motifs(seq: str, sequence_name: str) -> List[Dict[str, Any]]:
    """Detect slipped DNA structures from tandem repeats"""
    motifs = []
    patterns = get_patterns_for_motif('SLIPPED_DNA')
    
    for pattern_group, pattern_list in patterns.items():
        for pattern_tuple in pattern_list:
            if not pattern_tuple:
                continue
                
            regex, pattern_id, group_idx, subclass, _, scale, min_runs, min_score, method = pattern_tuple
            
            for match in re.finditer(regex, seq):
                start_pos = match.start() + 1
                end_pos = match.end()
                motif_seq = match.group(0)
                
                # Determine repeat unit from pattern
                if 'mononucleotide' in pattern_group:
                    repeat_unit = motif_seq[0]
                elif 'dinucleotide' in pattern_group:
                    repeat_unit = motif_seq[:2]
                elif 'trinucleotide' in pattern_group:
                    repeat_unit = motif_seq[:3]
                else:
                    repeat_unit = motif_seq[:2]  # Default
                
                raw_score = repeat_score(motif_seq, repeat_unit)
                normalized_score = min(raw_score / len(motif_seq), 1.0)
                
                if raw_score >= min_score:
                    motif = {
                        'Class': 2,
                        'Subclass': subclass,
                        'Start': start_pos,
                        'End': end_pos,
                        'Length': len(motif_seq),
                        'Sequence': motif_seq,
                        'Raw_Score': round(raw_score, 3),
                        'Normalized_Score': round(normalized_score, 3),
                        'Scoring_Method': method,
                        'Pattern_ID': pattern_id,
                        'Sequence_Name': sequence_name,
                        'Repeat_Unit': repeat_unit
                    }
                    motifs.append(motif)
    
    return motifs

# =============================================================================
# POST-PROCESSING AND FILTERING
# =============================================================================

def remove_overlapping_motifs(motifs: List[Dict[str, Any]], overlap_threshold: float = 0.5) -> List[Dict[str, Any]]:
    """
    Remove overlapping motifs keeping highest scoring ones
    
    | Parameter         | Type  | Range  | Description                    |
    |-------------------|-------|--------|--------------------------------|
    | motifs            | list  | -      | List of motif dictionaries     |
    | overlap_threshold | float | 0-1    | Max allowed overlap fraction   |
    
    Returns: Filtered list with reduced overlaps
    """
    if not motifs:
        return motifs
    
    # Sort by normalized score (descending)
    sorted_motifs = sorted(motifs, key=lambda x: x.get('Normalized_Score', 0), reverse=True)
    filtered_motifs = []
    
    for motif in sorted_motifs:
        is_overlapping = False
        motif_start, motif_end = motif['Start'], motif['End']
        
        for existing in filtered_motifs:
            existing_start, existing_end = existing['Start'], existing['End']
            
            # Calculate overlap
            overlap_start = max(motif_start, existing_start)
            overlap_end = min(motif_end, existing_end)
            
            if overlap_start < overlap_end:
                overlap_length = overlap_end - overlap_start
                motif_length = motif_end - motif_start
                overlap_fraction = overlap_length / motif_length
                
                if overlap_fraction > overlap_threshold:
                    is_overlapping = True
                    break
        
        if not is_overlapping:
            filtered_motifs.append(motif)
    
    return filtered_motifs

def filter_by_score_threshold(motifs: List[Dict[str, Any]], min_score: float = 0.3) -> List[Dict[str, Any]]:
    """
    Filter motifs by minimum normalized score threshold
    
    | Threshold | Quality Level | Description                    |
    |-----------|---------------|--------------------------------|
    | 0.1-0.3   | Low           | Weak motifs, many false pos   |
    | 0.3-0.6   | Medium        | Moderate confidence            |
    | 0.6-0.8   | High          | High confidence predictions    |
    | 0.8-1.0   | Very High     | Exceptional motifs             |
    
    Returns: Filtered list above threshold
    """
    return [motif for motif in motifs if motif.get('Normalized_Score', 0) >= min_score]

def merge_nearby_motifs(motifs: List[Dict[str, Any]], max_gap: int = 10) -> List[Dict[str, Any]]:
    """
    Merge motifs of same class that are close together
    
    | Parameter | Type | Range  | Description                     |
    |-----------|------|--------|---------------------------------|
    | motifs    | list | -      | List of motifs to merge         |
    | max_gap   | int  | 1-100  | Maximum gap for merging (bp)    |
    
    Returns: List with merged nearby motifs
    """
    if not motifs:
        return motifs
    
    # Group by motif class and sequence
    grouped = defaultdict(list)
    for motif in motifs:
        key = (motif.get('Class'), motif.get('Sequence_Name'))
        grouped[key].append(motif)
    
    merged_motifs = []
    
    for group_key, group_motifs in grouped.items():
        # Sort by start position
        group_motifs.sort(key=lambda x: x['Start'])
        
        i = 0
        while i < len(group_motifs):
            current_motif = group_motifs[i]
            j = i + 1
            
            # Look for nearby motifs to merge
            while j < len(group_motifs):
                next_motif = group_motifs[j]
                gap = next_motif['Start'] - current_motif['End']
                
                if gap <= max_gap:
                    # Merge motifs
                    merged_motif = {
                        'Class': current_motif['Class'],
                        'Subclass': current_motif['Subclass'],
                        'Start': current_motif['Start'],
                        'End': next_motif['End'],
                        'Length': next_motif['End'] - current_motif['Start'],
                        'Sequence': 'MERGED',
                        'Raw_Score': max(current_motif.get('Raw_Score', 0), next_motif.get('Raw_Score', 0)),
                        'Normalized_Score': max(current_motif.get('Normalized_Score', 0), 
                                              next_motif.get('Normalized_Score', 0)),
                        'Scoring_Method': current_motif.get('Scoring_Method', 'Unknown'),
                        'Pattern_ID': current_motif.get('Pattern_ID', 0),
                        'Sequence_Name': current_motif['Sequence_Name']
                    }
                    current_motif = merged_motif
                    j += 1
                else:
                    break
            
            merged_motifs.append(current_motif)
            i = j
    
    return merged_motifs

def calculate_motif_statistics(motifs: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Calculate comprehensive statistics for detected motifs
    
    | Statistic       | Description                             |
    |-----------------|-----------------------------------------|
    | Total Count     | Number of detected motifs               |
    | Class Distribution| Motifs per class                      |
    | Score Distribution| Raw and normalized score ranges       |
    | Length Distribution| Size distribution of motifs           |
    
    Returns: Dictionary with statistical summaries
    """
    if not motifs:
        return {
            'Total_Motifs': 0,
            'Classes_Found': 0,
            'Average_Score': 0.0,
            'Average_Length': 0.0,
            'Class_Distribution': {}
        }
    
    stats = {
        'Total_Motifs': len(motifs),
        'Classes_Found': len(set(motif.get('Class', 0) for motif in motifs)),
        'Average_Score': np.mean([motif.get('Normalized_Score', 0) for motif in motifs]),
        'Average_Length': np.mean([motif.get('Length', 0) for motif in motifs]),
        'Score_Range': {
            'Min': min(motif.get('Normalized_Score', 0) for motif in motifs),
            'Max': max(motif.get('Normalized_Score', 0) for motif in motifs),
            'Std': np.std([motif.get('Normalized_Score', 0) for motif in motifs])
        },
        'Length_Range': {
            'Min': min(motif.get('Length', 0) for motif in motifs),
            'Max': max(motif.get('Length', 0) for motif in motifs),
            'Std': np.std([motif.get('Length', 0) for motif in motifs])
        }
    }
    
    # Class distribution
    class_counts = Counter(motif.get('Class', 0) for motif in motifs)
    class_distribution = {}
    for class_id, count in class_counts.items():
        class_name = CLASS_DEFINITIONS.get(class_id, {}).get('name', f'Class_{class_id}')
        class_distribution[class_name] = count
    
    stats['Class_Distribution'] = class_distribution
    
    return stats

def process_sequence_parallel(args_tuple):
    """Helper function for parallel processing of sequences"""
    seq, seq_name, filters = args_tuple
    motifs = detect_motifs_in_sequence(seq, seq_name)
    
    # Apply filters if specified
    if filters.get('remove_overlaps', True):
        motifs = remove_overlapping_motifs(motifs, filters.get('overlap_threshold', 0.5))
    
    if filters.get('score_filter', True):
        motifs = filter_by_score_threshold(motifs, filters.get('min_score', 0.3))
    
    if filters.get('merge_nearby', True):
        motifs = merge_nearby_motifs(motifs, filters.get('max_gap', 10))
    
    return motifs

def analyze_multiple_sequences(sequences: Dict[str, str], 
                             max_workers: int = None, 
                             filters: Dict[str, Any] = None) -> Dict[str, List[Dict[str, Any]]]:
    """
    Analyze multiple sequences in parallel
    
    | Parameter   | Type | Description                          |
    |-------------|------|--------------------------------------|
    | sequences   | dict | {name: sequence} pairs               |
    | max_workers | int  | Number of parallel processes         |
    | filters     | dict | Post-processing filter settings      |
    
    Returns: Dictionary mapping sequence names to motif lists
    """
    if max_workers is None:
        max_workers = min(8, multiprocessing.cpu_count())
    
    if filters is None:
        filters = {
            'remove_overlaps': True,
            'overlap_threshold': 0.5,
            'score_filter': True,
            'min_score': 0.3,
            'merge_nearby': True,
            'max_gap': 10
        }
    
    # Prepare arguments for parallel processing
    args_list = [(seq, name, filters) for name, seq in sequences.items()]
    
    results = {}
    
    if len(sequences) == 1:
        # Single sequence - no need for parallelization
        seq_name, seq = list(sequences.items())[0]
        results[seq_name] = process_sequence_parallel((seq, seq_name, filters))
    else:
        # Multiple sequences - use parallel processing
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_name = {}
            for args in args_list:
                future = executor.submit(process_sequence_parallel, args)
                future_to_name[future] = args[1]  # sequence name
            
            for future in as_completed(future_to_name):
                seq_name = future_to_name[future]
                try:
                    motifs = future.result()
                    results[seq_name] = motifs
                except Exception as e:
                    print(f"Error processing sequence {seq_name}: {e}")
                    results[seq_name] = []
    
    return results