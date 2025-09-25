"""
Central Registry for Non-B DNA Pattern Recognition
==================================================

| Class | Motif Type    | Patterns | Description                           |
|-------|---------------|----------|---------------------------------------|
| 1     | Curved DNA    | 15       | A-tract mediated DNA curvature        |
| 2     | Slipped DNA   | 37       | Direct/tandem repeats                 |
| 3     | Cruciform     | 89       | Inverted repeats, four-way junctions  |
| 4     | R-loop        | 11       | RNA-DNA hybrids                       |
| 5     | Triplex       | 14       | Three-stranded DNA structures         |
| 6     | G-Quadruplex  | 16       | Four-stranded G-rich structures       |
| 7     | i-Motif       | 12       | C-rich quadruplex structures          |
| 8     | Z-DNA         | 13       | Left-handed double helix              |

Pattern Format: (regex, id, group, subclass, scorer, scale, runs, score, method)
"""

import re
import numpy as np

# =============================================================================
# G-QUADRUPLEX PATTERNS (Class 6) - G4Hunter Compatible
# =============================================================================
G_QUADRUPLEX_PATTERNS = {
    'canonical_g4': [
        (r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}", 1, 0, "Canonical_G4", None, 1.0, 4, 1.2, "G4Hunter"),
        (r"G{4,}[ACGT]{1,7}G{4,}[ACGT]{1,7}G{4,}[ACGT]{1,7}G{4,}", 2, 0, "Perfect_G4", None, 1.2, 4, 1.5, "G4Hunter"),
    ],
    'relaxed_g4': [
        (r"G{2,}[ACGT]{1,12}G{2,}[ACGT]{1,12}G{2,}[ACGT]{1,12}G{2,}", 3, 0, "Relaxed_G4", None, 0.8, 4, 0.8, "G4Hunter"),
        (r"G{3,}[ACGT]{8,12}G{3,}[ACGT]{8,12}G{3,}[ACGT]{8,12}G{3,}", 4, 0, "Long_Loop_G4", None, 0.8, 4, 0.8, "G4Hunter"),
    ],
    'bulged_g4': [
        (r"G{3,}[ACGT]{0,3}G{3,}[ACGT]{0,3}G{3,}[ACGT]{0,3}G{3,}", 5, 0, "Bulged_G4", None, 0.7, 4, 0.7, "G4Hunter"),
        (r"G{2,}[ACGT]{0,3}G{3,}[ACGT]{0,3}G{3,}[ACGT]{0,3}G{3,}", 6, 0, "Imperfect_Bulged_G4", None, 0.7, 4, 0.5, "G4Hunter"),
    ],
    'imperfect_g4': [
        (r"G{2,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}", 8, 0, "Imperfect_G4", None, 0.9, 4, 0.8, "G4Hunter"),
        (r"G{3,}[ACGT]{1,7}G{2,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}", 9, 0, "Imperfect_G4", None, 0.9, 4, 0.8, "G4Hunter"),
    ]
}

# =============================================================================
# I-MOTIF PATTERNS (Class 7) - pH-dependent C-rich structures  
# =============================================================================
I_MOTIF_PATTERNS = {
    'canonical_imotif': [
        (r"C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}", 20, 0, "Canonical_iMotif", None, 1.0, 4, 1.0, "iM_Hunter"),
        (r"C{4,}[ACGT]{1,7}C{4,}[ACGT]{1,7}C{4,}[ACGT]{1,7}C{4,}", 21, 0, "Perfect_iMotif", None, 1.2, 4, 1.2, "iM_Hunter"),
    ],
    'relaxed_imotif': [
        (r"C{2,}[ACGT]{1,12}C{2,}[ACGT]{1,12}C{2,}[ACGT]{1,12}C{2,}", 22, 0, "Relaxed_iMotif", None, 0.8, 4, 0.8, "iM_Hunter"),
        (r"C{3,}[ACGT]{8,12}C{3,}[ACGT]{8,12}C{3,}[ACGT]{8,12}C{3,}", 23, 0, "Long_Loop_iMotif", None, 0.8, 4, 0.8, "iM_Hunter"),
    ],
    'ac_motif': [
        (r"A{2,}C{2,}[ACGT]{1,7}A{2,}C{2,}[ACGT]{1,7}A{2,}C{2,}", 24, 0, "AC_Motif", None, 0.9, 3, 0.7, "AC_Score"),
        (r"C{2,}A{2,}[ACGT]{1,7}C{2,}A{2,}[ACGT]{1,7}C{2,}A{2,}", 25, 0, "CA_Motif", None, 0.9, 3, 0.7, "AC_Score"),
    ],
    'imperfect_imotif': [
        (r"C{2,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}", 26, 0, "Imperfect_iMotif", None, 0.9, 4, 0.8, "iM_Hunter"),
    ]
}

# =============================================================================
# Z-DNA PATTERNS (Class 8) - Left-handed helix structures
# =============================================================================
Z_DNA_PATTERNS = {
    'cg_zdna': [
        (r"(?:CG){6,}", 30, 0, "CG_repeat", None, 1.0, 6, 50.0, "ZSeeker"),
        (r"(?:GC){6,}", 31, 0, "GC_repeat", None, 1.0, 6, 50.0, "ZSeeker"),
        (r"(?:CG|GC){6,}", 32, 0, "Mixed_CG", None, 1.0, 6, 50.0, "ZSeeker"),
    ],
    'at_zdna': [
        (r"(?:AT){8,}", 33, 0, "AT_repeat", None, 0.7, 8, 30.0, "ZSeeker"),
        (r"(?:TA){8,}", 34, 0, "TA_repeat", None, 0.7, 8, 30.0, "ZSeeker"),
    ],
    'egz_dna': [
        (r"(?:CGG){4,}", 35, 0, "eGZ_CGG", None, 1.2, 4, 60.0, "eGZ_Score"),
        (r"(?:CCG){4,}", 36, 0, "eGZ_CCG", None, 1.2, 4, 60.0, "eGZ_Score"),
        (r"(?:GGC){4,}", 37, 0, "eGZ_GGC", None, 1.2, 4, 60.0, "eGZ_Score"),
    ]
}

# =============================================================================
# CURVED DNA PATTERNS (Class 1) - A-tract mediated curvature
# =============================================================================
CURVED_DNA_PATTERNS = {
    'a_tracts': [
        (r"A{4,}[^A]{1,6}A{4,}", 40, 0, "A_tract_spaced", None, 1.0, 2, 15.0, "Curvature"),
        (r"(?:A{4,}[^A]{1,6}){3,}A{4,}", 41, 0, "Phased_A_tracts", None, 1.2, 4, 20.0, "Curvature"),
        (r"A{6,}", 42, 0, "Long_A_tract", None, 0.8, 1, 12.0, "Curvature"),
    ],
    't_tracts': [
        (r"T{4,}[^T]{1,6}T{4,}", 43, 0, "T_tract_spaced", None, 1.0, 2, 15.0, "Curvature"),
        (r"(?:T{4,}[^T]{1,6}){3,}T{4,}", 44, 0, "Phased_T_tracts", None, 1.2, 4, 20.0, "Curvature"),
        (r"T{6,}", 45, 0, "Long_T_tract", None, 0.8, 1, 12.0, "Curvature"),
    ],
    'mixed_at_tracts': [
        (r"(?:[AT]{4,}[^AT]{1,6}){3,}[AT]{4,}", 46, 0, "Mixed_AT_tracts", None, 1.1, 4, 18.0, "Curvature"),
        (r"(?:A{3,}T{3,}|T{3,}A{3,})", 47, 0, "AT_junction", None, 0.9, 1, 10.0, "Curvature"),
    ]
}

# =============================================================================
# TRIPLEX PATTERNS (Class 5) - Three-stranded structures
# =============================================================================
TRIPLEX_PATTERNS = {
    'homopurine_tracts': [
        (r"(?:AG|GA){6,}", 50, 0, "Homopurine_AG", None, 1.0, 6, 25.0, "Triplex_Score"),
        (r"(?:AAG|AGA|GAA){4,}", 51, 0, "Homopurine_AAG", None, 1.0, 4, 25.0, "Triplex_Score"),
        (r"[AG]{12,}", 52, 0, "Homopurine_tract", None, 1.1, 1, 30.0, "Triplex_Score"),
    ],
    'homopyrimidine_tracts': [
        (r"(?:CT|TC){6,}", 53, 0, "Homopyrimidine_CT", None, 1.0, 6, 25.0, "Triplex_Score"),
        (r"(?:TTC|TCT|CTT){4,}", 54, 0, "Homopyrimidine_TTC", None, 1.0, 4, 25.0, "Triplex_Score"),
        (r"[CT]{12,}", 55, 0, "Homopyrimidine_tract", None, 1.1, 1, 30.0, "Triplex_Score"),
    ],
    'mirror_repeats': [
        (r"([AG]{6,})[ACGT]{10,50}(?=.*([CT]{6,}))", 56, 1, "Mirror_repeat", None, 1.2, 2, 35.0, "Mirror_Score"),
        (r"([CT]{6,})[ACGT]{10,50}(?=.*([AG]{6,}))", 57, 1, "Mirror_repeat_rev", None, 1.2, 2, 35.0, "Mirror_Score"),
    ]
}

# =============================================================================
# CRUCIFORM PATTERNS (Class 3) - Inverted repeats and palindromes
# =============================================================================
CRUCIFORM_PATTERNS = {
    'perfect_palindromes': [
        (r"([ACGT]{6,})[ACGT]{0,10}\1", 60, 0, "Perfect_palindrome", None, 1.5, 1, 40.0, "Palindrome"),
        (r"([ACGT]{8,})[ACGT]{0,10}\1", 61, 0, "Long_palindrome", None, 1.8, 1, 50.0, "Palindrome"),
        (r"([ACGT]{4,6})[ACGT]{0,8}\1", 62, 0, "Short_palindrome", None, 1.2, 1, 25.0, "Palindrome"),
    ],
    'at_rich_palindromes': [
        (r"([AT]{6,})[ACGT]{0,10}\1", 63, 0, "AT_palindrome", None, 1.3, 1, 35.0, "AT_Palindrome"),
        (r"([AT]{4,})[ACGT]{0,8}\1", 64, 0, "Short_AT_palindrome", None, 1.0, 1, 20.0, "AT_Palindrome"),
    ],
    'gc_rich_palindromes': [
        (r"([GC]{6,})[ACGT]{0,10}\1", 65, 0, "GC_palindrome", None, 1.4, 1, 38.0, "GC_Palindrome"),
        (r"([GC]{4,})[ACGT]{0,8}\1", 66, 0, "Short_GC_palindrome", None, 1.1, 1, 22.0, "GC_Palindrome"),
    ]
}

# =============================================================================
# R-LOOP PATTERNS (Class 4) - RNA-DNA hybrid structures
# =============================================================================
R_LOOP_PATTERNS = {
    'rlfs_model1': [
        (r"(?:G{3,}[ACGT]{1,12}){4,}", 70, 0, "RLFS_G_clusters", None, 1.0, 4, 30.0, "RLoop_Score"),
        (r"[GC]{15,}(?=.*[GC]{15,})", 71, 0, "RLFS_GC_rich", None, 1.1, 2, 35.0, "RLoop_Score"),
    ],
    'cpg_rlfs': [
        (r"(?:CG[ACGT]{1,5}){6,}", 72, 0, "CpG_RLFS", None, 1.2, 6, 40.0, "CpG_RLoop"),
        (r"(?:GC[ACGT]{1,5}){6,}", 73, 0, "GpC_RLFS", None, 1.2, 6, 40.0, "CpG_RLoop"),
    ],
    'gc_skew_rlfs': [
        (r"G{3,}[^G]{1,10}G{3,}[^G]{1,10}G{3,}", 74, 0, "GC_skew_positive", None, 1.0, 3, 25.0, "GC_Skew"),
        (r"C{3,}[^C]{1,10}C{3,}[^C]{1,10}C{3,}", 75, 0, "GC_skew_negative", None, 1.0, 3, 25.0, "GC_Skew"),
    ]
}

# =============================================================================
# SLIPPED DNA PATTERNS (Class 2) - Direct and tandem repeats
# =============================================================================
SLIPPED_DNA_PATTERNS = {
    'mononucleotide_repeats': [
        (r"A{12,}", 80, 0, "Poly_A", None, 0.8, 1, 12.0, "Repeat_Score"),
        (r"T{12,}", 81, 0, "Poly_T", None, 0.8, 1, 12.0, "Repeat_Score"),
        (r"G{12,}", 82, 0, "Poly_G", None, 1.0, 1, 12.0, "Repeat_Score"),
        (r"C{12,}", 83, 0, "Poly_C", None, 1.0, 1, 12.0, "Repeat_Score"),
    ],
    'dinucleotide_repeats': [
        (r"(?:AT){6,}", 84, 0, "AT_repeat", None, 1.0, 6, 18.0, "Dinuc_Score"),
        (r"(?:TA){6,}", 85, 0, "TA_repeat", None, 1.0, 6, 18.0, "Dinuc_Score"),
        (r"(?:GC){6,}", 86, 0, "GC_repeat", None, 1.1, 6, 20.0, "Dinuc_Score"),
        (r"(?:CG){6,}", 87, 0, "CG_repeat", None, 1.2, 6, 22.0, "Dinuc_Score"),
        (r"(?:AC){6,}", 88, 0, "AC_repeat", None, 0.9, 6, 16.0, "Dinuc_Score"),
        (r"(?:CA){6,}", 89, 0, "CA_repeat", None, 0.9, 6, 16.0, "Dinuc_Score"),
        (r"(?:AG){6,}", 90, 0, "AG_repeat", None, 0.9, 6, 16.0, "Dinuc_Score"),
        (r"(?:GA){6,}", 91, 0, "GA_repeat", None, 0.9, 6, 16.0, "Dinuc_Score"),
        (r"(?:GT){6,}", 92, 0, "GT_repeat", None, 0.9, 6, 16.0, "Dinuc_Score"),
        (r"(?:TG){6,}", 93, 0, "TG_repeat", None, 0.9, 6, 16.0, "Dinuc_Score"),
        (r"(?:TC){6,}", 94, 0, "TC_repeat", None, 0.9, 6, 16.0, "Dinuc_Score"),
        (r"(?:CT){6,}", 95, 0, "CT_repeat", None, 0.9, 6, 16.0, "Dinuc_Score"),
    ],
    'trinucleotide_repeats': [
        (r"(?:CAG){4,}", 96, 0, "CAG_repeat", None, 1.2, 4, 20.0, "Trinuc_Score"),
        (r"(?:CTG){4,}", 97, 0, "CTG_repeat", None, 1.2, 4, 20.0, "Trinuc_Score"),
        (r"(?:CGG){4,}", 98, 0, "CGG_repeat", None, 1.3, 4, 22.0, "Trinuc_Score"),
        (r"(?:CCG){4,}", 99, 0, "CCG_repeat", None, 1.3, 4, 22.0, "Trinuc_Score"),
        (r"(?:GAA){4,}", 100, 0, "GAA_repeat", None, 1.1, 4, 18.0, "Trinuc_Score"),
        (r"(?:TTC){4,}", 101, 0, "TTC_repeat", None, 1.1, 4, 18.0, "Trinuc_Score"),
        (r"(?:AAG){4,}", 102, 0, "AAG_repeat", None, 1.0, 4, 16.0, "Trinuc_Score"),
        (r"(?:CTT){4,}", 103, 0, "CTT_repeat", None, 1.0, 4, 16.0, "Trinuc_Score"),
        (r"(?:GCA){4,}", 104, 0, "GCA_repeat", None, 1.0, 4, 16.0, "Trinuc_Score"),
        (r"(?:TGC){4,}", 105, 0, "TGC_repeat", None, 1.0, 4, 16.0, "Trinuc_Score"),
    ]
}

# =============================================================================
# PATTERN AGGREGATION AND ACCESS FUNCTIONS
# =============================================================================
ALL_PATTERNS = {
    'G_QUADRUPLEX': G_QUADRUPLEX_PATTERNS,
    'I_MOTIF': I_MOTIF_PATTERNS,
    'Z_DNA': Z_DNA_PATTERNS,
    'CURVED_DNA': CURVED_DNA_PATTERNS,
    'TRIPLEX': TRIPLEX_PATTERNS,
    'CRUCIFORM': CRUCIFORM_PATTERNS,
    'R_LOOP': R_LOOP_PATTERNS,
    'SLIPPED_DNA': SLIPPED_DNA_PATTERNS
}

def get_patterns_for_motif(motif_type: str) -> dict:
    """
    Retrieve all patterns for a specific motif type
    
    | Parameter  | Type | Description                    |
    |------------|------|--------------------------------|
    | motif_type | str  | Motif class name (uppercase)   |
    
    Returns: Dictionary of pattern groups for the motif type
    """
    return ALL_PATTERNS.get(motif_type.upper(), {})

def get_all_hyperscan_patterns() -> list:
    """
    Extract all regex patterns for Hyperscan compilation
    
    Returns: List of (pattern, flags, id) tuples for Hyperscan
    """
    patterns = []
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        for pattern_group, pattern_list in pattern_dict.items():
            for pattern_tuple in pattern_list:
                if pattern_tuple:
                    regex, pattern_id = pattern_tuple[0], pattern_tuple[1]
                    patterns.append((regex, 0, pattern_id))  # flags=0 for basic matching
    return patterns

def get_pattern_info(pattern_id: int) -> dict:
    """
    Retrieve detailed information about a specific pattern
    
    | Field       | Description                           |
    |-------------|---------------------------------------|
    | motif_class | Parent motif class                    |
    | subclass    | Specific motif subclass               |
    | pattern     | Regular expression pattern            |
    | score_method| Scientific scoring algorithm          |
    | min_score   | Minimum threshold for reporting       |
    
    Returns: Dictionary with pattern metadata or None if not found
    """
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        for pattern_group, pattern_list in pattern_dict.items():
            for pattern_tuple in pattern_list:
                if pattern_tuple and pattern_tuple[1] == pattern_id:
                    return {
                        'motif_class': motif_class,
                        'pattern_group': pattern_group,
                        'subclass': pattern_tuple[3],
                        'pattern': pattern_tuple[0],
                        'score_method': pattern_tuple[8],
                        'min_score': pattern_tuple[7],
                        'scale_factor': pattern_tuple[5]
                    }
    return None

def validate_patterns() -> tuple:
    """
    Validate all regex patterns for syntax and Hyperscan compatibility
    
    | Validation Type   | Description                        |
    |-------------------|------------------------------------|
    | Syntax Check      | PCRE regex compilation test        |
    | ID Uniqueness     | Ensure all pattern IDs are unique |
    | Hyperscan Safety  | Check for incompatible constructs  |
    
    Returns: (is_valid: bool, error_messages: list)
    """
    errors = []
    pattern_ids = set()
    
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        for pattern_group, pattern_list in pattern_dict.items():
            for i, pattern_tuple in enumerate(pattern_list):
                if not pattern_tuple or len(pattern_tuple) != 9:
                    errors.append(f"{motif_class}.{pattern_group}[{i}]: Invalid tuple format")
                    continue
                
                regex, pattern_id = pattern_tuple[0], pattern_tuple[1]
                
                # Check pattern ID uniqueness
                if pattern_id in pattern_ids:
                    errors.append(f"{motif_class}.{pattern_group}[{i}]: Duplicate pattern ID {pattern_id}")
                pattern_ids.add(pattern_id)
                
                # Test regex compilation
                try:
                    re.compile(regex)
                except re.error as e:
                    errors.append(f"{motif_class}.{pattern_group}[{i}]: Regex error - {e}")
                
                # Basic Hyperscan compatibility checks
                if '\\b' in regex or '\\B' in regex:
                    errors.append(f"{motif_class}.{pattern_group}[{i}]: Word boundaries not supported in Hyperscan")
                if '(?=' in regex or '(?!' in regex:
                    errors.append(f"{motif_class}.{pattern_group}[{i}]: Lookahead assertions may not be supported")
    
    return len(errors) == 0, errors

# =============================================================================
# CLASSIFICATION CONFIGURATION
# =============================================================================
CLASS_DEFINITIONS = {
    1: {"name": "Curved DNA", "description": "A-tract mediated curvature", "color": "#FF6B6B"},
    2: {"name": "Slipped DNA", "description": "Direct/tandem repeats", "color": "#4ECDC4"}, 
    3: {"name": "Cruciform", "description": "Four-way junctions", "color": "#45B7D1"},
    4: {"name": "R-loop", "description": "RNA-DNA hybrids", "color": "#96CEB4"},
    5: {"name": "Triplex", "description": "Three-stranded DNA", "color": "#FECA57"},
    6: {"name": "G-Quadruplex", "description": "Four-stranded G-rich", "color": "#FF9FF3"},
    7: {"name": "i-Motif", "description": "C-rich quadruplex", "color": "#54A0FF"},
    8: {"name": "Z-DNA", "description": "Left-handed helix", "color": "#5F27CD"}
}

DEFAULT_SELECTED_CLASSES = [1, 2, 3, 4, 5, 6, 7, 8]
DEFAULT_SELECTED_SUBCLASSES = list(range(1, 25))

MOTIF_LENGTH_LIMITS = {
    'G-Quadruplex': {'min': 15, 'max': 200},
    'i-Motif': {'min': 15, 'max': 200}, 
    'Z-DNA': {'min': 12, 'max': 100},
    'Curved DNA': {'min': 8, 'max': 150},
    'Triplex': {'min': 12, 'max': 150},
    'Cruciform': {'min': 8, 'max': 200},
    'R-loop': {'min': 15, 'max': 300},
    'Slipped DNA': {'min': 6, 'max': 500}
}

def get_motif_limits(motif_class: str) -> dict:
    """
    Get length limits for a specific motif class
    
    | Parameter    | Type | Description              |
    |--------------|------|--------------------------|
    | motif_class  | str  | Motif class name         |
    
    Returns: Dictionary with 'min' and 'max' length limits
    """
    return MOTIF_LENGTH_LIMITS.get(motif_class, {'min': 6, 'max': 500})