"""
Central Regex Registry for NBDFinder - All Hyperscan-Safe Patterns
================================================================

This module contains all regex patterns used by NBDFinder motif detection modules,
organized by motif class for shared Hyperscan database compilation and maintenance.

Scientific References and Standards:
====================================
- G4Hunter: Bedrat et al. NAR 44(4):1746-1759 (2016) - G-quadruplex prediction algorithm
- Z-DNA Seeker: Ho et al. EMBO J 5:2737-2744 (1986); Rich & Zhang Nature 302:209-217 (1983)
- i-motif: Zeraati et al. Nat Chem 10:631-637 (2018) - C-rich quadruplex structures
- Triplex: Frank-Kamenetskii & Mirkin Annu Rev Biochem 64:65-95 (1995) - Triple helix DNA
- Cruciform: Lilley & Clegg Annu Rev Biophys Biomol Struct 22:299-328 (1993) - Inverted repeats
- Curved DNA: Crothers et al. Methods Enzymol 212:3-29 (1992); Marini et al. Cell 28:871 (1982)
- R-loops: Ginno et al. Mol Cell 45:511-522 (2012); Skourti-Stathaki et al. NSMB 18:1325 (2011)
- Slipped DNA: Wells et al. Trends Biochem Sci 30:437-444 (2005); Pearson et al. Nat Rev Genet 6:729 (2005)

Pattern Structure (9-tuple format):
==================================
Each pattern is a tuple: (regex_pattern, pattern_id, group_number, subclass_name, 
                         scoring_function, score_scale, min_runs, min_score, score_method)

Where:
- regex_pattern: PCRE-compatible regular expression (Hyperscan-safe)
- pattern_id: Unique integer identifier (1-999)
- group_number: Capture group index (usually 0 for full match)
- subclass_name: Scientific motif subclass name
- scoring_function: Optional Python function for advanced scoring (None for basic)
- score_scale: Scaling factor for motif importance (0.1-2.0)
- min_runs: Minimum number of pattern repeats required
- min_score: Minimum threshold score for reporting motif
- score_method: Scientific scoring algorithm name

Hyperscan Compatibility:
========================
All patterns are designed for Intel Hyperscan PCRE compatibility:
- No back-references or recursive patterns
- No lookahead/lookbehind assertions (except where marked)
- Fixed repetition counts preferred over unbounded quantifiers
- Character classes explicitly defined [ACGT] instead of \\w

Performance Optimization:
========================
Patterns are ordered by specificity and biological importance:
1. High-confidence canonical patterns first
2. Relaxed/variant patterns second  
3. Experimental/research patterns last
"""

import re
import numpy as np

# === G-QUADRUPLEX FAMILY PATTERNS (Class 6) ===
# Pattern Table: G4Hunter-compliant definitions (Bedrat et al. 2016, NAR 44:1746-1759)
# | Pattern Type      | Definition                    | Min Score | Scientific Ref                      | Loop Range |
# |-------------------|-------------------------------|-----------|-------------------------------------|------------|
# | Canonical G4      | G3+N1-7G3+N1-7G3+N1-7G3+     | 1.2       | Williamson 2005, Sen & Gilbert 1990| 1-7 bp     |
# | Relaxed G4        | G2+N1-12G2+N1-12G2+N1-12G2+  | 0.8       | Kikin et al. 2006, Huppert 2007    | 1-12 bp    |
# | Bulged G4         | G3+N0-3G3+N0-3G3+N0-3G3+     | 0.7       | Todd et al. 2005, Vorlíčková 2012  | 0-3 bp     |
# | Long-loop G4      | G3+N8-15G3+N8-15G3+N8-15G3+  | 0.6       | Bochman et al. 2012, EMBO J        | 8-15 bp    |
# | Two-quartet G4    | G2N1-7G2N1-7G2N1-7G2         | 0.5       | Adrian et al. 2014, NAR            | 1-7 bp     |
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
    'two_quartet_g4': [
        (r"G{2}[ACGT]{1,7}G{2}[ACGT]{1,7}G{2}[ACGT]{1,7}G{2}", 7, 0, "Two_Quartet_G4", None, 0.5, 4, 0.5, "G4Hunter"),
    ],
    'imperfect_g4': [
        (r"G{2,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}", 8, 0, "Imperfect_G4_pos1", None, 0.9, 4, 0.8, "G4Hunter"),
        (r"G{3,}[ACGT]{1,7}G{2,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}", 9, 0, "Imperfect_G4_pos2", None, 0.9, 4, 0.8, "G4Hunter"),
        (r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{2,}[ACGT]{1,7}G{3,}", 10, 0, "Imperfect_G4_pos3", None, 0.9, 4, 0.8, "G4Hunter"),
        (r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{2,}", 11, 0, "Imperfect_G4_pos4", None, 0.9, 4, 0.8, "G4Hunter"),
    ],
    'multimeric_g4': [
        (r"(?:G{3,}[ACGT]{1,12}){4,}", 12, 0, "Multimeric_G4", None, 1.2, 4, 1.0, "G4Hunter_Extended"),
        (r"(?:G{2,}[ACGT]{1,12}){5,}", 13, 0, "Extended_Multimeric_G4", None, 1.0, 5, 0.8, "G4Hunter_Extended"),
    ],
    'bipartite_g4': [
        (r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{30,100}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}", 14, 0, "Bipartite_G4", None, 0.9, 8, 0.6, "Bipartite_G4Hunter"),
    ],
    'g_triplex': [
        (r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}", 15, 0, "G_Triplex", None, 1.0, 3, 0.8, "G3_Score"),
        (r"G{4,}[ACGT]{1,7}G{4,}[ACGT]{1,7}G{4,}", 16, 0, "Perfect_G_Triplex", None, 1.2, 3, 1.0, "G3_Score"),
    ]
}

# === I-MOTIF FAMILY PATTERNS (Class 7) ===
# Pattern Table: C-rich quadruplex structures (Zeraati et al. 2018, Nat Chem 10:631-637)
# | Pattern Type      | Definition                    | pH Range  | Scientific Ref                     | Stability  |
# |-------------------|-------------------------------|-----------|-----------------------------------|------------|
# | Canonical iMotif  | C3+N1-7C3+N1-7C3+N1-7C3+     | pH 5.0    | Gehring et al. 1993, Leroy 1993  | High       |
# | Relaxed iMotif    | C2+N1-12C2+N1-12C2+N1-12C2+  | pH 5.5    | Zhou et al. 2010, Mergny 2009    | Medium     |
# | AC-motif          | A3N4-6C3N4-6C3N4-6C3          | pH 6.0    | Felsenfeld 1967, Jain 2019       | Low        |
# | Intercalated      | C3+N0-3C3+N0-3C3+N0-3C3+     | pH 4.5    | Dzatko et al. 2018, Kaiser 2017  | Very High  |
I_MOTIF_PATTERNS = {
    'canonical_imotif': [
        (r"C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}", 17, 0, "Canonical_iMotif", None, 1.0, 4, 1.2, "iM_Hunter"),
        (r"C{4,}[ACGT]{1,7}C{4,}[ACGT]{1,7}C{4,}[ACGT]{1,7}C{4,}", 18, 0, "Perfect_iMotif", None, 1.2, 4, 1.5, "iM_Hunter"),
    ],
    'relaxed_imotif': [
        (r"C{2,}[ACGT]{1,12}C{2,}[ACGT]{1,12}C{2,}[ACGT]{1,12}C{2,}", 19, 0, "Relaxed_iMotif", None, 0.8, 4, 0.8, "iM_Hunter"),
        (r"C{3,}[ACGT]{8,12}C{3,}[ACGT]{8,12}C{3,}[ACGT]{8,12}C{3,}", 20, 0, "Long_Loop_iMotif", None, 0.8, 4, 0.8, "iM_Hunter"),
    ],
    'intercalated_imotif': [
        (r"C{3,}[ACGT]{0,3}C{3,}[ACGT]{0,3}C{3,}[ACGT]{0,3}C{3,}", 21, 0, "Intercalated_iMotif", None, 1.3, 4, 1.3, "iM_Hunter"),
    ],
    'ac_motif': [
        (r"A{3,}[ACGT]{4,6}C{3,}[ACGT]{4,6}C{3,}[ACGT]{4,6}C{3,}", 22, 0, "AC_motif_type1", None, 0.8, 4, 0.6, "AC_Score"),
        (r"C{3,}[ACGT]{4,6}C{3,}[ACGT]{4,6}C{3,}[ACGT]{4,6}A{3,}", 23, 0, "AC_motif_type2", None, 0.8, 4, 0.6, "AC_Score"),
        (r"A{3,}[ACGT]{4,8}C{3,}[ACGT]{4,8}C{3,}[ACGT]{4,8}A{3,}", 24, 0, "Symmetric_AC_motif", None, 0.7, 4, 0.5, "AC_Score"),
    ],
    'imperfect_imotif': [
        (r"C{2,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}", 25, 0, "Imperfect_iMotif_pos1", None, 0.9, 4, 0.8, "iM_Hunter"),
        (r"C{3,}[ACGT]{1,7}C{2,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}", 26, 0, "Imperfect_iMotif_pos2", None, 0.9, 4, 0.8, "iM_Hunter"),
        (r"C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{2,}[ACGT]{1,7}C{3,}", 27, 0, "Imperfect_iMotif_pos3", None, 0.9, 4, 0.8, "iM_Hunter"),
        (r"C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{2,}", 28, 0, "Imperfect_iMotif_pos4", None, 0.9, 4, 0.8, "iM_Hunter"),
    ]
}

# === Z-DNA PATTERNS (Class 8) ===
# Pattern Table: Left-handed DNA structures (Ho et al. 1986, Rich & Zhang 2003)
# | Pattern Type      | Definition              | Z-Score   | Scientific Ref                     | Transition |
# |-------------------|-------------------------|-----------|-----------------------------------|------------|
# | CG Z-DNA          | (CG){6,} alternating    | >8.0      | Rich & Zhang 1981, Wang 1979     | High       |
# | AT Z-DNA          | (AT){8,} alternating    | >3.0      | Nordheim & Rich 1983, PNAS       | Medium     |
# | Mixed Z-DNA       | CG/AT mixed patterns    | >5.0      | Ellison et al. 1986, Cell        | Variable   |
# | eGZ DNA           | CGG expansion slippage  | >6.0      | Napierala et al. 1997, EMBOJ     | Disease    |
Z_DNA_PATTERNS = {
    'cg_zdna': [
        (r"(?:CG){6,}", 29, 0, "CG_Z_DNA", None, 1.0, 6, 8.0, "Z_Seeker"),
        (r"(?:GC){6,}", 30, 0, "GC_Z_DNA", None, 1.0, 6, 8.0, "Z_Seeker"),
        (r"[CG]{12,}", 31, 0, "CG_Rich_Z_DNA", None, 1.0, 12, 6.0, "Z_Seeker"),
    ],
    'at_zdna': [
        (r"(?:AT){8,}", 32, 0, "AT_Z_DNA", None, 0.8, 8, 3.0, "Z_Seeker"),
        (r"(?:TA){8,}", 33, 0, "TA_Z_DNA", None, 0.8, 8, 3.0, "Z_Seeker"),
    ],
    'mixed_zdna': [
        (r"(?:[CG]{2}[AT]{2}){4,}", 34, 0, "Mixed_Z_DNA_CGAT", None, 0.9, 4, 5.0, "Z_Seeker"),
        (r"(?:[AT]{2}[CG]{2}){4,}", 35, 0, "Mixed_Z_DNA_ATCG", None, 0.9, 4, 5.0, "Z_Seeker"),
        # Simplified pattern for general mixed Z-DNA
        (r"[CG]{6}[AT]{6}[CG]{6}", 36, 0, "General_Mixed_Z_DNA", None, 0.7, 12, 4.0, "Z_Seeker"),
    ],
    'egz_dna': [
        (r"(?:CGG){10,}", 37, 0, "eGZ_CGG_Expansion", None, 1.2, 10, 6.0, "Expansion_Score"),
        (r"(?:CCG){10,}", 38, 0, "eGZ_CCG_Expansion", None, 1.2, 10, 6.0, "Expansion_Score"),
        (r"G[CG]{8,}G", 39, 0, "eGZ_Extruded_G", None, 1.0, 8, 5.0, "Z_Seeker"),
    ],
    'z_dna_junction': [
        (r"[ACGT]{3,8}(?:CG){4,}[ACGT]{3,8}", 40, 0, "Z_B_Junction", None, 0.8, 4, 4.0, "Junction_Score"),
        (r"[ACGT]{3,8}(?:AT){6,}[ACGT]{3,8}", 41, 0, "Z_B_Junction_AT", None, 0.6, 6, 2.5, "Junction_Score"),
    ]
}

# === CURVED DNA PATTERNS (Class 1) ===
# Pattern Table: A-tract curvature and phasing (Crothers et al. 1992, Marini et al. 1982)
# | Pattern Type      | Definition           | Period    | Scientific Ref                        | Curvature  |
# |-------------------|---------------------|-----------|---------------------------------------|------------|
# | A-tract           | A4+ or T4+ runs     | N/A       | Marini et al. 1982, Cell 28:871      | High       |
# | Phased A-tract    | ~10 bp periodic     | 10.0-10.6 | Crothers et al. 1992, Meth Enzymol   | Very High  |
# | Mixed AT-tract    | Mixed A/T regions   | Variable  | Olson et al. 1998, PNAS 95:11163     | Medium     |
# | Global Array      | Multiple phased     | 10.0-11.0 | Yella & Bansal 2017, Sci Rep         | Maximum    |
CURVED_DNA_PATTERNS = {
    'a_tracts': [
        (r"A{4,}", 42, 0, "A_tract", None, 1.0, 4, 0.5, "Curvature_Calladine"),
        (r"A{7,}", 43, 0, "Long_A_tract", None, 1.2, 7, 0.8, "Curvature_Calladine"),
        (r"A{10,}", 44, 0, "Very_Long_A_tract", None, 1.5, 10, 1.2, "Curvature_Calladine"),
    ],
    't_tracts': [
        (r"T{4,}", 45, 0, "T_tract", None, 1.0, 4, 0.5, "Curvature_Calladine"),
        (r"T{7,}", 46, 0, "Long_T_tract", None, 1.2, 7, 0.8, "Curvature_Calladine"),
        (r"T{10,}", 47, 0, "Very_Long_T_tract", None, 1.5, 10, 1.2, "Curvature_Calladine"),
    ],
    'phased_a_tracts': [
        (r"A{3,}[ACGT]{6,14}A{3,}[ACGT]{6,14}A{3,}", 48, 0, "Phased_A_tract_3x", None, 1.3, 3, 1.0, "Phasing_Score"),
        (r"A{4,}[ACGT]{6,14}A{4,}[ACGT]{6,14}A{4,}[ACGT]{6,14}A{4,}", 49, 0, "Phased_A_tract_4x", None, 1.5, 4, 1.2, "Phasing_Score"),
    ],
    'phased_t_tracts': [
        (r"T{3,}[ACGT]{6,14}T{3,}[ACGT]{6,14}T{3,}", 50, 0, "Phased_T_tract_3x", None, 1.3, 3, 1.0, "Phasing_Score"),
        (r"T{4,}[ACGT]{6,14}T{4,}[ACGT]{6,14}T{4,}[ACGT]{6,14}T{4,}", 51, 0, "Phased_T_tract_4x", None, 1.5, 4, 1.2, "Phasing_Score"),
    ],
    'mixed_at_tracts': [
        (r"[AT]{6,}", 52, 0, "Mixed_AT_tract", None, 0.8, 6, 0.3, "AT_Content_Score"),
        (r"[AT]{10,}", 53, 0, "Long_Mixed_AT_tract", None, 1.0, 10, 0.5, "AT_Content_Score"),
    ],
    'global_arrays': [
        (r"(?:[AT]{3,}[ACGT]{7,13}){3,}", 54, 0, "Global_Array_candidate", None, 1.4, 3, 1.0, "Global_Array_Score"),
        (r"(?:A{3,}[ACGT]{8,12}){4,}", 55, 0, "Perfect_A_Array", None, 1.6, 4, 1.3, "Global_Array_Score"),
        (r"(?:T{3,}[ACGT]{8,12}){4,}", 56, 0, "Perfect_T_Array", None, 1.6, 4, 1.3, "Global_Array_Score"),
    ]
}

# === TRIPLEX PATTERNS (Class 5) ===
# Pattern Table: Triple-stranded DNA (Frank-Kamenetskii & Mirkin 1995, Annu Rev Biochem 64:65-95)
# | Pattern Type      | Definition              | Length    | Scientific Ref                        | Stability  |
# |-------------------|------------------------|-----------|---------------------------------------|------------|
# | Homopurine        | A/G tracts ≥15 bp     | 15+ bp    | Lyamichev et al. 1986, J Biomol Str  | High       |
# | Homopyrimidine    | C/T tracts ≥15 bp     | 15+ bp    | Moser & Dervan 1987, Science         | High       |
# | Purine-rich       | >70% A/G content      | 20+ bp    | Thuong & Hélène 1993, Angew Chem     | Medium     |
# | Mirror Repeats    | Palindromic Pu/Py     | 15+ bp    | Htun & Dahlberg 1989, Science        | Variable   |
# | Sticky DNA        | Homopyrimidine runs   | 15+ bp    | Kohwi & Kohwi-Shigematsu 1988 Cell   | Medium     |
TRIPLEX_PATTERNS = {
    'homopurine_tracts': [
        (r"[AG]{15,}", 57, 0, "Homopurine_Triplex", None, 1.0, 15, 0.8, "Triplex_Stability"),
        (r"[AG]{20,}", 58, 0, "Long_Homopurine_Triplex", None, 1.2, 20, 1.0, "Triplex_Stability"),
        (r"[AG]{30,}", 59, 0, "Very_Long_Homopurine", None, 1.5, 30, 1.2, "Triplex_Stability"),
    ],
    'homopyrimidine_tracts': [
        (r"[CT]{15,}", 60, 0, "Homopyrimidine_Triplex", None, 1.0, 15, 0.8, "Triplex_Stability"),
        (r"[CT]{20,}", 61, 0, "Long_Homopyrimidine_Triplex", None, 1.2, 20, 1.0, "Triplex_Stability"),
        (r"[CT]{30,}", 62, 0, "Very_Long_Homopyrimidine", None, 1.5, 30, 1.2, "Triplex_Stability"),
    ],
    'purine_rich_regions': [
        # Simplified pattern for purine-rich regions (70%+ AG content estimated)
        (r"[AG]{14}[ACGT]{0,6}[AG]{14}[ACGT]{0,6}[AG]{14}", 63, 0, "Purine_Rich_Triplex", None, 0.9, 20, 0.7, "Purine_Content"),
    ],
    'mirror_repeats': [
        (r"[AG]{8,}[ACGT]{0,10}[CT]{8,}", 64, 0, "Mirror_Repeat_Pu_Py", None, 1.1, 16, 0.8, "Mirror_Score"),
        (r"[CT]{8,}[ACGT]{0,10}[AG]{8,}", 65, 0, "Mirror_Repeat_Py_Pu", None, 1.1, 16, 0.8, "Mirror_Score"),
    ],
    'sticky_dna': [
        (r"[CT]{15,}", 66, 0, "Sticky_DNA", None, 1.0, 15, 0.8, "Sticky_Score"),
        (r"C{8,}[CT]{7,}", 67, 0, "C_Rich_Sticky_DNA", None, 1.2, 15, 1.0, "Sticky_Score"),
        (r"T{8,}[CT]{7,}", 68, 0, "T_Rich_Sticky_DNA", None, 1.2, 15, 1.0, "Sticky_Score"),
    ],
    'intermolecular_triplex': [
        (r"[AG]{12,}[ACGT]{10,50}[CT]{12,}", 69, 0, "Intermolecular_Triplex", None, 1.3, 2, 0.9, "Intermolecular_Score"),
    ],
    'intramolecular_triplex': [
        (r"[AG]{10,}[ACGT]{5,30}[CT]{10,}[ACGT]{5,30}[AG]{10,}", 70, 0, "Intramolecular_H_DNA", None, 1.4, 3, 1.0, "H_DNA_Score"),
    ]
}

# === CRUCIFORM PATTERNS (Class 3) ===
# Pattern Table: Inverted repeats and palindromes (Lilley & Clegg 1993, Annu Rev Biophys 22:299-328)
# | Pattern Type      | Definition                | Arm Len   | Scientific Ref                        | Stability  |
# |-------------------|---------------------------|-----------|---------------------------------------|------------|
# | Perfect Palindrome| Exact inverted repeats    | 6-20 bp   | Lilley & Clegg 1993, Annu Rev Biophy | Very High  |
# | Imperfect Palindrome| Minor mismatches allowed| 8-25 bp   | Pearson & Sinden 1996, Biochemistry  | High       |
# | Stem-loop         | IR with variable loop     | 4-15 bp   | Rouillard et al. 2003, NAR           | Medium     |
# | Long Inverted     | Extended palindromes      | 20+ bp    | Gordenin et al. 1993, PNAS           | Maximum    |
# | AT-rich Palindrome| A/T enriched cruciforms   | 6+ bp     | McClellan et al. 1990, J Biol Chem   | Variable   |
CRUCIFORM_PATTERNS = {
    'perfect_palindromes': [
        # Short perfect palindromes (6-10 bp arms)
        (r"[ACGT]{6}[ACGT]{0,10}[ACGT]{6}", 119, 0, "Perfect_Palindrome_Short", None, 1.2, 1, 0.8, "Palindrome_Score"),
        (r"[ACGT]{8}[ACGT]{0,15}[ACGT]{8}", 120, 0, "Perfect_Palindrome_Medium", None, 1.4, 1, 1.0, "Palindrome_Score"),
        (r"[ACGT]{10}[ACGT]{0,20}[ACGT]{10}", 121, 0, "Perfect_Palindrome_Long", None, 1.6, 1, 1.2, "Palindrome_Score"),
    ],
    'imperfect_palindromes': [
        # Allowing 1-2 mismatches in arms
        (r"[ACGT]{8,}[ACGT]{1,30}[ACGT]{8,}", 122, 0, "Imperfect_Palindrome_Candidate", None, 1.0, 1, 0.6, "Mismatch_Palindrome"),
        (r"[ACGT]{12,}[ACGT]{1,50}[ACGT]{12,}", 123, 0, "Long_Imperfect_Palindrome", None, 1.2, 1, 0.8, "Mismatch_Palindrome"),
    ],
    'at_rich_palindromes': [
        (r"[AT]{6,}[ACGT]{0,15}[AT]{6,}", 124, 0, "AT_Rich_Palindrome", None, 1.1, 1, 0.7, "AT_Palindrome_Score"),
        (r"[AT]{8,}[ACGT]{0,25}[AT]{8,}", 125, 0, "Long_AT_Palindrome", None, 1.3, 1, 0.9, "AT_Palindrome_Score"),
    ],
    'gc_rich_palindromes': [
        (r"[GC]{6,}[ACGT]{0,15}[GC]{6,}", 126, 0, "GC_Rich_Palindrome", None, 1.2, 1, 0.8, "GC_Palindrome_Score"),
        (r"[GC]{8,}[ACGT]{0,25}[GC]{8,}", 127, 0, "Long_GC_Palindrome", None, 1.4, 1, 1.0, "GC_Palindrome_Score"),
    ],
    'stem_loop_structures': [
        (r"[ACGT]{6}[ACGT]{3,15}[ACGT]{6}", 128, 0, "Stem_Loop_Short", None, 1.0, 1, 0.6, "Stem_Loop_Score"),
        (r"[ACGT]{8}[ACGT]{4,20}[ACGT]{8}", 129, 0, "Stem_Loop_Medium", None, 1.2, 1, 0.8, "Stem_Loop_Score"),
        (r"[ACGT]{10}[ACGT]{5,30}[ACGT]{10}", 130, 0, "Stem_Loop_Long", None, 1.4, 1, 1.0, "Stem_Loop_Score"),
    ],
    'direct_repeats': [
        # Pre-filter candidates for direct repeat analysis
        (r"[ACGT]{8,20}", 131, 0, "Direct_Repeat_Candidate", None, 0.8, 1, 0.4, "Repeat_Candidate"),
        (r"[ACGT]{15,30}", 132, 0, "Long_Direct_Repeat_Candidate", None, 1.0, 1, 0.6, "Repeat_Candidate"),
    ],
    'inverted_repeat_candidates': [],  # Generated dynamically
    'palindrome_candidates': []        # Generated dynamically
}

# === R-LOOP PATTERNS (Class 4) ===
# Pattern Table: R-loop forming sequences (Ginno et al. 2012, Skourti-Stathaki et al. 2011)
# | Pattern Type      | Definition                | GC%       | Scientific Ref                        | Stability  |
# |-------------------|---------------------------|-----------|---------------------------------------|------------|
# | RLFS Model 1      | G3+N1-10G3+(G3+N1-10)+   | >50%      | Ginno et al. 2012, Mol Cell          | High       |
# | RLFS Model 2      | G4+(G4+N1-10)+           | >60%      | Skourti-Stathaki et al. 2011, NSMB   | Very High  |
# | QmRLFS-RE         | G-quadruplex + G-rich    | >65%      | Wongsurawat et al. 2012, NAR         | Maximum    |
# | CpG RLFS          | CpG island associated    | >55%      | Santos-Pereira & Aguilera 2015, Cell | High       |
# | GC Skew RLFS      | G/C asymmetric regions   | >50%      | Ginno et al. 2013, Nature            | Medium     |
R_LOOP_PATTERNS = {
    'rlfs_model1': [
        (r"G{3,}[ACGT]{1,10}G{3,}(?:[ACGT]{1,10}G{3,}){1,}", 71, 0, "RLFS_Model1", None, 1.0, 2, 0.8, "QmRLFS"),
        (r"G{3,}[ACGT]{1,10}G{3,}(?:[ACGT]{1,10}G{3,}){2,}", 72, 0, "Extended_RLFS_Model1", None, 1.2, 3, 1.0, "QmRLFS"),
    ],
    'rlfs_model2': [
        (r"G{4,}(?:[ACGT]{1,10}G{4,}){1,}", 73, 0, "RLFS_Model2", None, 1.2, 2, 1.0, "QmRLFS"),
        (r"G{4,}(?:[ACGT]{1,10}G{4,}){2,}", 74, 0, "Extended_RLFS_Model2", None, 1.4, 3, 1.2, "QmRLFS"),
    ],
    'qmrlfs_re': [
        (r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{10,100}G{4,}(?:[ACGT]{1,10}G{4,}){1,}", 75, 0, "QmRLFS_RE", None, 1.5, 5, 1.3, "Combined_Score"),
    ],
    'cpg_rlfs': [
        (r"(?:CG[ACGT]{0,3}){5,}G{3,}[ACGT]{1,10}G{3,}", 76, 0, "CpG_RLFS", None, 1.1, 5, 0.9, "CpG_QmRLFS"),
        (r"(?:GC[ACGT]{0,3}){5,}G{3,}[ACGT]{1,10}G{3,}", 77, 0, "GC_Rich_RLFS", None, 1.1, 5, 0.9, "CpG_QmRLFS"),
    ],
    'gc_skew_rlfs': [
        (r"G{20,}[CT]{5,50}G{3,}[ACGT]{1,10}G{3,}", 78, 0, "GC_Skew_RLFS", None, 1.3, 3, 1.0, "GC_Skew_Score"),
        (r"C{20,}[AG]{5,50}C{3,}[ACGT]{1,10}C{3,}", 79, 0, "Complementary_GC_Skew", None, 1.3, 3, 1.0, "GC_Skew_Score"),
    ],
    'transcription_rlfs': [
        (r"G{3,}[ACGT]{1,15}G{3,}[ACGT]{1,15}G{3,}", 80, 0, "Transcription_RLFS", None, 1.0, 3, 0.8, "Transcription_Score"),
    ],
    'rnh1_target': [
        (r"[AG]{8,}G{3,}[ACGT]{1,10}G{3,}[CT]{8,}", 81, 0, "RNase_H1_Target", None, 1.2, 3, 1.0, "RNH1_Score"),
    ]
}

# === SLIPPED DNA PATTERNS (Class 2) ===
# Pattern Table: Short Tandem Repeats and Slipped Structures (Wells et al. 2005, Pearson et al. 2005)
# | Pattern Type      | Definition              | Disease   | Scientific Ref                        | Expansion  |
# |-------------------|------------------------|-----------|---------------------------------------|------------|
# | Dinucleotide STR  | (XY)n repeats         | Various   | Pearson et al. 2005, Nat Rev Genet   | Common     |
# | Trinucleotide STR | (XYZ)n repeats        | Severe    | Orr & Zoghbi 2007, Annu Rev Neurosci | Pathogenic |
# | Tetranucleotide   | (WXYZ)n repeats       | Rare      | Gatchel & Zoghbi 2005, Nat Rev Genet | Mild       |
# | Mononucleotide    | X10+ runs             | Cancer    | Boland & Goel 2010, Gastroenterology | MSI        |
# | Complex STR       | Mixed repeat patterns | Complex   | Gymrek et al. 2016, Nat Genet        | Variable   |
SLIPPED_DNA_PATTERNS = {
    'mononucleotide_repeats': [
        (r"A{8,}", 82, 0, "Poly_A_STR", None, 1.0, 8, 0.5, "STR_Score"),
        (r"T{8,}", 83, 0, "Poly_T_STR", None, 1.0, 8, 0.5, "STR_Score"),
        (r"G{8,}", 84, 0, "Poly_G_STR", None, 1.0, 8, 0.5, "STR_Score"),
        (r"C{8,}", 85, 0, "Poly_C_STR", None, 1.0, 8, 0.5, "STR_Score"),
        (r"A{12,}", 86, 0, "Long_Poly_A_STR", None, 1.2, 12, 0.8, "STR_Score"),
        (r"T{12,}", 87, 0, "Long_Poly_T_STR", None, 1.2, 12, 0.8, "STR_Score"),
    ],
    'dinucleotide_repeats': [
        (r"(?:AT){6,}", 88, 0, "AT_Dinucleotide_STR", None, 1.0, 6, 0.6, "STR_Score"),
        (r"(?:TA){6,}", 89, 0, "TA_Dinucleotide_STR", None, 1.0, 6, 0.6, "STR_Score"),
        (r"(?:GC){6,}", 90, 0, "GC_Dinucleotide_STR", None, 1.0, 6, 0.6, "STR_Score"),
        (r"(?:CG){6,}", 91, 0, "CG_Dinucleotide_STR", None, 1.0, 6, 0.6, "STR_Score"),
        (r"(?:AG){6,}", 92, 0, "AG_Dinucleotide_STR", None, 1.0, 6, 0.6, "STR_Score"),
        (r"(?:GA){6,}", 93, 0, "GA_Dinucleotide_STR", None, 1.0, 6, 0.6, "STR_Score"),
        (r"(?:CT){6,}", 94, 0, "CT_Dinucleotide_STR", None, 1.0, 6, 0.6, "STR_Score"),
        (r"(?:TC){6,}", 95, 0, "TC_Dinucleotide_STR", None, 1.0, 6, 0.6, "STR_Score"),
        (r"(?:AC){6,}", 96, 0, "AC_Dinucleotide_STR", None, 1.0, 6, 0.6, "STR_Score"),
        (r"(?:CA){6,}", 97, 0, "CA_Dinucleotide_STR", None, 1.0, 6, 0.6, "STR_Score"),
        (r"(?:GT){6,}", 98, 0, "GT_Dinucleotide_STR", None, 1.0, 6, 0.6, "STR_Score"),
        (r"(?:TG){6,}", 99, 0, "TG_Dinucleotide_STR", None, 1.0, 6, 0.6, "STR_Score"),
    ],
    'trinucleotide_repeats': [
        (r"(?:CAG){5,}", 100, 0, "CAG_Trinucleotide_STR", None, 1.2, 5, 0.8, "Disease_STR"),
        (r"(?:CTG){5,}", 101, 0, "CTG_Trinucleotide_STR", None, 1.2, 5, 0.8, "Disease_STR"),
        (r"(?:CGG){5,}", 102, 0, "CGG_Trinucleotide_STR", None, 1.2, 5, 0.8, "Disease_STR"),
        (r"(?:CCG){5,}", 103, 0, "CCG_Trinucleotide_STR", None, 1.2, 5, 0.8, "Disease_STR"),
        (r"(?:GAA){5,}", 104, 0, "GAA_Trinucleotide_STR", None, 1.2, 5, 0.8, "Disease_STR"),
        (r"(?:TTC){5,}", 105, 0, "TTC_Trinucleotide_STR", None, 1.2, 5, 0.8, "Disease_STR"),
        (r"(?:AAG){5,}", 106, 0, "AAG_Trinucleotide_STR", None, 1.0, 5, 0.7, "STR_Score"),
        (r"(?:CTT){5,}", 107, 0, "CTT_Trinucleotide_STR", None, 1.0, 5, 0.7, "STR_Score"),
        (r"(?:GCA){5,}", 108, 0, "GCA_Trinucleotide_STR", None, 1.0, 5, 0.7, "STR_Score"),
        (r"(?:TGC){5,}", 109, 0, "TGC_Trinucleotide_STR", None, 1.0, 5, 0.7, "STR_Score"),
    ],
    'tetranucleotide_repeats': [
        (r"(?:AAAG){4,}", 110, 0, "AAAG_Tetranucleotide_STR", None, 1.0, 4, 0.6, "STR_Score"),
        (r"(?:CTTT){4,}", 111, 0, "CTTT_Tetranucleotide_STR", None, 1.0, 4, 0.6, "STR_Score"),
        (r"(?:GATA){4,}", 112, 0, "GATA_Tetranucleotide_STR", None, 1.0, 4, 0.6, "STR_Score"),
        (r"(?:TATC){4,}", 113, 0, "TATC_Tetranucleotide_STR", None, 1.0, 4, 0.6, "STR_Score"),
    ],
    'pentanucleotide_repeats': [
        (r"(?:ATTCT){3,}", 114, 0, "ATTCT_Pentanucleotide_STR", None, 0.9, 3, 0.5, "STR_Score"),
        (r"(?:AGAAT){3,}", 115, 0, "AGAAT_Pentanucleotide_STR", None, 0.9, 3, 0.5, "STR_Score"),
    ],
    'hexanucleotide_repeats': [
        (r"(?:GGCCCC){3,}", 116, 0, "GGCCCC_Hexanucleotide_STR", None, 1.1, 3, 0.7, "Disease_STR"),
        (r"(?:GGGGCC){3,}", 117, 0, "GGGGCC_Hexanucleotide_STR", None, 1.1, 3, 0.7, "Disease_STR"),
    ],
    'complex_repeats': [
        # Simplified pattern for complex repeat detection - pre-filter for post-processing
        (r"[ACGT]{12,}", 118, 0, "Complex_Repeat_Candidate", None, 0.8, 3, 0.4, "Complex_STR"),
    ]
}

# === MASTER PATTERN REGISTRY ===
ALL_PATTERNS = {
    'g_quadruplex': G_QUADRUPLEX_PATTERNS,
    'i_motif': I_MOTIF_PATTERNS,
    'z_dna': Z_DNA_PATTERNS,
    'curved_dna': CURVED_DNA_PATTERNS,
    'triplex': TRIPLEX_PATTERNS,
    'cruciform': CRUCIFORM_PATTERNS,
    'r_loop': R_LOOP_PATTERNS,
    'slipped_dna': SLIPPED_DNA_PATTERNS,
}

# === UTILITY FUNCTIONS ===

def get_patterns_for_motif(motif_class: str) -> dict:
    """
    Get all patterns for a specific motif class.
    
    Args:
        motif_class: One of 'g_quadruplex', 'i_motif', 'z_dna', 'curved_dna', 
                    'triplex', 'cruciform', 'r_loop', 'slipped_dna'
    
    Returns:
        Dictionary of pattern categories for the motif class
    """
    return ALL_PATTERNS.get(motif_class, {})

def get_all_hyperscan_patterns() -> list:
    """
    Get all patterns suitable for Hyperscan compilation.
    
    Returns:
        List of (pattern, id) tuples for Hyperscan database compilation
    """
    all_patterns = []
    pattern_id = 0
    
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        for pattern_type, patterns in pattern_dict.items():
            for pattern_tuple in patterns:
                if pattern_tuple:  # Skip empty patterns
                    regex_pattern = pattern_tuple[0]
                    all_patterns.append((regex_pattern, pattern_id))
                    pattern_id += 1
    
    return all_patterns

def get_pattern_info(pattern_id: int) -> dict:
    """
    Get detailed information about a pattern by its ID.
    
    Args:
        pattern_id: Unique pattern identifier
        
    Returns:
        Dictionary with pattern details (motif_class, pattern_type, etc.)
    """
    current_id = 0
    
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        for pattern_type, patterns in pattern_dict.items():
            for pattern_tuple in patterns:
                if pattern_tuple and current_id == pattern_id:
                    return {
                        'motif_class': motif_class,
                        'pattern_type': pattern_type,
                        'regex': pattern_tuple[0],
                        'original_id': pattern_tuple[1],
                        'group_number': pattern_tuple[2],
                        'subclass': pattern_tuple[3],
                        'score_scale': pattern_tuple[5] if len(pattern_tuple) > 5 else 1.0,
                        'min_runs': pattern_tuple[6] if len(pattern_tuple) > 6 else 0,
                        'min_score': pattern_tuple[7] if len(pattern_tuple) > 7 else 0.0,
                        'score_method': pattern_tuple[8] if len(pattern_tuple) > 8 else "default"
                    }
                current_id += 1
    
    return {}

def generate_cruciform_patterns(min_arm_len: int = 6, max_arm_len: int = 20, 
                               max_loop_len: int = 50) -> None:
    """
    Generate dynamic cruciform patterns for Hyperscan pre-filtering.
    
    Scientific basis: Lilley & Clegg (1993) Annu Rev Biophys Biomol Struct 22:299-328
    
    Args:
        min_arm_len: Minimum arm length for palindromes/inverted repeats (default: 6 bp)
        max_arm_len: Maximum arm length (default: 20 bp, limited by Hyperscan performance)
        max_loop_len: Maximum loop length for inverted repeats (default: 50 bp)
    """
    # Clear existing dynamic patterns
    CRUCIFORM_PATTERNS['palindrome_candidates'] = []
    CRUCIFORM_PATTERNS['inverted_repeat_candidates'] = []
    
    pattern_id = 200  # Start with ID 200 to avoid conflicts with static patterns
    
    # Generate palindrome pre-filter patterns (exact palindromes detected post-processing)
    for arm_len in range(min_arm_len, min(max_arm_len + 1, 21)):
        # Perfect palindrome candidates (no loop)
        pattern = f'[ACGT]{{{2*arm_len}}}'
        CRUCIFORM_PATTERNS['palindrome_candidates'].append(
            (pattern, pattern_id, 0, f"Palindrome_Candidate_{arm_len}bp", None, 1.0, 1, 0.5, "Palindrome_Thermodynamics")
        )
        pattern_id += 1
    
    # Generate inverted repeat pre-filter patterns (with variable loops)
    for arm_len in range(min_arm_len, min(max_arm_len + 1, 16)):  # Reduced for performance
        for loop_len in [5, 10, 15, 20, 30, 50]:  # Specific loop lengths instead of all
            if loop_len <= max_loop_len:
                total_len = 2 * arm_len + loop_len
                if total_len <= 100:  # Reasonable upper limit for Hyperscan performance
                    pattern = f'[ACGT]{{{total_len}}}'
                    CRUCIFORM_PATTERNS['inverted_repeat_candidates'].append(
                        (pattern, pattern_id, 0, f"Inverted_Repeat_{arm_len}bp_loop{loop_len}", None, 
                         1.0, 1, 0.5, "IR_Thermodynamics")
                    )
                    pattern_id += 1

def validate_patterns() -> bool:
    """
    Validate all regex patterns for syntax correctness.
    
    Returns:
        True if all patterns are valid, False otherwise
    """
    try:
        for motif_class, pattern_dict in ALL_PATTERNS.items():
            for pattern_type, patterns in pattern_dict.items():
                for pattern_tuple in patterns:
                    if pattern_tuple:
                        regex_pattern = pattern_tuple[0]
                        re.compile(regex_pattern)
        return True
    except re.error as e:
        print(f"Invalid regex pattern found: {e}")
        return False

# Initialize cruciform patterns
generate_cruciform_patterns()

# Validate all patterns on import
if not validate_patterns():
    print("Warning: Some regex patterns in registry are invalid!")

# Export commonly used pattern collections
HYPERSCAN_SAFE_PATTERNS = get_all_hyperscan_patterns()

__all__ = [
    'ALL_PATTERNS',
    'G_QUADRUPLEX_PATTERNS',
    'I_MOTIF_PATTERNS', 
    'Z_DNA_PATTERNS',
    'CURVED_DNA_PATTERNS',
    'TRIPLEX_PATTERNS',
    'CRUCIFORM_PATTERNS',
    'R_LOOP_PATTERNS',
    'SLIPPED_DNA_PATTERNS',
    'get_patterns_for_motif',
    'get_all_hyperscan_patterns',
    'get_pattern_info',
    'generate_cruciform_patterns',
    'validate_patterns',
    'HYPERSCAN_SAFE_PATTERNS'
]