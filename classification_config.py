"""
Classification Configuration for NBDFinder
==========================================

This module defines scoring systems, length constraints, and normalization 
parameters for all Non-B DNA motif classes based on rigorous literature survey.

Scientific Annotations:
- Length constraints based on biological and computational expectations
- Values practical for typical eukaryotic genomes
- Scoring methods validated against experimental data
"""

import numpy as np
from typing import Dict, Tuple, Any

# Motif length constraints (S_min/S_max) table for NBDFinder
# Based on biological and code-derived expectations
# Scientific annotations per row; values practical for typical eukaryotic genomes
MOTIF_LENGTH_LIMITS = {
    "Curved_DNA":        {"S_min": 15,  "S_max": 200},   # A/T tracts ≥7bp; arrays up to ~100bp
    "Z-DNA":             {"S_min": 50,  "S_max": 200},   # Z-score threshold, long GC runs
    "eGZ":               {"S_min": 28,  "S_max": 200},   # (CGG)4 = 12bp, practical upper ~30x
    "Slipped_DNA_DR":    {"S_min": 20,  "S_max": 300},   # 2x10bp DR, up to 2x100bp
    "Slipped_DNA_STR":   {"S_min": 11,  "S_max": 100},   # 5x1bp, up to 50x2bp
    "R-Loop":            {"S_min": 140, "S_max": 400},   # RLFS+REZ, min 100bp, up to 400bp
    "Cruciform":         {"S_min": 23,  "S_max": 100},   # 2x10bp arm + 3bp spacer, max 2x50bp
    "Triplex":           {"S_min": 30,  "S_max": 400},   # 10bp arm + spacer, up to 2x100bp
    "A-philic_DNA":      {"S_min": 10,  "S_max": 50},    # Tetranucleotide windows, 10-20bp optimal
    "Sticky_DNA":        {"S_min": 236, "S_max": 1000},  # 59x GAA, up to ~333x
    "G4":                {"S_min": 13,  "S_max": 100},   # 4x3bp G run, up to long G4s
    "G-Triplex":         {"S_min": 28,  "S_max": 100},   # 3x3bp G runs
    "i-Motif":           {"S_min": 23,  "S_max": 100},   # 4x3bp C runs
    "AC-motif":          {"S_min": 21,  "S_max": 37},    # Shortest/longest consensus
}

# Scoring method configurations for each motif class
SCORING_METHODS = {
    "Curved_DNA": {
        "method": "curvature_raw",
        "description": "AT content weighted by tract length and spacing",
        "reference": "Olson et al. PNAS 1998",
        "base_score_range": (5.0, 200.0),
        "normalization_factor": "length_adjusted"
    },
    "Z-DNA": {
        "method": "z_seeker",
        "description": "Dinucleotide scoring with mismatch penalties",
        "reference": "Ho et al. Nucleic Acids Research 1986",
        "base_score_range": (50.0, 500.0),
        "normalization_factor": "gc_weighted"
    },
    "eGZ": {
        "method": "repeat_count_raw",
        "description": "CGG repeat number with G-bias weighting",
        "reference": "Usdin et al. Nature Reviews Genetics 2015",
        "base_score_range": (12.0, 300.0),
        "normalization_factor": "repeat_based"
    },
    "Curved_DNA": {
        "method": "enhanced_curvature",
        "description": "Dinucleotide bending angles with AT-tract periodicity",
        "reference": "Olson et al. PNAS 1998, Crothers model 1990",
        "base_score_range": (5.0, 200.0),
        "normalization_factor": "curvature_weighted"
    },
    "Slipped_DNA": {
        "method": "instability_based",
        "description": "Repeat instability scoring with AT-content and unit-specific factors",
        "reference": "McMurray Cell 2010, Ellegren Nature Reviews Genetics 2004",
        "base_score_range": (10.0, 600.0),
        "normalization_factor": "instability_weighted"
    },
    "Cruciform": {
        "method": "thermodynamic",
        "description": "Thermodynamic stability based on bp energy and loop penalties",
        "reference": "Lilley & Kemper Nature 1984, SantaLucia PNAS 1998",
        "base_score_range": (8.0, 150.0),
        "normalization_factor": "thermo_stability"
    },
    "Triplex": {
        "method": "triplex_enhanced",
        "description": "Homopurine/pyrimidine content with pH-dependent stability",
        "reference": "Frank-Kamenetskii & Mirkin Annual Review of Biochemistry 1995",
        "base_score_range": (30.0, 800.0),
        "normalization_factor": "homogeneity_weighted"
    },
    "A-philic_DNA": {
        "method": "tetranucleotide_log2_odds",
        "description": "Tetranucleotide log2 odds scoring with strong count thresholding",
        "reference": "Vinogradov Bioinformatics 2003, Bolshoy et al. PNAS 1991",
        "base_score_range": (10.0, 50.0),
        "normalization_factor": "log2_odds_weighted"
    },
    "Sticky_DNA": {
        "method": "sticky_enhanced", 
        "description": "GAA/TTC triplex potential with superlinear repeat scaling",
        "reference": "Sakamoto et al. Journal of Biological Chemistry 1999",
        "base_score_range": (200.0, 2000.0),
        "normalization_factor": "triplex_potential"
    },
    "i-Motif": {
        "method": "g4hunter_adapted",
        "description": "G4Hunter algorithm adapted for C-richness (+1 for C, -1 for G)",
        "reference": "Adapted from Bedrat et al. Nucleic Acids Research 2016",
        "base_score_range": (1.0, 180.0),
        "normalization_factor": "imotif_hunter_score"
    },
    "AC-motif": {
        "method": "g4hunter_adapted",
        "description": "G4Hunter-style scoring with C-richness and A-tract weighting",
        "reference": "Adapted from Bedrat et al. Nucleic Acids Research 2016",
        "base_score_range": (5.0, 100.0),
        "normalization_factor": "ac_motif_score"
    }
}

# Conservation analysis parameters
CONSERVATION_CONFIG = {
    "n_shuffles": 100,
    "pseudocount": 1e-6,
    "kmer_size": 6,
    "min_enrichment_threshold": 1.5,
    "significant_pvalue_threshold": 0.05
}

def get_motif_limits(motif_class: str, subclass: str = None) -> Tuple[int, int]:
    """
    Get length limits for a motif class/subclass.
    
    Args:
        motif_class: Main motif class
        subclass: Motif subclass (if applicable)
    
    Returns:
        Tuple of (S_min, S_max) length limits
    """
    # Handle special cases for subclasses
    if motif_class == "Z-DNA" and subclass == "eGZ (Extruded-G)":
        limits = MOTIF_LENGTH_LIMITS.get("eGZ", MOTIF_LENGTH_LIMITS["Z-DNA"])
    elif motif_class == "Slipped_DNA":
        if subclass and "STR" in subclass:
            limits = MOTIF_LENGTH_LIMITS.get("Slipped_DNA_STR", MOTIF_LENGTH_LIMITS["Slipped_DNA_DR"])
        else:
            limits = MOTIF_LENGTH_LIMITS.get("Slipped_DNA_DR", {"S_min": 20, "S_max": 300})
    elif motif_class in ["G4", "Relaxed_G4", "Bulged_G4", "Bipartite_G4", "Multimeric_G4"]:
        limits = MOTIF_LENGTH_LIMITS.get("G4", {"S_min": 13, "S_max": 100})
    elif motif_class == "AC-Motif":
        limits = MOTIF_LENGTH_LIMITS.get("AC-motif", {"S_min": 21, "S_max": 37})
    elif motif_class == "A-philic DNA":
        limits = MOTIF_LENGTH_LIMITS.get("A-philic_DNA", {"S_min": 10, "S_max": 50})
    else:
        # Try exact match first
        limits = MOTIF_LENGTH_LIMITS.get(motif_class)
        if limits is None:
            # Fallback to generic limits
            limits = {"S_min": 10, "S_max": 200}
    
    return limits["S_min"], limits["S_max"]

def normalize_score(actual_score: float, motif_length: int, motif_class: str, 
                   subclass: str = None) -> float:
    """
    Normalize motif score based on theoretical minimum to maximum score range.
    
    Args:
        actual_score: Raw motif score
        motif_length: Length of the motif
        motif_class: Main motif class
        subclass: Motif subclass
    
    Returns:
        Normalized score (0-1 scale) where 0=low stability, 1=high stability
    """
    s_min, s_max = get_motif_limits(motif_class, subclass)
    
    # Get scoring method info
    score_info = SCORING_METHODS.get(motif_class, SCORING_METHODS.get("Curved_DNA"))
    theoretical_min, theoretical_max = score_info["base_score_range"]
    
    # Adjust theoretical range based on motif length
    # Longer motifs generally have higher potential scores
    length_factor = max(0.5, min(2.0, motif_length / ((s_min + s_max) / 2)))
    adjusted_max = theoretical_max * length_factor
    adjusted_min = theoretical_min  # Keep minimum as baseline
    
    # Normalize to 0-1 scale using min-to-max normalization
    if adjusted_max <= adjusted_min:
        return 0.5  # Default for edge cases
    
    normalized = (actual_score - adjusted_min) / (adjusted_max - adjusted_min)
    
    # Ensure score is within bounds
    normalized = max(0.0, min(1.0, normalized))
    
    return round(normalized, 3)

def calculate_enrichment_score(observed_count: int, shuffled_counts: list, 
                             pseudocount: float = None) -> Tuple[float, float]:
    """
    Calculate enrichment score and p-value from shuffled distribution.
    
    Args:
        observed_count: Observed motif count
        shuffled_counts: List of counts from shuffled sequences
        pseudocount: Small value to avoid log(0)
    
    Returns:
        Tuple of (log2_enrichment, p_value)
    """
    if pseudocount is None:
        pseudocount = CONSERVATION_CONFIG["pseudocount"]
    
    if not shuffled_counts:
        return 0.0, 1.0
    
    mean_shuffled = np.mean(shuffled_counts)
    
    # Calculate log2 enrichment with pseudocounts
    log2_enrichment = np.log2((observed_count + pseudocount) / (mean_shuffled + pseudocount))
    
    # Calculate p-value as fraction of shuffled counts >= observed
    p_value = sum(1 for count in shuffled_counts if count >= observed_count) / len(shuffled_counts)
    
    return round(log2_enrichment, 3), round(p_value, 4)

def classify_conservation(log2_enrichment: float, p_value: float) -> str:
    """
    Classify motif conservation based on enrichment and significance.
    
    Args:
        log2_enrichment: Log2 enrichment score
        p_value: Statistical p-value
    
    Returns:
        Conservation classification string
    """
    threshold_enrich = CONSERVATION_CONFIG["min_enrichment_threshold"]
    threshold_pval = CONSERVATION_CONFIG["significant_pvalue_threshold"]
    
    if log2_enrichment >= np.log2(threshold_enrich) and p_value <= threshold_pval:
        return "Highly_Conserved"
    elif log2_enrichment >= np.log2(threshold_enrich * 0.75) and p_value <= threshold_pval * 2:
        return "Moderately_Conserved" 
    elif log2_enrichment > 0:
        return "Weakly_Conserved"
    elif log2_enrichment < -np.log2(threshold_enrich):
        return "Depleted"
    else:
        return "Neutral"

# Export main configuration elements
__all__ = [
    'MOTIF_LENGTH_LIMITS',
    'SCORING_METHODS', 
    'CONSERVATION_CONFIG',
    'get_motif_limits',
    'normalize_score',
    'calculate_enrichment_score',
    'classify_conservation'
]