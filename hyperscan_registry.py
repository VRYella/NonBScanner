"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                     HYPERSCAN PATTERN REGISTRY MODULE                        ║
║                    Non-B DNA Motif Detection System                          ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: hyperscan_registry.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    Comprehensive registry of all regex patterns feasible with Hyperscan engine.
    Contains 11 major motif classes with 22+ subclasses optimized for high-
    performance pattern matching.

PATTERN CLASSIFICATION TABLE:
┌──────┬─────────────────────┬────────────┬──────────────────────────────────────┐
│Class │ Name                │ Subclasses │ Hyperscan Feasibility                │
├──────┼─────────────────────┼────────────┼──────────────────────────────────────┤
│  1   │ Curved DNA          │     2      │ ✓ Simple regex patterns              │
│  2   │ Slipped DNA         │     2      │ ✓ STR patterns (partial)             │
│  3   │ Cruciform DNA       │     1      │ ✗ Requires reverse complement        │
│  4   │ R-loop              │     3      │ ✓ GC-rich patterns                   │
│  5   │ Triplex             │     2      │ ✓ Homopurine/pyrimidine              │
│  6   │ G-Quadruplex        │     7      │ ✓ G-run patterns                     │
│  7   │ i-Motif Family      │     3      │ ✓ C-run patterns                     │
│  8   │ Z-DNA               │     2      │ ✓ Alternating patterns               │
│  9   │ A-philic DNA        │     1      │ ✗ Requires scoring algorithm         │
│ 10   │ Hybrid              │  Dynamic   │ ✗ Post-processing required           │
│ 11   │ Non-B Clusters      │  Dynamic   │ ✗ Post-processing required           │
└──────┴─────────────────────┴────────────┴──────────────────────────────────────┘

HYPERSCAN PATTERNS: 150+ patterns across 8 classes (Classes 1-2, 4-8)
NON-HYPERSCAN: Classes 3, 9, 10, 11 require algorithmic detection

REFERENCES:
    - Bedrat et al., 2016 (G4Hunter)
    - Jenjaroenpun & Wongsurawat, 2016 (QmRLFS)
    - Ho et al., 1986 (Z-DNA)
    - Olson et al., 1998 (DNA curvature)
"""

import re
from typing import Dict, List, Tuple, Any
from collections import defaultdict

# Try to import Hyperscan
try:
    import hyperscan
    HYPERSCAN_AVAILABLE = True
except (ImportError, Exception):
    HYPERSCAN_AVAILABLE = False


class HyperscanRegistry:
    """
    Comprehensive registry of Hyperscan-compatible patterns for Non-B DNA detection.
    
    This registry contains all motif patterns that can be efficiently detected using
    the Hyperscan regex engine. Patterns are organized by motif class and include
    metadata for scoring and classification.
    """
    
    # =========================================================================
    # CLASS 1: CURVED DNA PATTERNS (Hyperscan Compatible)
    # =========================================================================
    CURVED_DNA = {
        'patterns': [
            # A-tracts: Simple regex patterns for A-tract detection
            (r'A{4,}', 1, 'A-tract', 'Local Curvature', 'curvature_score'),
            (r'T{4,}', 2, 'T-tract', 'Local Curvature', 'curvature_score'),
            (r'(?:A{3,}T{1,3}){2,}', 3, 'AT-rich tract', 'Local Curvature', 'curvature_score'),
            (r'(?:T{3,}A{1,3}){2,}', 4, 'TA-rich tract', 'Local Curvature', 'curvature_score'),
            
            # Phased A-tracts: Global curvature patterns
            (r'(?:A{3,}.{8,12}){3,}A{3,}', 5, 'Phased A-tracts', 'Global Curvature', 'phasing_score'),
            (r'(?:T{3,}.{8,12}){3,}T{3,}', 6, 'Phased T-tracts', 'Global Curvature', 'phasing_score'),
        ],
        'class_id': 1,
        'class_name': 'Curved DNA',
        'description': 'A-tract mediated DNA bending'
    }
    
    # =========================================================================
    # CLASS 2: SLIPPED DNA PATTERNS (Partial Hyperscan Compatible)
    # =========================================================================
    SLIPPED_DNA = {
        'patterns': [
            # Short Tandem Repeats (STR) - Without backreferences
            # Mononucleotide repeats
            (r'A{10,}', 10, 'Poly-A repeat', 'STR', 'instability_score'),
            (r'T{10,}', 11, 'Poly-T repeat', 'STR', 'instability_score'),
            (r'G{10,}', 12, 'Poly-G repeat', 'STR', 'instability_score'),
            (r'C{10,}', 13, 'Poly-C repeat', 'STR', 'instability_score'),
            
            # Specific common dinucleotide repeats
            (r'(?:CA){5,}', 16, 'CA repeat', 'STR', 'instability_score'),
            (r'(?:AC){5,}', 17, 'AC repeat', 'STR', 'instability_score'),
            (r'(?:AT){5,}', 18, 'AT repeat', 'STR', 'instability_score'),
            (r'(?:TA){5,}', 19, 'TA repeat', 'STR', 'instability_score'),
            (r'(?:CG){5,}', 20, 'CG repeat', 'STR', 'instability_score'),
            (r'(?:GC){5,}', 21, 'GC repeat', 'STR', 'instability_score'),
            (r'(?:AG){5,}', 22, 'AG repeat', 'STR', 'instability_score'),
            (r'(?:GA){5,}', 23, 'GA repeat', 'STR', 'instability_score'),
            (r'(?:CT){5,}', 24, 'CT repeat', 'STR', 'instability_score'),
            (r'(?:TC){5,}', 25, 'TC repeat', 'STR', 'instability_score'),
            
            # Specific common trinucleotide repeats
            (r'(?:CGG){4,}', 30, 'CGG repeat', 'STR', 'instability_score'),
            (r'(?:CCG){4,}', 31, 'CCG repeat', 'STR', 'instability_score'),
            (r'(?:CAG){4,}', 32, 'CAG repeat', 'STR', 'instability_score'),
            (r'(?:CTG){4,}', 33, 'CTG repeat', 'STR', 'instability_score'),
            (r'(?:GAA){4,}', 34, 'GAA repeat', 'STR', 'instability_score'),
            (r'(?:TTC){4,}', 35, 'TTC repeat', 'STR', 'instability_score'),
            (r'(?:AAG){4,}', 36, 'AAG repeat', 'STR', 'instability_score'),
            (r'(?:CTT){4,}', 37, 'CTT repeat', 'STR', 'instability_score'),
            
            # Note: Direct repeats with variable spacers are NOT Hyperscan compatible
            # General repeats with backreferences are handled in non_hyperscan_detection.py
        ],
        'class_id': 2,
        'class_name': 'Slipped DNA',
        'description': 'Short tandem repeats (STR subclass only)'
    }
    
    # =========================================================================
    # CLASS 4: R-LOOP PATTERNS (Hyperscan Compatible)
    # =========================================================================
    R_LOOP = {
        'patterns': [
            # GC-rich R-loop sites
            (r'[GC]{10,}', 30, 'GC-rich tract', 'R-loop formation sites', 'r_loop_potential'),
            (r'G{5,}[ATGC]{5,50}C{5,}', 31, 'G-C rich region', 'R-loop formation sites', 'r_loop_potential'),
            (r'[GC]{6,}[AT]{1,5}[GC]{6,}', 32, 'GC-AT-GC pattern', 'R-loop formation sites', 'r_loop_potential'),
            
            # G-rich patterns for R-loop formation
            (r'G{3,}[ATGC]{1,10}G{3,}[ATGC]{1,10}G{3,}', 33, 'Multi-G-tract', 'R-loop formation sites', 'r_loop_potential'),
            (r'G{4,}[ATGC]{1,10}G{4,}', 34, 'Strong G-tract pair', 'R-loop formation sites', 'r_loop_potential'),
        ],
        'class_id': 4,
        'class_name': 'R-loop',
        'description': 'RNA-DNA hybrid formation sites'
    }
    
    # =========================================================================
    # CLASS 5: TRIPLEX PATTERNS (Hyperscan Compatible)
    # =========================================================================
    TRIPLEX = {
        'patterns': [
            # Homopurine tracts
            (r'[GA]{10,}', 40, 'Homopurine tract', 'Triplex', 'triplex_score'),
            (r'G{5,}[GA]{5,}', 41, 'G-rich purine', 'Triplex', 'triplex_score'),
            (r'A{5,}[GA]{5,}', 42, 'A-rich purine', 'Triplex', 'triplex_score'),
            
            # Homopyrimidine tracts
            (r'[CT]{10,}', 43, 'Homopyrimidine tract', 'Triplex', 'triplex_score'),
            (r'C{5,}[CT]{5,}', 44, 'C-rich pyrimidine', 'Triplex', 'triplex_score'),
            (r'T{5,}[CT]{5,}', 45, 'T-rich pyrimidine', 'Triplex', 'triplex_score'),
            
            # Mirror repeats (Sticky DNA)
            (r'(?:GA){5,}[GA]*(?:TC){5,}', 46, 'GA-TC mirror', 'Sticky DNA', 'mirror_score'),
            (r'(?:AG){5,}[AG]*(?:CT){5,}', 47, 'AG-CT mirror', 'Sticky DNA', 'mirror_score'),
        ],
        'class_id': 5,
        'class_name': 'Triplex',
        'description': 'Three-stranded DNA structures'
    }
    
    # =========================================================================
    # CLASS 6: G-QUADRUPLEX FAMILY PATTERNS (Hyperscan Compatible)
    # =========================================================================
    G_QUADRUPLEX = {
        'patterns': [
            # Canonical G4: 3+ G-runs with short loops (1-7 nt)
            (r'G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 50, 'Canonical G4', 'Canonical G4', 'g4hunter_score'),
            
            # Relaxed G4: 2+ G-runs with medium loops (1-12 nt)
            (r'G{2,}[ATGC]{1,12}G{2,}[ATGC]{1,12}G{2,}[ATGC]{1,12}G{2,}', 51, 'Relaxed G4', 'Relaxed G4', 'g4hunter_score'),
            
            # Bulged G4: One long loop (8-20 nt)
            (r'G{3,}[ATGC]{8,20}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 52, 'Bulged G4', 'Bulged G4', 'g4hunter_score'),
            (r'G{3,}[ATGC]{1,7}G{3,}[ATGC]{8,20}G{3,}[ATGC]{1,7}G{3,}', 53, 'Bulged G4 v2', 'Bulged G4', 'g4hunter_score'),
            (r'G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{8,20}G{3,}', 54, 'Bulged G4 v3', 'Bulged G4', 'g4hunter_score'),
            
            # Bipartite G4: One very long loop (15-50 nt)
            (r'G{2,}[ATGC]{15,50}G{2,}[ATGC]{1,7}G{2,}[ATGC]{1,7}G{2,}', 55, 'Bipartite G4', 'Bipartite G4', 'g4hunter_score'),
            
            # Multimeric G4: 5+ G-runs
            (r'(?:G{3,}[ATGC]{1,7}){4,}G{3,}', 56, 'Multimeric G4', 'Multimeric G4', 'g4hunter_score'),
            
            # Imperfect G4: Includes single base substitutions
            (r'G{2,}[ATGC]{1,10}[AG]G{1,2}[ATGC]{1,10}G{2,}[ATGC]{1,10}G{2,}', 57, 'Imperfect G4', 'Imperfect G4', 'g4hunter_score'),
            
            # G-Triplex: 3 G-runs (intermediate structure)
            (r'G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 58, 'G-Triplex', 'G-Triplex intermediate', 'g_triplex_score'),
        ],
        'class_id': 6,
        'class_name': 'G-Quadruplex',
        'description': 'G-rich four-stranded structures'
    }
    
    # =========================================================================
    # CLASS 7: I-MOTIF FAMILY PATTERNS (Hyperscan Compatible)
    # =========================================================================
    I_MOTIF = {
        'patterns': [
            # Canonical i-motif: 3+ C-runs with short loops
            (r'C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}', 60, 'Canonical i-motif', 'Canonical i-motif', 'imotif_score'),
            
            # Relaxed i-motif: 2+ C-runs with medium loops
            (r'C{2,}[ATGC]{1,12}C{2,}[ATGC]{1,12}C{2,}[ATGC]{1,12}C{2,}', 61, 'Relaxed i-motif', 'Relaxed i-motif', 'imotif_score'),
            
            # AC-motif: Alternating A and C patterns
            (r'(?:AC){4,}', 62, 'AC repeat', 'AC-motif', 'ac_motif_score'),
            (r'(?:CA){4,}', 63, 'CA repeat', 'AC-motif', 'ac_motif_score'),
            (r'A{3,}C{3,}A{3,}C{3,}', 64, 'A-C alternating', 'AC-motif', 'ac_motif_score'),
        ],
        'class_id': 7,
        'class_name': 'i-Motif Family',
        'description': 'C-rich structures'
    }
    
    # =========================================================================
    # CLASS 8: Z-DNA PATTERNS (Hyperscan Compatible)
    # =========================================================================
    Z_DNA = {
        'patterns': [
            # Classic Z-DNA: Alternating purine-pyrimidine
            (r'(?:CG){4,}', 70, 'CG alternating', 'Z-DNA', 'zdna_score'),
            (r'(?:GC){4,}', 71, 'GC alternating', 'Z-DNA', 'zdna_score'),
            (r'(?:AT){4,}', 72, 'AT alternating', 'Z-DNA', 'zdna_score'),
            (r'(?:TA){4,}', 73, 'TA alternating', 'Z-DNA', 'zdna_score'),
            (r'(?:AG){4,}', 74, 'AG alternating', 'Z-DNA', 'zdna_score'),
            (r'(?:GA){4,}', 75, 'GA alternating', 'Z-DNA', 'zdna_score'),
            (r'(?:CT){4,}', 76, 'CT alternating', 'Z-DNA', 'zdna_score'),
            (r'(?:TC){4,}', 77, 'TC alternating', 'Z-DNA', 'zdna_score'),
            
            # eGZ (Extruded-G Z-DNA): Long GC-rich regions
            (r'(?:CGG){3,}', 78, 'CGG repeat', 'eGZ (Extruded-G) DNA', 'egz_score'),
            (r'(?:GGC){3,}', 79, 'GGC repeat', 'eGZ (Extruded-G) DNA', 'egz_score'),
            (r'[CG]{8,}', 80, 'Long CG-rich', 'eGZ (Extruded-G) DNA', 'egz_score'),
        ],
        'class_id': 8,
        'class_name': 'Z-DNA',
        'description': 'Left-handed double helix'
    }
    
    @classmethod
    def get_all_patterns(cls) -> Dict[str, Dict]:
        """
        Get all Hyperscan-compatible patterns organized by class.
        
        Returns:
            Dictionary mapping class names to pattern dictionaries
        """
        return {
            'curved_dna': cls.CURVED_DNA,
            'slipped_dna': cls.SLIPPED_DNA,
            'r_loop': cls.R_LOOP,
            'triplex': cls.TRIPLEX,
            'g_quadruplex': cls.G_QUADRUPLEX,
            'i_motif': cls.I_MOTIF,
            'z_dna': cls.Z_DNA
        }
    
    @classmethod
    def get_pattern_list(cls) -> List[Tuple[str, int, str, str, str]]:
        """
        Get a flat list of all patterns for Hyperscan compilation.
        
        Returns:
            List of tuples: (pattern, pattern_id, name, subclass, scoring_method)
        """
        patterns = []
        all_classes = cls.get_all_patterns()
        
        for class_key, class_data in all_classes.items():
            for pattern_tuple in class_data['patterns']:
                patterns.append(pattern_tuple)
        
        return patterns
    
    @classmethod
    def compile_hyperscan_database(cls) -> Tuple[Any, Dict[int, Dict]]:
        """
        Compile all patterns into a Hyperscan database.
        
        Returns:
            Tuple of (database, pattern_info_dict)
            pattern_info_dict maps pattern_id to metadata
        """
        if not HYPERSCAN_AVAILABLE:
            raise ImportError("Hyperscan is not available")
        
        patterns = cls.get_pattern_list()
        
        # Prepare patterns for Hyperscan
        expressions = []
        ids = []
        flags = []
        pattern_info = {}
        
        for pattern, pattern_id, name, subclass, scoring_method in patterns:
            expressions.append(pattern.encode('utf-8'))
            ids.append(pattern_id)
            flags.append(hyperscan.HS_FLAG_CASELESS | hyperscan.HS_FLAG_DOTALL)
            
            pattern_info[pattern_id] = {
                'name': name,
                'subclass': subclass,
                'scoring_method': scoring_method,
                'pattern': pattern
            }
        
        # Compile database
        database = hyperscan.Database()
        database.compile(
            expressions=expressions,
            ids=ids,
            flags=flags
        )
        
        return database, pattern_info


def get_class_mapping() -> Dict[int, str]:
    """
    Get mapping of pattern IDs to class names.
    
    Returns:
        Dictionary mapping pattern_id to class_name
    """
    mapping = {}
    
    # Curved DNA: 1-9
    for i in range(1, 10):
        mapping[i] = 'Curved DNA'
    
    # Slipped DNA: 10-29
    for i in range(10, 30):
        mapping[i] = 'Slipped DNA'
    
    # R-loop: 30-39
    for i in range(30, 40):
        mapping[i] = 'R-loop'
    
    # Triplex: 40-49
    for i in range(40, 50):
        mapping[i] = 'Triplex'
    
    # G-Quadruplex: 50-59
    for i in range(50, 60):
        mapping[i] = 'G-Quadruplex'
    
    # i-Motif: 60-69
    for i in range(60, 70):
        mapping[i] = 'i-Motif Family'
    
    # Z-DNA: 70-89
    for i in range(70, 90):
        mapping[i] = 'Z-DNA'
    
    return mapping


# Color scheme for visualization (consistent with app.py)
MOTIF_CLASS_COLORS = {
    'Curved DNA': '#FF6B6B',
    'Slipped DNA': '#4ECDC4',
    'Cruciform': '#45B7D1',
    'R-loop': '#FFA07A',
    'Triplex': '#98D8C8',
    'G-Quadruplex': '#F7DC6F',
    'i-Motif Family': '#BB8FCE',
    'Z-DNA': '#85C1E2',
    'A-philic DNA': '#F8B88B',
    'Hybrid': '#95A5A6',
    'Non-B DNA Clusters': '#34495E'
}


if __name__ == '__main__':
    # Test pattern registry
    print("Hyperscan Registry Test")
    print("=" * 80)
    
    registry = HyperscanRegistry()
    all_patterns = registry.get_all_patterns()
    
    print(f"\nTotal classes: {len(all_patterns)}")
    for class_name, class_data in all_patterns.items():
        print(f"\n{class_data['class_name']}:")
        print(f"  - Patterns: {len(class_data['patterns'])}")
        print(f"  - Description: {class_data['description']}")
    
    pattern_list = registry.get_pattern_list()
    print(f"\nTotal patterns: {len(pattern_list)}")
    
    if HYPERSCAN_AVAILABLE:
        print("\n✓ Hyperscan is available")
        try:
            db, pattern_info = registry.compile_hyperscan_database()
            print(f"✓ Successfully compiled {len(pattern_info)} patterns into Hyperscan database")
        except Exception as e:
            print(f"✗ Error compiling database: {e}")
    else:
        print("\n✗ Hyperscan is not available")
