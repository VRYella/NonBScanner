"""
NBDFinder Detector Modules
=========================

Individual detector modules for each of the 11 major Non-B DNA motif classes.
Each module provides specialized detection algorithms optimized for specific
structural features and biological contexts.

Detector Classes:
- Class 01: Curved DNA (A-tract mediated curvature)
- Class 02: Slipped DNA (Direct/tandem repeats forming slipped structures)
- Class 03: Cruciform (Inverted repeats forming four-way junctions)
- Class 04: R-Loop (RNA-DNA hybrids with displaced ssDNA)
- Class 05: Triplex (Three-stranded DNA structures)
- Class 06: G4 Family (G-quadruplex structures and variants)
- Class 07: i-Motif (C-rich structures complementary to G4)
- Class 08: Z-DNA (Left-handed double helix)
- Class 09: A-philic DNA (A-rich structural motifs)
- Class 10: Hybrid (Overlapping/composite motifs)
- Class 11: Cluster (High-density motif regions)
"""

# Import main detection functions from each class
try:
    from .class01_curved import find_curved_DNA as detect_curved_dna
except ImportError:
    detect_curved_dna = None

try:
    from .class02_slipped import find_slipped_dna as detect_slipped_dna
except ImportError:
    detect_slipped_dna = None

try:
    from .class03_cruciform import find_cruciform as detect_cruciform
except ImportError:
    detect_cruciform = None

try:
    from .class04_rloop import find_r_loop as detect_rloop
except ImportError:
    detect_rloop = None

try:
    from .class05_triplex import find_triplex as detect_triplex
except ImportError:
    detect_triplex = None

try:
    from .class06_g4_family import find_g_quadruplex as detect_g4_family
except ImportError:
    detect_g4_family = None

try:
    from .class07_imotif import find_i_motif as detect_imotif
except ImportError:
    detect_imotif = None

try:
    from .class08_zdna import find_z_dna as detect_zdna
except ImportError:
    detect_zdna = None

try:
    from .class09_a_philic import find_a_philic as detect_a_philic
except ImportError:
    detect_a_philic = None

try:
    from .class10_hybrid import find_hybrid as detect_hybrid
except ImportError:
    detect_hybrid = None

try:
    from .class11_cluster import find_cluster as detect_cluster
except ImportError:
    detect_cluster = None

# Detector registry mapping class names to functions
DETECTOR_REGISTRY = {
    'Curved_DNA': detect_curved_dna,           # Class 1
    'Slipped_DNA': detect_slipped_dna,         # Class 2
    'Cruciform': detect_cruciform,             # Class 3
    'R-Loop': detect_rloop,                    # Class 4
    'Triplex': detect_triplex,                 # Class 5
    'G-Quadruplex': detect_g4_family,          # Class 6
    'i-Motif': detect_imotif,                  # Class 7
    'Z-DNA': detect_zdna,                      # Class 8
    'A-philic_DNA': detect_a_philic,           # Class 9 (NEW)
    'Hybrid': detect_hybrid,                   # Class 10 (renumbered)
    'Cluster': detect_cluster,                 # Class 11 (renumbered)
}

# Filter out None values (failed imports)
DETECTOR_REGISTRY = {k: v for k, v in DETECTOR_REGISTRY.items() if v is not None}

def get_available_detectors() -> list:
    """Get list of available detector class names."""
    return list(DETECTOR_REGISTRY.keys())

def get_detector_function(motif_class: str):
    """Get detector function for a specific motif class."""
    return DETECTOR_REGISTRY.get(motif_class)

def detect_motif_class(sequence: str, motif_class: str, sequence_name: str = "sequence"):
    """
    Detect motifs for a specific class.
    
    Args:
        sequence: DNA sequence to analyze
        motif_class: Name of motif class to detect
        sequence_name: Name identifier for the sequence
        
    Returns:
        List of detected motifs for the class
    """
    detector_func = get_detector_function(motif_class)
    if detector_func is None:
        raise ValueError(f"No detector available for class: {motif_class}")
    
    try:
        return detector_func(sequence, sequence_name)
    except Exception as e:
        print(f"Error in {motif_class} detector: {e}")
        return []

# Class aliases for backward compatibility
Class01Curved = detect_curved_dna
Class02Slipped = detect_slipped_dna
Class03Cruciform = detect_cruciform
Class04RLoop = detect_rloop
Class05Triplex = detect_triplex
Class06G4Family = detect_g4_family
Class07IMotif = detect_imotif
Class08ZDna = detect_zdna
Class09APhilic = detect_a_philic
Class10Hybrid = detect_hybrid
Class11Cluster = detect_cluster

__all__ = [
    # Detection functions
    'detect_curved_dna',
    'detect_slipped_dna',
    'detect_cruciform',
    'detect_rloop',
    'detect_triplex',
    'detect_g4_family',
    'detect_imotif',
    'detect_zdna',
    'detect_a_philic',
    'detect_hybrid',
    'detect_cluster',
    
    # Registry and utilities
    'DETECTOR_REGISTRY',
    'get_available_detectors',
    'get_detector_function',
    'detect_motif_class',
    
    # Class aliases
    'Class01Curved',
    'Class02Slipped',
    'Class03Cruciform',
    'Class04RLoop',
    'Class05Triplex',
    'Class06G4Family',
    'Class07IMotif',
    'Class08ZDna',
    'Class09APhilic',
    'Class10Hybrid',
    'Class11Cluster',
]