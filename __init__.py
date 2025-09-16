"""
NBDFinder - Non-B DNA Structure Finder
=====================================

A comprehensive toolkit for detecting and analyzing non-B DNA structures 
in genomic sequences. Supports 10 major classes and 22 subclasses of 
non-B DNA motifs with scientifically validated scoring systems.

Main Features:
- Hyperscan-accelerated pattern matching
- Parallel motif detection
- Scientific scoring (G4Hunter, Z-DNA Seeker, etc.)
- Comprehensive visualization tools
- Multiple export formats (BED, CSV, Parquet, bedGraph)
- Command-line and web interfaces

Author: Dr. Venkata Rajesh Yella
"""

__version__ = "2.0.0"
__author__ = "Dr. Venkata Rajesh Yella"
__email__ = "vryella@example.com"

# Main API imports (optional, to avoid import loops)
try:
    from .orchestrators.all_motifs import detect_all_motifs
except ImportError:
    detect_all_motifs = None

from .core.regex_registry import ALL_PATTERNS as get_all_patterns, get_patterns_for_motif

try:
    from .nbdio.writers import export_to_bed, export_to_csv, export_to_parquet
except ImportError:
    export_to_bed = None
    export_to_csv = None  
    export_to_parquet = None

# Core classes
from .detectors import (
    Class01Curved, Class02Slipped, Class03Cruciform, Class04RLoop,
    Class05Triplex, Class06G4Family, Class07IMotif, Class08ZDna,
    Class09Hybrid, Class10Cluster
)

__all__ = [
    # Main API
    'detect_all_motifs',
    
    # Pattern registry
    'get_all_patterns',
    'get_patterns_for_motif',
    
    # Export functions
    'export_to_bed',
    'export_to_csv', 
    'export_to_parquet',
    
    # Detector classes
    'Class01Curved',
    'Class02Slipped',
    'Class03Cruciform',
    'Class04RLoop',
    'Class05Triplex',
    'Class06G4Family',
    'Class07IMotif',
    'Class08ZDna',
    'Class09Hybrid',
    'Class10Cluster',
    
    # Metadata
    '__version__',
    '__author__',
    '__email__'
]