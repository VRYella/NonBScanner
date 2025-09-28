"""
Modular Motif Detection System
=============================

Individual motif detector classes for each Non-B DNA motif type.
Split from centralized architecture for improved maintainability.
"""

from .base_detector import BaseMotifDetector
from .curved_dna_detector import CurvedDNADetector
from .slipped_dna_detector import SlippedDNADetector
from .cruciform_detector import CruciformDetector
from .r_loop_detector import RLoopDetector
from .triplex_detector import TriplexDetector
from .g_quadruplex_detector import GQuadruplexDetector
from .i_motif_detector import IMotifDetector
from .z_dna_detector import ZDNADetector
from .a_philic_detector import APhilicDetector

__all__ = [
    'BaseMotifDetector',
    'CurvedDNADetector',
    'SlippedDNADetector', 
    'CruciformDetector',
    'RLoopDetector',
    'TriplexDetector',
    'GQuadruplexDetector',
    'IMotifDetector',
    'ZDNADetector',
    'APhilicDetector'
]