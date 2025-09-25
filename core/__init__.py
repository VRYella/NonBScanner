"""
NBDFinder Core Module
===================

Core functionality for NBDFinder including:
- Regex pattern registry
- Hyperscan database management
- Vectorized scoring algorithms
- Sequence windowing and chunking
- Post-processing pipelines
- Utility functions
"""

from .regex_registry import (
    get_patterns_for_motif, get_pattern_info,
    get_all_hyperscan_patterns, ALL_PATTERNS
)

from .hyperscan_manager import HyperscanManager

from .scoring_simd import (
    g4hunter_score_vectorized, imotif_score_vectorized, 
    zdna_score_vectorized, score_sequence_region, batch_score_regions
)

from .windows import (
    SequenceWindow, SequenceWindower, GenomeWindower,
    calculate_window_statistics
)

from .postprocess import (
    remove_overlapping_motifs, merge_nearby_motifs, deduplicate_motifs,
    filter_by_score_threshold, filter_by_length_constraints,
    calculate_motif_statistics, apply_all_postprocessing
)

from .hs_dispatcher import (
    HyperscanDispatcher, get_global_dispatcher, register_all_detectors
)

from .utils import (
    parse_fasta, wrap, gc_content, reverse_complement, is_palindrome,
    calculate_tm, shuffle_sequence, get_basic_stats
)

__all__ = [
    # Pattern registry
    'ALL_PATTERNS',
    'get_patterns_for_motif', 
    'get_pattern_info',
    'get_all_hyperscan_patterns',
    
    # Hyperscan management
    'HyperscanManager',
    'HyperscanDispatcher',
    'get_global_dispatcher',
    'register_all_detectors',
    
    # Scoring algorithms
    'g4hunter_score_vectorized',
    'imotif_score_vectorized',
    'zdna_score_vectorized',
    'score_sequence_region',
    'batch_score_regions',
    
    # Windowing
    'SequenceWindow',
    'SequenceWindower',
    'GenomeWindower',
    'calculate_window_statistics',
    
    # Post-processing
    'remove_overlapping_motifs',
    'merge_nearby_motifs',
    'deduplicate_motifs',
    'filter_by_score_threshold',
    'filter_by_length_constraints',
    'calculate_motif_statistics',
    'apply_all_postprocessing',
    
    # Utilities
    'parse_fasta',
    'wrap',
    'gc_content',
    'reverse_complement',
    'is_palindrome',
    'calculate_tm',
    'shuffle_sequence',
    'get_basic_stats'
]