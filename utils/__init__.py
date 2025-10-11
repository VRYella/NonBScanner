"""
Utils Package - Helper Functions and Scanner Modules
====================================================

This package contains utility modules and scanner implementations for NBDScanner.

Modules:
- utils: General utility functions for sequence processing, I/O, and statistics
- visualization: Plotting and visualization functions
- advanced_visualizations: Advanced plotting capabilities
- motif_patterns: Pattern definitions and registry
- modular_scanner: Production scanner implementation
- nbdscanner: Legacy scanner wrapper with backward compatibility
"""

from .utils import (
    parse_fasta, gc_content, reverse_complement, wrap,
    get_basic_stats, export_to_bed, export_to_csv, export_to_json,
    validate_sequence, quality_check_motifs
)

from .visualization import (
    plot_motif_distribution, plot_coverage_map, plot_score_distribution,
    plot_length_distribution, plot_nested_pie_chart, save_all_plots,
    MOTIF_CLASS_COLORS
)

from .nbdscanner import (
    analyze_sequence, analyze_multiple_sequences,
    get_motif_classification_info, export_results_to_dataframe
)

__all__ = [
    # From utils
    'parse_fasta', 'gc_content', 'reverse_complement', 'wrap',
    'get_basic_stats', 'export_to_bed', 'export_to_csv', 'export_to_json',
    'validate_sequence', 'quality_check_motifs',
    # From visualization
    'plot_motif_distribution', 'plot_coverage_map', 'plot_score_distribution',
    'plot_length_distribution', 'plot_nested_pie_chart', 'save_all_plots',
    'MOTIF_CLASS_COLORS',
    # From nbdscanner
    'analyze_sequence', 'analyze_multiple_sequences',
    'get_motif_classification_info', 'export_results_to_dataframe'
]
