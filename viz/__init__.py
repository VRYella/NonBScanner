"""
NBDFinder Visualization Module
=============================

Visualization tools for NBDFinder results including:
- Static and interactive plots
- Genome browser track generation
- Track hubs for UCSC, IGV, JBrowse
"""

# Import visualization functions
try:
    from .plots import (
        generate_pseudodata, plot_motif_distribution, plot_motif_lengths,
        plot_score_distribution, create_motif_heatmap, create_interactive_plots
    )
except ImportError:
    # Fallback if plots module has import issues
    generate_pseudodata = None
    plot_motif_distribution = None
    plot_motif_lengths = None
    plot_score_distribution = None
    create_motif_heatmap = None
    create_interactive_plots = None

from .browser import TrackGenerator

__all__ = [
    # Plotting functions
    'generate_pseudodata',
    'plot_motif_distribution',
    'plot_motif_lengths', 
    'plot_score_distribution',
    'create_motif_heatmap',
    'create_interactive_plots',
    
    # Browser integration
    'TrackGenerator'
]