"""
NBDScanner Advanced Visualization Suite
========================================

Advanced static plotting functions for publication-quality motif analysis.
Implements 14 diverse visualization types with colorblind-friendly palettes.

Author: Dr. Venkata Rajesh Yella
License: MIT
Version: 2024.1
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.collections import LineCollection
import seaborn as sns
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Tuple
from collections import Counter, defaultdict
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# COLORBLIND-FRIENDLY PALETTES
# =============================================================================

# Okabe-Ito colorblind-safe palette (optimized for scientific publications)
COLORBLIND_PALETTE = [
    '#E69F00',  # Orange
    '#56B4E9',  # Sky Blue
    '#009E73',  # Bluish Green
    '#F0E442',  # Yellow
    '#0072B2',  # Blue
    '#D55E00',  # Vermillion
    '#CC79A7',  # Reddish Purple
]

# Class-specific colors using colorblind-safe palette
MOTIF_CLASS_COLORS_CB = {
    'Curved_DNA': '#E69F00',           # Orange
    'Slipped_DNA': '#F0E442',          # Yellow
    'Cruciform': '#56B4E9',            # Sky Blue
    'R-Loop': '#009E73',               # Bluish Green
    'Triplex': '#CC79A7',              # Reddish Purple
    'G-Quadruplex': '#D55E00',         # Vermillion
    'i-Motif': '#E69F00',              # Orange (lighter shade)
    'Z-DNA': '#0072B2',                # Blue
    'A-philic_DNA': '#F0E442',         # Yellow (lighter)
    'Hybrid': '#999999',               # Gray
    'Non-B_DNA_Clusters': '#666666'    # Dark Gray
}

def set_publication_style():
    """Set matplotlib style for publication-quality figures"""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'font.size': 10,
        'axes.titlesize': 12,
        'axes.labelsize': 11,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'figure.titlesize': 14,
        'axes.linewidth': 1.0,
        'xtick.major.width': 1.0,
        'ytick.major.width': 1.0,
        'grid.alpha': 0.3,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'savefig.dpi': 300,
        'savefig.format': 'png',
        'savefig.bbox': 'tight'
    })

# =============================================================================
# 1. GENOME LANDSCAPE TRACK
# =============================================================================

def plot_genome_landscape_track(motifs: List[Dict[str, Any]], 
                                sequence_length: int,
                                title: Optional[str] = None,
                                figsize: Tuple[int, int] = (16, 3)) -> plt.Figure:
    """
    Horizontal ruler with colored glyphs showing motif positions, lengths, and overlaps.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total sequence length
        title: Custom plot title
        figsize: Figure size (width, height)
    
    Returns:
        Matplotlib figure object
    """
    set_publication_style()
    
    fig, ax = plt.subplots(figsize=figsize)
    
    if not motifs:
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title or 'Genome Landscape Track')
        return fig
    
    # Group motifs by class for layering
    class_groups = defaultdict(list)
    for motif in motifs:
        class_groups[motif.get('Class', 'Unknown')].append(motif)
    
    # Draw each class on a different vertical level
    y_offset = 0
    class_positions = {}
    
    for i, (motif_class, class_motifs) in enumerate(sorted(class_groups.items())):
        color = MOTIF_CLASS_COLORS_CB.get(motif_class, COLORBLIND_PALETTE[i % len(COLORBLIND_PALETTE)])
        class_positions[motif_class] = y_offset
        
        for motif in class_motifs:
            start = motif.get('Start', 0)
            end = motif.get('End', 0)
            length = end - start
            
            # Draw as rectangle
            rect = patches.Rectangle((start, y_offset - 0.3), length, 0.6,
                                     linewidth=0.5, edgecolor='black',
                                     facecolor=color, alpha=0.7)
            ax.add_patch(rect)
        
        y_offset += 1
    
    # Add class labels
    for motif_class, y_pos in class_positions.items():
        ax.text(-sequence_length * 0.02, y_pos, motif_class, 
                ha='right', va='center', fontsize=9, style='italic')
    
    # Add baseline grid
    ax.axhline(y=-0.5, color='gray', linewidth=0.5, linestyle='--', alpha=0.3)
    ax.set_xlim(0, sequence_length)
    ax.set_ylim(-1, len(class_groups))
    ax.set_xlabel('Genomic Position (bp)', fontsize=11)
    ax.set_yticks([])
    ax.set_title(title or 'Genome Landscape Track', fontsize=14, fontweight='bold')
    
    # Add grid for major coordinates
    ax.grid(True, axis='x', alpha=0.3, linestyle=':')
    
    return fig

# =============================================================================
# 2. SLIDING WINDOW HEAT RIBBON
# =============================================================================

def plot_sliding_window_heat_ribbon(motifs: List[Dict[str, Any]],
                                   sequence_length: int,
                                   window_size: int = 1000,
                                   title: Optional[str] = None,
                                   figsize: Tuple[int, int] = (14, 4)) -> plt.Figure:
    """
    1D heatmap ribbon showing motif density with average score overlay.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total sequence length
        window_size: Window size for density calculation
        title: Custom plot title
        figsize: Figure size (width, height)
    
    Returns:
        Matplotlib figure object
    """
    set_publication_style()
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, 
                                    gridspec_kw={'height_ratios': [1, 3]})
    
    if not motifs:
        ax1.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax1.transAxes, fontsize=14)
        ax1.set_title(title or 'Sliding Window Heat Ribbon')
        ax2.axis('off')
        return fig
    
    # Calculate density and score for each window
    num_windows = max(1, sequence_length // window_size)
    densities = np.zeros(num_windows)
    avg_scores = np.zeros(num_windows)
    window_counts = np.zeros(num_windows)
    
    for motif in motifs:
        start = motif.get('Start', 0)
        end = motif.get('End', 0)
        score = motif.get('Score', motif.get('Normalized_Score', 0))
        
        # Find which windows this motif overlaps
        start_window = start // window_size
        end_window = min(end // window_size, num_windows - 1)
        
        for w in range(start_window, end_window + 1):
            densities[w] += 1
            avg_scores[w] += score
            window_counts[w] += 1
    
    # Calculate average scores
    avg_scores = np.divide(avg_scores, window_counts, 
                          out=np.zeros_like(avg_scores), 
                          where=window_counts != 0)
    
    # Plot density heatmap
    density_matrix = densities.reshape(1, -1)
    im = ax1.imshow(density_matrix, aspect='auto', cmap='YlOrRd',
                    extent=[0, sequence_length, 0, 1], interpolation='bilinear')
    ax1.set_yticks([])
    ax1.set_xlim(0, sequence_length)
    ax1.set_title('Motif Density Heatmap', fontsize=11)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax1, orientation='horizontal', pad=0.1, aspect=30)
    cbar.set_label('Motif Count', fontsize=9)
    
    # Plot average score line
    window_centers = np.arange(num_windows) * window_size + window_size / 2
    ax2.plot(window_centers, avg_scores, linewidth=2, color='#0072B2', label='Avg Score')
    ax2.fill_between(window_centers, 0, avg_scores, alpha=0.3, color='#0072B2')
    
    # Annotate peaks
    if len(avg_scores) > 0:
        peak_threshold = np.mean(avg_scores) + np.std(avg_scores)
        peaks = np.where(avg_scores > peak_threshold)[0]
        for peak in peaks[:3]:  # Annotate top 3 peaks
            ax2.annotate(f'{avg_scores[peak]:.2f}', 
                        xy=(window_centers[peak], avg_scores[peak]),
                        xytext=(0, 10), textcoords='offset points',
                        ha='center', fontsize=8,
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.5))
    
    ax2.set_xlabel('Genomic Position (bp)', fontsize=11)
    ax2.set_ylabel('Average Score', fontsize=11)
    ax2.set_xlim(0, sequence_length)
    ax2.grid(True, alpha=0.3, linestyle=':')
    ax2.legend(loc='upper right', fontsize=9)
    
    plt.suptitle(title or 'Sliding Window Analysis', fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    return fig

# =============================================================================
# 3. RIDGE PLOTS (JOYPLOTS) - LENGTH BY CLASS
# =============================================================================

def plot_ridge_plots_length_by_class(motifs: List[Dict[str, Any]],
                                    title: Optional[str] = None,
                                    figsize: Tuple[int, int] = (12, 8)) -> plt.Figure:
    """
    Stacked density ridges comparing length distributions across classes.
    
    Args:
        motifs: List of motif dictionaries
        title: Custom plot title
        figsize: Figure size (width, height)
    
    Returns:
        Matplotlib figure object
    """
    set_publication_style()
    
    fig, ax = plt.subplots(figsize=figsize)
    
    if not motifs:
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title or 'Ridge Plots: Length Distribution by Class')
        return fig
    
    # Group by class
    class_lengths = defaultdict(list)
    for motif in motifs:
        motif_class = motif.get('Class', 'Unknown')
        length = motif.get('Length', 0)
        if length > 0:
            class_lengths[motif_class].append(length)
    
    if not class_lengths:
        ax.text(0.5, 0.5, 'No length data available', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        return fig
    
    # Sort classes by median length
    sorted_classes = sorted(class_lengths.items(), 
                          key=lambda x: np.median(x[1]))
    
    # Plot ridge densities
    y_offset = 0
    overlap = 0.3  # Ridge overlap factor
    
    for i, (motif_class, lengths) in enumerate(sorted_classes):
        color = MOTIF_CLASS_COLORS_CB.get(motif_class, COLORBLIND_PALETTE[i % len(COLORBLIND_PALETTE)])
        
        # Calculate KDE
        from scipy import stats
        kde = stats.gaussian_kde(lengths)
        x_range = np.linspace(min(lengths), max(lengths), 200)
        density = kde(x_range)
        
        # Normalize density for plotting
        density_norm = density / density.max() * 0.8
        
        # Plot filled density
        ax.fill_between(x_range, y_offset, y_offset + density_norm,
                       color=color, alpha=0.6, linewidth=1.5, edgecolor='black')
        
        # Add median line
        median_val = np.median(lengths)
        median_density = kde(median_val) / kde(x_range).max() * 0.8
        ax.plot([median_val, median_val], [y_offset, y_offset + median_density],
               color='black', linewidth=1, linestyle='--', alpha=0.7)
        
        # Add class label
        ax.text(min(lengths) - 5, y_offset + 0.4, motif_class,
               ha='right', va='center', fontsize=10, style='italic')
        
        y_offset += 1 - overlap
    
    ax.set_xlabel('Motif Length (bp)', fontsize=12)
    ax.set_yticks([])
    ax.set_title(title or 'Ridge Plots: Length Distribution by Class', 
                fontsize=14, fontweight='bold')
    ax.spines['left'].set_visible(False)
    ax.grid(True, axis='x', alpha=0.3, linestyle=':')
    
    plt.tight_layout()
    return fig

# =============================================================================
# 4. SUNBURST/TREEMAP - HIERARCHICAL COMPOSITION
# =============================================================================

def plot_sunburst_treemap(motifs: List[Dict[str, Any]],
                         plot_type: str = 'sunburst',
                         title: Optional[str] = None,
                         figsize: Tuple[int, int] = (10, 10)) -> plt.Figure:
    """
    Hierarchical composition visualization (Class → Subclass → Method).
    
    Args:
        motifs: List of motif dictionaries
        plot_type: 'sunburst' or 'treemap'
        title: Custom plot title
        figsize: Figure size (width, height)
    
    Returns:
        Matplotlib figure object
    """
    set_publication_style()
    
    fig, ax = plt.subplots(figsize=figsize)
    
    if not motifs:
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title or f'{plot_type.title()} Chart')
        return fig
    
    # Build hierarchy: Class → Subclass
    hierarchy = defaultdict(lambda: defaultdict(int))
    for motif in motifs:
        motif_class = motif.get('Class', 'Unknown')
        subclass = motif.get('Subclass', 'Unknown')
        hierarchy[motif_class][subclass] += 1
    
    if plot_type == 'treemap':
        # Simple treemap using rectangles
        _plot_treemap(ax, hierarchy, MOTIF_CLASS_COLORS_CB)
    else:
        # Sunburst (concentric circles)
        _plot_sunburst(ax, hierarchy, MOTIF_CLASS_COLORS_CB)
    
    ax.set_title(title or f'{plot_type.title()}: Motif Hierarchy', 
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    return fig

def _plot_treemap(ax, hierarchy, color_map):
    """Helper to plot treemap"""
    # Calculate total
    total = sum(sum(subcounts.values()) for subcounts in hierarchy.values())
    
    # Squarify algorithm (simplified)
    x, y = 0, 0
    width, height = 1, 1
    
    class_sizes = [(cls, sum(subcounts.values())) 
                   for cls, subcounts in hierarchy.items()]
    class_sizes.sort(key=lambda x: x[1], reverse=True)
    
    for i, (motif_class, class_total) in enumerate(class_sizes):
        class_frac = class_total / total
        color = color_map.get(motif_class, COLORBLIND_PALETTE[i % len(COLORBLIND_PALETTE)])
        
        rect_width = class_frac * width
        rect = patches.Rectangle((x, y), rect_width, height,
                                linewidth=2, edgecolor='white',
                                facecolor=color, alpha=0.7)
        ax.add_patch(rect)
        
        # Add label
        if rect_width > 0.05:
            ax.text(x + rect_width/2, y + height/2, 
                   f"{motif_class}\n({class_total})",
                   ha='center', va='center', fontsize=9,
                   color='white', weight='bold',
                   bbox=dict(boxstyle='round', facecolor='black', alpha=0.3))
        
        x += rect_width
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

def _plot_sunburst(ax, hierarchy, color_map):
    """Helper to plot sunburst"""
    total = sum(sum(subcounts.values()) for subcounts in hierarchy.values())
    
    # Inner circle for classes
    theta_start = 0
    for i, (motif_class, subcounts) in enumerate(hierarchy.items()):
        class_total = sum(subcounts.values())
        theta_extent = (class_total / total) * 360
        color = color_map.get(motif_class, COLORBLIND_PALETTE[i % len(COLORBLIND_PALETTE)])
        
        wedge = patches.Wedge((0, 0), 0.5, theta_start, theta_start + theta_extent,
                             width=0.2, facecolor=color, edgecolor='white',
                             linewidth=2, alpha=0.8)
        ax.add_patch(wedge)
        
        # Add class label
        label_theta = np.radians(theta_start + theta_extent / 2)
        label_r = 0.35
        ax.text(label_r * np.cos(label_theta), label_r * np.sin(label_theta),
               motif_class.replace('_', '\n'), ha='center', va='center',
               fontsize=8, weight='bold', color='white')
        
        # Outer circle for subclasses
        sub_theta_start = theta_start
        for subclass, count in subcounts.items():
            sub_theta_extent = (count / class_total) * theta_extent
            lighter_color = plt.cm.colors.to_rgba(color)
            lighter_color = tuple(min(1, c + 0.2) for c in lighter_color[:3]) + (0.6,)
            
            sub_wedge = patches.Wedge((0, 0), 0.8, sub_theta_start, 
                                     sub_theta_start + sub_theta_extent,
                                     width=0.25, facecolor=lighter_color,
                                     edgecolor='white', linewidth=1)
            ax.add_patch(sub_wedge)
            
            sub_theta_start += sub_theta_extent
        
        theta_start += theta_extent
    
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_aspect('equal')
    ax.axis('off')

# =============================================================================
# 5. HEXBIN (START VS SCORE) WITH MARGINALS
# =============================================================================

def plot_hexbin_start_vs_score(motifs: List[Dict[str, Any]],
                               title: Optional[str] = None,
                               figsize: Tuple[int, int] = (12, 10)) -> plt.Figure:
    """
    2D density hexbin plot with marginal histograms showing positional score structure.
    
    Args:
        motifs: List of motif dictionaries
        title: Custom plot title
        figsize: Figure size (width, height)
    
    Returns:
        Matplotlib figure object
    """
    set_publication_style()
    
    # Create figure with gridspec for marginals
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(3, 3, hspace=0.02, wspace=0.02,
                          height_ratios=[1, 4, 0.3], width_ratios=[0.3, 4, 1])
    
    ax_main = fig.add_subplot(gs[1, 1])
    ax_top = fig.add_subplot(gs[0, 1], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1, 2], sharey=ax_main)
    
    if not motifs:
        ax_main.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                    transform=ax_main.transAxes, fontsize=14)
        ax_main.set_title(title or 'Hexbin: Start vs Score')
        return fig
    
    # Extract data
    starts = [m.get('Start', 0) for m in motifs]
    scores = [m.get('Score', m.get('Normalized_Score', 0)) for m in motifs]
    
    # Main hexbin plot
    hexbin = ax_main.hexbin(starts, scores, gridsize=30, cmap='YlOrRd', 
                           mincnt=1, edgecolors='none', alpha=0.8)
    ax_main.set_xlabel('Start Position (bp)', fontsize=11)
    ax_main.set_ylabel('Normalized Score', fontsize=11)
    ax_main.grid(True, alpha=0.3, linestyle=':')
    
    # Colorbar
    cbar = plt.colorbar(hexbin, cax=fig.add_subplot(gs[2, 1]), orientation='horizontal')
    cbar.set_label('Motif Count', fontsize=9)
    
    # Top marginal histogram
    ax_top.hist(starts, bins=50, color='#0072B2', alpha=0.7, edgecolor='black')
    ax_top.set_ylabel('Count', fontsize=9)
    ax_top.tick_params(labelbottom=False)
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    
    # Right marginal histogram
    ax_right.hist(scores, bins=50, orientation='horizontal', 
                 color='#E69F00', alpha=0.7, edgecolor='black')
    ax_right.set_xlabel('Count', fontsize=9)
    ax_right.tick_params(labelleft=False)
    ax_right.spines['top'].set_visible(False)
    ax_right.spines['right'].set_visible(False)
    
    plt.suptitle(title or 'Hexbin: Start Position vs Normalized Score', 
                fontsize=14, fontweight='bold')
    
    return fig

# =============================================================================
# EXPORT UTILITIES
# =============================================================================

def export_plot(fig: plt.Figure, filename: str, 
                formats: List[str] = ['png', 'svg'],
                dpi: int = 300) -> Dict[str, str]:
    """
    Export plot in multiple formats.
    
    Args:
        fig: Matplotlib figure
        filename: Base filename (without extension)
        formats: List of formats ('png', 'svg', 'pdf')
        dpi: DPI for raster formats
    
    Returns:
        Dictionary of {format: filepath}
    """
    saved_files = {}
    
    for fmt in formats:
        filepath = f"{filename}.{fmt}"
        if fmt == 'svg':
            fig.savefig(filepath, format='svg', bbox_inches='tight')
        elif fmt == 'pdf':
            fig.savefig(filepath, format='pdf', bbox_inches='tight')
        else:  # png or other raster
            fig.savefig(filepath, format=fmt, dpi=dpi, bbox_inches='tight')
        
        saved_files[fmt] = filepath
        print(f"✓ Saved {fmt.upper()}: {filepath}")
    
    return saved_files

# =============================================================================
# 6. UPSET PLOT - INTERSECTION ANALYSIS
# =============================================================================

def plot_upset_intersection(motifs: List[Dict[str, Any]],
                           title: Optional[str] = None,
                           figsize: Tuple[int, int] = (14, 8)) -> plt.Figure:
    """
    UpSet plot showing motif/hybrid overlaps (better than Venn diagrams).
    
    Args:
        motifs: List of motif dictionaries
        title: Custom plot title
        figsize: Figure size (width, height)
    
    Returns:
        Matplotlib figure object
    """
    set_publication_style()
    
    fig, (ax_matrix, ax_bars) = plt.subplots(2, 1, figsize=figsize,
                                              gridspec_kw={'height_ratios': [1, 2]})
    
    if not motifs:
        ax_matrix.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                      transform=ax_matrix.transAxes, fontsize=14)
        ax_matrix.set_title(title or 'UpSet Intersection Plot')
        ax_bars.axis('off')
        return fig
    
    # Find overlapping motifs
    from itertools import combinations
    
    # Get all unique classes
    all_classes = sorted(set(m.get('Class', 'Unknown') for m in motifs))
    
    # Find intersections
    intersections = defaultdict(list)
    
    for motif in motifs:
        # Check which classes this motif overlaps with
        motif_classes = set([motif.get('Class', 'Unknown')])
        
        # Check for overlaps (hybrids)
        if 'Component_Classes' in motif:
            motif_classes.update(motif['Component_Classes'])
        
        # Create frozenset key for intersection
        class_key = frozenset(motif_classes)
        intersections[class_key].append(motif)
    
    # Sort intersections by size
    sorted_intersections = sorted(intersections.items(), 
                                 key=lambda x: len(x[1]), reverse=True)[:15]  # Top 15
    
    # Plot bars
    counts = [len(motif_list) for _, motif_list in sorted_intersections]
    x_pos = np.arange(len(counts))
    bars = ax_bars.bar(x_pos, counts, color='#0072B2', alpha=0.7, edgecolor='black')
    ax_bars.set_ylabel('Intersection Size', fontsize=11)
    ax_bars.set_xticks([])
    ax_bars.grid(True, axis='y', alpha=0.3)
    ax_bars.spines['top'].set_visible(False)
    ax_bars.spines['right'].set_visible(False)
    
    # Annotate top bars
    for i, (bar, count) in enumerate(zip(bars, counts)):
        if i < 5:  # Annotate top 5
            ax_bars.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                        str(count), ha='center', va='bottom', fontsize=9)
    
    # Plot matrix
    matrix_data = np.zeros((len(all_classes), len(sorted_intersections)))
    
    for col, (class_set, _) in enumerate(sorted_intersections):
        for row, cls in enumerate(all_classes):
            if cls in class_set:
                matrix_data[row, col] = 1
    
    # Draw matrix as circles
    for row in range(len(all_classes)):
        for col in range(len(sorted_intersections)):
            if matrix_data[row, col] == 1:
                circle = patches.Circle((col, row), 0.3, 
                                       facecolor='black', edgecolor='none')
                ax_matrix.add_patch(circle)
                # Draw connecting lines for multi-class intersections
                if np.sum(matrix_data[:, col]) > 1:
                    ax_matrix.plot([col, col], 
                                 [np.where(matrix_data[:, col] == 1)[0].min(),
                                  np.where(matrix_data[:, col] == 1)[0].max()],
                                 'k-', linewidth=2, alpha=0.5)
    
    ax_matrix.set_yticks(range(len(all_classes)))
    ax_matrix.set_yticklabels(all_classes, fontsize=9)
    ax_matrix.set_xlim(-0.5, len(sorted_intersections) - 0.5)
    ax_matrix.set_ylim(-0.5, len(all_classes) - 0.5)
    ax_matrix.set_xticks([])
    ax_matrix.invert_yaxis()
    ax_matrix.spines['top'].set_visible(False)
    ax_matrix.spines['right'].set_visible(False)
    ax_matrix.spines['bottom'].set_visible(False)
    ax_matrix.set_ylabel('Motif Classes', fontsize=11)
    
    plt.suptitle(title or 'UpSet Plot: Motif Class Intersections', 
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    return fig

# =============================================================================
# 7. SCORE STACKED VIOLIN + BEESWARM
# =============================================================================

def plot_score_violin_beeswarm(motifs: List[Dict[str, Any]],
                              title: Optional[str] = None,
                              figsize: Tuple[int, int] = (14, 8)) -> plt.Figure:
    """
    Distribution shapes with individual points for outliers (per subclass).
    
    Args:
        motifs: List of motif dictionaries
        title: Custom plot title
        figsize: Figure size (width, height)
    
    Returns:
        Matplotlib figure object
    """
    set_publication_style()
    
    fig, ax = plt.subplots(figsize=figsize)
    
    if not motifs:
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title or 'Score Distribution: Violin + Beeswarm')
        return fig
    
    # Prepare data by class
    class_scores = defaultdict(list)
    for motif in motifs:
        motif_class = motif.get('Class', 'Unknown')
        score = motif.get('Score', motif.get('Normalized_Score', 0))
        if score > 0:
            class_scores[motif_class].append(score)
    
    if not class_scores:
        ax.text(0.5, 0.5, 'No score data available', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        return fig
    
    # Sort by median score
    sorted_classes = sorted(class_scores.items(), 
                          key=lambda x: np.median(x[1]), reverse=True)
    
    # Prepare data for violin plot
    data_for_violin = [scores for _, scores in sorted_classes]
    labels = [cls for cls, _ in sorted_classes]
    colors = [MOTIF_CLASS_COLORS_CB.get(cls, COLORBLIND_PALETTE[i % len(COLORBLIND_PALETTE)]) 
             for i, cls in enumerate(labels)]
    
    # Create violin plot
    parts = ax.violinplot(data_for_violin, positions=range(len(labels)),
                         widths=0.7, showmeans=True, showmedians=True)
    
    # Color violins
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_alpha(0.6)
        pc.set_edgecolor('black')
        pc.set_linewidth(1)
    
    # Style median and mean lines
    for key in ['cmeans', 'cmedians']:
        if key in parts:
            parts[key].set_edgecolor('black')
            parts[key].set_linewidth(2)
    
    # Add beeswarm points (jittered)
    for i, (cls, scores) in enumerate(sorted_classes):
        # Add jitter
        y_data = scores
        x_data = np.random.normal(i, 0.04, size=len(scores))
        
        # Identify top 1% as outliers
        threshold = np.percentile(scores, 99) if len(scores) > 10 else max(scores)
        
        for x, y in zip(x_data, y_data):
            if y >= threshold:
                ax.scatter(x, y, s=50, color='gold', edgecolor='black', 
                          linewidth=1.5, zorder=10, marker='*')
            else:
                ax.scatter(x, y, s=20, color=colors[i], edgecolor='black',
                          linewidth=0.5, alpha=0.6, zorder=5)
    
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Normalized Score', fontsize=11)
    ax.set_title(title or 'Score Distribution: Violin + Beeswarm (★ = top 1%)', 
                fontsize=14, fontweight='bold')
    ax.grid(True, axis='y', alpha=0.3, linestyle=':')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    return fig

# =============================================================================
# 8. CLUSTER HOTSPOT MAP
# =============================================================================

def plot_cluster_hotspot_map(motifs: List[Dict[str, Any]],
                            sequence_length: int,
                            title: Optional[str] = None,
                            figsize: Tuple[int, int] = (16, 6)) -> plt.Figure:
    """
    Bar chart showing cluster counts per region with representative sequences.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total sequence length
        title: Custom plot title
        figsize: Figure size (width, height)
    
    Returns:
        Matplotlib figure object
    """
    set_publication_style()
    
    fig, ax = plt.subplots(figsize=figsize)
    
    if not motifs:
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title or 'Cluster Hotspot Map')
        return fig
    
    # Filter for cluster motifs
    clusters = [m for m in motifs if m.get('Class') == 'Non-B_DNA_Clusters']
    
    if not clusters:
        ax.text(0.5, 0.5, 'No clusters found', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title or 'Cluster Hotspot Map')
        return fig
    
    # Divide sequence into regions
    region_size = max(1000, sequence_length // 10)
    n_regions = (sequence_length + region_size - 1) // region_size
    
    region_counts = defaultdict(list)
    for cluster in clusters:
        region_idx = cluster.get('Start', 0) // region_size
        region_counts[region_idx].append(cluster)
    
    # Plot bars
    regions = range(n_regions)
    counts = [len(region_counts.get(i, [])) for i in regions]
    region_labels = [f"{i*region_size//1000}-{(i+1)*region_size//1000}kb" 
                    for i in regions]
    
    bars = ax.bar(regions, counts, color='#0072B2', alpha=0.7, edgecolor='black')
    
    # Annotate top regions with cluster details
    top_regions = sorted(region_counts.items(), key=lambda x: len(x[1]), reverse=True)[:3]
    
    for region_idx, region_clusters in top_regions:
        if region_clusters:
            # Get top cluster
            top_cluster = max(region_clusters, key=lambda x: x.get('Score', 0))
            
            # Create annotation
            annotation_text = (f"Region {region_idx}:\n"
                             f"Motifs: {top_cluster.get('Motif_Count', 'N/A')}\n"
                             f"Score: {top_cluster.get('Score', 0):.2f}")
            
            ax.annotate(annotation_text, xy=(region_idx, counts[region_idx]),
                       xytext=(0, 20), textcoords='offset points',
                       ha='center', fontsize=8,
                       bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7),
                       arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    
    ax.set_xlabel('Genomic Region', fontsize=11)
    ax.set_ylabel('Cluster Count', fontsize=11)
    ax.set_xticks(regions)
    ax.set_xticklabels(region_labels, rotation=45, ha='right', fontsize=8)
    ax.set_title(title or 'Cluster Hotspot Map', fontsize=14, fontweight='bold')
    ax.grid(True, axis='y', alpha=0.3, linestyle=':')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    return fig

# =============================================================================
# TESTING
# =============================================================================

if __name__ == "__main__":
    print("Advanced Visualization Suite - NBDScanner")
    print("=" * 50)
    
    # Test with example data
    test_motifs = [
        {'Class': 'G-Quadruplex', 'Subclass': 'Canonical', 'Start': 10, 'End': 30, 
         'Length': 20, 'Score': 0.85, 'Normalized_Score': 0.75},
        {'Class': 'G-Quadruplex', 'Subclass': 'Relaxed', 'Start': 50, 'End': 68, 
         'Length': 18, 'Score': 0.72, 'Normalized_Score': 0.65},
        {'Class': 'Curved_DNA', 'Subclass': 'A-tract', 'Start': 90, 'End': 105, 
         'Length': 15, 'Score': 0.65, 'Normalized_Score': 0.55},
        {'Class': 'Z-DNA', 'Subclass': 'CG-alternating', 'Start': 130, 'End': 146, 
         'Length': 16, 'Score': 0.90, 'Normalized_Score': 0.82},
        {'Class': 'i-Motif', 'Subclass': 'Canonical', 'Start': 170, 'End': 190, 
         'Length': 20, 'Score': 0.78, 'Normalized_Score': 0.70}
    ]
    
    print(f"\nTesting with {len(test_motifs)} example motifs")
    
    # Test genome landscape
    try:
        fig1 = plot_genome_landscape_track(test_motifs, 200)
        plt.close(fig1)
        print("✓ Genome landscape track: PASS")
    except Exception as e:
        print(f"✗ Genome landscape track: FAIL - {e}")
    
    # Test heat ribbon
    try:
        fig2 = plot_sliding_window_heat_ribbon(test_motifs, 200)
        plt.close(fig2)
        print("✓ Sliding window heat ribbon: PASS")
    except Exception as e:
        print(f"✗ Sliding window heat ribbon: FAIL - {e}")
    
    # Test ridge plots
    try:
        fig3 = plot_ridge_plots_length_by_class(test_motifs)
        plt.close(fig3)
        print("✓ Ridge plots: PASS")
    except Exception as e:
        print(f"✗ Ridge plots: FAIL - {e}")
    
    print("\n✓ Testing completed")
