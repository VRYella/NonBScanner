"""
Optimized Visualization Module for NBDFinder
===========================================

Fast visualization generation with caching and optimized plotting.

Key optimizations:
- Cached visualization data preparation
- Optimized plot generation with minimal overhead
- Fast rendering options for large datasets
- Memory-efficient processing

Author: NBDFinder Optimization Team
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from matplotlib_venn import venn2
import random
import networkx as nx
from sklearn.manifold import TSNE
from functools import lru_cache
from typing import List, Dict, Any, Optional
import time
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Set matplotlib backend for better performance
plt.style.use('fast')
plt.rcParams['figure.max_open_warning'] = 0


@lru_cache(maxsize=32)
def _prepare_visualization_data(data_hash: str, motifs_tuple) -> pd.DataFrame:
    """Cached data preparation for visualizations."""
    # Convert back from tuple
    motifs_data = []
    for motif_tuple in motifs_tuple:
        motif_dict = dict(motif_tuple)
        motifs_data.append(motif_dict)
    
    return pd.DataFrame(motifs_data)


def motifs_to_dataframe(motifs: List[Dict[str, Any]]) -> pd.DataFrame:
    """Convert motifs to DataFrame with caching for repeated use."""
    if not motifs:
        return generate_pseudodata_optimized(50)  # Minimal demo data
    
    # Create hashable representation for caching
    data_hash = str(hash(str(sorted(motifs, key=lambda x: x.get('Start', 0)))))
    motifs_tuple = tuple(
        tuple(sorted(m.items())) for m in motifs
    )
    
    return _prepare_visualization_data(data_hash, motifs_tuple)


def generate_pseudodata_optimized(n_motifs=100):
    """Generate optimized pseudo data for fast demo visualizations."""
    np.random.seed(42)
    
    # Simplified class structure for faster generation
    motif_classes_subclasses = {
        'Curved_DNA': ['A-Tract', 'Intrinsic'],
        'Slipped_DNA': ['STR', 'Direct_Repeat'],
        'Cruciform_DNA': ['Inverted_Repeat', 'Hairpin'],
        'R-Loop': ['RLFS_m1', 'RLFS_m2'],
        'Triplex_DNA': ['Triplex', 'Sticky_DNA'],
        'G-Quadruplex': ['Canonical_G4', 'Relaxed_G4', 'Bulged_G4'],
        'i-Motif': ['Canonical_iMotif', 'Relaxed_iMotif'],
        'Z-DNA': ['Z-DNA', 'eGZ'],
        'Hybrid': ['Dynamic_Overlap'],
        'Cluster': ['Hotspot_Region']
    }

    data = []
    for i in range(n_motifs):
        cl = random.choice(list(motif_classes_subclasses.keys()))
        sub = random.choice(motif_classes_subclasses[cl])
        seq_len = np.random.randint(100, 1000)
        start = np.random.randint(1, seq_len - 20)
        end = start + np.random.randint(10, 50)
        score = np.abs(np.random.normal(loc=50.0, scale=20.0))
        
        data.append({
            'Class': cl, 'Subclass': sub, 'Start': start, 'End': end,
            'Actual_Score': score, 'Normalized_Score': min(score, 100.0),
            'Length': end-start, 'GC_Content': np.random.uniform(30, 80),
            'Sequence': 'ATGC'*3, 'Sequence_Name': f'Seq{np.random.randint(1,5)}',
            'Motif_ID': f'{cl}_{sub}_{i}', 'Scoring_Method': 'Optimized'
        })
    
    return pd.DataFrame(data)


def create_optimized_visualizations(motifs: List[Dict[str, Any]] = None, 
                                  save_plots: bool = False, 
                                  output_dir: str = './plots/',
                                  fast_mode: bool = True) -> Dict[str, Any]:
    """
    Create optimized visualizations with performance monitoring.
    
    Args:
        motifs: List of motif dictionaries
        save_plots: Whether to save plots to files
        output_dir: Directory for saved plots
        fast_mode: Use fast rendering options
        
    Returns:
        Dictionary with performance metrics and plot information
    """
    import os
    if save_plots and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Prepare data
    start_time = time.time()
    if motifs:
        df = motifs_to_dataframe(motifs)
    else:
        df = generate_pseudodata_optimized()
    
    data_prep_time = time.time() - start_time
    
    print(f"Creating optimized visualizations for {len(df)} motifs...")
    
    plot_times = {}
    created_plots = []
    
    try:
        # 1. Essential plots first (fast)
        plots_to_create = [
            ("motif_counts", plot_motif_counts_optimized, "Motif class counts"),
            ("score_distribution", plot_score_distribution_optimized, "Score distributions"),
            ("genomic_positions", plot_genomic_positions_optimized, "Genomic positions"),
            ("class_summary", plot_class_summary_optimized, "Class summary")
        ]
        
        if not fast_mode:
            # Add more complex plots in full mode
            plots_to_create.extend([
                ("stacked_distribution", plot_stacked_distribution_optimized, "Stacked distribution"),
                ("correlation_heatmap", plot_correlation_heatmap_optimized, "Feature correlations"),
                ("interactive_overview", create_interactive_overview, "Interactive overview"),
                ("network_analysis", plot_network_optimized, "Network analysis")
            ])
        
        for plot_name, plot_func, description in plots_to_create:
            try:
                plot_start = time.time()
                print(f"  Creating {description}...")
                
                plot_func(df, save_plots, output_dir)
                
                plot_time = time.time() - plot_start
                plot_times[plot_name] = plot_time
                created_plots.append(plot_name)
                
                if plot_time > 1.0:  # Log slow plots
                    print(f"    ⚠️  {description} took {plot_time:.2f}s")
                    
            except Exception as e:
                print(f"    ❌ Failed to create {description}: {e}")
        
        total_time = time.time() - start_time
        
        print(f"\n✓ Created {len(created_plots)} visualizations in {total_time:.2f}s")
        print(f"  Data preparation: {data_prep_time:.2f}s")
        print(f"  Plot generation: {total_time - data_prep_time:.2f}s")
        
        return {
            'success': True,
            'total_time': total_time,
            'data_prep_time': data_prep_time,
            'plot_times': plot_times,
            'created_plots': created_plots,
            'motif_count': len(df)
        }
        
    except Exception as e:
        print(f"❌ Visualization generation failed: {e}")
        import traceback
        traceback.print_exc()
        return {'success': False, 'error': str(e)}


def plot_motif_counts_optimized(df, save_plots=False, output_dir='./plots/'):
    """Optimized motif count visualization."""
    plt.figure(figsize=(10, 6))
    
    # Use value_counts for efficiency
    counts = df['Class'].value_counts()
    
    # Create bar plot
    ax = counts.plot(kind='bar', color='steelblue', alpha=0.8)
    plt.title("Motif Count per Class", fontsize=14, fontweight='bold')
    plt.xlabel("Motif Class", fontsize=12)
    plt.ylabel("Count", fontsize=12)
    plt.xticks(rotation=45, ha='right')
    
    # Add count labels on bars
    for i, v in enumerate(counts.values):
        ax.text(i, v + max(counts) * 0.01, str(v), ha='center', va='bottom')
    
    plt.tight_layout()
    
    if save_plots:
        plt.savefig(f"{output_dir}/motif_counts.png", dpi=150, bbox_inches='tight')
    
    plt.show()


def plot_score_distribution_optimized(df, save_plots=False, output_dir='./plots/'):
    """Optimized score distribution visualization."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Histogram of normalized scores
    ax1.hist(df['Normalized_Score'].dropna(), bins=20, color='lightblue', 
             alpha=0.7, edgecolor='black')
    ax1.set_title("Normalized Score Distribution")
    ax1.set_xlabel("Normalized Score")
    ax1.set_ylabel("Frequency")
    
    # Box plot by class (top 6 classes for clarity)
    top_classes = df['Class'].value_counts().head(6).index
    df_top = df[df['Class'].isin(top_classes)]
    
    df_top.boxplot(column='Normalized_Score', by='Class', ax=ax2)
    ax2.set_title("Score Distribution by Class (Top 6)")
    ax2.set_xlabel("Motif Class")
    ax2.set_ylabel("Normalized Score")
    
    plt.suptitle("")  # Remove default title
    plt.tight_layout()
    
    if save_plots:
        plt.savefig(f"{output_dir}/score_distribution.png", dpi=150, bbox_inches='tight')
    
    plt.show()


def plot_genomic_positions_optimized(df, save_plots=False, output_dir='./plots/'):
    """Optimized genomic position visualization."""
    plt.figure(figsize=(12, 6))
    
    # Sample data if too large for performance
    if len(df) > 1000:
        df_sample = df.sample(n=1000, random_state=42)
    else:
        df_sample = df
    
    # Scatter plot with class colors
    classes = df_sample['Class'].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(classes)))
    
    for i, cls in enumerate(classes):
        cls_data = df_sample[df_sample['Class'] == cls]
        plt.scatter(cls_data['Start'], cls_data['Length'], 
                   c=[colors[i]], label=cls, alpha=0.6, s=30)
    
    plt.xlabel("Genomic Position (Start)")
    plt.ylabel("Motif Length")
    plt.title("Motif Positions and Lengths")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    if save_plots:
        plt.savefig(f"{output_dir}/genomic_positions.png", dpi=150, bbox_inches='tight')
    
    plt.show()


def plot_class_summary_optimized(df, save_plots=False, output_dir='./plots/'):
    """Optimized class summary visualization."""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. Class counts (pie chart)
    class_counts = df['Class'].value_counts()
    ax1.pie(class_counts.values, labels=class_counts.index, autopct='%1.1f%%', 
            startangle=90)
    ax1.set_title("Motif Class Distribution")
    
    # 2. Average scores by class
    avg_scores = df.groupby('Class')['Normalized_Score'].mean().sort_values(ascending=True)
    avg_scores.plot(kind='barh', ax=ax2, color='lightcoral')
    ax2.set_title("Average Score by Class")
    ax2.set_xlabel("Average Normalized Score")
    
    # 3. Length distribution
    ax3.hist(df['Length'].dropna(), bins=15, color='lightgreen', alpha=0.7)
    ax3.set_title("Motif Length Distribution")
    ax3.set_xlabel("Length (bp)")
    ax3.set_ylabel("Frequency")
    
    # 4. GC content vs Score
    ax4.scatter(df['GC_Content'], df['Normalized_Score'], alpha=0.5, color='purple')
    ax4.set_title("GC Content vs Score")
    ax4.set_xlabel("GC Content (%)")
    ax4.set_ylabel("Normalized Score")
    
    plt.tight_layout()
    
    if save_plots:
        plt.savefig(f"{output_dir}/class_summary.png", dpi=150, bbox_inches='tight')
    
    plt.show()


def plot_stacked_distribution_optimized(df, save_plots=False, output_dir='./plots/'):
    """Optimized stacked distribution plot."""
    # Create pivot table
    pivot_data = df.groupby(['Class', 'Subclass']).size().unstack(fill_value=0)
    
    # Limit to top classes for readability
    if len(pivot_data) > 8:
        top_classes = df['Class'].value_counts().head(8).index
        pivot_data = pivot_data.loc[top_classes]
    
    plt.figure(figsize=(12, 6))
    pivot_data.plot(kind='bar', stacked=True, colormap='tab20', 
                   figsize=(12, 6), width=0.8)
    
    plt.title("Stacked Bar: Motif Subclass Distribution")
    plt.xlabel("Motif Class")
    plt.ylabel("Count")
    plt.xticks(rotation=45, ha='right')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    if save_plots:
        plt.savefig(f"{output_dir}/stacked_distribution.png", dpi=150, bbox_inches='tight')
    
    plt.show()


def plot_correlation_heatmap_optimized(df, save_plots=False, output_dir='./plots/'):
    """Optimized correlation heatmap."""
    # Select numeric columns
    numeric_cols = ['Start', 'End', 'Length', 'Normalized_Score', 'Actual_Score', 'GC_Content']
    available_cols = [col for col in numeric_cols if col in df.columns]
    
    if len(available_cols) < 2:
        print("  ⚠️  Insufficient numeric columns for correlation heatmap")
        return
    
    corr_data = df[available_cols].corr()
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(corr_data, annot=True, cmap='coolwarm', center=0, 
                square=True, cbar_kws={'shrink': 0.8})
    plt.title("Feature Correlation Heatmap")
    plt.tight_layout()
    
    if save_plots:
        plt.savefig(f"{output_dir}/correlation_heatmap.png", dpi=150, bbox_inches='tight')
    
    plt.show()


def create_interactive_overview(df, save_plots=False, output_dir='./plots/'):
    """Create interactive overview using Plotly."""
    try:
        # Sample for performance if dataset is large
        if len(df) > 2000:
            df_sample = df.sample(n=2000, random_state=42)
        else:
            df_sample = df
        
        fig = px.scatter(df_sample, 
                        x='Start', y='Normalized_Score',
                        color='Class', size='Length',
                        hover_data=['Subclass', 'GC_Content'],
                        title="Interactive Motif Overview")
        
        fig.update_layout(height=600)
        
        if save_plots:
            fig.write_html(f"{output_dir}/interactive_overview.html")
        
        fig.show()
        
    except Exception as e:
        print(f"  ⚠️  Interactive plot failed: {e}")


def plot_network_optimized(df, save_plots=False, output_dir='./plots/'):
    """Optimized network visualization."""
    try:
        # Create simplified network based on sequence co-occurrence
        seq_motifs = df.groupby('Sequence_Name')['Class'].apply(list).to_dict()
        
        # Build edge list
        edges = []
        for seq, classes in seq_motifs.items():
            unique_classes = list(set(classes))
            for i in range(len(unique_classes)):
                for j in range(i+1, len(unique_classes)):
                    edges.append((unique_classes[i], unique_classes[j]))
        
        if not edges:
            print("  ⚠️  No edges found for network plot")
            return
        
        # Count edge weights
        edge_counts = {}
        for edge in edges:
            edge_counts[edge] = edge_counts.get(edge, 0) + 1
        
        # Create network
        G = nx.Graph()
        for (n1, n2), weight in edge_counts.items():
            G.add_edge(n1, n2, weight=weight)
        
        # Plot
        plt.figure(figsize=(10, 8))
        pos = nx.spring_layout(G, k=1, iterations=50)
        
        # Draw nodes
        nx.draw_networkx_nodes(G, pos, node_color='lightblue', 
                              node_size=500, alpha=0.8)
        
        # Draw edges with weights
        nx.draw_networkx_edges(G, pos, alpha=0.5, width=1)
        
        # Draw labels
        nx.draw_networkx_labels(G, pos, font_size=8)
        
        plt.title("Motif Class Co-occurrence Network")
        plt.axis('off')
        plt.tight_layout()
        
        if save_plots:
            plt.savefig(f"{output_dir}/network_analysis.png", dpi=150, bbox_inches='tight')
        
        plt.show()
        
    except Exception as e:
        print(f"  ⚠️  Network plot failed: {e}")


def benchmark_visualizations(motifs: List[Dict[str, Any]] = None, iterations: int = 3):
    """Benchmark visualization performance."""
    print("Benchmarking visualization performance...")
    
    times = []
    for i in range(iterations):
        print(f"  Run {i+1}/{iterations}")
        result = create_optimized_visualizations(motifs, save_plots=False, fast_mode=True)
        if result['success']:
            times.append(result['total_time'])
    
    if times:
        avg_time = sum(times) / len(times)
        print(f"\nBenchmark Results:")
        print(f"  Average time: {avg_time:.2f}s")
        print(f"  Min time: {min(times):.2f}s") 
        print(f"  Max time: {max(times):.2f}s")
        print(f"  Performance: {'✓ Fast' if avg_time < 2 else '⚠️ Moderate' if avg_time < 5 else '❌ Slow'}")
    else:
        print("  ❌ Benchmark failed")


# Legacy compatibility functions
def create_all_visualizations(df=None, save_plots=False, output_dir='./plots/'):
    """Legacy compatibility wrapper."""
    if df is not None:
        motifs = df.to_dict('records')
    else:
        motifs = None
    
    return create_optimized_visualizations(motifs, save_plots, output_dir, fast_mode=False)