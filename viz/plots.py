"""
Advanced Non-B DNA Motif Visualization Suite
---------------------------------------------
Features: counts, distributions, locations, overlaps, interactions, density, networks, dimensionality reduction, interactive plots.

Author: Copilot Space, 2025
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


def generate_pseudodata(n_motifs=300):
    """Generate pseudo data for demonstration with all 22 subclasses"""
    np.random.seed(42)
    
    # Updated to match the exact 10 classes and 22 subclasses
    motif_classes_subclasses = {
        'Curved_DNA': ['Global_Array', 'Local_Tract'],
        'Slipped_DNA': ['Direct_Repeat', 'STR'],
        'Cruciform_DNA': ['Inverted_Repeat', 'Hairpin'],  # Added 2nd subclass
        'R-Loop': ['RLFS_m1', 'RLFS_m2'],
        'Triplex_DNA': ['Triplex', 'Sticky_DNA'],
        'G-Quadruplex': ['Canonical_G4', 'Relaxed_G4', 'Bulged_G4', 'Bipartite_G4', 
                        'Multimeric_G4', 'Imperfect_G4', 'G-Triplex_intermediate'],
        'i-Motif': ['Canonical_iMotif', 'Relaxed_iMotif', 'AC-motif'],
        'Z-DNA': ['Z-DNA', 'eGZ'],
        'Hybrid': ['Dynamic_Overlap'],
        'Cluster': ['Hotspot_Region']
    }

    data = []
    for _ in range(n_motifs):
        cl = random.choice(list(motif_classes_subclasses.keys()))
        sub = random.choice(motif_classes_subclasses[cl])
        seq_len = np.random.randint(50, 2000)
        start = np.random.randint(1, seq_len - 20)
        end = start + np.random.randint(10, 50)
        score = np.abs(np.random.normal(loc=3.0, scale=1.2 if 'G-Quadruplex' in cl else 1.0))
        data.append({
            'Class': cl, 'Subclass': sub, 'Start': start, 'End': end,
            'Actual_Score': score, 'Normalized_Score': min(score/5.0, 1.0),
            'Length': end-start, 'GC_Content': np.random.uniform(30, 80),
            'Sequence': 'A'*10, 'Sequence_Name': f'Seq{np.random.randint(1,10)}',
            'Motif_ID': f'{cl}_{sub}_{_}', 'Scoring_Method': 'Simulated'
        })
    return pd.DataFrame(data)


def plot_motif_counts(df):
    """Plot motif count per class"""
    plt.figure(figsize=(7,4))
    sns.countplot(data=df, x='Class', order=df['Class'].value_counts().index)
    plt.title("Motif Count per Class")
    plt.tight_layout()
    plt.show()


def plot_stacked_distribution(df):
    """Stacked bar chart of subclass distribution"""
    pivot = df.groupby(['Class','Subclass']).size().reset_index(name='Count')
    pivot_pivot = pivot.pivot(index='Class',columns='Subclass',values='Count').fillna(0)
    pivot_pivot.plot(kind='bar', stacked=True, figsize=(12,6), colormap='tab20')
    plt.title("Stacked Bar: Motif Subclass Distribution")
    plt.ylabel("Count")
    plt.xlabel("Motif Class")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()


def plot_pie_chart(df):
    """Pie/Donut chart of motif class proportions"""
    counts = df['Class'].value_counts()
    plt.figure(figsize=(6,6))
    plt.pie(counts, labels=counts.index, autopct='%1.1f%%', wedgeprops=dict(width=0.5))
    plt.title("Donut Chart: Motif Class Proportion")
    plt.show()


def plot_sunburst(df):
    """Sunburst chart using Plotly"""
    fig = px.sunburst(df, path=['Class','Subclass'], title="Sunburst: Class/Subclass")
    return fig


def plot_treemap(df):
    """Treemap using Plotly"""
    fig = px.treemap(df, path=['Class','Subclass'], title="Treemap: Class/Subclass")
    return fig


def plot_score_distributions(df):
    """Score distribution plots"""
    # Box plot
    plt.figure(figsize=(15,8))
    plt.subplot(2,2,1)
    sns.boxplot(data=df, x='Class', y='Actual_Score')
    plt.title("Score Distribution by Motif Class")
    plt.xticks(rotation=45)
    
    # Violin plot
    plt.subplot(2,2,2)
    sns.violinplot(data=df, x='Class', y='Actual_Score')
    plt.title("Violin Plot: Score Distribution")
    plt.xticks(rotation=45)
    
    # Histogram
    plt.subplot(2,2,3)
    plt.hist(df['Actual_Score'], bins=30, alpha=0.7, edgecolor='black')
    plt.title("Score Distribution Histogram")
    plt.xlabel("Actual Score")
    plt.ylabel("Frequency")
    
    # Subclass score comparison
    plt.subplot(2,2,4)
    top_subclasses = df['Subclass'].value_counts().head(8).index
    df_filtered = df[df['Subclass'].isin(top_subclasses)]
    sns.boxplot(data=df_filtered, x='Subclass', y='Actual_Score')
    plt.title("Score by Top Subclasses")
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.show()


def plot_motif_genomic_distribution(df):
    """Plot motif distribution along genomic coordinates"""
    plt.figure(figsize=(15,8))
    
    # Scatter plot of motif positions
    plt.subplot(2,2,1)
    colors = plt.cm.tab10(np.linspace(0, 1, len(df['Class'].unique())))
    class_colors = dict(zip(df['Class'].unique(), colors))
    
    for cls in df['Class'].unique():
        cls_data = df[df['Class'] == cls]
        plt.scatter(cls_data['Start'], cls_data['Length'], 
                   alpha=0.6, label=cls, color=class_colors[cls])
    plt.xlabel('Genomic Position (Start)')
    plt.ylabel('Motif Length')
    plt.title('Motif Length vs Genomic Position')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Density plot
    plt.subplot(2,2,2)
    plt.hist(df['Start'], bins=50, alpha=0.7, edgecolor='black')
    plt.xlabel('Genomic Position')
    plt.ylabel('Frequency')
    plt.title('Motif Position Distribution')
    
    # Length distribution by class
    plt.subplot(2,2,3)
    sns.boxplot(data=df, x='Class', y='Length')
    plt.xticks(rotation=45)
    plt.title('Motif Length by Class')
    
    # GC content vs Score
    plt.subplot(2,2,4)
    plt.scatter(df['GC_Content'], df['Actual_Score'], alpha=0.6)
    plt.xlabel('GC Content (%)')
    plt.ylabel('Actual Score')
    plt.title('Score vs GC Content')
    
    plt.tight_layout()
    plt.show()


def plot_class_subclass_heatmap(df):
    """Create a heatmap showing class-subclass relationships"""
    pivot_table = df.groupby(['Class', 'Subclass']).size().unstack(fill_value=0)
    
    plt.figure(figsize=(15,8))
    sns.heatmap(pivot_table, annot=True, fmt='d', cmap='RdYlGn', 
                cbar_kws={'label': 'Count'})
    plt.title('Class-Subclass Distribution Heatmap')
    plt.xlabel('Subclass')
    plt.ylabel('Class')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()


def plot_sequence_coverage(df):
    """Plot sequence coverage by motifs"""
    plt.figure(figsize=(12,6))
    
    # Calculate coverage per sequence
    coverage_data = []
    for seq_name in df['Sequence_Name'].unique():
        seq_data = df[df['Sequence_Name'] == seq_name]
        total_length = seq_data['End'].max() if not seq_data.empty else 0
        covered_positions = set()
        for _, row in seq_data.iterrows():
            covered_positions.update(range(int(row['Start']), int(row['End'])+1))
        coverage_pct = (len(covered_positions) / total_length * 100) if total_length else 0
        coverage_data.append({
            'Sequence': seq_name,
            'Coverage_Percent': coverage_pct,
            'Motif_Count': len(seq_data),
            'Total_Length': total_length
        })
    
    coverage_df = pd.DataFrame(coverage_data)
    
    plt.subplot(1,2,1)
    plt.bar(coverage_df['Sequence'], coverage_df['Coverage_Percent'])
    plt.xlabel('Sequence Name')
    plt.ylabel('Coverage (%)')
    plt.title('Motif Coverage by Sequence')
    plt.xticks(rotation=45)
    
    plt.subplot(1,2,2)
    plt.scatter(coverage_df['Total_Length'], coverage_df['Motif_Count'])
    plt.xlabel('Sequence Length')
    plt.ylabel('Motif Count')
    plt.title('Motif Count vs Sequence Length')
    
    plt.tight_layout()
    plt.show()


def create_interactive_motif_browser(df):
    """Create interactive plotly visualization for motif browsing"""
    # Ensure scores are positive for size parameter
    df_plot = df.copy()
    df_plot['Score_Size'] = np.maximum(df_plot['Actual_Score'], 0.1)  # Minimum size of 0.1
    
    # Interactive scatter plot
    fig = px.scatter(df_plot, x='Start', y='Length', color='Class', 
                    size='Score_Size', hover_data=['Subclass', 'GC_Content'],
                    title='Interactive Motif Browser')
    fig.update_layout(height=600)
    return fig


def plot_scoring_method_comparison(df):
    """Compare different scoring methods"""
    if 'Scoring_Method' not in df.columns:
        return
        
    plt.figure(figsize=(12,6))
    
    plt.subplot(1,2,1)
    sns.boxplot(data=df, x='Scoring_Method', y='Actual_Score')
    plt.title('Score Distribution by Scoring Method')
    plt.xticks(rotation=45)
    
    plt.subplot(1,2,2)
    scoring_counts = df['Scoring_Method'].value_counts()
    plt.pie(scoring_counts.values, labels=scoring_counts.index, autopct='%1.1f%%')
    plt.title('Usage of Scoring Methods')
    
    plt.tight_layout()
    plt.show()


def plot_cdf(df):
    """Cumulative distribution function"""
    scores_sorted = np.sort(df['Actual_Score'])
    cdf = np.arange(1, len(scores_sorted)+1)/len(scores_sorted)
    plt.figure(figsize=(6,4))
    plt.plot(scores_sorted, cdf)
    plt.xlabel('Actual Score')
    plt.ylabel('CDF')
    plt.title("CDF of Motif Scores")
    plt.tight_layout()
    plt.show()


def plot_motif_tracks(df, seq_length=2000):
    """Lollipop/Track plot"""
    motif_classes = df['Class'].unique()
    plt.figure(figsize=(10,3))
    for i, cl in enumerate(motif_classes):
        hits = df[df['Class']==cl]
        plt.hlines(i, 0, seq_length, color='gray', alpha=0.15)
        plt.scatter(hits['Start'], [i]*len(hits), label=cl, s=60, alpha=0.7)
    plt.yticks(range(len(motif_classes)), motif_classes)
    plt.xlabel("Sequence Position")
    plt.title("Lollipop/Track")
    plt.tight_layout()
    plt.show()


def plot_density_heatmap(df, seq_length=2000):
    """Motif density heatmap"""
    density = np.zeros(seq_length)
    for _, row in df.iterrows():
        start_idx = max(0, int(row['Start'])-1)
        end_idx = min(seq_length, int(row['End']))
        density[start_idx:end_idx] += 1
    
    plt.figure(figsize=(10,2))
    plt.imshow(density[np.newaxis,:], aspect='auto', cmap='RdYlGn', extent=[0,seq_length,0,1])
    plt.xlabel("Position")
    plt.yticks([])
    plt.title("Motif Density Heatmap")
    plt.tight_layout()
    plt.show()


def plot_cluster_density(df):
    """Cluster density per sequence"""
    density_data = df.groupby('Sequence_Name').size()
    plt.figure(figsize=(8,4))
    density_data.plot(kind='barh', color='teal', alpha=0.7)
    plt.xlabel("Motif Count")
    plt.title("Cluster Density per Sequence")
    plt.tight_layout()
    plt.show()


def plot_venn_diagram(df):
    """Venn diagram for motif overlaps"""
    # Example: G-Quadruplex vs Z-DNA
    g4_seqs = set(df[df['Class']=='G-Quadruplex']['Sequence_Name'])
    zdna_seqs = set(df[df['Class']=='Z-DNA']['Sequence_Name'])
    
    if len(g4_seqs) > 0 and len(zdna_seqs) > 0:
        plt.figure(figsize=(6,4))
        venn2([g4_seqs, zdna_seqs], set_labels=('G-Quadruplex','Z-DNA'))
        plt.title("Venn: Sequence Overlap")
        plt.show()


def plot_network_graph(df):
    """Network graph of motif interactions"""
    motif_classes = df['Class'].unique()
    # Generate random co-occurrences for demonstration
    pairs = [(random.choice(motif_classes), random.choice(motif_classes)) for _ in range(50)]
    mtrx = pd.crosstab(pd.DataFrame(pairs)[0], pd.DataFrame(pairs)[1])
    
    G = nx.from_pandas_adjacency(mtrx, create_using=nx.Graph())
    plt.figure(figsize=(8,6))
    nx.draw(G, with_labels=True, node_color='skyblue', node_size=2000, edge_color='gray')
    plt.title("Motif Interaction Network")
    plt.show()


def plot_gc_content_scatter(df):
    """GC content by motif position"""
    plt.figure(figsize=(8,3))
    plt.scatter(df['Start'], df['GC_Content'], alpha=0.5)
    plt.xlabel("Position")
    plt.ylabel("GC%")
    plt.title("GC Content by Motif Position")
    plt.tight_layout()
    plt.show()


def plot_tsne(df):
    """t-SNE dimensionality reduction"""
    features = df[['Actual_Score','Length','GC_Content']].fillna(0)
    tsne = TSNE(n_components=2, random_state=42)
    tsne_result = tsne.fit_transform(features)
    
    df_plot = df.copy()
    df_plot['TSNE1'] = tsne_result[:, 0]
    df_plot['TSNE2'] = tsne_result[:, 1]
    
    plt.figure(figsize=(8,6))
    sns.scatterplot(data=df_plot, x='TSNE1', y='TSNE2', hue='Class', palette='tab10')
    plt.title("t-SNE: Motif Feature Clustering")
    plt.tight_layout()
    plt.show()


def plot_manhattan(df):
    """Manhattan-like plot"""
    motif_classes = df['Class'].unique()
    plt.figure(figsize=(10,3))
    for cl in motif_classes:
        hits = df[df['Class']==cl]
        plt.scatter(hits['Start'], hits['Actual_Score'], label=cl, alpha=0.6)
    plt.xlabel("Position")
    plt.ylabel("Actual Score")
    plt.title("Manhattan Plot: Motif Scores")
    plt.legend(ncol=3)
    plt.tight_layout()
    plt.show()


def plot_interactive_track(df):
    """Interactive track plot using Plotly"""
    # Ensure scores are positive for size parameter
    df_plot = df.copy()
    df_plot['Score_Size'] = np.maximum(df_plot['Actual_Score'], 0.1)  # Minimum size of 0.1
    
    fig = px.scatter(df_plot, x='Start', y='Class', size='Score_Size', color='Subclass',
                    hover_data=['Sequence_Name','Length','GC_Content'], 
                    title="Interactive Motif Track")
    return fig


def create_all_visualizations(df=None, save_plots=False, output_dir='./plots/'):
    """Create all available visualizations for the motif data"""
    import os
    if save_plots and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    if df is None:
        df = generate_pseudodata()
    
    print("Creating comprehensive motif visualizations...")
    
    try:
        print("1. Motif counts...")
        plot_motif_counts(df)
        
        print("2. Stacked distribution...")
        plot_stacked_distribution(df)
        
        print("3. Pie chart...")
        plot_pie_chart(df)
        
        print("4. Sunburst chart...")
        plot_sunburst(df)
        
        print("5. Treemap...")
        plot_treemap(df)
        
        print("6. Score distributions...")
        plot_score_distributions(df)
        
        print("7. CDF plot...")
        plot_cdf(df)
        
        print("8. Genomic distribution...")
        plot_motif_genomic_distribution(df)
        
        print("9. Class-subclass heatmap...")
        plot_class_subclass_heatmap(df)
        
        print("10. Sequence coverage...")
        plot_sequence_coverage(df)
        
        print("11. Motif tracks...")
        plot_motif_tracks(df)
        
        print("12. Density heatmap...")
        plot_density_heatmap(df)
        
        print("13. Cluster density...")
        plot_cluster_density(df)
        
        print("14. Venn diagram...")
        plot_venn_diagram(df)
        
        print("15. Network graph...")
        plot_network_graph(df)
        
        print("16. GC content scatter...")
        plot_gc_content_scatter(df)
        
        print("17. t-SNE plot...")
        plot_tsne(df)
        
        print("18. Manhattan plot...")
        plot_manhattan(df)
        
        print("19. Interactive track...")
        plot_interactive_track(df)
        
        print("20. Interactive browser...")
        create_interactive_motif_browser(df)
        
        print("21. Scoring method comparison...")
        plot_scoring_method_comparison(df)
        
        print("✓ All visualizations created successfully!")
        
    except Exception as e:
        print(f"Error creating visualizations: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    # Demo with pseudodata
    demo_df = generate_pseudodata()
    print("Generated demo data with all 22 subclasses:")
    print(demo_df.head(10))
    print(f"\nUnique classes: {sorted(demo_df['Class'].unique())}")
    print(f"Unique subclasses: {sorted(demo_df['Subclass'].unique())}")
    print(f"Total subclasses count: {len(demo_df['Subclass'].unique())}")
    create_all_visualizations(demo_df)