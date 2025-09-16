"""
Enhanced NBDFinder Visualization Suite - Information-Type Based
==============================================================

This module provides comprehensive, automatically-generated visualizations organized by
information type rather than plot type, as required by the problem statement.

Key Features:
- Automatic generation of all plots (no user interaction required)
- Information-type organization (Coverage, Distribution, Sequence Analysis, etc.)
- Prominent motif coverage and non-B DNA density reporting
- Retains Intel Hyperscan approach and professional styling
- Rigorous testing and validation

Author: Enhanced by Copilot based on Dr. Venkata Rajesh Yella's original work
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import networkx as nx
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# Set professional styling
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

class InformationBasedVisualizer:
    """
    Information-type based visualization generator for Non-B DNA motifs.
    Automatically generates comprehensive analysis without user interaction.
    """
    
    def __init__(self, motifs_df, sequence_length, sequence_name="Sequence"):
        self.df = motifs_df
        self.seq_length = sequence_length
        self.seq_name = sequence_name
        self.coverage_stats = self._calculate_comprehensive_stats()
        
    def _calculate_comprehensive_stats(self):
        """Calculate comprehensive coverage and density statistics"""
        if self.df.empty:
            return {
                'motif_coverage_pct': 0,
                'non_b_dna_density': 0,
                'covered_positions': 0,
                'total_motifs': 0,
                'motif_density_per_kb': 0,
                'avg_motif_length': 0
            }
            
        # Calculate covered positions
        covered = set()
        for _, row in self.df.iterrows():
            covered.update(range(int(row['Start']), int(row['End']) + 1))
        
        coverage_pct = (len(covered) / self.seq_length * 100) if self.seq_length > 0 else 0
        total_motifs = len(self.df)
        motif_density_per_kb = (total_motifs / self.seq_length * 1000) if self.seq_length > 0 else 0
        avg_length = self.df['Length'].mean() if not self.df.empty else 0
        
        return {
            'motif_coverage_pct': round(coverage_pct, 2),
            'non_b_dna_density': round(motif_density_per_kb, 2),
            'covered_positions': len(covered),
            'total_motifs': total_motifs,
            'motif_density_per_kb': round(motif_density_per_kb, 2),
            'avg_motif_length': round(avg_length, 2)
        }
    
    def create_coverage_analysis(self):
        """
        INFORMATION TYPE: Coverage & Density Analysis
        Generate comprehensive coverage and density visualizations
        """
        plots = {}
        
        # Handle empty dataframe case
        if self.df.empty:
            fig, ax = plt.subplots(1, 1, figsize=(10, 6))
            ax.text(0.5, 0.5, f'No motifs found in {self.seq_name}\nCoverage: 0% | Density: 0 motifs/kb', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=16)
            ax.set_title('Coverage & Density Analysis - No Motifs Found', fontweight='bold')
            plots['coverage_analysis'] = fig
            plots['detailed_coverage_map'] = fig
            return plots
        
        # 1. Coverage Overview Dashboard
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))
        fig.suptitle(f'Coverage & Density Analysis - {self.seq_name}', fontsize=16, fontweight='bold')
        
        # Coverage percentage gauge
        coverage_pct = self.coverage_stats['motif_coverage_pct']
        ax1.pie([coverage_pct, 100-coverage_pct], labels=['Covered', 'Uncovered'], 
                autopct='%1.1f%%', startangle=90, colors=['#ff6b6b', '#e9e9e9'])
        ax1.set_title(f'Motif Coverage: {coverage_pct}%', fontweight='bold')
        
        # Density per class
        class_density = self.df.groupby('Class').size().sort_values(ascending=False)
        ax2.bar(range(len(class_density)), class_density.values, color='skyblue')
        ax2.set_xticks(range(len(class_density)))
        ax2.set_xticklabels(class_density.index, rotation=60, ha='right', fontsize=10)
        ax2.set_title(f'Non-B DNA Density by Class\n(Total: {self.coverage_stats["non_b_dna_density"]:.2f} motifs/kb)', fontweight='bold')
        ax2.set_ylabel('Count')
        
        # Coverage heatmap along sequence
        if not self.df.empty:
            bins = min(100, self.seq_length // 10)
            hist, edges = np.histogram(self.df['Start'], bins=bins, range=(0, self.seq_length))
            im = ax3.imshow([hist], aspect='auto', cmap='hot', extent=[0, self.seq_length, -0.5, 0.5])
            ax3.set_title('Motif Density Heatmap Along Sequence', fontweight='bold')
            ax3.set_xlabel('Sequence Position (bp)')
            ax3.set_yticks([])
            plt.colorbar(im, ax=ax3, label='Motif Count')
        else:
            ax3.text(0.5, 0.5, 'No motifs found', ha='center', va='center', transform=ax3.transAxes)
            ax3.set_title('Motif Density Heatmap Along Sequence', fontweight='bold')
        
        # Length distribution impact on coverage
        if not self.df.empty:
            ax4.scatter(self.df['Length'], self.df['Actual_Score'], alpha=0.6, c=self.df['Start'], cmap='viridis')
            ax4.set_xlabel('Motif Length (bp)')
            ax4.set_ylabel('Score')
            ax4.set_title('Length vs Score (colored by position)', fontweight='bold')
            cbar = plt.colorbar(ax4.collections[0], ax=ax4, label='Position')
        else:
            ax4.text(0.5, 0.5, 'No motifs found', ha='center', va='center', transform=ax4.transAxes)
            ax4.set_title('Length vs Score Analysis', fontweight='bold')
        
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.4, wspace=0.3)  # Add extra spacing
        plots['coverage_analysis'] = fig
        
        # 2. Detailed Coverage Map
        fig2, ax = plt.subplots(1, 1, figsize=(18, 10))
        
        if not self.df.empty:
            # Create detailed coverage visualization
            classes = self.df['Class'].unique()
            colors = plt.cm.Set3(np.linspace(0, 1, len(classes)))
            class_colors = dict(zip(classes, colors))
            
            for i, cls in enumerate(classes):
                cls_motifs = self.df[self.df['Class'] == cls]
                y_pos = i
                for _, motif in cls_motifs.iterrows():
                    width = motif['End'] - motif['Start']
                    rect = patches.Rectangle((motif['Start'], y_pos-0.4), width, 0.8, 
                                           facecolor=class_colors[cls], alpha=0.7, edgecolor='black')
                    ax.add_patch(rect)
            
            ax.set_xlim(0, self.seq_length)
            ax.set_ylim(-0.5, len(classes)-0.5)
            ax.set_yticks(range(len(classes)))
            ax.set_yticklabels(classes, fontsize=10)
            ax.set_xlabel('Sequence Position (bp)')
            ax.set_title(f'Detailed Motif Coverage Map - {self.seq_name}\nCoverage: {coverage_pct}% | Density: {self.coverage_stats["non_b_dna_density"]:.2f} motifs/kb', 
                        fontweight='bold', fontsize=14)
            
            # Add coverage statistics text
            stats_text = f"""
            Total Motifs: {self.coverage_stats['total_motifs']}
            Covered Positions: {self.coverage_stats['covered_positions']} bp
            Average Motif Length: {self.coverage_stats['avg_motif_length']:.1f} bp
            """
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10, 
                   verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        else:
            ax.text(0.5, 0.5, 'No motifs found for coverage analysis', ha='center', va='center', 
                   transform=ax.transAxes, fontsize=16)
            ax.set_title(f'Motif Coverage Map - {self.seq_name}', fontweight='bold')
        
        plt.tight_layout()
        plots['detailed_coverage_map'] = fig2
        
        return plots
    
    def create_distribution_analysis(self):
        """
        INFORMATION TYPE: Distribution Analysis  
        Generate comprehensive distribution visualizations
        """
        plots = {}
        
        if self.df.empty:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, f'No motifs found in {self.seq_name} for distribution analysis', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=16)
            ax.set_title('Distribution Analysis - No Motifs Found', fontweight='bold')
            plots['distribution_analysis'] = fig
            return plots
        
        # 1. Comprehensive Distribution Dashboard
        fig = plt.figure(figsize=(22, 18))
        gs = fig.add_gridspec(4, 3, hspace=0.4, wspace=0.4)
        
        # Class distribution
        ax1 = fig.add_subplot(gs[0, 0])
        class_counts = self.df['Class'].value_counts()
        ax1.pie(class_counts.values, labels=class_counts.index, autopct='%1.1f%%', startangle=90)
        ax1.set_title('Class Distribution', fontweight='bold')
        
        # Subclass distribution (top 10)
        ax2 = fig.add_subplot(gs[0, 1])
        subclass_counts = self.df['Subclass'].value_counts().head(10)
        ax2.barh(range(len(subclass_counts)), subclass_counts.values)
        ax2.set_yticks(range(len(subclass_counts)))
        ax2.set_yticklabels(subclass_counts.index, fontsize=9)
        ax2.set_title('Top 10 Subclass Distribution', fontweight='bold')
        ax2.set_xlabel('Count')
        
        # Length distribution
        ax3 = fig.add_subplot(gs[0, 2])
        ax3.hist(self.df['Length'], bins=30, alpha=0.7, edgecolor='black')
        ax3.set_xlabel('Motif Length (bp)')
        ax3.set_ylabel('Frequency')
        ax3.set_title('Length Distribution', fontweight='bold')
        ax3.axvline(self.df['Length'].mean(), color='red', linestyle='--', 
                   label=f'Mean: {self.df["Length"].mean():.1f}')
        ax3.legend()
        
        # Score distribution by class
        ax4 = fig.add_subplot(gs[1, :])
        sns.boxplot(data=self.df, x='Class', y='Actual_Score', ax=ax4)
        ax4.set_xticklabels(ax4.get_xticklabels(), rotation=60, ha='right', fontsize=10)
        ax4.set_title('Score Distribution by Class', fontweight='bold')
        ax4.set_ylabel('Actual Score')
        
        # Class-Subclass heatmap
        ax5 = fig.add_subplot(gs[2, :])
        class_subclass = pd.crosstab(self.df['Class'], self.df['Subclass'])
        sns.heatmap(class_subclass, annot=True, fmt='d', cmap='YlOrRd', ax=ax5)
        ax5.set_title('Class-Subclass Distribution Heatmap', fontweight='bold')
        ax5.set_xlabel('Subclass')
        ax5.set_ylabel('Class')
        
        # Position distribution analysis
        ax6 = fig.add_subplot(gs[3, 0])
        ax6.hist(self.df['Start'], bins=50, alpha=0.7, edgecolor='black')
        ax6.set_xlabel('Start Position')
        ax6.set_ylabel('Frequency')
        ax6.set_title('Position Distribution', fontweight='bold')
        
        # GC content distribution
        ax7 = fig.add_subplot(gs[3, 1])
        if 'GC_Content' in self.df.columns:
            ax7.hist(self.df['GC_Content'], bins=30, alpha=0.7, edgecolor='black')
            ax7.set_xlabel('GC Content (%)')
            ax7.set_ylabel('Frequency')
            ax7.set_title('GC Content Distribution', fontweight='bold')
        else:
            ax7.text(0.5, 0.5, 'GC Content\ndata not available', ha='center', va='center', 
                    transform=ax7.transAxes)
            ax7.set_title('GC Content Distribution', fontweight='bold')
        
        # Normalized score distribution
        ax8 = fig.add_subplot(gs[3, 2])
        ax8.hist(self.df['Normalized_Score'], bins=30, alpha=0.7, edgecolor='black')
        ax8.set_xlabel('Normalized Score')
        ax8.set_ylabel('Frequency')
        ax8.set_title('Normalized Score Distribution', fontweight='bold')
        
        fig.suptitle(f'Comprehensive Distribution Analysis - {self.seq_name}', fontsize=16, fontweight='bold')
        plots['distribution_analysis'] = fig
        
        return plots
    
    def create_sequence_analysis(self):
        """
        INFORMATION TYPE: Sequence Analysis
        Generate sequence-focused visualizations
        """
        plots = {}
        
        if self.df.empty:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, f'No motifs found in {self.seq_name} for sequence analysis', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=16)
            ax.set_title('Sequence Analysis - No Motifs Found', fontweight='bold')
            plots['sequence_analysis'] = fig
            return plots
        
        # 1. Motif Tracks and Overlaps
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(18, 14))
        
        # Track plot
        classes = self.df['Class'].unique()
        colors = plt.cm.tab20(np.linspace(0, 1, len(classes)))
        class_colors = dict(zip(classes, colors))
        
        for i, cls in enumerate(classes):
            cls_motifs = self.df[self.df['Class'] == cls]
            y_pos = i
            ax1.scatter(cls_motifs['Start'], [y_pos] * len(cls_motifs), 
                       c=[class_colors[cls]], s=60, alpha=0.7, label=cls)
        
        ax1.set_xlim(0, self.seq_length)
        ax1.set_yticks(range(len(classes)))
        ax1.set_yticklabels(classes, fontsize=10)
        ax1.set_xlabel('Sequence Position (bp)')
        ax1.set_title('Motif Track Plot - All Classes', fontweight='bold')
        ax1.grid(True, alpha=0.3)
        
        # Overlap analysis
        overlaps = []
        motifs_list = self.df.to_dict('records')
        for i, m1 in enumerate(motifs_list):
            for j, m2 in enumerate(motifs_list[i+1:], i+1):
                if (m1['Start'] <= m2['End'] and m2['Start'] <= m1['End']):
                    overlaps.append((m1['Class'], m2['Class']))
        
        if overlaps:
            overlap_df = pd.DataFrame(overlaps, columns=['Class1', 'Class2'])
            overlap_counts = overlap_df.groupby(['Class1', 'Class2']).size().reset_index(name='Count')
            
            # Create overlap matrix
            all_classes = list(self.df['Class'].unique())
            overlap_matrix = pd.DataFrame(0, index=all_classes, columns=all_classes)
            for _, row in overlap_counts.iterrows():
                overlap_matrix.loc[row['Class1'], row['Class2']] += row['Count']
                overlap_matrix.loc[row['Class2'], row['Class1']] += row['Count']
            
            sns.heatmap(overlap_matrix, annot=True, fmt='d', cmap='Blues', ax=ax2)
            ax2.set_title('Motif Overlap Analysis', fontweight='bold')
        else:
            ax2.text(0.5, 0.5, 'No overlapping motifs found', ha='center', va='center', 
                    transform=ax2.transAxes)
            ax2.set_title('Motif Overlap Analysis', fontweight='bold')
        
        # Clustering analysis
        if len(self.df) > 1:
            # Find clusters using distance threshold
            positions = self.df['Start'].values
            positions_sorted = np.sort(positions)
            gaps = np.diff(positions_sorted)
            threshold = np.mean(gaps) + 2 * np.std(gaps) if len(gaps) > 0 else 1000
            
            clusters = []
            current_cluster = [positions_sorted[0]]
            for i in range(1, len(positions_sorted)):
                if positions_sorted[i] - positions_sorted[i-1] <= threshold:
                    current_cluster.append(positions_sorted[i])
                else:
                    if len(current_cluster) > 1:
                        clusters.append(current_cluster)
                    current_cluster = [positions_sorted[i]]
            if len(current_cluster) > 1:
                clusters.append(current_cluster)
            
            # Plot clusters
            for i, cluster in enumerate(clusters):
                ax3.scatter(cluster, [i] * len(cluster), s=100, alpha=0.7, label=f'Cluster {i+1}')
            
            ax3.set_xlabel('Position')
            ax3.set_ylabel('Cluster')
            ax3.set_title(f'Motif Clustering Analysis (Found {len(clusters)} clusters)', fontweight='bold')
            if clusters:
                ax3.legend()
        else:
            ax3.text(0.5, 0.5, 'Insufficient motifs for clustering analysis', 
                    ha='center', va='center', transform=ax3.transAxes)
            ax3.set_title('Motif Clustering Analysis', fontweight='bold')
        
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.4)  # Add extra vertical spacing
        plots['sequence_analysis'] = fig
        
        return plots
    
    def create_comparative_analysis(self):
        """
        INFORMATION TYPE: Comparative Analysis
        Generate statistical and comparative visualizations
        """
        plots = {}
        
        if self.df.empty:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, f'No motifs found in {self.seq_name} for comparative analysis', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=16)
            ax.set_title('Comparative Analysis - No Motifs Found', fontweight='bold')
            plots['comparative_analysis'] = fig
            return plots
        
        # 1. Statistical Comparison Dashboard
        fig = plt.figure(figsize=(20, 16))
        gs = fig.add_gridspec(3, 3, hspace=0.4, wspace=0.4)
        
        # Score comparison by class
        ax1 = fig.add_subplot(gs[0, :])
        sns.violinplot(data=self.df, x='Class', y='Actual_Score', ax=ax1)
        ax1.set_xticklabels(ax1.get_xticklabels(), rotation=60, ha='right', fontsize=10)
        ax1.set_title('Score Distribution Comparison by Class', fontweight='bold')
        
        # Length vs Score scatter
        ax2 = fig.add_subplot(gs[1, 0])
        scatter = ax2.scatter(self.df['Length'], self.df['Actual_Score'], 
                             c=self.df['Start'], cmap='viridis', alpha=0.6)
        ax2.set_xlabel('Length')
        ax2.set_ylabel('Score')
        ax2.set_title('Length vs Score\n(colored by position)', fontweight='bold')
        plt.colorbar(scatter, ax=ax2, label='Position')
        
        # Class efficiency (score per bp)
        ax3 = fig.add_subplot(gs[1, 1])
        self.df['Score_per_bp'] = self.df['Actual_Score'] / self.df['Length']
        class_efficiency = self.df.groupby('Class')['Score_per_bp'].mean().sort_values(ascending=True)
        ax3.barh(range(len(class_efficiency)), class_efficiency.values)
        ax3.set_yticks(range(len(class_efficiency)))
        ax3.set_yticklabels(class_efficiency.index, fontsize=9)
        ax3.set_xlabel('Score per bp')
        ax3.set_title('Class Efficiency\n(Score per bp)', fontweight='bold')
        
        # Position preferences
        ax4 = fig.add_subplot(gs[1, 2])
        position_bins = pd.cut(self.df['Start'], bins=10)
        position_counts = position_bins.value_counts().sort_index()
        bin_centers = [interval.mid for interval in position_counts.index]
        ax4.bar(range(len(bin_centers)), position_counts.values)
        ax4.set_xticks(range(len(bin_centers)))
        ax4.set_xticklabels([f'{int(bc)}' for bc in bin_centers], rotation=60, fontsize=9)
        ax4.set_xlabel('Position (bin centers)')
        ax4.set_ylabel('Count')
        ax4.set_title('Position Preferences', fontweight='bold')
        
        # Statistical summary table
        ax5 = fig.add_subplot(gs[2, :])
        ax5.axis('tight')
        ax5.axis('off')
        
        # Create summary statistics
        summary_stats = []
        for cls in self.df['Class'].unique():
            cls_data = self.df[self.df['Class'] == cls]
            stats = {
                'Class': cls,
                'Count': len(cls_data),
                'Avg_Length': cls_data['Length'].mean(),
                'Avg_Score': cls_data['Actual_Score'].mean(),
                'Score_Std': cls_data['Actual_Score'].std(),
                'Min_Pos': cls_data['Start'].min(),
                'Max_Pos': cls_data['Start'].max(),
                'Coverage': (cls_data['Length'].sum() / self.seq_length * 100)
            }
            summary_stats.append(stats)
        
        stats_df = pd.DataFrame(summary_stats)
        stats_df = stats_df.round(2)
        
        table = ax5.table(cellText=stats_df.values, colLabels=stats_df.columns,
                         cellLoc='center', loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1.2, 1.5)
        ax5.set_title('Statistical Summary by Class', fontweight='bold', pad=20)
        
        fig.suptitle(f'Comparative Analysis - {self.seq_name}', fontsize=16, fontweight='bold')
        plots['comparative_analysis'] = fig
        
        return plots
    
    def create_advanced_analysis(self):
        """
        INFORMATION TYPE: Advanced Analysis
        Generate network, dimensionality reduction, and advanced visualizations
        """
        plots = {}
        
        if self.df.empty:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, f'No motifs found in {self.seq_name} for advanced analysis', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=16)
            ax.set_title('Advanced Analysis - No Motifs Found', fontweight='bold')
            plots['advanced_analysis'] = fig
            return plots
        
        # 1. Network and Dimensionality Analysis
        fig = plt.figure(figsize=(22, 14))
        gs = fig.add_gridspec(2, 3, hspace=0.4, wspace=0.4)
        
        # Network analysis of class co-occurrences
        ax1 = fig.add_subplot(gs[0, 0])
        
        # Create co-occurrence network
        classes = self.df['Class'].unique()
        co_occurrence = np.zeros((len(classes), len(classes)))
        class_to_idx = {cls: i for i, cls in enumerate(classes)}
        
        # Calculate co-occurrences based on proximity
        for i, m1 in self.df.iterrows():
            for j, m2 in self.df.iterrows():
                if i != j:
                    distance = abs(m1['Start'] - m2['Start'])
                    if distance < 1000:  # Within 1kb
                        idx1, idx2 = class_to_idx[m1['Class']], class_to_idx[m2['Class']]
                        co_occurrence[idx1, idx2] += 1
        
        # Create network graph
        G = nx.Graph()
        for i, cls1 in enumerate(classes):
            G.add_node(cls1)
            for j, cls2 in enumerate(classes):
                if i < j and co_occurrence[i, j] > 0:
                    G.add_edge(cls1, cls2, weight=co_occurrence[i, j])
        
        if G.edges():
            pos = nx.spring_layout(G)
            nx.draw(G, pos, ax=ax1, with_labels=True, node_color='lightblue', 
                   node_size=1000, font_size=8, font_weight='bold')
            ax1.set_title('Class Co-occurrence Network\n(within 1kb)', fontweight='bold')
        else:
            ax1.text(0.5, 0.5, 'No co-occurrences found', ha='center', va='center', 
                    transform=ax1.transAxes)
            ax1.set_title('Class Co-occurrence Network', fontweight='bold')
        
        # t-SNE analysis
        ax2 = fig.add_subplot(gs[0, 1])
        if len(self.df) > 2:
            # Prepare features for t-SNE
            features = ['Start', 'Length', 'Actual_Score']
            if 'GC_Content' in self.df.columns:
                features.append('GC_Content')
            
            X = self.df[features].fillna(0)
            if X.shape[1] >= 2 and len(X) > 3:
                scaler = StandardScaler()
                X_scaled = scaler.fit_transform(X)
                
                tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, len(X)-1))
                X_tsne = tsne.fit_transform(X_scaled)
                
                classes_for_color = self.df['Class'].values
                unique_classes = np.unique(classes_for_color)
                colors = plt.cm.tab10(np.linspace(0, 1, len(unique_classes)))
                
                for i, cls in enumerate(unique_classes):
                    mask = classes_for_color == cls
                    ax2.scatter(X_tsne[mask, 0], X_tsne[mask, 1], 
                              c=[colors[i]], label=cls, alpha=0.7)
                
                ax2.set_xlabel('t-SNE 1')
                ax2.set_ylabel('t-SNE 2')
                ax2.set_title('t-SNE Clustering\n(by features)', fontweight='bold')
                ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
            else:
                ax2.text(0.5, 0.5, 'Insufficient features\nfor t-SNE analysis', 
                        ha='center', va='center', transform=ax2.transAxes)
                ax2.set_title('t-SNE Clustering', fontweight='bold')
        else:
            ax2.text(0.5, 0.5, 'Insufficient data\nfor t-SNE analysis', 
                    ha='center', va='center', transform=ax2.transAxes)
            ax2.set_title('t-SNE Clustering', fontweight='bold')
        
        # Manhattan-style plot
        ax3 = fig.add_subplot(gs[0, 2])
        classes = self.df['Class'].unique()
        colors = plt.cm.tab10(np.linspace(0, 1, len(classes)))
        for i, cls in enumerate(classes):
            cls_data = self.df[self.df['Class'] == cls]
            ax3.scatter(cls_data['Start'], cls_data['Actual_Score'], 
                       c=[colors[i]], label=cls, alpha=0.6)
        
        ax3.set_xlabel('Position')
        ax3.set_ylabel('Score')
        ax3.set_title('Manhattan-style Plot\n(Score vs Position)', fontweight='bold')
        ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        
        # Advanced statistics heatmap
        ax4 = fig.add_subplot(gs[1, :])
        
        # Create advanced statistics matrix
        advanced_stats = []
        for cls in classes:
            cls_data = self.df[self.df['Class'] == cls]
            stats = {
                'Class': cls,
                'Count': len(cls_data),
                'Mean_Score': cls_data['Actual_Score'].mean(),
                'Std_Score': cls_data['Actual_Score'].std(),
                'Mean_Length': cls_data['Length'].mean(),
                'Std_Length': cls_data['Length'].std(),
                'Position_Spread': cls_data['Start'].max() - cls_data['Start'].min(),
                'Score_CV': cls_data['Actual_Score'].std() / cls_data['Actual_Score'].mean() if cls_data['Actual_Score'].mean() > 0 else 0
            }
            advanced_stats.append(stats)
        
        adv_df = pd.DataFrame(advanced_stats).set_index('Class')
        
        # Normalize for heatmap
        adv_df_norm = (adv_df - adv_df.min()) / (adv_df.max() - adv_df.min())
        
        sns.heatmap(adv_df_norm.T, annot=True, fmt='.2f', cmap='RdYlBu_r', ax=ax4)
        ax4.set_title('Advanced Statistics Heatmap (Normalized)', fontweight='bold')
        ax4.set_ylabel('Metrics')
        ax4.set_xlabel('Class')
        
        fig.suptitle(f'Advanced Analysis - {self.seq_name}', fontsize=16, fontweight='bold')
        plots['advanced_analysis'] = fig
        
        return plots
    
    def create_interactive_plots(self):
        """Generate interactive Plotly visualizations"""
        interactive_plots = {}
        
        if self.df.empty:
            return interactive_plots
        
        # 1. Interactive Motif Browser
        # Create a safe size column (ensure positive values for marker size)
        df_for_plot = self.df.copy()
        df_for_plot['Size_Safe'] = np.maximum(df_for_plot['Actual_Score'].abs(), 0.1)  # Ensure positive, minimum 0.1
        
        fig1 = px.scatter(df_for_plot, x='Start', y='Length', 
                         color='Class', size='Size_Safe',
                         hover_data=['Subclass', 'End', 'Normalized_Score', 'Actual_Score'],
                         title=f'Interactive Motif Browser - {self.seq_name}')
        fig1.update_layout(height=600)
        interactive_plots['interactive_browser'] = fig1
        
        # 2. Sunburst Chart
        fig2 = px.sunburst(self.df, path=['Class', 'Subclass'], 
                          title=f'Class/Subclass Hierarchy - {self.seq_name}')
        interactive_plots['sunburst'] = fig2
        
        return interactive_plots
    
    def generate_comprehensive_report(self):
        """
        Generate all visualizations automatically without user interaction.
        This is the main function that implements the problem statement requirements.
        """
        print(f"🎯 Generating comprehensive information-based visualizations for {self.seq_name}...")
        print(f"📊 Motif Coverage: {self.coverage_stats['motif_coverage_pct']}%")
        print(f"🧬 Non-B DNA Density: {self.coverage_stats['non_b_dna_density']:.2f} motifs/kb")
        print("=" * 80)
        
        all_plots = {}
        
        # Generate all information-type based visualizations
        all_plots.update(self.create_coverage_analysis())
        all_plots.update(self.create_distribution_analysis())
        all_plots.update(self.create_sequence_analysis())
        all_plots.update(self.create_comparative_analysis())
        all_plots.update(self.create_advanced_analysis())
        
        # Generate interactive plots
        interactive_plots = self.create_interactive_plots()
        
        print(f"✅ Generated {len(all_plots)} static plots and {len(interactive_plots)} interactive plots")
        print("📈 All visualizations organized by information type, not plot type")
        
        return all_plots, interactive_plots


def create_comprehensive_information_based_visualizations(motifs_df, sequence_length, sequence_name="Sequence"):
    """
    Main function to create comprehensive, information-based visualizations automatically.
    
    This function implements the key requirements from the problem statement:
    1. Generates all plots automatically without user interaction
    2. Organizes by information type rather than plot type  
    3. Reports motif coverage and non-B DNA density prominently
    4. Retains professional styling and approach
    
    Args:
        motifs_df: DataFrame with motif data
        sequence_length: Length of the analyzed sequence
        sequence_name: Name of the sequence
        
    Returns:
        tuple: (static_plots_dict, interactive_plots_dict, coverage_stats)
    """
    visualizer = InformationBasedVisualizer(motifs_df, sequence_length, sequence_name)
    static_plots, interactive_plots = visualizer.generate_comprehensive_report()
    
    return static_plots, interactive_plots, visualizer.coverage_stats


# For backward compatibility
def create_all_visualizations(df=None, save_plots=False, output_dir='./plots/'):
    """Backward compatibility wrapper"""
    if df is None:
        # Use pseudo data for testing
        from .visualization import generate_pseudodata
        df = generate_pseudodata()
        
    sequence_length = df['End'].max() if not df.empty else 1000
    sequence_name = df['Sequence_Name'].iloc[0] if not df.empty else "Test Sequence"
    
    return create_comprehensive_information_based_visualizations(df, sequence_length, sequence_name)