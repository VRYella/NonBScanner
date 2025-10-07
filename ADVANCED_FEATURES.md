# Advanced Features Documentation - NBDScanner

## Enhanced Hybrid and Cluster Detection

### Overview
The NBDScanner now includes significantly improved hybrid and cluster detection with the following enhancements:

### 1. Hybrid Motifs (Class 10)

**Key Improvements:**
- ✅ **Actual sequence extraction**: Hybrids now report the real DNA sequence instead of placeholder text
- ✅ **Longest non-overlapping regions**: Uses intelligent selection algorithm to avoid redundancy
- ✅ **Dual scoring**: Both raw and normalized scores (0-1 scale) provided
- ✅ **Component tracking**: Lists which motif classes participate in each hybrid

**Example Hybrid Output:**
```python
{
    'Class': 'Hybrid',
    'Subclass': 'R-Loop_Cruciform_Overlap',
    'Start': 1,
    'End': 77,
    'Length': 76,
    'Sequence': 'GGGGTTGGGGTTGGGGTTGGGGAAAAAAAAAAAAAAAAAAAAAAAAAGCG...',
    'Score': 0.193,                    # Raw score (average of components)
    'Normalized_Score': 0.0019,        # Normalized 0-1 scale
    'Component_Classes': ['R-Loop', 'Cruciform'],
    'Strand': '+',
    'Method': 'Hybrid_Detection'
}
```

**Normalization Formula:**
```
normalized = (raw_score / (raw_score + 50)) × length_factor × class_factor
where:
  length_factor = min(1.0, length / 100bp)  # Longer = more stable
  class_factor = min(1.0, num_classes / 3)   # More classes = more interesting
```

### 2. Cluster Motifs (Class 11)

**Key Improvements:**
- ✅ **Actual sequence extraction**: Clusters report the real DNA sequence of the region
- ✅ **Longest non-overlapping regions**: Intelligent selection to avoid overlap
- ✅ **Dual scoring**: Both raw density scores and normalized scores
- ✅ **Rich metadata**: Motif count, class diversity, and more

**Example Cluster Output:**
```python
{
    'Class': 'Non-B_DNA_Clusters',
    'Subclass': 'Mixed_Cluster_10_classes',
    'Start': 1,
    'End': 103,
    'Length': 102,
    'Sequence': 'GGGGTTGGGGTTGGGGTTGGGGAAAAAAAAAAAAAAAAAAAAAAAAAGCG...',
    'Score': 26.0,                     # Raw density score (motifs/100bp × 100)
    'Normalized_Score': 0.2063,        # Normalized 0-1 scale
    'Motif_Count': 26,                 # Number of motifs in cluster
    'Class_Diversity': 10,             # Number of distinct classes
    'Strand': '+',
    'Method': 'Cluster_Detection'
}
```

**Normalization Formula:**
```
normalized = (raw_score / (raw_score + 100)) × density_factor × diversity_factor
where:
  density_factor = min(1.0, motif_count / 10)      # 10 motifs = high density
  diversity_factor = min(1.0, class_diversity / 5)  # 5 classes = high diversity
```

## Score Normalization

### All Motif Classes
Every motif now includes both `Score` (raw) and `Normalized_Score` (0-1 scale):

**Class-Specific Parameters:**

| Class          | Max Typical Score | Length Weight | Notes                           |
|----------------|-------------------|---------------|---------------------------------|
| G-Quadruplex   | 3.0              | 0.7           | G4Hunter scores (-3 to +3)     |
| Curved DNA     | 1.0              | 0.5           | Curvature already normalized    |
| Z-DNA          | 100.0            | 0.6           | Z-DNA scores can be high        |
| Triplex        | 1.0              | 0.6           | Triplex potential (0-1)         |
| R-Loop         | 1.0              | 0.5           | R-loop potential                |
| i-Motif        | 2.0              | 0.6           | Similar to G4                   |
| Slipped DNA    | 1.0              | 0.8           | Instability-based               |
| Cruciform      | 1.0              | 0.7           | Structural stability            |
| A-philic       | 50.0             | 0.5           | A-philic scores                 |

**General Normalization:**
```
score_norm = raw_score / (raw_score + max_score)
length_factor = min(1.0, length / 50bp) × weight + (1 - weight)
normalized = score_norm × length_factor
```

## Advanced Visualization Suite

### Overview
Nine publication-quality static visualizations using colorblind-friendly palettes:

### 1. Genome Landscape Track
**Purpose**: Positional overview of every motif across the genome

**Features**:
- Horizontal ruler with colored glyphs (shape = Class, length = bar width)
- Motifs layered by class
- Baseline grid with major coordinate annotations
- Colorblind-safe Okabe-Ito palette

**Usage**:
```python
from advanced_visualizations import plot_genome_landscape_track

fig = plot_genome_landscape_track(motifs, sequence_length=10000)
fig.savefig('landscape.png', dpi=300, bbox_inches='tight')
```

### 2. Sliding Window Heat Ribbon
**Purpose**: Density + score trend across sequence

**Features**:
- Wide horizontal ribbon (color intensity = motif count)
- Overlaid line showing average score
- Annotated peaks with cluster IDs
- Viridis-like sequential palette

**Usage**:
```python
from advanced_visualizations import plot_sliding_window_heat_ribbon

fig = plot_sliding_window_heat_ribbon(motifs, sequence_length=10000, 
                                      window_size=1000)
fig.savefig('heat_ribbon.png', dpi=300)
```

### 3. Ridge Plots (Joyplots)
**Purpose**: Compare length distributions across motif classes

**Features**:
- Stacked density ridges ordered by median/mean
- Subtle color gradation
- Thin median lines for reference
- Clean, elegant presentation

**Usage**:
```python
from advanced_visualizations import plot_ridge_plots_length_by_class

fig = plot_ridge_plots_length_by_class(motifs)
fig.savefig('ridge_plots.png', dpi=300)
```

### 4. Sunburst / Treemap
**Purpose**: Hierarchical composition (Class → Subclass → Method)

**Features**:
- Two variants: sunburst (concentric rings) or treemap (nested rectangles)
- Sized by motif counts
- Rare subclasses grouped into "Other"
- Interactive color scheme

**Usage**:
```python
from advanced_visualizations import plot_sunburst_treemap

# Sunburst variant
fig1 = plot_sunburst_treemap(motifs, plot_type='sunburst')

# Treemap variant
fig2 = plot_sunburst_treemap(motifs, plot_type='treemap')
```

### 5. Hexbin with Marginals
**Purpose**: Show positional score structure with 2D density

**Features**:
- Central hexbin 2D density plot
- Marginal histograms (top and right)
- Log-scale option for skewed scores
- Wide format for legibility

**Usage**:
```python
from advanced_visualizations import plot_hexbin_start_vs_score

fig = plot_hexbin_start_vs_score(motifs)
fig.savefig('hexbin.png', dpi=300)
```

### 6. UpSet (Intersection) Plot
**Purpose**: Visualize common multi-motif overlaps (better than Venn)

**Features**:
- Set matrix below with vertical bars above
- Shows intersection counts
- Annotates top intersections
- Clean, publication-ready

**Usage**:
```python
from advanced_visualizations import plot_upset_intersection

fig = plot_upset_intersection(motifs)
fig.savefig('upset.png', dpi=300)
```

### 7. Score Violin + Beeswarm
**Purpose**: Show per-class score distributions with individual points

**Features**:
- Violins horizontally stacked
- Jittered points overlaid
- Top 1% points highlighted with gold/outline
- Individual outliers visible

**Usage**:
```python
from advanced_visualizations import plot_score_violin_beeswarm

fig = plot_score_violin_beeswarm(motifs)
fig.savefig('violin_beeswarm.png', dpi=300)
```

### 8. Cluster Hotspot Map
**Purpose**: Show cluster distribution and representative sequences

**Features**:
- Vertical bars for region counts
- Inset text boxes with top motif sequences/scores
- Regions ordered left→right by genomic coordinate
- Annotated top 3 hotspots

**Usage**:
```python
from advanced_visualizations import plot_cluster_hotspot_map

fig = plot_cluster_hotspot_map(motifs, sequence_length=10000)
fig.savefig('hotspot_map.png', dpi=300)
```

## Export Utilities

### Multi-Format Export
Export any plot in multiple formats simultaneously:

```python
from advanced_visualizations import export_plot

# Export in multiple formats
saved_files = export_plot(
    fig, 
    filename='my_analysis',
    formats=['png', 'svg', 'pdf'],
    dpi=300
)

# Returns: {'png': 'my_analysis.png', 'svg': 'my_analysis.svg', ...}
```

### Publication Guidelines

**For Print/Journal Submissions:**
- Use PNG @300 DPI or higher
- Keep figures under 10MB
- Use SVG for vector graphics

**For Presentations:**
- Use PNG @150 DPI
- Optimize file size
- Use landscape (16:9) aspect ratio

**Color Scheme:**
- Okabe-Ito colorblind-safe palette (7 colors)
- High contrast for B&W printing
- Consistent colors across all panels

## Complete Example

```python
from nbdscanner import MotifDetector
from advanced_visualizations import *
import matplotlib.pyplot as plt

# 1. Analyze sequence
detector = MotifDetector()
motifs = detector.analyze_sequence(my_sequence, "my_seq")

# 2. Check hybrids and clusters
hybrids = [m for m in motifs if m.get('Class') == 'Hybrid']
clusters = [m for m in motifs if m.get('Class') == 'Non-B_DNA_Clusters']

print(f"Found {len(hybrids)} hybrids and {len(clusters)} clusters")

for hybrid in hybrids:
    print(f"Hybrid: {hybrid['Subclass']}")
    print(f"  Sequence: {hybrid['Sequence'][:50]}...")
    print(f"  Raw Score: {hybrid['Score']:.3f}")
    print(f"  Normalized: {hybrid['Normalized_Score']:.4f}")

# 3. Generate all visualizations
viz_functions = [
    plot_genome_landscape_track,
    plot_sliding_window_heat_ribbon,
    plot_ridge_plots_length_by_class,
    lambda m: plot_sunburst_treemap(m, plot_type='sunburst'),
    plot_hexbin_start_vs_score,
    plot_upset_intersection,
    plot_score_violin_beeswarm,
    plot_cluster_hotspot_map,
]

for i, viz_func in enumerate(viz_functions):
    try:
        if 'sequence_length' in viz_func.__code__.co_varnames:
            fig = viz_func(motifs, len(my_sequence))
        else:
            fig = viz_func(motifs)
        
        export_plot(fig, f'figure_{i+1}', formats=['png', 'svg'])
        plt.close(fig)
    except Exception as e:
        print(f"Error in visualization {i+1}: {e}")

print("✓ Analysis complete!")
```

## Testing

Run the comprehensive test suite:

```bash
python test_advanced_features.py
```

This will test:
1. Enhanced hybrid detection
2. Enhanced cluster detection
3. Score normalization
4. All 9 advanced visualizations
5. Export functionality

## Performance Notes

- Hybrid detection: O(n²) where n = number of motifs
- Cluster detection: O(n × w) where w = number of windows
- Visualization: 1-5 seconds per plot for typical datasets
- Memory usage: ~50MB for 10,000 motifs

## References

**Scoring Methods:**
- G4Hunter: Bedrat et al. (2016) NAR
- QmRLFS: Jenjaroenpun & Wongsurawat (2016)
- Z-DNA: Ho et al. (1986)
- Curvature: Olson et al. (1998)

**Visualization Best Practices:**
- Okabe & Ito (2008) - Colorblind-safe palettes
- Tufte (2001) - The Visual Display of Quantitative Information
- Wilke (2019) - Fundamentals of Data Visualization

## Support

For issues or questions:
- GitHub: https://github.com/VRYella/NonBScanner
- Email: yvrajesh_bt@kluniversity.in
