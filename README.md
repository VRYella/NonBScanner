# 🧬 NonBScanner - Non-B DNA Motif Detection Suite

**Comprehensive detection and analysis of Non-B DNA structures with 10 major classes and 22+ specialized subclasses**

## 🎯 Overview

NonBScanner is a state-of-the-art bioinformatics tool for detecting and analyzing Non-B DNA motifs in genomic sequences. It combines high-performance optimized algorithms with scientific scoring methods to provide comprehensive analysis of structural DNA elements.

### 📊 Detection Coverage
- **10 Major Non-B DNA Classes** with comprehensive subclass analysis
- **22+ Specialized Subclasses** for detailed motif characterization
- **High-performance detection** (24,674 bp/s on 100K sequences)
- **Literature-validated** scoring algorithms
- **Normalized scoring** (0-1 scale) for cross-class comparison
- **✨ Enhanced hybrid/cluster detection** with actual sequence extraction
- **✨ 21+ advanced publication-quality visualizations** (colorblind-friendly)

### ⚡ Performance Highlights
- **100,000 bp in 4 seconds** with optimized detectors
- **24,674 bp/second** processing speed
- **Memory efficient**: ~5 MB for 100K sequences
- **Production ready**: Tested on large genomic datasets

## 🔬 Supported Motif Classes

| Class | Name | Subclasses | Key Features |
|-------|------|------------|--------------|
| **1** | Curved DNA | Global Curvature, Local Curvature | A-tract mediated curvature |
| **2** | Slipped DNA | Direct Repeat, STR | Tandem repeats, slipped structures |
| **3** | Cruciform | Palindromic Inverted Repeat | Four-way junctions |
| **4** | R-Loop | R-Loop Formation Site, QmRLFS-m1, QmRLFS-m2 | RNA-DNA hybrids, QmRLFS algorithm |
| **5** | Triplex | Mirror Repeat, Sticky DNA | Three-stranded structures |
| **6** | G-Quadruplex | Canonical, Bulged, Relaxed, Bipartite, Multimeric, Imperfect, G-Triplex | Four-stranded G-rich structures |
| **7** | i-Motif | Canonical, Extended, AC-Motif | C-rich structures |
| **8** | Z-DNA | Classic Z-DNA, eGZ | Left-handed double helix |
| **9** | Hybrid | Multi-class Overlap, Composite | Overlapping motifs |
| **10** | Cluster | Motif Hotspot, Mixed Cluster | High-density regions |

## 📚 Documentation

Comprehensive documentation is available:

- **[QUICK_START.md](QUICK_START.md)** - Fast installation and first analysis guide
- **[OPTIMIZED_SCANNER_ARCHITECTURE.md](OPTIMIZED_SCANNER_ARCHITECTURE.md)** - NEW! Optimized repeat scanner architecture
- **[HYPERSCAN_ARCHITECTURE.md](HYPERSCAN_ARCHITECTURE.md)** - Hyperscan integration details
- **[TOOL_DOCUMENTATION.md](TOOL_DOCUMENTATION.md)** - Complete technical documentation (nature paper-level writeup)
- **[VISUAL_FLOWCHARTS.md](VISUAL_FLOWCHARTS.md)** - 20+ interactive flowcharts and diagrams
- **[CODE_ORGANIZATION_SUMMARY.md](CODE_ORGANIZATION_SUMMARY.md)** - Code structure and organization
- **[OPTIMIZATION_SUMMARY.md](OPTIMIZATION_SUMMARY.md)** - Performance optimization details
- **[PERFORMANCE_OPTIMIZATION.md](PERFORMANCE_OPTIMIZATION.md)** - Detailed performance guide

## 🚀 Quick Start

### Web Interface
```bash
# Clone and setup
git clone https://github.com/VRYella/NonBScanner.git
cd NonBScanner
pip install -r requirements.txt

# Launch web interface
streamlit run app.py                    # Web interface on :8501
```

**New to NonBScanner?** Start with [QUICK_START.md](QUICK_START.md) for a step-by-step guide!

## 📱 User Interfaces

### **Streamlit Web Application** (`http://localhost:8501`)
- Interactive sequence upload and analysis
- Comprehensive visualization suite (21+ chart types)
- Real-time analysis with progress tracking
- Export capabilities (CSV, BED, BigWig)
- Documentation and tutorial sections

## 🛠️ Technical Features

### Performance
- **Optimized Algorithms**: Linear O(n) complexity for all major detectors
- **Seed-and-Extend K-mer Index**: Genome-scale efficient repeat detection
- **No Sequence Limits**: Handle sequences of any size without performance degradation
- **~280,000 bp/second**: Slipped DNA detector on 50kb sequences
- **~5,800 bp/second**: Overall rate including all detector types on 10kb sequences
- **Memory Efficient**: K-mer indexing with frequency limits for safety
- **Production Ready**: Rigorously tested on real genome sequences

### Architecture Highlights
- **Hybrid Detection Approach**:
  - **Optimized Python Scanner** (seed-and-extend k-mer index) for Slipped DNA, Cruciform, Triplex
  - **Hyperscan** (regex/pattern matching) for Z-DNA, G4, i-Motif, Curved DNA, A-Philic
  - **Algorithmic** (QmRLFS, GC-skew) for R-Loop detection
- See `OPTIMIZED_SCANNER_ARCHITECTURE.md` for detailed architecture documentation

### Performance Notes
- **Slipped DNA**: No size limits, O(n) complexity (was O(n²) with 50kb limit)
- **Cruciform**: No size limits, O(n) complexity (was O(n²) with 1kb limit)
- **Triplex**: No size limits, O(n) complexity with purine/pyrimidine filtering
- **All detectors**: Linear scaling validated on sequences up to 50kb+
- See `OPTIMIZED_SCANNER_ARCHITECTURE.md` for benchmarks

### Scientific Accuracy
- **Literature-Based**: Algorithms from peer-reviewed research
- **QmRLFS Integration**: Advanced R-loop detection with RIZ/REZ analysis
- **Validated Thresholds**: Biologically relevant cut-offs
- **Raw Score Reporting**: Direct algorithm outputs without normalization
- **Overlap Resolution**: Automatic removal of overlapping motifs within subclasses
- **Quality Control**: Built-in validation and error checking

### Export & Integration
- **Multiple Formats**: BED, BigWig, CSV, JSON export
- **Genome Browser**: UCSC/IGV compatible outputs
- **Batch Processing**: Multi-FASTA support

## 🔬 Scoring Algorithms

- **G4Hunter**: G-quadruplex prediction (Bedrat et al., 2016)
- **QmRLFS**: R-loop formation site detection with RIZ/REZ analysis (Jenjaroenpun & Wongsurawat, 2016)
- **Z-Seeker**: Z-DNA detection (Ho et al., 1986)
- **Enhanced Curvature**: A-tract curvature (Olson et al., 1998)
- **Instability Scoring**: Repeat instability analysis
- **Thermodynamic Models**: Cruciform stability prediction
- **Triplex Potential**: Three-strand formation scoring

## 📈 Visualization Suite

### Static Plots (Classic)
- Motif distribution analysis
- Coverage and density maps
- Length distribution analysis
- Sequence composition analysis
- Class/subclass comparisons

### Advanced Visualizations (NEW! 🎨)
**Publication-quality static plots with colorblind-friendly palettes:**

1. **Genome Landscape Track** - Horizontal ruler with colored glyphs showing motif positions
2. **Sliding Window Heat Ribbon** - 1D heatmap with density and score overlay
3. **Ridge Plots (Joyplots)** - Stacked density ridges for length distributions
4. **Sunburst/Treemap** - Hierarchical composition visualization
5. **Hexbin with Marginals** - 2D density plot with marginal histograms
6. **UpSet Plot** - Clear intersection visualization (better than Venn diagrams)
7. **Violin + Beeswarm** - Score distributions with individual data points
8. **Cluster Hotspot Map** - Regional cluster analysis with annotations

**Design Features:**
- ✅ Export as SVG (vector, editable) and PNG @300 DPI
- ✅ Okabe-Ito colorblind-safe palette (7 colors)
- ✅ Clean sans-serif typography
- ✅ Annotated peaks and top features
- ✅ Publication-ready layouts

### Interactive Visualizations
- Motif browser with zoom/pan
- Class hierarchy sunburst charts
- Position-based track plots
- Statistical correlation plots
- Network analysis graphs

## 🧪 Example Analysis

Use the web interface at `http://localhost:8501` to:
- Upload FASTA sequences or paste sequence data
- Analyze G-quadruplex, Z-DNA, R-loops, and other motifs
- Visualize results with comprehensive charts
- Export findings in BED, CSV, or JSON formats

## 🔗 Hybrid and Cluster Motif Separation (NEW)

NBDScanner now separates hybrid and cluster motifs from regular Non-B DNA motifs for cleaner analysis:

### What are Hybrid and Cluster Motifs?

- **Hybrid Motifs**: Regions where different Non-B DNA classes overlap (30-70% overlap)
  - Example: `R-Loop_Cruciform_Overlap` - indicates complex genomic regions
  - Shows interaction between different structural elements

- **Cluster Motifs**: High-density regions containing multiple Non-B DNA motifs from different classes
  - Example: `Mixed_Cluster_10_classes` - hotspots of Non-B DNA activity
  - Indicates regions with exceptional structural diversity

### How It Works

1. **Main Results**: Shows only primary Non-B DNA motifs (e.g., 52 motifs)
2. **Cluster/Hybrid Tab**: Dedicated visualization tab for hybrid/cluster motifs (e.g., 71 motifs)
3. **Downloads**: Export files contain only regular motifs for cleaner downstream analysis
4. **Clear Messaging**: Info indicators guide users to the appropriate location

### Benefits

- ✅ **Cleaner Results**: Main results focus on primary motifs
- ✅ **Focused Analysis**: Analyze regular motifs without complex overlaps
- ✅ **Advanced Access**: Hybrid/cluster data still available in dedicated tab
- ✅ **Better Downloads**: Export files contain clean datasets

## 📄 Citation

If you use NonBScanner in your research, please cite:

```
NonBScanner: Comprehensive Detection and Analysis of Non-B DNA Motifs
Dr. Venkata Rajesh Yella
GitHub: https://github.com/VRYella/NonBScanner
```

## 🤝 Contributing

We welcome contributions! Please see our contributing guidelines and submit pull requests for:
- New motif detection algorithms
- Additional visualization features
- Performance improvements
- Documentation enhancements

## 📞 Contact

**Dr. Venkata Rajesh Yella**
- Email: yvrajesh_bt@kluniversity.in
- GitHub: [@VRYella](https://github.com/VRYella)

## 📜 License

This project is licensed under the MIT License - see the LICENSE file for details.

---

*Developed for the scientific community to advance our understanding of Non-B DNA structures and their biological roles.*
