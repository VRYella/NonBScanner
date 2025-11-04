# NBDScanner Standalone Jupyter Notebook

## 📋 Overview

This Jupyter notebook provides a **complete standalone environment** for running NBDScanner, the comprehensive Non-B DNA motif detection system. It works independently of the Streamlit web interface and can be used for batch analysis, scripting, and integration into bioinformatics pipelines.

## 🚀 Quick Start

### Prerequisites
- Python 3.8 or higher
- Jupyter Notebook or JupyterLab
- Basic Python packages (numpy, pandas, matplotlib)

### Installation

1. **Clone the repository:**
```bash
git clone https://github.com/VRYella/NonBScanner.git
cd NonBScanner
```

2. **Install dependencies:**
```bash
pip install pandas numpy matplotlib seaborn biopython
```

3. **Launch Jupyter:**
```bash
jupyter notebook NBDScanner_Standalone.ipynb
```

## 📚 Features

### Comprehensive Detection
- **11 Major Motif Classes** with 22+ subclasses
- Detects: G-Quadruplex, i-Motif, Z-DNA, Curved DNA, R-loops, Cruciform, Triplex, Slipped DNA, A-philic DNA, Hybrid motifs, and Clusters

### Analysis Capabilities
- Single sequence analysis
- Batch processing of multiple sequences
- FASTA file support
- Customizable detection parameters
- Quality filtering options

### Visualization Suite
- Motif distribution charts
- Coverage maps
- Length distribution analysis
- Score distribution plots
- Nested pie charts (class-subclass hierarchy)

### Export Options
- **CSV** - For spreadsheet analysis
- **BED** - For genome browser visualization (UCSC, IGV)
- **JSON** - For programmatic access

## 🔬 Usage Examples

### Example 1: Analyze a Simple Sequence
```python
from scanner import analyze_sequence
from visualizations import plot_motif_distribution

# Define your sequence
sequence = "GGGTTAGGGTTAGGGTTAGGG..."

# Run analysis
motifs = analyze_sequence(sequence, "my_sequence")

# Visualize results
plot_motif_distribution(motifs, by='Class')
```

### Example 2: Batch Analysis from FASTA
```python
from utilities import parse_fasta
from scanner import analyze_sequence

# Load FASTA file
with open("sequences.fasta", "r") as f:
    fasta_content = f.read()
sequences = parse_fasta(fasta_content)

# Analyze all sequences
results = {}
for name, seq in sequences.items():
    motifs = analyze_sequence(seq, name)
    results[name] = motifs
```

### Example 3: Export Results
```python
from utilities import export_to_csv, export_to_bed, export_to_json

# Export to different formats
csv_data = export_to_csv(motifs)
bed_data = export_to_bed(motifs, "sequence_name")
json_data = export_to_json(motifs, pretty=True)

# Save to files
with open("results.csv", "w") as f:
    f.write(csv_data)
```

## 📊 Notebook Structure

The notebook is organized into clear sections:

1. **Introduction** - Overview and motif class table
2. **Installation & Setup** - Dependency installation and module imports
3. **Quick Start** - Simple examples to get started
4. **Motif Detection** - Detailed analysis workflows
5. **Visualization** - Comprehensive plotting options
6. **Export Results** - Data export in multiple formats
7. **Advanced Usage** - Batch processing and custom parameters

## 🔧 Technical Details

### Performance
- **Detection Speed**: 5,000-280,000 bp/second (varies by detector)
- **Memory Efficient**: ~5 MB for 100K sequences
- **Scalable**: Linear O(n) complexity for all detectors

### Architecture
The notebook uses the same detection engine as the web interface:
- **Modular detectors**: 9 specialized detector classes
- **Pattern-based detection**: Regex and Hyperscan acceleration
- **Algorithmic scoring**: Literature-based scoring methods
- **Overlap resolution**: Intelligent handling of overlapping motifs

### Supported Motif Classes

| # | Class | Subclasses | Detection Method |
|---|-------|-----------|-----------------|
| 1 | Curved DNA | Global, Local | A-tract phasing analysis |
| 2 | Slipped DNA | Direct Repeat, STR | K-mer index seed-and-extend |
| 3 | Cruciform | Inverted Repeats | Palindrome detection |
| 4 | R-Loop | Formation sites | GC-skew and QmRLFS |
| 5 | Triplex | Mirror Repeat, Sticky | Purine/pyrimidine content |
| 6 | G-Quadruplex | 7 subclasses | G4Hunter scoring |
| 7 | i-Motif | Canonical, Relaxed, AC | C-tract detection |
| 8 | Z-DNA | Classic, eGZ | 10-mer scoring table |
| 9 | A-philic | A-philic DNA | Tetranucleotide scoring |
| 10 | Hybrid | Overlaps | Multi-class intersection |
| 11 | Cluster | Hotspots | High-density detection |

## 📖 Documentation

### In-Notebook Documentation
The notebook includes:
- Inline code comments
- Markdown explanations
- Example outputs
- Usage tips and best practices

### Additional Resources
- **Main README**: [README.md](README.md)
- **API Documentation**: See function docstrings in source files
- **Web Interface**: Use `streamlit run app.py` for GUI

## 🧪 Testing

The notebook includes test cells to verify:
- Module imports
- Detector functionality
- Visualization generation
- Export capabilities

Run all cells to perform a complete system test.

## 🤝 Contributing

This is a read-only standalone notebook. For contributions to NBDScanner:
1. Fork the repository
2. Make changes to source files
3. Test with both notebook and web interface
4. Submit a pull request

## 📞 Contact

**Dr. Venkata Rajesh Yella**
- Email: yvrajesh_bt@kluniversity.in
- GitHub: [@VRYella](https://github.com/VRYella)
- Repository: [NonBScanner](https://github.com/VRYella/NonBScanner)

## 📜 License

MIT License - See [LICENSE](LICENSE) file for details.

## 🔖 Citation

If you use NBDScanner in your research, please cite:

```
NonBScanner: Comprehensive Detection and Analysis of Non-B DNA Motifs
Dr. Venkata Rajesh Yella
GitHub: https://github.com/VRYella/NonBScanner
Version: 2024.1
```

## 🎯 Tips for Best Results

1. **Sequence Quality**: Use only validated DNA sequences (ATGC)
2. **Length Considerations**: Works best on sequences 100 bp - 1 Mbp
3. **Score Thresholds**: Higher scores indicate higher confidence
4. **Hybrid/Cluster Analysis**: Analyze separately for complex regions
5. **Export Formats**: Choose based on downstream analysis needs

## 🐛 Troubleshooting

### Common Issues

**Problem**: Module import errors
**Solution**: Ensure you're running the notebook from the NonBScanner directory

**Problem**: No motifs detected
**Solution**: Check sequence validity and try longer sequences

**Problem**: Visualization not showing
**Solution**: Ensure matplotlib backend is properly configured

### Getting Help

1. Check the in-notebook documentation
2. Review the main README.md
3. Open an issue on GitHub
4. Contact the author directly

---

**Last Updated**: November 2024  
**Version**: 2024.1  
**Status**: Production Ready ✅
