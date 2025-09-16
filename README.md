# 🧬 NBDFinder - Non-B DNA Motif Detection Suite

**Comprehensive detection and analysis of Non-B DNA structures with 10 major classes and 22+ specialized subclasses**

## 🎯 Overview

NBDFinder is a state-of-the-art bioinformatics tool for detecting and analyzing Non-B DNA motifs in genomic sequences. It combines high-performance Hyperscan pattern matching with scientific scoring algorithms to provide comprehensive analysis of structural DNA elements.

### 📊 Detection Coverage
- **10 Major Non-B DNA Classes** with comprehensive subclass analysis
- **22+ Specialized Subclasses** for detailed motif characterization
- **Hyperscan-accelerated** detection (40x+ speed improvement)
- **Literature-validated** scoring algorithms
- **Normalized scoring** (0-1 scale) for cross-class comparison

## 🔬 Supported Motif Classes

| Class | Name | Subclasses | Key Features |
|-------|------|------------|--------------|
| **1** | Curved DNA | Global Curvature, Local Curvature | A-tract mediated curvature |
| **2** | Slipped DNA | Direct Repeat, STR | Tandem repeats, slipped structures |
| **3** | Cruciform | Palindromic Inverted Repeat | Four-way junctions |
| **4** | R-Loop | R-Loop Formation Site | RNA-DNA hybrids |
| **5** | Triplex | Mirror Repeat, Sticky DNA | Three-stranded structures |
| **6** | G-Quadruplex | Canonical, Bulged, Relaxed, Bipartite, Multimeric, Imperfect, G-Triplex | Four-stranded G-rich structures |
| **7** | i-Motif | Canonical, Extended, AC-Motif | C-rich structures |
| **8** | Z-DNA | Classic Z-DNA, eGZ | Left-handed double helix |
| **9** | Hybrid | Multi-class Overlap, Composite | Overlapping motifs |
| **10** | Cluster | Motif Hotspot, Mixed Cluster | High-density regions |

## 🚀 Quick Start

### Web Interface
```bash
# Clone and setup
git clone https://github.com/VRYella/NBDFinder.git
cd NBDFinder
pip install -r requirements.txt

# Launch services
./start_services.sh

# Or start individually
streamlit run app.py                    # Web interface on :8501
python api.py                          # REST API on :8000
```

### REST API Usage
```bash
# Health check
curl http://localhost:8000/api/v1/health

# Analyze sequence
curl -X POST http://localhost:8000/api/v1/analyze \
  -H "Content-Type: application/json" \
  -d '{"sequence": "GGGTTAGGGTTAGGGTTAGGG", "sequence_name": "test"}'

# Get comprehensive motif information
curl http://localhost:8000/api/v1/motif-info

# Get API statistics
curl http://localhost:8000/api/v1/stats
```

## 📱 User Interfaces

### 1. **Streamlit Web Application** (`http://localhost:8501`)
- Interactive sequence upload and analysis
- Comprehensive visualization suite (21+ chart types)
- Real-time analysis with progress tracking
- Export capabilities (CSV, BED, BigWig)
- Documentation and tutorial sections

### 2. **REST API** (`http://localhost:8000`)
- Programmatic access for automation
- JSON responses with comprehensive metadata
- Multiple analysis endpoints
- Automatic OpenAPI documentation at `/docs`
- CORS support for web applications

## 🛠️ Technical Features

### Performance
- **Hyperscan Integration**: Intel Hyperscan for ultra-fast pattern matching
- **Parallel Processing**: Multi-core analysis for large datasets
- **Optimized Algorithms**: Literature-validated scoring methods
- **Memory Efficient**: Streaming analysis for large files

### Scientific Accuracy
- **Literature-Based**: Algorithms from peer-reviewed research
- **Validated Thresholds**: Biologically relevant cut-offs
- **Comprehensive Scoring**: Class-specific normalization
- **Quality Control**: Built-in validation and error checking

### Export & Integration
- **Multiple Formats**: BED, BigWig, CSV, JSON export
- **Genome Browser**: UCSC/IGV compatible outputs
- **API Integration**: REST endpoints for pipelines
- **Batch Processing**: Multi-FASTA support

## 📚 API Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| `GET` | `/api/v1/health` | Health check and system status |
| `GET` | `/api/v1/classes` | List all motif classes and subclasses |
| `GET` | `/api/v1/classes/{class_id}` | Get detailed class information |
| `POST` | `/api/v1/analyze` | Comprehensive sequence analysis |
| `POST` | `/api/v1/analyze/{class_id}` | Class-specific analysis |
| `GET` | `/api/v1/motif-info` | Comprehensive motif information |
| `GET` | `/api/v1/stats` | API usage statistics |
| `GET` | `/docs` | Interactive API documentation |

## 🔬 Scoring Algorithms

- **G4Hunter**: G-quadruplex prediction (Bedrat et al., 2016)
- **Z-Seeker**: Z-DNA detection (Ho et al., 1986)
- **Enhanced Curvature**: A-tract curvature (Olson et al., 1998)
- **Instability Scoring**: Repeat instability analysis
- **Thermodynamic Models**: Cruciform stability prediction
- **Triplex Potential**: Three-strand formation scoring

## 📈 Visualization Suite

### Static Plots
- Motif distribution analysis
- Coverage and density maps
- Score distribution histograms
- Sequence composition analysis
- Class/subclass comparisons

### Interactive Visualizations
- Motif browser with zoom/pan
- Class hierarchy sunburst charts
- Position-based track plots
- Statistical correlation plots
- Network analysis graphs

## 🧪 Example Analysis

```python
# Python API usage
import requests

# Analyze a G-quadruplex sequence
response = requests.post('http://localhost:8000/api/v1/analyze', json={
    'sequence': 'GGGTTAGGGTTAGGGTTAGGG',
    'sequence_name': 'G4_example',
    'nonoverlap': False,
    'report_hotspots': True
})

results = response.json()
print(f"Found {results['total_motifs']} motifs")
print(f"Classes detected: {results['classes_detected']}")
print(f"Subclasses detected: {results['subclasses_detected']}")
```

## 📄 Citation

If you use NBDFinder in your research, please cite:

```
NBDFinder: Comprehensive Detection and Analysis of Non-B DNA Motifs
Dr. Venkata Rajesh Yella
GitHub: https://github.com/VRYella/NBDFinder
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
