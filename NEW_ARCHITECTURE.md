# NonBScanner - Consolidated Architecture

## 📁 New File Organization (5 Core Files)

NonBScanner has been reorganized from 29 Python files into **5 well-organized core files** for improved maintainability while preserving all functionality.

### Core Modules

#### 1. `app.py` (1,810 lines)
**Main Streamlit Web Application**
- Interactive UI for sequence analysis
- Real-time progress tracking
- Comprehensive visualization interface
- Export functionality (CSV, BED, JSON)
- Multi-FASTA support
- NCBI sequence fetching

**Usage:**
```bash
streamlit run app.py
```

#### 2. `detectors.py` (3,242 lines)
**All Motif Detection Classes**

Contains 10 detector classes for Non-B DNA structure detection:
- `BaseMotifDetector` - Abstract base class
- `CurvedDNADetector` - A-tract mediated curvature
- `ZDNADetector` - Left-handed helix and Z-DNA
- `APhilicDetector` - A-philic DNA sequences
- `SlippedDNADetector` - Direct repeats and STRs
- `CruciformDetector` - Palindromic inverted repeats
- `RLoopDetector` - R-loop formation sites
- `TriplexDetector` - Triplex and mirror repeats
- `GQuadruplexDetector` - G4 and variants
- `IMotifDetector` - i-Motif and AC-motif

**Usage:**
```python
from detectors import ZDNADetector, GQuadruplexDetector

# Create detector instance
detector = ZDNADetector()

# Detect motifs
motifs = detector.detect(sequence)
```

**Performance:**
- Hyperscan-based: 15,000-65,000 bp/s
- Algorithmic: 5,000-280,000 bp/s
- All detectors: O(n) linear complexity

#### 3. `scanner.py` (1,600 lines)
**Scanner Orchestration and Analysis**

Main analysis functions:
- `analyze_sequence()` - Single sequence analysis
- `analyze_multiple_sequences()` - Batch processing
- `detect_hybrids_and_clusters()` - Advanced overlap detection
- `export_results_to_dataframe()` - Export to pandas
- `get_motif_classification_info()` - Metadata retrieval

Also includes:
- Repeat detection algorithms
- Overlap resolution
- Hybrid and cluster detection
- Quality filtering

**Usage:**
```python
from scanner import analyze_sequence

# Analyze a DNA sequence
results = analyze_sequence(sequence, sequence_name="chr1")

# Results contain detected motifs with:
# - Class, Subclass, Start, End, Length
# - Score, GC content, Sequence
# - Pattern information
```

#### 4. `utilities.py` (2,076 lines)
**Utility Functions and Data Processing**

**Sequence Processing:**
- `parse_fasta()` - Parse FASTA format
- `gc_content()` - Calculate GC content
- `reverse_complement()` - Generate reverse complement
- `validate_sequence()` - Validate DNA sequence

**Pattern Loading:**
- `load_hyperscan_db()` - Load Hyperscan databases
- `load_regex_registry()` - Load pattern registries
- `get_motif_patterns()` - Get motif patterns

**Data Export:**
- `export_to_csv()` - CSV format export
- `export_to_bed()` - BED format for genome browsers
- `export_to_json()` - JSON format with metadata

**Statistics:**
- `get_basic_stats()` - Calculate sequence statistics
- `quality_check_motifs()` - Validate motif quality

**Usage:**
```python
from utilities import parse_fasta, export_to_bed, gc_content

# Parse FASTA file
sequences = parse_fasta(fasta_content)

# Calculate GC content
gc = gc_content(sequence)

# Export to BED format
bed_data = export_to_bed(motifs, sequence_name)
```

#### 5. `visualizations.py` (967 lines)
**Visualization and Plotting Functions**

**21+ Visualization Types:**
- Motif distribution charts
- Coverage maps and density plots
- Length distribution analysis
- Score distribution histograms
- Interactive plots (Plotly)
- Publication-quality figures

**Features:**
- Colorblind-friendly palettes (Okabe-Ito)
- SVG and PNG export (@300 DPI)
- Interactive zoom/pan capabilities
- Annotated peaks and features

**Usage:**
```python
from visualizations import (
    plot_motif_distribution, 
    plot_coverage_map,
    plot_length_distribution,
    MOTIF_CLASS_COLORS
)

# Create distribution plot
fig = plot_motif_distribution(motifs, by='Class')

# Create coverage map
fig = plot_coverage_map(motifs, sequence_length)
```

### Supporting Files

- **`registry/`** - Pattern registry files (9 JSON files)
- **`tests/`** - Test suite (6 test files)
- **`tools/`** - Utility scripts (pattern generation)
- **Documentation** - README, guides, architecture docs

## 🔄 Migration Guide

### Old Import Pattern → New Import Pattern

**Detectors:**
```python
# OLD
from motif_detection.z_dna_detector import ZDNADetector
from motif_detection.g_quadruplex_detector import GQuadruplexDetector

# NEW
from detectors import ZDNADetector, GQuadruplexDetector
```

**Scanner:**
```python
# OLD
from utils.nbdscanner import analyze_sequence
from utils.modular_scanner import ModularMotifDetector

# NEW
from scanner import analyze_sequence, ModularMotifDetector
```

**Utilities:**
```python
# OLD
from utils.utils import parse_fasta, gc_content
from utils.motif_patterns import get_motif_patterns

# NEW
from utilities import parse_fasta, gc_content, get_motif_patterns
```

**Visualizations:**
```python
# OLD
from utils.visualization import plot_motif_distribution

# NEW
from visualizations import plot_motif_distribution
```

## ✅ Benefits of Consolidation

### Code Organization
- **29 files → 5 files**: 80% reduction in file count
- Clear separation of concerns
- Easier navigation and understanding
- Reduced import complexity

### Maintainability
- Related functionality grouped together
- Single location for each major concern
- Easier to find and modify code
- Reduced cognitive load

### Development
- Faster imports (fewer files to parse)
- Easier debugging (less jumping between files)
- Better IDE navigation
- Simpler project structure

### Testing
- Clear module boundaries
- Easier to mock dependencies
- More focused unit tests
- Better integration testing

## 🧪 Testing

### Quick Test
```bash
# Test basic imports and functionality
python test_consolidated.py
```

### Full Test Suite
```bash
# Run all unit tests
python -m pytest tests/ -v
```

### Manual Testing
```bash
# Test web application
streamlit run app.py

# Test basic detector
python -c "from detectors import ZDNADetector; d = ZDNADetector(); print(d)"

# Test scanner
python -c "from scanner import analyze_sequence; print('Scanner OK')"
```

## 📊 File Statistics

| File | Lines | Size | Purpose |
|------|-------|------|---------|
| app.py | 1,810 | 90 KB | Web application |
| detectors.py | 3,242 | 145 KB | Motif detection |
| scanner.py | 1,600 | 67 KB | Analysis orchestration |
| utilities.py | 2,076 | 79 KB | Utilities & I/O |
| visualizations.py | 967 | 35 KB | Plotting functions |
| **Total** | **9,695** | **416 KB** | **5 core files** |

## 🔧 Architecture Preserved

✅ **All original functionality maintained:**
- All 10 detector classes intact
- All 11 motif classes supported
- 22+ subclass analysis
- Hyperscan integration preserved
- Performance characteristics unchanged
- Pattern registries compatible
- Export formats supported

✅ **No performance impact:**
- Same algorithms
- Same data structures
- Same optimization techniques
- Linear O(n) complexity maintained

## 📝 Development Workflow

### Adding a New Detector
1. Add class to `detectors.py`
2. Inherit from `BaseMotifDetector`
3. Implement required methods
4. Register in `scanner.py`
5. Add tests

### Modifying Utilities
1. Edit function in `utilities.py`
2. Update docstring
3. Add/update tests
4. Verify imports still work

### Adding Visualizations
1. Add function to `visualizations.py`
2. Follow existing patterns
3. Use MOTIF_CLASS_COLORS
4. Add to __all__ if public

## 🚀 Quick Start

```bash
# Clone repository
git clone https://github.com/VRYella/NonBScanner.git
cd NonBScanner

# Install dependencies
pip install -r requirements.txt

# Test installation
python test_consolidated.py

# Run web application
streamlit run app.py
```

## 📚 Documentation

- **README.md** - Main project documentation
- **CONSOLIDATION_SUMMARY.md** - Detailed consolidation information
- **ORGANIZATION.md** - Original organization documentation
- **QUICK_START.md** - Quick start guide
- **HYPERSCAN_ARCHITECTURE.md** - Hyperscan details
- **OPTIMIZED_SCANNER_ARCHITECTURE.md** - Performance details

## 🤝 Contributing

When contributing to NonBScanner:
1. Follow the 5-file structure
2. Keep related functionality together
3. Update imports in tests
4. Add tests for new features
5. Update documentation

## 📞 Support

**Dr. Venkata Rajesh Yella**
- Email: yvrajesh_bt@kluniversity.in
- GitHub: [@VRYella](https://github.com/VRYella)

---

*Consolidated architecture implemented: 2024*
*All functionality preserved • Performance maintained • Easier to maintain*
