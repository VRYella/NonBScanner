# NBDScanner Hyperscan Revamp - Architecture Documentation

## Overview

This document describes the revamped architecture of NBDScanner with integrated Hyperscan support for high-performance Non-B DNA motif detection.

## New Architecture Files

### Core Detection Modules

#### 1. `hyperscan_registry.py`
**Purpose**: Comprehensive registry of all Hyperscan-compatible regex patterns

**Contents**:
- 66 optimized regex patterns across 8 motif classes
- Pattern metadata (IDs, names, subclasses, scoring methods)
- Hyperscan database compilation
- Color schemes for visualization

**Hyperscan-Compatible Classes**:
1. Curved DNA (6 patterns)
2. Slipped DNA/STR (24 patterns - specific repeats)
3. R-loop (5 patterns)
4. Triplex (8 patterns)
5. G-Quadruplex (9 patterns)
6. i-Motif Family (5 patterns)
7. Z-DNA (11 patterns)

**Key Features**:
- No backreferences (Hyperscan limitation)
- Optimized pattern ordering
- Pattern ID mapping to classes
- Compile-once, use-many architecture

#### 2. `hyperscan_detection.py`
**Purpose**: High-performance pattern matching using Hyperscan engine

**Contents**:
- `HyperscanDetector`: Main detection class
- `HyperscanFallback`: Pure Python regex fallback
- `create_detector()`: Factory function

**Performance**:
- ~50,000 bp/s with Hyperscan
- Automatic fallback to regex if Hyperscan unavailable
- Batch processing support
- Performance statistics tracking

**Key Features**:
- Single database scan for all patterns
- Memory-efficient callback system
- Graceful error handling
- Pattern-to-class mapping

#### 3. `hyperscan_scoring_overlap_resolver.py`
**Purpose**: Scientific scoring and intelligent overlap resolution

**Contents**:
- `MotifScorer`: 12 literature-validated scoring algorithms
- `OverlapResolver`: Multiple overlap resolution strategies
- `score_and_resolve_motifs()`: Integrated function

**Scoring Algorithms**:
1. G4Hunter (G-quadruplex)
2. Curvature scoring (Curved DNA)
3. Phasing score (Global curvature)
4. Z-DNA alternating pattern score
5. eGZ (Extruded-G) score
6. Triplex potential score
7. Mirror repeat score
8. R-loop formation potential
9. i-Motif score
10. AC-motif score
11. Repeat instability score
12. G-Triplex score

**Overlap Strategies**:
- Remove overlaps within subclass (default)
- Keep highest score
- Keep longest
- Keep all
- Merge overlapping

**Normalization**: All scores normalized to [0, 1] range

#### 4. `non_hyperscan_detection.py`
**Purpose**: Algorithmic detection for complex motifs not feasible with regex

**Contents**:
- `CruciformDetector`: Inverted repeat finder
- `APhilicDetector`: Tetranucleotide scoring
- `DirectRepeatDetector`: Repeat finding with variable spacers
- `HybridMotifDetector`: Cross-class overlap analysis
- `ClusterDetector`: Density-based clustering

**Non-Hyperscan Classes**:
1. **Cruciform** (Class 3): Requires reverse complement computation
2. **A-philic DNA** (Class 9): Tetranucleotide log2-odds scoring
3. **Direct Repeats** (Class 2 subtype): Variable spacer repeats
4. **Hybrid Motifs** (Class 10): 30-70% overlap between classes
5. **Non-B DNA Clusters** (Class 11): Sliding window density analysis

**Performance Limits**:
- Cruciform: Limited to 10,000 bp sequences (O(n²) complexity)
- Direct Repeats: Skips sequences >50,000 bp
- Other methods: Linear time complexity

#### 5. `visualization.py`
**Purpose**: Comprehensive plotting and visualization suite

**Contents**:
- Static matplotlib plots
- Interactive Plotly visualizations
- Publication-quality styling
- Colorblind-friendly palettes

**Visualization Categories**:
1. Distribution plots (class/subclass)
2. Coverage maps (sequence-wide)
3. Score distributions
4. Length distributions
5. Nested pie charts
6. Interactive browsers

**Features**:
- Motif distribution over sequence
- Density heatmaps
- Multi-class comparisons
- Export to SVG/PNG

#### 6. `app.py`
**Purpose**: Streamlit web application interface

**Architecture**:
```
Input → Detection → Scoring → Overlap Resolution → Visualization → Export
```

**Features**:
- Multi-FASTA support
- Real-time progress tracking
- Interactive visualizations
- Export to CSV/BED/JSON
- NCBI sequence fetch
- Example datasets

**Tabs**:
1. Home: Overview and introduction
2. Upload & Analyze: Sequence input and analysis
3. Results: Visualization and statistics
4. Download: Export options
5. Documentation: Scientific references

## Integration Architecture

### Detection Pipeline

```
┌─────────────────────────────────────────────────────────────┐
│                    SEQUENCE INPUT                            │
│  (FASTA, paste, NCBI, examples)                             │
└────────────────────┬────────────────────────────────────────┘
                     │
        ┌────────────┴────────────┐
        │                         │
        ▼                         ▼
┌───────────────┐         ┌──────────────────┐
│   HYPERSCAN   │         │   ALGORITHMIC    │
│   DETECTION   │         │    DETECTION     │
│  (8 classes)  │         │  (4 classes)     │
└───────┬───────┘         └────────┬─────────┘
        │                          │
        │  ┌───────────────────┐   │
        └─►│    SCORING &      │◄──┘
           │ OVERLAP RESOLVER  │
           └─────────┬─────────┘
                     │
        ┌────────────┴────────────┐
        │                         │
        ▼                         ▼
┌───────────────┐         ┌──────────────────┐
│ VISUALIZATION │         │     EXPORT       │
│   & ANALYSIS  │         │ (CSV/BED/JSON)   │
└───────────────┘         └──────────────────┘
```

### Class Distribution

**Total**: 11 Major Classes, 22+ Subclasses

**Hyperscan-Compatible** (8/11 classes, ~73%):
- Curved DNA
- Slipped DNA (STR subclass)
- R-loop
- Triplex
- G-Quadruplex (7 subclasses)
- i-Motif Family (3 subclasses)
- Z-DNA (2 subclasses)

**Algorithmic** (3/11 classes, ~27%):
- Cruciform
- A-philic DNA
- Slipped DNA (Direct Repeat subclass)

**Post-Processing** (2/11 classes):
- Hybrid (from overlaps)
- Non-B DNA Clusters (from density)

## Performance Characteristics

### Speed
- **Hyperscan**: ~50,000 bp/s
- **Regex Fallback**: ~10,000 bp/s
- **Algorithmic**: Varies by method
  - A-philic: ~20,000 bp/s
  - Cruciform: ~500 bp/s (limited to 10kb)
  - Direct Repeats: ~1,000 bp/s (limited to 50kb)

### Memory
- Hyperscan database: ~5 MB (one-time compilation)
- Per-sequence overhead: Minimal
- Scalable to large genomes with streaming

### Accuracy
- Literature-validated patterns
- Scientific scoring algorithms
- Quality thresholds
- Overlap resolution

## File Organization

### Required Files (Root Directory)
```
.
├── app.py                                  # Streamlit web application
├── hyperscan_registry.py                   # Pattern registry (66 patterns)
├── hyperscan_detection.py                  # Hyperscan detection engine
├── hyperscan_scoring_overlap_resolver.py   # Scoring and overlap resolution
├── non_hyperscan_detection.py              # Algorithmic detection
├── visualization.py                         # Plotting suite
├── nbdcircle.JPG                           # Logo image
├── requirements.txt                        # Python dependencies
└── README.md                               # Project documentation
```

### Legacy Files (Kept for Compatibility)
```
utils/                    # Original utility modules (used by app.py)
├── __init__.py
├── utils.py              # Helper functions
├── nbdscanner.py         # Original detection (not used by new modules)
└── visualization.py      # Original plots (copied to root)
```

### Documentation Files
- README.md
- HYPERSCAN_REVAMP_README.md (this file)
- Various markdown documentation files

## Usage Examples

### Basic Detection
```python
from hyperscan_detection import create_detector
from hyperscan_scoring_overlap_resolver import score_and_resolve_motifs

# Create detector
detector = create_detector()

# Detect motifs
sequence = "AAAAAAGGGGTTAGGGTTAGGGTTAGGG"
raw_motifs = detector.detect_motifs(sequence, "MySequence")

# Score and resolve
motifs = score_and_resolve_motifs(raw_motifs, 'remove_within_subclass')

print(f"Found {len(motifs)} motifs")
```

### Full Analysis
```python
from hyperscan_detection import create_detector
from hyperscan_scoring_overlap_resolver import score_and_resolve_motifs
from non_hyperscan_detection import detect_non_hyperscan_motifs

def analyze_sequence(sequence, name="Sequence"):
    # Hyperscan detection
    detector = create_detector()
    hyper_motifs = detector.detect_motifs(sequence, name)
    scored_hyper = score_and_resolve_motifs(hyper_motifs)
    
    # Algorithmic detection
    algo_motifs = detect_non_hyperscan_motifs(sequence, name, scored_hyper)
    scored_algo = score_and_resolve_motifs(algo_motifs, 'keep_all')
    
    # Combine
    all_motifs = scored_hyper + scored_algo
    return sorted(all_motifs, key=lambda m: m['Start'])
```

### Visualization
```python
from visualization import plot_motif_distribution, plot_coverage_map

# Create plots
fig1 = plot_motif_distribution(motifs, by='Class')
fig2 = plot_coverage_map(motifs, len(sequence))

# Display or save
import matplotlib.pyplot as plt
plt.show()
```

## Testing

### Unit Tests
```bash
# Test individual modules
python3 hyperscan_registry.py        # Test pattern registry
python3 hyperscan_detection.py       # Test detection
python3 non_hyperscan_detection.py   # Test algorithmic detection
```

### Integration Test
```python
# Test full pipeline
from hyperscan_detection import create_detector
from hyperscan_scoring_overlap_resolver import score_and_resolve_motifs

test_seq = "AAAAAAGGGGTTAGGGTTAGGGTTAGGGCCCCTTCCCCTTCCCCTTCCCC"
detector = create_detector()
motifs = detector.detect_motifs(test_seq, "Test")
scored = score_and_resolve_motifs(motifs)
print(f"✓ Found {len(scored)} motifs")
```

### Web Application
```bash
streamlit run app.py
```

## Migration from Old Architecture

### What Changed
1. **Detection split** into Hyperscan and algorithmic methods
2. **Patterns optimized** for Hyperscan (no backreferences)
3. **Scoring centralized** in dedicated module
4. **Overlap resolution** improved with multiple strategies
5. **Visualization** copied to root with tabular header

### What Stayed
- app.py UI structure (with new header)
- utils/ folder (for backward compatibility)
- nbdcircle.JPG
- Configuration files
- Overall detection logic

### Breaking Changes
- New modules don't use utils/nbdscanner.py
- Pattern IDs reassigned
- Some subclass names standardized

## Dependencies

### Required
```
streamlit >= 1.28.0
numpy >= 1.21.0
pandas >= 1.3.0
matplotlib >= 3.5.0
seaborn >= 0.11.0
hyperscan >= 0.7.0  # For high-performance detection
```

### Optional
```
biopython >= 1.79   # For NCBI fetch
plotly >= 5.17.0    # For interactive plots
```

## Performance Tuning

### Hyperscan Optimization
- Compile database once at initialization
- Reuse database for multiple sequences
- Use appropriate pattern flags

### Memory Management
- Stream large files instead of loading all at once
- Clear results after processing
- Limit concurrent sequences

### Speed Optimization
- Use Hyperscan when available
- Skip algorithmic methods for large sequences when appropriate
- Parallelize multi-sequence analysis (future enhancement)

## Known Limitations

### Hyperscan
- No backreferences support
- Limited pattern complexity
- Platform-specific (x86_64 primarily)

### Algorithmic Methods
- Cruciform: O(n²) complexity, limited to 10kb
- Direct Repeats: Limited to 50kb sequences
- A-philic: Requires tetranucleotide table

### General
- Memory usage scales with sequence length
- Large pattern databases may slow compilation
- Some motifs may be missed if patterns incomplete

## Future Enhancements

### Planned
- [ ] Parallel sequence processing
- [ ] GPU acceleration for scoring
- [ ] Additional scoring algorithms
- [ ] More comprehensive pattern library
- [ ] REST API for programmatic access
- [ ] Cloud deployment support

### Under Consideration
- [ ] Real-time streaming analysis
- [ ] Integration with genome browsers
- [ ] Machine learning-based scoring
- [ ] Custom pattern upload
- [ ] Batch processing CLI tool

## References

### Hyperscan
- Intel Hyperscan: https://www.hyperscan.io/
- Python bindings: https://github.com/darvid/python-hyperscan

### Motif Detection
- Bedrat et al., 2016 (G4Hunter)
- Jenjaroenpun & Wongsurawat, 2016 (QmRLFS)
- Ho et al., 1986 (Z-DNA)
- Olson et al., 1998 (DNA curvature)
- Vinogradov, 2003 (A-philic tetranucleotides)

## Support

For issues, questions, or contributions:
- GitHub Issues: https://github.com/VRYella/NonBScanner/issues
- Email: yvrajesh_bt@kluniversity.in

## License

MIT License - See LICENSE file for details

## Authors

- Dr. Venkata Rajesh Yella
- Contributors: (see GitHub contributors page)

---

**Document Version**: 1.0
**Last Updated**: 2024-10-23
**Architecture Version**: 2024.1
