# NonBScanner Organization & Architecture

## Overview

NonBScanner is organized based on **hyperscan usage** for optimal performance. This document describes the structure and organization of the codebase.

## Repository Structure

### Core Application
- **`app.py`** - Main Streamlit web application
- **`requirements.txt`** - Python dependencies including hyperscan

### Documentation (Minimal Set)
- **`README.md`** - Main project documentation
- **`QUICK_START.md`** - Quick start guide
- **`HYPERSCAN_ARCHITECTURE.md`** - Hyperscan integration details
- **`OPTIMIZED_SCANNER_ARCHITECTURE.md`** - Optimized scanner architecture
- **`ORGANIZATION.md`** - This file

### Motif Detection Classes (`motif_detection/`)

Detectors are organized by their detection approach:

#### Hyperscan-Based Detectors (Fast Pattern Matching)
- **`z_dna_detector.py`** - Uses 126 pre-compiled 10-mer patterns
- **`a_philic_detector.py`** - Uses 208 pre-compiled 10-mer patterns
- **`g_quadruplex_detector.py`** - Uses 7 regex patterns
- **`i_motif_detector.py`** - Uses regex patterns
- **`curved_dna_detector.py`** - Uses 44 regex patterns for A-tract detection
- **`r_loop_detector.py`** - Uses regex patterns + QmRLFS algorithm

#### Algorithmic Detectors (Optimized Python)
- **`cruciform_detector.py`** - K-mer based inverted repeat detection (O(n))
- **`slipped_dna_detector.py`** - K-mer based direct repeat and STR detection (O(n))
- **`triplex_detector.py`** - K-mer based mirror repeat detection (O(n))

#### Base Class
- **`base_detector.py`** - Abstract base class for all detectors

### Utility Modules (`utils/`)

#### Hyperscan Integration
- **`load_hsdb.py`** - Load pre-compiled Hyperscan databases
- **`load_regex_registry.py`** - Load and compile regex pattern registries
- **`motif_patterns.py`** - Pattern definitions and Hyperscan compatibility checks

#### Scanner Modules
- **`modular_scanner.py`** - Main scanner coordinating all detectors
- **`nbdscanner.py`** - High-level analysis functions
- **`repeat_scanner.py`** - Optimized repeat detection algorithms

#### Support Modules
- **`utils.py`** - General utility functions (FASTA parsing, exports, etc.)
- **`visualization.py`** - Plotting and visualization functions
- **`canonicalize_motif.py`** - Motif canonicalization utilities

### Pattern Registries (`registry/`)

Each motif class has registry files in multiple formats:

**Format Priority:**
1. `.pkl` - Primary format (fastest to load)
2. `.json` - Fallback format (human-readable)

**Registry Files:**
- APhilic_registry.pkl/json
- Cruciform_registry.pkl/json
- CurvedDNA_registry.pkl/json
- G4_registry.pkl/json
- IMotif_registry.pkl/json
- RLoop_registry.pkl/json
- SlippedDNA_registry.pkl/json
- Triplex_registry.pkl/json
- ZDNA_registry.pkl/json

**Note:** `.txt` pattern files removed (redundant - patterns are in registries)

### Tests (`tests/`)

#### Essential Tests (Kept)
- **`test_all_motif_classes.py`** - Comprehensive test of all 10+ motif classes
  - Tests hyperscan database loading and usage
  - Tests all detector types
  - Tests hybrid/cluster detection
  - Tests scoring and overlap resolution
  
- **`test_pattern_registries.py`** - Validates pattern registries
  - Tests registry loading (pkl/json formats)
  - Tests hyperscan compilation
  - Tests registry metadata

- **`test_complete_pipeline.py`** - Full pipeline integration test
  - Tests detection → scoring → overlap resolution → visualization

- **`test_hyperscan_performance.py`** - NEW! Performance benchmarks
  - Tests hyperscan availability and DB loading
  - Measures detection speed and scalability
  - Provides motif class breakdown
  - Compares hyperscan vs algorithmic approaches

#### Additional Tests (Kept for Completeness)
- **`test_optimized_repeat_scanner.py`** - Tests optimized repeat algorithms
- **`test_real_genome_sequence.py`** - Tests with realistic sequences

### Tools (`tools/`)
- **`generate_pattern_registries.py`** - Script to generate registry files

## Hyperscan Integration

### What is Hyperscan?

Hyperscan is a high-performance regular expression matching library from Intel. NonBScanner uses it for:
- Ultra-fast pattern matching (10-100x faster than Python regex)
- Efficient 10-mer lookup for Z-DNA and A-philic detectors
- Pre-compiled databases for instant startup

### Hyperscan vs Algorithmic Approaches

| Approach | Classes | Method | Performance |
|----------|---------|--------|-------------|
| **Hyperscan** | Z-DNA, A-philic, G4, i-Motif, Curved DNA, R-Loop | Pre-compiled pattern databases | 15,000-65,000 bp/s |
| **Algorithmic** | Cruciform, Slipped DNA, Triplex | K-mer indexing + seed-and-extend | 5,000-280,000 bp/s |

### Loading Priority

1. Pre-compiled `.hsdb` files (if available)
2. `.pkl` registry files → compile to hyperscan
3. `.json` registry files → compile to hyperscan
4. Fallback to pure Python regex (if hyperscan unavailable)

## Files Removed in Organization

### Redundant Documentation (6 files)
- `IMPLEMENTATION_COMPLETE.md` - Development notes
- `IMPLEMENTATION_SUMMARY.md` - Development notes
- `PKL_REGISTRY_IMPLEMENTATION.md` - Development notes
- `QUICKSTART_REQUIREMENTS.md` - Development notes
- `REQUIREMENTS_IMPLEMENTATION.md` - Development notes
- `SUMMARY.md` - Duplicate summary

### Unused Code (2 files)
- `utils/advanced_visualizations.py` - Not imported anywhere
- `utils/detector_registry.py` - Not imported anywhere

### Demo Files (2 files)
- `demo_registries.py` - Demo script
- `test_requirements.py` - Test script

### Redundant Registry Files (9 files)
- `registry/*_patterns.txt` - Patterns already in .pkl/.json files

**Total: 19 files removed, ~3,200 lines deleted**

## Performance Characteristics

### Detection Rates (Measured)

| Sequence Length | Detection Rate | Notes |
|----------------|----------------|-------|
| 1,000 bp | ~65,000 bp/s | Fast hyperscan pattern matching |
| 5,000 bp | ~15,000 bp/s | Mixed hyperscan + algorithmic |
| 10,000 bp | ~4,000 bp/s | Cruciform becomes bottleneck |
| 25,000 bp | ~750 bp/s | Multiple overlap detection phases |

### Bottleneck Analysis

The main performance bottleneck for larger sequences is:
1. **Cruciform detection** - O(n) but with high constant factor
2. **Overlap resolution** - Must compare all detected motifs
3. **Hybrid/Cluster detection** - Secondary analysis phase

### Optimization Opportunities

For sequences > 10kb:
- Consider disabling Cruciform detection if not needed
- Use streaming analysis for genome-scale sequences
- Pre-filter regions of interest

## Testing Strategy

### Test Organization

1. **Unit Tests** - Individual detector functionality
2. **Integration Tests** - Full pipeline testing
3. **Performance Tests** - Speed and scalability benchmarks
4. **Registry Tests** - Pattern loading and validation

### Running Tests

```bash
# Run all tests
python -m pytest tests/ -v

# Run essential tests only
python -m pytest tests/test_all_motif_classes.py tests/test_pattern_registries.py -v

# Run performance benchmarks
python tests/test_hyperscan_performance.py

# Run specific test
python -m pytest tests/test_hyperscan_performance.py::test_hyperscan_availability -v
```

### Test Coverage

- **19 tests total** across 6 test files
- All tests passing
- Coverage includes:
  - All 10+ motif classes
  - Hyperscan integration
  - Pattern registries
  - Complete pipeline
  - Performance benchmarks

## Development Workflow

### Adding a New Motif Detector

1. Create detector class in `motif_detection/` inheriting from `BaseDetector`
2. Decide on detection method:
   - If simple patterns → use Hyperscan
   - If complex structure → use algorithmic approach
3. Add patterns to `registry/` directory
4. Register detector in `utils/modular_scanner.py`
5. Add tests in `tests/`

### Modifying Patterns

1. Edit pattern registry files in `registry/` (.pkl or .json)
2. Optionally pre-compile to `.hsdb` using hyperscan
3. Restart application to load new patterns

### Performance Tuning

1. Run `tests/test_hyperscan_performance.py` to establish baseline
2. Make changes
3. Re-run performance tests to measure improvement
4. Document changes in commit message

## Maintenance

### Regular Tasks

- **Update registries** - Add new validated patterns
- **Optimize algorithms** - Improve bottleneck detectors
- **Update dependencies** - Keep hyperscan and other libs current
- **Add tests** - Cover new functionality

### Code Health

- All imports verified (no unused modules)
- All tests passing (19/19)
- Documentation current
- Performance benchmarks established

## References

- **Hyperscan**: https://www.hyperscan.io/
- **Pattern Sources**: See individual detector files for references
- **Algorithms**: See `OPTIMIZED_SCANNER_ARCHITECTURE.md` for details

---

*Last Updated: 2024 - NonBScanner organized for hyperscan-based performance*
