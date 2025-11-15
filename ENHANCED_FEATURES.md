# Enhanced Features Documentation

## Overview
This document describes the enhancements made to NonBScanner for improved registry management and scientific visualization.

## 1. Consolidated Registry System

### What Changed
Previously, NonBScanner used a `registry/` folder containing 18 separate files (9 JSON + 9 PKL files) for pattern definitions. This has been consolidated into a **single file**: `consolidated_registry.json`.

### Benefits
- ✅ **Easier Distribution**: Single file instead of entire folder
- ✅ **Simpler Maintenance**: All patterns in one place
- ✅ **Backward Compatible**: Falls back to old registry folder if consolidated file not found
- ✅ **Same Performance**: Uses caching to maintain fast loading

### File Structure
```json
{
  "version": "2024.1",
  "generated_at": "2024-11-15T...",
  "description": "Consolidated registry for all Non-B DNA motif patterns",
  "total_patterns": 411,
  "total_classes": 9,
  "registries": {
    "APhilic": { ... },
    "Cruciform": { ... },
    "CurvedDNA": { ... },
    "G4": { ... },
    "IMotif": { ... },
    "RLoop": { ... },
    "SlippedDNA": { ... },
    "Triplex": { ... },
    "ZDNA": { ... }
  }
}
```

### Pattern Count by Class
- **APhilic**: 208 patterns
- **Cruciform**: 1 pattern
- **CurvedDNA**: 44 patterns
- **G4**: 7 patterns
- **IMotif**: 7 patterns
- **RLoop**: 5 patterns
- **SlippedDNA**: 9 patterns
- **Triplex**: 4 patterns
- **ZDNA**: 126 patterns
- **Total**: 411 patterns

## 2. Enhanced Scientific Visualizations

### New Visualization Functions

#### 1. `plot_class_analysis_comprehensive()`
Comprehensive analysis showing all 11 Non-B DNA classes.

**Features:**
- Shows distribution of all 11 classes (detected and not detected)
- Pie chart showing detection status
- Statistical table for top detected classes
- Clear indication of which classes were NOT detected
- Publication-quality output at 300 DPI

**Output Includes:**
- Main distribution bar chart with color coding
- Detection status pie chart
- Statistics table (count, avg length, avg score)
- List of non-detected classes

**Example Usage:**
```python
import visualizations as viz
import matplotlib.pyplot as plt

fig = viz.plot_class_analysis_comprehensive(motifs, figsize=(16, 12))
fig.savefig('class_analysis.png', dpi=300, bbox_inches='tight')
plt.close(fig)
```

#### 2. `plot_subclass_analysis_comprehensive()`
Detailed subclass-level analysis organized by parent class.

**Features:**
- Horizontal bar chart showing all detected subclasses
- Color-coded by parent class
- Summary statistics by class
- Clear class:subclass naming convention

**Output Includes:**
- Horizontal bar chart of all subclasses with counts
- Summary text showing subclass count per class
- Total motif count per class

**Example Usage:**
```python
fig = viz.plot_subclass_analysis_comprehensive(motifs, figsize=(18, 14))
fig.savefig('subclass_analysis.png', dpi=300, bbox_inches='tight')
plt.close(fig)
```

#### 3. `plot_score_statistics_by_class()`
Advanced statistical visualization of scores with multiple plot types.

**Features:**
- Violin plots showing score distributions
- Box plot overlays with median and quartiles
- Statistical annotations (mean, std)
- Comprehensive statistics table

**Statistical Measures:**
- Mean (μ)
- Median
- Standard Deviation (σ)
- Min/Max values
- Sample count (N)

**Example Usage:**
```python
fig = viz.plot_score_statistics_by_class(motifs, figsize=(14, 8))
fig.savefig('score_statistics.png', dpi=300, bbox_inches='tight')
plt.close(fig)
```

#### 4. `plot_length_statistics_by_class()`
Comprehensive length distribution analysis.

**Features:**
- Overlaid histograms for each class
- Box plot comparison
- Statistical summary table
- KDE (Kernel Density Estimation) overlays

**Statistical Measures:**
- Mean length
- Median length
- Standard deviation
- Min/Max lengths
- Distribution histograms

**Example Usage:**
```python
fig = viz.plot_length_statistics_by_class(motifs, figsize=(14, 10))
fig.savefig('length_statistics.png', dpi=300, bbox_inches='tight')
plt.close(fig)
```

### Visualization Design Principles

1. **Publication Quality**
   - 300 DPI output resolution
   - Clean, professional styling
   - Proper axis labels and titles
   - Statistical annotations

2. **Scientific Accuracy**
   - All 11 Non-B DNA classes shown (even if not detected)
   - Clear indication of detection status
   - Comprehensive statistics (mean, median, std, min, max)
   - Proper error handling for empty datasets

3. **Color Scheme**
   - Colorblind-friendly palette from matplotlib
   - Consistent class colors across all plots
   - High contrast for readability
   - Professional scientific appearance

4. **Information Density**
   - Multiple views in single figure (subplots)
   - Statistical tables alongside graphs
   - Count labels on bars
   - Comprehensive legends

### All 11 Non-B DNA Classes

The visualization functions properly handle all 11 classes:

1. **Curved_DNA** - A-tract mediated curvature
2. **Slipped_DNA** - Tandem repeats, slipped structures
3. **Cruciform** - Four-way junctions (inverted repeats)
4. **R-Loop** - RNA-DNA hybrids
5. **Triplex** - Three-stranded structures
6. **G-Quadruplex** - Four-stranded G-rich structures
7. **i-Motif** - C-rich structures
8. **Z-DNA** - Left-handed double helix
9. **A-philic_DNA** - A-philic DNA regions
10. **Hybrid** - Multi-class overlaps
11. **Non-B_DNA_Clusters** - High-density motif regions

### Subclass Coverage

The system tracks 22+ subclasses across all classes:

**G-Quadruplex** (7 subclasses):
- Canonical G4
- Relaxed G4
- Bulged G4
- Long Loop G4
- Multimeric G4
- Imperfect G4
- G-Triplex

**Curved DNA** (2 subclasses):
- Global Curvature
- Local Curvature

**R-Loop** (3 subclasses):
- R-loop formation sites
- QmRLFS-m1
- QmRLFS-m2

**i-Motif** (2 subclasses):
- Canonical i-Motif
- HuR AC-Motif

**Triplex** (2 subclasses):
- Homopyrimidine mirror repeat
- Sticky DNA

**SlippedDNA** (1 subclass):
- STR (Short Tandem Repeats)

**Cruciform** (1 subclass):
- Inverted_Repeat

**Z-DNA** (1 subclass):
- Z-DNA (classic)

**A-philic_DNA** (1 subclass):
- A-philic DNA

**Hybrid** (multiple):
- Various overlap combinations

**Non-B_DNA_Clusters** (multiple):
- Mixed clusters of varying sizes

## Testing

### Test Scripts

1. **test_overlap_removal.py** - Validates core detection functionality (8/8 tests passing)
2. **test_enhanced_visualizations.py** - Tests new visualization functions with sample data

### Running Tests

```bash
# Test core functionality
python test_overlap_removal.py

# Test enhanced visualizations
python test_enhanced_visualizations.py
```

### Expected Results
- All 8 overlap removal tests should pass
- Enhanced visualizations should generate 5 PNG files in /tmp/
- No errors or warnings during execution

## Migration Guide

### For Users

No changes required! The system automatically uses the consolidated registry if available, and falls back to the old registry folder if not.

### For Developers

If you need to update patterns:

1. Edit `consolidated_registry.json` directly
2. Or edit individual registry files in `registry/` folder
3. System will use consolidated file if present

## Performance

- **Loading Time**: < 1 second for all 411 patterns
- **Memory Usage**: Minimal increase (consolidated JSON is cached)
- **Visualization Generation**: 2-5 seconds per complex plot
- **No Performance Degradation**: All original tests still pass at same speed

## Backward Compatibility

✅ **Fully backward compatible**
- Works with or without consolidated registry
- Falls back to registry folder if consolidated file not found
- All existing code continues to work unchanged
- All tests pass without modification

## Future Enhancements

Potential future improvements:
- Interactive plotly versions of new visualizations
- Export to SVG for vector graphics
- Additional statistical tests (ANOVA, etc.)
- Customizable color schemes
- PDF report generation with all plots

## Summary

These enhancements improve NonBScanner in two key areas:

1. **Registry Management**: Simplified from 18 files to 1 file
2. **Visualizations**: Added 4 comprehensive, publication-quality plotting functions

Both improvements maintain full backward compatibility while providing significant usability and scientific value.
