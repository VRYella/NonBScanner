# Hybrid and Cluster Motif Separation - NBDScanner

## Overview

This document describes the implementation of the hybrid and cluster motif separation feature in NBDScanner. This feature addresses the user requirement to display hybrid and cluster motifs separately from regular Non-B DNA motifs.

## Problem Statement

The original issue requested:
- Don't show Hybrid as a separate class in main results
- Merge definitions of hybrid and non-B DNA clusters
- Show them only in a separate visualizations tab
- Exclude them from downloads and exports
- Keep the longest non-overlapping selection logic

## Implementation

### Changes Made

#### 1. Main Results Display (`app.py`)

**Filtering Logic:**
```python
# Filter motifs for main display (exclude hybrid and cluster)
filtered_motifs = [m for m in motifs if m.get('Class') not in ['Hybrid', 'Non-B_DNA_Clusters']]
hybrid_cluster_motifs = [m for m in motifs if m.get('Class') in ['Hybrid', 'Non-B_DNA_Clusters']]
```

**Statistics Update:**
- Total motifs count now shows only regular motifs (excluding hybrid/cluster)
- Coverage and density calculations use only filtered motifs
- Info message indicates how many hybrid/cluster motifs were detected

**Example:**
```
Total Motifs: 52 (regular Non-B DNA motifs)
ℹ️ 71 Hybrid/Cluster motifs detected. View them in the 'Cluster/Hybrid' tab below.
```

#### 2. Visualization Tabs

**New Tab Structure:**
- 📈 Distribution (regular motifs only)
- 🗺️ Coverage Map (regular motifs only)
- 📊 Statistics (regular motifs only)
- 🎯 Interactive (regular motifs only)
- 🔗 **Cluster/Hybrid** (NEW - dedicated tab)

**Cluster/Hybrid Tab Features:**
- Informational text explaining what hybrid and cluster motifs are
- Summary statistics (hybrid count, cluster count, average length)
- Detailed table with relevant columns
- Position map visualization
- Score distribution visualization

#### 3. Download/Export Exclusion

**Export Filtering:**
```python
# Separate hybrid/cluster from regular motifs
if export_motif.get('Class') in ['Hybrid', 'Non-B_DNA_Clusters']:
    hybrid_cluster_motifs_export.append(export_motif)
else:
    all_motifs.append(export_motif)
```

**User Feedback:**
```
ℹ️ 71 Hybrid/Cluster motifs are excluded from downloads. 
These are shown only in the Cluster/Hybrid visualization tab.

Showing first 10 of 52 total records (Hybrid/Cluster motifs excluded)
```

## User Experience

### Before
- Hybrid and cluster motifs were mixed with regular motifs in all views
- Downloads included all motif types
- No clear separation between simple and complex motifs
- Confusion about overlapping regions

### After
- Clean separation between regular and complex motifs
- Regular motifs in main results (52 motifs)
- Hybrid/cluster in dedicated tab (71 motifs)
- Clear messaging about where to find each type
- Downloads contain only regular motifs for cleaner analysis

## Technical Details

### Motif Classes Separated

1. **Hybrid Motifs** (`Class: 'Hybrid'`)
   - Regions where different non-B DNA classes overlap (30-70%)
   - Example: `R-Loop_Cruciform_Overlap`
   - Indicates complex genomic regions with multiple features

2. **Cluster Motifs** (`Class: 'Non-B_DNA_Clusters'`)
   - High-density regions with multiple motifs from different classes
   - Example: `Mixed_Cluster_10_classes`
   - Indicates hotspots of non-B DNA activity

### Data Flow

```
analyze_sequence()
    ↓
All motifs (123 total)
    ↓
Filtered motifs (52 regular) + Hybrid/Cluster (71)
    ↓
Main Display: 52 regular motifs
    ↓
Cluster/Hybrid Tab: 71 hybrid/cluster motifs
    ↓
Downloads: 52 regular motifs only
```

## Benefits

1. **Cleaner Results**: Main results show only primary Non-B DNA motifs
2. **Focused Analysis**: Users can analyze regular motifs without confusion
3. **Advanced Users**: Can still access hybrid/cluster data in dedicated tab
4. **Better Downloads**: Export files contain cleaner datasets for downstream analysis
5. **Clear Communication**: Users understand the separation through info messages

## Testing

The implementation was tested with:
- Example single FASTA sequence (516 bp)
- Multi-FASTA examples
- Advanced feature test suite
- UI verification through Streamlit interface

**Test Results:**
- ✅ 52 regular motifs displayed in main results
- ✅ 71 hybrid/cluster motifs in dedicated tab
- ✅ Downloads exclude hybrid/cluster motifs
- ✅ All visualizations use correct filtered data
- ✅ Info messages clearly communicate the separation

## Future Enhancements

Potential future improvements:
1. Optional toggle to include/exclude hybrid/cluster from downloads
2. Advanced filtering options in the Cluster/Hybrid tab
3. Comparative visualization showing regular vs hybrid/cluster
4. Export option specifically for hybrid/cluster motifs

## References

- `app.py`: Main UI implementation
- `nbdscanner.py`: Core detection logic (unchanged)
- `ADVANCED_FEATURES.md`: Documentation on hybrid/cluster detection
- `test_advanced_features.py`: Test suite
