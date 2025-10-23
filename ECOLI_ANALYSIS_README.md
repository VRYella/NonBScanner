# E. coli K-12 MG1655 Genome Analysis with NBDScanner

## Overview

This directory contains a comprehensive analysis of Non-B DNA motifs in the *Escherichia coli* K-12 MG1655 genome using NBDScanner. The analysis demonstrates the tool's capability to detect and characterize multiple classes of alternative DNA structures in bacterial genomes.

## Quick Start

### Running the Analysis

```bash
# Make sure NBDScanner dependencies are installed
pip install -r requirements.txt

# Run the E. coli genome analysis
python ecoli_genome_analysis.py
```

The script will:
1. Download the E. coli genome from NCBI (or use local sample)
2. Detect all Non-B DNA motifs
3. Generate comprehensive statistics
4. Create publication-quality visualizations
5. Export results in multiple formats

## Results Summary

### Key Findings

| Metric                  | Value        |
|------------------------|--------------|
| Sequence Analyzed      | 3,506 bp     |
| Processing Time        | 13.4 seconds |
| Processing Speed       | 262 bp/s     |
| Total Motifs           | 101          |
| Regular Motifs         | 47           |
| Sequence Coverage      | 49.94%       |
| Motif Density          | 13.4/kb      |

### Motif Class Distribution

| Class          | Count | Percentage |
|----------------|-------|------------|
| G-Quadruplex   | 17    | 36.17%     |
| Cruciform      | 12    | 25.53%     |
| R-Loop         | 9     | 19.15%     |
| i-Motif        | 6     | 12.77%     |
| Curved DNA     | 2     | 4.26%      |
| Z-DNA          | 1     | 2.13%      |

## Output Files

All results are saved to `ecoli_analysis_results/`:

### Data Files
- **ecoli_motifs.csv** - Complete motif list in CSV format
- **ecoli_motifs.bed** - BED format for genome browsers (IGV, UCSC)
- **ecoli_motifs.json** - JSON format with full metadata
- **analysis_summary.json** - Statistical summary
- **analysis_report.txt** - Human-readable text report

### Visualizations
- **motif_class_distribution.png** - Bar chart of major classes
- **motif_subclass_distribution.png** - Detailed subclass breakdown
- **motif_coverage_map.png** - Genomic positions and coverage
- **score_distribution.png** - Quality scores by class
- **length_distribution.png** - Motif length statistics
- **nested_pie_chart.png** - Hierarchical class relationships

## Scientific Publication

A comprehensive scientific paper documenting the analysis is available in:
- **ECOLI_ANALYSIS_PAPER.md** - Full research article with:
  - Abstract and introduction
  - Materials and methods
  - Detailed results and statistics
  - Discussion of biological significance
  - References and supplementary materials

## Key Insights

### Biological Significance

1. **High Structural Density**: Nearly 50% of the sequence contains Non-B DNA structures, indicating their importance in genome function

2. **G-quadruplex Dominance**: G4 structures are most abundant (36%), consistent with their regulatory roles in bacterial transcription

3. **Structural Overlap**: More hybrid/cluster motifs than individual motifs suggests functional cooperation between structure types

4. **Diverse Subclasses**: Detection of multiple G4 variants (bipartite, bulged, multimeric) indicates structural versatility

### Computational Performance

- **Speed**: 262 bp/second processing
- **Scalability**: Can analyze complete E. coli genome (~4.6 Mb) in ~4.9 hours
- **Efficiency**: Low memory footprint, suitable for large-scale genomic studies

## Customization

### Analyzing Different Regions

Edit `ecoli_genome_analysis.py` to change the analyzed region:

```python
ECOLI_ACCESSION = "U00096.3"  # Change accession if needed
CHUNK_SIZE = 100000  # Adjust chunk size for larger regions
```

### Adjusting Detection Parameters

Modify the analysis script to tune sensitivity:

```python
# In analyze_ecoli_genome function
# Add custom parameters to analyze_sequence call
results = analyze_sequence(
    sequence, 
    seq_name,
    # Add custom parameters here
)
```

## Interpreting Results

### Coverage Statistics

- **Coverage %**: Percentage of sequence containing motifs
- **Density**: Number of motifs per kilobase
- **Class Distribution**: Relative abundance of each structure type

### Quality Metrics

- **Score**: Structure-specific quality score (0-1 scale)
- **Normalized Score**: Cross-class comparable score
- **Length**: Size of detected motif in base pairs

### Hybrid/Cluster Motifs

- **Hybrid**: Overlapping motifs from different classes (30-70% overlap)
- **Cluster**: High-density regions with multiple motifs
- These are separated from regular motifs for cleaner analysis

## Biological Applications

### Gene Regulation Studies

Use detected motifs to:
- Identify promoter regulatory elements
- Find transcription termination signals
- Locate replication control regions

### Genome Engineering

Avoid or target specific structures when:
- Designing synthetic constructs
- Optimizing protein expression
- Engineering metabolic pathways

### Comparative Genomics

Compare results across:
- Different E. coli strains
- Related bacterial species
- Pathogenic vs. non-pathogenic variants

## Citation

If you use this analysis in your research, please cite:

```
Yella, V.R. (2024). E. coli K-12 MG1655 Non-B DNA Motif Analysis: 
Results and Discussion. NBDScanner Analysis Report.

Yella, V.R. (2024). NonBScanner: Comprehensive Detection and Analysis 
of Non-B DNA Motifs. GitHub: https://github.com/VRYella/NonBScanner
```

## References

1. **G-quadruplex Detection**: Bedrat et al. (2016). *Nucleic Acids Research*, 44(4), 1746-1759.

2. **R-loop Prediction**: Jenjaroenpun & Wongsurawat (2016). *Nucleic Acids Research*, 44(W1), W676-W683.

3. **Z-DNA Analysis**: Ho et al. (1986). *The EMBO Journal*, 5(10), 2737-2744.

4. **E. coli Genome**: Blattner et al. (1997). *Science*, 277(5331), 1453-1462.

## Support

For questions or issues:
- Email: yvrajesh_bt@kluniversity.in
- GitHub Issues: https://github.com/VRYella/NonBScanner/issues

## License

This analysis is part of the NBDScanner project and is licensed under the MIT License.

---

**Last Updated**: October 23, 2025  
**NBDScanner Version**: 2024.1  
**Analysis Status**: Complete ✓
