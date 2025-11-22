# Enhanced Statistics: Density and Enrichment Analysis

## Overview

This update adds comprehensive density and enrichment analysis to NonBScanner, implementing rigorous statistical methods for quantifying and validating motif prevalence.

## New Features

### 1. 📏 Density Calculations

Density metrics quantify the prevalence of non-B DNA motifs across a given genomic region.

#### A. Genomic Density (Coverage)

Reports the percentage of the total sequence length covered by predicted non-B DNA motifs.

**Formula:**
```
Genomic Density (σ_G) = (Total length in bp of all predicted motifs / Total length in bp of analyzed region) × 100
```

**Example:**
- If Z-DNA motifs collectively cover 30 bp of a 1000 bp sequence
- Genomic density = (30 / 1000) × 100 = 3.0%

#### B. Positional Density (Frequency)

Reports how often a motif occurs per unit of length, useful for comparing density between different-sized regions.

**Formula:**
```
Positional Density (λ) = Total count of predicted motifs / Total length (in kbp or Mbp) of analyzed region
```

**Example:**
- If you find 15 G-quadruplexes in a 50 kbp sequence
- Positional density = 15 / 50 = 0.3 motifs per kbp (or 300 motifs per Mbp)

### 2. ✨ Enrichment Analysis

Enrichment analysis determines if a motif occurs non-randomly (more or less often) compared to a background distribution generated through sequence shuffling.

#### A. Fold Enrichment

Quantifies the proportional increase (or decrease) in the motif's density relative to the background.

**Formula:**
```
Fold Enrichment = D_Observed / D_Background

where D = Motif Density = Total bp of motif / Total bp of region
```

**Calculation Steps:**
1. Calculate observed density in the original sequence
2. Generate background distribution through 100 sequence shuffles
3. Calculate mean background density
4. Compute fold enrichment ratio

**Interpretation:**
- Fold Enrichment = 2.5 → motif is 2.5× more frequent than expected
- Fold Enrichment = 0.5 → motif is depleted (half the expected frequency)
- Fold Enrichment = 1.0 → motif occurs at expected frequency

#### B. Statistical Significance

To confirm that observed enrichment is not due to random chance, p-values are calculated using shuffling-based permutation testing.

**Method:**
1. Shuffle the input sequence 100 times (preserving nucleotide composition)
2. Detect motifs in each shuffled sequence
3. Calculate density for each shuffle
4. P-value = proportion of shuffled sequences with density ≥ observed density

**Interpretation:**
- p < 0.001: Highly significant (***) 
- p < 0.01: Very significant (**)
- p < 0.05: Significant (*)
- p ≥ 0.05: Not significant (ns)

## Usage in the Web Application

### Accessing Density Metrics

1. Navigate to the **Statistics** tab in the Results section
2. Select the **Density Metrics** sub-tab
3. View:
   - Overall genomic and positional density
   - Per-class density breakdown in table format
   - Visualization comparing genomic vs positional density

### Running Enrichment Analysis

1. Navigate to the **Statistics** tab
2. Select the **Enrichment Analysis** sub-tab
3. Click the **"Run Enrichment Analysis (100 shuffles)"** button
4. Wait for analysis to complete (may take 1-5 minutes depending on sequence length)
5. View results:
   - Summary table with fold enrichment and p-values
   - Visualization panels showing:
     - Fold enrichment by class
     - Statistical significance (p-values)
     - Observed vs. background density comparison
   - Formatted summary table

## Implementation Details

### Functions Added

**In `utilities.py`:**
- `shuffle_sequence()`: Generates shuffled sequences preserving composition
- `calculate_genomic_density()`: Computes coverage percentage
- `calculate_positional_density()`: Computes motifs per unit length
- `calculate_enrichment_with_shuffling()`: Performs enrichment analysis
- `calculate_enhanced_statistics()`: Comprehensive statistics wrapper

**In `visualizations.py`:**
- `plot_density_comparison()`: Side-by-side density visualization
- `plot_enrichment_analysis()`: Multi-panel enrichment plots
- `plot_enrichment_summary_table()`: Formatted results table

**In `app.py`:**
- Enhanced Statistics tab with three sub-tabs:
  - Density Metrics
  - Enrichment Analysis
  - Distributions

### Performance Considerations

- Enrichment analysis with 100 shuffles typically takes 1-5 minutes
- Progress bar shows real-time status during analysis
- Results are cached in session state to avoid re-computation
- Density calculations are near-instantaneous

## Scientific Background

The enrichment analysis methodology is based on:

1. **Permutation Testing**: Non-parametric statistical method that doesn't assume a particular distribution
2. **Null Hypothesis**: Motifs are distributed randomly in the sequence
3. **Alternative Hypothesis**: Motifs show non-random enrichment/depletion
4. **Background Generation**: Sequence shuffling preserves nucleotide composition while disrupting motif patterns

This approach is commonly used in genomics for:
- ChIP-seq peak enrichment
- Motif over-representation analysis
- Regulatory element identification
- Functional genomics studies

## Validation

The implementation has been tested with:
- Unit tests for density calculations
- Unit tests for sequence shuffling
- End-to-end tests with real sequences
- Visualization tests for all plot types

All tests pass successfully, confirming:
- Accurate density calculations
- Proper nucleotide composition preservation during shuffling
- Correct enrichment and p-value computation
- Publication-quality visualizations with no text overlap

## Example Results

For a typical sequence analysis:

```
Genomic Density:
  Z-DNA: 4.0%
  G-Quadruplex: 4.5%
  i-Motif: 1.8%

Positional Density (motifs/kbp):
  Z-DNA: 2.0
  G-Quadruplex: 2.0
  i-Motif: 1.0

Enrichment Results:
  i-Motif: 
    - Fold Enrichment: 3.6
    - P-value: 0.001 (***)
    - Interpretation: Highly enriched
  
  Z-DNA:
    - Fold Enrichment: 2.0
    - P-value: 0.02 (*)
    - Interpretation: Significantly enriched
  
  G-Quadruplex:
    - Fold Enrichment: 1.5
    - P-value: 0.15 (ns)
    - Interpretation: Not significantly enriched
```

## References

These metrics are based on standard methods in genomics and bioinformatics:

1. Genomic density calculations: Standard coverage metrics
2. Fold enrichment: Commonly used in ChIP-seq and motif analysis
3. Permutation testing: Non-parametric statistical method
4. P-value calculation: Empirical distribution-based significance testing

---

**Author:** Dr. Venkata Rajesh Yella  
**Date:** 2024  
**Module Version:** 2024.1 Enhanced
