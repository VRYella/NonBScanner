# NonBScanner - Consolidated Non-B DNA Structure Detection Tool

A streamlined tool for detecting A-philic DNA and Z-DNA motifs in genomic sequences.

## Overview

This consolidated version of NonBScanner provides focused functionality for detecting two important classes of Non-B DNA structures:
- **A-philic DNA**: Sequences with high affinity for A-tract formation and protein binding
- **Z-DNA**: Left-handed double helical DNA structures with Z-forming potential

## Core Scripts

The tool consists of 5 main Python scripts with comprehensive tabular annotations:

1. **`aphilic_scanner.py`** - A-philic DNA detection using 10-mer log2 odds scoring
2. **`zdna_scanner.py`** - Z-DNA detection using 10-mer Z-forming potential scoring  
3. **`motif_utils.py`** - Common utilities for sequence processing and analysis
4. **`nbd_scanner.py`** - Main CLI interface for running analyses
5. **`visualization.py`** - Text-based visualization and reporting utilities

## Quick Start

### Basic Usage

```bash
# Analyze a single sequence for both motif types
python nbd_scanner.py --sequence "ACGCGCGCGCGCGCGC" --detect both

# Analyze FASTA file for A-philic DNA only
python nbd_scanner.py --fasta input.fa --detect aphilic --output results.tsv

# Analyze with filtering and summary report
python nbd_scanner.py --fasta input.fa --detect both --filter --format summary
```

### Output Formats

- **TSV**: Tab-separated values for downstream analysis
- **JSON**: Structured data format for programmatic access
- **Summary**: Human-readable statistical report

## Features

### A-philic DNA Detection
- Uses tetranucleotide log2 odds scoring table with 10-mer patterns
- Classifies regions as High/Moderate/Low confidence A-philic
- Scores based on A-tract formation potential and protein binding affinity

### Z-DNA Detection  
- Uses 10-mer Z-forming potential scoring (range: 50-63)
- Incorporates CG dinucleotide content analysis
- Classifies regions based on left-handed helix formation likelihood

### Analysis Features
- Quality filtering by score and length thresholds
- Overlapping region merging to reduce redundancy
- Comprehensive statistical analysis and reporting
- Text-based visualization of motif distributions

## Scientific Background

### A-philic DNA
A-philic DNA represents sequences with high affinity for specific protein binding and structural features that favor A-tract formation. Based on research by Bolshoy et al. (PNAS 1991) and Rohs et al. (Nature 2009).

### Z-DNA
Z-DNA is a left-handed double helical form discovered by Alexander Rich, forming under superhelical tension. Critical for gene regulation and chromatin organization (Rich & Zhang, Nature Reviews 2003).

## Requirements

- Python 3.7+
- NumPy (for numerical computations)
- Standard library modules: csv, json, argparse, pathlib

## Installation

```bash
# Clone repository
git clone https://github.com/VRYella/NonBScanner.git
cd NonBScanner

# Install dependencies
pip install numpy
```

## Examples

### Detect A-philic DNA in sequence
```bash
python nbd_scanner.py --sequence "CCCCCCCCCCAAAAAAAAAA" --detect aphilic --verbose
```

### Analyze FASTA file for Z-DNA with filtering
```bash
python nbd_scanner.py --fasta genome.fa --detect zdna --filter --min-score 55.0 --output zdna_hits.tsv
```

### Generate comprehensive analysis report
```bash
python nbd_scanner.py --fasta sequences.fa --detect both --format summary --merge --filter
```

## Output Fields

### A-philic DNA Results
- **Raw_Score**: Mean log2 odds score from tetranucleotide table
- **Classification**: High/Moderate/Low confidence A-philic
- **N_Tetranucleotides**: Number of scoring tetranucleotides
- **Strong_Count**: Count of high-scoring regions

### Z-DNA Results  
- **Raw_Score**: Mean Z-forming potential score (50-63)
- **CG_Content**: Percentage of CG dinucleotides
- **Classification**: High/Moderate/Low Z-DNA forming potential
- **N_Scoring_10mers**: Number of scoring 10-mer patterns

## Author

Dr. Venkata Rajesh Yella  
Integration: 2024 - NBDFinder Consolidated Implementation

## License

This project is part of the NBDFinder system for comprehensive Non-B DNA structure detection.