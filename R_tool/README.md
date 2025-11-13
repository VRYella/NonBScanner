# NonBScanner R Interface

**R Wrapper for NonBScanner - Non-B DNA Motif Detection Tool**

---

## 📋 Overview

The NonBScanner R interface provides seamless integration between R and the Python-based NonBScanner tool for detecting Non-B DNA motifs. This wrapper allows R users to leverage the comprehensive motif detection capabilities of NonBScanner directly from their R environment.

## 🎯 Features

- **11 Major Motif Classes** with 22+ subclasses
- **Easy R Integration** using reticulate
- **Comprehensive Analysis** of single or multiple sequences
- **FASTA File Support** for batch processing
- **Visualization Tools** using ggplot2
- **Export Capabilities** to CSV and other formats

## 📦 Installation

### Prerequisites

1. **R** (version 4.0+)
2. **Python** (version 3.8+)
3. **NonBScanner** Python tool installed

### Required R Packages

```r
# Install required R packages
install.packages("reticulate")
install.packages("ggplot2")  # Optional, for visualizations
```

### Python Dependencies

Make sure NonBScanner Python dependencies are installed:

```bash
pip install numpy pandas matplotlib seaborn plotly biopython scipy scikit-learn
```

### Setup

1. Clone the NonBScanner repository:
```bash
git clone https://github.com/VRYella/NonBScanner.git
cd NonBScanner/R_tool
```

2. In R, set the Python path (if needed):
```r
library(reticulate)
use_python("/path/to/python3")  # Optional: specify Python path
```

## 🚀 Quick Start

### Basic Usage

```r
# Load the NonBScanner R interface
source("nonbscanner.R")

# Initialize
init_nonbscanner()

# Analyze a sequence
sequence <- "GGGTTAGGGTTAGGGTTAGGG"
motifs <- analyze_sequence(sequence, "my_sequence")

# View results
print(motifs)
```

## 📚 Functions Reference

### Initialization

#### `init_nonbscanner()`

Initialize the NonBScanner Python module.

**Returns:** `TRUE` if successful

**Example:**
```r
init_nonbscanner()
```

---

### Analysis Functions

#### `analyze_sequence(sequence, sequence_name = "sequence")`

Analyze a single DNA sequence for Non-B DNA motifs.

**Parameters:**
- `sequence`: Character string containing DNA sequence
- `sequence_name`: Name identifier for the sequence (default: "sequence")

**Returns:** Data frame with detected motifs

**Example:**
```r
seq <- "GGGTTAGGGTTAGGGTTAGGG"
motifs <- analyze_sequence(seq, "test_seq")
```

---

#### `analyze_multiple_sequences(sequences, use_multiprocessing = FALSE)`

Analyze multiple DNA sequences.

**Parameters:**
- `sequences`: Named list of sequences
- `use_multiprocessing`: Logical, use parallel processing (default: FALSE)

**Returns:** Data frame with all detected motifs

**Example:**
```r
seqs <- list(
  seq1 = "GGGTTAGGGTTAGGGTTAGGG",
  seq2 = "AAAAATTTTAAAAATTTT"
)
motifs <- analyze_multiple_sequences(seqs)
```

---

#### `analyze_fasta(filename, use_multiprocessing = FALSE)`

Read and analyze sequences from a FASTA file.

**Parameters:**
- `filename`: Path to FASTA file
- `use_multiprocessing`: Logical, use parallel processing (default: FALSE)

**Returns:** Data frame with detected motifs

**Example:**
```r
motifs <- analyze_fasta("sequences.fasta")
```

---

### Utility Functions

#### `read_fasta(filename)`

Read sequences from a FASTA file.

**Parameters:**
- `filename`: Path to FASTA file

**Returns:** Named list of sequences

**Example:**
```r
sequences <- read_fasta("sequences.fasta")
```

---

#### `get_motif_info()`

Get information about supported motif classes and subclasses.

**Returns:** List with classification information

**Example:**
```r
info <- get_motif_info()
print(info$classification)
```

---

#### `get_summary_stats(motifs)`

Generate summary statistics for detected motifs.

**Parameters:**
- `motifs`: Data frame of detected motifs

**Returns:** Data frame with summary statistics

**Example:**
```r
summary <- get_summary_stats(motifs)
print(summary)
```

---

### Visualization Functions

#### `plot_motif_distribution(motifs, title = "Non-B DNA Motif Distribution")`

Create a bar plot of motif class distribution.

**Parameters:**
- `motifs`: Data frame of detected motifs
- `title`: Plot title

**Returns:** ggplot2 plot object

**Example:**
```r
p <- plot_motif_distribution(motifs)
print(p)
ggsave("distribution.png", p)
```

---

#### `plot_motif_positions(motifs, sequence_length = NULL, title = "Non-B DNA Motif Positions")`

Create a position plot showing motif locations.

**Parameters:**
- `motifs`: Data frame of detected motifs
- `sequence_length`: Length of the sequence (optional)
- `title`: Plot title

**Returns:** ggplot2 plot object

**Example:**
```r
p <- plot_motif_positions(motifs, sequence_length = 1000)
print(p)
ggsave("positions.png", p)
```

---

## 📊 Output Format

The motif detection results are returned as a data frame with the following columns:

| Column | Description |
|--------|-------------|
| `ID` | Unique motif identifier |
| `Sequence_Name` | Name of the analyzed sequence |
| `Class` | Motif class (e.g., G-Quadruplex, Z-DNA) |
| `Subclass` | Motif subclass |
| `Start` | Start position (1-based) |
| `End` | End position |
| `Length` | Motif length in base pairs |
| `Sequence` | Actual sequence of the motif |
| `Score` | Motif score |
| `Strand` | Strand information |
| `Method` | Detection method used |

## 🔬 Supported Motif Classes

1. **Curved DNA** (2 subclasses)
   - Global curvature
   - Local Curvature

2. **Slipped DNA** (2 subclasses)
   - Direct Repeat
   - STR (Short Tandem Repeat)

3. **Cruciform** (1 subclass)
   - Inverted Repeats

4. **R-Loop** (1 subclass)
   - R-loop formation sites

5. **Triplex** (2 subclasses)
   - Triplex
   - Sticky DNA

6. **G-Quadruplex** (7 subclasses)
   - Multimeric G4
   - Canonical G4
   - Relaxed G4
   - Bulged G4
   - Bipartite G4
   - Imperfect G4
   - G-Triplex intermediate

7. **i-Motif** (3 subclasses)
   - Canonical i-motif
   - Relaxed i-motif
   - AC-motif

8. **Z-DNA** (2 subclasses)
   - Z-DNA
   - eGZ (Extruded-G) DNA

9. **A-philic DNA** (1 subclass)

10. **Hybrid** (dynamic)
    - Multi-class overlaps

11. **Non-B DNA Clusters** (dynamic)
    - High-density regions

## 💡 Example Workflows

### Workflow 1: Single Sequence Analysis

```r
# Load and initialize
source("nonbscanner.R")
init_nonbscanner()

# Define sequence
my_seq <- "GGGTTAGGGTTAGGGTTAGGGAAAAATTTTAAAAATTTT"

# Analyze
motifs <- analyze_sequence(my_seq, "my_sequence")

# Get summary
summary <- get_summary_stats(motifs)
print(summary)

# Export results
write.csv(motifs, "results.csv", row.names = FALSE)
```

### Workflow 2: FASTA File Analysis with Visualization

```r
# Load and initialize
source("nonbscanner.R")
init_nonbscanner()

# Analyze FASTA file
motifs <- analyze_fasta("sequences.fasta")

# Create visualizations
library(ggplot2)
p1 <- plot_motif_distribution(motifs)
ggsave("distribution.png", p1, width = 10, height = 6)

p2 <- plot_motif_positions(motifs)
ggsave("positions.png", p2, width = 12, height = 8)

# Export by class
for (cls in unique(motifs$Class)) {
  if (cls != "NA") {
    class_motifs <- motifs[motifs$Class == cls, ]
    filename <- paste0(gsub(" ", "_", cls), ".csv")
    write.csv(class_motifs, filename, row.names = FALSE)
  }
}
```

### Workflow 3: Batch Processing

```r
# Load and initialize
source("nonbscanner.R")
init_nonbscanner()

# Define multiple sequences
sequences <- list(
  gene1 = "GGGTTAGGGTTAGGGTTAGGG...",
  gene2 = "AAAAATTTTAAAAATTTT...",
  gene3 = "CGCGCGCGCGCGCG..."
)

# Analyze all sequences
all_motifs <- analyze_multiple_sequences(sequences)

# Summary by sequence
for (seq_name in names(sequences)) {
  seq_motifs <- all_motifs[all_motifs$Sequence_Name == seq_name, ]
  cat(seq_name, ":", nrow(seq_motifs), "motifs\n")
}

# Export combined results
write.csv(all_motifs, "batch_results.csv", row.names = FALSE)
```

## 🔧 Troubleshooting

### Common Issues

1. **Python module not found**
   ```r
   # Set Python path explicitly
   library(reticulate)
   use_python("/usr/bin/python3")
   ```

2. **Scanner module not found**
   ```r
   # Add parent directory to Python path
   library(reticulate)
   py_run_string("import sys; sys.path.insert(0, '..')")
   init_nonbscanner()
   ```

3. **ggplot2 visualizations not working**
   ```r
   # Install ggplot2
   install.packages("ggplot2")
   ```

## 📖 Additional Resources

- **Main Repository**: https://github.com/VRYella/NonBScanner
- **Python Documentation**: See main README.md
- **Web Interface**: Run `streamlit run ../app.py`

## 👨‍💻 Author

**Dr. Venkata Rajesh Yella**
- Email: yvrajesh_bt@kluniversity.in
- GitHub: [@VRYella](https://github.com/VRYella)

## 📜 License

MIT License - see LICENSE file for details

## 🙏 Citation

If you use NonBScanner in your research, please cite:

```
NonBScanner: Comprehensive Detection and Analysis of Non-B DNA Motifs
Dr. Venkata Rajesh Yella
GitHub: https://github.com/VRYella/NonBScanner
```

---

**Version**: 1.0  
**Last Updated**: 2024
