# NonBScanner CSV Generator

**Shell Script for Batch CSV Generation from FASTA Files**

---

## 📋 Overview

The `generate_csv_output.sh` script provides a convenient command-line interface for generating CSV output files for all Non-B DNA motif classes and subclasses detected by NonBScanner.

## 🎯 Features

- **Batch Processing**: Analyze multiple sequences from FASTA files
- **Comprehensive CSV Output**: Generate separate CSV files for each motif class
- **Summary Statistics**: Automatic generation of summary reports
- **Easy to Use**: Simple command-line interface with helpful options
- **Flexible Output**: Customizable output directory and file prefixes

## 📦 Prerequisites

- **Bash** shell (Linux, macOS, WSL on Windows)
- **Python 3.8+** with NonBScanner installed
- **Required Python packages**: numpy, pandas, and other NonBScanner dependencies

## 🚀 Quick Start

### Basic Usage

```bash
./generate_csv_output.sh -i sequences.fasta
```

This will:
1. Read sequences from `sequences.fasta`
2. Detect all Non-B DNA motifs
3. Generate CSV files in the `output/` directory

### Specify Output Directory

```bash
./generate_csv_output.sh -i sequences.fasta -o results
```

### Custom Output Prefix

```bash
./generate_csv_output.sh -i sequences.fasta -o results -p my_analysis
```

### Verbose Mode

```bash
./generate_csv_output.sh -i sequences.fasta -v
```

## 📚 Command-Line Options

### Required Arguments

| Option | Description |
|--------|-------------|
| `-i, --input FILE` | Input FASTA file containing DNA sequences |

### Optional Arguments

| Option | Description | Default |
|--------|-------------|---------|
| `-o, --output DIR` | Output directory | `output` |
| `-p, --prefix PREFIX` | Output file prefix | `nonbscanner` |
| `-v, --verbose` | Enable verbose output | disabled |
| `-h, --help` | Display help message | - |

## 📊 Output Files

The script generates the following CSV files:

### Main Output Files

1. **`PREFIX_all_motifs.csv`**
   - Contains all detected motifs from all sequences
   - Includes complete motif information

2. **`PREFIX_summary.csv`**
   - Summary statistics
   - Total motifs, classes, subclasses
   - Counts per class

3. **`PREFIX_by_class.csv`**
   - Motif counts grouped by class

4. **`PREFIX_by_subclass.csv`**
   - Motif counts grouped by class and subclass

### Individual Class Files

Separate CSV files for each detected motif class:

- `PREFIX_Curved_DNA.csv`
- `PREFIX_Slipped_DNA.csv`
- `PREFIX_Cruciform.csv`
- `PREFIX_R_Loop.csv`
- `PREFIX_Triplex.csv`
- `PREFIX_G_Quadruplex.csv`
- `PREFIX_i_Motif.csv`
- `PREFIX_Z_DNA.csv`
- `PREFIX_A_philic_DNA.csv`
- `PREFIX_Hybrid.csv`
- `PREFIX_Non_B_DNA_Clusters.csv`

## 📝 CSV Output Format

Each CSV file contains the following columns:

| Column | Description |
|--------|-------------|
| `ID` | Unique motif identifier |
| `Sequence_Name` | Name of the analyzed sequence |
| `Source` | Source information |
| `Class` | Motif class |
| `Subclass` | Motif subclass |
| `Pattern_ID` | Pattern/annotation ID |
| `Start` | Start position (1-based) |
| `End` | End position |
| `Length` | Motif length in base pairs |
| `Sequence` | Actual sequence of the motif |
| `Method` | Detection method |
| `Score` | Motif score |
| ... | Additional class-specific fields |

## 💡 Examples

### Example 1: Basic Analysis

```bash
./generate_csv_output.sh -i example_all_motifs.fasta
```

Output:
```
Found 52 motifs
Generated files:
  output/nonbscanner_all_motifs.csv
  output/nonbscanner_summary.csv
  output/nonbscanner_by_class.csv
  output/nonbscanner_by_subclass.csv
  output/nonbscanner_G_Quadruplex.csv
  output/nonbscanner_Z_DNA.csv
  ...
```

### Example 2: Custom Output Location

```bash
./generate_csv_output.sh -i sequences.fasta -o my_results -p genome1
```

Output files:
- `my_results/genome1_all_motifs.csv`
- `my_results/genome1_summary.csv`
- etc.

### Example 3: Verbose Analysis

```bash
./generate_csv_output.sh -i sequences.fasta -v
```

Output shows:
- Sequences being processed
- Motif detection progress
- Detailed file generation info

## 🔬 Motif Classes Detected

The script detects **11 major classes** with **22+ subclasses**:

1. **Curved DNA** (2 subclasses)
2. **Slipped DNA** (2 subclasses)
3. **Cruciform** (1 subclass)
4. **R-Loop** (1 subclass)
5. **Triplex** (2 subclasses)
6. **G-Quadruplex** (7 subclasses)
7. **i-Motif** (3 subclasses)
8. **Z-DNA** (2 subclasses)
9. **A-philic DNA** (1 subclass)
10. **Hybrid** (dynamic overlaps)
11. **Non-B DNA Clusters** (dynamic clusters)

## 🔧 Troubleshooting

### Permission Denied

```bash
chmod +x generate_csv_output.sh
```

### Python Not Found

Ensure Python 3 is installed and in your PATH:
```bash
which python3
```

### Module Not Found

Install NonBScanner dependencies:
```bash
pip install -r requirements.txt
```

### Script Fails

Run with verbose mode to see detailed error messages:
```bash
./generate_csv_output.sh -i sequences.fasta -v
```

## 📖 Input File Format

The script accepts standard FASTA format:

```
>sequence1
GGGTTAGGGTTAGGGTTAGGGAAAAATTTTAAAAATTTT
>sequence2
CGCGCGCGCGCGCGCACACACACACACA
>sequence3
CCCCTAACCCCTAACCCCTAACCC
```

## 🎓 Advanced Usage

### Process Multiple Files

```bash
for file in *.fasta; do
    ./generate_csv_output.sh -i "$file" -o results -p "${file%.fasta}"
done
```

### Integrate with Pipeline

```bash
# Download sequences
curl -o sequences.fasta "http://example.com/sequences.fasta"

# Run analysis
./generate_csv_output.sh -i sequences.fasta -o results -v

# Post-process results
python analyze_results.py results/nonbscanner_all_motifs.csv
```

## 📊 Example Output Summary

```
╔════════════════════════════════════════════════════════════════╗
║       NonBScanner CSV Generator                                ║
║       Non-B DNA Motif Detection and CSV Export                 ║
╚════════════════════════════════════════════════════════════════╝

[INFO] Input file: sequences.fasta
[INFO] Output directory: output
[INFO] Output prefix: nonbscanner

============================================================
ANALYSIS SUMMARY
============================================================
Total sequences analyzed: 1
Total motifs detected: 52
Unique classes: 6
Unique subclasses: 10

Motifs by class:
  G-Quadruplex                        12 motifs
  Z-DNA                                8 motifs
  Curved_DNA                           7 motifs
  Slipped_DNA                          6 motifs
  i-Motif                              5 motifs
  Cruciform                            4 motifs
============================================================

[SUCCESS] Analysis complete!
[SUCCESS] Output files saved to: output
```

## 🤝 Integration with Other Tools

### Load Results in Python

```python
import pandas as pd

# Load all motifs
motifs = pd.read_csv("output/nonbscanner_all_motifs.csv")

# Load specific class
g4_motifs = pd.read_csv("output/nonbscanner_G_Quadruplex.csv")
```

### Load Results in R

```r
# Load all motifs
motifs <- read.csv("output/nonbscanner_all_motifs.csv")

# Load summary
summary <- read.csv("output/nonbscanner_summary.csv")
```

### Export to BED Format

```bash
# After running the script, convert CSV to BED
awk -F',' 'NR>1 {print $2"\t"$7"\t"$8"\t"$4"_"$5}' \
    output/nonbscanner_all_motifs.csv > motifs.bed
```

## 📚 Additional Resources

- **Main Documentation**: See main README.md
- **Web Interface**: Run `streamlit run app.py`
- **Jupyter Notebook**: See NonBScanner_Local.ipynb
- **R Interface**: See R_tool/README.md

## 👨‍💻 Author

**Dr. Venkata Rajesh Yella**
- Email: yvrajesh_bt@kluniversity.in
- GitHub: [@VRYella](https://github.com/VRYella)

## 📜 License

MIT License - see LICENSE file for details

---

**Version**: 1.0  
**Last Updated**: 2024
