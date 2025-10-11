# NBDScanner Quick Start Guide

## Welcome to NBDScanner! 🧬

This guide will get you up and running with NBDScanner in just a few minutes. NBDScanner is a powerful tool for detecting Non-B DNA structures in genomic sequences.

---

## Table of Contents

1. [Installation](#installation)
2. [Quick Launch](#quick-launch)
3. [First Analysis](#first-analysis)
4. [Understanding Results](#understanding-results)
5. [Common Workflows](#common-workflows)
6. [Troubleshooting](#troubleshooting)

---

## Installation

### Prerequisites

- **Python 3.8 or higher**
- **pip** (Python package manager)
- **Git** (for cloning the repository)

### Step 1: Clone the Repository

```bash
git clone https://github.com/VRYella/NonBScanner.git
cd NonBScanner
```

### Step 2: Create Virtual Environment (Recommended)

**On Linux/Mac:**
```bash
python3 -m venv venv
source venv/bin/activate
```

**On Windows:**
```bash
python -m venv venv
venv\Scripts\activate
```

### Step 3: Install Dependencies

```bash
pip install -r requirements.txt
```

**Expected Installation Time:** 2-3 minutes

### Step 4: Verify Installation

```bash
python -c "from utils.nbdscanner import analyze_sequence; print('✅ Installation successful!')"
```

If you see `✅ Installation successful!`, you're ready to go!

---

## Quick Launch

### Start the Web Interface

```bash
streamlit run app.py
```

**What happens:**
1. Streamlit server starts
2. Browser opens automatically at `http://localhost:8501`
3. NBDScanner interface loads

**First-Time Launch:**
- May take 10-15 seconds to load all modules
- You'll see a welcome screen with tabs

---

## First Analysis

### Option 1: Use Demo Data (Fastest Way to Start)

1. **Navigate to "Upload & Analyze" tab**
2. **Click "Load Demo Sequence" button**
3. **Click "Analyze Sequences" button**
4. **Wait 2-5 seconds**
5. **View results in "Results" tab**

### Option 2: Paste Your Own Sequence

1. **Go to "Upload & Analyze" tab**
2. **Click "Paste Sequence(s)" option**
3. **Paste your DNA sequence:**
   ```
   >MySequence
   AAAATTTTGGGGCCCCAAGGTTCCAAAATTTTGGGGCCCCATGCATGC
   ```
4. **Click "Analyze Sequences"**
5. **View results**

### Option 3: Upload FASTA File

1. **Go to "Upload & Analyze" tab**
2. **Click "Upload FASTA File" option**
3. **Select your `.fasta` or `.fa` file**
4. **Click "Analyze Sequences"**
5. **View results**

---

## Understanding Results

### Results Tab Overview

When analysis completes, you'll see:

#### 1. Summary Statistics (Top Banner)
```
┌─────────────────────────────────────────────────┐
│  Coverage: 15.2%  |  Density: 8.5 motifs/kb    │
│  Total Motifs: 42 |  Classes Detected: 6       │
└─────────────────────────────────────────────────┘
```

#### 2. Visualizations

**A. Motif Distribution Chart**
- Bar chart showing count of each motif class
- Color-coded by class
- Helps identify dominant motif types

**B. Coverage Map**
- Heatmap showing motif locations
- Genomic position on X-axis
- Motif types on Y-axis

**C. Score Distribution**
- Histogram of motif quality scores
- Higher scores = more confident predictions

**D. Nested Pie Chart**
- Outer ring: Major classes
- Inner ring: Subclasses
- Interactive hover tooltips

#### 3. Detailed Motif Table

| Sequence | Start | End | Class | Subclass | Score | Sequence |
|----------|-------|-----|-------|----------|-------|----------|
| seq1 | 0 | 10 | Curved DNA | A-tract | 0.85 | AAAAAAAAAA |
| seq1 | 20 | 35 | G-Quadruplex | Canonical | 0.92 | GGGATTGGGATTGGG |

**Columns Explained:**
- **Start/End:** Genomic coordinates (0-based)
- **Class:** Major motif category
- **Subclass:** Specific motif type
- **Score:** Confidence (0-1 scale, higher is better)
- **Sequence:** Actual DNA sequence of motif

---

## Common Workflows

### Workflow 1: Analyze Gene Promoter

**Goal:** Find regulatory elements in a gene promoter

```python
# If using Python API instead of web interface
from utils import analyze_sequence, parse_fasta

# Load promoter sequence
with open('TP53_promoter.fasta', 'r') as f:
    sequences = parse_fasta(f)

# Analyze
motifs = analyze_sequence(
    sequence=sequences[0]['seq'],
    sequence_name='TP53_promoter'
)

# Focus on G-quadruplexes (common in promoters)
g4_motifs = [m for m in motifs if m['Class'] == 'G-Quadruplex']
print(f"Found {len(g4_motifs)} G4 motifs")
```

### Workflow 2: Compare Multiple Genes

**Goal:** Batch analyze multiple sequences

1. **Prepare FASTA file with multiple sequences:**
   ```
   >Gene1_Promoter
   ATGCATGCAAAATTTTGGGGCCCCATGC
   >Gene2_Promoter
   GGGGTTTTAAAACCCCGGGGTTTT
   >Gene3_Promoter
   CCCCAAAATTTTGGGGCCCCATGC
   ```

2. **Upload to NBDScanner**
3. **Results show comparison across all sequences**
4. **Export to Excel for further analysis**

### Workflow 3: Export Results for Publication

**Goal:** Generate publication-quality figures

1. **Run analysis**
2. **Go to "Results" tab**
3. **Scroll to visualizations**
4. **Right-click on each plot → "Save Image"**
5. **Or click "Download All Plots" button**
6. **Plots saved as high-resolution PNG files**

### Workflow 4: Focus on Specific Motif Classes

**Goal:** Only detect G-Quadruplexes and i-Motifs

```python
from utils.modular_scanner import ModularMotifDetector

detector = ModularMotifDetector()

# Only use specific detectors
results = detector.detectors['g_quadruplex'].detect(sequence, "seq1")
results.extend(detector.detectors['i_motif'].detect(sequence, "seq1"))

print(f"Found {len(results)} G4/i-Motif structures")
```

---

## Downloading and Exporting

### Available Export Formats

#### 1. CSV Format (Universal)
```
Sequence,Start,End,Length,Class,Subclass,Score
seq1,0,10,10,Curved DNA,A-tract,0.85
seq1,20,35,15,G-Quadruplex,Canonical,0.92
```

**Use for:** Excel, statistical analysis, databases

#### 2. BED Format (Genomic)
```
chr1    100    110    Curved_DNA_A-tract    850    +
chr1    500    515    G4_Canonical          920    +
```

**Use for:** UCSC Genome Browser, IGV, bedtools

#### 3. JSON Format (Structured)
```json
{
  "sequence_name": "seq1",
  "motifs": [
    {"class": "G-Quadruplex", "start": 20, "score": 0.92}
  ]
}
```

**Use for:** APIs, web applications, databases

#### 4. Excel Format (Multi-sheet)
- Sheet 1: All motifs
- Sheet 2: Statistics summary
- Sheet 3: Class breakdown

**Use for:** Reports, presentations, sharing

---

## Interpreting Motif Classes

### Quick Reference

| Motif Class | What It Means | Biological Importance |
|-------------|---------------|----------------------|
| **Curved DNA** | DNA bending (A/T-tracts) | Nucleosome positioning, transcription |
| **G-Quadruplex** | Four-stranded G-rich structure | Telomeres, gene regulation |
| **i-Motif** | Four-stranded C-rich structure | pH sensing, gene regulation |
| **Z-DNA** | Left-handed helix | Transcription activation |
| **R-Loop** | RNA-DNA hybrid | Transcription, genome instability |
| **Triplex** | Three-stranded DNA | Gene regulation |
| **Cruciform** | Four-way junction | Recombination, chromosomal fragility |
| **Slipped DNA** | Tandem repeat expansion | Disease-causing mutations |
| **A-philic** | A-rich protein binding | Transcription factor sites |
| **Hybrid** | Multiple overlapping motifs | Complex regulation |
| **Cluster** | High motif density region | Fragile sites, hotspots |

---

## Troubleshooting

### Issue 1: "ModuleNotFoundError"

**Error:**
```
ModuleNotFoundError: No module named 'pandas'
```

**Solution:**
```bash
pip install -r requirements.txt
```

### Issue 2: "Streamlit won't start"

**Error:**
```
streamlit: command not found
```

**Solution:**
```bash
# Reinstall streamlit
pip install streamlit

# Or use python -m
python -m streamlit run app.py
```

### Issue 3: "Analysis is very slow"

**Problem:** Taking more than 1 minute for small sequences

**Solutions:**
1. **Check sequence size:**
   - Cruciform detector slow for >1KB sequences
   - Slipped DNA slow for >50KB sequences

2. **Use fast detectors only:**
   - Curved DNA, G-Quadruplex, i-Motif, Z-DNA are fast

3. **Install Hyperscan (optional):**
   ```bash
   pip install hyperscan
   ```

### Issue 4: "Invalid sequence error"

**Error:**
```
Invalid DNA characters found
```

**Solutions:**
1. **Check for non-ATGC characters**
   - Remove spaces, numbers, special characters
   - Keep only A, T, G, C (uppercase or lowercase)

2. **Use sequence cleaning:**
   ```python
   from utils import validate_sequence
   
   clean_seq = validate_sequence(raw_sequence)
   ```

### Issue 5: "No motifs detected"

**Problem:** Analysis completes but finds 0 motifs

**Possible Reasons:**
1. **Sequence too short** (< 20 bp)
2. **No recognizable patterns**
3. **Score threshold too high**

**Solutions:**
1. Use longer sequences (>100 bp recommended)
2. Try demo data to verify installation
3. Check sequence for homopolymers (AAAA, GGGG, etc.)

---

## Performance Tips

### For Best Performance:

1. **Sequence Size:**
   - Optimal: 1 KB - 100 KB
   - Avoid: > 1 MB without chunking

2. **Batch Processing:**
   - Process multiple small sequences rather than one huge sequence
   - Use FASTA with multiple entries

3. **Memory:**
   - Close other applications
   - Process in batches if analyzing > 100 sequences

4. **Export:**
   - Export to CSV for large datasets (faster than Excel)
   - Download plots individually if analyzing many sequences

---

## Next Steps

### Learn More:

1. **Read Full Documentation:** `TOOL_DOCUMENTATION.md`
2. **View Flowcharts:** `VISUAL_FLOWCHARTS.md`
3. **Check Examples:** Browse `examples/` directory (if available)
4. **Watch Tutorials:** Video guides (coming soon)

### Advanced Usage:

1. **Python API:**
   ```python
   from utils import analyze_sequence, export_to_csv
   
   results = analyze_sequence(sequence, "my_seq")
   export_to_csv(results, "output.csv")
   ```

2. **Custom Configuration:**
   - Modify detector parameters
   - Adjust score thresholds
   - Enable/disable specific detectors

3. **Integration:**
   - Use as part of genomics pipeline
   - Integrate with other tools
   - Automate batch processing

---

## Getting Help

### Support Resources:

1. **Documentation:** Check `README.md` and `TOOL_DOCUMENTATION.md`
2. **GitHub Issues:** Report bugs or request features
3. **Email:** raazbiochem@gmail.com
4. **Discussion Forum:** GitHub Discussions (coming soon)

### Before Asking for Help:

1. ✅ Check this Quick Start Guide
2. ✅ Review Troubleshooting section
3. ✅ Try demo data to verify installation
4. ✅ Check GitHub Issues for similar problems
5. ✅ Prepare minimal example showing the issue

---

## Example: Complete Analysis Workflow

### Step-by-Step Tutorial

Let's analyze a small test sequence from start to finish:

#### 1. Start NBDScanner
```bash
streamlit run app.py
```

#### 2. Go to "Upload & Analyze" Tab

#### 3. Paste This Test Sequence
```
>test_sequence
AAAATTTTGGGGCCCCAAAATTTTGGGGCCCCATGCGGGGTTTTGGGGTTTTGGGGTTTTGGGGCCCCAAAACCCCATGC
```

This sequence contains:
- A-tracts (Curved DNA)
- G-quadruplexes
- i-Motifs

#### 4. Click "Analyze Sequences"

#### 5. View Results in "Results" Tab

You should see:
- **~6-8 motifs detected**
- **Classes:** Curved DNA, G-Quadruplex, i-Motif
- **Coverage:** ~40-50%
- **Density:** ~75-100 motifs/kb

#### 6. Export Results

- Go to "Download" tab
- Click "Download CSV"
- Open in Excel to view

**Congratulations!** You've completed your first NBDScanner analysis! 🎉

---

## Summary Checklist

Before you start analyzing:

- [ ] Python 3.8+ installed
- [ ] Repository cloned
- [ ] Virtual environment created
- [ ] Dependencies installed (`pip install -r requirements.txt`)
- [ ] Installation verified
- [ ] Streamlit launches successfully
- [ ] Demo data works

Ready to analyze your sequences? **Happy scanning!** 🧬🔬

---

**Document Version:** 1.0  
**Last Updated:** October 2024  
**Author:** Dr. Venkata Rajesh Yella

