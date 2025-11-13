# Quick Start Guide for New NonBScanner Tools

This guide helps you quickly get started with the three new tools added to NonBScanner.

---

## 🎯 Choose Your Tool

| Tool | Best For | Quick Start |
|------|----------|-------------|
| **Jupyter Notebook** | Interactive exploration, learning, visualization | `jupyter notebook NonBScanner_Local.ipynb` |
| **Shell Script** | Batch processing, automation, pipelines | `./generate_csv_output.sh -i file.fasta` |
| **R Interface** | R workflows, statistical analysis, R plots | `source("R_tool/nonbscanner.R")` |

---

## 📓 Jupyter Notebook (Interactive)

### Quick Start

```bash
# Install Jupyter (if needed)
pip install jupyter

# Launch notebook
jupyter notebook NonBScanner_Local.ipynb
```

### First Steps in Notebook

1. **Run cell 1**: Install dependencies (if needed)
2. **Run cell 2**: Import modules
3. **Run cell 3**: Check available motif classes
4. **Run cell 4**: Analyze example sequence

### Common Tasks

**Analyze a sequence:**
```python
from scanner import analyze_sequence
motifs = analyze_sequence("GGGTTAGGGTTAGGGTTAGGG", "test")
```

**Read FASTA file:**
```python
sequences = read_fasta("sequences.fasta")
results = analyze_multiple_sequences(sequences)
```

**Export results:**
```python
from scanner import export_results_to_dataframe
df = export_results_to_dataframe(motifs)
df.to_csv("results.csv", index=False)
```

📖 **Full Documentation**: [JUPYTER_NOTEBOOK_LOCAL_README.md](JUPYTER_NOTEBOOK_LOCAL_README.md)

---

## 🔧 Shell Script (Command Line)

### Quick Start

```bash
# Make executable (first time only)
chmod +x generate_csv_output.sh

# Run analysis
./generate_csv_output.sh -i sequences.fasta
```

### Output Location

Results saved in `output/` directory:
- `nonbscanner_all_motifs.csv` - All detected motifs
- `nonbscanner_summary.csv` - Summary statistics
- `nonbscanner_G_Quadruplex.csv` - Individual class files
- ... (one CSV per detected class)

### Common Options

**Custom output directory:**
```bash
./generate_csv_output.sh -i sequences.fasta -o my_results
```

**Custom file prefix:**
```bash
./generate_csv_output.sh -i sequences.fasta -p genome1
```

**Verbose mode:**
```bash
./generate_csv_output.sh -i sequences.fasta -v
```

**Help:**
```bash
./generate_csv_output.sh -h
```

📖 **Full Documentation**: [CSV_GENERATOR_README.md](CSV_GENERATOR_README.md)

---

## 📊 R Interface (R Programming)

### Quick Start

```bash
# Navigate to R tool directory
cd R_tool
```

```r
# In R console

# Step 1: Install dependencies (first time only)
source("install.R")

# Step 2: Load NonBScanner
source("nonbscanner.R")

# Step 3: Initialize
init_nonbscanner()

# Step 4: Analyze
motifs <- analyze_sequence("GGGTTAGGGTTAGGGTTAGGG", "test")
```

### Common Tasks

**Analyze FASTA file:**
```r
motifs <- analyze_fasta("../sequences.fasta")
```

**Visualize results:**
```r
library(ggplot2)
p <- plot_motif_distribution(motifs)
print(p)
ggsave("distribution.png", p)
```

**Export to CSV:**
```r
write.csv(motifs, "results.csv", row.names = FALSE)
```

**Get summary:**
```r
summary <- get_summary_stats(motifs)
print(summary)
```

### Run Examples

```r
# Run all examples
source("example_usage.R")
```

📖 **Full Documentation**: [R_tool/README.md](R_tool/README.md)

---

## 🔬 Example Workflows

### Workflow 1: Quick Exploration (Jupyter)

1. Open `NonBScanner_Local.ipynb`
2. Paste your sequence in a cell
3. Run analysis cell
4. View inline visualizations
5. Export if needed

**Time**: 5 minutes

---

### Workflow 2: Batch Processing (Shell Script)

1. Prepare FASTA file with multiple sequences
2. Run: `./generate_csv_output.sh -i input.fasta -v`
3. Check `output/` directory for results
4. Each motif class gets its own CSV file

**Time**: 1 minute + processing time

---

### Workflow 3: R Analysis Pipeline (R)

```r
# 1. Load tools
source("R_tool/nonbscanner.R")
init_nonbscanner()

# 2. Read data
motifs <- analyze_fasta("sequences.fasta")

# 3. Analyze
summary <- get_summary_stats(motifs)

# 4. Visualize
p1 <- plot_motif_distribution(motifs)
p2 <- plot_motif_positions(motifs)

# 5. Export
write.csv(motifs, "results.csv", row.names = FALSE)
ggsave("plots.pdf", gridExtra::grid.arrange(p1, p2))
```

**Time**: 10 minutes

---

## 🎓 Learning Path

### Beginner
1. Start with **Jupyter Notebook**
2. Work through the example cells
3. Try your own sequences
4. Learn about different motif classes

### Intermediate
1. Use **Shell Script** for batch processing
2. Process multiple FASTA files
3. Understand CSV output format
4. Integrate with other tools

### Advanced
1. Use **R Interface** for custom analyses
2. Create publication-quality plots
3. Build automated pipelines
4. Integrate with statistical workflows

---

## 📋 Checklist: First Time Setup

### For All Tools
- [ ] Clone NonBScanner repository
- [ ] Install Python 3.8+
- [ ] Install dependencies: `pip install -r requirements.txt`

### For Jupyter
- [ ] Install Jupyter: `pip install jupyter`
- [ ] Test: `jupyter notebook NonBScanner_Local.ipynb`

### For Shell Script
- [ ] Make executable: `chmod +x generate_csv_output.sh`
- [ ] Test help: `./generate_csv_output.sh -h`

### For R
- [ ] Install R 4.0+
- [ ] Install reticulate: `install.packages("reticulate")`
- [ ] Run install script: `source("R_tool/install.R")`

---

## 🆘 Common Issues

### Issue: "Module not found"
**Solution**: Make sure you're in the NonBScanner directory and dependencies are installed

### Issue: "Permission denied" (Shell script)
**Solution**: Run `chmod +x generate_csv_output.sh`

### Issue: "Python not found" (R)
**Solution**: In R, use `reticulate::use_python("/path/to/python3")`

### Issue: "No motifs detected"
**Solution**: Check sequence format (DNA only, ATGC characters)

---

## 📚 Documentation Index

| File | Description |
|------|-------------|
| `JUPYTER_NOTEBOOK_LOCAL_README.md` | Complete Jupyter guide |
| `CSV_GENERATOR_README.md` | Shell script reference |
| `R_tool/README.md` | R interface documentation |
| `NEW_TOOLS_SUMMARY.md` | Implementation summary |
| This file | Quick start guide |

---

## 🎯 Next Steps

After getting started:

1. **Read full documentation** for your preferred tool
2. **Try example datasets** provided in the repository
3. **Customize** analysis parameters
4. **Integrate** with your workflow
5. **Share** results with collaborators

---

## 💡 Tips

- **Jupyter**: Use Shift+Enter to run cells
- **Shell**: Use `-v` flag to see detailed progress
- **R**: Check `example_usage.R` for comprehensive examples
- **All tools**: Start with small test sequences

---

## 🤝 Getting Help

- Check tool-specific documentation
- Review example files
- See main README.md
- Contact: Dr. Venkata Rajesh Yella (yvrajesh_bt@kluniversity.in)

---

**Happy Analyzing! 🧬**
