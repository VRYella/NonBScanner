# NonBScanner Local Jupyter Notebook

**Interactive Notebook for Non-B DNA Motif Detection**

---

## 📋 Overview

The `NonBScanner_Local.ipynb` Jupyter notebook provides an interactive, user-friendly interface for detecting and analyzing Non-B DNA motifs directly in your local environment. This notebook is perfect for researchers who want to explore their sequences interactively with immediate visualization and analysis capabilities.

## 🎯 Features

- **Interactive Analysis**: Run analyses cell-by-cell with immediate results
- **Built-in Visualizations**: Publication-ready plots and charts
- **Multiple Input Methods**: Paste sequences, load FASTA files, or batch process
- **Comprehensive Examples**: Pre-configured examples for all motif classes
- **Export Capabilities**: Save results to CSV, Excel, and other formats
- **Educational**: Includes detailed explanations and usage examples

## 📦 Prerequisites

### Software Requirements

- **Python 3.8+**
- **Jupyter Notebook** or **JupyterLab**

### Required Python Packages

```bash
pip install numpy pandas matplotlib seaborn plotly biopython scipy scikit-learn openpyxl
```

Or install from requirements:

```bash
pip install -r requirements.txt
```

### Installing Jupyter

If you don't have Jupyter installed:

```bash
pip install jupyter
# or
pip install jupyterlab
```

## 🚀 Quick Start

### Launch the Notebook

#### Option 1: Jupyter Notebook
```bash
cd NonBScanner
jupyter notebook NonBScanner_Local.ipynb
```

#### Option 2: JupyterLab
```bash
cd NonBScanner
jupyter lab NonBScanner_Local.ipynb
```

#### Option 3: VS Code
1. Open VS Code
2. Install the Jupyter extension
3. Open `NonBScanner_Local.ipynb`
4. Select Python kernel

### First Run

1. **Execute the first cell** to install dependencies (if needed)
2. **Run the import cell** to load NonBScanner
3. **Check available motif classes** to see what can be detected
4. **Try the example analysis** to verify everything works

## 📚 Notebook Structure

The notebook is organized into the following sections:

### 1. Installation
- Dependency installation commands
- Environment setup

### 2. Quick Start
- Import NonBScanner modules
- Initialize the tool
- Check available motif classes

### 3. Single Sequence Analysis
- Analyze individual sequences
- View detailed results
- Understand output format

### 4. FASTA File Analysis
- Read FASTA files
- Batch process multiple sequences
- Handle large datasets

### 5. Visualization
- Motif distribution plots
- Position maps
- Score distributions
- Class comparisons

### 6. Export Results
- Export to CSV
- Export to Excel (multi-sheet)
- Generate summary statistics
- Create custom reports

### 7. Advanced Usage
- Batch processing workflows
- Filter results by class
- Custom analysis pipelines

## 💡 Usage Examples

### Example 1: Quick Analysis

```python
# Import and initialize
from scanner import analyze_sequence
import pandas as pd

# Analyze a sequence
sequence = "GGGTTAGGGTTAGGGTTAGGG"
motifs = analyze_sequence(sequence, "test")

# View results
print(f"Found {len(motifs)} motifs")
for m in motifs:
    print(f"{m['Class']}: {m['Start']}-{m['End']}")
```

### Example 2: FASTA File Analysis

```python
# Define FASTA reader
def read_fasta(filename):
    sequences = {}
    current_name = None
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        if current_name:
            sequences[current_name] = ''.join(current_seq)
    
    return sequences

# Read and analyze
sequences = read_fasta("sequences.fasta")
results = analyze_multiple_sequences(sequences)
```

### Example 3: Visualization

```python
import matplotlib.pyplot as plt

# Count motifs by class
class_counts = {}
for m in motifs:
    cls = m['Class']
    class_counts[cls] = class_counts.get(cls, 0) + 1

# Create bar plot
plt.figure(figsize=(12, 6))
plt.bar(class_counts.keys(), class_counts.values())
plt.xlabel('Motif Class')
plt.ylabel('Count')
plt.title('Non-B DNA Motif Distribution')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()
```

### Example 4: Export to CSV

```python
from scanner import export_results_to_dataframe

# Convert to DataFrame
df = export_results_to_dataframe(motifs)

# Save to CSV
df.to_csv("results.csv", index=False)

# Save by class
for cls in df['Class'].unique():
    class_df = df[df['Class'] == cls]
    filename = f"{cls.replace(' ', '_')}.csv"
    class_df.to_csv(filename, index=False)
```

## 🔬 Motif Classes Detected

The notebook can detect **11 major classes** with **22+ subclasses**:

| Class | Subclasses | Example Detection |
|-------|------------|-------------------|
| Curved DNA | 2 | A-tract runs |
| Slipped DNA | 2 | Tandem repeats |
| Cruciform | 1 | Palindromes |
| R-Loop | 1 | GC-rich regions |
| Triplex | 2 | Purine/pyrimidine tracts |
| G-Quadruplex | 7 | G-rich structures |
| i-Motif | 3 | C-rich structures |
| Z-DNA | 2 | Alternating purines/pyrimidines |
| A-philic | 1 | Poly-A tracts |
| Hybrid | Dynamic | Overlapping motifs |
| Clusters | Dynamic | High-density regions |

## 📊 Output Format

Results are provided as Python dictionaries/DataFrames with fields:

- `ID`: Unique identifier
- `Sequence_Name`: Sequence identifier
- `Class`: Motif class
- `Subclass`: Specific subclass
- `Start`: Start position (1-based)
- `End`: End position
- `Length`: Length in base pairs
- `Sequence`: Actual sequence
- `Score`: Detection score
- Additional class-specific fields

## 🎨 Visualization Gallery

The notebook includes examples for:

1. **Bar Charts**: Motif distribution by class
2. **Position Maps**: Genome-wide motif locations
3. **Violin Plots**: Score distributions
4. **Heatmaps**: Density analysis
5. **Custom Plots**: User-defined visualizations

## 💾 Export Formats

### CSV Export
```python
df.to_csv("results.csv", index=False)
```

### Excel Export (Multi-Sheet)
```python
with pd.ExcelWriter("results.xlsx") as writer:
    df.to_excel(writer, sheet_name='All_Motifs', index=False)
    # Add sheets by class
    for cls in df['Class'].unique():
        class_df = df[df['Class'] == cls]
        class_df.to_excel(writer, sheet_name=cls[:31], index=False)
```

### JSON Export
```python
import json
with open("results.json", "w") as f:
    json.dump(motifs, f, indent=2)
```

## 🔧 Troubleshooting

### Kernel Not Found
```bash
# Install ipykernel
pip install ipykernel
python -m ipykernel install --user
```

### Module Import Errors
```bash
# Install in the same environment as Jupyter
pip install numpy pandas matplotlib seaborn
```

### Large Files
For very large FASTA files:
```python
# Process in chunks
sequences = read_fasta("large_file.fasta")
batch_size = 10

for i in range(0, len(sequences), batch_size):
    batch = dict(list(sequences.items())[i:i+batch_size])
    results = analyze_multiple_sequences(batch)
    # Process results
```

### Memory Issues
```python
# Use multiprocessing for large datasets
results = analyze_multiple_sequences(sequences, 
                                    use_multiprocessing=True)
```

## 🎓 Learning Resources

### Tutorials in the Notebook

1. **Basic Tutorial**: Step-by-step guide for beginners
2. **Advanced Tutorial**: Complex analysis workflows
3. **Visualization Tutorial**: Create publication-quality plots
4. **Export Tutorial**: Save and share results

### Example Datasets

The notebook includes examples using:
- `example_all_motifs.fasta`: Comprehensive test sequence
- Built-in test sequences
- Simulated data for practice

## 🔬 Research Applications

### Use Cases

1. **Genome Analysis**: Scan entire genomes for motifs
2. **Gene Studies**: Analyze promoter regions
3. **Structural Biology**: Identify non-canonical DNA
4. **Evolution Studies**: Compare motif distributions
5. **Drug Discovery**: Target identification

### Integration with Pipelines

```python
# Example: Integration with BioPython
from Bio import SeqIO

# Read GenBank file
for record in SeqIO.parse("genome.gb", "genbank"):
    sequence = str(record.seq)
    motifs = analyze_sequence(sequence, record.id)
    # Process motifs
```

## 📖 Additional Features

### Interactive Widgets (Optional)

For enhanced interactivity, install ipywidgets:

```bash
pip install ipywidgets
jupyter nbextension enable --py widgetsnbextension
```

Then use interactive controls:

```python
from ipywidgets import interact, widgets

@interact(
    motif_class=widgets.Dropdown(
        options=['All', 'G-Quadruplex', 'Z-DNA', 'Cruciform'],
        value='All'
    )
)
def filter_motifs(motif_class):
    if motif_class == 'All':
        display(df)
    else:
        display(df[df['Class'] == motif_class])
```

## 🚀 Performance Tips

1. **Use `use_multiprocessing=True`** for multiple sequences
2. **Process large files in batches**
3. **Clear output cells** to save memory
4. **Restart kernel** if needed

## 🤝 Contributing

Found a bug or have a suggestion? Please:
1. Open an issue on GitHub
2. Submit a pull request
3. Contact the author

## 📚 Related Tools

- **Shell Script**: `generate_csv_output.sh` for batch processing
- **R Interface**: See `R_tool/` directory
- **Web App**: Run `streamlit run app.py`

## 👨‍💻 Author

**Dr. Venkata Rajesh Yella**
- Email: yvrajesh_bt@kluniversity.in
- GitHub: [@VRYella](https://github.com/VRYella)

## 📜 License

MIT License - see LICENSE file for details

## 🙏 Citation

If you use this notebook in your research:

```
NonBScanner: Comprehensive Detection and Analysis of Non-B DNA Motifs
Dr. Venkata Rajesh Yella
GitHub: https://github.com/VRYella/NonBScanner
```

---

**Version**: 1.0  
**Last Updated**: 2024

**Happy Analyzing! 🧬**
