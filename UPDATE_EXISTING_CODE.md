# How to Update Existing Code

## For Users with Existing NonBScanner Code

If you have existing code using NonBScanner, here's how to add genome-scale support:

### Minimal Change (Recommended)

```python
# OLD CODE:
from scanner import analyze_sequence
motifs = analyze_sequence(sequence, "my_sequence")

# NEW CODE (automatic optimization):
from auto_scanner import analyze_sequence_auto as analyze_sequence
motifs = analyze_sequence(sequence, "my_sequence")
```

### Explicit Control

```python
# For genomes > 10MB
from genome_scale_scanner import analyze_genome_sequence
motifs = analyze_genome_sequence(large_sequence, "genome")

# For sequences < 10MB (standard)
from scanner import analyze_sequence  
motifs = analyze_sequence(sequence, "region")
```

### No Changes Needed

Your existing code will continue to work as-is. The new functionality is opt-in.

## Command-Line Usage

```bash
# New: Genome-scale analysis
python genome_scale_example.py genome.fasta -o results/genome

# Existing: Web interface still works
streamlit run app.py
```

## Compatibility

- ✅ Existing `scanner.py` unchanged
- ✅ Existing `detectors.py` unchanged  
- ✅ Web interface (`app.py`) unchanged
- ✅ All existing tests pass
- ✅ New features are additive only

## Migration Path

1. **No changes**: Continue using existing code
2. **Try new features**: Use `auto_scanner.py` or `genome_scale_scanner.py`
3. **Gradual adoption**: Use genome-scale for specific large sequences

## What's Different

The genome-scale scanner:
- Disables SlippedDNA and Cruciform (O(n²) complexity)
- Processes in chunks (for memory efficiency)
- Uses optimized regex patterns (faster)

If you need ALL motif classes, continue using standard scanner on smaller sequences.
