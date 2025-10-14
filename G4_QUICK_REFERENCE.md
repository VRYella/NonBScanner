# G4 Overlap Resolver - Quick Reference Card

## One-Line Description
Automated agent for detecting and resolving overlapping G-quadruplex motifs in DNA sequences.

## Quick Start

### CLI
```bash
# Basic usage
python g4_overlap_resolver.py --sequence "GGGTTAGGGTTAGGGTTAGGG"

# From FASTA file
python g4_overlap_resolver.py --fasta input.fasta --format table --output results.tsv

# Get help
python g4_overlap_resolver.py --help
```

### Python API
```python
from g4_overlap_resolver import G4OverlapResolver

resolver = G4OverlapResolver()
annotations = resolver.resolve_and_annotate(sequence, "my_seq")
print(resolver.format_as_json(annotations))
```

## Output Formats

| Format | Use Case | Example |
|--------|----------|---------|
| JSON | APIs, detailed analysis | `--format json` |
| Table | Excel, spreadsheets | `--format table` |
| BED | Genome browsers (UCSC, IGV) | `--format bed` |

## G4 Motif Classes (Priority Order)

1. **canonical_g4** - Standard 4-tract G4 (highest priority)
2. **relaxed_g4** - Longer loops (up to 12bp)
3. **multimeric_g4** - Multiple G4 units
4. **bulged_g4** - Large loop G4
5. **imperfect_g4** - G4-like with interruptions
6. **bipartite_g4** - Two-block structures
7. **g_triplex** - Three-tract structures

## Key Features

- **Overlap Resolution**: Greedy algorithm, O(n log n)
- **Scoring**: G4Hunter-based (Bedrat et al. 2016)
- **Priority**: Score > Class > Length
- **Output**: Non-overlapping, annotated motifs

## Common Commands

```bash
# JSON with pretty printing
python g4_overlap_resolver.py --sequence "GGGG..." --format json

# Compact JSON for APIs
python g4_overlap_resolver.py --sequence "GGGG..." --format json --compact

# Save to file
python g4_overlap_resolver.py --sequence "GGGG..." --output results.json

# BED format for genome browser
python g4_overlap_resolver.py --fasta seq.fasta --format bed --output track.bed

# Table for Excel
python g4_overlap_resolver.py --fasta seq.fasta --format table --output results.tsv
```

## Output Fields

### JSON/Table
- `sequence_name` - Input sequence identifier
- `class_name` - G4 motif class
- `pattern_id` - Pattern identifier
- `start` - 0-based start position
- `end` - 0-based end position (exclusive)
- `length` - Motif length in bp
- `score` - G4Hunter score
- `matched_seq` - DNA sequence
- `details` - Scoring breakdown

### BED
- chromosome, start, end, name, score, strand

## Scientific References

- Burge 2006, Todd 2005: Canonical G4 [web:67][web:68]
- Huppert 2005, Phan 2006: Relaxed G4 [web:75]
- Lim 2009, Adrian 2014: Bulged G4 [web:22]
- Bedrat et al. 2016: G4Hunter scoring [web:3]
- Guédin 2010, Frasson 2022: Multimeric G4 [web:42][web:73]

## Testing

```bash
# Run test suite
python test_g4_overlap_resolver.py

# Run examples
python example_g4_overlap_resolver.py
```

## Full Documentation

See: `G4_OVERLAP_RESOLVER_README.md`

## Support

For issues or questions, see the main NBDScanner repository.
