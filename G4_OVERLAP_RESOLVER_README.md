# G4 Overlap Resolution Agent

## Overview

The G4 Overlap Resolution Agent is an automated tool for detecting and resolving overlapping G-quadruplex (G4) motifs in DNA sequences. It implements a robust overlap resolution algorithm based on G4Hunter scoring and class priority, outputting non-overlapping, biologically meaningful G4 annotations.

## Scientific Background

G-quadruplexes (G4s) are non-canonical DNA secondary structures formed by guanine-rich sequences. They play crucial roles in:
- Transcriptional regulation
- Telomere maintenance
- DNA replication control
- Genome stability

The agent implements patterns and scoring from multiple scientific sources:
- **Canonical G4**: Burge et al. (2006), Todd et al. (2005) [web:67][web:68]
- **Relaxed G4**: Huppert & Balasubramanian (2005), Phan et al. (2006) [web:75]
- **Bulged G4**: Lim et al. (2009), Adrian et al. (2014) [web:22]
- **Multimeric G4**: Guédin et al. (2010), Frasson et al. (2022) [web:42][web:73]
- **G4Hunter Scoring**: Bedrat et al. (2016) [web:3]

## Installation

No additional installation required beyond the main NBDScanner dependencies. The agent uses the existing `motif_detection.g_quadruplex_detector` module.

```bash
# Ensure you're in the NonBScanner directory
cd NonBScanner
```

## Usage

### Command-Line Interface (CLI)

#### Basic Usage with DNA Sequence

```bash
python g4_overlap_resolver.py --sequence "GGGTTAGGGTTAGGGTTAGGG"
```

This outputs JSON format by default:
```json
{
  "version": "G4OverlapResolver-v1.0",
  "analysis_type": "G-Quadruplex_Overlap_Resolution",
  "total_motifs": 1,
  "motifs": [
    {
      "class_name": "canonical_g4",
      "start": 0,
      "end": 21,
      "score": 0.691429,
      "matched_seq": "GGGTTAGGGTTAGGGTTAGGG"
    }
  ]
}
```

#### Output Formats

**JSON Format** (default, with detailed scoring):
```bash
python g4_overlap_resolver.py --sequence "GGGGTTTTGGGGTTTTGGGG" --format json
```

**Tab-Delimited Table** (for Excel, spreadsheets):
```bash
python g4_overlap_resolver.py --sequence "GGGGTTTTGGGGTTTTGGGG" --format table
```

**BED Format** (for genome browsers like UCSC, IGV):
```bash
python g4_overlap_resolver.py --sequence "GGGGTTTTGGGGTTTTGGGG" --format bed
```

#### FASTA Input

For FASTA files:
```bash
python g4_overlap_resolver.py --fasta my_sequence.fasta --format table --output results.tsv
```

Example FASTA file:
```
>telomeric_repeat
GGGTTAGGGTTAGGGTTAGGG
>promoter_region
GGGGTTTTGGGGTTTTGGGGTTTTGGGG
```

#### Additional Options

```bash
# Specify sequence name
python g4_overlap_resolver.py --sequence "GGGTTAGGG" --name "my_sequence"

# Save output to file
python g4_overlap_resolver.py --sequence "GGGTTAGGG" --output results.json

# Compact JSON (no indentation)
python g4_overlap_resolver.py --sequence "GGGTTAGGG" --format json --compact

# Get help
python g4_overlap_resolver.py --help
```

### Python API

#### Basic API Usage

```python
from g4_overlap_resolver import G4OverlapResolver

# Initialize resolver
resolver = G4OverlapResolver()

# Analyze sequence
sequence = "GGGTTAGGGTTAGGGTTAGGG"
annotations = resolver.resolve_and_annotate(sequence, "my_sequence")

# annotations is a list of dictionaries with:
# - class_name: G4 motif class
# - start, end: 0-based positions
# - score: G4Hunter-based score
# - matched_seq: actual DNA sequence
# - details: scoring breakdown

for ann in annotations:
    print(f"{ann['class_name']}: {ann['start']}-{ann['end']}, score={ann['score']:.4f}")
```

#### Format Output

```python
# JSON format
json_output = resolver.format_as_json(annotations, pretty=True)
print(json_output)

# Table format
table_output = resolver.format_as_table(annotations)
print(table_output)

# BED format
bed_output = resolver.format_as_bed(annotations, "chr1")
print(bed_output)
```

#### Complete Workflow Example

```python
from g4_overlap_resolver import G4OverlapResolver
import json

# Initialize
resolver = G4OverlapResolver()

# Multiple sequences
sequences = {
    "telomeric": "GGGTTAGGGTTAGGGTTAGGG",
    "promoter": "GGGGTTTTGGGGTTTTGGGGTTTTGGGG",
    "complex": "GGGTTGGGTTGGGTTGGGAAAAAGGGGTTTTGGGG"
}

# Analyze all sequences
all_results = []
for name, seq in sequences.items():
    annotations = resolver.resolve_and_annotate(seq, name)
    all_results.extend(annotations)

# Save as JSON
with open('g4_analysis.json', 'w') as f:
    json.dump({
        'total_motifs': len(all_results),
        'motifs': all_results
    }, f, indent=2)

# Save as table
with open('g4_analysis.tsv', 'w') as f:
    f.write(resolver.format_as_table(all_results))
```

## Output Format Details

### JSON Output Structure

```json
{
  "version": "G4OverlapResolver-v1.0",
  "analysis_type": "G-Quadruplex_Overlap_Resolution",
  "description": "Non-overlapping G4 motifs resolved by score and class priority",
  "total_motifs": 1,
  "motifs": [
    {
      "sequence_name": "input_seq",
      "class_name": "canonical_g4",
      "pattern_id": "G4_6_1",
      "start": 0,
      "end": 21,
      "length": 21,
      "score": 0.691429,
      "matched_seq": "GGGTTAGGGTTAGGGTTAGGG",
      "details": {
        "n_g_tracts": 4,
        "total_g_len": 12,
        "gc_balance": 0.5714,
        "max_window_abs": 12.0,
        "normalized_window": 0.571429,
        "tract_bonus": 0.12,
        "gc_penalty": 0.0,
        "normalized_score": 0.691429,
        "region_score": 0.691429
      }
    }
  ]
}
```

### Table Output Format

Tab-delimited table with columns:
- `Sequence_Name`: Identifier for the input sequence
- `Class`: G4 motif class (canonical_g4, relaxed_g4, etc.)
- `Pattern_ID`: Specific pattern identifier
- `Start`: 0-based start position
- `End`: 0-based end position (exclusive)
- `Length`: Motif length in base pairs
- `Score`: G4Hunter-based score
- `Sequence`: Actual DNA sequence

Example:
```
Sequence_Name	Class	Pattern_ID	Start	End	Length	Score	Sequence
input_seq	canonical_g4	G4_6_1	0	28	28	0.7616	GGGGTTTTGGGGTTTTGGGGTTTTGGGG
```

### BED Format

Standard 6-column BED format for genome browsers:
```
chr1	0	28	canonical_g4_G4_6_1	76	+
```

Columns: chromosome, start, end, name, score, strand

## G4 Motif Classes

The agent detects and prioritizes multiple G4 classes:

1. **canonical_g4** (highest priority)
   - Standard G4 with 4 G-tracts
   - Loop lengths 1-7bp
   - Highest stability

2. **relaxed_g4**
   - G4 with longer loops (up to 12bp)
   - Still biologically relevant

3. **multimeric_g4**
   - Multiple G4 units
   - Higher-order structures

4. **bulged_g4**
   - G4 with bulge loops (8-40bp)
   - Alternative topologies

5. **imperfect_g4**
   - G4-like with interruptions
   - G/A substitutions or mismatches

6. **bipartite_g4** (lowest priority)
   - Two-block G4 structures
   - Very long loops

7. **g_triplex**
   - Three-tract G structures
   - G4 intermediates

## Overlap Resolution Algorithm

The agent uses a greedy algorithm for overlap resolution:

1. **Scoring**: All candidate G4 motifs are scored using G4Hunter algorithm
   - Window-based scoring (default 25bp)
   - G-tract bonuses for multiple runs
   - GC balance penalties

2. **Sorting**: Candidates sorted by:
   - Score (descending)
   - Class priority (canonical > relaxed > bulged > imperfect)
   - Length (descending)

3. **Greedy Selection**: 
   - Iterate through sorted candidates
   - Accept non-overlapping motifs
   - Skip conflicting regions

Time complexity: O(n log n) for sorting + O(n) for selection

## Testing

Run the test suite:

```bash
python test_g4_overlap_resolver.py
```

The test suite includes:
- Basic G4 detection tests
- Overlap resolution tests
- Multiple output format tests
- Edge case handling
- API usage tests

All 18 tests should pass.

## Integration with NBDScanner

The G4 Overlap Resolution Agent is designed to complement the main NBDScanner application:

```python
# Use within NBDScanner workflow
from g4_overlap_resolver import G4OverlapResolver
from utils.utils import parse_fasta

# Load sequences
with open('input.fasta', 'r') as f:
    sequences = parse_fasta(f)

# Analyze with G4 resolver
resolver = G4OverlapResolver()
for name, seq in sequences.items():
    annotations = resolver.resolve_and_annotate(seq, name)
    # Process annotations...
```

## Performance Considerations

- **Time Complexity**: O(n·m) for pattern matching + O(k log k) for sorting
  - n = sequence length
  - m = number of patterns
  - k = number of candidates
  
- **Memory**: O(k) for storing candidates

For typical sequences (<100kb), analysis completes in milliseconds.

## Troubleshooting

### No G4s Detected

If no G4 motifs are detected:
1. Check that sequence contains G-rich regions
2. Verify sequence format (A, T, G, C)
3. Try longer sequences (minimum ~15-20bp for G4s)

### Overlapping Annotations

If you see overlapping annotations, please report as a bug. The overlap resolution should guarantee non-overlapping output.

### FASTA Parsing Issues

Ensure FASTA files:
- Start with `>` header line
- Contain only DNA sequences (A, T, G, C)
- Have no extra characters or whitespace in sequences

## References

1. Burge, S., et al. (2006). "Quadruplex DNA: sequence, topology and structure." *Nucleic Acids Research* [web:67][web:68]
2. Huppert, J.L. & Balasubramanian, S. (2005). "Prevalence of quadruplexes in the human genome." *Nucleic Acids Research* [web:75]
3. Lim, K.W., et al. (2009). "Structure of the human telomere in K+ solution: a stable basket-type G-quadruplex with only two G-tetrad layers." *JACS* [web:22]
4. Bedrat, A., et al. (2016). "Re-evaluation of G-quadruplex propensity with G4Hunter." *Nucleic Acids Research* [web:3]
5. Guédin, A., et al. (2010). "How long is too long? Effects of loop size on G-quadruplex stability." *NAR* [web:42]
6. Frasson, I., et al. (2022). "Multimeric G-quadruplexes: a review." *NAR* [web:73]

## License

Part of NBDScanner project. See main repository LICENSE file.

## Contact

For issues or questions, please open an issue in the NBDScanner repository or contact Dr. Venkata Rajesh Yella.
