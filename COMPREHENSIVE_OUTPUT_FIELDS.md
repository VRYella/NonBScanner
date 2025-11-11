# Comprehensive Output Fields Documentation

## Overview

NonBScanner now provides comprehensive output fields for all detected motifs. All results include 26 standard fields, with additional detector-specific fields when applicable.

## Standard Output Fields

All motifs include these 26 fields (non-applicable fields show "NA"):

| # | Field Name | Description | Example |
|---|------------|-------------|---------|
| 1 | **ID** | Unique identifier for the motif | `test_seq_SLP_STR_6_1` |
| 2 | **Sequence_Name** | Sequence name or accession | `test_sequence` |
| 3 | **Source** | Source reference (genome, experiment, study) | `Wells 2005` |
| 4 | **Class** | Motif class | `Slipped_DNA` |
| 5 | **Subclass** | Motif subclass | `STR` |
| 6 | **Pattern_ID** | Pattern/annotation identifier | `SLP_STR_6` |
| 7 | **Start** | Start position (1-based) | `1` |
| 8 | **End** | End position | `18` |
| 9 | **Length** | Length in base pairs | `18` |
| 10 | **Sequence** | Motif sequence | `GGGTTAGGGTTAGGGTTA` |
| 11 | **Method** | Detection method used | `Slipped_DNA_detection` |
| 12 | **Score** | Motif score (0-1 or raw score) | `0.735` |
| 13 | **Repeat_Type** | Type of repeat/tract | `6-mer STR` |
| 14 | **Left_Arm** | Left arm sequence (for palindromes) | `ATATATAT` |
| 15 | **Right_Arm** | Right arm sequence (for palindromes) | `ATATATAT` |
| 16 | **Loop_Seq** | Loop sequence (for hairpins/cruciforms) | `TTTT` |
| 17 | **Arm_Length** | Length of arms in base pairs | `8` |
| 18 | **Loop_Length** | Length of loop in base pairs | `4` |
| 19 | **Stem_Length** | Stem length(s) for structured motifs | `NA` |
| 20 | **Unit_Length** | Length of repeat unit | `6` |
| 21 | **Number_Of_Copies** | Number of copies/repeats | `3` |
| 22 | **Spacer_Length** | Spacer length in base pairs | `4` |
| 23 | **Spacer_Sequence** | Spacer sequence | `TTTT` |
| 24 | **GC_Content** | GC content percentage | `50.0` |
| 25 | **Structural_Features** | Additional structural features | `Tract:A-tract` |
| 26 | **Strand** | Strand orientation | `+` |

## Field Population by Motif Class

### Slipped DNA (STRs and Direct Repeats)

**Short Tandem Repeats (STR):**
- `Repeat_Type`: Unit size (e.g., "6-mer STR")
- `Unit_Length`: Repeat unit size
- `Number_Of_Copies`: Number of repeat units
- `GC_Content`: GC percentage

**Direct Repeats:**
- `Unit_Length`: Length of repeated unit
- `Spacer_Length`: Length of spacer between repeats
- `Spacer_Sequence`: Actual spacer sequence
- `GC_Content`: GC percentage

### Cruciform (Inverted Repeats)

- `Left_Arm`: Left palindromic arm sequence
- `Right_Arm`: Right palindromic arm sequence
- `Loop_Seq`: Loop/spacer sequence
- `Arm_Length`: Length of palindromic arms
- `Loop_Length`: Length of loop
- `GC_Content`: GC percentage

### G-Quadruplex, Z-DNA, R-Loop, etc.

Standard fields only:
- `ID`, `Sequence_Name`, `Class`, `Subclass`, `Start`, `End`, `Length`
- `Sequence`, `Method`, `Score`, `GC_Content`, `Strand`
- Other fields show "NA"

## Export Formats

### DataFrame Export

```python
from scanner import analyze_sequence, export_results_to_dataframe

motifs = analyze_sequence(sequence, "my_sequence")
df = export_results_to_dataframe(motifs)
# df contains all 26 standard columns with "NA" for non-applicable fields
```

### CSV Export

```python
from utilities import export_to_csv

csv_data = export_to_csv(motifs)
# CSV includes all 26 standard columns plus any detector-specific fields
# Non-applicable fields contain "NA"
```

### JSON Export

```python
from utilities import export_to_json

json_data = export_to_json(motifs, pretty=True)
# JSON includes all fields from motif dictionaries
```

## Examples

### STR Output
```csv
ID,Sequence_Name,Source,Class,Subclass,Unit_Length,Number_Of_Copies,GC_Content,...
test_SLP_STR_6_1,test,Wells 2005,Slipped_DNA,STR,6,3,50.0,...
```

### Cruciform Output
```csv
ID,Sequence_Name,Class,Subclass,Left_Arm,Right_Arm,Loop_Seq,Arm_Length,Loop_Length,...
test_CRU_1,test,Cruciform,Inverted_Repeat,ATATATAT,ATATATAT,TTTT,8,4,...
```

### Direct Repeat Output
```csv
ID,Sequence_Name,Class,Subclass,Unit_Length,Spacer_Length,Spacer_Sequence,GC_Content,...
test_SLP_DIR_1,test,Slipped_DNA,Direct_Repeat,10,4,TTTT,0.0,...
```

## Notes

- All fields are guaranteed to be present in DataFrame and CSV exports
- Non-applicable fields contain the string "NA" (not NaN or None)
- Numeric fields may be integers or floats depending on the measurement
- Empty sequences (e.g., when spacer_length=0) are shown as "NA"
- Additional detector-specific fields may be included beyond the 26 standard fields

## Testing

Run the comprehensive output test:
```bash
python3 test_comprehensive_output.py
```

Run the example demonstration:
```bash
python3 example_comprehensive_output.py
```
