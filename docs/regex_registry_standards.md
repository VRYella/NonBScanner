# Regex Registry Scientific Standards Documentation

## Overview

The NonBScanner regex registry has been comprehensively updated to align with standard scientific definitions from major tools and literature. All patterns are now Hyperscan-compatible and include proper scientific annotations.

## Scientific Standards Applied

### G-Quadruplex Patterns (Class 6)
**Standard**: G4Hunter algorithm (Bedrat et al. 2016, NAR 44:1746-1759)

- **Canonical G4**: G3+N1-7G3+N1-7G3+N1-7G3+ (Williamson 2005)
- **Relaxed G4**: G2+N1-12G2+N1-12G2+N1-12G2+ (Kikin et al. 2006)
- **Bulged G4**: G3+N0-3G3+N0-3G3+N0-3G3+ (Todd et al. 2005)
- **Two-quartet G4**: G2N1-7G2N1-7G2N1-7G2 (Adrian et al. 2014)
- **Scoring**: G4Hunter score thresholds (1.2+ for canonical, 0.8+ for relaxed)

### I-Motif Patterns (Class 7) 
**Standard**: C-rich quadruplex structures (Zeraati et al. 2018, Nat Chem 10:631-637)

- **Canonical iMotif**: C3+N1-7C3+N1-7C3+N1-7C3+ (pH 5.0 optimal)
- **Relaxed iMotif**: C2+N1-12C2+N1-12C2+N1-12C2+ (pH 5.5)
- **AC-motif**: A3N4-6C3N4-6C3N4-6C3 (Felsenfeld 1967)
- **Intercalated**: C3+N0-3C3+N0-3C3+N0-3C3+ (very stable at pH 4.5)

### Z-DNA Patterns (Class 8)
**Standard**: Z-DNA seeker algorithm (Ho et al. 1986, Rich & Zhang 2003)

- **CG Z-DNA**: (CG)6+ alternating (Z-score >8.0)
- **AT Z-DNA**: (AT)8+ alternating (Z-score >3.0)
- **Mixed Z-DNA**: CG/AT mixed patterns (Z-score >5.0)
- **eGZ DNA**: CGG expansion slippage (Napierala et al. 1997)

### Curved DNA Patterns (Class 1)
**Standard**: A-tract curvature models (Crothers et al. 1992, Marini et al. 1982)

- **A-tracts**: A4+/T4+ runs with curvature scoring
- **Phased A-tracts**: ~10 bp periodic arrays (10.0-10.6 bp period)
- **Global Arrays**: Multiple phased elements (Yella & Bansal 2017)

### Triplex Patterns (Class 5)
**Standard**: Triple-stranded DNA (Frank-Kamenetskii & Mirkin 1995)

- **Homopurine**: [AG]15+ tracts (high stability)
- **Homopyrimidine**: [CT]15+ tracts (sticky DNA)
- **Mirror Repeats**: Palindromic Pu/Py sequences
- **H-DNA**: Intramolecular triplex structures

### Cruciform Patterns (Class 3)
**Standard**: Inverted repeats (Lilley & Clegg 1993)

- **Perfect Palindromes**: Exact inverted repeats (6-20 bp arms)
- **Imperfect Palindromes**: Minor mismatches allowed
- **Stem-loops**: Variable loop structures
- **AT/GC-rich**: Composition-specific cruciforms

### R-Loop Patterns (Class 4)
**Standard**: R-loop forming sequences (Ginno et al. 2012, Skourti-Stathaki et al. 2011)

- **RLFS Model 1**: G3+N1-10G3+(G3+N1-10)+ (>50% GC)
- **RLFS Model 2**: G4+(G4+N1-10)+ (>60% GC)
- **QmRLFS-RE**: Combined G4+G-rich regions
- **CpG RLFS**: CpG island-associated sequences

### Slipped DNA Patterns (Class 2)
**Standard**: STR definitions (Wells et al. 2005, Pearson et al. 2005)

- **Trinucleotide STRs**: Disease-associated repeats (CAG, CGG, GAA)
- **Dinucleotide STRs**: Common microsatellites
- **Mononucleotide**: MSI-associated poly-runs
- **Disease thresholds**: Based on clinical expansion data

## Pattern Annotation Format

Each pattern uses a 9-tuple structure:
```python
(regex_pattern, pattern_id, group_number, subclass_name, 
 scoring_function, score_scale, min_runs, min_score, score_method)
```

### Scoring Methods Reference

- **G4Hunter**: Bedrat et al. 2016 algorithm
- **iM_Hunter**: G4Hunter adapted for C-richness
- **Z_Seeker**: Ho et al. 1986 dinucleotide scoring
- **Curvature_Calladine**: Calladine & Drew curvature model
- **Triplex_Stability**: Frank-Kamenetskii stability calculations
- **STR_Score**: Expansion threshold-based scoring
- **Palindrome_Score**: Thermodynamic stability prediction

## Hyperscan Compatibility

All patterns are designed for Intel Hyperscan PCRE compatibility:

✅ **Compatible Features**:
- Character classes `[ACGT]` 
- Quantifiers `{n,m}`
- Non-capturing groups `(?:...)`
- Alternation `|`

❌ **Avoided Features**:
- Back-references `\1, \2`
- Word boundaries `\b, \B`
- Complex lookahead/lookbehind
- Recursive patterns

## Testing and Validation

The registry includes comprehensive validation:

- **Syntax Validation**: All patterns compile successfully
- **Scientific Validation**: Known sequences match expected patterns
- **Performance Testing**: Sub-millisecond compilation times
- **Hyperscan Testing**: Compatible with Intel Hyperscan library

## Usage in Motif Detection

Patterns are accessed through the centralized registry:

```python
from core.regex_registry import get_patterns_for_motif

# Get all G4 patterns
g4_patterns = get_patterns_for_motif('g_quadruplex')

# Get specific pattern info
pattern_info = get_pattern_info(pattern_id)
```

## Literature References

1. Bedrat, A. et al. (2016) G4Hunter. Nucleic Acids Research 44:1746-1759
2. Zeraati, M. et al. (2018) I-motif DNA structures. Nature Chemistry 10:631-637
3. Ho, P.S. et al. (1986) Z-DNA structure. EMBO Journal 5:2737-2744
4. Frank-Kamenetskii, M.D. & Mirkin, S.M. (1995) Triplex DNA. Annu Rev Biochem 64:65-95
5. Lilley, D.M.J. & Clegg, R.M. (1993) Cruciform structures. Annu Rev Biophys 22:299-328
6. Crothers, D.M. et al. (1992) Curved DNA. Methods Enzymol 212:3-29
7. Ginno, P.A. et al. (2012) R-loops. Molecular Cell 45:511-522
8. Wells, R.D. et al. (2005) Slipped-strand structures. Trends Biochem Sci 30:437-444