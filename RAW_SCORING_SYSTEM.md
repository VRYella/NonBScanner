# Raw Scoring System Documentation

## Overview

NBDFinder uses a **raw scoring system** that provides scientifically accurate scores based on published algorithms for each motif class. Unlike normalized scoring systems that artificially scale all scores to 0-1, raw scoring preserves the natural scale and meaning of each algorithm.

## Advantages of Raw Scoring

1. **Scientific Accuracy**: Scores retain their original meaning from published research
2. **No Information Loss**: Raw values allow proper comparison within the same motif class
3. **Interpretability**: Researchers can directly relate scores to literature values
4. **Performance**: No additional normalization calculations needed

## Score Ranges by Motif Class

Different motif classes naturally produce different score ranges because they use different algorithms:

| Motif Class | Algorithm | Typical Range | Reference |
|-------------|-----------|---------------|-----------|
| G-Quadruplex | G4Hunter | 0 - 3+ | Bedrat et al., 2016 |
| Z-DNA | Alternating purine-pyrimidine | 0 - 1000+ | Rich & Zhang, 2003 |
| A-philic DNA | A-tract density | 0 - 100+ | Satchwell et al., 1986 |
| i-Motif | C-richness score | 0 - 3+ | Similar to G4Hunter |
| Curved DNA | Curvature propensity | 0 - 1 | Gabrielian & Pongor, 1996 |
| Triplex | Purine/pyrimidine content | 0 - 1 | Frank-Kamenetskii, 1995 |
| R-Loop | GC-skew and density | 0 - 1 | Ginno et al., 2012 |
| Cruciform | Inverted repeat score | 0 - 1 | Lilley, 1980 |

## Quality Thresholds

Each motif class has appropriate minimum thresholds based on raw scores:

```python
min_scores = {
    'g_quadruplex': 0.5,    # G4Hunter score
    'curved_dna': 0.2,      # Low threshold for curvature
    'z_dna': 0.5,           # Raw alternating score
    'triplex': 0.6          # Purine/pyrimidine threshold
}
```

## Score Interpretation Examples

### G-Quadruplex (G4Hunter)
- **Score 0.5**: Weak G4-forming potential
- **Score 1.0**: Moderate G4-forming potential
- **Score 2.0+**: Strong G4-forming potential
- **Score 3.0+**: Very strong G4-forming potential

### Z-DNA (Alternating Pattern)
- **Score 0-100**: Low Z-DNA forming potential
- **Score 100-500**: Moderate Z-DNA forming potential
- **Score 500+**: High Z-DNA forming potential
- **Score 1000+**: Very high Z-DNA forming potential (long alternating runs)

### A-philic DNA
- **Score 0-10**: Low A-richness
- **Score 10-30**: Moderate A-richness
- **Score 30+**: High A-richness (strong A-philic binding)

## Comparison with Normalized Scoring

### Old Normalized System (Removed)
- All scores artificially scaled to 0-1
- Loss of scientific meaning
- Cannot compare across different datasets
- Required class-specific normalization parameters

### New Raw Scoring System
- Preserves original algorithm output
- Direct correspondence to literature
- Enables cross-dataset comparison
- Simpler and faster

## Usage

When analyzing motifs, use the `Score` field directly:

```python
from utils.modular_scanner import ModularMotifDetector

scanner = ModularMotifDetector()
motifs = scanner.analyze_sequence(sequence, "my_sequence")

for motif in motifs:
    motif_class = motif['Class']
    score = motif['Score']
    
    if motif_class == 'G-Quadruplex' and score > 2.0:
        print(f"Strong G4 at {motif['Start']}-{motif['End']}: score={score:.3f}")
    elif motif_class == 'Z-DNA' and score > 500:
        print(f"High Z-DNA potential at {motif['Start']}-{motif['End']}: score={score:.0f}")
```

## Backward Compatibility

The system maintains backward compatibility:
- Old code using `Normalized_Score` will still work (falls back to `Score`)
- Visualizations automatically use raw scores when available
- `canonicalize_motif()` utility handles both formats

## References

1. Bedrat, A., et al. (2016). "Re-evaluation of G-quadruplex propensity with G4Hunter." *Nucleic Acids Research*, 44(4), 1746-1759.
2. Rich, A., & Zhang, S. (2003). "Z-DNA: the long road to biological function." *Nature Reviews Genetics*, 4(7), 566-572.
3. Satchwell, S. C., et al. (1986). "Sequence periodicities in chicken nucleosome core DNA." *Journal of Molecular Biology*, 191(4), 659-675.

## Testing

Run the comprehensive test suite:

```bash
python test_raw_scoring.py
```

This validates:
- Raw scores are present
- No normalized scores in new detections
- Class-specific score ranges are correct
- No overlaps within same subclass
- Performance is maintained
