# Quick Reference: Hyperscan Integration Pattern

## Core Principles

1. **Scan first, score later** - Hyperscan finds matches, separate functions compute scores
2. **Keep raw per-hit contributions** - Aggregate all hits with per-base redistribution for k-mer detectors
3. **Guarantee non-overlapping** - Deterministic greedy selection by score

---

## Quick Start

### Check Hyperscan Availability
```python
from utils.motif_patterns import HYPERSCAN_AVAILABLE

print(f"Hyperscan available: {HYPERSCAN_AVAILABLE}")
```

### Run Single Detector
```python
from motif_detection.a_philic_detector import APhilicDetector

detector = APhilicDetector()
motifs = detector.detect_motifs(sequence, "my_sequence")
# Returns merged regions only (no individual 10-mers)
```

### Cross-Detector Resolution
```python
from utils.utils import merge_detector_results

# Run multiple detectors
results = {
    'a_philic': detector1.detect_motifs(seq, "seq1"),
    'z_dna': detector2.detect_motifs(seq, "seq1"),
}

# Get non-overlapping motifs (highest score wins)
merged = merge_detector_results(results, overlap_mode='strict')
```

---

## Detector Categories

### K-mer Table Detectors (Merging)
- **A-philic** (`a_philic_detector.py`) - 10-mer log2 odds table
- **Z-DNA** (`z_dna_detector.py`) - 10-mer score table

**Pattern:**
```
Find all 10-mers → Merge overlapping → Redistribute scores → Sum region
```

### Interval Detectors (Overlap Removal)
- **Cruciform** (`cruciform_detector.py`) - Inverted repeats
- **Curved DNA** (`curved_dna_detector.py`) - A-tracts
- **R-loop** (`r_loop_detector.py`) - GC-rich sites
- **G4** (`g_quadruplex_detector.py`) - G-quadruplex

**Pattern:**
```
Find all candidates → Score each → Remove overlaps (greedy by score)
```

---

## Key Functions

### Merging (K-mer Detectors)
```python
def _merge_matches(matches, merge_gap=0):
    """
    Combine overlapping 10-mers into regions.
    merge_gap=0: only overlapping/adjacent
    merge_gap>0: merge if within gap bases
    """
```

### Overlap Removal (Interval Detectors)
```python
def _remove_overlaps(candidates):
    """
    Greedy selection: sort by score (desc), 
    pick highest non-overlapping
    """
```

### Cross-Class Resolution
```python
def resolve_cross_class_overlaps(motifs, mode='strict'):
    """
    mode='strict': highest score wins
    mode='hybrid': keep all, mark overlaps
    """
```

---

## Testing

```bash
# Run all integration tests
python3 test_hyperscan_integration.py

# Expected: ALL TESTS PASSED ✓
```

Tests validate:
- HYPERSCAN_AVAILABLE defined ✓
- Scoring separate from scanning ✓
- K-mer merging works ✓
- Overlap removal works ✓
- Cross-class resolution works ✓

---

## Performance Tips

1. **Compile Hyperscan DB once**, reuse for many sequences
2. **Merging reduces load** - 1000 10-mers → 10 regions
3. **Use adaptive params** - detectors auto-adjust for large sequences

---

## File Locations

| Component | File |
|-----------|------|
| Hyperscan integration | `utils/motif_patterns.py` |
| Cross-class resolution | `utils/utils.py` |
| A-philic detector | `motif_detection/a_philic_detector.py` |
| Z-DNA detector | `motif_detection/z_dna_detector.py` |
| Cruciform detector | `motif_detection/cruciform_detector.py` |
| Tests | `test_hyperscan_integration.py` |
| Full docs | `HYPERSCAN_INTEGRATION.md` |

---

## Common Patterns

### Pattern 1: Single-class analysis with merging
```python
detector = APhilicDetector()
motifs = detector.detect_motifs(seq, "seq1")
# Automatically merged regions
```

### Pattern 2: Multi-class with strict overlap resolution
```python
results = {det_name: det.detect_motifs(seq) for det_name, det in detectors.items()}
merged = merge_detector_results(results, overlap_mode='strict')
# One motif per region (highest score)
```

### Pattern 3: Multi-class with hybrid annotation
```python
merged = merge_detector_results(results, overlap_mode='hybrid')
# All motifs preserved, overlaps can be analyzed
```

---

## Guarantees

✓ **Deterministic** - Same input → same output  
✓ **No duplicates** - Merging/removal prevents split reporting  
✓ **Scientifically accurate** - Literature-based scoring preserved  
✓ **Non-overlapping** - Greedy selection in strict mode  
✓ **Fast** - Hyperscan acceleration where available  

---

For detailed information, see `HYPERSCAN_INTEGRATION.md`
