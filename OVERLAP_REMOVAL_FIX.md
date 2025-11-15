# Overlap Removal Fix - Summary

## Problem
The NonBScanner was reporting multiple overlapping motifs of the same subclass, leading to redundant and confusing results.

## Example Issue
For a G-Quadruplex sequence, 9 different overlapping G4 subclass patterns were being reported:
- Canonical G4 (2 matches)
- Relaxed G4 (1 match)
- Long-loop G4 (1 match)
- Bulged G4 (1 match)
- Multimeric G4 (1 match)
- Imperfect G4 (1 match)
- G-Triplex (2 matches)

All overlapping at positions 1-47 in the same 47bp sequence.

## Solution
Modified GQuadruplexDetector and IMotifDetector to use their existing overlap resolution methods:

1. **GQuadruplexDetector.detect_motifs()**: Now calls `annotate_sequence()` which uses `_resolve_overlaps()`
2. **IMotifDetector.detect_motifs()**: Now calls `annotate_sequence()` which uses `_resolve_overlaps_greedy()`

Both methods:
- Sort motifs by score (descending), then by length (descending)
- Select the highest-scoring motif first
- Skip any motifs that overlap with already-selected motifs
- Result: Only non-overlapping, highest-quality motifs are reported

## Results

### Before Fix
```
Test Sequence: GGGTTAGGGTTAGGGTTAGGGAAAAAGGGTTAGGGTTAGGGTTAGGG (47 bp)

Found 9 overlapping G4 motifs:
  1. Canonical G4     pos 1-21  (len=21) score=0.691
  2. Canonical G4     pos 27-47 (len=21) score=0.691
  3. Relaxed G4       pos 1-41  (len=41) score=1.279
  4. Long-loop G4     pos 1-35  (len=35) score=1.008
  5. Bulged G4        pos 1-29  (len=29) score=0.766
  6. Multimeric G4    pos 1-47  (len=47) score=1.579  ← HIGHEST SCORE
  7. Imperfect G4     pos 1-29  (len=29) score=0.766
  8. G-Triplex        pos 1-15  (len=15) score=0.660
  9. G-Triplex        pos 19-35 (len=17) score=0.589
```

### After Fix
```
Test Sequence: GGGTTAGGGTTAGGGTTAGGGAAAAAGGGTTAGGGTTAGGGTTAGGG (47 bp)

Found 1 G4 motif:
  1. Multimeric G4    pos 1-47  (len=47) score=1.579
```

## Key Benefits

1. **Cleaner Results**: No redundant overlapping motifs
2. **Accurate**: Only the best (highest-scoring/longest) motif per region
3. **Consistent**: All detector classes handle overlaps uniformly
4. **Validated**: Comprehensive test coverage (8 unit tests, 5 integration tests)

## Test Coverage

### Unit Tests (test_overlap_removal.py)
- ✅ G-Quadruplex overlap removal
- ✅ G-Quadruplex non-overlapping motifs
- ✅ i-Motif overlap removal
- ✅ Curved DNA overlap removal
- ✅ Cruciform overlap removal
- ✅ Scanner-level overlap removal
- ✅ Highest score selection
- ✅ Longest motif selection

### Integration Tests
- ✅ G-Quadruplex (overlapping subclasses): 9 candidates → 1 final motif
- ✅ i-Motif (overlapping patterns): Properly resolved
- ✅ Curved DNA (multiple A-tracts): Non-overlapping within subclasses
- ✅ Slipped DNA (STR): Working correctly
- ✅ Mixed classes: Different classes can overlap as expected

### Regression Tests
- ✅ Example file (28 sequences): All analyzed correctly
- ✅ No regressions in motif detection
- ✅ CodeQL security scan: 0 vulnerabilities

## Technical Details

### Modified Files
- `detectors.py`:
  - Added `GQuadruplexDetector.detect_motifs()` (66 lines)
  - Updated `IMotifDetector.detect_motifs()` (77 lines)

### Overlap Resolution Logic
Both detectors use similar logic:
1. Find all candidate motifs
2. Score each candidate
3. Sort by score (desc), then length (desc)
4. Iterate through sorted list:
   - If motif doesn't overlap with any selected motifs, add it
   - Otherwise, skip it
5. Return non-overlapping set

### Selection Criteria
1. **Primary**: Highest score
2. **Secondary** (if scores equal): Longest length
3. **Tertiary** (for G4): Class priority (canonical > relaxed > long-loop > multimeric > bulged > imperfect > g_triplex)

## Verification

Run the test suite:
```bash
python3 test_overlap_removal.py
```

Expected output:
```
======================================================================
OVERLAP REMOVAL TEST SUITE
======================================================================
...
======================================================================
RESULTS: 8 passed, 0 failed
======================================================================
```

## Backward Compatibility

✅ Fully backward compatible:
- No changes to API or data structures
- Only affects internal motif selection logic
- Reduces output (fewer redundant motifs) but doesn't break existing code
- All other detectors (CurvedDNA, Cruciform, R-Loop, Z-DNA, etc.) continue to work as before
