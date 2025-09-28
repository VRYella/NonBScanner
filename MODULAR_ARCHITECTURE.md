# Modular Motif Detection Architecture

## Overview

The NonBScanner codebase has been successfully refactored from a centralized pattern registry approach to a **modular architecture** where each motif class has its own dedicated detector. This change addresses the requirements to:

1. ✅ Create a `motif_detection` folder  
2. ✅ Split the code so each motif has its own separate code file
3. ✅ Remove the central registry file and split patterns for each motif

## New Architecture

### Directory Structure
```
motif_detection/
├── __init__.py                 # Module interface
├── base_detector.py           # Abstract base class
├── curved_dna_detector.py     # A-tract mediated bending
├── slipped_dna_detector.py    # STR and direct repeats  
├── cruciform_detector.py      # Inverted repeats
├── r_loop_detector.py         # R-loop formation sites
├── triplex_detector.py        # Three-strand structures
├── g_quadruplex_detector.py   # G4 and G-triplex
├── i_motif_detector.py        # C-rich structures
├── z_dna_detector.py          # Left-handed helix
└── a_philic_detector.py       # A-rich tracts
```

### Key Components

#### BaseMotifDetector
- Abstract base class defining common interface
- Handles pattern compilation, motif detection, and scoring
- Provides quality thresholds and statistics methods

#### Individual Detectors
Each detector contains:
- **Patterns**: Specific to that motif class (no shared registry)
- **Scoring algorithms**: Motif-specific scoring methods (G4Hunter, curvature, instability, etc.)
- **Quality thresholds**: Class-specific confidence thresholds
- **Biological metadata**: References and significance descriptions

## Pattern Distribution

| Motif Class | Detector File | Patterns | Key Features |
|-------------|---------------|----------|--------------|
| Curved DNA | `curved_dna_detector.py` | 5 | A-tract curvature, phasing analysis |
| Slipped DNA | `slipped_dna_detector.py` | 8 | STR instability, direct repeats |
| Cruciform | `cruciform_detector.py` | 3 | Palindrome stability analysis |
| R-loop | `r_loop_detector.py` | 5 | GC skew, QmRLFS models |
| Triplex | `triplex_detector.py` | 5 | Triplex potential, sticky DNA |
| G-Quadruplex | `g_quadruplex_detector.py` | 13 | G4Hunter algorithm, 7 subclasses |
| i-Motif | `i_motif_detector.py` | 6 | C-tract analysis, pH stability |
| Z-DNA | `z_dna_detector.py` | 5 | Purine-pyrimidine alternation |
| A-philic | `a_philic_detector.py` | 4 | A-tract propensity |

**Total: 54 patterns across 9 specialized detectors**

## Usage

### New Modular Approach (Default)
```python
from nbdscanner import analyze_sequence

# Uses modular architecture by default
motifs = analyze_sequence("GGGTTAGGGTTAGGGTTAGGG", use_modular=True)
```

### Individual Detectors
```python
from motif_detection import GQuadruplexDetector

detector = GQuadruplexDetector()
motifs = detector.detect_motifs("GGGTTAGGGTTAGGGTTAGGG", "seq_name")
```

### Legacy Support
```python
from nbdscanner import analyze_sequence

# Use legacy centralized approach if needed
motifs = analyze_sequence("GGGTTAGGGTTAGGGTTAGGG", use_modular=False)
```

## Benefits

### 1. **Maintainability**
- Each motif class is self-contained
- Easy to modify or extend individual detectors
- Clear separation of concerns

### 2. **Extensibility**
- Add new motif classes without touching existing code
- Easy to implement class-specific improvements
- Modular testing and validation

### 3. **Performance**
- Improved G4Hunter algorithm
- Optimized scoring methods
- Class-specific thresholds

### 4. **Code Organization**
- No central registry file to maintain
- Patterns stored with their specific detectors
- Clear biological context and references

## Backward Compatibility

The refactor maintains **100% backward compatibility**:

- All existing APIs work unchanged
- Original `nbdscanner.py` still functional
- Applications (Streamlit, FastAPI) work without changes
- Pure Python scanner integration preserved

## Testing

Comprehensive test suite validates:
- ✅ Individual detector functionality
- ✅ Modular architecture integration  
- ✅ API compatibility between approaches
- ✅ Pattern coverage and scoring
- ✅ Application compatibility

Run tests:
```bash
python test_modular_architecture.py
python test_split_approach.py
```

## Migration Guide

### For Users
No changes required - the modular architecture is now the default.

### For Developers
To add a new motif class:

1. Create `new_motif_detector.py` in `motif_detection/`
2. Inherit from `BaseMotifDetector`
3. Implement required methods:
   - `get_motif_class_name()`
   - `get_patterns()`
   - `calculate_score()`
4. Add to `__init__.py` imports
5. Update `ModularMotifDetector` in `modular_scanner.py`

### Example New Detector
```python
from .base_detector import BaseMotifDetector

class NewMotifDetector(BaseMotifDetector):
    def get_motif_class_name(self) -> str:
        return "New_Motif"
    
    def get_patterns(self) -> Dict[str, List[Tuple]]:
        return {
            'pattern_group': [
                (r'PATTERN', 'ID', 'Name', 'Subclass', 10, 'scoring_method', 0.8, 'Description', 'Reference')
            ]
        }
    
    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        # Implement motif-specific scoring
        return score
```

This modular architecture successfully addresses all requirements while maintaining the robustness and functionality of the original system.