# Automatic Detector Discovery System

## Overview

The NonBScanner now includes an **automatic detector discovery system** that ensures all motif detectors are picked up for testing and demonstration without manual maintenance. This addresses the requirement to "ensure that the tool has ability to pick all motifs for testing/retesting."

## Key Features

✅ **Automatic Discovery**: All detector classes are automatically discovered from the `motif_detection` module  
✅ **Zero Maintenance**: Adding new detectors requires no changes to test scripts  
✅ **Centralized Registry**: Test sequences are managed in one location  
✅ **Backward Compatible**: All existing tests continue to work  
✅ **Type Safe**: Uses proper type hints and base class checking  

## Architecture

### Core Component: `utils/detector_registry.py`

This module provides the automatic discovery and test sequence management:

```python
from utils import detector_registry

# Get all detector classes automatically
detectors = detector_registry.get_all_detector_classes()
# Returns: {'APhilicDetector': <class>, 'ZDNADetector': <class>, ...}

# Get all detectors with their test sequences
for name, detector_class, test_sequences in detector_registry.get_all_detectors_with_test_sequences():
    detector = detector_class()
    for description, sequence in test_sequences:
        results = detector.detect_motifs(sequence, 'test')
```

### Updated Test Scripts

Both `test_all_motifs.py` and `demo_all_detectors.py` now use automatic discovery:

**Before** (Manual maintenance required):
```python
from motif_detection.a_philic_detector import APhilicDetector
from motif_detection.z_dna_detector import ZDNADetector
# ... 9 manual imports

# Manual test definition for each class
a_philic_tests = [...]
results = test_motif_class("A-philic DNA", APhilicDetector(), a_philic_tests)
# ... repeated for all 9 classes
```

**After** (Automatic discovery):
```python
from utils import detector_registry

# Automatically discover all detectors
for display_name, detector_class, test_sequences in detector_registry.get_all_detectors_with_test_sequences():
    detector = detector_class()
    results = test_motif_class(display_name, detector, test_sequences)
```

## Functions Available

### `get_all_detector_classes()`
Returns all detector classes found in the motif_detection module.

```python
detectors = detector_registry.get_all_detector_classes()
# Returns: Dict[str, Type[BaseMotifDetector]]
# Example: {'APhilicDetector': APhilicDetector, 'ZDNADetector': ZDNADetector, ...}
```

### `get_detector_display_name(detector_class)`
Gets a human-readable name for a detector.

```python
name = detector_registry.get_detector_display_name(APhilicDetector)
# Returns: "A-philic_DNA"
```

### `get_test_sequences_registry()`
Returns comprehensive test sequences for all detectors.

```python
test_sequences = detector_registry.get_test_sequences_registry()
# Returns: Dict[str, List[Tuple[str, str]]]
# Each detector has 3-4 test sequences
```

### `get_demo_sequences_registry()`
Returns focused demo sequences (shorter list of best examples).

```python
demo_sequences = detector_registry.get_demo_sequences_registry()
# Returns: Dict[str, List[Tuple[str, str]]]
# Each detector has 2-3 demo sequences
```

### `get_all_detectors_with_test_sequences()`
One-stop function that combines everything for testing.

```python
for display_name, detector_class, test_sequences in detector_registry.get_all_detectors_with_test_sequences():
    # display_name: Human-readable name like "A-philic_DNA"
    # detector_class: Class object that can be instantiated
    # test_sequences: List of (description, sequence) tuples
    pass
```

### `get_all_detectors_with_demo_sequences()`
Same as above but with demo sequences instead of full test suite.

## How It Works

### 1. Discovery Mechanism

The system uses Python's introspection to find all classes:

```python
import inspect
import importlib

# Import the motif_detection module
motif_detection = importlib.import_module('motif_detection')

# Find all subclasses of BaseMotifDetector
for name in dir(motif_detection):
    obj = getattr(motif_detection, name)
    if (inspect.isclass(obj) and 
        issubclass(obj, BaseMotifDetector) and 
        obj != BaseMotifDetector):
        # Found a detector!
        detector_classes[name] = obj
```

### 2. Test Sequence Registry

Test sequences are centrally managed in dictionaries:

```python
{
    'APhilicDetector': [
        ("Problem statement sequence", "AGGGGGGGGGAGGGGGGGGC"),
        ("Multiple G-runs", "AGGGGGGGGGTGGGGGGGGC"),
        ("C-rich (complement)", "CCCCCCCCCCCCCCCCCCCC"),
        ("Mixed GC-rich", "GGGGGCCCCCGGGGGCCCCC"),
    ],
    'ZDNADetector': [...],
    # ... etc
}
```

### 3. Automatic Integration

Test scripts now automatically:
1. Discover all detector classes
2. Match them with their test sequences
3. Run tests on all detectors
4. Report comprehensive results

## Benefits

### For Developers

✅ **Add new detectors easily**: Just create a new class inheriting from `BaseMotifDetector`  
✅ **No test script updates needed**: New detectors are automatically tested  
✅ **Centralized test data**: All test sequences in one file  
✅ **Easy maintenance**: Update sequences in one place  

### For Testing

✅ **Comprehensive coverage**: Ensures ALL detectors are tested  
✅ **No missed detectors**: Automatic discovery catches everything  
✅ **Consistent testing**: Same test framework for all detectors  
✅ **Easy verification**: Single command tests all detectors  

### For Documentation

✅ **Always up-to-date**: Detector count is dynamically determined  
✅ **Self-documenting**: Scripts show what detectors are available  
✅ **Clear output**: Shows which detectors were discovered  

## Usage Examples

### Running Comprehensive Tests

```bash
$ python3 test_all_motifs.py
================================================================================
COMPREHENSIVE NON-B DNA MOTIF DETECTOR TEST SUITE
================================================================================

Automatically discovered 9 detector classes
================================================================================

Testing A-philic_DNA
✓ PASS     Problem statement sequence
...

Total Classes: 9
Total Tests:   28
Total Passed:  20
Success Rate:  71.4%
```

### Running Demo Script

```bash
$ python3 demo_all_detectors.py
================================================================================
              NON-B DNA MOTIF DETECTOR COMPREHENSIVE VERIFICATION               
================================================================================

Automatically discovered 9 detector classes
================================================================================

All 9 Non-B DNA motif detection classes have been tested:
  1. ✓ A-philic_DNA
  2. ✓ Cruciform
  ...
```

### Using in Custom Scripts

```python
#!/usr/bin/env python3
import sys
import os
import importlib.util

# Load detector_registry
project_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, project_root)

spec = importlib.util.spec_from_file_location(
    "detector_registry",
    os.path.join(project_root, "utils", "detector_registry.py")
)
detector_registry = importlib.util.module_from_spec(spec)
spec.loader.exec_module(detector_registry)

# Use the registry
for name, detector_class, sequences in detector_registry.get_all_detectors_with_test_sequences():
    print(f"Testing {name}...")
    detector = detector_class()
    for desc, seq in sequences:
        motifs = detector.detect_motifs(seq, 'test')
        print(f"  {desc}: {len(motifs)} motifs found")
```

## Adding New Detectors

When adding a new detector class:

1. **Create the detector class** in `motif_detection/` inheriting from `BaseMotifDetector`
2. **Add to `__init__.py`** in the motif_detection module's `__all__` list
3. **Add test sequences** to `utils/detector_registry.py` in both registries:
   - `get_test_sequences_registry()` - Full test suite (3-4 sequences)
   - `get_demo_sequences_registry()` - Demo sequences (2-3 sequences)

That's it! The new detector will automatically be:
- Discovered by test scripts
- Included in comprehensive tests
- Shown in demo output
- Counted in statistics

## Testing the Registry

To verify the registry is working:

```bash
$ python3 utils/detector_registry.py

Discovered Detector Classes:
================================================================================
  APhilicDetector                → A-philic_DNA
  CruciformDetector              → Cruciform
  CurvedDNADetector              → Curved_DNA
  GQuadruplexDetector            → G-Quadruplex
  IMotifDetector                 → i-Motif
  RLoopDetector                  → R-Loop
  SlippedDNADetector             → Slipped_DNA
  TriplexDetector                → Triplex
  ZDNADetector                   → Z-DNA

Total detectors discovered: 9
```

## Current Status

✅ **All 9 detectors are automatically discovered**  
✅ **28 test sequences defined across all detectors**  
✅ **71.4% test pass rate (20/28 tests)**  
✅ **100% backward compatible with existing tests**  
✅ **Both test_all_motifs.py and demo_all_detectors.py updated**  

## Future Enhancements

Possible future improvements:

1. **Auto-generate test sequences**: Use AI/heuristics to suggest test sequences
2. **Dynamic sequence validation**: Validate that test sequences match detector patterns
3. **Performance benchmarking**: Automatically benchmark each detector
4. **Coverage reporting**: Track which patterns in each detector are tested
5. **Test generation**: Auto-generate edge cases and boundary tests

## Troubleshooting

### Import Issues

If you get `ModuleNotFoundError`, ensure the project root is in your Python path:

```python
import sys
import os
project_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, project_root)
```

### Registry Not Finding Detectors

1. Check that the detector class inherits from `BaseMotifDetector`
2. Verify the class is exported in `motif_detection/__init__.py`
3. Run `python3 utils/detector_registry.py` to see what's discovered

### Test Sequences Not Matching

Update sequences in `utils/detector_registry.py`:
- Use `get_test_sequences_registry()` for comprehensive testing
- Use `get_demo_sequences_registry()` for demonstrations

## Conclusion

The automatic detector discovery system ensures that **all motifs are picked up for testing and retesting** without manual intervention. This satisfies the issue requirement and provides a maintainable, scalable solution for the NonBScanner project.
