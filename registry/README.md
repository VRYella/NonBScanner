# Hyperscan Registry System

This directory contains precompiled pattern registries for high-performance motif detection using Hyperscan (when available).

## Overview

The registry system provides:
- **Automatic discovery** of all detector classes with pattern tables
- **Multiple output formats**: plain text, pickle, JSON, and optional Hyperscan binary DBs
- **Deterministic pattern ordering** for stable IDs across regenerations
- **Pattern validation** ensuring only valid DNA characters (A/C/G/T/N)
- **Sharding support** for large pattern sets (>50k patterns)
- **Graceful fallbacks** when Hyperscan is not available

## Directory Structure

```
registry/
├── <CLASS>_patterns.txt      # Plain text list of patterns (one per line)
├── <CLASS>_registry.pkl       # Pickle serialization of registry (fast loading)
├── <CLASS>_registry.json      # JSON serialization (human-readable)
└── <CLASS>.hsdb               # Optional: compiled Hyperscan database
```

For large pattern sets (>50k), the system creates sharded DBs:
```
registry/
├── <CLASS>_shard_0.hsdb
├── <CLASS>_shard_1.hsdb
└── ...
```

## Supported Detector Classes

The system automatically discovers all detector classes with pattern tables.
Currently supported:

### 10-mer Exact Match Registries (Hyperscan)
- **ZDNA**: 126 Z-DNA 10-mer motifs (from `ZDNADetector.TENMER_SCORE`)
- **APhilic**: 208 A-philic 10-mer motifs (from `APhilicDetector.TENMER_LOG2`)

### Regex Pattern Registries (Hyperscan)
- **CurvedDNA**: 44 curved DNA patterns (2 local + 42 global curvature patterns)
- **G4**: 7 G-quadruplex pattern classes (canonical, relaxed, long-loop, bulged, multimeric, imperfect, g-triplex)
- **IMotif**: 7 i-motif patterns (1 canonical i-motif + 6 HUR AC-motif patterns)

Additional detectors can be added by implementing pattern tables with any of these attributes:
- `TENMER_SCORE`, `TENMER_LOG2`, `TENMER_TABLE`
- `PATTERN_TABLE`, `PATTERNS`, `PATTERN_DICT`

## Registry Contents

Each registry contains:
- **id**: Sequential integer ID for each pattern (deterministic, based on sorted pattern order)
- **tenmer**: The pattern sequence (uppercased)
- **score**: The scoring value for the pattern (float)
- **metadata**: Source module, attribute name, generation timestamp, and notes

### Registry Format (JSON example)

```json
{
  "class": "ZDNA",
  "generated_at": "2024-10-16T12:00:00Z",
  "n_patterns": 126,
  "patterns": [
    {
      "id": 0,
      "tenmer": "AACGCGCGCG",
      "score": 50.25
    },
    ...
  ],
  "meta": {
    "source_module": "motif_detection.z_dna_detector",
    "source_attribute": "TENMER_SCORE",
    "source": "ZDNADetector.TENMER_SCORE",
    "detector_class": "ZDNADetector"
  }
}
```

### Automatic Loading (Recommended)

The registries are automatically loaded by `ModularMotifDetector` when initialized:

```python
from utils.modular_scanner import ModularMotifDetector

# Registries are loaded automatically from default location
scanner = ModularMotifDetector()

# Or specify a custom registry directory
scanner = ModularMotifDetector(registry_dir='/path/to/custom/registry')
```

### Manual Loading

You can also load registries manually:

**For 10-mer registries (A-philic, Z-DNA):**
```python
from utils.motif_patterns import get_hs_db_for_class, get_pattern_registry

# Load compiled Hyperscan DB (or None if Hyperscan not available)
db, id_to_ten, id_to_score = get_hs_db_for_class('ZDNA', 'registry')

# Load registry metadata
registry = get_pattern_registry('ZDNA', 'registry')
print(f"Loaded {registry['n_patterns']} patterns")
```

**For regex pattern registries (CurvedDNA, G4, IMotif):**
```python
from utils.load_regex_registry import scan_with_registry, get_cached_registry

# Scan a sequence using a registry
matches = scan_with_registry('G4', sequence)
for start, end, pattern_id, subclass in matches:
    print(f"Match at {start}-{end}: {subclass}")

# Or load the registry manually
db, id_to_pattern, id_to_subclass, id_to_score = get_cached_registry('G4', 'registry')
```

### Environment Variable

Set the `NBD_REGISTRY_DIR` environment variable to use a custom registry location:

```bash
export NBD_REGISTRY_DIR=/path/to/custom/registry
python app.py
```

## Usage

## Regenerating Registries

To regenerate the registries (e.g., after updating detector pattern tables):

**For all registries (10-mer and regex):**
```bash
# Regenerate 10-mer registries (A-philic, Z-DNA)
python tools/generate_all_registries.py --out registry

# Regenerate regex pattern registries (CurvedDNA, G4, IMotif)
python tools/generate_pattern_registries.py

# Or regenerate both types
python tools/generate_all_registries.py --out registry && python tools/generate_pattern_registries.py
```

**For 10-mer registries only:**
```bash
# Regenerate all 10-mer registries
python tools/generate_all_registries.py --out registry

# Force regeneration even if files exist
python tools/generate_all_registries.py --out registry --force

# Use verbose logging
python tools/generate_all_registries.py --out registry --verbose

# Use smaller shard size for memory-constrained systems
python tools/generate_all_registries.py --out registry --shard-size 30000

# Discover from different package
python tools/generate_all_registries.py --out registry --pkg custom_detectors
```

The generator automatically discovers all detector classes with pattern tables in the specified package.

## Hyperscan Binary DBs (.hsdb)

If Hyperscan is installed with serialization support, the generator will also create `.hsdb` files containing precompiled Hyperscan databases. These are excluded from version control (via `.gitignore`) and will be generated at runtime if needed.

For pattern sets exceeding the shard size (default: 50,000), multiple sharded DBs are created:
- `<CLASS>_shard_0.hsdb`
- `<CLASS>_shard_1.hsdb`
- etc.

This prevents memory issues during compilation and allows parallel scanning.

## Pattern Validation

The generator validates all patterns to ensure they contain only valid DNA characters (A, C, G, T, N). Patterns with invalid characters are logged with warnings but still included in the registry for debugging purposes.

## Fallback Behavior

The system gracefully handles missing Hyperscan:
1. If Hyperscan is not installed, the loader returns `None` for the DB object
2. Detectors automatically fall back to pure-Python matching
3. All pattern and score information is still available via `id_to_ten` and `id_to_score`

## Testing

Run the test suites to verify the registry system:

```bash
# Test registry generation
python tests/test_registry_generation.py

# Test scanner integration
python tests/test_modular_scanner_registry_integration.py
```

These tests verify:
- Detector discovery and pattern table extraction
- Pattern validation and normalization  
- Deterministic ID assignment
- Registry file generation (txt, pkl, json)
- Loader compatibility
- Scanner integration and automatic registry discovery
- Fallback behavior when registries are missing
- Environment variable override (NBD_REGISTRY_DIR)

## Performance

When Hyperscan is available:
- **Pattern matching**: 10-100x faster than pure Python regex
- **Memory usage**: Minimal (patterns compiled to optimized state machine)
- **Startup time**: ~10ms per class for DB compilation (cached in memory)
- **Sharding**: Prevents memory issues for large pattern sets

Without Hyperscan:
- Falls back to pure Python exact string matching
- Still functional, just slower
- No additional dependencies required

## Adding New Detectors

To add support for a new detector class:

1. **Define pattern table in your detector**:
   ```python
   class MyDetector(BaseMotifDetector):
       TENMER_SCORE = {
           "AAAAAAAAAA": 10.0,
           "CCCCCCCCCC": 20.0,
           # ... more patterns
       }
   ```

2. **Run the generator** (auto-discovers new detectors):
   ```bash
   python tools/generate_all_registries.py --out registry --force
   ```

3. **Update scanner detector map** (if needed for HS_DB_INFO):
   ```python
   # In utils/modular_scanner.py, add to detector_map:
   detector_map = {
       "ZDNA": ("z_dna", ZDNADetector),
       "APhilic": ("a_philic", APhilicDetector),
       "MyClass": ("my_class", MyDetector)  # Add this
   }
   ```

The system will automatically discover and process your detector's patterns.
