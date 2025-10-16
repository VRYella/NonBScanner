# Hyperscan Registry System

This directory contains precompiled pattern registries for high-performance motif detection using Hyperscan (when available).

## Directory Structure

```
registry/
├── ZDNA_patterns.txt      # Plain text list of Z-DNA 10-mer patterns
├── ZDNA_registry.pkl      # Pickle serialization of Z-DNA registry
├── ZDNA_registry.json     # JSON serialization of Z-DNA registry
├── APhilic_patterns.txt   # Plain text list of A-philic 10-mer patterns
├── APhilic_registry.pkl   # Pickle serialization of A-philic registry
└── APhilic_registry.json  # JSON serialization of A-philic registry
```

## Registry Contents

Each registry contains:
- **id**: Sequential integer ID for each pattern (for Hyperscan)
- **tenmer**: The 10-mer sequence pattern
- **score**: The scoring value for the pattern
- **metadata**: Source and generation information

### Z-DNA Registry
- **Patterns**: 126 Z-DNA 10-mer motifs
- **Source**: `ZDNADetector.TENMER_SCORE`
- **Score Range**: 50.0 - 63.0

### A-philic Registry
- **Patterns**: 208 A-philic 10-mer motifs
- **Source**: `APhilicDetector.TENMER_LOG2`
- **Score Type**: log2 odds ratios
- **Score Range**: varies (log2 values)

## Usage

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

```python
from utils.motif_patterns import get_hs_db_for_class, get_pattern_registry

# Load compiled Hyperscan DB (or None if Hyperscan not available)
db, id_to_ten, id_to_score = get_hs_db_for_class('ZDNA', 'registry')

# Load registry metadata
registry = get_pattern_registry('ZDNA', 'registry')
print(f"Loaded {registry['n_patterns']} patterns")
```

### Environment Variable

Set the `NBD_REGISTRY_DIR` environment variable to use a custom registry location:

```bash
export NBD_REGISTRY_DIR=/path/to/custom/registry
python app.py
```

## Regenerating Registries

To regenerate the registries (e.g., after updating detector pattern tables):

```bash
# Regenerate all registries
python tools/generate_class_hsdb.py --out registry

# Regenerate only Z-DNA
python tools/generate_class_hsdb.py --out registry --skip-aphilic

# Regenerate only A-philic
python tools/generate_class_hsdb.py --out registry --skip-zdna
```

## Hyperscan Binary DBs (.hsdb)

If Hyperscan is installed with serialization support, the generator will also create `.hsdb` files containing precompiled Hyperscan databases. These are excluded from version control (via `.gitignore`) and will be generated at runtime if needed.

## Fallback Behavior

The system gracefully handles missing Hyperscan:
1. If Hyperscan is not installed, the loader returns `None` for the DB object
2. Detectors automatically fall back to pure-Python matching
3. All pattern and score information is still available via `id_to_ten` and `id_to_score`

## Testing

Run the integration test suite to verify the registry system:

```bash
python test_registry_integration.py
```

This will test:
- Registry generation
- Registry loading
- Integration with motif_patterns module
- Integration with modular_scanner
- Detector class attributes

## Performance

When Hyperscan is available:
- **Pattern matching**: 10-100x faster than pure Python regex
- **Memory usage**: Minimal (patterns compiled to optimized state machine)
- **Startup time**: ~10ms per class for DB compilation (cached in memory)

Without Hyperscan:
- Falls back to pure Python exact string matching
- Still faster than regex for exact 10-mer matching
- No additional dependencies required
