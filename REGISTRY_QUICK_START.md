# Quick Reference: Registry System

## Generate Registries

```bash
# Generate all registries (auto-discovers all detector classes)
python tools/generate_all_registries.py --out registry

# Force regeneration
python tools/generate_all_registries.py --out registry --force

# Verbose logging
python tools/generate_all_registries.py --out registry --verbose
```

## Load Registries in Python

```python
from utils.motif_patterns import get_hs_db_for_class, get_pattern_registry

# Load compiled DB or None if Hyperscan unavailable
db, id_to_ten, id_to_score = get_hs_db_for_class('ZDNA', 'registry')
print(f"Loaded {len(id_to_ten)} patterns")

# Load registry metadata
registry = get_pattern_registry('ZDNA', 'registry')
print(f"Source: {registry['meta']['source']}")
```

## Use with Scanner

```python
from utils.modular_scanner import ModularMotifDetector

# Scanner automatically loads registries from 'registry' directory
scanner = ModularMotifDetector()

# Or specify custom directory
scanner = ModularMotifDetector(registry_dir='/path/to/registry')

# Analyze sequence (uses preloaded registries for fast scanning)
results = scanner.analyze_sequence(sequence, "seq_name")
```

## Run Tests

```bash
# Test registry generation
python tests/test_registry_generation.py

# Test scanner integration
python tests/test_modular_scanner_registry_integration.py
```

## Add New Detector

1. Create detector with pattern table:
   ```python
   class MyDetector(BaseMotifDetector):
       TENMER_SCORE = {
           "AAAAAAAAAA": 10.0,
           # ... more patterns
       }
   ```

2. Generate registries (auto-discovers):
   ```bash
   python tools/generate_all_registries.py --out registry --force
   ```

3. Done! System automatically handles the rest.

## Environment Variable

```bash
# Set custom registry directory
export NBD_REGISTRY_DIR=/custom/path/to/registry
```

## File Structure

```
registry/
├── ZDNA_patterns.txt       # Plain text (one pattern per line)
├── ZDNA_registry.pkl       # Binary (fast loading)
├── ZDNA_registry.json      # JSON (human-readable)
└── ZDNA.hsdb              # Optional Hyperscan DB
```

## Key Features

- ✅ Automatic discovery of all detector classes
- ✅ Pattern validation (A/C/G/T/N only)
- ✅ Deterministic IDs (sorted patterns)
- ✅ Multiple formats (txt, pkl, json, hsdb)
- ✅ Graceful fallbacks (works without Hyperscan)
- ✅ Sharding for large pattern sets (>50k)
- ✅ In-memory caching for performance
- ✅ Comprehensive error handling
