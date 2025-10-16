# Tools Directory

This directory contains utility scripts for managing the NonBScanner system.

## Available Tools

### generate_all_registries.py (Recommended)

**Comprehensive registry generator** for all motif detector classes with automatic discovery.

**Purpose:**
- Automatically discovers all detector classes with pattern tables in motif_detection/
- Creates registry files (patterns.txt, registry.pkl, registry.json)
- Compiles Hyperscan databases (.hsdb) if Hyperscan is available
- Supports sharding for large pattern sets (>50k patterns)
- Validates patterns for DNA alphabet compliance

**Usage:**

```bash
# Generate all registries (default output: ./registry)
python tools/generate_all_registries.py

# Specify custom output directory
python tools/generate_all_registries.py --out /path/to/output

# Force regeneration even if files exist
python tools/generate_all_registries.py --force

# Use verbose logging for debugging
python tools/generate_all_registries.py --verbose

# Use smaller shard size for memory-constrained systems
python tools/generate_all_registries.py --shard-size 30000

# Discover from different package
python tools/generate_all_registries.py --pkg custom_detectors
```

**Auto-Discovery:**
The generator automatically finds detector classes by:
1. Scanning all modules in the specified package (default: motif_detection)
2. Looking for classes ending with "Detector"
3. Checking for pattern table attributes:
   - `TENMER_SCORE`, `TENMER_LOG2`, `TENMER_TABLE`
   - `PATTERN_TABLE`, `PATTERNS`, `PATTERN_DICT`
4. Extracting and validating patterns
5. Generating registries for each discovered class

**Output Files (per class):**
- `<CLASS>_patterns.txt` - Plain text list of patterns (one per line)
- `<CLASS>_registry.pkl` - Pickle serialization of registry (fast loading)
- `<CLASS>_registry.json` - JSON serialization of registry (human-readable)
- `<CLASS>.hsdb` - Compiled Hyperscan database (only if Hyperscan supports serialization)
- `<CLASS>_shard_N.hsdb` - Sharded DBs for large pattern sets (N=0,1,2,...)

**Registry Structure:**

Each registry contains:
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
    "detector_class": "ZDNADetector",
    "note": "Additional metadata if applicable"
  }
}
```

**Features:**
- ✓ Deterministic pattern ordering (case-insensitive alphabetical sort)
- ✓ Sequential ID assignment starting from 0
- ✓ Pattern validation (A/C/G/T/N only)
- ✓ Score normalization (converted to float)
- ✓ Comprehensive logging and error handling
- ✓ Graceful fallback when Hyperscan unavailable
- ✓ Sharding for memory-constrained compilation

**Dependencies:**
- Python 3.7+
- Access to motif_detection module
- Optional: hyperscan (for .hsdb compilation)

**Example Output:**

```
======================================================================
NON-B DNA SCANNER - REGISTRY GENERATOR
======================================================================
Output directory: registry
Package: motif_detection
Shard size: 50000
Hyperscan available: True

Discovering detector modules in 'motif_detection' package...
  ✓ APhilic: 208 patterns from APhilicDetector.TENMER_LOG2
  ✓ ZDNA: 126 patterns from ZDNADetector.TENMER_SCORE
Discovery complete: found 2 detector classes with pattern tables

======================================================================
GENERATING REGISTRIES
======================================================================
Generating registry for APhilic...
  ✓ Saved patterns text: registry/APhilic_patterns.txt
  ✓ Saved registry: registry/APhilic_registry.pkl, registry/APhilic_registry.json
  ✓ Saved Hyperscan DB: registry/APhilic.hsdb
✓ Registry generation complete for APhilic

Generating registry for ZDNA...
  ✓ Saved patterns text: registry/ZDNA_patterns.txt
  ✓ Saved registry: registry/ZDNA_registry.pkl, registry/ZDNA_registry.json
  ✓ Saved Hyperscan DB: registry/ZDNA.hsdb
✓ Registry generation complete for ZDNA

======================================================================
COMPLETE: Generated 2/2 registries
======================================================================
```

---

### generate_class_hsdb.py (Legacy - Deprecated)

**Legacy generator** for Z-DNA and A-philic 10-mers only.

**⚠️ Deprecated:** Use `generate_all_registries.py` instead for automatic discovery of all detector classes.

This legacy tool manually specifies Z-DNA and A-philic detectors and requires updates when adding new detector classes.

**Usage:**
```bash
python tools/generate_class_hsdb.py --out registry
```

---

## Adding New Detector Classes

To add support for a new detector class with pattern tables:

1. **Create your detector class with a pattern table:**
   
   ```python
   # In motif_detection/my_detector.py
   from .base_detector import BaseMotifDetector
   
   class MyDetector(BaseMotifDetector):
       """My custom motif detector"""
       
       # Define pattern table - use any of these names:
       # TENMER_SCORE, TENMER_LOG2, TENMER_TABLE,
       # PATTERN_TABLE, PATTERNS, or PATTERN_DICT
       TENMER_SCORE = {
           "AAAAAAAAAA": 10.0,
           "CCCCCCCCCC": 20.0,
           "GGGGGGGGGG": 30.0,
           # ... more patterns
       }
       
       def detect_motifs(self, sequence, sequence_name=""):
           # Implementation
           pass
   ```

2. **Run the generator** (automatically discovers your new detector):
   
   ```bash
   python tools/generate_all_registries.py --out registry --force
   ```
   
   The generator will automatically:
   - Discover your `MyDetector` class
   - Extract the `TENMER_SCORE` table
   - Generate all registry files
   - Create Hyperscan DB if available

3. **Update scanner detector map** (optional, for HS_DB_INFO attribute):
   
   If you want the scanner to set `HS_DB_INFO` on your detector class:
   
   ```python
   # In utils/modular_scanner.py, update detector_map in _preload_detector_dbs():
   detector_map = {
       "ZDNA": ("z_dna", ZDNADetector),
       "APhilic": ("a_philic", APhilicDetector),
       "My": ("my_class", MyDetector),  # Add this line
   }
   ```
   
   Replace "My" with your class name (without "Detector" suffix) to match the registry filename.

4. **Verify with tests:**
   
   ```bash
   python tests/test_registry_generation.py
   python tests/test_modular_scanner_registry_integration.py
   ```

That's it! The system handles everything else automatically.

## Testing

Test the registry generation system:

```bash
# Test registry generation (discovery, validation, file creation)
python tests/test_registry_generation.py

# Test scanner integration (automatic loading, fallbacks)
python tests/test_modular_scanner_registry_integration.py
```

## Troubleshooting

**No detectors discovered**
- Check that your detector class name ends with "Detector"
- Verify the pattern table attribute uses one of the recognized names
- Ensure the table is a non-empty dict
- Run with `--verbose` to see detailed discovery logs

**ImportError: Cannot import detector classes**
- Make sure you're running from the repository root
- Or set PYTHONPATH: `export PYTHONPATH=/path/to/NonBScanner`

**Hyperscan not available**
- The script will skip .hsdb generation but still create text/json/pkl registries
- These are sufficient for the system to work (with pure-Python fallback)
- Install Hyperscan if you want compiled databases: `pip install hyperscan`

**Pattern validation warnings**
- Patterns with invalid characters (not A/C/G/T/N) are logged
- These patterns are still included in the registry for debugging
- Review warnings to ensure patterns are correct

**Memory issues during compilation**
- Use a smaller `--shard-size` (e.g., `--shard-size 30000`)
- This creates multiple smaller DB files instead of one large one

## Development Notes

## Development Notes

**Pattern ID Assignment:**
- IDs are assigned sequentially (0, 1, 2, ...)
- Patterns are sorted alphabetically (case-insensitive) before ID assignment
- This ensures consistent IDs across regenerations
- Important for maintaining compatibility with saved Hyperscan DBs

**Pattern Validation:**
- Only A, C, G, T, and N are considered valid DNA characters
- Invalid patterns generate warnings but are still included
- Patterns are automatically uppercased for consistency

**Hyperscan Compilation:**
- Patterns are compiled as literal strings (exact match only)
- IDs correspond to array indices in the registry
- Case-insensitive matching is used
- Sharding prevents memory issues for large pattern sets

**Registry Format:**
- Pickle format is preferred for loading (fastest)
- JSON format is for human inspection and debugging
- Text format can be used for external tools or verification
- All three formats contain the same data

**Caching:**
- `utils/motif_patterns.py` caches loaded registries in memory
- `utils/load_hsdb.py` handles DB loading and caching
- Scanner maintains hsdb_map for quick access
- Reload by restarting the Python process
