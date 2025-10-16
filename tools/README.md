# Tools Directory

This directory contains utility scripts for managing the NonBScanner system.

## Available Tools

### generate_class_hsdb.py

Generator for creating Hyperscan pattern registries and compiled databases.

**Purpose:**
- Extracts 10-mer patterns from detector classes (Z-DNA, A-philic)
- Creates registry files (patterns.txt, registry.pkl, registry.json)
- Compiles Hyperscan databases (.hsdb) if Hyperscan is available

**Usage:**

```bash
# Generate all registries (default output: ./registry)
python tools/generate_class_hsdb.py

# Specify custom output directory
python tools/generate_class_hsdb.py --out /path/to/output

# Generate only Z-DNA registry
python tools/generate_class_hsdb.py --skip-aphilic

# Generate only A-philic registry
python tools/generate_class_hsdb.py --skip-zdna
```

**Output Files (per class):**
- `<CLASS>_patterns.txt` - Plain text list of 10-mer patterns (one per line)
- `<CLASS>_registry.pkl` - Pickle serialization of registry (fast loading)
- `<CLASS>_registry.json` - JSON serialization of registry (human-readable)
- `<CLASS>.hsdb` - Compiled Hyperscan database (only if Hyperscan supports serialization)

**Registry Structure:**

Each registry contains:
```json
{
  "class": "ZDNA",
  "generated_at": "2024-10-16T10:48:23.456Z",
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
    "source": "ZDNADetector.TENMER_SCORE",
    "note": "..."
  }
}
```

**Dependencies:**
- Python 3.6+
- Access to motif_detection module
- Optional: hyperscan (for .hsdb compilation)

**Example Output:**

```
[STEP] Building Z-DNA registry from ZDNADetector.TENMER_SCORE
[OK] saved patterns text -> registry/ZDNA_patterns.txt
[OK] saved registry (pickle/json) -> registry/ZDNA_registry.pkl, registry/ZDNA_registry.json
[OK] saved serialized Hyperscan DB -> registry/ZDNA.hsdb
[STEP] Building A-philic registry from APhilicDetector.TENMER_LOG2
[OK] saved patterns text -> registry/APhilic_patterns.txt
[OK] saved registry (pickle/json) -> registry/APhilic_registry.pkl, registry/APhilic_registry.json
[OK] saved serialized Hyperscan DB -> registry/APhilic.hsdb
[DONE] registry generation complete.
```

## Adding New Detector Classes

To add support for a new detector class:

1. **Update the generator script:**
   ```python
   # In generate_class_hsdb.py main():
   
   # Import the new detector
   from motif_detection.new_detector import NewDetector
   
   # Add generation step
   if not args.skip_newclass:
       print("[STEP] Building NewClass registry from NewDetector.PATTERN_TABLE")
       table = getattr(NewDetector, "PATTERN_TABLE", None)
       if table:
           generate_class_db(outdir, "NewClass", table, 
                           extra_meta={"source": "NewDetector.PATTERN_TABLE"})
   ```

2. **Add command-line flag:**
   ```python
   p.add_argument("--skip-newclass", action="store_true", 
                  help="skip NewClass generation")
   ```

3. **Update modular_scanner.py:**
   ```python
   # In ModularMotifDetector._preload_detector_dbs():
   candidate_classes = ["ZDNA", "APhilic", "NewClass"]
   
   # In detector_map:
   detector_map = {
       "ZDNA": ("z_dna", ZDNADetector),
       "APhilic": ("a_philic", APhilicDetector),
       "NewClass": ("new_class", NewDetector)
   }
   ```

## Troubleshooting

**ImportError: Cannot import detector classes**
- Make sure you're running from the repository root
- Or set PYTHONPATH: `export PYTHONPATH=/path/to/NonBScanner`

**Hyperscan not available**
- The script will skip .hsdb generation but still create text/json/pkl registries
- These are sufficient for the system to work (with pure-Python fallback)

**Detector table not found**
- Check that the detector class has the expected attribute (e.g., TENMER_SCORE)
- Verify the import path is correct

## Development Notes

**Pattern ID Assignment:**
- IDs are assigned sequentially (0, 1, 2, ...)
- Patterns are sorted alphabetically before ID assignment (for determinism)
- This ensures consistent IDs across regenerations

**Hyperscan Compilation:**
- Patterns are compiled as literal strings (exact match only)
- IDs correspond to array indices in the registry
- Case-insensitive matching is enabled by default

**Registry Format:**
- Pickle format is preferred for loading (fastest)
- JSON format is for human inspection and debugging
- Text format can be used for external tools or verification
