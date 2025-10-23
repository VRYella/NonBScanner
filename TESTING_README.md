# Testing and Verification Suite

This directory contains test scripts for the NonBScanner Non-B DNA motif detection system.

## Available Tests

The `tests/` directory contains the following test scripts:

### 1. test_registry_generation.py
- **Purpose**: Test pattern registry generation system
- **Tests**: Pattern table extraction, validation, and registry file generation
- **Run**: `python3 tests/test_registry_generation.py`

### 2. test_modular_scanner_registry_integration.py  
- **Purpose**: Test scanner integration with pattern registries
- **Tests**: Detector discovery, registry loading, and scanner functionality
- **Run**: `python3 tests/test_modular_scanner_registry_integration.py`

## Quick Start

To run the tests:

```bash
# Run registry generation tests
python3 tests/test_registry_generation.py

# Run scanner integration tests
python3 tests/test_modular_scanner_registry_integration.py
```

## Test Coverage

The test suite covers:

- ✓ Automatic detector discovery
- ✓ Pattern table extraction and validation
- ✓ Registry file generation (txt, pkl, json formats)
- ✓ Scanner initialization with registries
- ✓ Motif detection across all detector classes
- ✓ Fallback behavior when Hyperscan is unavailable

## Supported Detector Classes

The tests verify functionality for all 9 motif detection classes:

1. ✓ A-philic DNA    - A-rich protein binding sites
2. ✓ Z-DNA           - Left-handed double helix
3. ✓ i-Motif         - C-rich structures
4. ✓ G-Quadruplex    - Four-stranded G-rich structures
5. ✓ Triplex         - Three-stranded DNA
6. ✓ Cruciform       - Inverted repeats
7. ✓ Slipped DNA     - Tandem repeats (STRs)
8. ✓ Curved DNA      - A-tract mediated bending
9. ✓ R-Loop          - RNA-DNA hybrid sites

## Regenerating Pattern Registries

If you need to regenerate pattern registries (e.g., after updating detector pattern tables):

```bash
# Regenerate all registries
python tools/generate_all_registries.py --out registry

# Force regeneration even if files exist
python tools/generate_all_registries.py --out registry --force
```

See `registry/README.md` for more details on the registry system.
