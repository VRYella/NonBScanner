# NBDFinder Archives

This directory contains archived files that are no longer needed for the main NBDFinder functionality but are preserved for historical reference.

## Directory Structure

### `documentation/`
Process documentation files from the refactor:
- `ARCHITECTURE_REVAMP_COMPLETE.md` - Complete architecture revamp summary
- `HYPERSCAN_UPGRADE_SUMMARY.md` - Hyperscan performance upgrade documentation
- `MOTIFS_RESTRUCTURE.md` - Original motifs restructure documentation
- `REVAMP_COMPLETION_SUMMARY.md` - Revamp completion summary
- `VISUALIZATION_ENHANCEMENT_SUMMARY.md` - Visualization features enhancement

### `tests/`
Redundant test files moved to archives:
- `test_revamp_features.py` - REST API and integration tests (superseded by test_app_integration.py)
- `test_enhanced_features.py` - Enhanced features tests (functionality now in core tests)
- `test_ui_improvements.py` - UI improvement tests (functionality integrated into main app)
- `test_hyperscan_summary.py` - Simple hyperscan test (superseded by test_hyperscan_performance.py)
- `test_enhanced_visualizations_comprehensive.py` - Comprehensive visualization tests (superseded by test_enhanced_visualizations.py)
- `test_performance_optimization.py` - Performance optimization tests (functionality in core tests)

### Root level
- `motifs.py` - Original motifs implementation (replaced by all_motifs_refactored.py and modular architecture)

## What Remains Active

### Core Functionality
- `app.py` - Streamlit web interface
- `api.py` - REST API server  
- `all_motifs_refactored.py` - Unified orchestrator
- `regex_registry.py` - Central pattern registry
- `curated_examples.py` - NCBI test sequences

### Essential Tests
- `test_curated_examples.py` - Core functionality validation
- `test_app_integration.py` - App integration testing
- `test_comprehensive_motifs.py` - Comprehensive motif detection
- `test_enhanced_visualizations.py` - Visualization features
- `test_hyperscan_performance.py` - Performance validation
- `test_motifs.py` - Basic motif tests
- `test_orchestrator.py` - Orchestrator testing

### Modules
- `motifs/` - Modular motif detection implementations
- `utils.py` - Utility functions
- `classification_config.py` - Configuration and classification
- Other supporting modules

The refactor successfully consolidated 247 regex patterns into a central registry, implemented parallel processing, and preserved all scientific accuracy while improving maintainability.