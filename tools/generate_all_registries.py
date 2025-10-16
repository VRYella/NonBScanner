#!/usr/bin/env python3
"""
Generate Hyperscan DBs & registries for ALL motif detector classes.

This script discovers all detector modules in motif_detection/ package,
extracts pattern tables, and generates registries with optional Hyperscan
database compilation.

Outputs (per class):
  - <outdir>/<CLASS>_patterns.txt (one pattern per line)
  - <outdir>/<CLASS>_registry.pkl (pickle with pattern mapping + metadata)
  - <outdir>/<CLASS>_registry.json (JSON for human inspection)
  - <outdir>/<CLASS>.hsdb (optional: serialized Hyperscan DB)

For large pattern sets (>50k), supports sharding into multiple .hsdb files.

Usage:
  python tools/generate_all_registries.py --out registry
  python tools/generate_all_registries.py --out registry --force
  python tools/generate_all_registries.py --out registry --shard-size 30000
"""
import os
import sys
import argparse
import pickle
import json
import re
import logging
import importlib
import pkgutil
from datetime import datetime, timezone
from typing import Dict, List, Tuple, Any, Optional

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s'
)
logger = logging.getLogger("nbdscanner.registry")

# Adjust path to import from parent directory
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Attempt to import Hyperscan
try:
    import hyperscan
    _HYPERSCAN_AVAILABLE = True
except Exception:
    hyperscan = None
    _HYPERSCAN_AVAILABLE = False

# DNA alphabet for validation
DNA_ALPHABET = set('ACGTN')


def discover_detector_modules(package_name: str = "motif_detection") -> Dict[str, Any]:
    """
    Discover all modules in the specified package and find detector classes.
    
    Returns:
        Dict mapping class name -> (module, detector_class, pattern_table, table_attr_name)
    """
    logger.info(f"Discovering detector modules in '{package_name}' package...")
    
    # Import the package
    try:
        package = importlib.import_module(package_name)
    except ImportError as e:
        logger.error(f"Cannot import package '{package_name}': {e}")
        return {}
    
    discovered = {}
    
    # Iterate through all modules in the package
    package_path = package.__path__
    for importer, modname, ispkg in pkgutil.iter_modules(package_path):
        if ispkg or modname.startswith('_'):
            continue
        
        full_modname = f"{package_name}.{modname}"
        try:
            module = importlib.import_module(full_modname)
        except Exception as e:
            logger.warning(f"Failed to import {full_modname}: {e}")
            continue
        
        # Look for detector classes and pattern tables
        for attr_name in dir(module):
            if attr_name.startswith('_'):
                continue
            
            attr = getattr(module, attr_name, None)
            
            # Check if this is a detector class (name ends with 'Detector')
            if isinstance(attr, type) and attr_name.endswith('Detector'):
                # Look for pattern table attributes
                pattern_attrs = [
                    'TENMER_SCORE', 'TENMER_LOG2', 'TENMER_TABLE',
                    'PATTERN_TABLE', 'PATTERNS', 'PATTERN_DICT'
                ]
                
                for pat_attr in pattern_attrs:
                    if hasattr(attr, pat_attr):
                        table = getattr(attr, pat_attr)
                        if isinstance(table, dict) and len(table) > 0:
                            # Found a pattern table!
                            class_name = attr_name.replace('Detector', '')
                            discovered[class_name] = (module, attr, table, pat_attr)
                            logger.info(f"  ✓ {class_name}: {len(table)} patterns from {attr_name}.{pat_attr}")
                            break
    
    # Also check module-level pattern tables
    for importer, modname, ispkg in pkgutil.iter_modules(package_path):
        if ispkg or modname.startswith('_'):
            continue
        
        full_modname = f"{package_name}.{modname}"
        try:
            module = importlib.import_module(full_modname)
        except:
            continue
        
        # Check for module-level pattern tables
        pattern_attrs = [
            'TENMER_SCORE', 'TENMER_LOG2', 'TENMER_TABLE',
            'PATTERN_TABLE', 'PATTERNS', 'PATTERN_DICT'
        ]
        
        for pat_attr in pattern_attrs:
            if hasattr(module, pat_attr):
                table = getattr(module, pat_attr)
                if isinstance(table, dict) and len(table) > 0:
                    # Use module name as class name
                    class_name = modname.replace('_detector', '').replace('_', '').title()
                    if class_name not in discovered:
                        discovered[class_name] = (module, None, table, pat_attr)
                        logger.info(f"  ✓ {class_name}: {len(table)} patterns from module-level {pat_attr}")
                        break
    
    logger.info(f"Discovery complete: found {len(discovered)} detector classes with pattern tables")
    return discovered


def validate_pattern(pattern: str) -> Tuple[bool, Optional[str]]:
    """
    Validate that a pattern contains only valid DNA characters.
    
    Returns:
        (is_valid, warning_message)
    """
    pattern_upper = pattern.upper()
    invalid_chars = set(pattern_upper) - DNA_ALPHABET
    
    if invalid_chars:
        return False, f"Invalid characters: {invalid_chars}"
    
    return True, None


def normalize_pattern_table(table: Dict[str, Any]) -> List[Tuple[int, str, float]]:
    """
    Convert a pattern table to normalized format: [(id, pattern, score), ...]
    
    - Patterns are sorted lexicographically for deterministic IDs
    - IDs are sequential integers starting at 0
    - Scores are converted to float
    - Patterns are uppercased
    - Invalid patterns are logged but included with warning
    
    Returns:
        List of (id, pattern_str, score_float) tuples
    """
    items = []
    warnings = []
    
    # Sort patterns for determinism (uppercase first for case-insensitive sort)
    sorted_patterns = sorted(table.keys(), key=lambda x: x.upper())
    
    for idx, pattern in enumerate(sorted_patterns):
        pattern_upper = pattern.upper()
        score = float(table[pattern]) if table[pattern] is not None else 0.0
        
        # Validate pattern
        is_valid, warning = validate_pattern(pattern_upper)
        if not is_valid:
            warnings.append(f"Pattern '{pattern}' (ID {idx}): {warning}")
        
        items.append((idx, pattern_upper, score))
    
    if warnings:
        logger.warning(f"Pattern validation warnings ({len(warnings)} patterns):")
        for w in warnings[:10]:  # Show first 10
            logger.warning(f"  {w}")
        if len(warnings) > 10:
            logger.warning(f"  ... and {len(warnings) - 10} more")
    
    return items


def save_pattern_text(outdir: str, class_name: str, items: List[Tuple[int, str, float]]):
    """Save patterns as plain text (one per line)."""
    path = os.path.join(outdir, f"{class_name}_patterns.txt")
    with open(path, "w") as fh:
        for _, pattern, _ in items:
            fh.write(pattern + "\n")
    logger.info(f"  ✓ Saved patterns text: {path}")


def save_registry(outdir: str, class_name: str, items: List[Tuple[int, str, float]], 
                  extra_meta: Dict[str, Any]):
    """Save registry as both pickle and JSON."""
    registry = {
        "class": class_name,
        "generated_at": datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z'),
        "n_patterns": len(items),
        "patterns": [
            {"id": i, "tenmer": pattern, "score": score} 
            for (i, pattern, score) in items
        ],
        "meta": extra_meta,
    }
    
    # Save pickle
    pkl_path = os.path.join(outdir, f"{class_name}_registry.pkl")
    with open(pkl_path, "wb") as fh:
        pickle.dump(registry, fh)
    
    # Save JSON
    json_path = os.path.join(outdir, f"{class_name}_registry.json")
    with open(json_path, "w") as fh:
        json.dump(registry, fh, indent=2)
    
    logger.info(f"  ✓ Saved registry: {pkl_path}, {json_path}")


def compile_hyperscan_db(outdir: str, class_name: str, items: List[Tuple[int, str, float]], 
                         shard_size: int = 50000) -> bool:
    """
    Compile Hyperscan DB for given patterns.
    
    If number of patterns > shard_size, creates multiple sharded DBs:
      <class_name>_shard_0.hsdb, <class_name>_shard_1.hsdb, etc.
    
    Otherwise creates single DB: <class_name>.hsdb
    
    Args:
        outdir: Output directory
        class_name: Class name for output files
        items: List of (id, pattern, score) tuples
        shard_size: Maximum patterns per shard (default 50000)
    
    Returns:
        True if compilation succeeded
    """
    if not _HYPERSCAN_AVAILABLE:
        logger.info(f"  ⊘ Hyperscan not available — skipping DB compile for {class_name}")
        return False
    
    n_patterns = len(items)
    
    # Determine if sharding is needed
    if n_patterns > shard_size:
        logger.info(f"  ⚠ {n_patterns} patterns exceeds shard size {shard_size}, creating shards...")
        n_shards = (n_patterns + shard_size - 1) // shard_size
        
        for shard_idx in range(n_shards):
            start_idx = shard_idx * shard_size
            end_idx = min(start_idx + shard_size, n_patterns)
            shard_items = items[start_idx:end_idx]
            
            success = _compile_single_db(
                outdir, 
                f"{class_name}_shard_{shard_idx}", 
                shard_items
            )
            if not success:
                logger.warning(f"  ✗ Failed to compile shard {shard_idx} for {class_name}")
        
        return True
    else:
        # Single DB
        return _compile_single_db(outdir, class_name, items)


def _compile_single_db(outdir: str, db_name: str, items: List[Tuple[int, str, float]]) -> bool:
    """Helper to compile a single Hyperscan DB."""
    expressions = []
    ids = []
    
    for (i, pattern, score) in items:
        # For 10-mers and short patterns, use as literal patterns
        # Hyperscan requires patterns as bytes
        expressions.append(pattern.encode("ascii"))
        ids.append(i)
    
    try:
        db = hyperscan.Database()
        logger.debug(f"    Compiling {len(expressions)} patterns for {db_name}...")
        db.compile(expressions=expressions, ids=ids, elements=len(expressions))
        
        # Try to serialize if supported
        try:
            serial = db.serialize()
            hsdb_path = os.path.join(outdir, f"{db_name}.hsdb")
            with open(hsdb_path, "wb") as fh:
                fh.write(serial)
            logger.info(f"  ✓ Saved Hyperscan DB: {hsdb_path}")
            return True
        except Exception as e:
            logger.warning(f"  ⚠ DB compiled but serialize() failed: {e}")
            logger.info(f"    Runtime will compile from patterns.txt instead")
            return True  # compilation succeeded even if serialization failed
            
    except Exception as e:
        logger.error(f"  ✗ Hyperscan compilation failed for {db_name}: {e}")
        return False


def generate_class_registry(outdir: str, class_name: str, table: Dict[str, Any],
                            extra_meta: Dict[str, Any], shard_size: int = 50000):
    """
    Generate all registry files for a single detector class.
    
    Args:
        outdir: Output directory
        class_name: Name of the detector class
        table: Pattern table (pattern_str -> score)
        extra_meta: Metadata to include in registry
        shard_size: Maximum patterns per Hyperscan DB shard
    """
    logger.info(f"Generating registry for {class_name}...")
    
    # Normalize patterns
    items = normalize_pattern_table(table)
    
    # Save text file
    save_pattern_text(outdir, class_name, items)
    
    # Save registry files
    save_registry(outdir, class_name, items, extra_meta)
    
    # Compile Hyperscan DB
    compile_hyperscan_db(outdir, class_name, items, shard_size)
    
    logger.info(f"✓ Registry generation complete for {class_name}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate Hyperscan registries for all motif detector classes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate all registries (default output: ./registry)
  python tools/generate_all_registries.py
  
  # Specify custom output directory
  python tools/generate_all_registries.py --out /path/to/registry
  
  # Force regeneration even if files exist
  python tools/generate_all_registries.py --force
  
  # Use smaller shard size for memory-constrained systems
  python tools/generate_all_registries.py --shard-size 30000
  
  # Discover from different package
  python tools/generate_all_registries.py --pkg custom_detectors
        """
    )
    
    parser.add_argument(
        "--out", "-o", 
        default="registry",
        help="Output directory for registries and DBs (default: registry)"
    )
    parser.add_argument(
        "--pkg",
        default="motif_detection",
        help="Package name to scan for detectors (default: motif_detection)"
    )
    parser.add_argument(
        "--force", "-f",
        action="store_true",
        help="Force regeneration even if registry files exist"
    )
    parser.add_argument(
        "--shard-size",
        type=int,
        default=50000,
        help="Maximum patterns per Hyperscan DB shard (default: 50000)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose debug logging"
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Create output directory
    os.makedirs(args.out, exist_ok=True)
    
    logger.info("="*70)
    logger.info("NON-B DNA SCANNER - REGISTRY GENERATOR")
    logger.info("="*70)
    logger.info(f"Output directory: {args.out}")
    logger.info(f"Package: {args.pkg}")
    logger.info(f"Shard size: {args.shard_size}")
    logger.info(f"Hyperscan available: {_HYPERSCAN_AVAILABLE}")
    logger.info("")
    
    # Discover all detector classes
    discovered = discover_detector_modules(args.pkg)
    
    if not discovered:
        logger.warning("No detector classes with pattern tables found!")
        logger.info("Make sure detector classes have one of these attributes:")
        logger.info("  TENMER_SCORE, TENMER_LOG2, TENMER_TABLE,")
        logger.info("  PATTERN_TABLE, PATTERNS, PATTERN_DICT")
        return 1
    
    logger.info("")
    logger.info("="*70)
    logger.info("GENERATING REGISTRIES")
    logger.info("="*70)
    
    # Generate registry for each discovered class
    success_count = 0
    for class_name, (module, detector_cls, table, attr_name) in discovered.items():
        try:
            # Check if files exist and skip if not forcing
            pkl_path = os.path.join(args.out, f"{class_name}_registry.pkl")
            if os.path.exists(pkl_path) and not args.force:
                logger.info(f"⊘ Skipping {class_name} (registry exists, use --force to regenerate)")
                continue
            
            # Prepare metadata
            source = f"{detector_cls.__name__}.{attr_name}" if detector_cls else f"{module.__name__}.{attr_name}"
            meta = {
                "source_module": module.__name__,
                "source_attribute": attr_name,
                "source": source,
                "detector_class": detector_cls.__name__ if detector_cls else None,
            }
            
            # Add notes for specific attribute types
            if "LOG2" in attr_name:
                meta["note"] = "Values are log2 odds ratios"
            
            # Generate registry
            generate_class_registry(args.out, class_name, table, meta, args.shard_size)
            success_count += 1
            logger.info("")
            
        except Exception as e:
            logger.error(f"✗ Failed to generate registry for {class_name}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    logger.info("="*70)
    logger.info(f"COMPLETE: Generated {success_count}/{len(discovered)} registries")
    logger.info("="*70)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
