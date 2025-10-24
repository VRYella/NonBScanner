"""
Helper to load regex pattern registries and compile them with Hyperscan.

This module provides utilities for loading and using regex-based pattern registries
(CurvedDNA, G4, IMotif) similar to how the 10-mer registries work for A-philic and Z-DNA.

API:
    load_registry_for_class(class_name: str, registry_dir: str) 
        -> (db, id_to_pattern, id_to_subclass, id_to_score)
      - db: hyperscan.Database instance (or None if hyperscan not available)
      - id_to_pattern: dict mapping id -> regex pattern (string)
      - id_to_subclass: dict mapping id -> subclass name (string)
      - id_to_score: dict mapping id -> score (float)
"""

import os
import pickle
import json
import logging

logger = logging.getLogger(__name__)

try:
    import hyperscan
    _HYPERSCAN_AVAILABLE = True
except Exception:
    hyperscan = None
    _HYPERSCAN_AVAILABLE = False


def _load_registry(registry_dir: str, class_name: str):
    """Load registry from .pkl or .json file."""
    pkl_path = os.path.join(registry_dir, f"{class_name}_registry.pkl")
    json_path = os.path.join(registry_dir, f"{class_name}_registry.json")
    
    if os.path.isfile(pkl_path):
        with open(pkl_path, "rb") as fh:
            return pickle.load(fh)
    if os.path.isfile(json_path):
        with open(json_path, "r") as fh:
            return json.load(fh)
    
    raise FileNotFoundError(f"No registry found for {class_name} in {registry_dir}")


def load_registry_for_class(class_name: str, registry_dir: str = "registry"):
    """
    Load regex pattern registry and compile with Hyperscan.
    
    Returns (db, id_to_pattern, id_to_subclass, id_to_score).
    - db is a hyperscan.Database instance (or None if hyperscan not available).
    - id_to_pattern: dict mapping integer id -> regex pattern string
    - id_to_subclass: dict mapping integer id -> subclass name
    - id_to_score: dict mapping integer id -> score value
    """
    # Load registry data
    reg = _load_registry(registry_dir, class_name)
    patterns = reg.get("patterns", [])
    
    # Extract mappings
    id_to_pattern = {int(p["id"]): p["pattern"] for p in patterns}
    id_to_subclass = {int(p["id"]): p.get("subclass", "unknown") for p in patterns}
    id_to_score = {int(p["id"]): float(p.get("score", 0.0)) for p in patterns}
    
    # Compile Hyperscan database if available
    db = None
    if _HYPERSCAN_AVAILABLE:
        try:
            # Try to load pre-compiled DB first
            hsdb_path = os.path.join(registry_dir, f"{class_name}.hsdb")
            if os.path.isfile(hsdb_path):
                try:
                    db = hyperscan.Database()
                    with open(hsdb_path, "rb") as fh:
                        raw = fh.read()
                    db.deserialize(raw)
                    logger.info(f"Loaded serialized Hyperscan DB for {class_name}")
                except Exception as e:
                    logger.warning(f"Failed to deserialize {class_name}.hsdb: {e}")
                    db = None
            
            # If no pre-compiled DB, compile from patterns
            if db is None:
                expressions = []
                ids = []
                flags = []
                
                for i in sorted(id_to_pattern.keys()):
                    expressions.append(id_to_pattern[i].encode("ascii"))
                    ids.append(i)
                    # Use CASELESS and DOTALL flags for DNA matching
                    flags.append(hyperscan.HS_FLAG_CASELESS | hyperscan.HS_FLAG_DOTALL)
                
                db = hyperscan.Database()
                db.compile(
                    expressions=expressions,
                    ids=ids,
                    elements=len(expressions),
                    flags=flags
                )
                logger.info(f"Compiled Hyperscan DB for {class_name} from {len(expressions)} patterns")
                
        except Exception as e:
            logger.error(f"Hyperscan compilation failed for {class_name}: {e}")
            db = None
    else:
        logger.debug(f"Hyperscan not available for {class_name}; use pure-Python matching")
    
    return db, id_to_pattern, id_to_subclass, id_to_score


# Cache compiled databases to avoid recompilation
_REGISTRY_CACHE = {}


def get_cached_registry(class_name: str, registry_dir: str = "registry"):
    """
    Get cached registry or load and compile it.
    Returns (db, id_to_pattern, id_to_subclass, id_to_score).
    """
    cache_key = f"{registry_dir}/{class_name}"
    
    if cache_key not in _REGISTRY_CACHE:
        _REGISTRY_CACHE[cache_key] = load_registry_for_class(class_name, registry_dir)
    
    return _REGISTRY_CACHE[cache_key]


def scan_with_registry(class_name: str, sequence: str, registry_dir: str = "registry"):
    """
    Scan a sequence using a registry's compiled Hyperscan DB.
    
    Returns list of (start, end, pattern_id, subclass) matches.
    Falls back to pure-Python regex matching if Hyperscan not available.
    """
    db, id_to_pattern, id_to_subclass, id_to_score = get_cached_registry(class_name, registry_dir)
    matches = []
    
    if db is not None and _HYPERSCAN_AVAILABLE:
        # Use Hyperscan for fast matching
        def on_match(pattern_id, start, end, flags, context):
            subclass = id_to_subclass.get(pattern_id, "unknown")
            matches.append((start, end, pattern_id, subclass))
        
        try:
            db.scan(sequence.encode(), match_event_handler=on_match)
        except Exception as e:
            logger.error(f"Hyperscan scan failed for {class_name}: {e}")
            return []
    else:
        # Fallback to pure-Python regex matching
        import re
        for pattern_id in sorted(id_to_pattern.keys()):
            pattern = id_to_pattern[pattern_id]
            subclass = id_to_subclass[pattern_id]
            
            try:
                compiled_pattern = re.compile(pattern, re.IGNORECASE)
                for match in compiled_pattern.finditer(sequence):
                    matches.append((match.start(), match.end(), pattern_id, subclass))
            except re.error as e:
                logger.warning(f"Invalid regex pattern {pattern_id} in {class_name}: {e}")
                continue
    
    # Sort by start position
    matches.sort(key=lambda x: x[0])
    return matches
