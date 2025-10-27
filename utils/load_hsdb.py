"""
Helper to load Hyperscan DB or compile from patterns when needed.

API:
    load_db_for_class(class_name: str, registry_dir: str) -> (db, id_to_pattern, id_to_score)
      - db: hyperscan.Database instance (or None if hyperscan not available)
      - id_to_pattern: dict mapping id -> pattern string (tenmer for 10-mer, regex for others)
      - id_to_score: dict mapping id -> score (float)
      
    Supports both:
      - 10-mer patterns (ZDNA, APhilic): patterns have 'tenmer' key
      - Regex patterns (G4, IMotif, etc.): patterns have 'pattern' key
"""

import os
import pickle
import logging

logger = logging.getLogger(__name__)

try:
    import hyperscan
    _HYPERSCAN_AVAILABLE = True
except Exception:
    hyperscan = None
    _HYPERSCAN_AVAILABLE = False


def _load_registry(registry_dir: str, class_name: str):
    pkl_path = os.path.join(registry_dir, f"{class_name}_registry.pkl")
    json_path = os.path.join(registry_dir, f"{class_name}_registry.json")
    if os.path.isfile(pkl_path):
        with open(pkl_path, "rb") as fh:
            return pickle.load(fh)
    if os.path.isfile(json_path):
        import json
        with open(json_path, "r") as fh:
            return json.load(fh)
    raise FileNotFoundError(f"No registry found for {class_name} in {registry_dir}")


def load_db_for_class(class_name: str, registry_dir: str = "registry"):
    """
    Returns (db, id_to_pattern, id_to_score).
    db is a hyperscan.Database instance (or None if hyperscan not available).
    id_to_pattern and id_to_score are dicts mapping integer id -> pattern/tenmer / score.
    
    Handles both 10-mer patterns (with 'tenmer' key) and regex patterns (with 'pattern' key).
    """
    reg = _load_registry(registry_dir, class_name)
    patterns = reg.get("patterns", [])
    
    # Handle both 10-mer registries (tenmer) and regex registries (pattern)
    id_to_pattern = {}
    for p in patterns:
        pattern_id = int(p["id"])
        # Try 'tenmer' first (for 10-mer patterns like ZDNA, APhilic)
        if "tenmer" in p:
            id_to_pattern[pattern_id] = p["tenmer"]
        # Fall back to 'pattern' (for regex patterns like G4, IMotif, etc.)
        elif "pattern" in p:
            id_to_pattern[pattern_id] = p["pattern"]
        else:
            logger.warning(f"Pattern {pattern_id} in {class_name} has neither 'tenmer' nor 'pattern' key")
    
    id_to_score = {int(p["id"]): float(p.get("score", 0.0)) for p in patterns}

    db = None
    hsdb_path = os.path.join(registry_dir, f"{class_name}.hsdb")
    if _HYPERSCAN_AVAILABLE:
        try:
            db = hyperscan.Database()
            if os.path.isfile(hsdb_path):
                # Many hyperscan Python bindings provide deserialize()
                with open(hsdb_path, "rb") as fh:
                    raw = fh.read()
                try:
                    db.deserialize(raw)
                    logger.info(f"Loaded serialized DB for {class_name} from {hsdb_path}")
                except Exception:
                    # Fallback: compile from patterns
                    expressions = [id_to_pattern[i].encode("ascii") for i in sorted(id_to_pattern.keys())]
                    ids = sorted(id_to_pattern.keys())
                    db.compile(expressions=expressions, ids=ids, elements=len(expressions))
                    logger.warning(f"deserialize() failed; compiled DB for {class_name} from patterns")
            else:
                # No serialized DB; compile from patterns
                expressions = [id_to_pattern[i].encode("ascii") for i in sorted(id_to_pattern.keys())]
                ids = sorted(id_to_pattern.keys())
                db.compile(expressions=expressions, ids=ids, elements=len(expressions))
                logger.info(f"Compiled DB for {class_name} from patterns in {registry_dir}")
        except Exception as e:
            logger.error(f"Hyperscan operations failed: {e}")
            db = None
    else:
        logger.debug("Hyperscan not installed; returning None DB (use pure-Python matcher).")

    return db, id_to_pattern, id_to_score
