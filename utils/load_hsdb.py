"""
Helper to load Hyperscan DB or compile from patterns when needed.

API:
    load_db_for_class(class_name: str, registry_dir: str) -> (db, id_to_ten, id_to_score)
      - db: hyperscan.Database instance (or None if hyperscan not available)
      - id_to_ten: dict mapping id -> tenmer (string)
      - id_to_score: dict mapping id -> score (float)
"""

import os
import pickle

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
    Returns (db, id_to_ten, id_to_score).
    db is a hyperscan.Database instance (or None if hyperscan not available).
    id_to_ten and id_to_score are dicts mapping integer id -> tenmer / score.
    """
    reg = _load_registry(registry_dir, class_name)
    patterns = reg.get("patterns") or reg.get("patterns", [])
    id_to_ten = {int(p["id"]): p["tenmer"] for p in patterns}
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
                    print(f"[OK] Loaded serialized DB for {class_name} from {hsdb_path}")
                except Exception:
                    # Fallback: compile from patterns
                    expressions = [id_to_ten[i].encode("ascii") for i in sorted(id_to_ten.keys())]
                    ids = sorted(id_to_ten.keys())
                    db.compile(expressions=expressions, ids=ids, elements=len(expressions))
                    print(f"[WARN] deserialize() failed; compiled DB for {class_name} from patterns")
            else:
                # No serialized DB; compile from patterns
                expressions = [id_to_ten[i].encode("ascii") for i in sorted(id_to_ten.keys())]
                ids = sorted(id_to_ten.keys())
                db.compile(expressions=expressions, ids=ids, elements=len(expressions))
                print(f"[OK] Compiled DB for {class_name} from patterns in {registry_dir}")
        except Exception as e:
            print(f"[ERROR] Hyperscan operations failed: {e}")
            db = None
    else:
        print("[SKIP] Hyperscan not installed; returning None DB (use pure-Python matcher).")

    return db, id_to_ten, id_to_score
