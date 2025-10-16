#!/usr/bin/env python3
"""
Generate Hyperscan DBs & registries for Z-DNA and A-philic 10-mers.

Outputs (per class):
  - <outdir>/<CLASS>.hsdb         (serialized Hyperscan DB, if supported)
  - <outdir>/<CLASS>_patterns.txt (one 10-mer per line)
  - <outdir>/<CLASS>_registry.pkl (pickle with id->tenmer->score mapping + metadata)

Usage:
  python tools/generate_class_hsdb.py --out registry_dir
"""
import os
import sys
import argparse
import pickle
import json
from datetime import datetime

# Adjust path to import from parent directory
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Attempt to import Hyperscan
try:
    import hyperscan
    _HYPERSCAN_AVAILABLE = True
except Exception:
    hyperscan = None
    _HYPERSCAN_AVAILABLE = False

# Import detector classes
try:
    from motif_detection.z_dna_detector import ZDNADetector
    from motif_detection.a_philic_detector import APhilicDetector
except Exception as e:
    print(f"[ERROR] Cannot import detector classes: {e}")
    print("Make sure you're running from the repository root or PYTHONPATH is set correctly.")
    sys.exit(1)


def build_registry_from_table(table: dict) -> list:
    """
    Convert a mapping {10mer: score} -> list of records [(id, tenmer, score), ...]
    IDs are sequential integers starting at 0 (suitable for Hyperscan).
    """
    items = []
    # canonical order: sorted by tenmer for determinism
    for idx, ten in enumerate(sorted(table.keys())):
        score = float(table[ten])
        items.append((idx, ten.upper(), score))
    return items


def save_pattern_text(outdir: str, cls: str, items: list):
    path = os.path.join(outdir, f"{cls}_patterns.txt")
    with open(path, "w") as fh:
        for _, ten, _ in items:
            fh.write(ten + "\n")
    print(f"[OK] saved patterns text -> {path}")


def save_registry(outdir: str, cls: str, items: list, extra_meta: dict):
    registry = {
        "class": cls,
        "generated_at": datetime.utcnow().isoformat() + "Z",
        "n_patterns": len(items),
        "patterns": [{"id": i, "tenmer": ten, "score": score} for (i, ten, score) in items],
        "meta": extra_meta,
    }
    pkl_path = os.path.join(outdir, f"{cls}_registry.pkl")
    with open(pkl_path, "wb") as fh:
        pickle.dump(registry, fh)
    json_path = os.path.join(outdir, f"{cls}_registry.json")
    with open(json_path, "w") as fh:
        json.dump(registry, fh, indent=2)
    print(f"[OK] saved registry (pickle/json) -> {pkl_path}, {json_path}")


def compile_hyperscan_db(outdir: str, cls: str, items: list):
    """
    Compile hyperscan DB for given items: list of (id, tenmer, score).
    Saves serialized DB to outdir/<cls>.hsdb if hyperscan supports serialize().
    """
    if not _HYPERSCAN_AVAILABLE:
        print(f"[SKIP] Hyperscan not available — skipping binary DB compile for {cls}.")
        return False

    # Build expressions and ids (as bytes)
    expressions = []
    ids = []
    for (i, ten, score) in items:
        expressions.append(ten.encode("ascii"))
        ids.append(i)

    db = hyperscan.Database()
    print(f"[INFO] compiling Hyperscan DB for {cls} ({len(expressions)} patterns)...")
    try:
        db.compile(expressions=expressions, ids=ids, elements=len(expressions))
    except Exception as e:
        print(f"[ERROR] Hyperscan compilation failed for {cls}: {e}")
        return False

    # Try to serialize binary DB if supported
    try:
        serial = db.serialize()
        hsdb_path = os.path.join(outdir, f"{cls}.hsdb")
        with open(hsdb_path, "wb") as fh:
            fh.write(serial)
        print(f"[OK] saved serialized Hyperscan DB -> {hsdb_path}")
        return True
    except Exception as e:
        print(f"[WARN] Hyperscan DB compiled but serialize() not supported or failed: {e}")
        return True  # compile succeeded even if we couldn't write binary


def generate_class_db(outdir: str, cls: str, table: dict, extra_meta: dict = None):
    os.makedirs(outdir, exist_ok=True)
    extra_meta = extra_meta or {}
    items = build_registry_from_table(table)
    save_pattern_text(outdir, cls, items)
    save_registry(outdir, cls, items, extra_meta)
    compile_hyperscan_db(outdir, cls, items)
    return items


def main():
    p = argparse.ArgumentParser(description="Generate Hyperscan DBs & registries for Z-DNA and A-philic 10-mers.")
    p.add_argument("--out", "-o", default="registry", help="output directory for registries and DBs")
    p.add_argument("--skip-aphilic", action="store_true", help="skip A-philic generation")
    p.add_argument("--skip-zdna", action="store_true", help="skip Z-DNA generation")
    args = p.parse_args()

    outdir = args.out
    os.makedirs(outdir, exist_ok=True)

    # Z-DNA
    if not args.skip_zdna:
        print("[STEP] Building Z-DNA registry from ZDNADetector.TENMER_SCORE")
        zd_table = getattr(ZDNADetector, "TENMER_SCORE", None)
        if not zd_table:
            print("[ERROR] ZDNADetector.TENMER_SCORE not found in your module. Check import paths.")
        else:
            generate_class_db(outdir, "ZDNA", zd_table, extra_meta={"source": "ZDNADetector.TENMER_SCORE"})

    # A-philic
    if not args.skip_aphilic:
        print("[STEP] Building A-philic registry from APhilicDetector.TENMER_LOG2")
        aph_table = getattr(APhilicDetector, "TENMER_LOG2", None)
        if not aph_table:
            print("[ERROR] APhilicDetector.TENMER_LOG2 not found in your module. Check import paths.")
        else:
            # table might be log2 odds rather than direct 'score' — store the value as-is for post-filtering/scoring
            generate_class_db(outdir, "APhilic", aph_table, extra_meta={"source": "APhilicDetector.TENMER_LOG2", "note": "values are log2 odds"})

    print("[DONE] registry generation complete.")


if __name__ == "__main__":
    main()
