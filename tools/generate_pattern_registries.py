#!/usr/bin/env python3
"""
Generate pattern registries for regex-based motif detectors.

This script creates registry files (.txt, .pkl, .json) for regex patterns
used by curved DNA, G4, i-motif, and HUR AC-motif detectors.

Unlike the 10-mer registries (A-philic, Z-DNA), these registries contain
regex patterns that can be compiled with Hyperscan for high-performance matching.
"""

import os
import json
import pickle
import logging
from datetime import datetime
from typing import Dict, List, Tuple, Any

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Try to import hyperscan
try:
    import hyperscan
    HYPERSCAN_AVAILABLE = True
    logger.info("Hyperscan is available")
except ImportError:
    HYPERSCAN_AVAILABLE = False
    logger.warning("Hyperscan not available - will generate registries without .hsdb files")


def create_registry_dir(output_dir: str):
    """Create registry directory if it doesn't exist."""
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Registry directory: {output_dir}")


def save_registry_files(class_name: str, output_dir: str, registry_data: Dict[str, Any]):
    """
    Save registry data in multiple formats (.txt, .pkl, .json).
    
    Args:
        class_name: Name of the detector class (e.g., 'CurvedDNA', 'G4')
        output_dir: Directory to save files
        registry_data: Registry dictionary containing patterns and metadata
    """
    base_path = os.path.join(output_dir, class_name)
    
    # Save plain text list of patterns
    txt_path = f"{base_path}_patterns.txt"
    with open(txt_path, 'w') as f:
        for pattern in registry_data['patterns']:
            f.write(f"{pattern['pattern']}\n")
    logger.info(f"Saved {txt_path}")
    
    # Save pickle
    pkl_path = f"{base_path}_registry.pkl"
    with open(pkl_path, 'wb') as f:
        pickle.dump(registry_data, f)
    logger.info(f"Saved {pkl_path}")
    
    # Save JSON
    json_path = f"{base_path}_registry.json"
    with open(json_path, 'w') as f:
        json.dump(registry_data, f, indent=2)
    logger.info(f"Saved {json_path}")


def compile_hyperscan_db(class_name: str, output_dir: str, patterns: List[Dict[str, Any]]):
    """
    Compile and serialize Hyperscan database for patterns.
    
    Args:
        class_name: Name of the detector class
        output_dir: Directory to save .hsdb file
        patterns: List of pattern dictionaries
    """
    if not HYPERSCAN_AVAILABLE:
        logger.info(f"Skipping Hyperscan DB compilation for {class_name} (not available)")
        return
    
    try:
        # Prepare expressions for compilation
        expressions = []
        ids = []
        flags = []
        
        for p in patterns:
            expressions.append(p['pattern'].encode('utf-8'))
            ids.append(p['id'])
            # Use CASELESS and DOTALL flags for DNA pattern matching
            flags.append(hyperscan.HS_FLAG_CASELESS | hyperscan.HS_FLAG_DOTALL)
        
        # Compile database
        db = hyperscan.Database()
        db.compile(
            expressions=expressions,
            ids=ids,
            elements=len(expressions),
            flags=flags
        )
        
        # Serialize if supported
        try:
            serialized = db.serialize()
            hsdb_path = os.path.join(output_dir, f"{class_name}.hsdb")
            with open(hsdb_path, 'wb') as f:
                f.write(serialized)
            logger.info(f"Saved Hyperscan DB: {hsdb_path}")
        except AttributeError:
            logger.warning(f"Hyperscan serialization not supported - skipping .hsdb file for {class_name}")
            
    except Exception as e:
        logger.error(f"Failed to compile Hyperscan DB for {class_name}: {e}")


def generate_curved_dna_registry(output_dir: str):
    """Generate registry for Curved DNA patterns (local + global curvature)."""
    logger.info("Generating Curved DNA registry...")
    
    # Patterns from problem statement
    LOCAL_CURVED_PATTERNS = [
        r"A{7,}",      # covers 7,8,9,... (002)
        r"T{7,}",      # covers 7,8,9,... (003)
    ]
    
    GLOBAL_CURVED_PATTERNS = [
        # ----- Global curvature: A-phased repeats -----
        # 3-tract APRs (center spacing 9–11)
        r"(?:A{3}[ACGT]{6,8}A{3}[ACGT]{6,8}A{3})",                   # 008
        r"(?:A{4}[ACGT]{5,7}A{4}[ACGT]{5,7}A{4})",                   # 009
        r"(?:A{5}[ACGT]{4,6}A{5}[ACGT]{4,6}A{5})",                   # 010
        r"(?:A{6}[ACGT]{3,5}A{6}[ACGT]{3,5}A{6})",                   # 011
        r"(?:A{7}[ACGT]{2,4}A{7}[ACGT]{2,4}A{7})",                   # 012
        r"(?:A{8}[ACGT]{1,3}A{8}[ACGT]{1,3}A{8})",                   # 013
        r"(?:A{9}[ACGT]{0,2}A{9}[ACGT]{0,2}A{9})",                   # 014

        # 4-tract APRs
        r"(?:A{3}[ACGT]{6,8}A{3}[ACGT]{6,8}A{3}[ACGT]{6,8}A{3})",    # 015
        r"(?:A{4}[ACGT]{5,7}A{4}[ACGT]{5,7}A{4}[ACGT]{5,7}A{4})",    # 016
        r"(?:A{5}[ACGT]{4,6}A{5}[ACGT]{4,6}A{5}[ACGT]{4,6}A{5})",    # 017
        r"(?:A{6}[ACGT]{3,5}A{6}[ACGT]{3,5}A{6}[ACGT]{3,5}A{6})",    # 018
        r"(?:A{7}[ACGT]{2,4}A{7}[ACGT]{2,4}A{7}[ACGT]{2,4}A{7})",    # 019
        r"(?:A{8}[ACGT]{1,3}A{8}[ACGT]{1,3}A{8}[ACGT]{1,3}A{8})",    # 020
        r"(?:A{9}[ACGT]{0,2}A{9}[ACGT]{0,2}A{9}[ACGT]{0,2}A{9})",    # 021

        # 5-tract APRs
        r"(?:A{3}[ACGT]{6,8}A{3}[ACGT]{6,8}A{3}[ACGT]{6,8}A{3}[ACGT]{6,8}A{3})",  # 022
        r"(?:A{4}[ACGT]{5,7}A{4}[ACGT]{5,7}A{4}[ACGT]{5,7}A{4}[ACGT]{5,7}A{4})",  # 023
        r"(?:A{5}[ACGT]{4,6}A{5}[ACGT]{4,6}A{5}[ACGT]{4,6}A{5}[ACGT]{4,6}A{5})",  # 024
        r"(?:A{6}[ACGT]{3,5}A{6}[ACGT]{3,5}A{6}[ACGT]{3,5}A{6}[ACGT]{3,5}A{6})",  # 025
        r"(?:A{7}[ACGT]{2,4}A{7}[ACGT]{2,4}A{7}[ACGT]{2,4}A{7}[ACGT]{2,4}A{7})",  # 026
        r"(?:A{8}[ACGT]{1,3}A{8}[ACGT]{1,3}A{8}[ACGT]{1,3}A{8}[ACGT]{1,3}A{8})",  # 027
        r"(?:A{9}[ACGT]{0,2}A{9}[ACGT]{0,2}A{9}[ACGT]{0,2}A{9}[ACGT]{0,2}A{9})",  # 028

        # ----- Global curvature: T-phased repeats (mirror analogs) -----
        # 3-tract T-tracts
        r"(?:T{3}[ACGT]{6,8}T{3}[ACGT]{6,8}T{3})",                   # 029
        r"(?:T{4}[ACGT]{5,7}T{4}[ACGT]{5,7}T{4})",                   # 030
        r"(?:T{5}[ACGT]{4,6}T{5}[ACGT]{4,6}T{5})",                   # 031
        r"(?:T{6}[ACGT]{3,5}T{6}[ACGT]{3,5}T{6})",                   # 032
        r"(?:T{7}[ACGT]{2,4}T{7}[ACGT]{2,4}T{7})",                   # 033
        r"(?:T{8}[ACGT]{1,3}T{8}[ACGT]{1,3}T{8})",                   # 034
        r"(?:T{9}[ACGT]{0,2}T{9}[ACGT]{0,2}T{9})",                   # 035

        # 4-tract T-tracts
        r"(?:T{3}[ACGT]{6,8}T{3}[ACGT]{6,8}T{3}[ACGT]{6,8}T{3})",    # 036
        r"(?:T{4}[ACGT]{5,7}T{4}[ACGT]{5,7}T{4}[ACGT]{5,7}T{4})",    # 037
        r"(?:T{5}[ACGT]{4,6}T{5}[ACGT]{4,6}T{5}[ACGT]{4,6}T{5})",    # 038
        r"(?:T{6}[ACGT]{3,5}T{6}[ACGT]{3,5}T{6}[ACGT]{3,5}T{6})",    # 039
        r"(?:T{7}[ACGT]{2,4}T{7}[ACGT]{2,4}T{7}[ACGT]{2,4}T{7})",    # 040
        r"(?:T{8}[ACGT]{1,3}T{8}[ACGT]{1,3}T{8}[ACGT]{1,3}T{8})",    # 041
        r"(?:T{9}[ACGT]{0,2}T{9}[ACGT]{0,2}T{9}[ACGT]{0,2}T{9})",    # 042

        # 5-tract T-tracts
        r"(?:T{3}[ACGT]{6,8}T{3}[ACGT]{6,8}T{3}[ACGT]{6,8}T{3}[ACGT]{6,8}T{3})",  # 043
        r"(?:T{4}[ACGT]{5,7}T{4}[ACGT]{5,7}T{4}[ACGT]{5,7}T{4}[ACGT]{5,7}T{4})",  # 044
        r"(?:T{5}[ACGT]{4,6}T{5}[ACGT]{4,6}T{5}[ACGT]{4,6}T{5}[ACGT]{4,6}T{5})",  # 045
        r"(?:T{6}[ACGT]{3,5}T{6}[ACGT]{3,5}T{6}[ACGT]{3,5}T{6}[ACGT]{3,5}T{6})",  # 046
        r"(?:T{7}[ACGT]{2,4}T{7}[ACGT]{2,4}T{7}[ACGT]{2,4}T{7}[ACGT]{2,4}T{7})",  # 047
        r"(?:T{8}[ACGT]{1,3}T{8}[ACGT]{1,3}T{8}[ACGT]{1,3}T{8}[ACGT]{1,3}T{8})",  # 048
        r"(?:T{9}[ACGT]{0,2}T{9}[ACGT]{0,2}T{9}[ACGT]{0,2}T{9}[ACGT]{0,2}T{9})",  # 049
    ]
    
    # Combine all patterns with metadata
    patterns = []
    pattern_id = 0
    
    # Add local patterns
    for pattern in LOCAL_CURVED_PATTERNS:
        patterns.append({
            'id': pattern_id,
            'pattern': pattern,
            'subclass': 'Local Curvature',
            'score': 0.95
        })
        pattern_id += 1
    
    # Add global patterns
    for pattern in GLOBAL_CURVED_PATTERNS:
        patterns.append({
            'id': pattern_id,
            'pattern': pattern,
            'subclass': 'Global Curvature',
            'score': 0.90
        })
        pattern_id += 1
    
    # Create registry data
    registry_data = {
        'class': 'CurvedDNA',
        'generated_at': datetime.utcnow().isoformat() + 'Z',
        'n_patterns': len(patterns),
        'patterns': patterns,
        'meta': {
            'source': 'Problem statement LOCAL_CURVED_PATTERNS and GLOBAL_CURVED_PATTERNS',
            'detector_class': 'CurvedDNADetector',
            'pattern_type': 'regex'
        }
    }
    
    # Save files
    save_registry_files('CurvedDNA', output_dir, registry_data)
    compile_hyperscan_db('CurvedDNA', output_dir, patterns)
    
    logger.info(f"Generated Curved DNA registry with {len(patterns)} patterns")


def generate_g4_registry(output_dir: str):
    """Generate registry for G-Quadruplex patterns."""
    logger.info("Generating G-Quadruplex registry...")
    
    # Patterns from problem statement (G4_HS_PATTERNS)
    G4_PATTERNS = [
        (0, "canonical_g4",    r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}"),
        (1, "relaxed_g4",      r"G{2,}[ACGT]{1,12}G{2,}[ACGT]{1,12}G{2,}[ACGT]{1,12}G{2,}"),
        (2, "long_loop_g4",    r"G{3,}[ACGT]{8,15}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}"),
        (3, "bulged_g4",       r"(?:G{2,3}[ACGT]{0,3}G{1,3})[ACGT]{1,7}G{2,3}[ACGT]{1,7}G{2,3}[ACGT]{1,7}G{2,3}"),
        (4, "multimeric_g4",   r"(?:G{3,}[ACGT]{1,7}){4,}G{3,}"),
        (5, "imperfect_g4",    r"G{2,}[ACGT]{1,10}[AG]G{1,3}[ACGT]{1,10}G{2,}[ACGT]{1,10}G{2,}"),
        (6, "g_triplex",       r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}"),
    ]
    
    patterns = []
    for pattern_id, subclass, pattern in G4_PATTERNS:
        patterns.append({
            'id': pattern_id,
            'pattern': pattern,
            'subclass': subclass,
            'score': 0.90  # Default score
        })
    
    # Create registry data
    registry_data = {
        'class': 'G4',
        'generated_at': datetime.utcnow().isoformat() + 'Z',
        'n_patterns': len(patterns),
        'patterns': patterns,
        'meta': {
            'source': 'Problem statement G4_HS_PATTERNS',
            'detector_class': 'GQuadruplexDetector',
            'pattern_type': 'regex'
        }
    }
    
    # Save files
    save_registry_files('G4', output_dir, registry_data)
    compile_hyperscan_db('G4', output_dir, patterns)
    
    logger.info(f"Generated G-Quadruplex registry with {len(patterns)} patterns")


def generate_imotif_registry(output_dir: str):
    """Generate registry for i-Motif patterns."""
    logger.info("Generating i-Motif registry...")
    
    # Patterns from problem statement (IMOTIF_PATTERNS + HUR_AC_PATTERNS)
    IMOTIF_PATTERNS = [
        r"C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}"
    ]
    
    HUR_AC_PATTERNS = [
        r"A{3}[ACGT]{4}C{3}[ACGT]{4}C{3}[ACGT]{4}C{3}",
        r"C{3}[ACGT]{4}C{3}[ACGT]{4}C{3}[ACGT]{4}A{3}",
        r"A{3}[ACGT]{5}C{3}[ACGT]{5}C{3}[ACGT]{5}C{3}",
        r"C{3}[ACGT]{5}C{3}[ACGT]{5}C{3}[ACGT]{5}A{3}",
        r"A{3}[ACGT]{6}C{3}[ACGT]{6}C{3}[ACGT]{6}C{3}",
        r"C{3}[ACGT]{6}C{3}[ACGT]{6}C{3}[ACGT]{6}A{3}",
    ]
    
    patterns = []
    pattern_id = 0
    
    # Add canonical i-motif pattern
    for pattern in IMOTIF_PATTERNS:
        patterns.append({
            'id': pattern_id,
            'pattern': pattern,
            'subclass': 'canonical_imotif',
            'score': 0.95
        })
        pattern_id += 1
    
    # Add HUR AC-motif patterns
    for pattern in HUR_AC_PATTERNS:
        patterns.append({
            'id': pattern_id,
            'pattern': pattern,
            'subclass': 'hur_ac_motif',
            'score': 0.85
        })
        pattern_id += 1
    
    # Create registry data
    registry_data = {
        'class': 'IMotif',
        'generated_at': datetime.utcnow().isoformat() + 'Z',
        'n_patterns': len(patterns),
        'patterns': patterns,
        'meta': {
            'source': 'Problem statement IMOTIF_PATTERNS and HUR_AC_PATTERNS',
            'detector_class': 'IMotifDetector',
            'pattern_type': 'regex'
        }
    }
    
    # Save files
    save_registry_files('IMotif', output_dir, registry_data)
    compile_hyperscan_db('IMotif', output_dir, patterns)
    
    logger.info(f"Generated i-Motif registry with {len(patterns)} patterns")


def main():
    """Main entry point for registry generation."""
    # Default output directory
    output_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'registry')
    
    # Create registry directory
    create_registry_dir(output_dir)
    
    # Generate all registries
    logger.info("="*60)
    logger.info("GENERATING PATTERN REGISTRIES")
    logger.info("="*60)
    
    generate_curved_dna_registry(output_dir)
    generate_g4_registry(output_dir)
    generate_imotif_registry(output_dir)
    
    logger.info("="*60)
    logger.info("REGISTRY GENERATION COMPLETE")
    logger.info("="*60)
    logger.info(f"All registries saved to: {output_dir}")


if __name__ == '__main__':
    main()
