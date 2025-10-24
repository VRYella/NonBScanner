#!/usr/bin/env python3
"""
Demonstration of Pattern Registry System

This script shows how to use the pattern registries (CurvedDNA, G4, IMotif)
for high-performance motif detection with Hyperscan.
"""

import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from utils.load_regex_registry import scan_with_registry, get_cached_registry


def demo_simple_scanning():
    """Demonstrate simple pattern scanning with registries."""
    print("="*70)
    print("SIMPLE PATTERN SCANNING DEMO")
    print("="*70)
    
    # Example: Scan for G-Quadruplex patterns
    print("\n1. Scanning for G-Quadruplex (G4) patterns")
    print("-" * 70)
    
    sequence = "GGGTTAGGGTTAGGGTTAGGG"  # Human telomeric repeat
    print(f"Sequence: {sequence}")
    
    matches = scan_with_registry('G4', sequence)
    print(f"\nFound {len(matches)} G4 pattern matches:")
    
    for start, end, pattern_id, subclass in matches:
        matched_seq = sequence[start:end]
        print(f"  • {subclass:20s} (ID {pattern_id}): pos {start:2d}-{end:2d}  '{matched_seq}'")
    
    # Example: Scan for Curved DNA patterns
    print("\n2. Scanning for Curved DNA patterns")
    print("-" * 70)
    
    sequence = "AAAAAAAAAAAAAATTTTTTTTTTT"  # Long A-tract and T-tract
    print(f"Sequence: {sequence}")
    
    matches = scan_with_registry('CurvedDNA', sequence)
    print(f"\nFound {len(matches)} Curved DNA pattern matches:")
    
    # Group by subclass
    by_subclass = {}
    for start, end, pattern_id, subclass in matches:
        if subclass not in by_subclass:
            by_subclass[subclass] = []
        by_subclass[subclass].append((start, end, pattern_id))
    
    for subclass, subclass_matches in by_subclass.items():
        print(f"\n  {subclass}:")
        for start, end, pattern_id in subclass_matches[:3]:  # Show first 3
            matched_seq = sequence[start:end]
            print(f"    • Pattern {pattern_id}: pos {start:2d}-{end:2d}  '{matched_seq}'")
    
    # Example: Scan for i-Motif patterns
    print("\n3. Scanning for i-Motif patterns")
    print("-" * 70)
    
    sequence = "CCCTAACCCTAACCCTAACCC"  # C-rich i-motif forming sequence
    print(f"Sequence: {sequence}")
    
    matches = scan_with_registry('IMotif', sequence)
    print(f"\nFound {len(matches)} i-Motif pattern matches:")
    
    for start, end, pattern_id, subclass in matches:
        matched_seq = sequence[start:end]
        print(f"  • {subclass:20s} (ID {pattern_id}): pos {start:2d}-{end:2d}  '{matched_seq}'")


def demo_registry_information():
    """Demonstrate accessing registry metadata."""
    print("\n" + "="*70)
    print("REGISTRY METADATA DEMO")
    print("="*70)
    
    registries = ['G4', 'CurvedDNA', 'IMotif']
    
    for reg_name in registries:
        print(f"\n{reg_name} Registry Information:")
        print("-" * 70)
        
        db, id_to_pattern, id_to_subclass, id_to_score = get_cached_registry(reg_name)
        
        print(f"Total patterns: {len(id_to_pattern)}")
        print(f"Hyperscan support: {'Yes' if db is not None else 'No (using regex fallback)'}")
        
        # Show pattern counts by subclass
        subclass_counts = {}
        for subclass in id_to_subclass.values():
            subclass_counts[subclass] = subclass_counts.get(subclass, 0) + 1
        
        print(f"\nPattern breakdown by subclass:")
        for subclass, count in sorted(subclass_counts.items()):
            print(f"  • {subclass:30s}: {count:2d} patterns")
        
        # Show a few example patterns
        print(f"\nExample patterns:")
        for i, (pattern_id, pattern) in enumerate(sorted(id_to_pattern.items())[:3]):
            subclass = id_to_subclass[pattern_id]
            score = id_to_score[pattern_id]
            print(f"  [{pattern_id}] {subclass:20s} (score: {score:.2f})")
            print(f"      Pattern: {pattern}")


def demo_performance_comparison():
    """Demonstrate performance benefits of Hyperscan."""
    print("\n" + "="*70)
    print("PERFORMANCE COMPARISON")
    print("="*70)
    
    import time
    import re
    
    # Create a longer test sequence
    test_seq = "GGGTTAGGGTTAGGGTTAGGG" * 100  # 2100 bp
    
    print(f"\nTest sequence length: {len(test_seq)} bp")
    print("\nComparing G4 pattern matching performance:")
    print("-" * 70)
    
    # Method 1: Hyperscan-based (registry)
    start = time.time()
    matches_hs = scan_with_registry('G4', test_seq)
    time_hs = time.time() - start
    
    print(f"Hyperscan (registry):  {time_hs*1000:.2f} ms - {len(matches_hs)} matches")
    
    # Method 2: Pure Python regex (for comparison)
    start = time.time()
    pattern = r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}"
    compiled = re.compile(pattern, re.IGNORECASE)
    matches_re = list(compiled.finditer(test_seq))
    time_re = time.time() - start
    
    print(f"Pure Python regex:     {time_re*1000:.2f} ms - {len(matches_re)} matches")
    
    if time_re > 0:
        speedup = time_re / time_hs if time_hs > 0 else float('inf')
        print(f"\nSpeedup: {speedup:.1f}x faster with Hyperscan")


def main():
    """Run all demonstrations."""
    print("\n" + "="*70)
    print("PATTERN REGISTRY SYSTEM - DEMONSTRATION")
    print("="*70)
    print("\nThis demo shows how to use the pattern registries for")
    print("high-performance motif detection in DNA sequences.")
    
    try:
        demo_simple_scanning()
        demo_registry_information()
        demo_performance_comparison()
        
        print("\n" + "="*70)
        print("DEMO COMPLETE")
        print("="*70)
        print("\nThe pattern registry system provides:")
        print("  ✓ Fast Hyperscan-based pattern matching")
        print("  ✓ Easy access to 58+ curated Non-B DNA patterns")
        print("  ✓ Automatic fallback to Python regex if needed")
        print("  ✓ Organized by motif class and subclass")
        print("\nFor more information, see:")
        print("  • registry/README.md - Registry documentation")
        print("  • utils/load_regex_registry.py - Registry loader API")
        print("  • tools/generate_pattern_registries.py - Registry generator")
        print("="*70)
        
    except Exception as e:
        print(f"\n✗ Demo failed with error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
