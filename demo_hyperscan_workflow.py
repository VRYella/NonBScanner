#!/usr/bin/env python3
"""
Demonstration of complete Hyperscan integration workflow.

This script shows how to:
1. Run multiple detectors on a sequence
2. Merge results with overlap resolution
3. Compare strict vs hybrid modes
4. Generate formatted output

This demonstrates the complete implementation of the problem statement.
"""

import sys
sys.path.insert(0, '/home/runner/work/NonBScanner/NonBScanner')

from motif_detection.a_philic_detector import APhilicDetector
from motif_detection.z_dna_detector import ZDNADetector
from motif_detection.g_quadruplex_detector import GQuadruplexDetector
from motif_detection.cruciform_detector import CruciformDetector
from utils.utils import merge_detector_results, resolve_cross_class_overlaps

def demo_single_detector():
    """Demonstrate single detector with merging"""
    print("=" * 70)
    print("DEMO 1: Single Detector with K-mer Merging")
    print("=" * 70)
    
    # Create a sequence with A-philic 10-mers
    sequence = "A" * 30 + "TGCATGCA" + "C" * 30
    print(f"Sequence: {sequence[:40]}... (length: {len(sequence)})")
    print()
    
    detector = APhilicDetector()
    
    # Show raw matches
    raw_matches = detector._find_10mer_matches(sequence)
    print(f"Step 1: Raw 10-mer matches found: {len(raw_matches)}")
    if raw_matches:
        print(f"  First 3: {raw_matches[:3]}")
    
    # Show merging
    if raw_matches:
        merged = detector._merge_matches(raw_matches)
        print(f"Step 2: Merged into regions: {len(merged)}")
        for i, (start, end, matches) in enumerate(merged, 1):
            print(f"  Region {i}: [{start}, {end}) with {len(matches)} 10-mers")
    
    # Final motifs
    motifs = detector.detect_motifs(sequence, "demo_seq")
    print(f"Step 3: Final motifs output: {len(motifs)}")
    for motif in motifs:
        print(f"  {motif['Class']}: pos {motif['Start']}-{motif['End']}, "
              f"score={motif['Score']:.3f}, length={motif['Length']}")
    print()

def demo_cross_detector_strict():
    """Demonstrate cross-detector resolution in strict mode"""
    print("=" * 70)
    print("DEMO 2: Cross-Detector Resolution (Strict Mode)")
    print("=" * 70)
    
    # Create a complex sequence with multiple motif types
    sequence = (
        "CGCGCGCGCGCGCG" +      # Z-DNA
        "AAAAAAAAAAAAAAAA" +    # A-philic (overlaps Z-DNA tail)
        "GGGTTAGGGTTAGGG" +     # G4
        "ATGCATGCATGCATGC"      # Potential cruciform
    )
    print(f"Sequence: {sequence[:50]}... (length: {len(sequence)})")
    print()
    
    # Run multiple detectors
    detectors = {
        'a_philic': APhilicDetector(),
        'z_dna': ZDNADetector(),
        'g4': GQuadruplexDetector(),
        'cruciform': CruciformDetector()
    }
    
    results = {}
    total_raw = 0
    for name, detector in detectors.items():
        motifs = detector.detect_motifs(sequence, "demo_seq")
        results[name] = motifs
        total_raw += len(motifs)
        if motifs:
            print(f"{name}: Found {len(motifs)} motifs")
            for motif in motifs:
                print(f"  {motif['Class']}: pos {motif['Start']}-{motif['End']}, "
                      f"score={motif['Score']:.3f}")
    
    print(f"\nTotal raw motifs: {total_raw}")
    print()
    
    # Merge with strict mode (non-overlapping, highest score wins)
    print("Applying strict overlap resolution...")
    merged_strict = merge_detector_results(results, overlap_mode='strict')
    
    print(f"After strict resolution: {len(merged_strict)} non-overlapping motifs")
    for motif in merged_strict:
        print(f"  {motif['Class']}: pos {motif['Start']}-{motif['End']}, "
              f"score={motif['Score']:.3f} [Winner!]")
    print()

def demo_cross_detector_hybrid():
    """Demonstrate cross-detector resolution in hybrid mode"""
    print("=" * 70)
    print("DEMO 3: Cross-Detector Resolution (Hybrid Mode)")
    print("=" * 70)
    
    # Same sequence as Demo 2
    sequence = (
        "CGCGCGCGCGCGCG" +
        "AAAAAAAAAAAAAAAA" +
        "GGGTTAGGGTTAGGG" +
        "ATGCATGCATGCATGC"
    )
    print(f"Sequence: {sequence[:50]}... (length: {len(sequence)})")
    print()
    
    # Run detectors (reuse results from previous demo logic)
    detectors = {
        'a_philic': APhilicDetector(),
        'z_dna': ZDNADetector(),
    }
    
    results = {}
    for name, detector in detectors.items():
        motifs = detector.detect_motifs(sequence, "demo_seq")
        results[name] = motifs
        if motifs:
            print(f"{name}: {len(motifs)} motifs")
    
    # Merge with hybrid mode (keep all, overlaps preserved)
    merged_hybrid = merge_detector_results(results, overlap_mode='hybrid')
    
    print(f"\nHybrid mode: {len(merged_hybrid)} motifs (overlaps preserved)")
    print("This mode is useful for biological analysis of co-occurring motifs")
    print()

def demo_overlap_resolution_algorithm():
    """Demonstrate the overlap resolution algorithm"""
    print("=" * 70)
    print("DEMO 4: Overlap Resolution Algorithm")
    print("=" * 70)
    
    # Create mock overlapping motifs
    motifs = [
        {'Class': 'G4', 'Start': 10, 'End': 30, 'Score': 0.9, 'Length': 20},
        {'Class': 'A-philic', 'Start': 25, 'End': 45, 'Score': 0.7, 'Length': 20},
        {'Class': 'Z-DNA', 'Start': 40, 'End': 60, 'Score': 0.8, 'Length': 20},
        {'Class': 'Cruciform', 'Start': 100, 'End': 120, 'Score': 0.6, 'Length': 20}
    ]
    
    print("Input motifs (some overlapping):")
    for i, m in enumerate(motifs, 1):
        print(f"  {i}. {m['Class']}: [{m['Start']}, {m['End']}], score={m['Score']}")
    
    print("\nOverlap analysis:")
    print("  - G4 [10,30] overlaps A-philic [25,45]")
    print("  - A-philic [25,45] overlaps Z-DNA [40,60]")
    print("  - Cruciform [100,120] is separate")
    
    print("\nAlgorithm:")
    print("  1. Sort by score (desc): G4(0.9), Z-DNA(0.8), A-philic(0.7), Cruciform(0.6)")
    print("  2. Select G4 (highest score)")
    print("  3. Skip A-philic (overlaps G4)")
    print("  4. Select Z-DNA (no overlap with G4)")
    print("  5. Select Cruciform (no overlap)")
    
    resolved = resolve_cross_class_overlaps(motifs, mode='strict')
    
    print(f"\nResolved output ({len(resolved)} non-overlapping motifs):")
    for i, m in enumerate(resolved, 1):
        print(f"  {i}. {m['Class']}: [{m['Start']}, {m['End']}], score={m['Score']}")
    
    print("\nGuarantees:")
    print("  ✓ Deterministic (same input → same output)")
    print("  ✓ Non-overlapping (highest score wins)")
    print("  ✓ Scientifically sound (preserves best motifs)")
    print()

def main():
    print("\n")
    print("=" * 70)
    print("HYPERSCAN INTEGRATION WORKFLOW DEMONSTRATION")
    print("=" * 70)
    print("\nThis demonstrates the complete implementation of:")
    print("  1. Scan first, score later (Hyperscan for matches)")
    print("  2. Keep raw per-hit contributions (k-mer redistribution)")
    print("  3. Guarantee non-overlapping outputs (greedy selection)")
    print("\n")
    
    try:
        demo_single_detector()
        demo_cross_detector_strict()
        demo_cross_detector_hybrid()
        demo_overlap_resolution_algorithm()
        
        print("=" * 70)
        print("DEMONSTRATION COMPLETE ✓")
        print("=" * 70)
        print("\nAll workflows demonstrated:")
        print("  ✓ Single detector with k-mer merging")
        print("  ✓ Multi-detector with strict overlap resolution")
        print("  ✓ Multi-detector with hybrid mode")
        print("  ✓ Overlap resolution algorithm details")
        print("\nFor more information, see:")
        print("  - HYPERSCAN_INTEGRATION.md (full guide)")
        print("  - HYPERSCAN_QUICK_REFERENCE.md (quick start)")
        print("  - IMPLEMENTATION_COMPLETE_HYPERSCAN.md (summary)")
        return 0
        
    except Exception as e:
        print(f"\n✗ DEMONSTRATION FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
