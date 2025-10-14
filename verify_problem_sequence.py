#!/usr/bin/env python3
"""
Verification script specifically for the problem statement sequence.

This script demonstrates that the sequence AGGGGGGGGGAGGGGGGGGC
IS correctly detected as an A-philic DNA motif.
"""

print("=" * 80)
print("VERIFICATION: A-philic DNA Detection for Problem Statement Sequence")
print("=" * 80)
print()

# The sequence from the problem statement
test_sequence = "AGGGGGGGGGAGGGGGGGGC"
print(f"Test Sequence: {test_sequence}")
print(f"Length: {len(test_sequence)} bp")
print()

# ============================================================================
# Test 1: Direct detector test
# ============================================================================
print("-" * 80)
print("TEST 1: Direct A-philic Detector Test")
print("-" * 80)

from motif_detection.a_philic_detector import APhilicDetector

detector = APhilicDetector()

# Find 10-mer matches
matches = detector._find_10mer_matches(test_sequence)
print(f"\n✓ Found {len(matches)} 10-mer matches:")
for i, (start, tenmer, log2) in enumerate(matches, 1):
    print(f"  {i}. Position {start:2d}: {tenmer} (log2 score: {log2:.3f})")

# Check annotations (merged regions)
annotations = detector.annotate_sequence(test_sequence)
print(f"\n✓ Found {len(annotations)} merged region(s):")
for i, ann in enumerate(annotations, 1):
    print(f"  {i}. Region [{ann['start']}-{ann['end']}): "
          f"length={ann['length']} bp, "
          f"score={ann['sum_log2']:.3f}, "
          f"n_10mers={ann['n_10mers']}")

# Check detect_motifs output
motifs = detector.detect_motifs(test_sequence, "problem_sequence")
print(f"\n✓ detect_motifs() returned {len(motifs)} motif(s):")
for i, motif in enumerate(motifs, 1):
    print(f"  {i}. ID: {motif['ID']}")
    print(f"     Class: {motif['Class']}")
    print(f"     Subclass: {motif['Subclass']}")
    print(f"     Position: [{motif['Start']}-{motif['End']}]")
    print(f"     Sequence: {motif['Sequence']}")
    print(f"     Score: {motif['Score']:.3f}")

# ============================================================================
# Test 2: Full application test via modular_scanner
# ============================================================================
print()
print("-" * 80)
print("TEST 2: Full Application Test (via modular_scanner)")
print("-" * 80)

from utils.modular_scanner import analyze_sequence
from collections import defaultdict

all_motifs = analyze_sequence(test_sequence, "problem_sequence")
print(f"\n✓ Total motifs detected: {len(all_motifs)}")

# Group by class
by_class = defaultdict(list)
for motif in all_motifs:
    by_class[motif['Class']].append(motif)

print(f"\n✓ Motifs by class:")
for motif_class in sorted(by_class.keys()):
    class_motifs = by_class[motif_class]
    print(f"  - {motif_class}: {len(class_motifs)} motif(s)")
    
    # Show A-philic details
    if motif_class == "A-philic_DNA":
        for motif in class_motifs:
            print(f"    → {motif['Subclass']} [{motif['Start']}-{motif['End']}] "
                  f"Score: {motif['Score']:.3f}")

# ============================================================================
# Summary
# ============================================================================
print()
print("=" * 80)
print("VERIFICATION RESULT")
print("=" * 80)

if motifs:
    print("\n✓✓✓ SUCCESS ✓✓✓")
    print()
    print(f"The sequence '{test_sequence}' IS correctly detected")
    print(f"as an A-philic DNA motif with score {motifs[0]['Score']:.3f}")
    print()
    print("The A-philic DNA detector is working correctly!")
else:
    print("\n✗✗✗ FAILURE ✗✗✗")
    print()
    print(f"The sequence '{test_sequence}' was NOT detected.")
    print("There may be an issue with the detector.")

print("=" * 80)
