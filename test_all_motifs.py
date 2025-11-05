#!/usr/bin/env python3
"""
Test Script: Comprehensive Motif Detection Validation
=====================================================
This script creates an example DNA sequence that contains all motif classes
and subclasses, then validates that all are detected by the scanner.

Author: Auto-generated for NonBScanner validation
"""

from scanner import analyze_sequence, get_motif_classification_info
from collections import defaultdict
import json


def create_comprehensive_test_sequence():
    """
    Create a DNA sequence containing motifs from all 11 classes and their subclasses.
    Each motif pattern is carefully designed to trigger specific detectors.
    """
    
    # Build sequence by concatenating motif-specific patterns
    sequence_parts = []
    
    # Class 1: Curved DNA (Global curvature, Local Curvature)
    # A-tracts cause DNA bending
    sequence_parts.append("AAAAA")  # A-tract (Local curvature)
    sequence_parts.append("TTTTT")  # T-tract (Local curvature)
    sequence_parts.append("AAAATGCAAAATGCAAAATGCAAAA")  # Phased A-tracts (Global curvature)
    
    # Class 2: Slipped DNA (Direct Repeat, STR)
    sequence_parts.append("ATCGATCGATCGATCG")  # Direct repeat pattern
    sequence_parts.append("CACACACACACA")  # CA repeat (STR)
    sequence_parts.append("CGGCGGCGGCGG")  # CGG repeat (STR)
    
    # Class 3: Cruciform DNA (Inverted Repeats/Palindromes)
    # Inverted repeats form hairpin/cruciform structures
    sequence_parts.append("ATCGATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCGATAT")  # Inverted repeat with loop
    sequence_parts.append("GCTAGCNNNNNNNNNNNGCTAGC")  # Palindrome
    
    # Class 4: R-loop (R-loop formation sites)
    # GC-rich regions with AT spacers
    sequence_parts.append("GGGGGCCCCCATTTGGGGGCCCCC")  # GC-rich R-loop site
    sequence_parts.append("GGGGGATCCCCC")  # G-C rich region
    
    # Class 5: Triplex (Triplex, Sticky DNA)
    # Homopurine or homopyrimidine tracts
    sequence_parts.append("GGGGGAAAAAGGGGG")  # Homopurine tract (Triplex)
    sequence_parts.append("CCCCCTTTTTCCCCC")  # Homopyrimidine tract (Triplex)
    sequence_parts.append("GAGAGAGAGAGAGATCTCTCTCTCTCTC")  # Mirror repeat (Sticky DNA)
    
    # Class 6: G-Quadruplex Family (7 subclasses)
    # Multimeric G4
    sequence_parts.append("GGGAGGGAGGGAGGGAGGGAGGGA")  # Multimeric G4
    # Canonical G4
    sequence_parts.append("GGGATGGGTAGGGTGGGG")  # Canonical G4 (G3+N1-7G3+N1-7G3+N1-7G3+)
    # Relaxed G4
    sequence_parts.append("GGATCGGATCGGATCGG")  # Relaxed G4 (G2+N1-12G2+N1-12G2+N1-12G2+)
    # Bulged G4
    sequence_parts.append("GGGAAAAAAAAAAGGGATGGGTAGGG")  # Bulged G4 (with longer loop)
    # Bipartite G4
    sequence_parts.append("GGATCGATCGATCGATCGATCGATCGATCGGATAGGATGG")  # Bipartite G4 (long spacer)
    # Imperfect G4
    sequence_parts.append("GGATCAGGATCGGATCGG")  # Imperfect G4
    # G-Triplex
    sequence_parts.append("GGGATGGGAGGGT")  # G-Triplex intermediate
    
    # Class 7: i-Motif Family (3 subclasses)
    # Canonical i-motif
    sequence_parts.append("CCCATCCCGTCCCACCCC")  # Canonical i-motif (C3+N1-7C3+N1-7C3+N1-7C3+)
    # Relaxed i-motif
    sequence_parts.append("CCATGCCATGCCATGCC")  # Relaxed i-motif (C2+N1-12C2+N1-12C2+N1-12C2+)
    # AC-motif
    sequence_parts.append("ACACACACACAC")  # AC-motif
    sequence_parts.append("CACACACACACACA")  # CA variant
    
    # Class 8: Z-DNA (2 subclasses)
    # Z-DNA (alternating purine-pyrimidine)
    sequence_parts.append("CGCGCGCGCGCG")  # CG alternating (Z-DNA)
    sequence_parts.append("GCGCGCGCGCGC")  # GC alternating (Z-DNA)
    sequence_parts.append("ATATATATATATAT")  # AT alternating (Z-DNA)
    sequence_parts.append("TATATATATATATAT")  # TA alternating (Z-DNA)
    # eGZ DNA (CG-rich regions)
    sequence_parts.append("CCCGGGCCCGGG")  # CG-rich region (eGZ)
    sequence_parts.append("GGGGGGCCCCCC")  # G/C-rich region (eGZ)
    
    # Class 9: A-philic DNA
    sequence_parts.append("AAAAAAAAAA")  # Poly-A tract (A-philic DNA)
    sequence_parts.append("AAATAAAAATAAAAT")  # A-rich region (A-philic DNA)
    
    # Add spacers between major sections to avoid unwanted overlaps
    spacer = "NNNNN"
    full_sequence = spacer.join(sequence_parts)
    
    # Replace N with balanced bases to make a valid DNA sequence
    n_count = full_sequence.count('N')
    bases = ['A', 'C', 'G', 'T'] * (n_count // 4 + 1)
    n_index = 0
    result = []
    for char in full_sequence:
        if char == 'N':
            result.append(bases[n_index % len(bases)])
            n_index += 1
        else:
            result.append(char)
    full_sequence = ''.join(result)
    
    return full_sequence


def analyze_and_report(sequence, sequence_name="comprehensive_test"):
    """Analyze sequence and generate detailed report"""
    
    print("=" * 80)
    print(f"COMPREHENSIVE MOTIF DETECTION TEST")
    print("=" * 80)
    print(f"Sequence length: {len(sequence)} bp")
    print(f"Sequence name: {sequence_name}")
    print()
    
    # Get expected classification
    classification_info = get_motif_classification_info()
    print(f"Expected classes: {classification_info['total_classes']}")
    print(f"Expected subclasses: {classification_info['total_subclasses']}")
    print()
    
    # Analyze sequence
    print("Running motif detection...")
    motifs = analyze_sequence(sequence, sequence_name)
    print(f"Total motifs detected: {len(motifs)}")
    print()
    
    # Group by class and subclass
    class_motifs = defaultdict(list)
    subclass_motifs = defaultdict(list)
    
    for motif in motifs:
        motif_class = motif.get('Class', 'Unknown')
        motif_subclass = motif.get('Subclass', 'Unknown')
        class_motifs[motif_class].append(motif)
        subclass_key = f"{motif_class} :: {motif_subclass}"
        subclass_motifs[subclass_key].append(motif)
    
    # Report detected classes
    print("=" * 80)
    print("DETECTED CLASSES")
    print("=" * 80)
    detected_classes = sorted(class_motifs.keys())
    for i, class_name in enumerate(detected_classes, 1):
        count = len(class_motifs[class_name])
        print(f"{i:2}. {class_name:25} - {count:3} motifs")
    print()
    
    # Report detected subclasses
    print("=" * 80)
    print("DETECTED SUBCLASSES")
    print("=" * 80)
    detected_subclasses = sorted(subclass_motifs.keys())
    for i, subclass_key in enumerate(detected_subclasses, 1):
        count = len(subclass_motifs[subclass_key])
        print(f"{i:2}. {subclass_key:45} - {count:3} motifs")
    print()
    
    # Check coverage against expected classes
    print("=" * 80)
    print("COVERAGE ANALYSIS")
    print("=" * 80)
    
    expected_classes = {
        1: "Curved_DNA",
        2: "Slipped_DNA",
        3: "Cruciform",
        4: "R-Loop",
        5: "Triplex",
        6: "G-Quadruplex",
        7: "i-Motif",
        8: "Z-DNA",
        9: "A-philic_DNA",
        10: "Hybrid",
        11: "Non-B_DNA_Clusters"
    }
    
    missing_classes = []
    for class_num, class_name in expected_classes.items():
        if class_name not in detected_classes:
            missing_classes.append(f"{class_num}. {class_name}")
    
    if missing_classes:
        print("⚠ MISSING CLASSES:")
        for missing in missing_classes:
            print(f"  - {missing}")
    else:
        print("✓ ALL 11 CLASSES DETECTED!")
    print()
    
    # Check expected subclasses
    expected_subclasses = {
        "Curved_DNA": ["Global curvature", "Local Curvature"],
        "Slipped_DNA": ["Direct Repeat", "STR"],
        "Cruciform": ["Inverted Repeats"],
        "R-Loop": ["R-loop formation sites"],
        "Triplex": ["Triplex", "Sticky DNA"],
        "G-Quadruplex": ["Multimeric G4", "Canonical G4", "Relaxed G4", 
                         "Bulged G4", "Bipartite G4", "Imperfect G4", 
                         "G-Triplex intermediate"],
        "i-Motif": ["Canonical i-motif", "Relaxed i-motif", "AC-motif"],
        "Z-DNA": ["Z-DNA", "eGZ (Extruded-G) DNA"],
        "A-philic_DNA": ["A-philic DNA"]
    }
    
    missing_subclasses = []
    for class_name, expected_subs in expected_subclasses.items():
        for sub in expected_subs:
            subclass_key = f"{class_name} :: {sub}"
            found = any(subclass_key in detected_key for detected_key in detected_subclasses)
            if not found:
                missing_subclasses.append(subclass_key)
    
    if missing_subclasses:
        print("⚠ MISSING SUBCLASSES:")
        for missing in missing_subclasses:
            print(f"  - {missing}")
    else:
        print("✓ ALL PRIMARY SUBCLASSES DETECTED!")
    print()
    
    # Detailed motif listing
    print("=" * 80)
    print("DETAILED MOTIF LISTING")
    print("=" * 80)
    for i, motif in enumerate(motifs, 1):
        print(f"{i:3}. {motif.get('Class', 'N/A'):20} | {motif.get('Subclass', 'N/A'):30} | "
              f"Pos: {motif.get('Start', 0):4}-{motif.get('End', 0):4} | "
              f"Len: {motif.get('Length', 0):3} | Score: {motif.get('Score', 0):.3f}")
    
    print()
    print("=" * 80)
    
    # Return results for programmatic use
    return {
        'sequence': sequence,
        'motifs': motifs,
        'detected_classes': detected_classes,
        'detected_subclasses': detected_subclasses,
        'missing_classes': missing_classes,
        'missing_subclasses': missing_subclasses,
        'total_motifs': len(motifs),
        'class_count': len(detected_classes),
        'subclass_count': len(detected_subclasses)
    }


def save_example_sequence(sequence, filename="example_all_motifs.fasta"):
    """Save the example sequence to a FASTA file"""
    with open(filename, 'w') as f:
        f.write(">comprehensive_test_sequence | All motif classes and subclasses\n")
        # Write sequence in 80 character lines
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + "\n")
    print(f"Example sequence saved to: {filename}")


if __name__ == "__main__":
    # Create comprehensive test sequence
    print("Creating comprehensive test sequence...")
    sequence = create_comprehensive_test_sequence()
    
    # Analyze and report
    results = analyze_and_report(sequence)
    
    # Save example sequence
    save_example_sequence(sequence)
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Sequence length: {len(sequence)} bp")
    print(f"Total motifs detected: {results['total_motifs']}")
    print(f"Classes detected: {results['class_count']}/11")
    print(f"Subclasses detected: {results['subclass_count']}")
    
    if not results['missing_classes'] and not results['missing_subclasses']:
        print("\n✓✓✓ SUCCESS: All classes and primary subclasses detected! ✓✓✓")
    else:
        print(f"\n⚠ INCOMPLETE: {len(results['missing_classes'])} classes and "
              f"{len(results['missing_subclasses'])} subclasses missing")
        print("\nRecommendation: Adjust sequence to include missing motif patterns")
    
    print("=" * 80)
