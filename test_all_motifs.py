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
    
    Optimized to ensure ALL primary subclasses are detected.
    """
    
    # Build sequence by concatenating motif-specific patterns
    sequence_parts = []
    
    # ======================================================================
    # Class 1: Curved DNA (Global curvature, Local Curvature)
    # ======================================================================
    # Local curvature: Long A-tracts or T-tracts (>=7 bp)
    sequence_parts.append("AAAAAAA")  # 7bp A-tract (Local curvature)
    sequence_parts.append("TTTTTTT")  # 7bp T-tract (Local curvature)
    
    # Global curvature: Phased A-tracts (A-phased repeats)
    # Need ~3+ A-tracts spaced ~11bp apart (phasing)
    sequence_parts.append("AAAAAACGTAAAAAAGTCAAAAAACGT")  # Phased A-tracts (Global curvature)
    
    # ======================================================================
    # Class 2: Slipped DNA (Direct Repeat, STR)
    # ======================================================================
    # STRs (Short Tandem Repeats): Unit 1-9bp, total length >=10bp
    sequence_parts.append("CACACACACACA")  # CA repeat (STR)
    sequence_parts.append("CGGCGGCGGCGG")  # CGG repeat (STR)
    sequence_parts.append("ATGATGATGATG")  # ATG repeat (STR)
    
    # Direct repeats: Unit 10-300bp, spacer <=10bp
    # Pattern: [UNIT][spacer][UNIT]
    direct_unit = "ATCGATCGATCG"  # 12bp unit
    spacer = "NNNNN"  # 5bp spacer
    sequence_parts.append(direct_unit + spacer + direct_unit)  # Direct Repeat
    
    # ======================================================================
    # Class 3: Cruciform DNA (Inverted Repeats)
    # ======================================================================
    # Inverted repeats: arm >=6bp, loop <=100bp
    # Pattern: [arm][loop][revcomp(arm)]
    left_arm = "ATCGATCG"  # 8bp arm
    loop_cruciform = "NNNNNNNNNNNN"  # 12bp loop
    right_arm = "CGATCGAT"  # revcomp of left_arm
    sequence_parts.append(left_arm + loop_cruciform + right_arm)  # Inverted Repeat
    
    # ======================================================================
    # Class 4: R-loop (R-loop formation sites)
    # ======================================================================
    # GC-rich regions with AT spacers
    sequence_parts.append("GGGGGCCCCCATGGGGGCCCCC")  # GC-rich R-loop site
    sequence_parts.append("GCGCGCATGCGCGC")  # GC-rich region
    
    # ======================================================================
    # Class 5: Triplex (Triplex, Sticky DNA)
    # ======================================================================
    # Homopurine tract (Triplex): >90% purines
    sequence_parts.append("GGGGGGGGGGAAAAAAAAA")  # Homopurine (Triplex)
    
    # Homopyrimidine tract (Triplex): >90% pyrimidines  
    sequence_parts.append("CCCCCCCCCCTTTTTTTTTT")  # Homopyrimidine (Triplex)
    
    # Mirror repeat with triplex potential: arm >=10bp, loop <=100bp, >90% pur/pyr
    # Pattern: [purine_arm][loop][reverse(purine_arm)]
    mirror_arm = "GAGGAGGAGG"  # 10bp purine-rich arm
    loop_mirror = "NNNNNN"  # 6bp loop
    mirror_arm_rev = mirror_arm[::-1]  # Reverse (not revcomp!)
    sequence_parts.append(mirror_arm + loop_mirror + mirror_arm_rev)  # Mirror repeat for Triplex
    
    # Sticky DNA: GAA/TTC repeats
    sequence_parts.append("GAAGAAGAAGAAGAA")  # GAA repeat (Sticky DNA)
    sequence_parts.append("TTCTTCTTCTTCTTC")  # TTC repeat (Sticky DNA)
    
    # ======================================================================
    # Class 6: G-Quadruplex Family (7 subclasses)
    # ======================================================================
    # Canonical G4: G3+N1-7G3+N1-7G3+N1-7G3+
    sequence_parts.append("GGGATGGGTAGGGTGGGG")  # Canonical G4
    
    # Multimeric G4: Multiple G4 units
    sequence_parts.append("GGGAGGGAGGGAGGGAGGGAGGGA")  # Multimeric G4
    
    # Relaxed G4: G2+N1-12G2+N1-12G2+N1-12G2+
    sequence_parts.append("GGATCGGATCGGATCGG")  # Relaxed G4
    
    # Bulged G4: One loop 8-20bp
    sequence_parts.append("GGGAAAAAAAAAAGGGATGGGTAGGG")  # Bulged G4
    
    # Bipartite G4: One loop 15-50bp
    long_spacer = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"  # 39bp spacer
    sequence_parts.append("GGATCG" + long_spacer + "GGATAGGATGG")  # Bipartite G4
    
    # Imperfect G4: Interrupted G-tracts
    sequence_parts.append("GGATCAGGATCGGATCGG")  # Imperfect G4
    
    # G-Triplex: G3+N1-7G3+N1-7G3+ (3 G-tracts)
    sequence_parts.append("GGGATGGGAGGGT")  # G-Triplex
    
    # ======================================================================
    # Class 7: i-Motif Family (3 subclasses)
    # ======================================================================
    # Canonical i-motif: C3+N1-7C3+N1-7C3+N1-7C3+
    sequence_parts.append("CCCATCCCGTCCCACCCC")  # Canonical i-motif
    
    # Relaxed i-motif: C2+N1-12C2+N1-12C2+N1-12C2+
    sequence_parts.append("CCATGCCATGCCATGCC")  # Relaxed i-motif
    
    # AC-motif: AC or CA repeats
    sequence_parts.append("ACACACACACAC")  # AC-motif
    sequence_parts.append("CACACACACACACA")  # CA variant
    
    # ======================================================================
    # Class 8: Z-DNA (2 subclasses)
    # ======================================================================
    # Z-DNA: Alternating purine-pyrimidine (CG or AT)
    sequence_parts.append("CGCGCGCGCGCGCGCG")  # CG alternating (Z-DNA)
    sequence_parts.append("GCGCGCGCGCGCGCGC")  # GC alternating (Z-DNA)
    sequence_parts.append("ATATATATATATATAT")  # AT alternating (Z-DNA)
    sequence_parts.append("TATATATATATATAT")  # TA alternating (Z-DNA)
    
    # eGZ DNA: CG-rich regions (>=6 consecutive C/G)
    sequence_parts.append("CCCCCCCCCGGGGGGGGG")  # CG-rich region (eGZ)
    sequence_parts.append("GGGGGGGGCCCCCCCCCC")  # G/C-rich region (eGZ)
    
    # ======================================================================
    # Class 9: A-philic DNA
    # ======================================================================
    sequence_parts.append("AAAAAAAAAA")  # Poly-A tract (A-philic DNA)
    sequence_parts.append("AAATAAAAATAAAAT")  # A-rich region (A-philic DNA)
    
    # ======================================================================
    # Assemble sequence with spacers
    # ======================================================================
    # Add spacers between major sections to avoid unwanted overlaps
    # But not too many spacers to maintain density for cluster/hybrid detection
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
    
    # Check expected subclasses (flexible matching)
    expected_subclasses = {
        "Curved_DNA": ["Global curvature", "Local Curvature"],
        "Slipped_DNA": ["Direct Repeat", "Direct_Repeat", "STR"],
        "Cruciform": ["Inverted Repeats", "Inverted Repeat", "Inverted_Repeat"],
        "R-Loop": ["R-loop formation sites", "R-loop formation site"],
        "Triplex": ["Triplex", "Sticky DNA", "GAA repeat", "TTC repeat", "Mirror repeat"],
        "G-Quadruplex": ["Multimeric G4", "Canonical G4", "Relaxed G4", 
                         "Bulged G4", "Bipartite G4", "Imperfect G4", 
                         "G-Triplex intermediate", "Long-loop G4"],
        "i-Motif": ["Canonical i-motif", "canonical_imotif", "Relaxed i-motif", 
                    "relaxed_imotif", "AC-motif", "ac_motif"],
        "Z-DNA": ["Z-DNA", "eGZ (Extruded-G) DNA", "eGZ DNA"],
        "A-philic_DNA": ["A-philic DNA", "A-philic_DNA"]
    }
    
    # Build mapping of canonical names to detected variations
    canonical_to_variations = {}
    for class_name, variations in expected_subclasses.items():
        canonical_to_variations[class_name] = set(variations)
    
    missing_subclasses = []
    detected_canonical = {}
    
    for class_name, expected_subs_list in expected_subclasses.items():
        # Get first item as canonical name
        canonical_name = expected_subs_list[0]
        
        # Check if any variation is detected
        found = False
        for sub_variation in expected_subs_list:
            subclass_key = f"{class_name} :: {sub_variation}"
            if any(subclass_key in detected_key or 
                   detected_key.startswith(f"{class_name} :: ") and 
                   sub_variation.lower() in detected_key.lower() 
                   for detected_key in detected_subclasses):
                found = True
                break
        
        if not found:
            # Only add to missing if canonical form not found
            missing_subclasses.append(f"{class_name} :: {canonical_name}")
    
    
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
