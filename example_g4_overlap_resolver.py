#!/usr/bin/env python3
"""
Example Usage of G4 Overlap Resolution Agent
============================================

This script demonstrates various use cases for the G4 Overlap Resolution Agent,
including:
1. Basic single sequence analysis
2. Multiple sequence analysis
3. Different output formats
4. Overlapping motif resolution

Author: Dr. Venkata Rajesh Yella
"""

from g4_overlap_resolver import G4OverlapResolver


def example_1_basic_usage():
    """Example 1: Basic single sequence analysis."""
    print("\n" + "="*80)
    print("Example 1: Basic Single Sequence Analysis")
    print("="*80)
    
    resolver = G4OverlapResolver()
    
    # Human telomeric repeat sequence (Zahler et al. 1991)
    sequence = "GGGTTAGGGTTAGGGTTAGGG"
    print(f"\nSequence: {sequence}")
    print(f"Length: {len(sequence)} bp")
    
    # Analyze sequence
    annotations = resolver.resolve_and_annotate(sequence, "telomeric_repeat")
    
    print(f"\nFound {len(annotations)} G4 motif(s):")
    for i, ann in enumerate(annotations, 1):
        print(f"\n  Motif {i}:")
        print(f"    Class: {ann['class_name']}")
        print(f"    Position: {ann['start']}-{ann['end']}")
        print(f"    Length: {ann['length']} bp")
        print(f"    Score: {ann['score']:.4f}")
        print(f"    Sequence: {ann['matched_seq']}")


def example_2_multiple_sequences():
    """Example 2: Analyzing multiple sequences."""
    print("\n" + "="*80)
    print("Example 2: Multiple Sequence Analysis")
    print("="*80)
    
    resolver = G4OverlapResolver()
    
    # Multiple test sequences with different G4 types
    sequences = {
        "telomeric": "GGGTTAGGGTTAGGGTTAGGG",
        "canonical": "GGGGTTTTGGGGTTTTGGGGTTTTGGGG",
        "relaxed": "GGGTTTTTTTTTTTGGGTTGGGTTGGG",
        "bulged": "GGGGAAAAAAAAAAAAAAAAAAAGGGGTTGGGTTGGGG"
    }
    
    print("\nAnalyzing sequences:")
    all_results = []
    
    for name, seq in sequences.items():
        print(f"\n  {name}:")
        print(f"    Sequence: {seq}")
        annotations = resolver.resolve_and_annotate(seq, name)
        print(f"    Found {len(annotations)} motif(s)")
        
        for ann in annotations:
            print(f"      - {ann['class_name']}: score={ann['score']:.4f}")
            all_results.append(ann)
    
    print(f"\nTotal motifs across all sequences: {len(all_results)}")


def example_3_output_formats():
    """Example 3: Different output formats."""
    print("\n" + "="*80)
    print("Example 3: Different Output Formats")
    print("="*80)
    
    resolver = G4OverlapResolver()
    sequence = "GGGGTTTTGGGGTTTTGGGGTTTTGGGG"
    annotations = resolver.resolve_and_annotate(sequence, "example_seq")
    
    # JSON format
    print("\nJSON Format:")
    print("-" * 40)
    json_output = resolver.format_as_json(annotations, pretty=True)
    print(json_output[:300] + "..." if len(json_output) > 300 else json_output)
    
    # Table format
    print("\nTable Format:")
    print("-" * 40)
    table_output = resolver.format_as_table(annotations)
    print(table_output)
    
    # BED format
    print("\nBED Format:")
    print("-" * 40)
    bed_output = resolver.format_as_bed(annotations, "chr1")
    print(bed_output)


def example_4_overlapping_motifs():
    """Example 4: Demonstrating overlap resolution."""
    print("\n" + "="*80)
    print("Example 4: Overlapping Motif Resolution")
    print("="*80)
    
    resolver = G4OverlapResolver()
    
    # Sequence with potential overlapping G4 motifs
    # Multiple patterns can match the same region
    sequence = "GGGTTGGGTTGGGTTGGGAAAAAGGGGTTTTGGGGTTTTGGGGTTTTGGGG"
    print(f"\nSequence: {sequence}")
    print(f"Length: {len(sequence)} bp")
    print("\nThis sequence contains multiple G-tracts that can form")
    print("overlapping G4 motifs of different classes.")
    
    annotations = resolver.resolve_and_annotate(sequence, "overlapping_test")
    
    print(f"\nAfter overlap resolution, found {len(annotations)} non-overlapping motif(s):")
    for i, ann in enumerate(annotations, 1):
        print(f"\n  Motif {i}:")
        print(f"    Class: {ann['class_name']}")
        print(f"    Position: {ann['start']}-{ann['end']}")
        print(f"    Score: {ann['score']:.4f}")
        print(f"    Sequence: {ann['matched_seq']}")
    
    # Verify no overlaps
    print("\nVerifying no overlaps:")
    for i in range(len(annotations)):
        for j in range(i + 1, len(annotations)):
            ann_i = annotations[i]
            ann_j = annotations[j]
            gap = ann_j['start'] - ann_i['end']
            print(f"  Gap between motif {i+1} and {j+1}: {gap} bp")


def example_5_scoring_details():
    """Example 5: Understanding G4Hunter scoring."""
    print("\n" + "="*80)
    print("Example 5: G4Hunter Scoring Details")
    print("="*80)
    
    resolver = G4OverlapResolver()
    
    # Compare different G4 sequences
    sequences = {
        "high_density": "GGGGTTTTGGGGTTTTGGGGTTTTGGGG",  # Dense G-tracts
        "telomeric": "GGGTTAGGGTTAGGGTTAGGG",            # Standard telomeric
        "relaxed": "GGTTTTTTTTGGGTTGGGTTGGG",            # Longer loops
    }
    
    print("\nComparing G4Hunter scores for different sequences:")
    
    for name, seq in sequences.items():
        annotations = resolver.resolve_and_annotate(seq, name)
        
        if annotations:
            ann = annotations[0]  # First motif
            details = ann['details']
            
            print(f"\n  {name}:")
            print(f"    Sequence: {seq}")
            print(f"    Final Score: {ann['score']:.4f}")
            print(f"    Details:")
            print(f"      - G-tracts: {details['n_g_tracts']}")
            print(f"      - Total G length: {details['total_g_len']}")
            print(f"      - GC balance: {details['gc_balance']:.4f}")
            print(f"      - Normalized window score: {details['normalized_window']:.4f}")
            print(f"      - Tract bonus: {details['tract_bonus']:.4f}")
            print(f"      - GC penalty: {details['gc_penalty']:.4f}")


def example_6_priority_demonstration():
    """Example 6: Demonstrating class priority in overlap resolution."""
    print("\n" + "="*80)
    print("Example 6: Class Priority in Overlap Resolution")
    print("="*80)
    
    resolver = G4OverlapResolver()
    
    print("\nClass Priority Order (highest to lowest):")
    priorities = [
        "canonical_g4",
        "relaxed_g4",
        "multimeric_g4",
        "bulged_g4",
        "imperfect_g4",
        "bipartite_g4",
        "g_triplex"
    ]
    
    for i, cls in enumerate(priorities, 1):
        print(f"  {i}. {cls}")
    
    print("\nWhen overlapping motifs have similar scores,")
    print("the agent selects the higher-priority class.")
    
    # Sequence that can match multiple classes
    sequence = "GGGGTTTTGGGGTTTTGGGGTTTTGGGG"
    annotations = resolver.resolve_and_annotate(sequence, "priority_demo")
    
    print(f"\nSequence: {sequence}")
    print(f"Detected {len(annotations)} motif(s) after resolution:")
    
    for ann in annotations:
        priority_rank = priorities.index(ann['class_name']) + 1 if ann['class_name'] in priorities else "N/A"
        print(f"  - {ann['class_name']} (priority rank: {priority_rank})")
        print(f"    Score: {ann['score']:.4f}")


def main():
    """Run all examples."""
    print("\n" + "="*80)
    print("G4 OVERLAP RESOLUTION AGENT - EXAMPLE DEMONSTRATIONS")
    print("="*80)
    print("\nThis script demonstrates various capabilities of the")
    print("G4 Overlap Resolution Agent.\n")
    
    # Run examples
    example_1_basic_usage()
    example_2_multiple_sequences()
    example_3_output_formats()
    example_4_overlapping_motifs()
    example_5_scoring_details()
    example_6_priority_demonstration()
    
    print("\n" + "="*80)
    print("EXAMPLES COMPLETE")
    print("="*80)
    print("\nFor more information, see:")
    print("  - G4_OVERLAP_RESOLVER_README.md")
    print("  - test_g4_overlap_resolver.py")
    print("  - python g4_overlap_resolver.py --help")
    print()


if __name__ == '__main__':
    main()
