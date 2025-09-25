#!/usr/bin/env python3
"""
Example usage of NonBScanner consolidated modules
Shows how to use the individual scanner modules programmatically
"""

from aphilic_scanner import find_a_philic_dna, calculate_aphilic_score
from zdna_scanner import find_z_dna, calculate_zdna_score
from motif_utils import gc_content, reverse_complement, merge_overlapping_motifs
from visualization import create_text_visualization, generate_statistics_report

def example_basic_usage():
    """Basic example showing individual scanner usage"""
    
    print("=== NonBScanner Example Usage ===\n")
    
    # Test sequences
    aphilic_seq = "CCCCCCCCCCCCCCAAAAAAAAAATTTTTTTTT"
    zdna_seq = "ACGCGCGCGCGCGCGCGCGCGC"
    
    print(f"A-philic test sequence: {aphilic_seq}")
    print(f"Z-DNA test sequence: {zdna_seq}")
    print()
    
    # A-philic DNA detection
    print("--- A-philic DNA Detection ---")
    aphilic_motifs = find_a_philic_dna(aphilic_seq, "test_aphilic")
    print(f"Found {len(aphilic_motifs)} A-philic motifs:")
    
    for motif in aphilic_motifs[:5]:  # Show first 5
        score = motif.get('Raw_Score', 0)
        conf = motif.get('Confidence', 'Unknown')
        start = motif.get('Start', 0)
        end = motif.get('End', 0)
        print(f"  {start:3d}-{end:3d}: {score:.3f} ({conf})")
    
    if len(aphilic_motifs) > 5:
        print(f"  ... and {len(aphilic_motifs) - 5} more")
    print()
    
    # Z-DNA detection
    print("--- Z-DNA Detection ---")
    zdna_motifs = find_z_dna(zdna_seq, "test_zdna")
    print(f"Found {len(zdna_motifs)} Z-DNA motifs:")
    
    for motif in zdna_motifs[:5]:  # Show first 5
        score = motif.get('Raw_Score', 0)
        cg_content = motif.get('CG_Content', 0)
        conf = motif.get('Confidence', 'Unknown')
        start = motif.get('Start', 0)
        end = motif.get('End', 0)
        print(f"  {start:3d}-{end:3d}: {score:.1f} (CG: {cg_content:.1f}%, {conf})")
    
    if len(zdna_motifs) > 5:
        print(f"  ... and {len(zdna_motifs) - 5} more")
    print()


def example_scoring_functions():
    """Example showing direct use of scoring functions"""
    
    print("--- Direct Scoring Examples ---")
    
    # A-philic scoring
    test_seq = "CCCCCCCCCC"
    aphilic_score = calculate_aphilic_score(test_seq)
    print(f"A-philic score for '{test_seq}': {aphilic_score:.3f}")
    
    # Z-DNA scoring  
    test_seq = "CGCGCGCGCG"
    zdna_score = calculate_zdna_score(test_seq)
    cg_pct = gc_content(test_seq)
    print(f"Z-DNA score for '{test_seq}': {zdna_score:.1f}")
    print(f"GC content: {cg_pct:.1f}%")
    print()


def example_utilities():
    """Example showing utility functions"""
    
    print("--- Utility Functions ---")
    
    seq = "ATCGATCG"
    print(f"Original sequence: {seq}")
    print(f"Reverse complement: {reverse_complement(seq)}")
    print(f"GC content: {gc_content(seq):.1f}%")
    print()


def example_visualization():
    """Example showing visualization capabilities"""
    
    print("--- Visualization Example ---")
    
    # Create some sample motifs
    sample_motifs = [
        {
            'Class': 'A-philic DNA',
            'Start': 10,
            'End': 25,
            'Raw_Score': 2.5,
            'Confidence': 'High'
        },
        {
            'Class': 'Z-DNA', 
            'Start': 50,
            'End': 65,
            'Raw_Score': 58.0,
            'Confidence': 'Moderate'
        },
        {
            'Class': 'Z-DNA',
            'Start': 80, 
            'End': 95,
            'Raw_Score': 62.0,
            'Confidence': 'High'
        }
    ]
    
    # Text visualization
    viz = create_text_visualization(sample_motifs, 100, width=50)
    print(viz)
    print()
    
    # Statistics report
    stats = generate_statistics_report(sample_motifs)
    print(stats)


def example_combined_analysis():
    """Example showing combined analysis workflow"""
    
    print("\n=== Combined Analysis Workflow ===\n")
    
    # Sample sequence with both motif types
    sequence = "NNNCCCCCCCCCCCCACGCGCGCGCGCGCAAAAAAANNNN"
    print(f"Test sequence: {sequence}")
    print()
    
    # Detect both motif types
    aphilic_motifs = find_a_philic_dna(sequence, "combined_test")
    zdna_motifs = find_z_dna(sequence, "combined_test")
    
    # Combine results
    all_motifs = aphilic_motifs + zdna_motifs
    
    print(f"Total motifs found: {len(all_motifs)}")
    print(f"  A-philic: {len(aphilic_motifs)}")
    print(f"  Z-DNA: {len(zdna_motifs)}")
    print()
    
    # Merge overlapping regions
    merged_motifs = merge_overlapping_motifs(all_motifs, max_distance=5)
    print(f"After merging overlapping regions: {len(merged_motifs)}")
    print()
    
    # Show merged results
    for i, motif in enumerate(merged_motifs, 1):
        motif_class = motif.get('Class', 'Unknown')
        start = motif.get('Start', 0)
        end = motif.get('End', 0)
        score = motif.get('Raw_Score', 0)
        conf = motif.get('Confidence', 'Unknown')
        print(f"  {i}. {motif_class}: {start}-{end} (score: {score:.2f}, {conf})")


if __name__ == "__main__":
    example_basic_usage()
    example_scoring_functions() 
    example_utilities()
    example_visualization()
    example_combined_analysis()
    print("\n=== Example Complete ===")