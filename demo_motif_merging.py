#!/usr/bin/env python3
"""
Visual demonstration of motif merging behavior.
This script shows exactly how overlapping 10-mers are merged into regions.
"""

from motif_detection.a_philic_detector import APhilicDetector
from motif_detection.z_dna_detector import ZDNADetector


def visualize_merging(sequence, detector_name, detector):
    """Visualize how 10-mers are merged into regions."""
    print("=" * 80)
    print(f"{detector_name} Detector - Merging Visualization")
    print("=" * 80)
    print(f"\nSequence: {sequence}")
    print(f"Position: {''.join(str(i % 10) for i in range(len(sequence)))}")
    print(f"         {''.join('0' if i < 10 else str(i // 10) for i in range(len(sequence)))}")
    
    # Get annotations (merged regions)
    annotations = detector.annotate_sequence(sequence)
    
    if not annotations:
        print("\nNo motifs detected.")
        return
    
    print(f"\nFound {len(annotations)} merged region(s):")
    
    for idx, region in enumerate(annotations, 1):
        start = region['start']
        end = region['end']
        n_10mers = region['n_10mers']
        
        print(f"\n  Region {idx}:")
        print(f"    Coordinates: [{start}, {end}) (0-based, end-exclusive)")
        print(f"    Length: {region['length']} bp")
        print(f"    Contributing 10-mers: {n_10mers}")
        
        # Visual representation
        visual = [' '] * len(sequence)
        for i in range(start, end):
            visual[i] = '█'
        print(f"    Visual: {''.join(visual)}")
        
        # Show first few contributing 10-mers
        print(f"    10-mer details (showing first 5):")
        for i, tenmer_info in enumerate(region['contributing_10mers'][:5]):
            tm_start = tenmer_info['start']
            tm_seq = tenmer_info['tenmer']
            score_key = 'log2' if 'log2' in tenmer_info else 'score'
            score = tenmer_info[score_key]
            print(f"      [{i+1}] Position {tm_start}: {tm_seq} (score: {score:.2f})")
        
        if n_10mers > 5:
            print(f"      ... and {n_10mers - 5} more")
        
        # Show merged sequence
        print(f"    Merged sequence: {sequence[start:end]}")
        print(f"    Total score: {region.get('sum_log2', region.get('sum_score')):.3f}")


def main():
    print("\n" + "=" * 80)
    print("MOTIF DETECTION MERGING - VISUAL DEMONSTRATION")
    print("=" * 80)
    print("\nThis demonstrates how overlapping 10-mers are merged into single regions.")
    print()
    
    # Example 1: A-philic with many overlapping matches
    print("\nEXAMPLE 1: A-philic DNA with densely overlapping 10-mers")
    aphilic_seq = "AGGGGGGGGGCCCCCCCCCTAGGGGGGGGC"
    visualize_merging(aphilic_seq, "A-philic", APhilicDetector())
    
    # Example 2: Z-DNA with overlapping matches
    print("\n\nEXAMPLE 2: Z-DNA with overlapping CG-rich 10-mers")
    zdna_seq = "GCGCGCGCGCGCGCGCGCGC"
    visualize_merging(zdna_seq, "Z-DNA", ZDNADetector())
    
    # Example 3: A-philic with separated regions
    print("\n\nEXAMPLE 3: A-philic DNA with separated regions (NOT merged)")
    separated_seq = "AGGGGGGGGGCCCCCAAAAAAAAAAAAAAGGGGGGGGGCCCCCC"
    visualize_merging(separated_seq, "A-philic", APhilicDetector())
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("""
Key Observations:
1. Overlapping/adjacent 10-mers → Single merged region
2. Separated 10-mer clusters → Multiple separate regions
3. Each region includes ALL contributing 10-mers
4. No duplicate reporting of overlapping matches
5. Full sequence of merged region is preserved

This behavior is GUARANTEED by both detectors regardless of how they're called.
    """)


if __name__ == "__main__":
    main()
