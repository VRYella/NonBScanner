#!/usr/bin/env python3
"""
Complete Pipeline Test: Detection → Scoring → Overlap Resolution → Visualization
================================================================================

This test verifies the complete pipeline from sequence input to visualization output.
Tests the architecture described in the problem statement:
- Hyperscan database for feasible motifs
- Scoring algorithms
- Non-overlap resolution
- Ready for visualization
"""

import os
import sys

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from utils.modular_scanner import ModularMotifDetector
import pandas as pd


def test_complete_pipeline():
    """Test the complete detection pipeline"""
    print("="*80)
    print("COMPLETE PIPELINE TEST")
    print("="*80)
    
    # Create test sequence with multiple motif types
    test_sequence = (
        # G-quadruplex region
        "GGGTTAGGGTTAGGGTTAGGG" +
        # Spacer
        "ATCGATCG" +
        # Curved DNA (A-tract)
        "AAAAAAAAAAAA" +
        # Spacer  
        "GCGCGC" +
        # i-Motif
        "CCCTAACCCTAACCCTAACCC" +
        # Spacer
        "TATATATA" +
        # Z-DNA
        "CGCGCGCGCGCGCGCGCGCG" +
        # Spacer
        "GGGGGG" +
        # Slipped DNA (repeat)
        "CACACACACACACA" +
        # Spacer
        "TTTTTTTT" +
        # Potential cruciform-forming sequence (inverted repeat)
        "ATCGATCGAAAAAATCGATCGAT"
    )
    
    print(f"\nTest Sequence Length: {len(test_sequence)} bp")
    print(f"Test Sequence: {test_sequence[:80]}...")
    
    # Initialize detector
    print("\n[1] DETECTION PHASE")
    print("-" * 80)
    detector = ModularMotifDetector()
    
    # Check Hyperscan databases loaded
    if hasattr(detector, 'hsdb_map'):
        print(f"✓ Loaded {len(detector.hsdb_map)} Hyperscan databases")
        for cls_name, info in detector.hsdb_map.items():
            n_patterns = len(info.get('id_to_pattern', {}))
            print(f"  - {cls_name}: {n_patterns} patterns")
    
    # Run detection
    print("\n[2] SCANNING SEQUENCE")
    print("-" * 80)
    motifs = detector.analyze_sequence(test_sequence, "test_pipeline")
    
    print(f"✓ Detected {len(motifs)} total motifs")
    
    # Analyze by class
    by_class = {}
    for motif in motifs:
        cls = motif['Class']
        if cls not in by_class:
            by_class[cls] = []
        by_class[cls].append(motif)
    
    print(f"✓ Found motifs in {len(by_class)} classes:")
    for cls, class_motifs in sorted(by_class.items()):
        print(f"  - {cls}: {len(class_motifs)} motifs")
    
    # Test scoring
    print("\n[3] SCORING VERIFICATION")
    print("-" * 80)
    
    for motif in motifs:
        score = motif.get('Score', 0)
        seq = motif.get('Sequence', '')
        
        # Verify score is reasonable
        if score < 0:
            print(f"✗ INVALID: Negative score for {motif['Class']}/{motif['Subclass']}")
        elif 'Cluster' not in motif['Class']:
            # Regular motifs should have scores in reasonable range
            if score > 10000:  # Z-DNA can have high scores
                print(f"  ⚠️  High score: {motif['Class']}/{motif['Subclass']} = {score:.2f}")
            else:
                print(f"  ✓ {motif['Class']}/{motif['Subclass']}: score={score:.3f}")
    
    # Test overlap resolution
    print("\n[4] OVERLAP RESOLUTION CHECK")
    print("-" * 80)
    
    has_overlap = False
    checked_pairs = 0
    
    # Group by class/subclass and check for overlaps
    class_subclass_groups = {}
    for motif in motifs:
        key = f"{motif['Class']}-{motif['Subclass']}"
        if key not in class_subclass_groups:
            class_subclass_groups[key] = []
        class_subclass_groups[key].append(motif)
    
    for key, group_motifs in class_subclass_groups.items():
        for i, m1 in enumerate(group_motifs):
            for m2 in group_motifs[i+1:]:
                checked_pairs += 1
                # Check if they overlap
                start1, end1 = m1['Start'], m1['End']
                start2, end2 = m2['Start'], m2['End']
                
                if not (end1 <= start2 or end2 <= start1):
                    has_overlap = True
                    overlap_len = min(end1, end2) - max(start1, start2)
                    print(f"  ✗ OVERLAP in {key}:")
                    print(f"    Motif 1: {start1}-{end1} (len={end1-start1})")
                    print(f"    Motif 2: {start2}-{end2} (len={end2-start2})")
                    print(f"    Overlap: {overlap_len} bp")
    
    if not has_overlap:
        print(f"  ✓ No overlaps detected (checked {checked_pairs} pairs)")
    
    # Test hybrid and cluster detection
    print("\n[5] HYBRID/CLUSTER DETECTION")
    print("-" * 80)
    
    hybrid_motifs = [m for m in motifs if m['Class'] == 'Hybrid']
    cluster_motifs = [m for m in motifs if 'cluster' in m['Class'].lower()]
    regular_motifs = [m for m in motifs if m['Class'] != 'Hybrid' and 'cluster' not in m['Class'].lower()]
    
    print(f"  Regular motifs: {len(regular_motifs)}")
    print(f"  Hybrid motifs: {len(hybrid_motifs)}")
    print(f"  Cluster motifs: {len(cluster_motifs)}")
    
    if hybrid_motifs:
        print(f"\n  Hybrid details:")
        for h in hybrid_motifs[:3]:
            print(f"    {h['Subclass']}: {h['Start']}-{h['End']} (score: {h['Score']:.3f})")
    
    if cluster_motifs:
        print(f"\n  Cluster details:")
        for c in cluster_motifs[:3]:
            print(f"    {c['Subclass']}: {c['Start']}-{c['End']} (score: {c['Score']:.3f})")
    
    # Prepare for visualization (convert to DataFrame)
    print("\n[6] OUTPUT FORMAT FOR VISUALIZATION")
    print("-" * 80)
    
    # Convert to DataFrame (standard format for visualization)
    df = pd.DataFrame(motifs)
    
    print(f"  ✓ Created DataFrame with {len(df)} rows")
    print(f"  ✓ Columns: {list(df.columns)}")
    
    # Show sample rows
    print(f"\n  Sample output (first 5 motifs):")
    for i, row in df.head(5).iterrows():
        print(f"    {row['Class']:20s} {row['Subclass']:25s} "
              f"{row['Start']:4d}-{row['End']:4d} "
              f"score={row['Score']:.3f}")
    
    # Export capability check
    print("\n[7] EXPORT CAPABILITY CHECK")
    print("-" * 80)
    
    # CSV export
    csv_data = df.to_csv(index=False)
    print(f"  ✓ CSV export: {len(csv_data)} bytes")
    
    # BED format export
    bed_lines = []
    for _, row in df.iterrows():
        # BED format: chr start end name score strand
        # BED scores must be integers 0-1000
        bed_score = min(1000, max(0, int(row['Score'] * 10)))  # Scale and clamp to valid range
        bed_line = f"chr1\t{row['Start']}\t{row['End']}\t{row['Class']}_{row['Subclass']}\t{bed_score}\t+"
        bed_lines.append(bed_line)
    
    print(f"  ✓ BED export: {len(bed_lines)} features")
    print(f"    Sample: {bed_lines[0]}")
    
    # Summary statistics
    print("\n[8] SUMMARY STATISTICS")
    print("-" * 80)
    
    print(f"  Total sequence length: {len(test_sequence)} bp")
    print(f"  Total motifs detected: {len(motifs)}")
    print(f"  Motif classes: {len(by_class)}")
    print(f"  Average motif length: {df['Length'].mean():.1f} bp")
    print(f"  Average score: {df['Score'].mean():.3f}")
    print(f"  Coverage: {df['Length'].sum() / len(test_sequence) * 100:.1f}%")
    
    # Pipeline validation
    print("\n[9] PIPELINE VALIDATION")
    print("-" * 80)
    
    validations = {
        'Hyperscan databases loaded': hasattr(detector, 'hsdb_map') and len(detector.hsdb_map) > 0,
        'Motifs detected': len(motifs) > 0,
        'Multiple classes detected': len(by_class) > 1,
        'Scores are valid': all(m['Score'] >= 0 for m in motifs),
        'No overlaps within class': not has_overlap,
        'Output format ready': len(df) > 0,
        'Export capability': len(csv_data) > 0 and len(bed_lines) > 0
    }
    
    all_passed = True
    for check, passed in validations.items():
        status = "✓" if passed else "✗"
        print(f"  {status} {check}")
        if not passed:
            all_passed = False
    
    return all_passed


def main():
    """Run complete pipeline test"""
    print("\n" + "="*80)
    print("TESTING COMPLETE DETECTION PIPELINE")
    print("Detection → Scoring → Overlap Resolution → Output → Visualization")
    print("="*80)
    
    try:
        success = test_complete_pipeline()
        
        print("\n" + "="*80)
        if success:
            print("✓ COMPLETE PIPELINE TEST PASSED")
            print("="*80)
            print("\nPipeline is ready for production use:")
            print("  ✓ Hyperscan database acceleration active")
            print("  ✓ Scoring algorithms functional")
            print("  ✓ Overlap resolution working correctly")
            print("  ✓ Output format ready for visualization")
            print("="*80)
            return 0
        else:
            print("✗ COMPLETE PIPELINE TEST FAILED")
            print("="*80)
            return 1
            
    except Exception as e:
        print(f"\n✗ TEST FAILED WITH EXCEPTION: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
