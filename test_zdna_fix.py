#!/usr/bin/env python3
"""
Test script to verify Z-DNA and A-philic detector fixes for proper merging of overlapping 10-mers.

This test specifically addresses the issue where GCGCGCGCGCGCGCGCGCGC (20-mer) was being 
reported as only the first 10 nucleotides instead of the complete 20-mer as a single Z-DNA region.
"""

from detectors import ZDNADetector, APhilicDetector


def test_zdna_20mer_merging():
    """Test that a 20-mer Z-DNA sequence is reported as a single merged region."""
    print("=" * 80)
    print("TEST: Z-DNA 20-mer Merging")
    print("=" * 80)
    
    # The example sequence from the issue
    test_sequence = 'GCGCGCGCGCGCGCGCGCGC'  # 20-mer
    
    detector = ZDNADetector()
    
    # Test 1: Verify hyperscan finds all 11 matches at correct positions
    print("\n1. Testing hyperscan match positions...")
    matches_hs = detector._hs_find_matches(test_sequence.upper())
    assert len(matches_hs) == 11, f"Expected 11 matches, got {len(matches_hs)}"
    
    # Verify positions are 0-10 (not all at position 0)
    positions = [m[0] for m in matches_hs]
    expected_positions = list(range(11))
    assert positions == expected_positions, f"Expected positions {expected_positions}, got {positions}"
    print(f"   ✓ Found {len(matches_hs)} matches at positions {positions}")
    
    # Test 2: Verify annotations merge into single region
    print("\n2. Testing annotation merging...")
    annotations = detector.annotate_sequence(test_sequence)
    assert len(annotations) == 1, f"Expected 1 merged annotation, got {len(annotations)}"
    
    ann = annotations[0]
    assert ann['start'] == 0, f"Expected start=0, got {ann['start']}"
    assert ann['end'] == 20, f"Expected end=20, got {ann['end']}"
    assert ann['length'] == 20, f"Expected length=20, got {ann['length']}"
    assert ann['n_10mers'] == 11, f"Expected 11 contributing 10-mers, got {ann['n_10mers']}"
    print(f"   ✓ Single merged region: positions {ann['start']}-{ann['end']}, length {ann['length']}")
    print(f"   ✓ Contributing 10-mers: {ann['n_10mers']}")
    
    # Test 3: Verify motif detection reports the full 20-mer
    print("\n3. Testing motif detection...")
    motifs = detector.detect_motifs(test_sequence, 'test_seq')
    assert len(motifs) == 1, f"Expected 1 motif, got {len(motifs)}"
    
    motif = motifs[0]
    assert motif['Length'] == 20, f"Expected length=20, got {motif['Length']}"
    assert motif['Sequence'] == test_sequence, f"Expected sequence={test_sequence}, got {motif['Sequence']}"
    assert motif['Contributing_10mers'] == 11, f"Expected 11 contributing 10-mers, got {motif['Contributing_10mers']}"
    print(f"   ✓ Detected motif with correct length: {motif['Length']}")
    print(f"   ✓ Detected motif sequence: {motif['Sequence']}")
    print(f"   ✓ Contributing 10-mers: {motif['Contributing_10mers']}")
    
    print("\n" + "=" * 80)
    print("✓ Z-DNA 20-mer test PASSED!")
    print("=" * 80)


def test_aphilic_20mer_merging():
    """Test that a 20-mer A-philic sequence is reported as a single merged region."""
    print("\n" + "=" * 80)
    print("TEST: A-philic 20-mer Merging")
    print("=" * 80)
    
    # Test sequence with repeated 10-mer pattern
    test_sequence = 'GGGGGGGGGGGGGGGGGGGG'  # 20 Gs
    
    detector = APhilicDetector()
    
    # Test 1: Verify hyperscan finds all 11 matches at correct positions
    print("\n1. Testing hyperscan match positions...")
    matches_hs = detector._hs_find_matches(test_sequence.upper())
    assert len(matches_hs) == 11, f"Expected 11 matches, got {len(matches_hs)}"
    
    # Verify positions are 0-10
    positions = [m[0] for m in matches_hs]
    expected_positions = list(range(11))
    assert positions == expected_positions, f"Expected positions {expected_positions}, got {positions}"
    print(f"   ✓ Found {len(matches_hs)} matches at positions {positions}")
    
    # Test 2: Verify annotations merge into single region
    print("\n2. Testing annotation merging...")
    annotations = detector.annotate_sequence(test_sequence)
    assert len(annotations) == 1, f"Expected 1 merged annotation, got {len(annotations)}"
    
    ann = annotations[0]
    assert ann['start'] == 0, f"Expected start=0, got {ann['start']}"
    assert ann['end'] == 20, f"Expected end=20, got {ann['end']}"
    assert ann['length'] == 20, f"Expected length=20, got {ann['length']}"
    assert ann['n_10mers'] == 11, f"Expected 11 contributing 10-mers, got {ann['n_10mers']}"
    print(f"   ✓ Single merged region: positions {ann['start']}-{ann['end']}, length {ann['length']}")
    print(f"   ✓ Contributing 10-mers: {ann['n_10mers']}")
    
    print("\n" + "=" * 80)
    print("✓ A-philic 20-mer test PASSED!")
    print("=" * 80)


def test_zdna_30mer_merging():
    """Test that a 30-mer Z-DNA sequence is properly merged."""
    print("\n" + "=" * 80)
    print("TEST: Z-DNA 30-mer Merging")
    print("=" * 80)
    
    # 30-mer alternating CG sequence
    test_sequence = 'CGCGCGCGCGCGCGCGCGCGCGCGCGCGCG'  # 30-mer
    
    detector = ZDNADetector()
    
    print("\n1. Testing motif detection...")
    motifs = detector.detect_motifs(test_sequence, 'test_seq')
    assert len(motifs) == 1, f"Expected 1 motif, got {len(motifs)}"
    
    motif = motifs[0]
    assert motif['Length'] == 30, f"Expected length=30, got {motif['Length']}"
    assert motif['Sequence'] == test_sequence, f"Expected sequence={test_sequence}, got {motif['Sequence']}"
    print(f"   ✓ Detected motif with correct length: {motif['Length']}")
    print(f"   ✓ Contributing 10-mers: {motif['Contributing_10mers']}")
    
    print("\n" + "=" * 80)
    print("✓ Z-DNA 30-mer test PASSED!")
    print("=" * 80)


if __name__ == '__main__':
    test_zdna_20mer_merging()
    test_aphilic_20mer_merging()
    test_zdna_30mer_merging()
    
    print("\n" + "=" * 80)
    print("ALL TESTS PASSED! ✓✓✓")
    print("=" * 80)
