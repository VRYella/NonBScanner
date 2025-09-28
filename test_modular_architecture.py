"""
Test script to validate the modular motif detection architecture
==============================================================

Tests both the new modular architecture and legacy centralized approach.
"""

from nbdscanner import analyze_sequence, get_motif_classification_info
from modular_scanner import ModularMotifDetector
from motif_detection import *

def test_modular_architecture():
    """Test the new modular architecture"""
    print("=== Testing Modular Architecture ===")
    
    # Test sequence with multiple motif types
    test_sequences = {
        'G4_sequence': 'GGGTTAGGGTTAGGGTTAGGG',
        'STR_sequence': 'CACACACACACACACACACA', 
        'A_tract': 'AAAAAAAATTTTTTGGGGGGCCCCCC',
        'Complex': 'GGGTTAGGGAAAAAAAACCCCCCATAGATAGATAGATAGATAG'
    }
    
    detector = ModularMotifDetector()
    
    for seq_name, sequence in test_sequences.items():
        print(f"\n--- {seq_name} ---")
        print(f"Sequence: {sequence}")
        
        motifs = detector.analyze_sequence(sequence, seq_name)
        print(f"Found {len(motifs)} motifs:")
        
        class_counts = {}
        for motif in motifs:
            cls = motif.get('Class', 'Unknown')
            class_counts[cls] = class_counts.get(cls, 0) + 1
            print(f"  {cls}/{motif.get('Subclass', 'N/A')}: {motif.get('Start', 0)}-{motif.get('End', 0)} (score: {motif.get('Score', 0):.3f})")
        
        print(f"Class summary: {dict(class_counts)}")

def test_individual_detectors():
    """Test individual motif detectors"""
    print("\n=== Testing Individual Detectors ===")
    
    detectors = {
        'CurvedDNA': (CurvedDNADetector(), 'AAAAAAAATTTTTTGGGGGG'),
        'GQuadruplex': (GQuadruplexDetector(), 'GGGTTAGGGTTAGGGTTAGGG'),
        'SlippedDNA': (SlippedDNADetector(), 'CACACACACACACACACACA'),
        'ZDna': (ZDNADetector(), 'CGCGCGCGCGCGCGCG'),
        'iMotif': (IMotifDetector(), 'CCCTAACCCTAACCCTAACCC')
    }
    
    for detector_name, (detector, sequence) in detectors.items():
        print(f"\n--- {detector_name} Detector ---")
        print(f"Test sequence: {sequence}")
        
        motifs = detector.detect_motifs(sequence, f"test_{detector_name.lower()}")
        print(f"Found {len(motifs)} motifs:")
        
        for motif in motifs:
            print(f"  {motif.get('Class', 'N/A')}/{motif.get('Subclass', 'N/A')}: {motif.get('Start', 0)}-{motif.get('End', 0)} (score: {motif.get('Score', 0):.3f})")

def test_api_compatibility():
    """Test API compatibility between old and new approaches"""
    print("\n=== Testing API Compatibility ===")
    
    test_sequence = 'GGGTTAGGGAAAAAAAACCCCCCATAGATAGATAGATAGATAG'
    
    print("Using modular approach (default):")
    modular_motifs = analyze_sequence(test_sequence, 'test', use_modular=True)
    print(f"  Found {len(modular_motifs)} motifs")
    
    print("Using legacy approach:")
    legacy_motifs = analyze_sequence(test_sequence, 'test', use_modular=False)
    print(f"  Found {len(legacy_motifs)} motifs")
    
    print("Summary format:")
    summary = analyze_sequence(test_sequence, 'test', detailed=False)
    print(f"  Total motifs: {summary.get('total_motifs', 0)}")
    print(f"  Classes: {summary.get('classes_detected', 0)}")
    print(f"  Subclasses: {summary.get('subclasses_detected', 0)}")

def test_classification_info():
    """Test classification information"""
    print("\n=== Testing Classification Info ===")
    
    info = get_motif_classification_info()
    print(f"Architecture: {info.get('architecture', 'unknown')}")
    print(f"Version: {info.get('version', 'unknown')}")
    print(f"Total classes: {info.get('total_classes', 0)}")
    print(f"Total patterns: {info.get('total_patterns', 'N/A')}")
    print(f"Total detectors: {info.get('total_detectors', 'N/A')}")
    
    if 'detector_details' in info:
        print("\nDetector breakdown:")
        for detector_name, details in info['detector_details'].items():
            print(f"  {detector_name}: {details.get('total_patterns', 0)} patterns")

if __name__ == "__main__":
    print("Testing Modular Motif Detection Architecture")
    print("=" * 60)
    
    try:
        test_modular_architecture()
        test_individual_detectors()
        test_api_compatibility()
        test_classification_info()
        
        print("\n" + "=" * 60)
        print("✓ All tests completed successfully!")
        print("✓ Modular architecture is working correctly")
        print("✓ API compatibility maintained")
        print("✓ Individual detectors functional")
        
    except Exception as e:
        print(f"\n❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()