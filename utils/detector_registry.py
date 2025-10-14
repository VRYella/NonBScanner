#!/usr/bin/env python3
"""
Detector Registry and Test Sequence Manager
============================================

This module provides utilities to automatically discover all detector classes
and maintain test sequences for comprehensive testing.

Features:
- Automatic discovery of all detector classes from motif_detection module
- Test sequence registry for each detector class
- Easy integration for testing and demonstration scripts
"""

import inspect
import importlib
from typing import Dict, List, Tuple, Type, Any
from motif_detection.base_detector import BaseMotifDetector


def get_all_detector_classes() -> Dict[str, Type[BaseMotifDetector]]:
    """
    Automatically discover all detector classes from the motif_detection module.
    
    Returns:
        Dictionary mapping detector class names to their class objects
    """
    # Import the motif_detection module
    motif_detection = importlib.import_module('motif_detection')
    
    # Get all classes that are subclasses of BaseMotifDetector (excluding the base itself)
    detector_classes = {}
    
    for name in dir(motif_detection):
        obj = getattr(motif_detection, name)
        
        # Check if it's a class, subclass of BaseMotifDetector, and not the base class itself
        if (inspect.isclass(obj) and 
            issubclass(obj, BaseMotifDetector) and 
            obj != BaseMotifDetector):
            detector_classes[name] = obj
    
    return detector_classes


def get_detector_display_name(detector_class: Type[BaseMotifDetector]) -> str:
    """
    Get a human-readable display name for a detector class.
    
    Args:
        detector_class: The detector class
        
    Returns:
        Display name (e.g., "A-philic DNA", "Z-DNA", etc.)
    """
    try:
        # Try to instantiate and get the motif class name
        instance = detector_class()
        return instance.get_motif_class_name()
    except Exception:
        # Fallback to class name formatting
        name = detector_class.__name__.replace('Detector', '')
        # Convert camelCase to readable format
        import re
        name = re.sub(r'([a-z])([A-Z])', r'\1 \2', name)
        return name


def get_test_sequences_registry() -> Dict[str, List[Tuple[str, str]]]:
    """
    Registry of test sequences for each detector class.
    
    Returns:
        Dictionary mapping detector class names to lists of (description, sequence) tuples
    """
    return {
        'APhilicDetector': [
            ("Problem statement sequence", "AGGGGGGGGGAGGGGGGGGC"),
            ("Multiple G-runs", "AGGGGGGGGGTGGGGGGGGC"),
            ("C-rich (complement)", "CCCCCCCCCCCCCCCCCCCC"),
            ("Mixed GC-rich", "GGGGGCCCCCGGGGGCCCCC"),
        ],
        'ZDNADetector': [
            ("CG alternating repeats", "CGCGCGCGCGCGCG"),
            ("Longer CG repeat", "CGCGCGCGCGCGCGCGCGCG"),
            ("Embedded CG repeat", "ATCGCGCGCGCGCGAT"),
        ],
        'IMotifDetector': [
            ("C-rich telomeric", "CCCTAACCCTAACCCTAACCCT"),
            ("C-tract region", "CCCCTTCCCCTTCCCC"),
            ("Multiple C-runs", "CCCAACCCAACCCAACCC"),
        ],
        'GQuadruplexDetector': [
            ("G4 telomeric", "GGGTTAGGGTTAGGGTTAGGG"),
            ("Canonical G4", "GGGGTTTTGGGGTTTTGGGG"),
            ("Long loop G4", "GGGGAAAAAAAGGGGAAAAAAAGGGG"),
        ],
        'TriplexDetector': [
            ("Polypurine-Polypyrimidine", "AGAGAGAGAGAGAGAGAGAGAGAGAGAG"),
            ("GA repeat", "GAGAGAGAGAGAGAGAGAGAGA"),
            ("Purine tract", "AAAAGGGGAAAAGGGGAAAA"),
        ],
        'CruciformDetector': [
            ("Inverted repeat", "AAAATTTTAAAATTTT"),
            ("Palindrome", "ATCGATCGATCGATCG"),
            ("Long inverted repeat", "AAAAATTTTGGGGGCCCCCCCCCCCCCCGGGGG"),
        ],
        'SlippedDNADetector': [
            ("CAG repeat (Huntington's)", "CAGCAGCAGCAGCAGCAG"),
            ("CGG repeat (Fragile X)", "CGGCGGCGGCGGCGGCGG"),
            ("GAA repeat", "GAAGAAGAAGAAGAAGAA"),
        ],
        'CurvedDNADetector': [
            ("A-tract with phasing", "AAAAAAAATAAAAAAA"),
            ("Multiple A-tracts", "AAAAATTTTTAAAAATTTTTAAAAA"),
            ("AT-rich region", "AAAAAATTTTTTAAAAAATTTTTT"),
        ],
        'RLoopDetector': [
            ("G-rich skew", "GGGGGAAAAAGGGGGAAAAAGGGGGAAAAA"),
            ("Long G-skew", "GGGGGGGGAAAAAAAAAAGGGGGGGGAAAAAAAAAA"),
            ("Multiple G clusters", "GGGGGTTGGGGGTTGGGGG"),
        ],
    }


def get_demo_sequences_registry() -> Dict[str, List[Tuple[str, str]]]:
    """
    Registry of demonstration sequences for each detector class.
    Uses shorter list focused on best examples.
    
    Returns:
        Dictionary mapping detector class names to lists of (description, sequence) tuples
    """
    return {
        'APhilicDetector': [
            ("Problem statement sequence (G-rich)", "AGGGGGGGGGAGGGGGGGGC"),
            ("Poly-C (complement of A-philic)", "CCCCCCCCCCCCCCCCCCCC"),
        ],
        'ZDNADetector': [
            ("CG alternating repeats (classic Z-DNA)", "CGCGCGCGCGCGCG"),
            ("Embedded CG repeat", "ATCGCGCGCGCGCGAT"),
        ],
        'IMotifDetector': [
            ("C-rich telomeric sequence", "CCCTAACCCTAACCCTAACCCT"),
            ("Multiple C-runs with A spacers", "CCCAACCCAACCCAACCC"),
        ],
        'GQuadruplexDetector': [
            ("Telomeric G4 (human telomere)", "GGGTTAGGGTTAGGGTTAGGG"),
            ("G4 with A spacers", "GGGAAAGGGAAAGGGAAAGGG"),
        ],
        'TriplexDetector': [
            ("Polypurine-polypyrimidine (AG repeat)", "AGAGAGAGAGAGAGAGAGAGAGAGAGAG"),
            ("GA repeat (triplex forming)", "GAGAGAGAGAGAGAGAGAGAGA"),
        ],
        'CruciformDetector': [
            ("Inverted repeat (AT-rich)", "AAAATTTTAAAATTTT"),
            ("Palindromic sequence", "ATCGATCGATCGATCG"),
        ],
        'SlippedDNADetector': [
            ("CAG repeat (Huntington's disease)", "CAGCAGCAGCAGCAGCAG"),
            ("CGG repeat (Fragile X syndrome)", "CGGCGGCGGCGGCGGCGG"),
            ("GAA repeat (Friedreich's ataxia)", "GAAGAAGAAGAAGAAGAA"),
        ],
        'CurvedDNADetector': [
            ("A-tract with phasing", "AAAAAAAATAAAAAAA"),
            ("Long A-tract (intrinsic curvature)", "AAAAAAAAAAAAAAAAAAAA"),
        ],
        'RLoopDetector': [
            ("G-rich with A spacers (R-loop prone)", "GGGGGAAAAAGGGGGAAAAAGGGGGAAAAA"),
            ("Long G-runs with skew", "GGGGGGGGAAAAAAAAAAGGGGGGGGAAAAAAAAAA"),
        ],
    }


def get_all_detectors_with_test_sequences() -> List[Tuple[str, Type[BaseMotifDetector], List[Tuple[str, str]]]]:
    """
    Get all detectors with their test sequences in a convenient format.
    
    Returns:
        List of tuples: (display_name, detector_class, test_sequences)
    """
    detector_classes = get_all_detector_classes()
    test_sequences = get_test_sequences_registry()
    
    result = []
    for class_name, detector_class in sorted(detector_classes.items()):
        display_name = get_detector_display_name(detector_class)
        sequences = test_sequences.get(class_name, [])
        result.append((display_name, detector_class, sequences))
    
    return result


def get_all_detectors_with_demo_sequences() -> List[Tuple[str, Type[BaseMotifDetector], List[Tuple[str, str]]]]:
    """
    Get all detectors with their demo sequences in a convenient format.
    
    Returns:
        List of tuples: (display_name, detector_class, demo_sequences)
    """
    detector_classes = get_all_detector_classes()
    demo_sequences = get_demo_sequences_registry()
    
    result = []
    for class_name, detector_class in sorted(detector_classes.items()):
        display_name = get_detector_display_name(detector_class)
        sequences = demo_sequences.get(class_name, [])
        result.append((display_name, detector_class, sequences))
    
    return result


if __name__ == "__main__":
    # Quick test of the registry
    print("Discovered Detector Classes:")
    print("=" * 80)
    
    detectors = get_all_detector_classes()
    for class_name, detector_class in sorted(detectors.items()):
        display_name = get_detector_display_name(detector_class)
        print(f"  {class_name:30s} → {display_name}")
    
    print(f"\nTotal detectors discovered: {len(detectors)}")
    
    # Show test sequence counts
    print("\n\nTest Sequences per Detector:")
    print("=" * 80)
    test_seqs = get_test_sequences_registry()
    for class_name, sequences in sorted(test_seqs.items()):
        print(f"  {class_name:30s} → {len(sequences)} test sequences")
