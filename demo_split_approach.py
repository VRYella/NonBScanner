#!/usr/bin/env python3
"""
Demo script showing the split approach implementation:
1. Hyperscan-based detection (all classes except Slipped DNA, Cruciform DNA, Triplex DNA)  
2. Non-hyperscan pure Python detection (Slipped DNA, Cruciform DNA, Triplex DNA)

Usage:
    python demo_split_approach.py
"""

import sys
import os

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from nbdscanner import MotifDetector
from nonb_pure_python import scan_sequence, scan_sequence_tsv_output

def demo_split_approach():
    """Demonstrate the split approach with test sequences"""
    
    # Test sequences targeting different motif types
    test_sequences = {
        "STR_test": "CACACACACACACAACACACACACACA",  # STRs (CA repeats)
        "Direct_repeat": "ATCGATCGATCGAAATCGATCGATCG",  # Direct repeats with spacer
        "Cruciform": "AAAGCTTTAGCGATCGCTAAAGCTTT",   # Palindromic sequence
        "Mirror_triplex": "AGGAGGAGGAGGATTAGGAGGAGGAG", # Mirror repeat with spacer
        "G4_sequence": "GGGTTAGGGTTAGGGTTAGGGTACC",    # G-quadruplex
        "Complex": "GGGTTAGGGAAAAAAAACCCCCCATAGATAGATAGATAGATAG"  # Mixed motifs
    }
    
    print("=== NBDScanner Split Approach Demo ===\n")
    
    # Initialize detector
    detector = MotifDetector()
    
    print("1. Testing Pure Python Scanner (Slipped/Cruciform/Triplex only):\n")
    
    for name, seq in test_sequences.items():
        print(f"--- {name} ---")
        print(f"Sequence: {seq}")
        
        # Pure Python scanner results
        pure_results = scan_sequence(seq, name)
        
        if pure_results:
            print("Pure Python Scanner Results:")
            for motif in pure_results:
                print(f"  {motif['Class']}/{motif['Subclass']}: {motif['Start']}-{motif['End']} (score: {motif['Score']:.3f})")
        else:
            print("  No pure Python motifs found")
        
        print()
    
    print("\n2. Testing Full NBDScanner (Split Approach - Hyperscan + Pure Python):\n")
    
    for name, seq in test_sequences.items():
        print(f"--- {name} ---")
        
        # Full NBDScanner results
        full_results = detector.analyze_sequence(seq, name)
        
        if full_results:
            print("NBDScanner Results (Split Approach):")
            hyperscan_count = 0
            pure_python_count = 0
            
            for motif in full_results:
                method = motif.get('Method', 'Unknown')
                if 'Pure_Python' in method:
                    pure_python_count += 1
                    marker = "[PP]"
                else:
                    hyperscan_count += 1
                    marker = "[HS]"
                    
                print(f"  {marker} {motif['Class']}/{motif['Subclass']}: {motif['Start']}-{motif['End']} (score: {motif['Score']:.3f})")
            
            print(f"  Summary: {hyperscan_count} Hyperscan, {pure_python_count} Pure Python")
        else:
            print("  No motifs found")
        
        print()

def demo_tsv_output():
    """Demo the original TSV output format"""
    print("\n3. Testing Pure Python Scanner TSV Output:\n")
    
    test_seq = "CACACACACACACAGGGTTAGGGAAAAAAAACCCCCCATAGATAGATAGATAGATAG"
    print(f"Test sequence: {test_seq}\n")
    
    print("TSV Output Format:")
    scan_sequence_tsv_output(test_seq, "demo_sequence", verbose=True)

def print_approach_info():
    """Print information about the split approach"""
    print("\n=== Split Approach Information ===\n")
    
    print("HYPERSCAN-BASED DETECTION:")
    print("  • Curved DNA")
    print("  • R-loop")  
    print("  • G-Quadruplex Family")
    print("  • i-Motif Family")
    print("  • Z-DNA")
    print("  • A-philic DNA")
    print("  • Hybrid (overlaps)")
    print("  • Non-B DNA Clusters")
    
    print("\nPURE PYTHON DETECTION:")
    print("  • Slipped DNA (STRs, Direct Repeats)")
    print("  • Cruciform DNA (Inverted Repeats)")
    print("  • Triplex DNA (Mirror Repeats)")
    
    print("\nPure Python Features:")
    print("  • 2-bit sequence encoding")
    print("  • Double rolling hash (64-bit)")
    print("  • K-mer seeding for efficiency")
    print("  • Binary search for maximal extensions")
    print("  • Normalized scoring (0-1 scale)")
    print("  • O(1) substring equality checks")

if __name__ == "__main__":
    print_approach_info()
    demo_split_approach()
    demo_tsv_output()
    
    print("\n=== Demo Complete ===")
    print("\nThe split approach successfully:")
    print("1. Uses pure Python for Slipped/Cruciform/Triplex DNA detection")
    print("2. Uses Hyperscan for all other Non-B DNA motif classes")
    print("3. Integrates seamlessly with the existing NBDScanner API")
    print("4. Maintains backward compatibility")