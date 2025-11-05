#!/usr/bin/env python3
"""
Simple usage example for the comprehensive test sequence
========================================================

This script demonstrates how to use the example_all_motifs.fasta file
with NonBScanner to detect all motif classes.
"""

from scanner import analyze_sequence

def main():
    # Read the example FASTA file
    print("Reading example_all_motifs.fasta...")
    with open('example_all_motifs.fasta', 'r') as f:
        lines = f.readlines()
        header = lines[0].strip()
        sequence = ''.join(line.strip() for line in lines[1:])
    
    print(f"Header: {header}")
    print(f"Sequence length: {len(sequence)} bp\n")
    
    # Analyze the sequence
    print("Analyzing sequence...")
    motifs = analyze_sequence(sequence, "example")
    
    # Summary
    classes = {}
    for motif in motifs:
        cls = motif.get('Class', 'Unknown')
        if cls not in classes:
            classes[cls] = []
        classes[cls].append(motif)
    
    print(f"\n{'='*60}")
    print("ANALYSIS RESULTS")
    print('='*60)
    print(f"Total motifs detected: {len(motifs)}")
    print(f"Unique classes detected: {len(classes)}\n")
    
    # Show class breakdown
    print("Class Distribution:")
    print("-" * 60)
    for cls in sorted(classes.keys()):
        print(f"  {cls:25} {len(classes[cls]):3} motifs")
    
    # Show first few motifs from each class
    print(f"\n{'='*60}")
    print("SAMPLE MOTIFS (first 2 from each class)")
    print('='*60)
    for cls in sorted(classes.keys()):
        print(f"\n{cls}:")
        for motif in classes[cls][:2]:  # Show first 2
            print(f"  Pos: {motif.get('Start', 0):4}-{motif.get('End', 0):4} | "
                  f"Subclass: {motif.get('Subclass', 'N/A'):30} | "
                  f"Score: {motif.get('Score', 0):.3f}")
    
    print(f"\n{'='*60}")
    print(f"✓ Example demonstrates detection of all {len(classes)} motif classes!")
    print('='*60)

if __name__ == "__main__":
    main()
