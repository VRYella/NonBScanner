#!/usr/bin/env python3
"""
Generate realistic test pathogenic genomes for NonBScanner analysis.
These sequences mimic the characteristics of real pathogenic genomes.
"""

import os
import random

# Set random seed for reproducibility
random.seed(42)

DATA_DIR = "data"
os.makedirs(DATA_DIR, exist_ok=True)

def generate_sequence(length, gc_content=0.5, pattern_freq=None):
    """Generate a DNA sequence with specified characteristics."""
    sequence = []
    
    # Basic nucleotide frequencies based on GC content
    at_prob = (1 - gc_content) / 2
    gc_prob = gc_content / 2
    
    for i in range(length):
        rand = random.random()
        if rand < at_prob:
            sequence.append('A')
        elif rand < 2 * at_prob:
            sequence.append('T')
        elif rand < 2 * at_prob + gc_prob:
            sequence.append('G')
        else:
            sequence.append('C')
    
    return ''.join(sequence)

def insert_motifs(sequence, motif_type="random"):
    """Insert specific Non-B DNA motifs into sequence."""
    seq_list = list(sequence)
    length = len(seq_list)
    
    # Insert G-quadruplex motifs
    for _ in range(max(1, length // 5000)):
        pos = random.randint(0, length - 30)
        g4_motif = "GGGGTTAGGGTTAGGGTTAGGG"  # Classic G4 motif
        for i, base in enumerate(g4_motif):
            if pos + i < length:
                seq_list[pos + i] = base
    
    # Insert Z-DNA prone sequences
    for _ in range(max(1, length // 3000)):
        pos = random.randint(0, length - 20)
        z_motif = "CGCGCGCGCGCG"  # Z-DNA forming sequence
        for i, base in enumerate(z_motif):
            if pos + i < length:
                seq_list[pos + i] = base
    
    # Insert A-tracts (curved DNA)
    for _ in range(max(1, length // 2000)):
        pos = random.randint(0, length - 15)
        a_tract = "AAAAAAAA"
        for i, base in enumerate(a_tract):
            if pos + i < length:
                seq_list[pos + i] = base
    
    # Insert inverted repeats (cruciform-forming)
    for _ in range(max(1, length // 4000)):
        pos = random.randint(0, length - 30)
        palindrome = "ATCGATCGATCGATCGAT"
        for i, base in enumerate(palindrome):
            if pos + i < length:
                seq_list[pos + i] = base
    
    return ''.join(seq_list)

# Pathogenic genome specifications
GENOMES = {
    "SARS-CoV-2": {
        "length": 29903,
        "gc_content": 0.38,
        "description": "Severe acute respiratory syndrome coronavirus 2",
        "type": "RNA virus (presented as DNA equivalent)"
    },
    "Hepatitis_B_Virus": {
        "length": 3215,
        "gc_content": 0.48,
        "description": "Hepatitis B virus genotype D",
        "type": "DNA virus"
    },
    "Human_Papillomavirus_16": {
        "length": 7906,
        "gc_content": 0.42,
        "description": "Human papillomavirus type 16",
        "type": "DNA virus"
    },
    "Ebola_Virus": {
        "length": 18959,
        "gc_content": 0.41,
        "description": "Zaire ebolavirus",
        "type": "RNA virus (presented as DNA equivalent)"
    },
    "Influenza_A_H1N1": {
        "length": 1778,
        "gc_content": 0.46,
        "description": "Influenza A virus H1N1 (HA segment)",
        "type": "RNA virus (presented as DNA equivalent)"
    }
}

def main():
    """Generate all test genomes."""
    print("=" * 80)
    print("GENERATING TEST PATHOGENIC GENOMES FOR NONBSCANNER ANALYSIS")
    print("=" * 80)
    print()
    
    for name, info in GENOMES.items():
        print(f"Generating {name}...")
        print(f"  Description: {info['description']}")
        print(f"  Length: {info['length']:,} bp")
        print(f"  GC Content: {info['gc_content']*100:.1f}%")
        print(f"  Type: {info['type']}")
        
        # Generate sequence
        sequence = generate_sequence(info['length'], info['gc_content'])
        
        # Insert Non-B DNA motifs
        sequence = insert_motifs(sequence)
        
        # Calculate actual GC content
        gc_count = sequence.count('G') + sequence.count('C')
        actual_gc = gc_count / len(sequence)
        
        # Write to FASTA file
        output_file = os.path.join(DATA_DIR, f"{name}.fasta")
        with open(output_file, 'w') as f:
            f.write(f">{name} | {info['description']}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + '\n')
        
        print(f"  Actual GC: {actual_gc*100:.1f}%")
        print(f"  ✓ Saved to {output_file}")
        print()
    
    print("=" * 80)
    print(f"Generation complete: {len(GENOMES)} genomes created")
    print("=" * 80)

if __name__ == "__main__":
    main()
