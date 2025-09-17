#!/usr/bin/env python3
"""
Example usage of scan_aphilic_hyperscan.py

This example demonstrates how to use the A-philic DNA scanner
with Intel Hyperscan acceleration.
"""

import scan_aphilic_hyperscan as sahp

def main():
    print("A-philic DNA Scanner with Intel Hyperscan")
    print("=" * 45)
    
    # Example sequences
    sequences = {
        "High A-philic content": (
            "NNNNN" +
            "AGGGGGGGGGCCCCTGGGGGCCCAAGGG" +
            "NNNNN"
        ),
        "Mixed sequence": (
            "ATCGATCGAGGGCCCCGGGGACCCCCCATGC"
        ),
        "Low A-philic content": (
            "ATCGATCGATCGATCGATCGATCG"
        )
    }
    
    print(f"Hyperscan available: {sahp.USE_HYPERSCAN}")
    print()
    
    for name, seq in sequences.items():
        print(f"Analyzing: {name}")
        print(f"Sequence: {seq}")
        print(f"Length: {len(seq)} bp")
        
        # Analyze the sequence
        output_file = f"{name.lower().replace(' ', '_')}_results.tsv"
        result_file = sahp.scan_sequence(
            seq, 
            sahp.TET_LOG2, 
            min_window=10, 
            max_window=15, 
            step=1, 
            out_tsv=output_file
        )
        
        # Count results
        with open(result_file, 'r') as f:
            lines = f.readlines()
        
        window_count = len(lines) - 1  # exclude header
        
        if window_count > 0:
            # Count by classification
            content = ''.join(lines)
            high_count = content.count("A_high_confidence")
            moderate_count = content.count("A_moderate")
            
            print(f"Results: {window_count} A-philic windows found")
            print(f"  - High confidence: {high_count}")
            print(f"  - Moderate confidence: {moderate_count}")
            print(f"  - Output saved to: {result_file}")
            
            # Show first few results
            if len(lines) > 1:
                print("  - First result:")
                headers = lines[0].strip().split('\t')
                first_result = lines[1].strip().split('\t')
                for h, v in zip(headers[:7], first_result[:7]):  # Show key fields
                    print(f"    {h}: {v}")
        else:
            print("Results: No A-philic windows found")
        
        print()

if __name__ == "__main__":
    main()