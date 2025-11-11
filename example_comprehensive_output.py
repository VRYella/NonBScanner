#!/usr/bin/env python3
"""
Example script demonstrating comprehensive output fields in NonBScanner.

This script shows how to use NonBScanner to get comprehensive analysis results
with all required fields including structural details for different motif types.
"""

from scanner import analyze_sequence, export_results_to_dataframe
from utilities import export_to_csv
import pandas as pd

# Set pandas display options for better output
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 50)

def main():
    print("=" * 80)
    print("NonBScanner - Comprehensive Output Fields Example")
    print("=" * 80)
    
    # Example sequences with different motif types
    examples = {
        'G-Quadruplex': 'GGGTTAGGGTTAGGGTTAGGG',
        'STR (Short Tandem Repeat)': 'CACACACACACACACACA',
        'Direct Repeat with Spacer': 'AAAAAAAAAATTTTAAAAAAAAAA',
        'Cruciform (Inverted Repeat)': 'ATATATATTTTTATATAT',
        'Z-DNA': 'CGCGCGCGCGCGCG',
    }
    
    for motif_type, sequence in examples.items():
        print(f"\n{motif_type}")
        print("-" * 80)
        print(f"Sequence: {sequence}")
        
        # Analyze sequence
        motifs = analyze_sequence(sequence, motif_type.replace(' ', '_'))
        
        if not motifs:
            print("No motifs detected")
            continue
        
        # Export to DataFrame
        df = export_results_to_dataframe(motifs)
        
        # Show relevant fields based on motif type
        print(f"\nDetected {len(motifs)} motif(s):")
        
        if 'Cruciform' in motif_type:
            # Show cruciform-specific fields
            cols = ['ID', 'Class', 'Subclass', 'Start', 'End', 'Length',
                   'Left_Arm', 'Right_Arm', 'Loop_Seq', 'Arm_Length', 'Loop_Length', 
                   'GC_Content', 'Score']
            print(df[cols].to_string(index=False))
        
        elif 'Direct' in motif_type:
            # Show direct repeat-specific fields
            cols = ['ID', 'Class', 'Subclass', 'Start', 'End', 'Length',
                   'Unit_Length', 'Spacer_Length', 'Spacer_Sequence', 
                   'GC_Content', 'Score']
            direct_df = df[df['Subclass'].str.contains('Direct', na=False)]
            if len(direct_df) > 0:
                print(direct_df[cols].to_string(index=False))
        
        elif 'STR' in motif_type:
            # Show STR-specific fields
            cols = ['ID', 'Class', 'Subclass', 'Start', 'End', 'Length',
                   'Repeat_Type', 'Unit_Length', 'Number_Of_Copies', 
                   'GC_Content', 'Score']
            print(df[cols].head(3).to_string(index=False))
        
        else:
            # Show standard fields for other motifs
            cols = ['ID', 'Class', 'Subclass', 'Start', 'End', 'Length',
                   'Sequence', 'Score', 'Method']
            print(df[cols].head(3).to_string(index=False))
    
    # Demonstrate CSV export
    print("\n" + "=" * 80)
    print("CSV Export Example")
    print("=" * 80)
    
    # Use G-Quadruplex sequence for demo
    seq = 'GGGTTAGGGTTAGGGTTAGGG'
    motifs = analyze_sequence(seq, 'csv_demo')
    
    if motifs:
        csv_output = export_to_csv(motifs[:1])  # Just first motif
        lines = csv_output.split('\n')
        
        print("\nCSV Header (comprehensive fields):")
        header_fields = lines[0].split(',')
        for i, field in enumerate(header_fields[:26], 1):  # Show first 26 comprehensive fields
            print(f"  {i:2d}. {field}")
        
        print(f"\n... and {len(header_fields) - 26} additional detector-specific fields")
        
        print("\nSample CSV row:")
        if len(lines) > 1:
            values = lines[1].split(',')
            for field, value in zip(header_fields[:10], values[:10]):  # Show first 10
                display_value = value if value != 'NA' and value != '' else 'NA'
                print(f"  {field:20} : {display_value}")
            print("  ...")
    
    print("\n" + "=" * 80)
    print("All outputs include comprehensive fields with 'NA' for non-applicable data")
    print("=" * 80)

if __name__ == "__main__":
    main()
