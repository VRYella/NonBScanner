#!/usr/bin/env python3
"""
Test script to verify comprehensive output fields are present in NonBScanner results.

This script tests that all required fields from the problem statement are included
in the output/results with "NA" for non-applicable fields.

Required Fields:
1. ID
2. Sequence Name (or Accession)
3. Source (e.g., genome, experiment, study)
4. Motif Class
5. Motif Subclass
6. Pattern/Annotation ID
7. Start Position
8. End Position
9. Length (bp)
10. Sequence
11. Detection Method
12. Motif Score
13. Repeat/Tract Type
14. Left Arm Sequence
15. Right Arm Sequence
16. Loop Sequence
17. Arm Length
18. Loop Length
19. Stem Length(s)
20. Unit/Repeat Length
21. Number of Copies/Repeats
22. Spacer Length
23. Spacer Sequence
24. GC Content (%)
25. Structural Features (e.g., Tract Type, Curvature Score)
26. Strand
"""

from scanner import analyze_sequence, export_results_to_dataframe
from utilities import export_to_csv
import sys

def test_comprehensive_fields():
    """Test that all comprehensive fields are present in outputs"""
    
    # Required fields as per problem statement
    required_fields = [
        'ID',
        'Sequence_Name',  # Sequence Name (or Accession)
        'Source',  # Source (e.g., genome, experiment, study)
        'Class',  # Motif Class
        'Subclass',  # Motif Subclass
        'Pattern_ID',  # Pattern/Annotation ID
        'Start',  # Start Position
        'End',  # End Position
        'Length',  # Length (bp)
        'Sequence',  # Sequence
        'Method',  # Detection Method
        'Score',  # Motif Score
        'Repeat_Type',  # Repeat/Tract Type
        'Left_Arm',  # Left Arm Sequence
        'Right_Arm',  # Right Arm Sequence
        'Loop_Seq',  # Loop Sequence
        'Arm_Length',  # Arm Length
        'Loop_Length',  # Loop Length
        'Stem_Length',  # Stem Length(s)
        'Unit_Length',  # Unit/Repeat Length
        'Number_Of_Copies',  # Number of Copies/Repeats
        'Spacer_Length',  # Spacer Length
        'Spacer_Sequence',  # Spacer Sequence
        'GC_Content',  # GC Content (%)
        'Structural_Features',  # Structural Features
        'Strand'  # Strand information
    ]
    
    # Test sequence with multiple motif types
    test_seq = 'GGGTTAGGGTTAGGGTTAGGGAAAAAAAAAATTTTAAAAAAAAAACGCGCGCGCGCGATATATATTTTTATATAT'
    
    print("Testing NonBScanner Comprehensive Output Fields")
    print("=" * 70)
    
    # Analyze sequence
    print(f"\nAnalyzing test sequence ({len(test_seq)} bp)...")
    motifs = analyze_sequence(test_seq, 'test_sequence')
    print(f"✓ Detected {len(motifs)} motifs")
    
    # Test DataFrame export
    print("\nTesting DataFrame export...")
    df = export_results_to_dataframe(motifs)
    print(f"✓ DataFrame created with shape: {df.shape}")
    
    # Check for required fields
    print("\nVerifying required fields:")
    missing_fields = []
    for field in required_fields:
        if field in df.columns:
            print(f"  ✓ {field}")
        else:
            print(f"  ✗ {field} - MISSING!")
            missing_fields.append(field)
    
    if missing_fields:
        print(f"\n✗ FAILED: {len(missing_fields)} fields missing: {missing_fields}")
        return False
    
    print(f"\n✓ All {len(required_fields)} required fields present!")
    
    # Check for NaN values (should all be "NA" strings)
    nan_count = df.isna().sum().sum()
    na_string_count = (df == 'NA').sum().sum()
    print(f"\nNA value handling:")
    print(f"  NaN values: {nan_count} (should be 0)")
    print(f"  'NA' string values: {na_string_count}")
    
    if nan_count > 0:
        print(f"  ✗ FAILED: Found {nan_count} NaN values (should use 'NA' strings)")
        return False
    else:
        print(f"  ✓ No NaN values found (all use 'NA' strings)")
    
    # Test CSV export
    print("\nTesting CSV export...")
    csv_output = export_to_csv(motifs)
    csv_lines = csv_output.split('\n')
    csv_header = csv_lines[0].split(',')
    
    # Check that all required fields are in CSV header
    csv_missing = []
    for field in required_fields:
        if field not in csv_header:
            csv_missing.append(field)
    
    if csv_missing:
        print(f"  ✗ FAILED: {len(csv_missing)} fields missing from CSV: {csv_missing}")
        return False
    else:
        print(f"  ✓ All required fields present in CSV export")
    
    # Sample output
    print("\nSample output (first motif):")
    print("-" * 70)
    if len(motifs) > 0:
        sample = df.iloc[0]
        for field in required_fields[:10]:  # Show first 10 fields
            value = sample[field]
            print(f"  {field:20} : {value}")
        print("  ...")
    
    print("\n" + "=" * 70)
    print("✓ ALL TESTS PASSED!")
    print("=" * 70)
    return True

if __name__ == "__main__":
    try:
        success = test_comprehensive_fields()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
