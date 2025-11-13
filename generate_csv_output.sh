#!/bin/bash

###############################################################################
# NonBScanner CSV Generator
# 
# Description: Shell script to generate CSV output for all motif classes and 
#              subclasses from FASTA input files
#
# Author: Dr. Venkata Rajesh Yella
# Version: 1.0
# License: MIT
###############################################################################

# Default values
INPUT_FILE=""
OUTPUT_DIR="output"
OUTPUT_PREFIX="nonbscanner"
VERBOSE=false

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored messages
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Function to display usage information
usage() {
    cat << EOF
Usage: $0 -i INPUT_FILE [OPTIONS]

Generate CSV output for all Non-B DNA motif classes and subclasses.

Required Arguments:
    -i, --input FILE        Input FASTA file containing DNA sequences

Optional Arguments:
    -o, --output DIR        Output directory (default: output)
    -p, --prefix PREFIX     Output file prefix (default: nonbscanner)
    -v, --verbose           Verbose output
    -h, --help              Display this help message

Examples:
    # Basic usage
    $0 -i sequences.fasta

    # Specify output directory and prefix
    $0 -i sequences.fasta -o results -p my_analysis

    # Verbose mode
    $0 -i sequences.fasta -v

Output Files:
    PREFIX_all_motifs.csv              - All detected motifs
    PREFIX_summary.csv                 - Summary statistics
    PREFIX_by_class.csv                - Motifs grouped by class
    PREFIX_by_subclass.csv             - Motifs grouped by subclass
    PREFIX_CLASS_NAME.csv              - Individual CSV for each class

Motif Classes Detected:
    1.  Curved DNA (2 subclasses)
    2.  Slipped DNA (2 subclasses)
    3.  Cruciform (1 subclass)
    4.  R-Loop (1 subclass)
    5.  Triplex (2 subclasses)
    6.  G-Quadruplex (7 subclasses)
    7.  i-Motif (3 subclasses)
    8.  Z-DNA (2 subclasses)
    9.  A-philic DNA (1 subclass)
    10. Hybrid (dynamic)
    11. Non-B DNA Clusters (dynamic)

EOF
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--prefix)
            OUTPUT_PREFIX="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            print_error "Unknown option: $1"
            usage
            ;;
    esac
done

# Check if input file is provided
if [ -z "$INPUT_FILE" ]; then
    print_error "Input file is required"
    usage
fi

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    print_error "Input file not found: $INPUT_FILE"
    exit 1
fi

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    print_error "python3 is not installed or not in PATH"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Print banner
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║       NonBScanner CSV Generator                                ║"
echo "║       Non-B DNA Motif Detection and CSV Export                 ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

print_info "Input file: $INPUT_FILE"
print_info "Output directory: $OUTPUT_DIR"
print_info "Output prefix: $OUTPUT_PREFIX"
echo ""

# Create Python script for analysis
PYTHON_SCRIPT=$(cat << 'PYTHON_EOF'
import sys
import os
from scanner import analyze_sequence, analyze_multiple_sequences, export_results_to_dataframe
import pandas as pd

def read_fasta(filename):
    """Read sequences from a FASTA file"""
    sequences = {}
    current_name = None
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        if current_name:
            sequences[current_name] = ''.join(current_seq)
    
    return sequences

def main(input_file, output_dir, output_prefix, verbose):
    # Read sequences
    if verbose:
        print(f"Reading sequences from {input_file}...")
    sequences = read_fasta(input_file)
    
    if verbose:
        print(f"Found {len(sequences)} sequence(s)")
        for name, seq in sequences.items():
            print(f"  {name}: {len(seq)} bp")
    
    # Analyze sequences
    if verbose:
        print("\nAnalyzing sequences...")
    results = analyze_multiple_sequences(sequences, use_multiprocessing=False)
    
    # Combine all results
    all_motifs = []
    for seq_name, motifs in results.items():
        all_motifs.extend(motifs)
    
    if verbose:
        print(f"Total motifs detected: {len(all_motifs)}")
    
    if len(all_motifs) == 0:
        print("WARNING: No motifs detected")
        return
    
    # Convert to DataFrame
    df = export_results_to_dataframe(all_motifs)
    
    # Save all motifs
    all_motifs_file = os.path.join(output_dir, f"{output_prefix}_all_motifs.csv")
    df.to_csv(all_motifs_file, index=False)
    if verbose:
        print(f"Saved: {all_motifs_file}")
    
    # Generate summary statistics
    summary = {
        'Total_Motifs': len(all_motifs),
        'Unique_Classes': len(df['Class'].unique()),
        'Unique_Subclasses': len(df['Subclass'].unique()),
        'Total_Sequences': len(sequences),
        'Total_Length': sum(len(seq) for seq in sequences.values())
    }
    
    # Count by class
    class_counts = df['Class'].value_counts().to_dict()
    for cls, count in class_counts.items():
        summary[f'Class_{cls}'] = count
    
    summary_df = pd.DataFrame([summary])
    summary_file = os.path.join(output_dir, f"{output_prefix}_summary.csv")
    summary_df.to_csv(summary_file, index=False)
    if verbose:
        print(f"Saved: {summary_file}")
    
    # Group by class
    by_class = df.groupby('Class').size().reset_index(name='Count')
    by_class_file = os.path.join(output_dir, f"{output_prefix}_by_class.csv")
    by_class.to_csv(by_class_file, index=False)
    if verbose:
        print(f"Saved: {by_class_file}")
    
    # Group by subclass
    by_subclass = df.groupby(['Class', 'Subclass']).size().reset_index(name='Count')
    by_subclass_file = os.path.join(output_dir, f"{output_prefix}_by_subclass.csv")
    by_subclass.to_csv(by_subclass_file, index=False)
    if verbose:
        print(f"Saved: {by_subclass_file}")
    
    # Save individual CSV files for each class
    for cls in df['Class'].unique():
        if cls != 'NA':
            class_df = df[df['Class'] == cls]
            class_name = cls.replace('/', '_').replace(' ', '_').replace('-', '_')
            class_file = os.path.join(output_dir, f"{output_prefix}_{class_name}.csv")
            class_df.to_csv(class_file, index=False)
            if verbose:
                print(f"Saved: {class_file} ({len(class_df)} motifs)")
    
    # Print summary to stdout
    print("\n" + "="*60)
    print("ANALYSIS SUMMARY")
    print("="*60)
    print(f"Total sequences analyzed: {len(sequences)}")
    print(f"Total motifs detected: {len(all_motifs)}")
    print(f"Unique classes: {len(df['Class'].unique())}")
    print(f"Unique subclasses: {len(df['Subclass'].unique())}")
    print("\nMotifs by class:")
    for cls, count in sorted(class_counts.items()):
        print(f"  {cls:30} {count:5} motifs")
    print("="*60)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py INPUT_FILE OUTPUT_DIR OUTPUT_PREFIX [VERBOSE]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    output_prefix = sys.argv[3]
    verbose = len(sys.argv) > 4 and sys.argv[4].lower() == 'true'
    
    main(input_file, output_dir, output_prefix, verbose)
PYTHON_EOF
)

# Run Python analysis
TEMP_SCRIPT="/tmp/nonbscanner_analysis_$$.py"
echo "$PYTHON_SCRIPT" > "$TEMP_SCRIPT"

if [ "$VERBOSE" = true ]; then
    print_info "Running analysis..."
    python3 "$TEMP_SCRIPT" "$INPUT_FILE" "$OUTPUT_DIR" "$OUTPUT_PREFIX" "true"
    EXIT_CODE=$?
else
    python3 "$TEMP_SCRIPT" "$INPUT_FILE" "$OUTPUT_DIR" "$OUTPUT_PREFIX" "false"
    EXIT_CODE=$?
fi

# Clean up temporary script
rm -f "$TEMP_SCRIPT"

# Check exit code
if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    print_success "Analysis complete!"
    print_success "Output files saved to: $OUTPUT_DIR"
    echo ""
    print_info "Generated files:"
    ls -lh "$OUTPUT_DIR"/${OUTPUT_PREFIX}*.csv | awk '{print "  " $9 " (" $5 ")"}'
else
    echo ""
    print_error "Analysis failed with exit code $EXIT_CODE"
    exit $EXIT_CODE
fi

exit 0
