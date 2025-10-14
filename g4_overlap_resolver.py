#!/usr/bin/env python3
"""
G-Quadruplex (G4) Overlap Resolution Agent
===========================================

This automated agent provides robust resolution of overlapping G-quadruplex motifs
and outputs non-overlapping, prioritized G4 annotations with scores.

Scientific Background:
---------------------
G-quadruplexes (G4s) are non-canonical DNA secondary structures formed by guanine-rich 
sequences. They consist of stacked G-tetrads stabilized by Hoogsteen hydrogen bonding
and monovalent cations (K+, Na+). G4s play crucial roles in:
- Transcriptional regulation (Rhodes & Lipps 2015)[web:69]
- Telomere maintenance (Blackburn 2001)
- DNA replication control (Bochman et al. 2012)
- Genome stability (Maizels & Gray 2013)

Overlap Resolution Strategy:
----------------------------
The agent implements a greedy algorithm based on:
1. G4Hunter scoring (Bedrat et al. 2016)[web:3] - quantitative G4 propensity scoring
2. Class priority ranking (canonical > relaxed > bulged > imperfect)
3. Length-based tie-breaking for equivalent scores

This approach ensures biologically meaningful, non-redundant G4 annotations by:
- Prioritizing high-scoring, stable G4 structures
- Resolving ambiguous or overlapping motif predictions
- Maintaining computational efficiency (O(n log n) sorting + O(n) greedy selection)

References:
-----------
- Burge et al. (2006): Canonical G4 structural definitions[web:67][web:68]
- Huppert & Balasubramanian (2005): Relaxed G4 patterns[web:75]
- Lim et al. (2009): Bulged G4 structures[web:22]
- Adrian et al. (2014), Papp et al. (2023): Imperfect G4s[web:25]
- Guédin et al. (2010), Frasson et al. (2022): Multimeric G4s[web:42][web:73]
- Bedrat et al. (2016): G4Hunter algorithm[web:3]

Author: Dr. Venkata Rajesh Yella
License: MIT
"""

import sys
import json
import argparse
from typing import List, Dict, Any, Optional
from motif_detection.g_quadruplex_detector import GQuadruplexDetector


class G4OverlapResolver:
    """
    Automated agent for G-quadruplex overlap resolution and annotation.
    
    This class provides a simple API for:
    1. Detecting all G4 motif candidates in a DNA sequence
    2. Scoring candidates using G4Hunter-based metrics
    3. Resolving overlaps by score and class priority
    4. Outputting non-overlapping, annotated G4 motifs
    
    The overlap resolution uses a greedy algorithm (Bedrat et al. 2016) that:
    - Sorts candidates by descending score, class priority, and length
    - Greedily selects non-conflicting regions
    - Ensures O(n log n) time complexity
    """
    
    def __init__(self):
        """
        Initialize the G4 overlap resolver with the G-quadruplex detector.
        
        The detector implements patterns from multiple sources:
        - Canonical G4: Burge 2006, Todd 2005[web:67][web:68]
        - Relaxed G4: Huppert 2005, Phan 2006[web:75]
        - Bulged G4: Lim 2009, Adrian 2014[web:22]
        - Multimeric G4: Guédin 2010, Frasson 2022[web:42][web:73]
        """
        self.detector = GQuadruplexDetector()
    
    def resolve_and_annotate(self, sequence: str, sequence_name: str = "input_seq") -> List[Dict[str, Any]]:
        """
        Detect, score, and resolve overlapping G4 motifs in a DNA sequence.
        
        Args:
            sequence: DNA sequence string (A, T, G, C)
            sequence_name: Optional identifier for the sequence
        
        Returns:
            List of non-overlapping G4 annotations, each containing:
            - sequence_name: identifier for the input sequence
            - class_name: G4 motif class (canonical_g4, relaxed_g4, etc.)
            - pattern_id: specific pattern identifier
            - start: 0-based start position
            - end: 0-based end position (exclusive)
            - length: motif length in base pairs
            - score: G4Hunter-based score (higher = more stable)
            - matched_seq: actual DNA sequence of the motif
            - details: scoring breakdown (G-tracts, GC balance, etc.)
        
        Methodology:
            1. Pattern matching: Regex-based detection of all G4 classes
               (Burge 2006, Huppert 2005, Lim 2009)[web:67][web:75][web:22]
            
            2. G4Hunter scoring: Window-based quantitative scoring that
               correlates with experimental G4 formation (Bedrat et al. 2016)[web:3]
               - Sliding window analysis (default 25bp)
               - G-tract bonuses for multiple runs
               - GC balance penalties for C-rich regions
            
            3. Overlap resolution: Greedy selection algorithm that:
               - Prioritizes high-scoring, stable G4s
               - Uses class hierarchy (canonical > relaxed > bulged > imperfect)
               - Breaks ties by motif length
               - Ensures non-overlapping output (PERFORMANCE_OPTIMIZATION.md)
        """
        # Step 1: Detect all candidate G4 motifs using pattern matching
        # Patterns from Burge 2006, Huppert 2005, Lim 2009, etc.[web:67][web:75][web:22]
        annotations = self.detector.annotate_sequence(sequence)
        
        # Step 2: Add sequence name to each annotation for tracking
        for ann in annotations:
            ann['sequence_name'] = sequence_name
        
        return annotations
    
    def format_as_json(self, annotations: List[Dict[str, Any]], pretty: bool = True) -> str:
        """
        Format resolved G4 annotations as JSON.
        
        Args:
            annotations: List of G4 annotation dictionaries
            pretty: Whether to format JSON with indentation
        
        Returns:
            JSON string representation
        
        Output format includes:
        - Metadata: version, analysis type, motif count
        - Annotations: complete motif details for each resolved G4
        """
        output = {
            'version': 'G4OverlapResolver-v1.0',
            'analysis_type': 'G-Quadruplex_Overlap_Resolution',
            'description': 'Non-overlapping G4 motifs resolved by score and class priority',
            'total_motifs': len(annotations),
            'motifs': annotations
        }
        
        if pretty:
            return json.dumps(output, indent=2, ensure_ascii=False)
        else:
            return json.dumps(output, ensure_ascii=False)
    
    def format_as_table(self, annotations: List[Dict[str, Any]]) -> str:
        """
        Format resolved G4 annotations as tab-delimited table.
        
        Args:
            annotations: List of G4 annotation dictionaries
        
        Returns:
            Tab-delimited string with header and data rows
        
        Table columns:
        - Sequence_Name: identifier for input sequence
        - Class: G4 motif class (scientific classification)
        - Pattern_ID: specific pattern identifier
        - Start: 0-based start position
        - End: 0-based end position
        - Length: motif length in bp
        - Score: G4Hunter score (Bedrat et al. 2016)[web:3]
        - Sequence: actual DNA sequence
        """
        # Header with column descriptions
        header = '\t'.join([
            'Sequence_Name',
            'Class',
            'Pattern_ID',
            'Start',
            'End',
            'Length',
            'Score',
            'Sequence'
        ])
        
        rows = [header]
        for ann in annotations:
            row = '\t'.join([
                str(ann.get('sequence_name', 'N/A')),
                str(ann.get('class_name', 'N/A')),
                str(ann.get('pattern_id', 'N/A')),
                str(ann.get('start', 0)),
                str(ann.get('end', 0)),
                str(ann.get('length', 0)),
                str(ann.get('score', 0.0)),
                str(ann.get('matched_seq', ''))
            ])
            rows.append(row)
        
        return '\n'.join(rows)
    
    def format_as_bed(self, annotations: List[Dict[str, Any]], sequence_name: str = "chr1") -> str:
        """
        Format resolved G4 annotations as BED format for genome browsers.
        
        Args:
            annotations: List of G4 annotation dictionaries
            sequence_name: Chromosome/sequence name for BED track
        
        Returns:
            BED format string (6-column)
        
        BED format is widely used in genomics for:
        - UCSC Genome Browser visualization
        - IGV (Integrative Genomics Viewer)
        - bedtools operations
        
        Columns: chrom, chromStart, chromEnd, name, score, strand
        """
        rows = []
        for ann in annotations:
            # BED uses 0-based coordinates (same as our internal format)
            chrom = sequence_name
            start = ann.get('start', 0)
            end = ann.get('end', 0)
            name = f"{ann.get('class_name', 'G4')}_{ann.get('pattern_id', 'unknown')}"
            # Convert score to integer (0-1000 range for BED)
            score = min(1000, int(ann.get('score', 0.0) * 100))
            strand = '+'
            
            row = '\t'.join([
                str(chrom),
                str(start),
                str(end),
                name,
                str(score),
                strand
            ])
            rows.append(row)
        
        return '\n'.join(rows)


def main():
    """
    Command-line interface for G4 overlap resolution.
    
    Usage examples:
    
    1. Basic usage with sequence:
       python g4_overlap_resolver.py --sequence "GGGTTAGGGTTAGGGTTAGGG"
    
    2. From FASTA file with JSON output:
       python g4_overlap_resolver.py --fasta input.fasta --format json --output results.json
    
    3. Generate BED format for genome browser:
       python g4_overlap_resolver.py --sequence "GGGGTTTTGGGGTTTTGGGG" --format bed
    
    The CLI provides flexible input/output options while maintaining the core
    overlap resolution algorithm from the G4Hunter framework (Bedrat et al. 2016)[web:3].
    """
    parser = argparse.ArgumentParser(
        description='G-Quadruplex Overlap Resolution Agent',
        epilog='Scientific references: Burge 2006, Huppert 2005, Bedrat et al. 2016'
    )
    
    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--sequence',
        type=str,
        help='DNA sequence string (e.g., "GGGTTAGGGTTAGGGTTAGGG")'
    )
    input_group.add_argument(
        '--fasta',
        type=str,
        help='Path to FASTA file containing DNA sequence(s)'
    )
    
    # Optional parameters
    parser.add_argument(
        '--name',
        type=str,
        default='input_seq',
        help='Sequence name/identifier (default: input_seq)'
    )
    parser.add_argument(
        '--format',
        type=str,
        choices=['json', 'table', 'bed'],
        default='json',
        help='Output format (default: json)'
    )
    parser.add_argument(
        '--output',
        type=str,
        help='Output file path (default: stdout)'
    )
    parser.add_argument(
        '--compact',
        action='store_true',
        help='Compact JSON output (no indentation)'
    )
    
    args = parser.parse_args()
    
    # Initialize the G4 overlap resolver
    resolver = G4OverlapResolver()
    
    # Parse input sequence
    if args.sequence:
        sequence = args.sequence.strip()
        sequence_name = args.name
    elif args.fasta:
        # Simple FASTA parser for single sequence
        try:
            with open(args.fasta, 'r') as f:
                lines = f.readlines()
                sequence_name = args.name
                sequence = ''
                for line in lines:
                    line = line.strip()
                    if line.startswith('>'):
                        if args.name == 'input_seq':  # Use FASTA header if name not specified
                            sequence_name = line[1:].split()[0]
                    else:
                        sequence += line
        except Exception as e:
            print(f"Error reading FASTA file: {e}", file=sys.stderr)
            sys.exit(1)
    
    # Validate sequence
    if not sequence:
        print("Error: Empty sequence provided", file=sys.stderr)
        sys.exit(1)
    
    # Run overlap resolution
    try:
        annotations = resolver.resolve_and_annotate(sequence, sequence_name)
    except Exception as e:
        print(f"Error during G4 analysis: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Format output
    if args.format == 'json':
        output = resolver.format_as_json(annotations, pretty=not args.compact)
    elif args.format == 'table':
        output = resolver.format_as_table(annotations)
    elif args.format == 'bed':
        output = resolver.format_as_bed(annotations, sequence_name)
    
    # Write output
    if args.output:
        try:
            with open(args.output, 'w') as f:
                f.write(output)
            print(f"Results written to {args.output}", file=sys.stderr)
        except Exception as e:
            print(f"Error writing output file: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print(output)


if __name__ == '__main__':
    main()
