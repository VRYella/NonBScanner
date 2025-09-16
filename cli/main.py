"""
NBDFinder Command Line Interface
===============================

Main CLI for NBDFinder with comprehensive motif detection capabilities.
"""

import argparse
import sys
import time
from pathlib import Path
from typing import List, Optional, Dict, Any
import json

from ..orchestrators import detect_all_motifs, StreamingOrchestrator
from ..io import FastaReader, export_to_bed, export_to_csv, export_to_gff3
from ..io.schemas import AnalysisConfig, ExportConfig
from ..core import get_basic_stats

def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="NBDFinder - Non-B DNA Structure Detection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis
  nbdfinder run --fasta sequence.fa --output results.bed
  
  # Comprehensive analysis with all formats
  nbdfinder run --fasta genome.fa --output-dir results/ --format bed,csv,gff
  
  # Streaming analysis for large files
  nbdfinder run --fasta large_genome.fa --streaming --chunk-size 1000000
  
  # Filter by motif classes
  nbdfinder run --fasta seq.fa --classes G-Quadruplex,i-Motif --min-score 0.5
  
  # Slice FASTA regions
  nbdfinder slice --fasta genome.fa --regions chr1:1000-2000,chr2:5000-6000
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Main run command
    run_parser = subparsers.add_parser('run', help='Run motif detection analysis')
    _add_run_arguments(run_parser)
    
    # FASTA slicing command
    slice_parser = subparsers.add_parser('slice', help='Extract FASTA regions')
    _add_slice_arguments(slice_parser)
    
    # Version command
    version_parser = subparsers.add_parser('version', help='Show version information')
    
    args = parser.parse_args()
    
    if args.command == 'run':
        return run_analysis(args)
    elif args.command == 'slice':
        return slice_fasta(args)
    elif args.command == 'version':
        return show_version()
    else:
        parser.print_help()
        return 1

def _add_run_arguments(parser):
    """Add arguments for run command."""
    # Input/Output
    parser.add_argument('--fasta', '-f', required=True, type=Path,
                       help='Input FASTA file')
    parser.add_argument('--output', '-o', type=Path,
                       help='Output file (extension determines format)')
    parser.add_argument('--output-dir', '-d', type=Path,
                       help='Output directory for multiple files')
    parser.add_argument('--format', default='bed',
                       help='Output format(s): bed,csv,gff,json (comma-separated)')
    
    # Analysis parameters
    parser.add_argument('--classes', type=str,
                       help='Motif classes to detect (comma-separated)')
    parser.add_argument('--min-score', type=float, default=0.1,
                       help='Minimum score threshold')
    parser.add_argument('--min-length', type=int, default=10,
                       help='Minimum motif length')
    parser.add_argument('--max-length', type=int, default=1000,
                       help='Maximum motif length')
    
    # Processing options
    parser.add_argument('--streaming', action='store_true',
                       help='Use streaming mode for large files')
    parser.add_argument('--chunk-size', type=int, default=100000,
                       help='Chunk size for streaming mode')
    parser.add_argument('--workers', type=int,
                       help='Number of parallel workers')
    parser.add_argument('--memory-limit', type=int, default=8,
                       help='Memory limit in GB')
    
    # Filtering options
    parser.add_argument('--no-overlaps', action='store_true',
                       help='Remove overlapping motifs')
    parser.add_argument('--merge-nearby', action='store_true',
                       help='Merge nearby motifs')
    parser.add_argument('--merge-distance', type=int, default=10,
                       help='Maximum distance for merging motifs')
    
    # Output options
    parser.add_argument('--include-sequence', action='store_true', default=True,
                       help='Include sequences in output')
    parser.add_argument('--score-type', choices=['raw', 'normalized'], default='normalized',
                       help='Score type for output')
    parser.add_argument('--precision', type=int, default=3,
                       help='Decimal precision for scores')
    
    # Verbosity
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Quiet mode')

def _add_slice_arguments(parser):
    """Add arguments for slice command."""
    parser.add_argument('--fasta', '-f', required=True, type=Path,
                       help='Input FASTA file')
    parser.add_argument('--regions', '-r', required=True,
                       help='Regions to extract (format: chr:start-end)')
    parser.add_argument('--output', '-o', type=Path,
                       help='Output FASTA file')
    parser.add_argument('--format', choices=['fasta', 'bed'], default='fasta',
                       help='Output format')

def run_analysis(args) -> int:
    """Run motif detection analysis."""
    try:
        # Validate input
        if not args.fasta.exists():
            print(f"Error: FASTA file not found: {args.fasta}")
            return 1
        
        # Create analysis configuration
        config = AnalysisConfig(
            chunk_size=args.chunk_size,
            max_workers=args.workers,
            min_score_threshold=args.min_score,
            min_motif_length=args.min_length,
            max_motif_length=args.max_length,
            remove_overlaps=args.no_overlaps,
            merge_nearby=args.merge_nearby,
            merge_distance=args.merge_distance
        )
        
        # Filter enabled classes
        if args.classes:
            enabled_classes = [cls.strip() for cls in args.classes.split(',')]
            config.enabled_classes = enabled_classes
        
        if not args.quiet:
            print(f"Analyzing FASTA file: {args.fasta}")
            if args.streaming:
                print("Using streaming mode for large file processing")
        
        start_time = time.time()
        
        if args.streaming:
            # Use streaming orchestrator
            orchestrator = StreamingOrchestrator(config)
            
            if args.verbose:
                def progress_callback(info):
                    stage = info.get('stage', 'unknown')
                    if stage == 'processing':
                        percent = info.get('progress_percent', 0)
                        print(f"\rProgress: {percent:.1f}%", end='', flush=True)
                    elif stage == 'complete':
                        print(f"\nCompleted: {info.get('total_motifs', 0)} motifs found")
                
                orchestrator.add_progress_callback(progress_callback)
            
            # Process sequences
            all_results = []
            with FastaReader(args.fasta) as reader:
                for seq_name in reader.get_sequence_names():
                    sequence = reader.get_sequence(seq_name)
                    result = orchestrator.analyze_sequence_stream(sequence, seq_name)
                    all_results.append(result)
            
            # Combine all motifs
            all_motifs = []
            for result in all_results:
                all_motifs.extend(result.motifs)
        
        else:
            # Standard analysis
            all_motifs = []
            with FastaReader(args.fasta) as reader:
                for seq_name, sequence in reader.iterate_sequences():
                    if args.verbose:
                        print(f"Processing sequence: {seq_name}")
                    
                    motifs = detect_all_motifs(sequence, seq_name)
                    all_motifs.extend(motifs)
        
        processing_time = time.time() - start_time
        
        if not args.quiet:
            print(f"Analysis completed in {processing_time:.2f} seconds")
            print(f"Found {len(all_motifs)} motifs")
        
        # Export results
        output_formats = [fmt.strip() for fmt in args.format.split(',')]
        
        for fmt in output_formats:
            if args.output_dir:
                output_path = args.output_dir / f"nbdfinder_results.{fmt}"
                args.output_dir.mkdir(exist_ok=True)
            elif args.output:
                output_path = args.output
            else:
                output_path = Path(f"nbdfinder_results.{fmt}")
            
            _export_results(all_motifs, fmt, output_path, args)
            
            if not args.quiet:
                print(f"Results exported to: {output_path}")
        
        return 0
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1

def _export_results(motifs: List[Dict[str, Any]], format_type: str, 
                   output_path: Path, args) -> None:
    """Export results in specified format."""
    if format_type == 'bed':
        content = export_to_bed(
            motifs,
            score_type=args.score_type,
            include_subclass=True
        )
    elif format_type == 'csv':
        content = export_to_csv(
            motifs,
            include_sequence=args.include_sequence,
            precision=args.precision
        )
    elif format_type == 'gff':
        content = export_to_gff3(motifs)
    elif format_type == 'json':
        content = json.dumps(motifs, indent=2)
    else:
        raise ValueError(f"Unsupported format: {format_type}")
    
    with open(output_path, 'w') as f:
        f.write(content)

def slice_fasta(args) -> int:
    """Extract FASTA regions."""
    try:
        if not args.fasta.exists():
            print(f"Error: FASTA file not found: {args.fasta}")
            return 1
        
        # Parse regions
        regions = []
        for region_str in args.regions.split(','):
            region_str = region_str.strip()
            if ':' in region_str and '-' in region_str:
                chrom, coords = region_str.split(':')
                start, end = map(int, coords.split('-'))
                regions.append((chrom, start, end))
            else:
                print(f"Invalid region format: {region_str}")
                return 1
        
        # Extract sequences
        extracted = {}
        with FastaReader(args.fasta) as reader:
            for chrom, start, end in regions:
                try:
                    sequence = reader.get_sequence_slice(chrom, start-1, end)  # Convert to 0-based
                    region_id = f"{chrom}:{start}-{end}"
                    extracted[region_id] = sequence
                except KeyError:
                    print(f"Warning: Sequence '{chrom}' not found")
        
        # Output results
        if args.output:
            output_path = args.output
        else:
            output_path = Path("extracted_regions.fasta")
        
        if args.format == 'fasta':
            with open(output_path, 'w') as f:
                for region_id, sequence in extracted.items():
                    f.write(f">{region_id}\n")
                    # Write with line breaks
                    for i in range(0, len(sequence), 80):
                        f.write(sequence[i:i+80] + '\n')
        
        elif args.format == 'bed':
            with open(output_path, 'w') as f:
                for region_id, sequence in extracted.items():
                    chrom, coords = region_id.split(':')
                    start, end = map(int, coords.split('-'))
                    f.write(f"{chrom}\t{start-1}\t{end}\t{region_id}\t0\t+\n")
        
        print(f"Extracted {len(extracted)} regions to: {output_path}")
        return 0
        
    except Exception as e:
        print(f"Error during extraction: {e}")
        return 1

def show_version() -> int:
    """Show version information."""
    from .. import __version__, __author__
    print(f"NBDFinder version {__version__}")
    print(f"Author: {__author__}")
    return 0

if __name__ == '__main__':
    sys.exit(main())