#!/usr/bin/env python3
"""
Basic visualization and output formatting for NonBScanner results
- Simple text-based visualization of motif positions
- Summary statistics and distribution plots
- Export utilities for downstream analysis
- Basic graphical output for motif density and scores

This module provides lightweight visualization capabilities for NonBScanner
results without requiring heavy dependencies like matplotlib or plotly.

Author: Dr. Venkata Rajesh Yella
Integration: 2024 - NBDFinder Visualization
"""

import sys
import os
from typing import List, Dict, Any, Optional, Tuple
from collections import defaultdict
import math


def create_text_visualization(motifs: List[Dict[str, Any]], 
                            sequence_length: int,
                            width: int = 80) -> str:
    """
    Create text-based visualization of motif positions along sequence.
    
    | Parameter       | Type | Description                    | Default | Range      |
    |-----------------|------|--------------------------------|---------|------------|
    | motifs          | list | List of motif dictionaries    | -       | any        |
    | sequence_length | int  | Total length of sequence       | -       | >0         |
    | width           | int  | Width of text display          | 80      | 20-200     |
    | return          | str  | Text visualization string      | -       | any        |
    """
    if not motifs or sequence_length <= 0:
        return "No motifs to visualize"
    
    # Create scale
    scale_factor = sequence_length / width
    visualization = []
    
    # Header with scale
    visualization.append("Motif Distribution Visualization")
    visualization.append("=" * width)
    
    # Scale markers
    scale_line = ""
    number_line = ""
    for i in range(width + 1):
        pos = int(i * scale_factor) + 1
        if i % 10 == 0:
            scale_line += "|"
            number_line += str(pos).ljust(10)[:10]
        else:
            scale_line += "-"
            number_line += " "
    
    visualization.append(number_line[:width])
    visualization.append(scale_line[:width])
    
    # Group motifs by class
    motif_groups = defaultdict(list)
    for motif in motifs:
        motif_class = motif.get('Class', 'Unknown')
        motif_groups[motif_class].append(motif)
    
    # Symbol mapping for different motif types
    symbols = {
        'A-philic DNA': 'A',
        'Z-DNA': 'Z',
        'Unknown': '?'
    }
    
    # Create visualization lines for each motif type
    for motif_class, class_motifs in motif_groups.items():
        line = [' '] * width
        symbol = symbols.get(motif_class, motif_class[0].upper())
        
        for motif in class_motifs:
            start = motif.get('Start', 1)
            end = motif.get('End', start)
            
            # Convert to display coordinates
            start_pos = min(width - 1, int((start - 1) / scale_factor))
            end_pos = min(width - 1, int((end - 1) / scale_factor))
            
            # Fill in the motif region
            for pos in range(start_pos, end_pos + 1):
                if pos < width:
                    line[pos] = symbol
        
        line_str = ''.join(line)
        visualization.append(f"{motif_class:12}: {line_str}")
    
    # Legend
    visualization.append("")
    visualization.append("Legend:")
    for motif_class, symbol in symbols.items():
        if motif_class in motif_groups:
            count = len(motif_groups[motif_class])
            visualization.append(f"  {symbol} = {motif_class} ({count} motifs)")
    
    return '\n'.join(visualization)


def generate_statistics_report(motifs: List[Dict[str, Any]], 
                             sequences: Optional[List[Dict]] = None) -> str:
    """
    Generate detailed statistics report for detected motifs.
    
    | Parameter | Type | Description                    | Default | Range      |
    |-----------|------|--------------------------------|---------|------------|
    | motifs    | list | List of motif dictionaries    | -       | any        |
    | sequences | list | List of sequence info dicts   | None    | any        |
    | return    | str  | Formatted statistics report    | -       | any        |
    """
    if not motifs:
        return "No motifs found for statistics"
    
    report = []
    report.append("NonBScanner Statistics Report")
    report.append("=" * 50)
    report.append("")
    
    # Overall statistics
    total_motifs = len(motifs)
    report.append(f"Total motifs detected: {total_motifs}")
    
    # Motif class distribution
    class_counts = defaultdict(int)
    subclass_counts = defaultdict(int)
    
    for motif in motifs:
        motif_class = motif.get('Class', 'Unknown')
        subclass = motif.get('Subclass', 'Unknown')
        class_counts[motif_class] += 1
        subclass_counts[f"{motif_class}.{subclass}"] += 1
    
    report.append("")
    report.append("Motif Class Distribution:")
    report.append("-" * 25)
    for motif_class, count in sorted(class_counts.items()):
        percentage = (count / total_motifs) * 100
        report.append(f"  {motif_class}: {count} ({percentage:.1f}%)")
    
    # Score statistics
    scores = [motif.get('Raw_Score', motif.get('Score', 0)) for motif in motifs]
    scores = [s for s in scores if isinstance(s, (int, float))]
    
    if scores:
        report.append("")
        report.append("Score Statistics:")
        report.append("-" * 17)
        report.append(f"  Mean score: {sum(scores)/len(scores):.2f}")
        report.append(f"  Min score:  {min(scores):.2f}")
        report.append(f"  Max score:  {max(scores):.2f}")
        
        # Score distribution (simple histogram)
        report.append("")
        report.append("Score Distribution:")
        score_bins = create_histogram(scores, bins=10)
        for bin_range, count in score_bins:
            bar = '#' * min(50, int(count * 50 / max(1, max(c for _, c in score_bins))))
            report.append(f"  {bin_range[0]:6.1f}-{bin_range[1]:6.1f}: {bar} ({count})")
    
    # Length statistics
    lengths = [motif.get('Length', 0) for motif in motifs]
    lengths = [l for l in lengths if isinstance(l, int) and l > 0]
    
    if lengths:
        report.append("")
        report.append("Length Statistics:")
        report.append("-" * 18)
        report.append(f"  Mean length: {sum(lengths)/len(lengths):.1f} bp")
        report.append(f"  Min length:  {min(lengths)} bp")
        report.append(f"  Max length:  {max(lengths)} bp")
    
    # Confidence distribution
    confidence_counts = defaultdict(int)
    for motif in motifs:
        confidence = motif.get('Confidence', 'Unknown')
        confidence_counts[confidence] += 1
    
    if confidence_counts:
        report.append("")
        report.append("Confidence Distribution:")
        report.append("-" * 24)
        for confidence, count in sorted(confidence_counts.items()):
            percentage = (count / total_motifs) * 100
            report.append(f"  {confidence}: {count} ({percentage:.1f}%)")
    
    # Per-sequence statistics if available
    if sequences:
        report.append("")
        report.append("Per-Sequence Statistics:")
        report.append("-" * 25)
        
        for seq_info in sequences:
            seq_name = seq_info.get('name', 'Unknown')
            seq_length = len(seq_info.get('sequence', ''))
            seq_motifs = [m for m in motifs if m.get('Sequence_Name', '') == seq_name]
            
            if seq_motifs:
                density = (len(seq_motifs) * 1000.0) / max(1, seq_length)
                report.append(f"  {seq_name}:")
                report.append(f"    Length: {seq_length:,} bp")
                report.append(f"    Motifs: {len(seq_motifs)}")
                report.append(f"    Density: {density:.2f} motifs/kb")
    
    return '\n'.join(report)


def create_histogram(values: List[float], bins: int = 10) -> List[Tuple[Tuple[float, float], int]]:
    """
    Create histogram bins for numeric values.
    
    | Parameter | Type | Description                    | Default | Range      |
    |-----------|------|--------------------------------|---------|------------|
    | values    | list | List of numeric values         | -       | any        |
    | bins      | int  | Number of histogram bins       | 10      | 1-100      |
    | return    | list | List of (bin_range, count) tuples | []  | any        |
    """
    if not values:
        return []
    
    min_val = min(values)
    max_val = max(values)
    
    if min_val == max_val:
        return [((min_val, max_val), len(values))]
    
    bin_width = (max_val - min_val) / bins
    histogram = []
    
    for i in range(bins):
        bin_start = min_val + i * bin_width
        bin_end = min_val + (i + 1) * bin_width
        
        count = sum(1 for v in values if bin_start <= v < bin_end)
        
        # Handle the last bin to include max value
        if i == bins - 1:
            count = sum(1 for v in values if bin_start <= v <= bin_end)
        
        histogram.append(((bin_start, bin_end), count))
    
    return histogram


def create_motif_summary_table(motifs: List[Dict[str, Any]]) -> str:
    """
    Create formatted table summary of top motifs.
    
    | Parameter | Type | Description                    | Range      |
    |-----------|------|--------------------------------|------------|
    | motifs    | list | List of motif dictionaries    | any        |
    | return    | str  | Formatted table string         | any        |
    """
    if not motifs:
        return "No motifs to summarize"
    
    # Sort by score (descending)
    sorted_motifs = sorted(motifs, 
                          key=lambda x: x.get('Raw_Score', x.get('Score', 0)), 
                          reverse=True)
    
    # Take top 20 motifs
    top_motifs = sorted_motifs[:20]
    
    table = []
    table.append("Top Detected Motifs")
    table.append("=" * 80)
    table.append("")
    
    # Header
    header = f"{'Seq':12} {'Class':10} {'Start':8} {'End':8} {'Length':8} {'Score':8} {'Conf':10}"
    table.append(header)
    table.append("-" * 80)
    
    # Data rows
    for motif in top_motifs:
        seq_name = motif.get('Sequence_Name', 'Unknown')[:12]
        motif_class = motif.get('Class', 'Unknown')[:10]
        start = motif.get('Start', 0)
        end = motif.get('End', 0)
        length = motif.get('Length', 0)
        score = motif.get('Raw_Score', motif.get('Score', 0))
        confidence = motif.get('Confidence', 'Unknown')[:10]
        
        row = f"{seq_name:12} {motif_class:10} {start:8d} {end:8d} {length:8d} {score:8.2f} {confidence:10}"
        table.append(row)
    
    if len(motifs) > 20:
        table.append("")
        table.append(f"... and {len(motifs) - 20} more motifs")
    
    return '\n'.join(table)


def export_for_genome_browser(motifs: List[Dict[str, Any]], 
                            output_file: str, 
                            format_type: str = 'bed') -> str:
    """
    Export motifs in genome browser compatible format.
    
    | Parameter   | Type | Description                    | Default | Range      |
    |-------------|------|--------------------------------|---------|------------|
    | motifs      | list | List of motif dictionaries    | -       | any        |
    | output_file | str  | Output file path               | -       | any        |
    | format_type | str  | Export format type             | 'bed'   | bed/gff    |
    | return      | str  | Status message                 | -       | any        |
    """
    if not motifs:
        return "No motifs to export"
    
    try:
        with open(output_file, 'w') as f:
            if format_type.lower() == 'bed':
                # BED format: chrom start end name score strand
                f.write("track name=\"NonBScanner\" description=\"A-philic and Z-DNA motifs\"\n")
                
                for motif in motifs:
                    chrom = motif.get('Sequence_Name', 'chr1')
                    start = motif.get('Start', 1) - 1  # BED uses 0-based coordinates
                    end = motif.get('End', start + 1)
                    name = f"{motif.get('Class', 'Unknown')}_{motif.get('Subclass', '')}"
                    score = int(motif.get('Raw_Score', motif.get('Score', 0)) * 10)  # Scale for BED
                    strand = "."
                    
                    f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
            
            elif format_type.lower() == 'gff':
                # GFF3 format
                f.write("##gff-version 3\n")
                
                for i, motif in enumerate(motifs, 1):
                    seqid = motif.get('Sequence_Name', 'sequence')
                    source = "NonBScanner"
                    feature_type = motif.get('Class', 'motif').replace(' ', '_')
                    start = motif.get('Start', 1)
                    end = motif.get('End', start)
                    score = motif.get('Raw_Score', motif.get('Score', '.'))
                    strand = "."
                    phase = "."
                    
                    # Attributes
                    attributes = []
                    attributes.append(f"ID=motif_{i}")
                    attributes.append(f"Name={motif.get('Class', 'Unknown')}")
                    if 'Subclass' in motif:
                        attributes.append(f"subclass={motif['Subclass']}")
                    if 'Classification' in motif:
                        attributes.append(f"classification={motif['Classification']}")
                    
                    attr_str = ";".join(attributes)
                    
                    f.write(f"{seqid}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attr_str}\n")
        
        return f"Exported {len(motifs)} motifs to {output_file} in {format_type.upper()} format"
    
    except Exception as e:
        return f"Export failed: {e}"


def create_density_plot(motifs: List[Dict[str, Any]], 
                       sequence_length: int,
                       window_size: int = 1000) -> str:
    """
    Create text-based density plot showing motif frequency along sequence.
    
    | Parameter       | Type | Description                    | Default | Range      |
    |-----------------|------|--------------------------------|---------|------------|
    | motifs          | list | List of motif dictionaries    | -       | any        |
    | sequence_length | int  | Total sequence length          | -       | >0         |
    | window_size     | int  | Window size for density calc   | 1000    | 100-10000  |
    | return          | str  | Text density plot              | -       | any        |
    """
    if not motifs or sequence_length <= 0:
        return "No data for density plot"
    
    # Calculate density in windows
    num_windows = max(1, sequence_length // window_size)
    densities = [0] * num_windows
    
    for motif in motifs:
        start = motif.get('Start', 1)
        window_idx = min(num_windows - 1, (start - 1) // window_size)
        densities[window_idx] += 1
    
    # Create plot
    plot_lines = []
    plot_lines.append("Motif Density Plot")
    plot_lines.append("=" * 50)
    plot_lines.append(f"Window size: {window_size:,} bp")
    plot_lines.append("")
    
    max_density = max(densities) if densities else 1
    
    for i, density in enumerate(densities):
        start_pos = i * window_size + 1
        end_pos = min(sequence_length, (i + 1) * window_size)
        
        # Scale bar length
        bar_length = int((density / max_density) * 40) if max_density > 0 else 0
        bar = '#' * bar_length
        
        plot_lines.append(f"{start_pos:8d}-{end_pos:8d}: {bar} ({density})")
    
    return '\n'.join(plot_lines)


# === Example usage and testing ===
if __name__ == "__main__":
    # Test visualization functions with sample data
    sample_motifs = [
        {
            'Sequence_Name': 'test_seq',
            'Class': 'A-philic DNA',
            'Subclass': 'A-philic',
            'Start': 100,
            'End': 120,
            'Length': 21,
            'Raw_Score': 2.5,
            'Score': 2.5,
            'Confidence': 'High'
        },
        {
            'Sequence_Name': 'test_seq',
            'Class': 'Z-DNA', 
            'Subclass': 'Z-DNA',
            'Start': 200,
            'End': 215,
            'Length': 16,
            'Raw_Score': 55.0,
            'Score': 55.0,
            'Confidence': 'Moderate'
        }
    ]
    
    print("Text Visualization:")
    print(create_text_visualization(sample_motifs, 500))
    print("\n" + "="*50 + "\n")
    
    print("Statistics Report:")
    print(generate_statistics_report(sample_motifs))
    print("\n" + "="*50 + "\n")
    
    print("Motif Summary Table:")
    print(create_motif_summary_table(sample_motifs))
    print("\n" + "="*50 + "\n")
    
    print("Density Plot:")
    print(create_density_plot(sample_motifs, 500, 100))