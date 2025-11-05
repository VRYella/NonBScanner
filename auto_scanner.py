"""
Backward-compatible wrapper for genome-scale analysis
======================================================

This module provides a drop-in replacement for analyze_sequence that
automatically uses genome-scale optimization for large sequences.
"""

from typing import List, Dict, Any


def analyze_sequence_auto(sequence: str, sequence_name: str = "sequence", 
                          auto_optimize: bool = True) -> List[Dict[str, Any]]:
    """
    Automatically choose between standard and genome-scale analysis
    
    Args:
        sequence: DNA sequence to analyze
        sequence_name: Name identifier
        auto_optimize: If True, automatically use genome-scale for large sequences
        
    Returns:
        List of detected motifs
    """
    seq_len = len(sequence)
    
    # Use genome-scale for sequences > 10MB
    if auto_optimize and seq_len > 10_000_000:
        print(f"Auto-selecting genome-scale scanner for {seq_len:,} bp sequence")
        from genome_scale_scanner import analyze_genome_sequence
        return analyze_genome_sequence(
            sequence, 
            sequence_name,
            enable_hybrid_cluster=False
        )
    else:
        # Use standard scanner
        from scanner import analyze_sequence
        return analyze_sequence(sequence, sequence_name)


# Alias for convenience
analyze = analyze_sequence_auto
