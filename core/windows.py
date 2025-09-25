"""
Sequence Windowing and Chunking for Large-Scale Analysis
========================================================

Handles chunking of large genomic sequences (100kb default) with overlap
management and merge operations for seamless motif detection across boundaries.
"""

import numpy as np
from typing import List, Tuple, Iterator, Dict, Any, Optional
from dataclasses import dataclass

@dataclass
class SequenceWindow:
    """Represents a sequence window with metadata."""
    sequence: str
    start_pos: int
    end_pos: int
    chunk_id: int
    overlap_start: bool = False
    overlap_end: bool = False

class SequenceWindower:
    """
    Manages chunking of large sequences with overlap handling.
    """
    
    def __init__(self, chunk_size: int = 100000, overlap_size: int = 1000):
        """
        Initialize windower with chunking parameters.
        
        Args:
            chunk_size: Size of each chunk in base pairs
            overlap_size: Overlap between adjacent chunks
        """
        self.chunk_size = chunk_size
        self.overlap_size = overlap_size
        
    def create_windows(self, sequence: str, sequence_name: str = "sequence") -> List[SequenceWindow]:
        """
        Create overlapping windows from a large sequence.
        
        Args:
            sequence: Input DNA sequence
            sequence_name: Name identifier for the sequence
            
        Returns:
            List of SequenceWindow objects
        """
        seq_len = len(sequence)
        windows = []
        
        if seq_len <= self.chunk_size:
            # Single window for short sequences
            windows.append(SequenceWindow(
                sequence=sequence,
                start_pos=0,
                end_pos=seq_len,
                chunk_id=0
            ))
            return windows
        
        chunk_id = 0
        start = 0
        
        while start < seq_len:
            end = min(start + self.chunk_size, seq_len)
            
            # Determine overlap flags
            overlap_start = start > 0
            overlap_end = end < seq_len
            
            # Extract sequence with overlaps
            window_start = max(0, start - (self.overlap_size if overlap_start else 0))
            window_end = min(seq_len, end + (self.overlap_size if overlap_end else 0))
            
            window_seq = sequence[window_start:window_end]
            
            windows.append(SequenceWindow(
                sequence=window_seq,
                start_pos=window_start,
                end_pos=window_end,
                chunk_id=chunk_id,
                overlap_start=overlap_start,
                overlap_end=overlap_end
            ))
            
            chunk_id += 1
            start = end - self.overlap_size  # Step with overlap
            
        return windows
    
    def merge_overlapping_motifs(self, motif_lists: List[List[Dict[str, Any]]], 
                                windows: List[SequenceWindow]) -> List[Dict[str, Any]]:
        """
        Merge motifs from overlapping windows, handling duplicates at boundaries.
        
        Args:
            motif_lists: List of motif lists from each window
            windows: Corresponding window objects
            
        Returns:
            Merged and deduplicated motif list
        """
        all_motifs = []
        
        for motifs, window in zip(motif_lists, windows):
            for motif in motifs:
                # Adjust positions to global coordinates
                global_start = motif['Start'] + window.start_pos
                global_end = motif['End'] + window.start_pos
                
                # Create adjusted motif
                adjusted_motif = motif.copy()
                adjusted_motif['Start'] = global_start
                adjusted_motif['End'] = global_end
                adjusted_motif['Chunk_ID'] = window.chunk_id
                
                all_motifs.append(adjusted_motif)
        
        # Remove duplicates at overlap boundaries
        return self._deduplicate_boundary_motifs(all_motifs)
    
    def _deduplicate_boundary_motifs(self, motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Remove duplicate motifs found in overlapping regions.
        
        Args:
            motifs: List of motifs with global coordinates
            
        Returns:
            Deduplicated motif list
        """
        if not motifs:
            return []
        
        # Sort by position for efficient comparison
        sorted_motifs = sorted(motifs, key=lambda x: (x['Start'], x['End']))
        deduplicated = [sorted_motifs[0]]
        
        for current in sorted_motifs[1:]:
            last = deduplicated[-1]
            
            # Check for duplicates (same class, overlapping positions)
            if (current['Class'] == last['Class'] and
                current['Subclass'] == last['Subclass'] and
                self._positions_overlap(current, last, tolerance=5)):
                
                # Keep the motif with higher score
                if current.get('Score', 0) > last.get('Score', 0):
                    deduplicated[-1] = current
                # If scores are equal, keep the first one (skip current)
                    
            else:
                deduplicated.append(current)
        
        return deduplicated
    
    def _positions_overlap(self, motif1: Dict[str, Any], motif2: Dict[str, Any], 
                          tolerance: int = 5) -> bool:
        """
        Check if two motifs overlap within a tolerance.
        
        Args:
            motif1, motif2: Motif dictionaries with Start/End positions
            tolerance: Maximum difference in positions to consider overlap
            
        Returns:
            True if motifs overlap within tolerance
        """
        start_diff = abs(motif1['Start'] - motif2['Start'])
        end_diff = abs(motif1['End'] - motif2['End'])
        
        return start_diff <= tolerance and end_diff <= tolerance

class GenomeWindower(SequenceWindower):
    """
    Extended windower for whole-genome analysis with chromosome handling.
    """
    
    def __init__(self, chunk_size: int = 1000000, overlap_size: int = 10000):
        """
        Initialize genome windower with larger default chunks.
        
        Args:
            chunk_size: Size of each chunk (default 1Mb)
            overlap_size: Overlap between chunks (default 10kb)
        """
        super().__init__(chunk_size, overlap_size)
    
    def create_chromosome_windows(self, fasta_dict: Dict[str, str]) -> Dict[str, List[SequenceWindow]]:
        """
        Create windows for multiple chromosomes/contigs.
        
        Args:
            fasta_dict: Dictionary mapping chromosome names to sequences
            
        Returns:
            Dictionary mapping chromosome names to window lists
        """
        chromosome_windows = {}
        
        for chrom_name, sequence in fasta_dict.items():
            windows = self.create_windows(sequence, chrom_name)
            
            # Add chromosome information to windows
            for window in windows:
                window.chromosome = chrom_name
                
            chromosome_windows[chrom_name] = windows
            
        return chromosome_windows

def calculate_window_statistics(window: SequenceWindow) -> Dict[str, Any]:
    """
    Calculate basic statistics for a sequence window.
    
    Args:
        window: SequenceWindow object
        
    Returns:
        Dictionary of window statistics
    """
    sequence = window.sequence.upper()
    length = len(sequence)
    
    if length == 0:
        return {'Length': 0, 'GC_Content': 0, 'N_Content': 0}
    
    gc_count = sequence.count('G') + sequence.count('C')
    n_count = sequence.count('N')
    
    return {
        'Window_ID': window.chunk_id,
        'Start_Position': window.start_pos,
        'End_Position': window.end_pos,
        'Length': length,
        'GC_Content': round(100 * gc_count / length, 2),
        'N_Content': round(100 * n_count / length, 2),
        'Has_Overlap_Start': window.overlap_start,
        'Has_Overlap_End': window.overlap_end
    }

__all__ = [
    'SequenceWindow',
    'SequenceWindower', 
    'GenomeWindower',
    'calculate_window_statistics'
]