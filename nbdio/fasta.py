"""
FASTA File Handling with Memory Mapping
======================================

Efficient FASTA file reading and sequence slicing utilities using
memory mapping for large genomic files.
"""

import mmap
import os
from typing import Dict, List, Tuple, Iterator, Optional, Union
from pathlib import Path
import re

class FastaReader:
    """
    Memory-mapped FASTA file reader for efficient large file handling.
    """
    
    def __init__(self, fasta_path: Union[str, Path]):
        """
        Initialize FASTA reader with memory mapping.
        
        Args:
            fasta_path: Path to FASTA file
        """
        self.fasta_path = Path(fasta_path)
        self.file_handle = None
        self.mmap_obj = None
        self.sequences = {}
        self.sequence_indices = {}
        
        if not self.fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
        
        self._open_and_index()
    
    def _open_and_index(self):
        """Open file with memory mapping and create sequence index."""
        self.file_handle = open(self.fasta_path, 'rb')
        self.mmap_obj = mmap.mmap(self.file_handle.fileno(), 0, access=mmap.ACCESS_READ)
        
        # Index sequences for fast access
        self._build_sequence_index()
    
    def _build_sequence_index(self):
        """Build index of sequence positions in the memory-mapped file."""
        current_name = None
        current_start = 0
        
        # Find all header lines
        header_pattern = re.compile(rb'^>([^\s]+)', re.MULTILINE)
        headers = []
        
        for match in header_pattern.finditer(self.mmap_obj):
            header_name = match.group(1).decode('ascii')
            header_pos = match.start()
            headers.append((header_name, header_pos))
        
        # Calculate sequence regions
        for i, (name, header_pos) in enumerate(headers):
            # Find start of sequence (after header line)
            seq_start = self.mmap_obj.find(b'\n', header_pos) + 1
            
            # Find end of sequence (before next header or EOF)
            if i + 1 < len(headers):
                seq_end = headers[i + 1][1]
            else:
                seq_end = len(self.mmap_obj)
            
            self.sequence_indices[name] = (seq_start, seq_end)
    
    def get_sequence_names(self) -> List[str]:
        """Get list of all sequence names in the FASTA file."""
        return list(self.sequence_indices.keys())
    
    def get_sequence(self, name: str) -> str:
        """
        Get complete sequence by name.
        
        Args:
            name: Sequence name/identifier
            
        Returns:
            Complete DNA sequence as string
        """
        if name not in self.sequence_indices:
            raise KeyError(f"Sequence '{name}' not found in FASTA file")
        
        start, end = self.sequence_indices[name]
        raw_seq = self.mmap_obj[start:end]
        
        # Remove newlines and convert to string
        clean_seq = raw_seq.replace(b'\n', b'').replace(b'\r', b'').decode('ascii')
        return clean_seq.upper()
    
    def get_sequence_slice(self, name: str, start: int, end: int) -> str:
        """
        Get a slice of sequence efficiently.
        
        Args:
            name: Sequence name/identifier
            start: Start position (0-based)
            end: End position (0-based, exclusive)
            
        Returns:
            Sequence slice as string
        """
        full_sequence = self.get_sequence(name)
        return full_sequence[start:end]
    
    def get_sequence_length(self, name: str) -> int:
        """
        Get length of sequence without loading it completely.
        
        Args:
            name: Sequence name/identifier
            
        Returns:
            Sequence length in base pairs
        """
        if name not in self.sequence_indices:
            raise KeyError(f"Sequence '{name}' not found in FASTA file")
        
        start, end = self.sequence_indices[name]
        raw_seq = self.mmap_obj[start:end]
        
        # Count non-newline characters
        newline_count = raw_seq.count(b'\n') + raw_seq.count(b'\r')
        return (end - start) - newline_count
    
    def iterate_sequences(self) -> Iterator[Tuple[str, str]]:
        """
        Iterate over all sequences in the file.
        
        Yields:
            Tuples of (sequence_name, sequence_string)
        """
        for name in self.sequence_indices:
            yield name, self.get_sequence(name)
    
    def get_all_sequences(self) -> Dict[str, str]:
        """
        Load all sequences into a dictionary.
        
        Returns:
            Dictionary mapping sequence names to sequences
        """
        return {name: self.get_sequence(name) for name in self.sequence_indices}
    
    def create_windows(self, name: str, window_size: int = 100000, 
                      overlap: int = 1000) -> Iterator[Tuple[int, int, str]]:
        """
        Create sliding windows over a sequence.
        
        Args:
            name: Sequence name
            window_size: Size of each window
            overlap: Overlap between windows
            
        Yields:
            Tuples of (start_pos, end_pos, window_sequence)
        """
        seq_length = self.get_sequence_length(name)
        
        start = 0
        while start < seq_length:
            end = min(start + window_size, seq_length)
            window_seq = self.get_sequence_slice(name, start, end)
            
            yield start, end, window_seq
            
            if end >= seq_length:
                break
                
            start = end - overlap
    
    def close(self):
        """Close memory-mapped file and file handle."""
        if self.mmap_obj:
            self.mmap_obj.close()
        if self.file_handle:
            self.file_handle.close()
    
    def __enter__(self):
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()

def parse_fasta_string(fasta_content: str) -> Dict[str, str]:
    """
    Parse FASTA content from string.
    
    Args:
        fasta_content: FASTA formatted string
        
    Returns:
        Dictionary mapping sequence names to sequences
    """
    sequences = {}
    current_name = None
    current_seq = []
    
    for line in fasta_content.strip().split('\n'):
        line = line.strip()
        
        if line.startswith('>'):
            # Save previous sequence
            if current_name is not None:
                sequences[current_name] = ''.join(current_seq).upper()
            
            # Start new sequence
            current_name = line[1:].split()[0]  # Take first word after >
            current_seq = []
            
        elif current_name is not None:
            # Add to current sequence
            current_seq.append(line.replace(' ', '').replace('U', 'T'))
    
    # Save last sequence
    if current_name is not None:
        sequences[current_name] = ''.join(current_seq).upper()
    
    return sequences

def write_fasta(sequences: Dict[str, str], output_path: Union[str, Path], 
               line_width: int = 80):
    """
    Write sequences to FASTA file.
    
    Args:
        sequences: Dictionary mapping names to sequences
        output_path: Output file path
        line_width: Maximum line width for sequence lines
    """
    with open(output_path, 'w') as f:
        for name, sequence in sequences.items():
            f.write(f'>{name}\n')
            
            # Write sequence with line breaks
            for i in range(0, len(sequence), line_width):
                f.write(sequence[i:i + line_width] + '\n')

def extract_sequence_regions(fasta_path: Union[str, Path], 
                           regions: List[Tuple[str, int, int]]) -> Dict[str, str]:
    """
    Extract specific regions from FASTA file efficiently.
    
    Args:
        fasta_path: Path to FASTA file
        regions: List of (sequence_name, start, end) tuples
        
    Returns:
        Dictionary mapping region identifiers to extracted sequences
    """
    extracted = {}
    
    with FastaReader(fasta_path) as reader:
        for seq_name, start, end in regions:
            try:
                region_seq = reader.get_sequence_slice(seq_name, start, end)
                region_id = f"{seq_name}:{start}-{end}"
                extracted[region_id] = region_seq
            except KeyError:
                print(f"Warning: Sequence '{seq_name}' not found")
    
    return extracted

__all__ = [
    'FastaReader',
    'parse_fasta_string',
    'write_fasta',
    'extract_sequence_regions'
]