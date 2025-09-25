"""
Streaming Orchestrator for Large-Scale Genomic Analysis
======================================================

Handles chunked processing of large genomic sequences with memory-efficient
streaming, progress tracking, and result aggregation.
"""

import time
import gc
from typing import Iterator, List, Dict, Any, Optional, Callable, Tuple
from pathlib import Path
import threading
import queue
from concurrent.futures import ThreadPoolExecutor, as_completed

from ..core.windows import SequenceWindower, GenomeWindower, SequenceWindow
from ..core.postprocess import apply_all_postprocessing
from ..io.fasta import FastaReader
from ..io.schemas import AnalysisResults, SequenceInfo, AnalysisConfig, PerformanceMetrics
from .all_motifs import detect_all_motifs

class StreamingOrchestrator:
    """
    High-performance streaming orchestrator for large-scale motif detection.
    """
    
    def __init__(self, config: Optional[AnalysisConfig] = None):
        """
        Initialize streaming orchestrator.
        
        Args:
            config: Analysis configuration (uses defaults if None)
        """
        self.config = config or AnalysisConfig()
        self.windower = SequenceWindower(
            chunk_size=self.config.chunk_size,
            overlap_size=self.config.overlap_size
        )
        self._progress_callbacks = []
        self._cancelled = False
        
    def add_progress_callback(self, callback: Callable[[Dict[str, Any]], None]):
        """Add callback function for progress updates."""
        self._progress_callbacks.append(callback)
    
    def _notify_progress(self, progress_info: Dict[str, Any]):
        """Notify all registered progress callbacks."""
        for callback in self._progress_callbacks:
            try:
                callback(progress_info)
            except Exception as e:
                print(f"Progress callback error: {e}")
    
    def cancel(self):
        """Cancel ongoing analysis."""
        self._cancelled = True
    
    def analyze_sequence_stream(self, sequence: str, sequence_name: str = "sequence") -> AnalysisResults:
        """
        Analyze sequence using streaming approach.
        
        Args:
            sequence: Input DNA sequence
            sequence_name: Sequence identifier
            
        Returns:
            Complete analysis results
        """
        start_time = time.time()
        
        # Create sequence info
        seq_info = SequenceInfo(
            name=sequence_name,
            length=len(sequence),
            gc_content=self._calculate_gc_content(sequence)
        )
        
        # Create windows
        windows = self.windower.create_windows(sequence, sequence_name)
        total_windows = len(windows)
        
        self._notify_progress({
            'stage': 'initialization',
            'total_windows': total_windows,
            'sequence_length': len(sequence)
        })
        
        # Process windows
        all_motifs = []
        processed_windows = 0
        
        for window in windows:
            if self._cancelled:
                break
                
            # Detect motifs in window
            window_motifs = self._process_window(window, sequence)
            all_motifs.extend(window_motifs)
            
            processed_windows += 1
            
            # Progress update
            self._notify_progress({
                'stage': 'processing',
                'windows_processed': processed_windows,
                'total_windows': total_windows,
                'progress_percent': 100 * processed_windows / total_windows,
                'current_motifs': len(all_motifs)
            })
            
            # Memory management
            if processed_windows % 10 == 0:
                gc.collect()
        
        # Post-processing
        self._notify_progress({'stage': 'postprocessing'})
        
        postprocess_config = {
            'min_score_threshold': self.config.min_score_threshold,
            'length_limits': self._get_length_limits(),
            'remove_overlaps': self.config.remove_overlaps,
            'merge_nearby': self.config.merge_nearby,
            'merge_distance': self.config.merge_distance
        }
        
        final_motifs, processing_stats = apply_all_postprocessing(all_motifs, postprocess_config)
        
        # Calculate final statistics
        processing_time = time.time() - start_time
        
        # Update sequence info with results
        seq_info.total_motifs = len(final_motifs)
        seq_info.motif_coverage = self._calculate_coverage(final_motifs, len(sequence))
        
        # Create results
        results = AnalysisResults(
            analysis_id=f"{sequence_name}_{int(time.time())}",
            sequence_info=seq_info,
            config=self.config,
            motifs=final_motifs,
            total_motifs=len(final_motifs),
            class_distribution=processing_stats.get('class_distribution', {}),
            coverage_statistics=processing_stats,
            processing_time_seconds=processing_time
        )
        
        self._notify_progress({
            'stage': 'complete',
            'total_motifs': len(final_motifs),
            'processing_time': processing_time
        })
        
        return results
    
    def analyze_fasta_file_stream(self, fasta_path: Path) -> Iterator[AnalysisResults]:
        """
        Stream analysis of FASTA file with multiple sequences.
        
        Args:
            fasta_path: Path to FASTA file
            
        Yields:
            AnalysisResults for each sequence
        """
        with FastaReader(fasta_path) as reader:
            sequence_names = reader.get_sequence_names()
            total_sequences = len(sequence_names)
            
            self._notify_progress({
                'stage': 'file_analysis_start',
                'total_sequences': total_sequences,
                'file_path': str(fasta_path)
            })
            
            for i, seq_name in enumerate(sequence_names):
                if self._cancelled:
                    break
                
                self._notify_progress({
                    'stage': 'sequence_start',
                    'sequence_index': i + 1,
                    'total_sequences': total_sequences,
                    'sequence_name': seq_name
                })
                
                try:
                    # Load sequence
                    sequence = reader.get_sequence(seq_name)
                    
                    # Analyze sequence
                    results = self.analyze_sequence_stream(sequence, seq_name)
                    
                    yield results
                    
                except Exception as e:
                    print(f"Error analyzing sequence {seq_name}: {e}")
                    continue
            
            self._notify_progress({'stage': 'file_analysis_complete'})
    
    def analyze_genomic_windows(self, fasta_path: Path, 
                              window_size: int = 1000000,
                              overlap: int = 10000) -> Iterator[Tuple[str, int, int, List[Dict[str, Any]]]]:
        """
        Stream analysis of genomic windows for very large sequences.
        
        Args:
            fasta_path: Path to FASTA file
            window_size: Size of genomic windows (default 1Mb)
            overlap: Overlap between windows
            
        Yields:
            Tuples of (sequence_name, start_pos, end_pos, motifs)
        """
        genome_windower = GenomeWindower(window_size, overlap)
        
        with FastaReader(fasta_path) as reader:
            for seq_name in reader.get_sequence_names():
                if self._cancelled:
                    break
                
                # Create windows for this chromosome/sequence
                for start, end, window_seq in reader.create_windows(
                    seq_name, window_size, overlap):
                    
                    if self._cancelled:
                        break
                    
                    # Detect motifs in window
                    motifs = detect_all_motifs(window_seq, f"{seq_name}:{start}-{end}")
                    
                    # Adjust coordinates to global positions
                    for motif in motifs:
                        motif['Start'] += start
                        motif['End'] += start
                        motif['Chromosome'] = seq_name
                    
                    yield seq_name, start, end, motifs
    
    def _process_window(self, window: SequenceWindow, full_sequence: str) -> List[Dict[str, Any]]:
        """Process a single sequence window."""
        # Extract window sequence with proper coordinates
        window_seq = window.sequence
        
        # Detect motifs
        motifs = detect_all_motifs(window_seq, f"window_{window.chunk_id}")
        
        # Adjust coordinates to global positions
        for motif in motifs:
            motif['Start'] += window.start_pos
            motif['End'] += window.start_pos
            motif['Window_ID'] = window.chunk_id
        
        return motifs
    
    def _calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content percentage."""
        if not sequence:
            return 0.0
        
        gc_count = sequence.upper().count('G') + sequence.upper().count('C')
        return 100 * gc_count / len(sequence)
    
    def _calculate_coverage(self, motifs: List[Dict[str, Any]], sequence_length: int) -> float:
        """Calculate percentage of sequence covered by motifs."""
        if not motifs or sequence_length == 0:
            return 0.0
        
        covered_bases = sum(motif.get('Length', 0) for motif in motifs)
        return 100 * covered_bases / sequence_length
    
    def _get_length_limits(self) -> Dict[str, Tuple[int, int]]:
        """Get length limits for motif classes."""
        # This could be loaded from config files
        return {
            'Curved_DNA': (15, 200),
            'Slipped_DNA': (11, 300),
            'Cruciform': (23, 100),
            'R-Loop': (140, 400),
            'Triplex': (30, 400),
            'G-Quadruplex': (13, 100),
            'i-Motif': (23, 100),
            'Z-DNA': (50, 200),
            'Hybrid': (10, 500),
            'Cluster': (100, 10000)
        }

class BatchStreamingOrchestrator(StreamingOrchestrator):
    """
    Orchestrator for batch processing multiple files.
    """
    
    def __init__(self, config: Optional[AnalysisConfig] = None, max_workers: int = None):
        """
        Initialize batch orchestrator.
        
        Args:
            config: Analysis configuration
            max_workers: Maximum number of parallel workers
        """
        super().__init__(config)
        self.max_workers = max_workers or self.config.max_workers
    
    def analyze_multiple_files(self, fasta_files: List[Path]) -> Iterator[Tuple[Path, AnalysisResults]]:
        """
        Analyze multiple FASTA files in parallel.
        
        Args:
            fasta_files: List of FASTA file paths
            
        Yields:
            Tuples of (file_path, analysis_results)
        """
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all files for processing
            future_to_file = {
                executor.submit(self._analyze_single_file, fasta_file): fasta_file
                for fasta_file in fasta_files
            }
            
            # Yield results as they complete
            for future in as_completed(future_to_file):
                if self._cancelled:
                    break
                
                fasta_file = future_to_file[future]
                try:
                    results = future.result()
                    for result in results:
                        yield fasta_file, result
                except Exception as e:
                    print(f"Error processing {fasta_file}: {e}")
    
    def _analyze_single_file(self, fasta_path: Path) -> List[AnalysisResults]:
        """Analyze a single FASTA file and return all results."""
        results = []
        for result in self.analyze_fasta_file_stream(fasta_path):
            results.append(result)
        return results

__all__ = [
    'StreamingOrchestrator',
    'BatchStreamingOrchestrator'
]