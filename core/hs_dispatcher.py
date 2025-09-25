"""
Hyperscan Dispatcher for One-Pass Motif Detection
================================================

Routes Hyperscan pattern matches to appropriate class-specific modules
for detailed analysis and scoring. Provides high-performance pattern
matching with intelligent dispatch to motif detectors.
"""

import hyperscan
from typing import List, Dict, Any, Callable, Tuple, Optional
from collections import defaultdict
import threading
import time

from .hyperscan_manager import HyperscanManager
from .regex_registry import get_pattern_info, get_all_hyperscan_patterns

class HyperscanDispatcher:
    """
    High-performance dispatcher for routing Hyperscan matches to motif detectors.
    """
    
    def __init__(self):
        """Initialize dispatcher with pattern registry and detector mapping."""
        self.hs_manager = HyperscanManager()
        self.pattern_registry = {}
        self.class_detectors = {}
        self.database = None
        self._lock = threading.Lock()
        
        # Initialize pattern registry
        self._build_pattern_registry()
        
    def _build_pattern_registry(self):
        """Build mapping from pattern IDs to motif classes and detectors."""
        patterns = get_all_hyperscan_patterns()
        
        for pattern, pattern_id in patterns:
            pattern_info = get_pattern_info(pattern_id)
            if pattern_info:
                self.pattern_registry[pattern_id] = {
                    'pattern': pattern,
                    'motif_class': pattern_info.get('motif_class', 'Unknown'),
                    'subclass': pattern_info.get('subclass', ''),
                    'scoring_function': pattern_info.get('scoring_function', None),
                    'score_method': pattern_info.get('score_method', 'default')
                }
    
    def register_detector(self, motif_class: str, detector_function: Callable):
        """
        Register a detector function for a specific motif class.
        
        Args:
            motif_class: Name of the motif class
            detector_function: Function to call for detailed analysis
        """
        self.class_detectors[motif_class] = detector_function
    
    def compile_database(self, force_recompile: bool = False) -> bool:
        """
        Compile Hyperscan database from all registered patterns.
        
        Args:
            force_recompile: Force recompilation even if cached
            
        Returns:
            True if compilation successful
        """
        with self._lock:
            try:
                patterns = get_all_hyperscan_patterns()
                
                if not patterns:
                    raise ValueError("No patterns available for compilation")
                
                # Use HyperscanManager for compilation
                self.database = self.hs_manager.compile_database(
                    patterns, 
                    flags=hyperscan.HS_FLAG_CASELESS | hyperscan.HS_FLAG_DOTALL
                )
                
                return self.database is not None
                
            except Exception as e:
                print(f"Error compiling Hyperscan database: {e}")
                return False
    
    def scan_sequence(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """
        Perform one-pass scan with intelligent dispatch to detectors.
        
        Args:
            sequence: DNA sequence to scan
            sequence_name: Name identifier for the sequence
            
        Returns:
            List of detected motifs with detailed analysis
        """
        if not self.database:
            if not self.compile_database():
                return []
        
        # Collect raw matches
        raw_matches = []
        
        def match_handler(pattern_id: int, start: int, end: int, flags: int, context=None):
            """Handle Hyperscan match events."""
            raw_matches.append({
                'pattern_id': pattern_id,
                'start': start,
                'end': end,
                'length': end - start,
                'sequence': sequence[start:end]
            })
        
        try:
            # Perform Hyperscan scan
            hyperscan.hs_scan(
                self.database,
                sequence.encode(),
                match_handler,
                None
            )
            
        except Exception as e:
            print(f"Error during Hyperscan scan: {e}")
            return []
        
        # Dispatch matches to appropriate detectors
        return self._dispatch_matches(raw_matches, sequence, sequence_name)
    
    def _dispatch_matches(self, raw_matches: List[Dict[str, Any]], 
                         sequence: str, sequence_name: str) -> List[Dict[str, Any]]:
        """
        Dispatch raw Hyperscan matches to class-specific detectors.
        
        Args:
            raw_matches: Raw pattern matches from Hyperscan
            sequence: Full sequence context
            sequence_name: Sequence identifier
            
        Returns:
            List of fully analyzed motifs
        """
        # Group matches by motif class
        class_matches = defaultdict(list)
        
        for match in raw_matches:
            pattern_id = match['pattern_id']
            pattern_info = self.pattern_registry.get(pattern_id, {})
            motif_class = pattern_info.get('motif_class', 'Unknown')
            
            # Add pattern metadata to match
            enhanced_match = match.copy()
            enhanced_match.update(pattern_info)
            
            class_matches[motif_class].append(enhanced_match)
        
        # Dispatch to class-specific detectors
        all_motifs = []
        
        for motif_class, matches in class_matches.items():
            if motif_class in self.class_detectors:
                try:
                    # Call class-specific detector
                    detector = self.class_detectors[motif_class]
                    detailed_motifs = detector(matches, sequence, sequence_name)
                    
                    if detailed_motifs:
                        all_motifs.extend(detailed_motifs)
                        
                except Exception as e:
                    print(f"Error in {motif_class} detector: {e}")
                    # Fallback to basic motif creation
                    fallback_motifs = self._create_fallback_motifs(matches, motif_class)
                    all_motifs.extend(fallback_motifs)
            else:
                # No specific detector - create basic motifs
                fallback_motifs = self._create_fallback_motifs(matches, motif_class)
                all_motifs.extend(fallback_motifs)
        
        return all_motifs
    
    def _create_fallback_motifs(self, matches: List[Dict[str, Any]], 
                               motif_class: str) -> List[Dict[str, Any]]:
        """
        Create basic motif records when no specific detector is available.
        
        Args:
            matches: Raw matches for the class
            motif_class: Motif class name
            
        Returns:
            List of basic motif dictionaries
        """
        motifs = []
        
        for match in matches:
            motif = {
                'Class': motif_class,
                'Subclass': match.get('subclass', 'Generic'),
                'Start': match['start'] + 1,  # Convert to 1-based
                'End': match['end'],
                'Length': match['length'],
                'Sequence': match['sequence'],
                'Score': 1.0,  # Default score
                'Strand': '+',
                'Method': 'Hyperscan_Pattern',
                'Pattern_ID': match['pattern_id']
            }
            
            motifs.append(motif)
        
        return motifs
    
    def get_performance_stats(self) -> Dict[str, Any]:
        """
        Get performance statistics for the dispatcher.
        
        Returns:
            Dictionary of performance metrics
        """
        return {
            'total_patterns': len(self.pattern_registry),
            'registered_detectors': len(self.class_detectors),
            'database_compiled': self.database is not None,
            'pattern_classes': list(set(
                info.get('motif_class', 'Unknown') 
                for info in self.pattern_registry.values()
            ))
        }
    
    def batch_scan_sequences(self, sequences: Dict[str, str]) -> Dict[str, List[Dict[str, Any]]]:
        """
        Scan multiple sequences efficiently.
        
        Args:
            sequences: Dictionary mapping sequence names to sequences
            
        Returns:
            Dictionary mapping sequence names to motif lists
        """
        results = {}
        
        for seq_name, sequence in sequences.items():
            start_time = time.time()
            motifs = self.scan_sequence(sequence, seq_name)
            scan_time = time.time() - start_time
            
            results[seq_name] = {
                'motifs': motifs,
                'scan_time_seconds': scan_time,
                'sequence_length': len(sequence),
                'motifs_per_kb': len(motifs) / (len(sequence) / 1000) if sequence else 0
            }
        
        return results

# Global dispatcher instance
_global_dispatcher = None
_dispatcher_lock = threading.Lock()

def get_global_dispatcher() -> HyperscanDispatcher:
    """Get or create the global dispatcher instance."""
    global _global_dispatcher
    
    if _global_dispatcher is None:
        with _dispatcher_lock:
            if _global_dispatcher is None:
                _global_dispatcher = HyperscanDispatcher()
    
    return _global_dispatcher

def register_all_detectors(dispatcher: HyperscanDispatcher = None):
    """
    Register all available detector functions with the dispatcher.
    
    Args:
        dispatcher: Dispatcher instance (uses global if None)
    """
    if dispatcher is None:
        dispatcher = get_global_dispatcher()
    
    try:
        # Import detector functions (these would be implemented in detector modules)
        from ..detectors.class01_curved import detect_curved_dna
        from ..detectors.class02_slipped import detect_slipped_dna  
        from ..detectors.class03_cruciform import detect_cruciform
        from ..detectors.class04_rloop import detect_rloop
        from ..detectors.class05_triplex import detect_triplex
        from ..detectors.class06_g4_family import detect_g4_family
        from ..detectors.class07_imotif import detect_imotif
        from ..detectors.class08_zdna import detect_zdna
        from ..detectors.class09_hybrid import detect_hybrid
        from ..detectors.class10_cluster import detect_cluster
        
        # Register detectors
        dispatcher.register_detector('Curved_DNA', detect_curved_dna)
        dispatcher.register_detector('Slipped_DNA', detect_slipped_dna)
        dispatcher.register_detector('Cruciform', detect_cruciform)
        dispatcher.register_detector('R-Loop', detect_rloop)
        dispatcher.register_detector('Triplex', detect_triplex)
        dispatcher.register_detector('G-Quadruplex', detect_g4_family)
        dispatcher.register_detector('i-Motif', detect_imotif)
        dispatcher.register_detector('Z-DNA', detect_zdna)
        dispatcher.register_detector('Hybrid', detect_hybrid)
        dispatcher.register_detector('Cluster', detect_cluster)
        
    except ImportError as e:
        print(f"Warning: Could not import all detectors: {e}")

__all__ = [
    'HyperscanDispatcher',
    'get_global_dispatcher',
    'register_all_detectors'
]