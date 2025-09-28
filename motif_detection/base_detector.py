"""
Base Detector Class for Modular Motif Detection
==============================================

Abstract base class defining common interface for all motif detectors.
"""

import re
from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple, Optional


class BaseMotifDetector(ABC):
    """Abstract base class for all motif detectors"""
    
    def __init__(self):
        self.patterns = self.get_patterns()
        self.compiled_patterns = self._compile_patterns()
    
    @abstractmethod
    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return patterns specific to this motif class"""
        pass
    
    @abstractmethod  
    def get_motif_class_name(self) -> str:
        """Return the motif class name"""
        pass
    
    @abstractmethod
    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        """Calculate motif-specific score"""
        pass
    
    def _compile_patterns(self) -> Dict[str, List[Tuple]]:
        """Compile all regex patterns for performance"""
        compiled_patterns = {}
        
        for pattern_group, patterns in self.patterns.items():
            compiled_group = []
            for pattern_info in patterns:
                pattern, pattern_id, name, subclass = pattern_info[:4]
                try:
                    compiled_pattern = re.compile(pattern, re.IGNORECASE)
                    compiled_group.append((compiled_pattern, pattern_id, name, subclass, pattern_info))
                except re.error as e:
                    print(f"Warning: Invalid pattern {pattern}: {e}")
                    continue
            compiled_patterns[pattern_group] = compiled_group
        
        return compiled_patterns
    
    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Main detection method"""
        sequence = sequence.upper().strip()
        motifs = []
        
        for pattern_group, compiled_patterns in self.compiled_patterns.items():
            for compiled_pattern, pattern_id, name, subclass, full_pattern_info in compiled_patterns:
                for match in compiled_pattern.finditer(sequence):
                    start, end = match.span()
                    motif_seq = sequence[start:end]
                    
                    # Calculate motif-specific score
                    score = self.calculate_score(motif_seq, full_pattern_info)
                    
                    # Apply quality thresholds
                    if self.passes_quality_threshold(motif_seq, score, full_pattern_info):
                        motifs.append({
                            'ID': f"{sequence_name}_{pattern_id}_{start+1}",
                            'Sequence_Name': sequence_name,
                            'Class': self.get_motif_class_name(),
                            'Subclass': subclass,
                            'Start': start + 1,  # 1-based coordinates
                            'End': end,
                            'Length': len(motif_seq),
                            'Sequence': motif_seq,
                            'Score': round(score, 3),
                            'Strand': '+',
                            'Method': f'{self.get_motif_class_name()}_detection',
                            'Pattern_ID': pattern_id
                        })
        
        return motifs
    
    def passes_quality_threshold(self, sequence: str, score: float, pattern_info: Tuple) -> bool:
        """Apply quality thresholds - can be overridden by subclasses"""
        # Default threshold from pattern info if available
        if len(pattern_info) > 6:
            min_threshold = pattern_info[6]  # confidence/threshold from pattern
            return score >= min_threshold
        
        # Default minimum score threshold
        return score >= 0.5
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get detector statistics"""
        total_patterns = sum(len(patterns) for patterns in self.patterns.values())
        return {
            'motif_class': self.get_motif_class_name(),
            'total_patterns': total_patterns,
            'pattern_groups': list(self.patterns.keys()),
            'patterns_by_group': {k: len(v) for k, v in self.patterns.items()}
        }