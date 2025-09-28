#!/usr/bin/env python3
"""
QmRLFS-finder Integration for NBDScanner
========================================

R-loop Forming Sequence (RLFS) detection based on QmRLFS-finder v1.5 algorithm
by Piroon Jenjaroenpun and Thidathip Wongsurawat

This module integrates the QmRLFS logic for detecting R-loop forming regions
in DNA sequences, including RIZ (R-loop Initiating Zone) and REZ (R-loop 
Extending Zone) detection.

References:
- QmRLFS-finder v1.5 (2016) 
- Piroon Jenjaroenpun, Thidathip Wongsurawat
- Email: piroonj@bii.a-star.edu.sg
"""

import re
import math
from typing import List, Dict, Any, Tuple, Optional, Union
from collections import defaultdict

class QmRLFSDetector:
    """
    QmRLFS-based R-loop formation site detector
    
    Implements the QmRLFS-finder algorithm with RIZ and REZ detection
    """
    
    def __init__(self, 
                 models: str = "m1,m2",
                 min_perc_g_riz: float = 50.0,
                 num_linker: int = 50,
                 window_step: int = 100,
                 max_length_rez: int = 2000,
                 min_perc_g_rez: float = 40.0,
                 quick_mode: bool = False):
        """
        Initialize QmRLFS detector with parameters
        
        Args:
            models: Comma-separated model names (m1, m2)
            min_perc_g_riz: Minimum G percentage for RIZ
            num_linker: Number of linker positions to check
            window_step: Window step size for REZ search
            max_length_rez: Maximum REZ length
            min_perc_g_rez: Minimum G percentage for REZ
            quick_mode: Enable quick mode (shorter search)
        """
        self.models = {}
        self.min_perc_g_riz = min_perc_g_riz
        self.num_linker = num_linker
        self.window_step = window_step
        self.max_length_rez = max_length_rez
        self.min_perc_g_rez = min_perc_g_rez
        self.quick_mode = quick_mode
        
        # Initialize regex models
        input_models = set(models.split(","))
        if "m1" in input_models:
            self.models['m1'] = re.compile(r"G{3,}[ATCGU]{1,10}?G{3,}(?:[ATCGU]{1,10}?G{3,}){1,}?")
        if "m2" in input_models:
            self.models['m2'] = re.compile(r"G{4,}(?:[ATCGU]{1,10}?G{4,}){1,}?")
    
    def percent_g(self, seq: str) -> float:
        """Calculate percentage of G bases in sequence"""
        if not seq:
            return 0.0
        return (seq.count("G") / len(seq)) * 100.0
    
    def riz_search(self, seq: str, model: str) -> List[Dict[str, Any]]:
        """
        Search for RIZ (R-loop Initiating Zone) regions
        
        Args:
            seq: DNA sequence
            model: Model name (m1 or m2)
            
        Returns:
            List of RIZ region dictionaries
        """
        result_list = []
        
        if model not in self.models:
            return result_list
        
        for result in self.models[model].finditer(seq):
            start = result.start()
            end = result.end()
            riz_seq = result.group(0)
            perc_g = self.percent_g(riz_seq)
            
            if perc_g >= self.min_perc_g_riz:
                dict_result = {
                    "start": start,
                    "end": end,
                    "length": end - start,
                    "G": riz_seq.count("G"),
                    "G3s": riz_seq.count("GGG"),
                    "G4s": riz_seq.count("GGGG"),
                    "perc_g": perc_g,
                    "model": model,
                    "seq": riz_seq
                }
                result_list.append(dict_result)
        
        return result_list
    
    def rez_search(self, riz_end: int, seq: str, seqlength: int) -> Dict[str, Any]:
        """
        Search for REZ (R-loop Extending Zone) regions
        
        Args:
            riz_end: End position of RIZ
            seq: Sequence after RIZ end
            seqlength: Total sequence length
            
        Returns:
            Dictionary with REZ information or empty dict
        """
        if not seq or len(seq) < 20:
            return {}
            
        max_length = 0
        best_start = 0
        best_end = 0
        
        # Simplified REZ search - look for G-rich or C-rich extending regions
        search_limit = min(len(seq), self.num_linker) if not self.quick_mode else min(len(seq), self.num_linker // 2)
        
        for start_pos in range(search_limit):
            # Start with minimum required length
            min_length = max(20, self.window_step // 5)  # Minimum 20 bp
            max_search_len = min(len(seq) - start_pos, self.max_length_rez)
            
            for end_pos in range(start_pos + min_length, start_pos + max_search_len + 1, 5):
                if end_pos > len(seq):
                    break
                    
                candidate_seq = seq[start_pos:end_pos]
                if not candidate_seq:
                    continue
                
                # Calculate G% and C% for REZ evaluation
                perc_g = self.percent_g(candidate_seq)
                perc_c = (candidate_seq.count("C") / len(candidate_seq)) * 100.0
                gc_content = perc_g + perc_c
                
                # REZ criteria: either G-rich, C-rich, or high GC content
                is_valid_rez = (perc_g >= self.min_perc_g_rez or 
                               perc_c >= self.min_perc_g_rez or
                               gc_content >= self.min_perc_g_rez * 1.2)
                
                # Additional checks
                within_sequence_bounds = (riz_end + end_pos <= seqlength)
                longer_than_current = len(candidate_seq) > max_length
                
                if is_valid_rez and within_sequence_bounds and longer_than_current:
                    max_length = len(candidate_seq)
                    best_start = start_pos
                    best_end = end_pos
            
            # In quick mode, stop after finding first valid REZ
            if self.quick_mode and max_length > 0:
                break
        
        # Create REZ region dictionary if found
        if max_length > 0:
            return self.get_rez_element(riz_end, best_start, best_end, seq)
        
        return {}
    
    def get_rez_element(self, riz_end: int, rez_start: int, rez_end: int, seq: str) -> Dict[str, Any]:
        """Create REZ element dictionary"""
        rez_seq = seq[rez_start:rez_end]
        return {
            "start": riz_end + rez_start,
            "end": riz_end + rez_end,
            "length": rez_end - rez_start,
            "G": rez_seq.count("G"),
            "G3s": rez_seq.count("GGG"),
            "G4s": rez_seq.count("GGGG"),
            "perc_g": self.percent_g(rez_seq),
            "seq": rez_seq
        }
    
    def search_rlfs(self, sequence: str, strand: str = "+", sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """
        Search for R-loop forming sequences (RLFS) in DNA sequence
        
        Args:
            sequence: DNA sequence to analyze
            strand: Strand orientation (+ or -)
            sequence_name: Name of the sequence
            
        Returns:
            List of RLFS dictionaries with RIZ and REZ information
        """
        results = []
        seq = sequence.upper()
        seqlength = len(seq)
        
        for model in self.models:
            # Find RIZ regions
            riz_list = self.riz_search(seq, model)
            
            if riz_list:
                for riz in riz_list:
                    # Find REZ for each RIZ
                    rez = self.rez_search(riz["end"], seq[riz["end"]:], seqlength)
                    
                    if rez:
                        # Calculate additional metrics
                        riz_rez_sorted = sorted([riz['start'], riz['end'], rez['start'], rez['end']])
                        linker_length = riz_rez_sorted[2] - riz_rez_sorted[1]
                        linker_seq = seq[riz_rez_sorted[1]:riz_rez_sorted[2]]
                        
                        # Adjust coordinates for negative strand
                        if strand == "-":
                            RIZ_start = seqlength - riz["end"]
                            RIZ_end = seqlength - riz["start"]
                            REZ_start = seqlength - rez["end"]
                            REZ_end = seqlength - rez["start"]
                            riz["start"] = RIZ_start
                            riz["end"] = RIZ_end
                            rez["start"] = REZ_start
                            rez["end"] = REZ_end
                        
                        # Create RLFS result
                        rlfs_result = {
                            "sequence_name": sequence_name,
                            "model": model,
                            "strand": strand,
                            "location": f"{sequence_name}:{min(riz['start'], rez['start'])+1}-{max(riz['end'], rez['end'])}",
                            
                            # RIZ information
                            "start_RIZ": riz["start"],
                            "end_RIZ": riz["end"],
                            "length_RIZ": riz["length"],
                            "G_RIZ": riz["G"],
                            "3Gs_RIZ": riz["G3s"],
                            "4Gs_RIZ": riz["G4s"],
                            "perc_G_RIZ": riz["perc_g"],
                            "sequence_RIZ": riz["seq"],
                            
                            # Linker information
                            "linker_length": linker_length,
                            "linker_seq": linker_seq,
                            
                            # REZ information
                            "start_REZ": rez["start"],
                            "end_REZ": rez["end"],
                            "length_REZ": rez["length"],
                            "G_REZ": rez["G"],
                            "3Gs_REZ": rez["G3s"],
                            "4Gs_REZ": rez["G4s"],
                            "perc_G_REZ": rez["perc_g"],
                            "sequence_REZ": rez["seq"],
                            
                            # Overall metrics
                            "total_length": max(riz['end'], rez['end']) - min(riz['start'], rez['start']),
                            "qmrlfs_score": self.calculate_qmrlfs_score(riz, rez, linker_length)
                        }
                        
                        results.append(rlfs_result)
        
        return results
    
    def calculate_qmrlfs_score(self, riz: Dict[str, Any], rez: Dict[str, Any], linker_length: int) -> float:
        """
        Calculate QmRLFS score based on RIZ, REZ, and linker characteristics
        
        Args:
            riz: RIZ region dictionary
            rez: REZ region dictionary  
            linker_length: Length of linker region
            
        Returns:
            QmRLFS score (0.0-1.0)
        """
        # RIZ scoring (G content and G-tract density)
        riz_score = (riz["perc_g"] / 100.0) * 0.4
        riz_g_density = (riz["G3s"] + riz["G4s"] * 1.5) / max(riz["length"], 1)
        riz_score += min(riz_g_density, 1.0) * 0.3
        
        # REZ scoring (G content and length)
        rez_score = (rez["perc_g"] / 100.0) * 0.2
        rez_length_score = min(rez["length"] / 500.0, 1.0) * 0.1
        
        # Linker penalty (shorter linkers are better)
        linker_penalty = max(0, (linker_length - 50) / 200.0)
        
        total_score = riz_score + rez_score + rez_length_score - linker_penalty
        return max(0.0, min(1.0, total_score))
    
    def analyze_sequence(self, sequence: str, sequence_name: str = "sequence", 
                        analyze_both_strands: bool = True) -> List[Dict[str, Any]]:
        """
        Analyze sequence for R-loop forming sites on both strands
        
        Args:
            sequence: DNA sequence to analyze
            sequence_name: Name of the sequence
            analyze_both_strands: Whether to analyze both strands
            
        Returns:
            List of all RLFS found
        """
        all_results = []
        
        # Positive strand
        positive_results = self.search_rlfs(sequence, "+", sequence_name)
        all_results.extend(positive_results)
        
        # Negative strand (reverse complement)
        if analyze_both_strands:
            # Simple reverse complement
            complement_map = str.maketrans("ATGC", "TACG")
            reverse_complement = sequence.translate(complement_map)[::-1]
            
            negative_results = self.search_rlfs(reverse_complement, "-", sequence_name)
            all_results.extend(negative_results)
        
        return all_results
    
    @staticmethod
    def qmrlfs_score_function(sequence: str) -> float:
        """
        Simplified QmRLFS scoring function for integration with NBDScanner
        
        Args:
            sequence: DNA sequence
            
        Returns:
            QmRLFS-based score (0.0-1.0)
        """
        if len(sequence) < 20:
            return 0.0
        
        detector = QmRLFSDetector(quick_mode=True)
        results = detector.analyze_sequence(sequence, analyze_both_strands=False)
        
        if not results:
            return 0.0
        
        # Return the highest scoring RLFS
        max_score = max(result["qmrlfs_score"] for result in results)
        return max_score


# Convenience functions for NBDScanner integration
def detect_r_loop_qmrlfs(sequence: str, sequence_name: str = "sequence", 
                        models: str = "m1,m2", quick_mode: bool = False) -> List[Dict[str, Any]]:
    """
    Detect R-loop forming sites using QmRLFS algorithm
    
    Args:
        sequence: DNA sequence to analyze
        sequence_name: Name of the sequence
        models: Models to use (m1, m2, or m1,m2)
        quick_mode: Enable quick mode for faster analysis
        
    Returns:
        List of R-loop forming sites
    """
    detector = QmRLFSDetector(models=models, quick_mode=quick_mode)
    return detector.analyze_sequence(sequence, sequence_name)


def qmrlfs_to_nbdscanner_format(qmrlfs_results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Convert QmRLFS results to NBDScanner motif format
    
    Args:
        qmrlfs_results: Results from QmRLFS detector
        
    Returns:
        List of motifs in NBDScanner format
    """
    motifs = []
    
    for result in qmrlfs_results:
        motif = {
            'ID': f"{result['sequence_name']}_rlfs_{len(motifs)+1}",
            'Sequence_Name': result['sequence_name'],
            'Class': 'R-loop',
            'Subclass': f'QmRLFS-{result["model"]}',
            'Start': min(result['start_RIZ'], result['start_REZ']),
            'End': max(result['end_RIZ'], result['end_REZ']),
            'Length': result['total_length'],
            'Strand': result['strand'],
            'Raw_Score': result['qmrlfs_score'],
            'Normalized_Score': result['qmrlfs_score'],
            'Sequence': result['sequence_RIZ'] + result['linker_seq'] + result['sequence_REZ'],
            'Method': 'QmRLFS_Pure_Python',
            'Details': {
                'model': result['model'],
                'RIZ_start': result['start_RIZ'],
                'RIZ_end': result['end_RIZ'],
                'RIZ_G_percent': result['perc_G_RIZ'],
                'REZ_start': result['start_REZ'],
                'REZ_end': result['end_REZ'],
                'REZ_G_percent': result['perc_G_REZ'],
                'linker_length': result['linker_length']
            }
        }
        motifs.append(motif)
    
    return motifs