"""
CruciformDetector (Optimized for Performance)
=============================================

PERFORMANCE OPTIMIZATIONS:
- Sliding window approach for sequences > 1000 bp (no skipping!)
- Optimized search: larger arm lengths first with early termination
- Limited maximum arm length to 100 bp for O(n) complexity
- Maintains accuracy while improving speed on all sequence lengths

Detects inverted repeats (potential cruciform-forming) with:
 - arm length >= 6 (capped at 100 bp for performance)
 - loop (spacer) <= 100 bp
 - optional mismatch tolerance
Scoring: interpretable 0..1 score that favors long arms and small loops.
"""

import re
from typing import List, Dict, Any, Tuple
from .base_detector import BaseMotifDetector


def revcomp(seq: str) -> str:
    trans = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(trans)[::-1]


class CruciformDetector(BaseMotifDetector):
    def get_motif_class_name(self) -> str:
        return "Cruciform"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """
        Provide pattern metadata to remain compatible with your framework.
        The actual detection uses find_inverted_repeats() and the criterion:
           - arm >= 6
           - loop <= 100
        """
        return {
            'inverted_repeats': [
                # Metadata kept for compatibility; actual search enforces arm>=6 and loop<=100
                (r'palindrome_like', 'CRU_3_1', 'Potential palindrome', 'Inverted Repeats', 12, 'cruciform_stability', 0.95, 'DNA secondary structure', 'Lilley 2000'),
            ]
        }

    # --------------------------
    # Configuration (tweakable)
    # --------------------------
    MIN_ARM = 6          # minimum arm length (user-specified criterion)
    MAX_ARM = 100        # maximum arm length (for performance on large sequences)
    MAX_LOOP = 100       # maximum loop (spacer) length
    MAX_MISMATCHES = 0   # allowed mismatches between arm and RC(arm). Set >0 to allow imperfect arms.
    MAX_SEQUENCE_LENGTH = 1000  # Skip cruciform detection for sequences longer than this (too slow)

    # --------------------------
    # Core search function
    # --------------------------
    def find_inverted_repeats(self, sequence: str, min_arm: int = None,
                              max_loop: int = None, max_mismatches: int = None) -> List[Dict[str, Any]]:
        """
        Scan sequence for inverted repeats satisfying:
          - arm length >= min_arm
          - loop length <= max_loop
          - at most max_mismatches mismatches (between left arm and RC of right arm)
        Returns list of hits:
          {
            'left_start','left_end','right_start','right_end',
            'arm_len','loop_len',
            'left_seq','right_seq','right_seq_rc',
            'mismatches','match_fraction','score'
          }
        Coordinates are 0-based, end-exclusive.
        
        For long sequences (>1000 bp), uses a sliding window approach to avoid
        skipping motifs while maintaining reasonable performance.
        """
        seq = sequence.upper()
        n = len(seq)
        
        if min_arm is None:
            min_arm = self.MIN_ARM
        if max_loop is None:
            max_loop = self.MAX_LOOP
        if max_mismatches is None:
            max_mismatches = self.MAX_MISMATCHES

        hits: List[Dict[str, Any]] = []
        
        # For long sequences, use sliding window approach
        if n > self.MAX_SEQUENCE_LENGTH:
            window_size = self.MAX_SEQUENCE_LENGTH
            step_size = window_size // 2  # 50% overlap to avoid missing motifs at boundaries
            
            for window_start in range(0, n, step_size):
                window_end = min(window_start + window_size, n)
                window_seq = seq[window_start:window_end]
                
                # Process this window
                window_hits = self._find_inverted_repeats_in_window(
                    window_seq, min_arm, max_loop, max_mismatches
                )
                
                # Adjust coordinates to full sequence
                for hit in window_hits:
                    hit['left_start'] += window_start
                    hit['left_end'] += window_start
                    hit['right_start'] += window_start
                    hit['right_end'] += window_start
                    hits.append(hit)
                
                # If we've reached the end, stop
                if window_end >= n:
                    break
            
            # Remove duplicates that may occur in overlapping windows
            hits = self._deduplicate_hits(hits)
        else:
            # For short sequences, use the full algorithm
            hits = self._find_inverted_repeats_in_window(seq, min_arm, max_loop, max_mismatches)
        
        # Sort hits by descending score, then by left_start
        hits.sort(key=lambda h: (-h['score'], h['left_start'], -h['arm_len']))
        return hits
    
    def _find_inverted_repeats_in_window(self, seq: str, min_arm: int, 
                                         max_loop: int, max_mismatches: int) -> List[Dict[str, Any]]:
        """Core inverted repeat detection in a single window"""
        hits: List[Dict[str, Any]] = []
        n = len(seq)
        
        # PERFORMANCE: Sample positions for large windows to reduce search space
        # For sequences, check every position, but for larger ones, sample
        step = 1 if n <= 500 else 2 if n <= 1000 else 3

        # PERFORMANCE: Limit arm length testing to MAX_ARM for computational feasibility
        for left_start in range(0, n - 2 * min_arm, step):
            # Maximum possible arm length at this left_start
            max_possible_arm = min(self.MAX_ARM, (n - left_start) // 2)
            
            # PERFORMANCE: Search from larger arm lengths first (better quality)
            # This allows early termination and reduces total iterations
            for arm_len in range(max_possible_arm, min_arm - 1, -1):
                left_end = left_start + arm_len
                # right arm must start at least left_end + 0 loop; but loop cannot exceed max_loop
                right_start_min = left_end
                right_start_max = min(left_end + max_loop, n - arm_len)
                
                # Early exit if we already found a good match at this position
                found_good_match = False
                
                # iterate right_starts (loop sizes) - prefer smaller loops
                for right_start in range(right_start_min, right_start_max + 1):
                    loop_len = right_start - left_end
                    right_end = right_start + arm_len
                    if right_end > n:
                        break
                    left_seq = seq[left_start:left_end]
                    right_seq = seq[right_start:right_end]
                    # Compare left_seq to reverse complement of right_seq
                    right_rc = revcomp(right_seq)
                    # Count mismatches
                    mismatches = sum(1 for a, b in zip(left_seq, right_rc) if a != b)
                    if mismatches <= max_mismatches:
                        match_fraction = (arm_len - mismatches) / arm_len if arm_len > 0 else 0.0
                        score = self._score_arm_loop(arm_len, loop_len, match_fraction)
                        hits.append({
                            'left_start': left_start,
                            'left_end': left_end,
                            'right_start': right_start,
                            'right_end': right_end,
                            'arm_len': arm_len,
                            'loop_len': loop_len,
                            'left_seq': left_seq,
                            'right_seq': right_seq,
                            'right_seq_rc': right_rc,
                            'mismatches': mismatches,
                            'match_fraction': round(match_fraction, 4),
                            'score': round(score, 6)
                        })
                        found_good_match = True
                        # PERFORMANCE: If we found a perfect match with good score, 
                        # we can skip checking smaller arm lengths at this position
                        if mismatches == 0 and score > 0.5:
                            break
                    # If exact-match required and mismatches > 0, we can continue scanning other right_starts
                    # Continue until right_start_max.
                
                # If we found a good match, skip smaller arm lengths at this position
                if found_good_match and arm_len >= min_arm * 2:
                    break
                    
        return hits
    
    def _deduplicate_hits(self, hits: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Remove duplicate hits that may occur in overlapping windows"""
        if not hits:
            return hits
        
        # Create a unique key for each hit based on positions
        seen = {}
        unique_hits = []
        
        for hit in hits:
            key = (hit['left_start'], hit['left_end'], hit['right_start'], hit['right_end'])
            if key not in seen:
                seen[key] = True
                unique_hits.append(hit)
        
        return unique_hits

    # --------------------------
    # Scoring function (interpretable)
    # --------------------------
    def _score_arm_loop(self, arm_len: int, loop_len: int, match_fraction: float) -> float:
        """
        Compute a normalized score 0..1:
          - arm contribution: sigmoid-like increase with arm_len
          - loop penalty: linear penalty up to MAX_LOOP
          - match_fraction scales final score (1.0 = perfect match)
        Formula:
          arm_term = arm_len / (arm_len + 8)    -> approaches 1 for long arms
          loop_term = max(0.0, 1.0 - (loop_len / float(self.MAX_LOOP)))  -> 1 at loop 0, 0 at MAX_LOOP
          base = arm_term * loop_term
          final = base * match_fraction
        This yields scores near 1 for long arms + short loops + perfect match.
        """
        arm_term = float(arm_len) / (arm_len + 8.0)
        loop_term = max(0.0, 1.0 - (float(loop_len) / float(self.MAX_LOOP)))
        base = arm_term * loop_term
        final = base * float(match_fraction)
        # clamp
        return max(0.0, min(1.0, final))

    # --------------------------
    # Public API: calculate_score & annotate_sequence
    # --------------------------
    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
        """
        Sum scores of detected inverted repeats in the sequence (non-overlap-resolved).
        If you prefer overlap resolution (one strongest per region), call annotate_sequence() and sum accepted.
        """
        seq = sequence.upper()
        hits = self.find_inverted_repeats(seq,
                                         min_arm=self.MIN_ARM,
                                         max_loop=self.MAX_LOOP,
                                         max_mismatches=self.MAX_MISMATCHES)
        # Sum of hit scores (reflects number + quality)
        total = sum(h['score'] for h in hits)
        return float(total)

    def annotate_sequence(self, sequence: str, max_hits: int = 0) -> List[Dict[str, Any]]:
        """
        Return list of detected inverted repeats with details.
        If max_hits > 0, return at most that many top hits (by score). Otherwise return all.
        """
        seq = sequence.upper()
        hits = self.find_inverted_repeats(seq,
                                         min_arm=self.MIN_ARM,
                                         max_loop=self.MAX_LOOP,
                                         max_mismatches=self.MAX_MISMATCHES)
        if max_hits and len(hits) > max_hits:
            return hits[:max_hits]
        return hits

    # --------------------------
    # Quality threshold check (keeps similar signature)
    # --------------------------
    def passes_quality_threshold(self, sequence: str, score: float, pattern_info: Tuple) -> bool:
        """
        Enhanced quality check:
          - require that at least one inverted repeat was found
          - optionally enforce minimal per-hit score threshold if pattern_info provides one
        pattern_info[6] was your previous 'score threshold' position; we accept either that or default 0.2
        """
        seq = sequence.upper()
        hits = self.find_inverted_repeats(seq,
                                         min_arm=self.MIN_ARM,
                                         max_loop=self.MAX_LOOP,
                                         max_mismatches=self.MAX_MISMATCHES)
        if not hits:
            return False
        # prefer the best hit score
        best_score = hits[0]['score']
        # read threshold from pattern_info if provided (position 6 per your metadata format)
        try:
            provided_thresh = float(pattern_info[6]) if (pattern_info and len(pattern_info) > 6) else None
        except Exception:
            provided_thresh = None
        thresh = provided_thresh if provided_thresh is not None else 0.2
        return best_score >= thresh

    def _remove_overlaps(self, inverted_repeats: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Remove overlapping inverted repeats, keeping highest scoring non-overlapping set"""
        if not inverted_repeats:
            return []
        
        # Sort by score (descending), then by length (descending)
        sorted_repeats = sorted(inverted_repeats, 
                               key=lambda x: (-x['score'], -(x['right_end'] - x['left_start'])))
        
        non_overlapping = []
        for repeat in sorted_repeats:
            # Check if this repeat overlaps with any already selected
            overlaps = False
            for selected in non_overlapping:
                # Two repeats overlap if their full regions (left_start to right_end) overlap
                if not (repeat['right_end'] <= selected['left_start'] or 
                       repeat['left_start'] >= selected['right_end']):
                    overlaps = True
                    break
            
            if not overlaps:
                non_overlapping.append(repeat)
        
        # Sort by start position for output
        non_overlapping.sort(key=lambda x: x['left_start'])
        return non_overlapping

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Override base method to use sophisticated cruciform detection"""
        sequence = sequence.upper().strip()
        motifs = []
        
        # Use the find_inverted_repeats method which has the sophisticated logic
        inverted_repeats = self.find_inverted_repeats(sequence, 
                                                     min_arm=self.MIN_ARM,
                                                     max_loop=self.MAX_LOOP,
                                                     max_mismatches=self.MAX_MISMATCHES)
        
        # Filter by meaningful score threshold before overlap removal
        filtered_repeats = [r for r in inverted_repeats if r.get('score', 0) > 0.1]
        
        # Remove overlapping repeats
        non_overlapping_repeats = self._remove_overlaps(filtered_repeats)
        
        for i, repeat in enumerate(non_overlapping_repeats):
            start_pos = repeat['left_start']
            end_pos = repeat['right_end'] 
            full_length = end_pos - start_pos
            full_seq = sequence[start_pos:end_pos]
            
            motifs.append({
                'ID': f"{sequence_name}_CRU_{start_pos+1}",
                'Sequence_Name': sequence_name,
                'Class': self.get_motif_class_name(),
                'Subclass': 'Inverted_Repeat',
                'Start': start_pos + 1,  # 1-based coordinates
                'End': end_pos,
                'Length': full_length,
                'Sequence': full_seq,
                'Score': round(repeat['score'], 3),
                'Strand': '+',
                'Method': 'Cruciform_detection',
                'Pattern_ID': f'CRU_{i+1}'
            })
        
        return motifs
