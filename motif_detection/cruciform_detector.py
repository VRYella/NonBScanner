"""
CruciformDetector (strict inverted-repeat search)

Detects inverted repeats (potential cruciform-forming) with:
 - arm length >= 6 (no explicit upper cutoff)
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
    MAX_LOOP = 100       # maximum loop (spacer) length
    MAX_MISMATCHES = 0   # allowed mismatches between arm and RC(arm). Set >0 to allow imperfect arms.

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

        # To avoid excessive work on very long sequences, we limit max arm tested to length that fits
        # left_start from 0..n- (2*min_arm + 0) ; left_arm_len can be up to (n - left_start - min_arm - min_loop)
        # We'll iterate left start and arm_len and loop_len — complexity O(n * possible_arm_lengths * max_loop).
        # This is conservative but acceptable for typical genomic windows. For very long sequences you may
        # want to chunk the input.
        for left_start in range(0, n - 2 * min_arm):
            # Maximum possible arm length at this left_start given minimal right arm of min_arm and loop 0
            max_possible_arm = n - left_start - min_arm
            # But we must leave room for loop and right arm; cap arm length to reasonable bound
            # Here arm_len ranges from min_arm to max_possible_arm//2 (ensures room for right arm)
            arm_max_bound = max(min_arm, (n - left_start) // 2)
            for arm_len in range(min_arm, arm_max_bound + 1):
                left_end = left_start + arm_len
                # right arm must start at least left_end + 0 loop; but loop cannot exceed max_loop
                right_start_min = left_end
                right_start_max = min(left_end + max_loop, n - arm_len)
                # iterate right_starts (loop sizes)
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
                    # If exact-match required and mismatches > 0, we can continue scanning other right_starts
                    # Continue until right_start_max.
        # Sort hits by descending score, then by left_start
        hits.sort(key=lambda h: (-h['score'], h['left_start'], -h['arm_len']))
        return hits

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
