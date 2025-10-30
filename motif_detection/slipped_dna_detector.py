"""
Slipped DNA Motif Detector (Optimized for Performance)
------------------------------------------------------
PERFORMANCE OPTIMIZATIONS:
- Uses optimized seed-and-extend k-mer index approach from repeat_scanner
- Efficient linear scanning for STRs
- Direct repeats detection with k-mer seeding (no catastrophic backtracking!)
- O(N) complexity for most operations
- Safe-guards to avoid explosion on highly-repetitive seeds

Detects and annotates complete repeat regions, following:
- STRs: Unit size 1–9 bp, ≥10 bp in span, non-overlapping, match full region
- Direct repeats: Unit length 10-300 bp, spacer <= 10 bp

References:
Wells 2005, Schlötterer 2000, Weber 1989, Verkerk 1991
"""

import re
from typing import List, Dict, Any, Tuple

try:
    from .base_detector import BaseMotifDetector
except ImportError:
    class BaseMotifDetector: pass

# Import optimized repeat scanner
try:
    from utils.repeat_scanner import find_direct_repeats, find_strs
except ImportError:
    # Fallback if import fails
    find_direct_repeats = None
    find_strs = None

class SlippedDNADetector(BaseMotifDetector):
    """Detector for slipped DNA motifs: captures full repeat regions[web:79]"""

    def get_motif_class_name(self) -> str:
        return "Slipped_DNA"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        # All detection now done via optimized repeat_scanner
        # Keep patterns for metadata/compatibility but don't use for regex matching
        return {
            "short_tandem_repeats": [],
            "direct_repeats": []
        }

    def find_direct_repeats_fast(self, seq: str, used: List[bool]) -> List[Dict[str, Any]]:
        """
        Fast algorithmic detection of direct repeats without catastrophic backtracking.
        Scans for unit 10-30 bp (limited for performance), spacer ≤10 bp.
        """
        regions = []
        n = len(seq)
        
        # PERFORMANCE: Adaptive parameters based on sequence length
        if n > 100000:
            max_unit = min(20, n // 3)  # Reduced unit size
            step_size = max(5, n // 20000)
            max_spacer = 5  # Reduced spacer size
        elif n > 50000:
            max_unit = min(25, n // 3)
            step_size = 4
            max_spacer = 8
        elif n >= 10000:
            max_unit = min(30, n // 3)
            step_size = 3
            max_spacer = 10
        else:
            max_unit = min(30, n // 3)
            step_size = 1
            max_spacer = 10
        
        for unit_len in range(10, max_unit + 1):
            for i in range(0, n - 2 * unit_len, step_size):
                unit = seq[i:i + unit_len]
                if 'N' in unit:
                    continue
                    
                # Check for repeat with 0 to max_spacer bp spacer
                for spacer_len in range(max_spacer + 1):
                    j = i + unit_len + spacer_len
                    if j + unit_len > n:
                        break
                    
                    if seq[j:j + unit_len] == unit:
                        start = i
                        end = j + unit_len
                        
                        # Check if already used
                        if any(used[start:end]):
                            continue
                        
                        # Mark as used
                        for k in range(start, end):
                            used[k] = True
                        
                        regions.append({
                            'class_name': 'Direct_Repeat',
                            'pattern_id': 'SLP_DIR_1',
                            'start': start,
                            'end': end,
                            'length': end - start,
                            'score': min(unit_len / 30.0, 0.95),  # Simple fast score
                            'matched_seq': seq[start:end],
                            'details': {
                                'unit_length': unit_len,
                                'spacer_length': spacer_len,
                                'repeat_type': f'Direct repeat ({unit_len} bp unit, {spacer_len} bp spacer)',
                                'source': 'Wells 2005'
                            }
                        })
                        # Found a match, skip to next position
                        break
        
        return regions

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        seq = sequence.upper()
        regions = []

        # Use optimized repeat_scanner if available
        if find_strs and find_direct_repeats:
            # STRs (unit 1–9 bp)
            str_results = find_strs(seq, min_u=1, max_u=9, min_total=10)
            for str_rec in str_results:
                regions.append({
                    'class_name': 'STR',
                    'pattern_id': f'SLP_STR_{str_rec["Unit_Length"]}',
                    'start': str_rec['Start'] - 1,  # Convert to 0-based
                    'end': str_rec['End'],
                    'length': str_rec['Length'],
                    'score': self._instability_score(str_rec['Sequence']),
                    'matched_seq': str_rec['Sequence'],
                    'details': {
                        'unit_length': str_rec['Unit_Length'],
                        'repeat_units': str_rec['Copies'],
                        'repeat_type': f"{str_rec['Unit_Length']}-mer STR",
                        'source': 'Wells 2005'
                    }
                })

            # Direct repeats
            direct_results = find_direct_repeats(seq, min_unit=10, max_unit=300, max_spacer=10)
            for dr_rec in direct_results:
                regions.append({
                    'class_name': 'Direct_Repeat',
                    'pattern_id': 'SLP_DIR_1',
                    'start': dr_rec['Start'] - 1,  # Convert to 0-based
                    'end': dr_rec['End'],
                    'length': dr_rec['Length'],
                    'score': min(dr_rec['Unit_Length'] / 30.0, 0.95),
                    'matched_seq': dr_rec['Sequence'],
                    'details': {
                        'unit_length': dr_rec['Unit_Length'],
                        'spacer_length': dr_rec['Spacer'],
                        'repeat_type': f'Direct repeat ({dr_rec["Unit_Length"]} bp unit, {dr_rec["Spacer"]} bp spacer)',
                        'source': 'Wells 2005'
                    }
                })
        else:
            # Fallback to old implementation if imports fail
            used = [False] * len(seq)
            pat_groups = self.get_patterns()

            # STRs (unit 1–9 bp) - fallback regex
            for k in range(1, 10):
                regex = rf"((?:[ATGC]{{{k}}}){{3,}})"
                for m in re.finditer(regex, seq):
                    s, e = m.span()
                    if (e - s) < 10 or any(used[s:e]):
                        continue
                    for i in range(s, e):
                        used[i] = True
                    n_units = (e - s) // k
                    regions.append({
                        'class_name': 'STR',
                        'pattern_id': f'SLP_STR_{k}',
                        'start': s,
                        'end': e,
                        'length': e - s,
                        'score': self._instability_score(seq[s:e]),
                        'matched_seq': seq[s:e],
                        'details': {
                            'unit_length': k,
                            'repeat_units': n_units,
                            'repeat_type': f'{k}-mer STR',
                            'source': 'Wells 2005'
                        }
                    })

            # Direct repeats - fallback
            direct_regions = self.find_direct_repeats_fast(seq, used)
            regions.extend(direct_regions)
        
        # Sort by start
        regions.sort(key=lambda r: r['start'])
        return regions

    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        scoring_method = pattern_info[5] if len(pattern_info) > 5 else 'instability_score'
        if scoring_method == 'instability_score':
            return self._instability_score(sequence)
        elif scoring_method == 'repeat_score':
            return self._repeat_score(sequence)
        else:
            return 0.0

    def _instability_score(self, sequence: str) -> float:
        """Score based on unit count and size (normalized)
        
        Optimized version: Since sequence is already a matched repeat region,
        we can compute the score more efficiently.
        """
        N = len(sequence)
        if N < 10:
            return 0.0
        
        max_instability = 0
        # PERFORMANCE: Limit the search range for large sequences
        # Only check first 100 bp to determine the repeat unit
        search_len = min(N, 100)
        
        for unit_length in range(1, min(10, search_len // 3 + 1)):
            for i in range(min(search_len - unit_length * 3 + 1, 20)):  # Only check first 20 positions
                unit = sequence[i:i + unit_length]
                if 'N' in unit:
                    continue
                count = 1
                pos = i + unit_length
                # Only count up to reasonable limit
                max_check = min(N, i + unit_length * 50)  # Limit to 50 repeats
                while pos + unit_length <= max_check and sequence[pos:pos + unit_length] == unit:
                    count += 1
                    pos += unit_length
                length = count * unit_length
                if count >= 3 and length >= 10:
                    instability = count * (unit_length ** 0.5)
                    max_instability = max(max_instability, instability)
                    # If we found a good score, we can stop early
                    if max_instability > 8:
                        return 1.0
        return min(max_instability / 10, 1.0)

    def _repeat_score(self, sequence: str) -> float:
        """Score for direct repeats (unit size, short spacer normalized)"""
        N = len(sequence)
        max_score = 0
        for unit_size in range(10, min(301, N // 2 + 1)):
            for i in range(N - 2 * unit_size - 10 + 1):
                unit = sequence[i:i + unit_size]
                for spacer_size in range(0, min(11, N - (i + 2 * unit_size))):
                    j = i + unit_size + spacer_size
                    if j + unit_size > N:
                        break
                    if sequence[j:j + unit_size] == unit:
                        score = unit_size / (1 + spacer_size)
                        max_score = max(max_score, score)
        return min(max_score / 50, 1.0)
