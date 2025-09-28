"""
Slipped DNA Motif Detector (Full Region Capture)
------------------------------------------------
Detects and annotates complete repeat regions, following:
- STRs: Unit size 1–9 bp, ≥10 bp in span, non-overlapping, match full region[web:79][web:78]
- Direct repeats: Repeat unit 10–300 bp, max 10 bp spacer, region = repeat1 + spacer + repeat2, non-overlapping[web:79]

References:
Wells 2005, Schlötterer 2000, Weber 1989, Verkerk 1991[web:79][web:78]
"""

import re
from typing import List, Dict, Any, Tuple

try:
    from .base_detector import BaseMotifDetector
except ImportError:
    class BaseMotifDetector: pass

class SlippedDNADetector(BaseMotifDetector):
    """Detector for slipped DNA motifs: captures full repeat regions[web:79]"""

    def get_motif_class_name(self) -> str:
        return "Slipped_DNA"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        patterns = []
        # STRs: 1-9 nt units, ≥10 bp span
        for k in range(1, 10):
            pat = (
                rf"((?:[ATGC]{{{k}}}){{3,}})",   # At least 3 consecutive unit matches
                f"SLP_STR_{k}",
                f"{k}-mer STR",
                "STR",
                max(10, 3 * k),
                "instability_score",
                0.90 if k == 1 else 0.85 if k == 2 else 0.80,
                "Tandem repeat region (full span)",
                "Wells 2005"
            )
            patterns.append(pat)
        # Direct repeats: unit 10–300 nt, spacer ≤10 bp
        direct_pat = (
            r"([ATGC]{10,300})([ATGC]{0,10})\1",
            "SLP_DIR_1",
            "Direct repeat (10–300 bp units, ≤10bp spacer, full region)",
            "Direct_Repeat",
            20,
            "repeat_score",
            0.85,
            "Direct repeat (full region including spacer)",
            "Wells 2005"
        )
        return {
            "short_tandem_repeats": patterns,
            "direct_repeats": [direct_pat]
        }

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        seq = sequence.upper()
        regions = []
        pat_groups = self.get_patterns()
        used = [False] * len(seq)  # Mark used bases for overlap resolution

        # STRs (unit 1–9 bp)
        for patinfo in pat_groups["short_tandem_repeats"]:
            regex = patinfo[0]
            unit_len = int(patinfo[1].split("_")[-1])
            for m in re.finditer(regex, seq):
                s, e = m.span()
                repeat_region = seq[s:e]
                # Confirm repeat length, only accept ≥10 bp
                if (e - s) < 10:
                    continue
                # Avoid overlapping matches
                if any(used[s:e]):
                    continue
                # Mark as used
                for i in range(s, e):
                    used[i] = True
                # Count unit number (sanity check)
                n_units = (e - s) // unit_len
                regions.append({
                    'class_name': patinfo[3],
                    'pattern_id': patinfo[1],
                    'start': s,
                    'end': e,
                    'length': e - s,
                    'score': self._instability_score(repeat_region),
                    'matched_seq': repeat_region,
                    'details': {
                        'unit_length': unit_len,
                        'repeat_units': n_units,
                        'repeat_type': patinfo[2],
                        'source': patinfo[8]
                    }
                })

        # Direct repeats (unit 10–300 bp + ≤10bp spacer)
        for patinfo in pat_groups["direct_repeats"]:
            regex = patinfo[0]
            for m in re.finditer(regex, seq):
                s, e = m.span()
                if (e - s) < 20:
                    continue
                if any(used[s:e]):
                    continue
                for i in range(s, e):
                    used[i] = True
                unit = m.group(1)
                spacer = m.group(2)
                regions.append({
                    'class_name': patinfo[3],
                    'pattern_id': patinfo[1],
                    'start': s,
                    'end': e,
                    'length': e - s,
                    'score': self._repeat_score(seq[s:e]),
                    'matched_seq': seq[s:e],
                    'details': {
                        'unit_length': len(unit),
                        'spacer_length': len(spacer),
                        'repeat_type': patinfo[2],
                        'source': patinfo[8]
                    }
                })
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
        """Score based on unit count and size (normalized)"""
        N = len(sequence)
        max_instability = 0
        for unit_length in range(1, 10):
            for i in range(N - unit_length * 3 + 1):
                unit = sequence[i:i + unit_length]
                if 'N' in unit:
                    continue
                count = 1
                pos = i + unit_length
                while pos + unit_length <= N and sequence[pos:pos + unit_length] == unit:
                    count += 1
                    pos += unit_length
                length = count * unit_length
                if count >= 3 and length >= 10:
                    instability = count * (unit_length ** 0.5)
                    max_instability = max(max_instability, instability)
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
