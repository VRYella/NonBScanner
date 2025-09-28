"""
Z-DNA Motif Detector (10-mer table)
==================================

Detects Z-DNA-like 10-mer motifs using a provided motif -> score table.

Behavior:
 - Use Hyperscan (if available) for very fast matching of the 10-mer list.
 - Fallback to a pure-Python exact matcher if Hyperscan isn't installed.
 - Merge overlapping/adjacent 10-mer matches into contiguous regions.
 - Redistribute each 10-mer's score evenly across its 10 bases (score/10 per base).
 - Region sum_score = sum of per-base contributions inside merged region.
 - calculate_score returns total sum_score across merged regions.
"""

import re
from typing import List, Dict, Any, Tuple
from .base_detector import BaseMotifDetector

# optional hyperscan
try:
    import hyperscan
    _HYPERSCAN_AVAILABLE = True
except Exception:
    _HYPERSCAN_AVAILABLE = False


class ZDNADetector(BaseMotifDetector):
    """Detector for Z-DNA motifs using a 10-mer scoring table."""

    # -------------------------
    # Full provided 10-mer -> score table (paste as-is)
    # -------------------------
    TENMER_SCORE: Dict[str, float] = {
        "AACGCGCGCG": 50.25,
        "ATGCGCGCGC": 51.25,
        "ATCGCGCGCG": 50.0,
        "AGCGCGCGCA": 50.25,
        "AGCGCGCGCG": 56.0,
        "ACGGCGCGCG": 50.25,
        "ACGCGGCGCG": 50.25,
        "ACGCGCGGCG": 50.25,
        "ACGCGCGCGA": 50.25,
        "ACGCGCGCGT": 51.5,
        "ACGCGCGCGG": 50.25,
        "ACGCGCGCGC": 57.25,
        "ACGCGCGCCG": 50.25,
        "ACGCGCCGCG": 50.25,
        "ACGCCGCGCG": 50.25,
        "ACCGCGCGCG": 50.25,
        "TAGCGCGCGC": 50.0,
        "TACGCGCGCG": 51.25,
        "TTGCGCGCGC": 50.25,
        "TGGCGCGCGC": 50.25,
        "TGCGGCGCGC": 50.25,
        "TGCGCGGCGC": 50.25,
        "TGCGCGCGGC": 50.25,
        "TGCGCGCGCA": 51.5,
        "TGCGCGCGCT": 50.25,
        "TGCGCGCGCG": 57.25,
        "TGCGCGCGCC": 50.25,
        "TGCGCGCCGC": 50.25,
        "TGCGCCGCGC": 50.25,
        "TGCCGCGCGC": 50.25,
        "TCGCGCGCGT": 50.25,
        "TCGCGCGCGC": 56.0,
        "GACGCGCGCG": 50.25,
        "GTGCGCGCGC": 51.5,
        "GTCGCGCGCG": 50.25,
        "GGCGCGCGCA": 50.25,
        "GGCGCGCGCG": 56.0,
        "GCAGCGCGCG": 50.25,
        "GCACGCGCGC": 51.5,
        "GCTGCGCGCG": 50.25,
        "GCGACGCGCG": 50.25,
        "GCGTGCGCGC": 51.5,
        "GCGTCGCGCG": 50.25,
        "GCGGCGCGCA": 50.25,
        "GCGGCGCGCG": 56.0,
        "GCGCAGCGCG": 50.25,
        "GCGCACGCGC": 51.5,
        "GCGCTGCGCG": 50.25,
        "GCGCGACGCG": 50.25,
        "GCGCGTGCGC": 51.5,
        "GCGCGTCGCG": 50.25,
        "GCGCGGCGCA": 50.25,
        "GCGCGGCGCG": 56.0,
        "GCGCGCAGCG": 50.25,
        "GCGCGCACGC": 51.5,
        "GCGCGCTGCG": 50.25,
        "GCGCGCGACG": 50.25,
        "GCGCGCGTGC": 51.5,
        "GCGCGCGTCG": 50.25,
        "GCGCGCGGCA": 50.25,
        "GCGCGCGGCG": 56.0,
        "GCGCGCGCAA": 50.25,
        "GCGCGCGCAT": 51.25,
        "GCGCGCGCAG": 50.25,
        "GCGCGCGCAC": 51.5,
        "GCGCGCGCTA": 50.0,
        "GCGCGCGCTG": 50.25,
        "GCGCGCGCGA": 56.0,
        "GCGCGCGCGT": 57.25,
        "GCGCGCGCGG": 56.0,
        "GCGCGCGCGC": 63.0,
        "GCGCGCGCCA": 50.25,
        "GCGCGCGCCG": 56.0,
        "GCGCGCCGCA": 50.25,
        "GCGCGCCGCG": 56.0,
        "GCGCCGCGCA": 50.25,
        "GCGCCGCGCG": 56.0,
        "GCCGCGCGCA": 50.25,
        "GCCGCGCGCG": 56.0,
        "CAGCGCGCGC": 50.25,
        "CACGCGCGCG": 51.5,
        "CTGCGCGCGC": 50.25,
        "CGACGCGCGC": 50.25,
        "CGTGCGCGCG": 51.5,
        "CGTCGCGCGC": 50.25,
        "CGGCGCGCGT": 50.25,
        "CGGCGCGCGC": 56.0,
        "CGCAGCGCGC": 50.25,
        "CGCACGCGCG": 51.5,
        "CGCTGCGCGC": 50.25,
        "CGCGACGCGC": 50.25,
        "CGCGTGCGCG": 51.5,
        "CGCGTCGCGC": 50.25,
        "CGCGGCGCGT": 50.25,
        "CGCGGCGCGC": 56.0,
        "CGCGCAGCGC": 50.25,
        "CGCGCACGCG": 51.5,
        "CGCGCTGCGC": 50.25,
        "CGCGCGACGC": 50.25,
        "CGCGCGTGCG": 51.5,
        "CGCGCGTCGC": 50.25,
        "CGCGCGGCGT": 50.25,
        "CGCGCGGCGC": 56.0,
        "CGCGCGCAGC": 50.25,
        "CGCGCGCACG": 51.5,
        "CGCGCGCTGC": 50.25,
        "CGCGCGCGAT": 50.0,
        "CGCGCGCGAC": 50.25,
        "CGCGCGCGTA": 51.25,
        "CGCGCGCGTT": 50.25,
        "CGCGCGCGTG": 51.5,
        "CGCGCGCGTC": 50.25,
        "CGCGCGCGGT": 50.25,
        "CGCGCGCGGC": 56.0,
        "CGCGCGCGCA": 57.25,
        "CGCGCGCGCT": 56.0,
        "CGCGCGCGCG": 63.0,
        "CGCGCGCGCC": 56.0,
        "CGCGCGCCGT": 50.25,
        "CGCGCGCCGC": 56.0,
        "CGCGCCGCGT": 50.25,
        "CGCGCCGCGC": 56.0,
        "CGCCGCGCGT": 50.25,
        "CGCCGCGCGC": 56.0,
        "CCGCGCGCGT": 50.25,
        "CCGCGCGCGC": 56.0,
    }

    def get_motif_class_name(self) -> str:
        return "Z-DNA"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """
        Keep compatibility with framework: return a representative pattern entry.
        The actual matching uses the TENMER_SCORE table inside calculate_score/annotate_sequence.
        """
        return {
            "z_dna_10mers": [
                (r"", "ZDN_10MER", "Z-DNA 10-mer table", "Z-DNA", 10, "z_dna_10mer_score", 0.9, "Z-DNA 10mer motif", "user_table"),
            ]
        }

    # -------------------------
    # Public API
    # -------------------------
    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        """
        Return total sum_score across merged Z-like regions in sequence.
        sum_score is computed by redistributing each matched 10-mer's score equally across its 10 bases
        and summing per-base contributions inside merged regions.
        """
        seq = sequence.upper()
        merged = self._find_and_merge_10mer_matches(seq)
        if not merged:
            return 0.0
        contrib = self._build_per_base_contrib(seq)
        total = 0.0
        for s, e in merged:
            total += sum(contrib[s:e])
        return float(total)

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        """
        Return list of merged region annotations:
          - start, end (0-based, end-exclusive)
          - length
          - sum_score (sum of per-base contributions)
          - mean_score_per10mer (mean of 10-mer scores contributing)
          - n_10mers
          - contributing_10mers: list of dicts {tenmer, start, score}
        """
        seq = sequence.upper()
        matches = self._find_10mer_matches(seq)
        if not matches:
            return []
        merged = self._merge_matches(matches)
        contrib = self._build_per_base_contrib(seq)
        annotations = []
        for (s, e, region_matches) in merged:
            sum_score = sum(contrib[s:e])
            n10 = len(region_matches)
            mean10 = (sum(m[2] for m in region_matches) / n10) if n10 > 0 else 0.0
            ann = {
                "start": s,
                "end": e,
                "length": e - s,
                "sum_score": round(sum_score, 6),
                "mean_score_per10mer": round(mean10, 6),
                "n_10mers": n10,
                "contributing_10mers": [{"tenmer": m[1], "start": m[0], "score": m[2]} for m in region_matches]
            }
            annotations.append(ann)
        return annotations

    # -------------------------
    # Core helpers
    # -------------------------
    def _find_10mer_matches(self, seq: str) -> List[Tuple[int, str, float]]:
        """Return list of (start, tenmer, score)."""
        if _HYPERSCAN_AVAILABLE:
            try:
                return self._hs_find_matches(seq)
            except Exception:
                return self._py_find_matches(seq)
        else:
            return self._py_find_matches(seq)

    def _hs_find_matches(self, seq: str) -> List[Tuple[int, str, float]]:
        """Hyperscan-based matching."""
        expressions = []
        ids = []
        id_to_ten = {}
        id_to_score = {}
        for idx, (ten, score) in enumerate(self.TENMER_SCORE.items()):
            expressions.append(ten.encode())
            ids.append(idx)
            id_to_ten[idx] = ten
            id_to_score[idx] = float(score)
        db = hyperscan.Database()
        db.compile(expressions=expressions, ids=ids, elements=len(expressions))
        matches: List[Tuple[int, str, float]] = []

        def on_match(id, start, end, flags, context):
            matches.append((start, id_to_ten[id], id_to_score[id]))

        db.scan(seq.encode(), match_event_handler=on_match)
        matches.sort(key=lambda x: x[0])
        return matches

    def _py_find_matches(self, seq: str) -> List[Tuple[int, str, float]]:
        """Pure-Python exact search (overlapping matches allowed)."""
        n = len(seq)
        matches: List[Tuple[int, str, float]] = []
        for i in range(0, n - 10 + 1):
            ten = seq[i:i + 10]
            score = self.TENMER_SCORE.get(ten)
            if score is not None:
                matches.append((i, ten, float(score)))
        return matches

    def _merge_matches(self, matches: List[Tuple[int, str, float]],
                       merge_gap: int = 0) -> List[Tuple[int, int, List[Tuple[int, str, float]]]]:
        """
        Merge matches into regions. merge_gap controls how far apart matches can be and still merge.
        Returns list of (start, end, matches_in_region). end is exclusive.
        """
        if not matches:
            return []
        merged = []
        cur_start = matches[0][0]
        cur_end = matches[0][0] + 10
        cur_matches = [matches[0]]
        for m in matches[1:]:
            s = m[0]
            m_end = s + 10
            if s <= cur_end + merge_gap:
                cur_end = max(cur_end, m_end)
                cur_matches.append(m)
            else:
                merged.append((cur_start, cur_end, cur_matches))
                cur_start, cur_end = s, m_end
                cur_matches = [m]
        merged.append((cur_start, cur_end, cur_matches))
        return merged

    def _find_and_merge_10mer_matches(self, seq: str, merge_gap: int = 0) -> List[Tuple[int, int]]:
        matches = self._find_10mer_matches(seq)
        merged = self._merge_matches(matches, merge_gap=merge_gap)
        return [(s, e) for (s, e, _) in merged]

    def _build_per_base_contrib(self, seq: str) -> List[float]:
        """
        Build per-base contribution array. For each matched 10-mer (start j, score S),
        add S/10 to bases j..j+9.
        """
        n = len(seq)
        contrib = [0.0] * n
        matches = self._find_10mer_matches(seq)
        for (start, ten, score) in matches:
            per_base = float(score) / 10.0
            for k in range(start, min(start + 10, n)):
                contrib[k] += per_base
        return contrib
