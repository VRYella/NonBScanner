"""
G-Quadruplex Detector (G4Hunter-based, overlap-resolved)

- Uses the regex pattern groups present in get_patterns() to find candidate regions.
- Scores candidates using a G4Hunter-like algorithm (per-base G=+1, C=-1, others 0,
  sliding window aggregation + tract bonus), producing a numeric score in 0..1.
- Resolves overlaps between classes by accepting candidates in descending score order
  (higher-score regions get precedence). The 'g_triplex' class is given lowest
  priority (ties broken by class priority list).
- Exposes:
    - get_patterns() (keeps original pattern metadata)
    - calculate_score(sequence, pattern_info): returns total sum of accepted region scores
    - annotate_sequence(sequence): returns detailed region annotations (class, start/end, score, matched_text)
- No external dependencies (pure Python). Replace with Hyperscan-based matching if desired.

Notes:
- Tuning knobs at top of file: WINDOW_SIZE (G4Hunter window) and CLASS_PRIORITY (ordering).
- This module expects BaseMotifDetector to be present in same package.
"""

import re
from typing import List, Dict, Any, Tuple
from .base_detector import BaseMotifDetector

# -----------------------
# Tuning parameters
# -----------------------
WINDOW_SIZE_DEFAULT = 25  # G4Hunter sliding window size used to compute per-region score
MIN_REGION_LEN = 8        # minimum length for candidate evaluation (can be tuned)
CLASS_PRIORITY = [
    "canonical_g4",
    "relaxed_g4",
    "multimeric_g4",
    "bulged_g4",
    "imperfect_g4",
    "bipartite_g4",
    "g_triplex",  # lowest priority: triplex gets least preference
]


class GQuadruplexDetector(BaseMotifDetector):
    """Detector for G-quadruplex DNA motifs using G4Hunter scoring and overlap resolution."""

    def get_motif_class_name(self) -> str:
        return "G-Quadruplex"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return G-quadruplex DNA regex patterns (keeps existing metadata)."""
        return {
            'canonical_g4': [
                (r'G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G4_6_1', 'Canonical G4', 'Canonical G4', 15, 'g4hunter_score', 0.95, 'Stable G4 structures', 'Burge 2006'),
                (r'G{4,}[ATGC]{1,5}G{4,}[ATGC]{1,5}G{4,}[ATGC]{1,5}G{4,}', 'G4_6_2', 'High-density G4', 'Canonical G4', 16, 'g4hunter_score', 0.98, 'Very stable G4', 'Todd 2005'),
            ],
            'relaxed_g4': [
                (r'G{2,}[ATGC]{1,12}G{2,}[ATGC]{1,12}G{2,}[ATGC]{1,12}G{2,}', 'G4_6_3', 'Relaxed G4', 'Relaxed G4', 12, 'g4hunter_score', 0.80, 'Potential G4 structures', 'Huppert 2005'),
                (r'G{3,}[ATGC]{8,15}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G4_6_4', 'Long-loop G4', 'Relaxed G4', 18, 'g4hunter_score', 0.75, 'Alternative G4 topology', 'Phan 2006'),
            ],
            'bulged_g4': [
                (r'G{3,}[ATGC]{8,25}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G4_6_5', 'Bulged G4', 'Bulged G4', 20, 'g4hunter_score', 0.85, 'G4 with bulge loops', 'Lim 2009'),
                (r'G{2,}[ATGC]{15,40}G{2,}[ATGC]{1,7}G{2,}[ATGC]{1,7}G{2,}', 'G4_6_6', 'Large bulge G4', 'Bulged G4', 25, 'g4hunter_score', 0.70, 'Extended bulge G4', 'Adrian 2014'),
            ],
            'bipartite_g4': [
                (r'G{2,}[ATGC]{15,70}G{2,}[ATGC]{1,7}G{2,}[ATGC]{1,7}G{2,}', 'G4_6_7', 'Bipartite G4', 'Bipartite G4', 30, 'g4hunter_score', 0.75, 'Two-block G4 structures', 'Guédin 2010'),
            ],
            'multimeric_g4': [
                (r'(?:G{3,}[ATGC]{1,7}){4,}G{3,}', 'G4_6_8', 'Multimeric G4', 'Multimeric G4', 25, 'g4hunter_score', 0.90, 'Multiple G4 units', 'Phan 2007'),
                (r'(?:G{2,}[ATGC]{1,10}){5,}G{2,}', 'G4_6_9', 'Extended multimeric G4', 'Multimeric G4', 30, 'g4hunter_score', 0.85, 'Long G4 arrays', 'Maizels 2006'),
            ],
            'imperfect_g4': [
                (r'G{2,}[ATGC]{1,10}[AG]G{1,3}[ATGC]{1,10}G{2,}[ATGC]{1,10}G{2,}', 'G4_6_10', 'Imperfect G4', 'Imperfect G4', 15, 'g4hunter_score', 0.65, 'G4-like with interruptions', 'Kuryavyi 2010'),
                (r'G{3,}[ATGC]{1,7}[AG]{1,2}G{2,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G4_6_11', 'G-rich imperfect', 'Imperfect G4', 18, 'g4hunter_score', 0.70, 'Interrupted G-tracts', 'Webba da Silva 2007'),
            ],
            'g_triplex': [
                (r'G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G4_6_12', 'G-Triplex', 'G-Triplex intermediate', 12, 'g_triplex_score', 0.80, 'Three G-tract structures', 'Sen 1988'),
                (r'G{4,}[ATGC]{1,5}G{4,}[ATGC]{1,5}G{4,}', 'G4_6_13', 'High-density G-triplex', 'G-Triplex intermediate', 14, 'g_triplex_score', 0.85, 'Stable three-tract G-structure', 'Williamson 1989'),
            ]
        }

    # -------------------------
    # Public API
    # -------------------------
    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
        """
        Scan sequence for all pattern groups, score candidates, resolve overlaps by score
        and return total sum of accepted region scores (sum of region scores, not normalized).
        """
        seq = sequence.upper()
        candidates = self._find_all_candidates(seq)
        if not candidates:
            return 0.0
        # score candidates
        scored = [self._score_candidate(c, seq) for c in candidates]
        # resolve overlaps by selecting highest-scoring candidates first (respecting class priority tie-break)
        accepted = self._resolve_overlaps(scored)
        # total score sum
        total = sum(a['score'] for a in accepted)
        return float(total)

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        """
        Return list of accepted region annotations after overlap resolution.
        Each annotation includes:
          - class_name (pattern group)
          - pattern_id (metadata)
          - start, end (0-based, end-exclusive)
          - length
          - score (region G4Hunter-derived score 0..1 scaled by length; we return raw region score)
          - matched_seq
          - details: g_tracts, n_g_tracts, gc_balance, max_window_score, normalized_window_score
        """
        seq = sequence.upper()
        candidates = self._find_all_candidates(seq)
        scored = [self._score_candidate(c, seq) for c in candidates]
        accepted = self._resolve_overlaps(scored)
        anns = []
        for a in accepted:
            ann = {
                'class_name': a['class_name'],
                'pattern_id': a['pattern_id'],
                'start': a['start'],
                'end': a['end'],
                'length': a['end'] - a['start'],
                'score': round(a['score'], 6),
                'matched_seq': seq[a['start']:a['end']],
                'details': a['details']
            }
            anns.append(ann)
        return anns

    # -------------------------
    # Candidate finding & scoring
    # -------------------------
    def _find_all_candidates(self, seq: str) -> List[Dict[str, Any]]:
        """
        Run through all pattern groups and regex patterns, collect candidate regions.
        Returns list of dicts:
          { 'class_name', 'pattern_id', 'start', 'end', 'match_text' }
        """
        patt_groups = self.get_patterns()
        candidates = []
        for class_name, patterns in patt_groups.items():
            for pat in patterns:
                regex = pat[0]
                pattern_id = pat[1] if len(pat) > 1 else f"{class_name}_pat"
                for m in re.finditer(regex, seq):
                    s, e = m.start(), m.end()
                    # enforce minimum length
                    if (e - s) < MIN_REGION_LEN:
                        continue
                    candidates.append({
                        'class_name': class_name,
                        'pattern_id': pattern_id,
                        'start': s,
                        'end': e,
                        'match_text': seq[s:e]
                    })
        return candidates

    def _score_candidate(self, candidate: Dict[str, Any], seq: str, window_size: int = WINDOW_SIZE_DEFAULT) -> Dict[str, Any]:
        """
        Compute a G4Hunter-like score for the candidate region.
        Returns candidate dict extended with:
          - score: raw region score (sum of per-base redistributed window contributions)
          - details: breakdown dict
        Scoring approach:
          - For candidate region R, compute sliding-window sums where G=+1, C=-1, others=0.
          - max_window_abs = max absolute window sum across windows inside R (windows limited to window_size or region length)
          - normalized_window = max_window_abs / window_size
          - tract bonus: number and length of G-tracts (G{2,}) inside R increases score modestly
          - final normalized score = clamp(normalized_window + tract_bonus, 0.0, 1.0)
          - region_score (raw) = normalized_score * (region_length / window_size)  <-- this makes longer, dense regions accumulate larger raw scores
        """
        s = candidate['start']
        e = candidate['end']
        region = seq[s:e]
        L = len(region)
        # per-base values
        vals = []
        for ch in region:
            if ch == 'G':
                vals.append(1)
            elif ch == 'C':
                vals.append(-1)
            else:
                vals.append(0)
        # sliding window sums
        ws = min(window_size, L)
        if ws <= 0:
            max_abs = 0
        else:
            # compute first window sum
            cur = sum(vals[0:ws])
            max_abs = abs(cur)
            for i in range(1, L - ws + 1):
                cur += vals[i + ws - 1] - vals[i - 1]
                if abs(cur) > max_abs:
                    max_abs = abs(cur)
        normalized_window = (max_abs / ws) if ws > 0 else 0.0

        # G-tracts
        g_tracts = re.findall(r'G{2,}', region)
        n_g = len(g_tracts)
        total_g_len = sum(len(t) for t in g_tracts)
        # tract bonus: proportional to number and average tract length
        tract_bonus = 0.0
        if n_g >= 3:
            # modest bonus: 0.08 per tract beyond 2, scaled by mean tract length / 4
            tract_bonus = min(0.5, 0.08 * (n_g - 2) * ( (total_g_len / n_g) / 4.0 ))
        # small GC-balance penalty if many Cs dominate
        total_c = region.count('C')
        total_g = region.count('G')
        gc_balance = (total_g - total_c) / (L if L>0 else 1)
        gc_penalty = 0.0
        if gc_balance < -0.3:
            gc_penalty = 0.2  # heavy C-rich region penalized
        elif gc_balance < -0.1:
            gc_penalty = 0.1

        normalized_score = max(0.0, min(1.0, normalized_window + tract_bonus - gc_penalty))

        # raw region score scales with region length relative to window (so longer dense regions accumulate)
        region_score = normalized_score * (L / float(ws))

        details = {
            'n_g_tracts': n_g,
            'total_g_len': total_g_len,
            'gc_balance': round(gc_balance, 4),
            'max_window_abs': float(max_abs),
            'normalized_window': round(normalized_window, 6),
            'tract_bonus': round(tract_bonus, 6),
            'gc_penalty': round(gc_penalty, 6),
            'normalized_score': round(normalized_score, 6),
            'region_score': round(region_score, 6)
        }
        out = candidate.copy()
        out['score'] = float(region_score)
        out['details'] = details
        return out

    # -------------------------
    # Overlap resolution
    # -------------------------
    def _resolve_overlaps(self, scored_candidates: List[Dict[str, Any]], merge_gap: int = 0) -> List[Dict[str, Any]]:
        """
        Resolve overlaps between scored candidates.
        Strategy:
          - Sort candidates by (score desc, class_priority asc, region length desc)
          - Greedily accept candidate if it does not overlap any previously accepted region (respecting merge_gap)
          - If it overlaps, skip (so higher-score region keeps the space)
        Returns list of accepted candidates (same dicts with score/details).
        """
        if not scored_candidates:
            return []
        # attach class priority index
        def class_prio_idx(class_name):
            try:
                return CLASS_PRIORITY.index(class_name)
            except ValueError:
                return len(CLASS_PRIORITY)
        # sort
        scored_sorted = sorted(
            scored_candidates,
            key=lambda x: (-x['score'], class_prio_idx(x['class_name']), -(x['end'] - x['start']))
        )
        accepted = []
        occupied = []  # list of (s,e) accepted intervals
        for cand in scored_sorted:
            s, e = cand['start'], cand['end']
            conflict = False
            for (as_, ae) in occupied:
                # if overlap when extended by merge_gap -> conflict
                if not (e <= as_ - merge_gap or s >= ae + merge_gap):
                    conflict = True
                    break
            if not conflict:
                accepted.append(cand)
                occupied.append((s, e))
        # optionally sort accepted by position
        accepted.sort(key=lambda x: x['start'])
        return accepted
