import re
from typing import List, Dict, Any, Tuple

try:
    from .base_detector import BaseMotifDetector
except ImportError:
    class BaseMotifDetector: pass

# Helper: reverse complement
def _rc(seq: str) -> str:
    trans = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(trans)[::-1]

VALIDATED_SEQS = [
    ("IM_VAL_001", "CCCCTCCCCTCCCCTCCCC", "Validated i-motif sequence 1", "Gehring 1993"),
    ("IM_VAL_002", "CCCCACCCCACCCCACCCC", "Validated i-motif sequence 2", "Leroy 1995"),
]

MIN_REGION_LEN = 10
CLASS_PRIORITIES = {'canonical_imotif': 1, 'relaxed_imotif': 2, 'ac_motif_hur': 3}

def _class_prio_idx(class_name: str) -> int:
    return CLASS_PRIORITIES.get(class_name, 999)

class IMotifDetector(BaseMotifDetector):
    """Detector for i-motif DNA structures (updated: Hur et al. 2021, Benabou 2014)[web:92][web:98][web:100]"""

    def get_motif_class_name(self) -> str:
        return "i-Motif"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        return {
            'canonical_imotif': [
                (r'C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}', 'IM_7_1', 'Canonical i-motif', 'canonical_imotif', 15, 'imotif_score', 0.95, 'pH-dependent C-rich structure', 'Gehring 1993'),
                (r'C{4,}[ATGC]{1,8}C{4,}[ATGC]{1,8}C{4,}[ATGC]{1,8}C{4,}', 'IM_7_2', 'High-density i-motif', 'canonical_imotif', 16, 'imotif_score', 0.98, 'Stable i-motif', 'Leroy 1995'),
            ],
            'relaxed_imotif': [
                (r'C{2,}[ATGC]{1,15}C{2,}[ATGC]{1,15}C{2,}[ATGC]{1,15}C{2,}', 'IM_7_3', 'Relaxed i-motif', 'relaxed_imotif', 12, 'imotif_score', 0.80, 'Potential i-motif structures', 'Mergny 1995'),
            ]
        }

    def find_validated_matches(self, sequence: str, check_revcomp: bool = False) -> List[Dict[str, Any]]:
        seq = sequence.upper()
        out = []
        for vid, vseq, desc, cite in VALIDATED_SEQS:
            idx = seq.find(vseq)
            if idx >= 0:
                out.append({'id': vid, 'seq': vseq, 'start': idx, 'end': idx+len(vseq), 'strand': '+', 'desc': desc, 'cite': cite})
            elif check_revcomp:
                rc = _rc(vseq)
                idx2 = seq.find(rc)
                if idx2 >= 0:
                    out.append({'id': vid, 'seq': vseq, 'start': idx2, 'end': idx2+len(vseq), 'strand': '-', 'desc': desc, 'cite': cite})
        return out

    def find_hur_ac_candidates(self, sequence: str, scan_rc: bool = True) -> List[Dict[str, Any]]:
        seq = sequence.upper()
        candidates = []

        def _matches_hur_ac(target, strand):
            for nlink in (4, 5, 6):
                # A at start, or A at end
                pat1 = r"A{3}[ACGT]{%d}C{3}[ACGT]{%d}C{3}[ACGT]{%d}C{3}" % (nlink, nlink, nlink)
                pat2 = r"C{3}[ACGT]{%d}C{3}[ACGT]{%d}C{3}[ACGT]{%d}A{3}" % (nlink, nlink, nlink)
                for pat in (pat1, pat2):
                    for m in re.finditer(pat, target):
                        s, e = m.span()
                        matched = m.group(0).upper()
                        candidates.append({
                            'start': s if strand == '+' else len(seq) - e,
                            'end': e if strand == '+' else len(seq) - s,
                            'strand': strand,
                            'linker': nlink,
                            'pattern': pat,
                            'matched_seq': matched,
                            'loose_mode': True,
                            'high_confidence': (nlink == 4 or nlink == 5)
                        })

        _matches_hur_ac(seq, '+')
        if scan_rc:
            _matches_hur_ac(_rc(seq), '-')
        candidates.sort(key=lambda x: x['start'])
        return candidates

    def _find_regex_candidates(self, sequence: str) -> List[Dict[str, Any]]:
        seq = sequence.upper()
        patterns = self.get_patterns()
        out = []
        for class_name, pats in patterns.items():
            for patt in pats:
                regex = patt[0]
                pid = patt[1] if len(patt) > 1 else f"{class_name}_pat"
                for m in re.finditer(regex, seq, flags=re.I):
                    s, e = m.start(), m.end()
                    if (e - s) < MIN_REGION_LEN:
                        continue
                    out.append({
                        'class_name': class_name,
                        'pattern_id': pid,
                        'start': s,
                        'end': e,
                        'matched_seq': seq[s:e]
                    })
        return out

    def _score_imotif_candidate(self, matched_seq: str) -> float:
        region = matched_seq.upper()
        L = len(region)
        if L < 12:
            return 0.0
        c_tracts = [m.group(0) for m in re.finditer(r"C{2,}", region)]
        if len(c_tracts) < 3:
            return 0.0
        total_c = sum(len(t) for t in c_tracts)
        c_density = total_c / L
        tract_bonus = min(0.4, 0.12 * (len(c_tracts) - 2))
        score = max(0.0, min(1.0, c_density + tract_bonus))
        return float(score)

    def _score_hur_ac_candidate(self, matched_seq: str, linker: int, high_confidence: bool) -> float:
        r = matched_seq.upper()
        L = len(r)
        ac_count = r.count('A') + r.count('C')
        ac_frac = ac_count / L if L > 0 else 0.0
        a_tracts = [len(m.group(0)) for m in re.finditer(r"A{2,}", r)]
        c_tracts = [len(m.group(0)) for m in re.finditer(r"C{2,}", r)]
        tract_score = 0.0
        if any(x >= 3 for x in a_tracts) and sum(1 for x in c_tracts if x >= 3) >= 3:
            tract_score = 0.5
        base = min(0.6, ac_frac * 0.8)
        linker_boost = 0.25 if high_confidence else (0.12 if linker == 6 else 0.0)
        return max(0.0, min(1.0, base + tract_score + linker_boost))

    def _resolve_overlaps_greedy(self, scored: List[Dict[str, Any]], merge_gap: int = 0) -> List[Dict[str, Any]]:
        if not scored:
            return []
        scored_sorted = sorted(scored, key=lambda x: (-x['score'], _class_prio_idx(x.get('class_name','')), -(x['end']-x['start'])))
        accepted = []
        occupied = []
        for cand in scored_sorted:
            s, e = cand['start'], cand['end']
            conflict = False
            for (as_, ae) in occupied:
                if not (e <= as_ - merge_gap or s >= ae + merge_gap):
                    conflict = True
                    break
            if not conflict:
                accepted.append(cand)
                occupied.append((s, e))
        accepted.sort(key=lambda x: x['start'])
        return accepted

    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
        seq = sequence.upper()
        validated = self.find_validated_matches(seq, check_revcomp=False)
        if validated:
            return 0.99
        hur_cands = self.find_hur_ac_candidates(seq, scan_rc=True)
        hur_scored = [dict(
            class_name='ac_motif_hur',
            pattern_id=h['pattern'],
            start=h['start'],
            end=h['end'],
            matched_seq=h['matched_seq'],
            linker=h['linker'],
            high_confidence=h['high_confidence'],
            score=self._score_hur_ac_candidate(h['matched_seq'], h['linker'], h['high_confidence']),
            details=h
        ) for h in hur_cands]
        regex_cands = self._find_regex_candidates(seq)
        regex_scored = [dict(
            class_name=r['class_name'],
            pattern_id=r['pattern_id'],
            start=r['start'],
            end=r['end'],
            matched_seq=r['matched_seq'],
            score=self._score_imotif_candidate(r['matched_seq']),
            details={}
        ) for r in regex_cands]
        combined = hur_scored + regex_scored
        accepted = self._resolve_overlaps_greedy(combined, merge_gap=0)
        total = float(sum(a['score'] * max(1, (a['end']-a['start'])/10.0) for a in accepted))
        return total

    def annotate_sequence(self, sequence: str) -> Dict[str, Any]:
        seq = sequence.upper()
        res = {}
        res['validated_matches'] = self.find_validated_matches(seq, check_revcomp=True)
        hur_cands = self.find_hur_ac_candidates(seq, scan_rc=True)
        for h in hur_cands:
            h['score'] = self._score_hur_ac_candidate(h['matched_seq'], h['linker'], h['high_confidence'])
        res['hur_candidates'] = hur_cands
        regex_cands = self._find_regex_candidates(seq)
        for r in regex_cands:
            r['score'] = self._score_imotif_candidate(r['matched_seq'])
        res['regex_matches'] = regex_cands
        combined = [dict(class_name='ac_motif_hur', start=h['start'], end=h['end'], score=h['score'], details=h) for h in hur_cands]
        combined += [dict(class_name=r['class_name'], start=r['start'], end=r['end'], score=r['score'], details=r) for r in regex_cands]
        res['accepted'] = self._resolve_overlaps_greedy(combined, merge_gap=0)
        return res
