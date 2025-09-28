def get_motif_class_name(self) -> str:
    return "i-Motif"

def get_patterns(self) -> Dict[str, List[Tuple]]:
    """
    Provide canonical pattern metadata (kept for compatibility). Note:
    actual Hur-style ac motif detection is implemented in find_hur_ac_candidates().
    """
    return {
        'canonical_imotif': [
            (r'C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}', 'IM_7_1', 'Canonical i-motif', 'Canonical i-motif', 15, 'imotif_score', 0.95, 'pH-dependent C-rich structure', 'Gehring 1993'),
            (r'C{4,}[ATGC]{1,5}C{4,}[ATGC]{1,5}C{4,}[ATGC]{1,5}C{4,}', 'IM_7_2', 'High-density i-motif', 'Canonical i-motif', 16, 'imotif_score', 0.98, 'Stable i-motif', 'Leroy 1995'),
        ],
        'relaxed_imotif': [
            (r'C{2,}[ATGC]{1,12}C{2,}[ATGC]{1,12}C{2,}[ATGC]{1,12}C{2,}', 'IM_7_3', 'Relaxed i-motif', 'Relaxed i-motif', 12, 'imotif_score', 0.80, 'Potential i-motif structures', 'Mergny 1995'),
        ],
        # keep older ac heuristics (optional), but Hur detection is separate
        'ac_motif_heuristic': [
            (r'(?:AC){6,}', 'IM_7_h1', 'AC-motif (heuristic)', 'AC-motif', 10, 'imotif_score', 0.6, 'heuristic alternating AC', 'heuristic'),
        ],
    }

# ---------------- validated-sequence utilities ----------------
def find_validated_matches(self, sequence: str, check_revcomp: bool = False) -> List[Dict[str, Any]]:
    """Return exact validated-seq substring matches (forward by default)."""
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

# ---------------- Hur-specific AC candidate finder ----------------
def find_hur_ac_candidates(self, sequence: str, scan_rc: bool = True) -> List[Dict[str, Any]]:
    """
    Find candidates that match Hur et al.'s AC motif schemes:
      - primary scheme: A{3} N{6} C{3} N{6} C{3} N{6} C{3}  (Hur reported 2151 hits genome-wide)
      - narrower searches: linker length = 4 or 5 (Hur called these high-confidence)
    Returns list of dicts:
      {'start','end','strand','linker','pattern_regex','matched_seq','high_confidence'(bool)}
    Citation: Hur et al. NAR 2021 (PMC). :contentReference[oaicite:2]{index=2}
    """
    seq = sequence.upper()
    candidates = []

    def _search_with_L(L: int, target_seq: str, strand: str):
        # pattern: A3 N{L} C3 N{L} C3 N{L} C3
        pat = re.compile(r"A{3}[ACGT]{" + str(L) + r"}C{3}[ACGT]{" + str(L) + r"}C{3}[ACGT]{" + str(L) + r"}C{3}", re.I)
        for m in pat.finditer(target_seq):
            s, e = m.start(), m.end()
            matched = m.group(0).upper()
            candidates.append({
                'start': s if strand == '+' else len(seq) - e,
                'end': e if strand == '+' else len(seq) - s,
                'strand': strand,
                'linker': L,
                'pattern': f"A3N{L}C3N{L}C3N{L}C3",
                'matched_seq': matched,
                'high_confidence': (L in (4,5))
            })

    # search forward strand
    for L in HUR_LINKERS:
        _search_with_L(L, seq, '+')

    # optionally search reverse complement (report coordinates on forward reference)
    if scan_rc:
        rseq = _rc(seq)
        for L in HUR_LINKERS:
            _search_with_L(L, rseq, '-')

    # sort by start
    candidates.sort(key=lambda x: x['start'])
    return candidates

# ---------------- general candidate enumeration (regex groups) ----------------
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

# ---------------- scoring functions ----------------
def _score_imotif_candidate(self, matched_seq: str) -> float:
    """
    Simple normalized heuristic for i-motif-like (C-tract based).
    Returns 0..1. (Conservative, largely unchanged from earlier heuristics.)
    """
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
    """
    Scoring rule for Hur AC-candidates:
     - baseline: check A/C fraction and presence of A-tracts and C-tracts of length >=3
     - high_confidence (L==4 or L==5) => boost score (Hur reported many validated promo hits for L=4/5)
     - return value 0..1 (normalized)
    Implementation is conservative — it marks high L=4/5 as stronger.
    """
    r = matched_seq.upper()
    L = len(r)
    # A/C density
    ac_count = r.count('A') + r.count('C')
    ac_frac = ac_count / L if L>0 else 0.0
    # tract lengths
    a_tracts = [len(m.group(0)) for m in re.finditer(r"A{2,}", r)]
    c_tracts = [len(m.group(0)) for m in re.finditer(r"C{2,}", r)]
    # require at least 2 A-tracts and 2 C-tracts (Hur's motif has A3 + multiple C3)
    tract_score = 0.0
    if any(x >= 3 for x in a_tracts) and sum(1 for x in c_tracts if x>=3) >= 3:
        tract_score = 0.5
    # base from ac_frac
    base = min(0.6, ac_frac * 0.8)
    # linker boost for L in (4,5)
    linker_boost = 0.25 if high_confidence else (0.12 if linker==6 else 0.0)
    # final normalized
    score = max(0.0, min(1.0, base + tract_score + linker_boost))
    return float(score)

# ---------------- overlap-resolve and public API ----------------
def _resolve_overlaps_greedy(self, scored: List[Dict[str, Any]], merge_gap: int = 0) -> List[Dict[str, Any]]:
    """
    Greedy resolution: sort by (-score, class priority, length desc) and accept non-overlapping
    """
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
            occupied.append((s,e))
    accepted.sort(key=lambda x: x['start'])
    return accepted

def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
    """
    Main API: returns a total raw sum score across accepted regions (validated > Hur AC > regex matches).
    Also ensures exact validated sequences are prioritized (returns high score quickly).
    """
    seq = sequence.upper()

    # 1) exact validated sequence check (forward only by default) => high confidence
    validated = self.find_validated_matches(seq, check_revcomp=False)
    if validated:
        # return a very high normalized score (caller can also inspect validated hits via find_validated_matches)
        return 0.99

    # 2) collect Hur AC candidates
    hur_cands = self.find_hur_ac_candidates(seq, scan_rc=True)
    hur_scored = []
    for h in hur_cands:
        score = self._score_hur_ac_candidate(h['matched_seq'], h['linker'], h['high_confidence'])
        hur_scored.append({
            'class_name': 'ac_motif_hur',
            'pattern_id': h['pattern'],
            'start': h['start'],
            'end': h['end'],
            'matched_seq': h['matched_seq'],
            'linker': h['linker'],
            'high_confidence': h['high_confidence'],
            'score': score,
            'details': h
        })

    # 3) collect regex motif candidates (canonical/relaxed/heuristic)
    regex_cands = self._find_regex_candidates(seq)
    regex_scored = []
    for r in regex_cands:
        score = self._score_imotif_candidate(r['matched_seq'])
        regex_scored.append({
            'class_name': r['class_name'],
            'pattern_id': r['pattern_id'],
            'start': r['start'],
            'end': r['end'],
            'matched_seq': r['matched_seq'],
            'score': score,
            'details': {}
        })

    # 4) combine all and resolve overlaps
    combined = hur_scored + regex_scored
    accepted = self._resolve_overlaps_greedy(combined, merge_gap=0)

    # 5) total score: sum of accepted region scores (so longer/more/stronger regions increase total)
    total = float(sum(a['score'] * max(1, (a['end']-a['start'])/10.0) for a in accepted))
    return total

def annotate_sequence(self, sequence: str) -> Dict[str, Any]:
    """
    Returns a dictionary containing:
     - validated_matches (exact)
     - hur_candidates (raw)
     - regex_matches (raw)
     - accepted (post-overlap resolution, with scores and details)
    """
    seq = sequence.upper()
    res = {}
    res['validated_matches'] = self.find_validated_matches(seq, check_revcomp=True)
    hur_cands = self.find_hur_ac_candidates(seq, scan_rc=True)
    # score hur
    for h in hur_cands:
        h['score'] = self._score_hur_ac_candidate(h['matched_seq'], h['linker'], h['high_confidence'])
    res['hur_candidates'] = hur_cands
    # regex
    regex_cands = self._find_regex_candidates(seq)
    for r in regex_cands:
        r['score'] = self._score_imotif_candidate(r['matched_seq'])
    res['regex_matches'] = regex_cands
    # combine & accept
    combined = []
    for h in hur_cands:
        combined.append({'class_name': 'ac_motif_hur', 'start': h['start'], 'end': h['end'], 'score': h['score'], 'details': h})
    for r in regex_cands:
        combined.append({'class_name': r['class_name'], 'start': r['start'], 'end': r['end'], 'score': r['score'], 'details': r})
    res['accepted'] = self._resolve_overlaps_greedy(combined, merge_gap=0)
    return res
