import re
try:
    import re2 as re2
    RE = re2
except Exception:
    RE = re
import ahocorasick
import numpy as np
from numba import njit

G4_PATTERN = re.compile(r'G{2,}[ATGC]{1,7}G{2,}[ATGC]{1,7}G{2,}[ATGC]{1,7}G{2,}', flags=re.IGNORECASE)
A = ahocorasick.Automaton()
for seed in ("GG", "GGG", "GGGG"):
    A.add_word(seed, seed)
A.make_automaton()

def find_g4_candidates(seq: str, seed_window=50):
    seen = set()
    for end_idx, seed in A.iter(seq):
        start_win = max(0, end_idx - seed_window)
        end_win = min(len(seq), end_idx + seed_window)
        key = (start_win, end_win)
        if key in seen:
            continue
        seen.add(key)
        region = seq[start_win:end_win]
        for m in G4_PATTERN.finditer(region):
            s = start_win + m.start()
            e = start_win + m.end()
            yield {'start': s, 'end': e, 'matched_seq': m.group()}

def extract_g4_features(mseq: str):
    import re as _re
    parts = _re.findall(r'(G+)|([^G]+)', mseq, flags=_re.IGNORECASE)
    g_tracts, loop_lengths = [], []
    for gpart, loop in parts:
        if gpart:
            g_tracts.append(len(gpart))
        else:
            loop_lengths.append(len(loop))
    return {
        'length': len(mseq),
        'g_tracts': g_tracts,
        'loop_lengths': loop_lengths,
        'mean_tract': float(np.mean(g_tracts)) if g_tracts else 0.0,
        'mean_loop': float(np.mean(loop_lengths)) if loop_lengths else 0.0
    }

@njit
def _score_g4_vectorized(mean_tracts, mean_loops, w1, w2):
    n = mean_tracts.shape[0]
    out = np.empty(n, dtype=np.float64)
    for i in range(n):
        raw = w1*(mean_tracts[i] - 3.0) - w2*(mean_loops[i] - 3.0)
        if raw >= 0:
            ex = np.exp(-raw)
            out[i] = 1.0 / (1.0 + ex)
        else:
            ex = np.exp(raw)
            out[i] = ex / (1.0 + ex)
    return out

def score_g4_candidates(feats_list):
    n = len(feats_list)
    mean_tracts = np.zeros(n, dtype=np.float64)
    mean_loops = np.zeros(n, dtype=np.float64)
    for i, f in enumerate(feats_list):
        mean_tracts[i] = f['mean_tract']
        mean_loops[i] = f['mean_loop'] if f['mean_loop'] > 0 else 999.0
    return _score_g4_vectorized(mean_tracts, mean_loops, 1.2, 0.6)

def fast_scan_and_score_g4(seq: str, name: str):
    cands = list(find_g4_candidates(seq))
    feats = [extract_g4_features(c['matched_seq']) for c in cands]
    if not feats:
        return []
    scores = score_g4_candidates(feats)
    return [
        {
            'Class': 'G-Quadruplex',
            'Subclass': 'G4_Canonical',
            'Start': int(c['start']),
            'End': int(c['end']),
            'Length': f['length'],
            'matched_seq': c['matched_seq'],
            'Score': float(sc),
            'Actual_Score': float(sc),
            'Sequence_Name': name
        }
        for c, f, sc in zip(cands, feats, scores)
    ]
