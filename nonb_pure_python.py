#!/usr/bin/env python3
"""
| **Parameter/Function**           | **Purpose/Range**                                              | **Typical Value**       | **Relevant Detector/Scoring**             | **Notes/Comment**                                                             |
|----------------------------------|----------------------------------------------------------------|------------------------|-------------------------------------------|-------------------------------------------------------------------------------|
| `KMER_DIRECT`                    | Direct repeat seed k-mer size                                  | 12                     | Direct Repeat                            | Seed for initial match, unit >= 10                                          |
| `KMER_INVERT_SEED`               | Inverted/mirror repeat seed size                               | 8                      | Cruciform, Mirror/Triplex                 | Seed for palindrome/mirror seeding                                           |
| `MIN_DIRECT_UNIT`                | Min direct repeat unit length                                  | 10                     | Direct Repeat                            | Lower bound for repeat unit                                                  |
| `MAX_DIRECT_UNIT`                | Max direct repeat unit length                                  | 300                    | Direct Repeat                            | Upper bound for repeat unit                                                  |
| `MAX_DIRECT_SPACER`              | Max allowed spacer between direct repeats                      | 10                     | Direct Repeat                            | Max gap between repeat units                                                 |
| `STR_MIN_UNIT`                   | Min STR unit length                                            | 1                      | STR                                      | Shortest repeat unit for STR detection                                       |
| `STR_MAX_UNIT`                   | Max STR unit length                                            | 9                      | STR                                      | Longest repeat unit for STR detection                                        |
| `STR_MIN_TOTAL_LEN`              | Min total STR length                                           | 10                     | STR                                      | Ensures STRs are significant                                                 |
| `MIN_CRUCIFORM_ARM`              | Min arm length for cruciform/inverted repeat                   | 6                      | Cruciform/Inverted                       | Shortest palindrome arm detected                                             |
| `MAX_CRUCIFORM_SPACER`           | Max allowed spacer in cruciform/inverted repeats               | 100                    | Cruciform/Inverted                       | Max gap between palindrome arms                                              |
| `MIN_TRIPLEX_ARM`                | Min arm length for triplex/mirror repeat                       | 10                     | Mirror/Triplex                           | Shortest arm for triplex detection                                           |
| `MAX_TRIPLEX_SPACER`             | Max allowed spacer in triplex/mirror repeats                   | 100                    | Mirror/Triplex                           | Max gap between arms                                                         |
| `TRIPLEX_PURITY`                 | Min purine/pyrimidine purity for triplex/mirror                | 0.90                   | Mirror/Triplex                           | Fraction of A/G or C/T required in arm                                      |
| **score_direct**                 | Scores direct repeats by unit len, repeat count, and spacer    | -                      | Direct Repeat                            | Penalizes by spacer, saturates via logistic curve                            |
| **score_STR**                    | Scores STRs by total length and unit length                    | -                      | STR                                      | Favors longer repeats                                                        |
| **score_cruciform**              | Scores cruciform/inverted by arm length and spacer             | -                      | Cruciform/Inverted                       | Prefers long arms, penalizes large spacer                                    |
| **score_triplex**                | Scores triplex/mirror by arm length, purity, and spacer        | -                      | Mirror/Triplex                           | Penalizes by spacer, boosts for purity above threshold          

"""

from collections import defaultdict
import math
import sys
import itertools
import bisect

# --------------------------
# Parameters (tuneable)
# --------------------------
KMER_DIRECT = 12        # k-mer seed size for direct repeats (>=10 unit)
KMER_INVERT_SEED = 8    # k-mer seed size for inverted/mirror seeding
MIN_DIRECT_UNIT = 10
MAX_DIRECT_UNIT = 300
MAX_DIRECT_SPACER = 10

STR_MIN_UNIT = 1
STR_MAX_UNIT = 9
STR_MIN_TOTAL_LEN = 10

MIN_CRUCIFORM_ARM = 6
MAX_CRUCIFORM_SPACER = 100

MIN_TRIPLEX_ARM = 10
MAX_TRIPLEX_SPACER = 100
TRIPLEX_PURITY = 0.90

# Rolling hash params (64-bit-ish)
MASK64 = (1<<64) - 1
BASE1 = 1315423911  # odd random-ish
BASE2 = 1000003

# --------------------------
# Utilities
# --------------------------
ENC = {ord('A'):0, ord('C'):1, ord('G'):2, ord('T'):3,
       ord('a'):0, ord('c'):1, ord('g'):2, ord('t'):3}

RCMAP = str.maketrans("ACGTacgt", "TGCAtgca")
def rc_str(s): return s.translate(RCMAP)[::-1]
def rev_str(s): return s[::-1]

def seq_to_int_kmer(seq):
    """k-mer to int with 2-bit packing."""
    x = 0
    for ch in seq:
        x = (x << 2) | (ENC.get(ord(ch), 0) & 3)
    return x

def build_kmer_index(seq, k):
    """Return dict: kmer_int -> sorted list of start positions."""
    idx = defaultdict(list)
    n = len(seq)
    if n < k: return idx
    # rolling kmer via 2-bit
    mask = (1 << (2*k)) - 1
    x = seq_to_int_kmer(seq[:k])
    idx[x].append(0)
    for i in range(1, n-k+1):
        x = ((x << 2) & mask) | (ENC.get(ord(seq[i+k-1]), 0) & 3)
        idx[x].append(i)
    return idx

def build_double_hashes(seq):
    """Return (h1, p1, h2, p2) arrays where h[i] is prefix hash up to i."""
    n = len(seq)
    h1 = [0]*(n+1); p1 = [1]*(n+1)
    h2 = [0]*(n+1); p2 = [1]*(n+1)
    for i,ch in enumerate(seq):
        v = (ENC.get(ord(ch),0) + 1)  # avoid zero
        h1[i+1] = ((h1[i] * BASE1) + v) & MASK64
        p1[i+1] = (p1[i] * BASE1) & MASK64
        h2[i+1] = ((h2[i] * BASE2) + v) & MASK64
        p2[i+1] = (p2[i] * BASE2) & MASK64
    return h1, p1, h2, p2

def subhash(h, p, l, r):
    """Return hash of seq[l:r] using precomputed prefix arrays (mod 2^64)."""
    return (h[r] - (h[l] * p[r-l] & MASK64)) & MASK64

# --------------------------
# Scoring helpers
# --------------------------
def score_direct(unit_len, repeat_count, spacer, max_spacer=MAX_DIRECT_SPACER):
    # raw_score ~ repeat_count * unit_len penalized by spacer
    raw = repeat_count * unit_len
    spacer_pen = max(0.0, 1.0 - (spacer / (max_spacer+1)))
    raw_adj = raw * spacer_pen
    # normalized to 0-1 roughly: map to 0..1 using logistic-ish / saturation
    norm = raw_adj / (raw_adj + 100.0)  # tweak denominator to change scale
    return raw_adj, norm

def score_STR(total_len, unit_len):
    # prefer longer total length and shorter unit complexity somewhat
    raw = total_len
    norm = total_len / (total_len + 20.0)
    return raw, norm

def score_cruciform(arm_len, spacer, max_spacer=MAX_CRUCIFORM_SPACER):
    # prefer long arms and small spacers
    raw = arm_len
    spacer_pen = max(0.0, 1.0 - (spacer / (max_spacer + 1)))
    raw_adj = raw * spacer_pen
    norm = raw_adj / (raw_adj + 50.0)
    return raw_adj, norm

def score_triplex(arm_len, purity, spacer, max_spacer=MAX_TRIPLEX_SPACER):
    # purity should be > threshold; weight arm_len heavily
    purity_bonus = (purity - TRIPLEX_PURITY) if purity >= TRIPLEX_PURITY else (purity - TRIPLEX_PURITY)
    raw = arm_len * max(0.0, purity)
    spacer_pen = max(0.0, 1.0 - (spacer / (max_spacer + 1)))
    raw_adj = raw * spacer_pen
    norm = raw_adj / (raw_adj + 80.0)
    return raw_adj, norm

# --------------------------
# Detectors
# --------------------------

def detect_STRs(seq, min_unit=STR_MIN_UNIT, max_unit=STR_MAX_UNIT, min_total=STR_MIN_TOTAL_LEN):
    n = len(seq)
    results = []
    for p in range(min_unit, max_unit+1):
        i = 0
        while i < n:
            if i + p > n:
                break
            unit = seq[i:i+p]
            if len(unit) < p:
                break
            j = i + p
            copies = 1
            while j + p <= n and seq[j:j+p] == unit:
                copies += 1
                j += p
            total_len = copies * p
            if total_len >= min_total and copies >= 2:
                raw, norm = score_STR(total_len, p)
                results.append({
                    "class":"Slipped_DNA",
                    "subclass":"STR",
                    "start":i+1, "end":j, "length": total_len,
                    "sequence": seq[i:j],
                    "raw_score": raw, "normalized_score": norm,
                    "strand": "+",
                    "unit_len": p, "copies": copies})
                i = j
            else:
                i += 1
    return results

def detect_direct_repeats(seq, unit_min=MIN_DIRECT_UNIT, unit_max=MAX_DIRECT_UNIT, spacer_max=MAX_DIRECT_SPACER, k=KMER_DIRECT):
    n = len(seq)
    results = []
    if n < k:
        return results
    km_idx = build_kmer_index(seq, k)
    h1, p1, h2, p2 = build_double_hashes(seq)
    # iterate over kmers (seed): for each occurrence p and for candidate q we test unit lengths
    # To avoid huge loops, iterate positions and unit sizes smartly
    for kmer_int, pos_list in km_idx.items():
        # pos_list sorted ascending
        for ppos in pos_list:
            # prefer to iterate feasible unit lengths where q stays inside sequence
            max_unit_here = min(unit_max, n - ppos - k)  # ensure q+k in bounds for seed check
            for unit in range(unit_min, max_unit_here+1):
                q_base = ppos + unit
                # test possible spacers 0..spacer_max
                for spacer in range(0, spacer_max+1):
                    q = q_base + spacer
                    if q + k > n:
                        break
                    # quick k-mer equality test
                    # compute k-mer int at q (fast building not cached here - acceptable cost)
                    if seq_to_int_kmer(seq[q:q+k]) != kmer_int:
                        continue
                    # verify full unit equality using double hash
                    if ppos + unit > n or q + unit > n:
                        continue
                    h1a = subhash(h1, p1, ppos, ppos+unit)
                    h1b = subhash(h1, p1, q, q+unit)
                    h2a = subhash(h2, p2, ppos, ppos+unit)
                    h2b = subhash(h2, p2, q, q+unit)
                    if h1a == h1b and h2a == h2b:
                        # compute repeat count: see how many tandem copies (rare here because direct repeats may be non-tandem)
                        # For direct repeats we consider it a pair (ppos, q)
                        raw, norm = score_direct(unit, 2, spacer, spacer_max)
                        results.append({
                            "class":"Slipped_DNA",
                            "subclass":"Direct_Repeat",
                            "start": ppos+1, "end": ppos+unit,
                            "start2": q+1, "end2": q+unit,
                            "length": unit,
                            "sequence1": seq[ppos:ppos+unit],
                            "sequence2": seq[q:q+unit],
                            "spacer": spacer,
                            "raw_score": raw, "normalized_score": norm,
                            "strand": "+",
                        })
    return results

def detect_inverted_repeats(seq, min_arm=MIN_CRUCIFORM_ARM, max_spacer=MAX_CRUCIFORM_SPACER, seed_k=KMER_INVERT_SEED):
    n = len(seq)
    results = []
    if n < seed_k:
        return results
    # build k-mer index for forward sequence
    km_idx = build_kmer_index(seq, seed_k)
    # rolling hashes for seq and for rc(seq)
    h1, p1, h2, p2 = build_double_hashes(seq)
    rcseq = rc_str(seq)
    h1_rc, p1_rc, h2_rc, p2_rc = build_double_hashes(rcseq)
    # iterate left_end as rightmost index of left arm
    for left_end in range(min_arm-1, n):
        left_k_start = left_end - seed_k + 1
        if left_k_start < 0:
            continue
        seed = seq[left_k_start:left_k_start+seed_k]
        seed_rc = rc_str(seed)
        seed_rc_int = seq_to_int_kmer(seed_rc)
        # candidate right_start positions where k-mer equals seed_rc
        cand = km_idx.get(seed_rc_int, [])
        if not cand:
            continue
        # restrict candidates to spacer window
        rmin = left_end + 1
        rmax = min(n - min_arm, left_end + 1 + max_spacer)
        lo = bisect.bisect_left(cand, rmin)
        hi = bisect.bisect_right(cand, rmax)
        for q in cand[lo:hi]:
            # binary search maximal arm L
            loL = min_arm
            hiL = min(n - q, left_end + 1)  # arm cannot exceed boundaries
            best = 0
            while loL <= hiL:
                mid = (loL + hiL) // 2
                l1 = left_end - mid + 1
                r1 = left_end + 1
                pos_in_rc = n - (q + mid)
                if pos_in_rc < 0:
                    hiL = mid - 1
                    continue
                h1a = subhash(h1, p1, l1, r1)
                h1b = subhash(h1_rc, p1_rc, pos_in_rc, pos_in_rc + mid)
                h2a = subhash(h2, p2, l1, r1)
                h2b = subhash(h2_rc, p2_rc, pos_in_rc, pos_in_rc + mid)
                if h1a == h1b and h2a == h2b:
                    best = mid
                    loL = mid + 1
                else:
                    hiL = mid - 1
            if best >= min_arm:
                left_start = left_end - best + 1
                right_end = q + best - 1
                spacer = q - (left_end + 1)
                raw, norm = score_cruciform(best, spacer, max_spacer)
                results.append({
                    "class":"Cruciform",
                    "subclass":"Inverted_Repeat",
                    "left_start": left_start+1, "left_end": left_end+1,
                    "right_start": q+1, "right_end": right_end+1,
                    "arm_len": best, "spacer": spacer,
                    "sequence_left": seq[left_start:left_start+best],
                    "sequence_right": seq[q:q+best],
                    "raw_score": raw, "normalized_score": norm,
                    "strand":"+"
                })
    return results

def detect_mirror_triplex(seq, min_arm=MIN_TRIPLEX_ARM, max_spacer=MAX_TRIPLEX_SPACER, seed_k=KMER_INVERT_SEED, purity_thresh=TRIPLEX_PURITY):
    n = len(seq)
    results = []
    if n < seed_k:
        return results
    km_idx = build_kmer_index(seq, seed_k)
    h1, p1, h2, p2 = build_double_hashes(seq)
    revseq = rev_str(seq)
    h1_rev, p1_rev, h2_rev, p2_rev = build_double_hashes(revseq)
    for left_end in range(min_arm-1, n):
        left_k_start = left_end - seed_k + 1
        if left_k_start < 0:
            continue
        seed = seq[left_k_start:left_k_start+seed_k]
        seed_rev = rev_str(seed)
        seed_rev_int = seq_to_int_kmer(seed_rev)
        cand = km_idx.get(seed_rev_int, [])
        if not cand:
            continue
        rmin = left_end + 1
        rmax = min(n - min_arm, left_end + 1 + max_spacer)
        lo = bisect.bisect_left(cand, rmin)
        hi = bisect.bisect_right(cand, rmax)
        for q in cand[lo:hi]:
            loL = min_arm
            hiL = min(n - q, left_end + 1)
            best = 0
            while loL <= hiL:
                mid = (loL + hiL)//2
                l1 = left_end - mid + 1
                r1 = left_end + 1
                pos_in_rev = n - (q + mid)
                if pos_in_rev < 0:
                    hiL = mid - 1; continue
                h1a = subhash(h1, p1, l1, r1)
                h1b = subhash(h1_rev, p1_rev, pos_in_rev, pos_in_rev + mid)
                h2a = subhash(h2, p2, l1, r1)
                h2b = subhash(h2_rev, p2_rev, pos_in_rev, pos_in_rev + mid)
                if h1a == h1b and h2a == h2b:
                    best = mid
                    loL = mid + 1
                else:
                    hiL = mid - 1
            if best >= min_arm:
                left_start = left_end - best + 1
                right_end = q + best - 1
                spacer = q - (left_end + 1)
                arm_seq = seq[left_start:left_start+best]
                pur = sum(1 for c in arm_seq if c in 'AGag') / best
                pyr = sum(1 for c in arm_seq if c in 'CTct') / best
                purity = max(pur, pyr)
                if purity >= purity_thresh:
                    raw, norm = score_triplex(best, purity, spacer, max_spacer)
                    results.append({
                        "class":"Triplex",
                        "subclass":"Mirror_Repeat",
                        "left_start": left_start+1, "left_end": left_end+1,
                        "right_start": q+1, "right_end": right_end+1,
                        "arm_len": best, "spacer": spacer,
                        "purity": purity,
                        "sequence_arm": arm_seq,
                        "raw_score": raw, "normalized_score": norm,
                        "strand":"+"
                    })
    return results

# --------------------------
# Combined scan entrypoints
# --------------------------

def scan_sequence(seq, seq_name="sequence", verbose=False):
    """
    Run all detectors and return structured results compatible with NBDScanner format.
    """
    seq = seq.strip()
    n = len(seq)
    
    # 1) STRs
    strs = detect_STRs(seq)
    # 2) Direct repeats
    directs = detect_direct_repeats(seq)
    # 3) Cruciform/inverted
    cruciforms = detect_inverted_repeats(seq)
    # 4) Mirror/triplex
    triplexes = detect_mirror_triplex(seq)

    # Convert to NBDScanner-compatible format
    all_motifs = []
    motif_id = 1
    
    # STRs
    for r in strs:
        motif = {
            'ID': f"{seq_name}_motif_{motif_id}",
            'Sequence_Name': seq_name,
            'Class': r["class"],
            'Subclass': r["subclass"],
            'Start': r["start"],
            'End': r["end"],
            'Length': r["length"],
            'Sequence': r["sequence"],
            'Score': r["normalized_score"],
            'Strand': r["strand"],
            'Method': 'Pure_Python_STR',
            'Details': f"unit_len={r['unit_len']};copies={r['copies']}"
        }
        all_motifs.append(motif)
        motif_id += 1
    
    # Direct repeats
    for r in directs:
        motif = {
            'ID': f"{seq_name}_motif_{motif_id}",
            'Sequence_Name': seq_name,
            'Class': r["class"],
            'Subclass': r["subclass"],
            'Start': r["start"],
            'End': r["end2"],
            'Length': r["length"],
            'Sequence': r["sequence1"] + "|" + r["sequence2"],
            'Score': r["normalized_score"],
            'Strand': r["strand"],
            'Method': 'Pure_Python_Direct',
            'Details': f"unit_len={r['length']};left={r['start']}:{r['end']};right={r['start2']}:{r['end2']};spacer={r['spacer']}"
        }
        all_motifs.append(motif)
        motif_id += 1
    
    # Cruciforms
    for r in cruciforms:
        motif = {
            'ID': f"{seq_name}_motif_{motif_id}",
            'Sequence_Name': seq_name,
            'Class': r["class"],
            'Subclass': r["subclass"],
            'Start': r['left_start'],
            'End': r['right_end'],
            'Length': r['arm_len'],
            'Sequence': r["sequence_left"] + "|" + r["sequence_right"],
            'Score': r["normalized_score"],
            'Strand': r["strand"],
            'Method': 'Pure_Python_Cruciform',
            'Details': f"arm_len={r['arm_len']};left={r['left_start']}:{r['left_end']};right={r['right_start']}:{r['right_end']};spacer={r['spacer']}"
        }
        all_motifs.append(motif)
        motif_id += 1
    
    # Triplexes
    for r in triplexes:
        motif = {
            'ID': f"{seq_name}_motif_{motif_id}",
            'Sequence_Name': seq_name,
            'Class': r["class"],
            'Subclass': r["subclass"],
            'Start': r['left_start'],
            'End': r['right_end'],
            'Length': r['arm_len'],
            'Sequence': r["sequence_arm"],
            'Score': r["normalized_score"],
            'Strand': r["strand"],
            'Method': 'Pure_Python_Triplex',
            'Details': f"arm_len={r['arm_len']};left={r['left_start']}:{r['left_end']};right={r['right_start']}:{r['right_end']};spacer={r['spacer']};purity={r['purity']:.3f}"
        }
        all_motifs.append(motif)
        motif_id += 1

    if verbose:
        print(f"# Detected: STRs={len(strs)} directs={len(directs)} cruciforms={len(cruciforms)} triplexes={len(triplexes)}", file=sys.stderr)

    return all_motifs

def scan_sequence_tsv_output(seq, seq_name="sequence", verbose=False):
    """
    Run all detectors and print TSV to stdout (original format).
    """
    seq = seq.strip()
    n = len(seq)
    
    # 1) STRs
    strs = detect_STRs(seq)
    # 2) Direct repeats
    directs = detect_direct_repeats(seq)
    # 3) Cruciform/inverted
    cruciforms = detect_inverted_repeats(seq)
    # 4) Mirror/triplex
    triplexes = detect_mirror_triplex(seq)

    # assemble output rows
    # header
    header = ["motif_id","sequence_name","class","subclass","start","end","length","sequence","raw_score","normalized_score","strand","details"]
    print("\t".join(header))
    mid = 1
    # STRs
    for r in strs:
        details = f"unit_len={r['unit_len']};copies={r['copies']}"
        row = [f"motif_{mid}", seq_name, r["class"], r["subclass"], str(r["start"]), str(r["end"]),
               str(r["length"]), r["sequence"], f"{r['raw_score']:.3f}", f"{r['normalized_score']:.3f}", r["strand"], details]
        print("\t".join(row)); mid += 1
    # Direct repeats
    for r in directs:
        details = f"unit_len={r['length']};left={r['start']}:{r['end']};right={r['start2']}:{r['end2']};spacer={r['spacer']}"
        seqfrag = r["sequence1"] + "|" + r["sequence2"]
        row = [f"motif_{mid}", seq_name, "Slipped_DNA", "Direct_Repeat", str(r["start"]), str(r["end2"]), str(r["length"]),
               seqfrag, f"{r['raw_score']:.3f}", f"{r['normalized_score']:.3f}", r["strand"], details]
        print("\t".join(row)); mid += 1
    # Cruciforms
    for r in cruciforms:
        details = f"arm_len={r['arm_len']};left={r['left_start']}:{r['left_end']};right={r['right_start']}:{r['right_end']};spacer={r['spacer']}"
        seqfrag = r["sequence_left"] + "|" + r["sequence_right"]
        row = [f"motif_{mid}", seq_name, "Cruciform", "Inverted_Repeat", str(r['left_start']), str(r['right_end']), str(r['arm_len']),
               seqfrag, f"{r['raw_score']:.3f}", f"{r['normalized_score']:.3f}", r["strand"], details]
        print("\t".join(row)); mid += 1
    # Triplexes
    for r in triplexes:
        details = f"arm_len={r['arm_len']};left={r['left_start']}:{r['left_end']};right={r['right_start']}:{r['right_end']};spacer={r['spacer']};purity={r['purity']:.3f}"
        row = [f"motif_{mid}", seq_name, "Triplex", "Mirror_Repeat", str(r['left_start']), str(r['right_end']), str(r['arm_len']),
               r["sequence_arm"], f"{r['raw_score']:.3f}", f"{r['normalized_score']:.3f}", r["strand"], details]
        print("\t".join(row)); mid += 1

    if verbose:
        print(f"# Detected: STRs={len(strs)} directs={len(directs)} cruciforms={len(cruciforms)} triplexes={len(triplexes)}", file=sys.stderr)

# --------------------------
# FASTA convenience
# --------------------------
def scan_fasta(path):
    with open(path) as f:
        name = None; seq_chunks=[]
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seq = "".join(seq_chunks)
                    scan_sequence_tsv_output(seq, seq_name=name)
                name = line[1:].strip().split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if name is not None:
            seq = "".join(seq_chunks)
            scan_sequence_tsv_output(seq, seq_name=name)

# --------------------------
# If run as script
# --------------------------
if __name__ == "__main__":
    # quick demo: if a filename passed -> scan FASTA, else read stdin sequence
    if len(sys.argv) > 1:
        scan_fasta(sys.argv[1])
    else:
        seq = sys.stdin.read().strip()
        if not seq:
            print("Usage: python nonb_pure_python.py [input.fasta]  OR cat seq | python nonb_pure_python.py", file=sys.stderr)
            sys.exit(1)
        scan_sequence_tsv_output(seq, seq_name="stdin_seq", verbose=True)
