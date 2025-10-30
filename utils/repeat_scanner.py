#!/usr/bin/env python3
"""
Optimized genome-scale Python scanner for:
 - Direct repeats (unit length 10..300 bp, spacer <= 10 bp)
 - Inverted repeats (arm >= 6 bp, loop <= 100 bp)
 - Mirror repeats (arm >= 10 bp, loop <= 100 bp)
 - STRs (unit size 1..9 bp, total repeated length >= 10 bp)

Design:
 - Seed-and-extend using k-mer indices (dict: kmer -> list(positions))
 - Rolling-hash-free verification: uses direct slice compare for correctness
 - Candidate pruning: only consider seed position pairs with delta <= max_unit + max_spacer
 - Tuned defaults and safe-guards to avoid explosion on highly-repetitive seeds

Usage:
    from utils.repeat_scanner import find_direct_repeats, find_inverted_repeats, find_mirror_repeats, find_strs

Requirements:
    Python 3.8+
    (pure Python; no external C deps required)
"""

from __future__ import annotations
import math
from collections import defaultdict
from typing import List, Dict, Tuple

# -------------------------
# Parameters (user constraints)
# -------------------------
DIRECT_MIN_UNIT = 10
DIRECT_MAX_UNIT = 300
DIRECT_MAX_SPACER = 10

INVERTED_MIN_ARM = 6
INVERTED_MAX_LOOP = 100

MIRROR_MIN_ARM = 10
MIRROR_MAX_LOOP = 100

STR_MIN_UNIT = 1
STR_MAX_UNIT = 9
STR_MIN_TOTAL = 10

# Seeds (k) used for k-mer indexing (trade-off: larger k -> fewer hits)
K_DIRECT = 10   # equal to DIRECT_MIN_UNIT (good seed)
K_INVERTED = 6  # equal to INVERTED_MIN_ARM (min arm)
K_MIRROR = 10   # equal to MIRROR_MIN_ARM

# Safety thresholds
MAX_POSITIONS_PER_KMER = 10000  # skip extremely frequent k-mers to avoid blow-up

# -------------------------
# Utilities
# -------------------------
_RC_TRANS = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s: str) -> str:
    return s.translate(_RC_TRANS)[::-1]


# -------------------------
# K-mer index builder
# -------------------------
def build_kmer_index(seq: str, k: int) -> Dict[str, List[int]]:
    """
    Build dictionary mapping k-mer -> list(positions). Skips k-mers containing non-ACGT.
    """
    n = len(seq)
    idx: Dict[str, List[int]] = defaultdict(list)
    for i in range(0, n - k + 1):
        kmer = seq[i:i + k]
        if 'N' in kmer or any(ch not in "ACGT" for ch in kmer):
            continue
        lst = idx[kmer]
        if len(lst) <= MAX_POSITIONS_PER_KMER:
            lst.append(i)
    # Remove too-frequent k-mers entirely
    to_delete = [kmer for kmer, lst in idx.items() if len(lst) > MAX_POSITIONS_PER_KMER]
    for kmer in to_delete:
        del idx[kmer]
    return idx


# -------------------------
# Direct Repeats (seed-and-extend)
# -------------------------
def find_direct_repeats(seq: str, 
                       min_unit: int = DIRECT_MIN_UNIT, 
                       max_unit: int = DIRECT_MAX_UNIT,
                       max_spacer: int = DIRECT_MAX_SPACER) -> List[Dict]:
    """
    Use k-mer index (k == len(kmer) >= min_unit) to find direct repeats:
    - For each kmer positions list, iterate pairs (i, j) with j > i and delta <= max_unit + max_spacer
    - For each pair compute delta = j - i; for s in 0..max_spacer, L = delta - s; if L in [min_unit,max_unit]
      verify seq[i:i+L] == seq[j:j+L] and record.
    """
    n = len(seq)
    idx = build_kmer_index(seq, K_DIRECT)
    results = []
    max_delta = max_unit + max_spacer
    
    for kmer, poses in idx.items():
        poses_sorted = poses  # already in increasing order
        m = len(poses_sorted)
        if m < 2:
            continue
        # for each pair (i_pos, j_pos) with proximity constraint
        for a in range(m):
            i = poses_sorted[a]
            b = a + 1
            while b < m:
                j = poses_sorted[b]
                delta = j - i
                if delta > max_delta:
                    break  # further j will only be larger
                # for each possible spacer s
                # s = delta - L  => L = delta - s
                # s in [0, max_spacer] -> L in [delta - max_spacer, delta]
                L_min = max(min_unit, delta - max_spacer)
                L_max = min(max_unit, delta)
                if L_min <= L_max:
                    # test candidate lengths L in this small interval; prefer larger L first (maximal)
                    for L in range(L_max, L_min - 1, -1):
                        j_start = j
                        if j_start + L > n or i + L > n:
                            continue
                        if seq[i:i + L] == seq[j_start:j_start + L]:
                            rec = {
                                'Class': 'Direct_Repeat',
                                'Subclass': f'Direct_L{L}',
                                'Start': i + 1,
                                'End': j_start + L,
                                'Length': (j_start + L) - i,
                                'Unit_Length': L,
                                'Spacer': j_start - (i + L),
                                'Left_Pos': i + 1,
                                'Right_Pos': j_start + 1,
                                'Unit_Seq': seq[i:i + L],
                                'Sequence': seq[i:j_start + L]
                            }
                            results.append(rec)
                            break
                b += 1
    # dedupe: prefer maximal unit_length per (Left_Pos, Spacer)
    dedup: Dict[Tuple[int, int], Dict] = {}
    for rec in results:
        key = (rec['Left_Pos'], rec['Spacer'])
        if key not in dedup or rec['Unit_Length'] > dedup[key]['Unit_Length']:
            dedup[key] = rec
    return list(dedup.values())


# -------------------------
# Inverted Repeats (Cruciform)
# -------------------------
def find_inverted_repeats(seq: str, 
                         min_arm: int = INVERTED_MIN_ARM, 
                         max_loop: int = INVERTED_MAX_LOOP) -> List[Dict]:
    """
    For each k-mer (of size >= min_arm) and its reverse complement locations:
    iterate pairs (i in idx[kmer], j in idx[rc_kmer]) with j > i, delta = j - i,
    possible arm lengths: arm = delta - loop for loop in 0..max_loop, require arm >= min_arm.
    Verify equality: seq[i:i+arm] == revcomp(seq[j:j+arm])
    """
    n = len(seq)
    idx = build_kmer_index(seq, K_INVERTED)
    results = []
    
    # Precompute reverse complement mapping of keys present
    keys = list(idx.keys())
    rc_map = {}
    for kmer in keys:
        rc = revcomp(kmer)
        if rc in idx:
            rc_map[kmer] = idx[rc]
    
    # For each kmer and its rc positions, test pairs
    for kmer, rc_positions in rc_map.items():
        left_positions = idx[kmer]
        # iterate positions pairs
        for i in left_positions:
            for j in rc_positions:
                if j <= i:
                    continue
                delta = j - i
                # arm = delta - loop must be >= min_arm and <= available space
                # loop in 0..max_loop -> arm in [delta - max_loop, delta]
                arm_min = max(min_arm, delta - max_loop)
                arm_max = delta
                # cap arm_max by sequence boundary
                arm_max = min(arm_max, n - j, n - i)
                if arm_min > arm_max:
                    continue
                # try maximal arms first
                for arm in range(arm_max, arm_min - 1, -1):
                    if i + arm > n or j + arm > n:
                        continue
                    left_sub = seq[i:i + arm]
                    right_sub = seq[j:j + arm]
                    if left_sub == revcomp(right_sub):
                        rec = {
                            'Class': 'Inverted_Repeat',
                            'Subclass': f'Inverted_arm_{arm}',
                            'Start': i + 1,
                            'End': j + arm,
                            'Length': (j + arm) - i,
                            'Left_Start': i + 1,
                            'Right_Start': j + 1,
                            'Arm_Length': arm,
                            'Loop': delta - arm,
                            'Left_Arm': left_sub,
                            'Right_Arm': right_sub,
                            'Sequence': seq[i:j + arm]
                        }
                        results.append(rec)
                        break
    
    # dedupe: keep maximal arm per (Left_Start, Loop)
    dedup: Dict[Tuple[int, int], Dict] = {}
    for rec in results:
        key = (rec['Left_Start'], rec['Loop'])
        if key not in dedup or rec['Arm_Length'] > dedup[key]['Arm_Length']:
            dedup[key] = rec
    return list(dedup.values())


# -------------------------
# Mirror Repeats (Triplex DNA component)
# -------------------------
def find_mirror_repeats(seq: str, 
                       min_arm: int = MIRROR_MIN_ARM, 
                       max_loop: int = MIRROR_MAX_LOOP,
                       purine_pyrimidine_threshold: float = 0.9) -> List[Dict]:
    """
    Mirror repeats: left arm matches reverse (not complement) of right arm.
    Same pattern as inverted but compare seq[i:i+arm] == seq[j:j+arm][::-1]
    
    For Triplex DNA, we also filter for >90% purine or pyrimidine content in arms.
    """
    n = len(seq)
    idx = build_kmer_index(seq, K_MIRROR)
    results = []
    
    keys = list(idx.keys())
    # For mirror, the reverse of substring is simply reversed string
    rev_map = {}
    for kmer in keys:
        rev = kmer[::-1]
        if rev in idx:
            rev_map[kmer] = idx[rev]
    
    for kmer, rev_positions in rev_map.items():
        left_positions = idx[kmer]
        for i in left_positions:
            for j in rev_positions:
                if j <= i:
                    continue
                delta = j - i
                arm_min = max(min_arm, delta - max_loop)
                arm_max = min(delta, n - j, n - i)
                if arm_min > arm_max:
                    continue
                for arm in range(arm_max, arm_min - 1, -1):
                    if i + arm > n or j + arm > n:
                        continue
                    left_arm = seq[i:i + arm]
                    right_arm = seq[j:j + arm]
                    if left_arm == right_arm[::-1]:
                        # Check purine/pyrimidine content for Triplex DNA
                        combined_arms = left_arm + right_arm
                        purine_count = sum(1 for b in combined_arms if b in 'AG')
                        pyrimidine_count = sum(1 for b in combined_arms if b in 'CT')
                        total_bases = len(combined_arms)
                        
                        purine_fraction = purine_count / total_bases if total_bases > 0 else 0
                        pyrimidine_fraction = pyrimidine_count / total_bases if total_bases > 0 else 0
                        
                        is_triplex = (purine_fraction >= purine_pyrimidine_threshold or 
                                     pyrimidine_fraction >= purine_pyrimidine_threshold)
                        
                        rec = {
                            'Class': 'Mirror_Repeat',
                            'Subclass': f'Mirror_arm_{arm}',
                            'Start': i + 1,
                            'End': j + arm,
                            'Length': (j + arm) - i,
                            'Left_Start': i + 1,
                            'Right_Start': j + 1,
                            'Arm_Length': arm,
                            'Loop': delta - arm,
                            'Left_Arm': left_arm,
                            'Right_Arm': right_arm,
                            'Sequence': seq[i:j + arm],
                            'Is_Triplex': is_triplex,
                            'Purine_Fraction': round(purine_fraction, 3),
                            'Pyrimidine_Fraction': round(pyrimidine_fraction, 3)
                        }
                        results.append(rec)
                        break
    
    # dedupe maximal arms
    dedup: Dict[Tuple[int, int], Dict] = {}
    for rec in results:
        key = (rec['Left_Start'], rec['Loop'])
        if key not in dedup or rec['Arm_Length'] > dedup[key]['Arm_Length']:
            dedup[key] = rec
    return list(dedup.values())


# -------------------------
# STRs (Short Tandem Repeats)
# -------------------------
def find_strs(seq: str, 
              min_u: int = STR_MIN_UNIT, 
              max_u: int = STR_MAX_UNIT, 
              min_total: int = STR_MIN_TOTAL) -> List[Dict]:
    """
    Greedy detection of perfect STRs (tandem repeats).
    For each unit size k in 1..9, slide and count consecutive copies.
    """
    n = len(seq)
    results = []
    for k in range(min_u, max_u + 1):
        i = 0
        while i <= n - k:
            unit = seq[i:i + k]
            # require at least one repeat following to qualify quickly
            if i + k >= n or seq[i + k:i + 2 * k] != unit:
                i += 1
                continue
            # count copies
            copies = 1
            j = i + k
            while j + k <= n and seq[j:j + k] == unit:
                copies += 1
                j += k
            total_len = copies * k
            if total_len >= min_total:
                results.append({
                    'Class': 'STR',
                    'Subclass': f'unit_{k}',
                    'Start': i + 1,
                    'End': j,
                    'Unit_Length': k,
                    'Copies': copies,
                    'Length': total_len,
                    'Unit_Seq': unit,
                    'Sequence': seq[i:j]
                })
                i = j  # skip to end of tandem block
            else:
                i += 1
    return results
