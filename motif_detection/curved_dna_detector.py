"""
CurvedDNADetector

Detects:
 - Global curvature (A-phased repeats, APRs): >= 3 A-tract centers phased ~11 bp apart
 - Local curvature: long A-tracts (>=7) or T-tracts (>=7)

Implements A-tract detection logic similar to the provided C code: within AT-rich windows,
computes the longest A/AnTn run and the longest T-only run, uses difference (maxATlen - maxTlen)
to decide bona fide A-tracts and reports the tract center.
"""

import re
from typing import List, Dict, Any, Tuple
from .base_detector import BaseMotifDetector


def revcomp(seq: str) -> str:
    trans = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(trans)[::-1]


class CurvedDNADetector(BaseMotifDetector):
    def get_motif_class_name(self) -> str:
        return "Curved_DNA"

    # ---------- Parameters you can tune ----------
    MIN_AT_TRACT = 3         # minimum A-tract length (for global APR detection)
    MAX_AT_WINDOW = None     # None => no hard upper limit on AT window used to search AnTn patterns
    PHASING_CENTER_SPACING = 11.0  # ideal center-to-center spacing in bp for APR phasing
    PHASING_TOL_LOW = 9.9    # lower tolerance
    PHASING_TOL_HIGH = 11.1  # upper tolerance
    MIN_APR_TRACTS = 3       # at least this many A-tract centers to call an APR
    LOCAL_LONG_TRACT = 7     # local curvature: A>=7 or T>=7
    # --------------------------------------------

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """
        Comprehensive curved DNA patterns including:
        - Local curvature: long A/T tracts (>=7)
        - Global curvature: A-phased and T-phased repeats (APRs)
        Based on problem statement specifications.
        """
        return {
            # Local curvature patterns
            'local_curved': [
                (r'A{7,}', 'CRV_002', 'Long A-tract', 'Local Curvature', 7, 'curvature_score', 0.95, 'A-tract curvature', 'Olson 1998'),
                (r'T{7,}', 'CRV_003', 'Long T-tract', 'Local Curvature', 7, 'curvature_score', 0.95, 'T-tract curvature', 'Olson 1998'),
            ],
            
            # Global curvature: 3-tract A-phased repeats (APRs)
            'global_curved_a_3tract': [
                (r'(?:A{3}[ACGT]{6,8}A{3}[ACGT]{6,8}A{3})', 'CRV_008', 'A3-APR', 'Global Curvature', 20, 'phasing_score', 0.90, '3-tract APR', 'Koo 1986'),
                (r'(?:A{4}[ACGT]{5,7}A{4}[ACGT]{5,7}A{4})', 'CRV_009', 'A4-APR', 'Global Curvature', 20, 'phasing_score', 0.90, '3-tract APR', 'Koo 1986'),
                (r'(?:A{5}[ACGT]{4,6}A{5}[ACGT]{4,6}A{5})', 'CRV_010', 'A5-APR', 'Global Curvature', 20, 'phasing_score', 0.90, '3-tract APR', 'Koo 1986'),
                (r'(?:A{6}[ACGT]{3,5}A{6}[ACGT]{3,5}A{6})', 'CRV_011', 'A6-APR', 'Global Curvature', 20, 'phasing_score', 0.90, '3-tract APR', 'Koo 1986'),
                (r'(?:A{7}[ACGT]{2,4}A{7}[ACGT]{2,4}A{7})', 'CRV_012', 'A7-APR', 'Global Curvature', 20, 'phasing_score', 0.90, '3-tract APR', 'Koo 1986'),
                (r'(?:A{8}[ACGT]{1,3}A{8}[ACGT]{1,3}A{8})', 'CRV_013', 'A8-APR', 'Global Curvature', 20, 'phasing_score', 0.90, '3-tract APR', 'Koo 1986'),
                (r'(?:A{9}[ACGT]{0,2}A{9}[ACGT]{0,2}A{9})', 'CRV_014', 'A9-APR', 'Global Curvature', 20, 'phasing_score', 0.90, '3-tract APR', 'Koo 1986'),
            ],
            
            # Global curvature: 4-tract A-phased repeats
            'global_curved_a_4tract': [
                (r'(?:A{3}[ACGT]{6,8}A{3}[ACGT]{6,8}A{3}[ACGT]{6,8}A{3})', 'CRV_015', 'A3-APR-4', 'Global Curvature', 25, 'phasing_score', 0.92, '4-tract APR', 'Koo 1986'),
                (r'(?:A{4}[ACGT]{5,7}A{4}[ACGT]{5,7}A{4}[ACGT]{5,7}A{4})', 'CRV_016', 'A4-APR-4', 'Global Curvature', 25, 'phasing_score', 0.92, '4-tract APR', 'Koo 1986'),
                (r'(?:A{5}[ACGT]{4,6}A{5}[ACGT]{4,6}A{5}[ACGT]{4,6}A{5})', 'CRV_017', 'A5-APR-4', 'Global Curvature', 25, 'phasing_score', 0.92, '4-tract APR', 'Koo 1986'),
                (r'(?:A{6}[ACGT]{3,5}A{6}[ACGT]{3,5}A{6}[ACGT]{3,5}A{6})', 'CRV_018', 'A6-APR-4', 'Global Curvature', 25, 'phasing_score', 0.92, '4-tract APR', 'Koo 1986'),
                (r'(?:A{7}[ACGT]{2,4}A{7}[ACGT]{2,4}A{7}[ACGT]{2,4}A{7})', 'CRV_019', 'A7-APR-4', 'Global Curvature', 25, 'phasing_score', 0.92, '4-tract APR', 'Koo 1986'),
                (r'(?:A{8}[ACGT]{1,3}A{8}[ACGT]{1,3}A{8}[ACGT]{1,3}A{8})', 'CRV_020', 'A8-APR-4', 'Global Curvature', 25, 'phasing_score', 0.92, '4-tract APR', 'Koo 1986'),
                (r'(?:A{9}[ACGT]{0,2}A{9}[ACGT]{0,2}A{9}[ACGT]{0,2}A{9})', 'CRV_021', 'A9-APR-4', 'Global Curvature', 25, 'phasing_score', 0.92, '4-tract APR', 'Koo 1986'),
            ],
            
            # Global curvature: 5-tract A-phased repeats
            'global_curved_a_5tract': [
                (r'(?:A{3}[ACGT]{6,8}A{3}[ACGT]{6,8}A{3}[ACGT]{6,8}A{3}[ACGT]{6,8}A{3})', 'CRV_022', 'A3-APR-5', 'Global Curvature', 30, 'phasing_score', 0.95, '5-tract APR', 'Koo 1986'),
                (r'(?:A{4}[ACGT]{5,7}A{4}[ACGT]{5,7}A{4}[ACGT]{5,7}A{4}[ACGT]{5,7}A{4})', 'CRV_023', 'A4-APR-5', 'Global Curvature', 30, 'phasing_score', 0.95, '5-tract APR', 'Koo 1986'),
                (r'(?:A{5}[ACGT]{4,6}A{5}[ACGT]{4,6}A{5}[ACGT]{4,6}A{5}[ACGT]{4,6}A{5})', 'CRV_024', 'A5-APR-5', 'Global Curvature', 30, 'phasing_score', 0.95, '5-tract APR', 'Koo 1986'),
                (r'(?:A{6}[ACGT]{3,5}A{6}[ACGT]{3,5}A{6}[ACGT]{3,5}A{6}[ACGT]{3,5}A{6})', 'CRV_025', 'A6-APR-5', 'Global Curvature', 30, 'phasing_score', 0.95, '5-tract APR', 'Koo 1986'),
                (r'(?:A{7}[ACGT]{2,4}A{7}[ACGT]{2,4}A{7}[ACGT]{2,4}A{7}[ACGT]{2,4}A{7})', 'CRV_026', 'A7-APR-5', 'Global Curvature', 30, 'phasing_score', 0.95, '5-tract APR', 'Koo 1986'),
                (r'(?:A{8}[ACGT]{1,3}A{8}[ACGT]{1,3}A{8}[ACGT]{1,3}A{8}[ACGT]{1,3}A{8})', 'CRV_027', 'A8-APR-5', 'Global Curvature', 30, 'phasing_score', 0.95, '5-tract APR', 'Koo 1986'),
                (r'(?:A{9}[ACGT]{0,2}A{9}[ACGT]{0,2}A{9}[ACGT]{0,2}A{9}[ACGT]{0,2}A{9})', 'CRV_028', 'A9-APR-5', 'Global Curvature', 30, 'phasing_score', 0.95, '5-tract APR', 'Koo 1986'),
            ],
            
            # Global curvature: 3-tract T-phased repeats
            'global_curved_t_3tract': [
                (r'(?:T{3}[ACGT]{6,8}T{3}[ACGT]{6,8}T{3})', 'CRV_029', 'T3-TPR', 'Global Curvature', 20, 'phasing_score', 0.90, '3-tract TPR', 'Koo 1986'),
                (r'(?:T{4}[ACGT]{5,7}T{4}[ACGT]{5,7}T{4})', 'CRV_030', 'T4-TPR', 'Global Curvature', 20, 'phasing_score', 0.90, '3-tract TPR', 'Koo 1986'),
                (r'(?:T{5}[ACGT]{4,6}T{5}[ACGT]{4,6}T{5})', 'CRV_031', 'T5-TPR', 'Global Curvature', 20, 'phasing_score', 0.90, '3-tract TPR', 'Koo 1986'),
                (r'(?:T{6}[ACGT]{3,5}T{6}[ACGT]{3,5}T{6})', 'CRV_032', 'T6-TPR', 'Global Curvature', 20, 'phasing_score', 0.90, '3-tract TPR', 'Koo 1986'),
                (r'(?:T{7}[ACGT]{2,4}T{7}[ACGT]{2,4}T{7})', 'CRV_033', 'T7-TPR', 'Global Curvature', 20, 'phasing_score', 0.90, '3-tract TPR', 'Koo 1986'),
                (r'(?:T{8}[ACGT]{1,3}T{8}[ACGT]{1,3}T{8})', 'CRV_034', 'T8-TPR', 'Global Curvature', 20, 'phasing_score', 0.90, '3-tract TPR', 'Koo 1986'),
                (r'(?:T{9}[ACGT]{0,2}T{9}[ACGT]{0,2}T{9})', 'CRV_035', 'T9-TPR', 'Global Curvature', 20, 'phasing_score', 0.90, '3-tract TPR', 'Koo 1986'),
            ],
            
            # Global curvature: 4-tract T-phased repeats
            'global_curved_t_4tract': [
                (r'(?:T{3}[ACGT]{6,8}T{3}[ACGT]{6,8}T{3}[ACGT]{6,8}T{3})', 'CRV_036', 'T3-TPR-4', 'Global Curvature', 25, 'phasing_score', 0.92, '4-tract TPR', 'Koo 1986'),
                (r'(?:T{4}[ACGT]{5,7}T{4}[ACGT]{5,7}T{4}[ACGT]{5,7}T{4})', 'CRV_037', 'T4-TPR-4', 'Global Curvature', 25, 'phasing_score', 0.92, '4-tract TPR', 'Koo 1986'),
                (r'(?:T{5}[ACGT]{4,6}T{5}[ACGT]{4,6}T{5}[ACGT]{4,6}T{5})', 'CRV_038', 'T5-TPR-4', 'Global Curvature', 25, 'phasing_score', 0.92, '4-tract TPR', 'Koo 1986'),
                (r'(?:T{6}[ACGT]{3,5}T{6}[ACGT]{3,5}T{6}[ACGT]{3,5}T{6})', 'CRV_039', 'T6-TPR-4', 'Global Curvature', 25, 'phasing_score', 0.92, '4-tract TPR', 'Koo 1986'),
                (r'(?:T{7}[ACGT]{2,4}T{7}[ACGT]{2,4}T{7}[ACGT]{2,4}T{7})', 'CRV_040', 'T7-TPR-4', 'Global Curvature', 25, 'phasing_score', 0.92, '4-tract TPR', 'Koo 1986'),
                (r'(?:T{8}[ACGT]{1,3}T{8}[ACGT]{1,3}T{8}[ACGT]{1,3}T{8})', 'CRV_041', 'T8-TPR-4', 'Global Curvature', 25, 'phasing_score', 0.92, '4-tract TPR', 'Koo 1986'),
                (r'(?:T{9}[ACGT]{0,2}T{9}[ACGT]{0,2}T{9}[ACGT]{0,2}T{9})', 'CRV_042', 'T9-TPR-4', 'Global Curvature', 25, 'phasing_score', 0.92, '4-tract TPR', 'Koo 1986'),
            ],
            
            # Global curvature: 5-tract T-phased repeats
            'global_curved_t_5tract': [
                (r'(?:T{3}[ACGT]{6,8}T{3}[ACGT]{6,8}T{3}[ACGT]{6,8}T{3}[ACGT]{6,8}T{3})', 'CRV_043', 'T3-TPR-5', 'Global Curvature', 30, 'phasing_score', 0.95, '5-tract TPR', 'Koo 1986'),
                (r'(?:T{4}[ACGT]{5,7}T{4}[ACGT]{5,7}T{4}[ACGT]{5,7}T{4}[ACGT]{5,7}T{4})', 'CRV_044', 'T4-TPR-5', 'Global Curvature', 30, 'phasing_score', 0.95, '5-tract TPR', 'Koo 1986'),
                (r'(?:T{5}[ACGT]{4,6}T{5}[ACGT]{4,6}T{5}[ACGT]{4,6}T{5}[ACGT]{4,6}T{5})', 'CRV_045', 'T5-TPR-5', 'Global Curvature', 30, 'phasing_score', 0.95, '5-tract TPR', 'Koo 1986'),
                (r'(?:T{6}[ACGT]{3,5}T{6}[ACGT]{3,5}T{6}[ACGT]{3,5}T{6}[ACGT]{3,5}T{6})', 'CRV_046', 'T6-TPR-5', 'Global Curvature', 30, 'phasing_score', 0.95, '5-tract TPR', 'Koo 1986'),
                (r'(?:T{7}[ACGT]{2,4}T{7}[ACGT]{2,4}T{7}[ACGT]{2,4}T{7}[ACGT]{2,4}T{7})', 'CRV_047', 'T7-TPR-5', 'Global Curvature', 30, 'phasing_score', 0.95, '5-tract TPR', 'Koo 1986'),
                (r'(?:T{8}[ACGT]{1,3}T{8}[ACGT]{1,3}T{8}[ACGT]{1,3}T{8}[ACGT]{1,3}T{8})', 'CRV_048', 'T8-TPR-5', 'Global Curvature', 30, 'phasing_score', 0.95, '5-tract TPR', 'Koo 1986'),
                (r'(?:T{9}[ACGT]{0,2}T{9}[ACGT]{0,2}T{9}[ACGT]{0,2}T{9}[ACGT]{0,2}T{9})', 'CRV_049', 'T9-TPR-5', 'Global Curvature', 30, 'phasing_score', 0.95, '5-tract TPR', 'Koo 1986'),
            ]
        }

    def _remove_overlaps(self, motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Remove overlapping motifs within the same subclass.
        Allows overlaps between different subclasses (Global vs Local curvature).
        """
        if not motifs:
            return []
        
        from collections import defaultdict
        
        # Group by subclass
        groups = defaultdict(list)
        for motif in motifs:
            subclass = motif.get('Subclass', 'unknown')
            groups[subclass].append(motif)
        
        non_overlapping = []
        
        # Process each subclass separately
        for subclass, group_motifs in groups.items():
            # Sort by score (descending), then by length (descending)
            sorted_motifs = sorted(group_motifs, 
                                  key=lambda x: (-x.get('Score', 0), -x.get('Length', 0)))
            
            selected = []
            for motif in sorted_motifs:
                # Check if this motif overlaps with any already selected in this subclass
                overlaps = False
                for selected_motif in selected:
                    if not (motif['End'] <= selected_motif['Start'] or 
                           motif['Start'] >= selected_motif['End']):
                        overlaps = True
                        break
                
                if not overlaps:
                    selected.append(motif)
            
            non_overlapping.extend(selected)
        
        # Sort by start position for output
        non_overlapping.sort(key=lambda x: x['Start'])
        return non_overlapping

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Override base method to use sophisticated curved DNA detection"""
        sequence = sequence.upper().strip()
        motifs = []
        
        # Use the sophisticated annotation method
        annotation = self.annotate_sequence(sequence)
        
        # Extract APR (A-phased repeat) motifs
        for i, apr in enumerate(annotation.get('aprs', [])):
            if apr.get('score', 0) > 0.1:  # Lower threshold for sensitivity
                start_pos = int(min(apr['center_positions'])) - 10  # Estimate start
                end_pos = int(max(apr['center_positions'])) + 10    # Estimate end
                start_pos = max(0, start_pos)
                end_pos = min(len(sequence), end_pos)
                
                motifs.append({
                    'ID': f"{sequence_name}_CRV_APR_{start_pos+1}",
                    'Sequence_Name': sequence_name,
                    'Class': self.get_motif_class_name(),
                    'Subclass': 'Global Curvature',
                    'Start': start_pos + 1,  # 1-based coordinates
                    'End': end_pos,
                    'Length': end_pos - start_pos,
                    'Sequence': sequence[start_pos:end_pos],
                    'Score': round(apr.get('score', 0), 3),
                    'Strand': '+',
                    'Method': 'Curved_DNA_detection',
                    'Pattern_ID': f'CRV_APR_{i+1}'
                })
        
        # Extract long tract motifs
        for i, tract in enumerate(annotation.get('long_tracts', [])):
            if tract.get('score', 0) > 0.1:  # Lower threshold for sensitivity
                start_pos = tract['start']
                end_pos = tract['end']
                
                motifs.append({
                    'ID': f"{sequence_name}_CRV_TRACT_{start_pos+1}",
                    'Sequence_Name': sequence_name,
                    'Class': self.get_motif_class_name(),
                    'Subclass': 'Local Curvature',
                    'Start': start_pos + 1,  # 1-based coordinates
                    'End': end_pos,
                    'Length': end_pos - start_pos,
                    'Sequence': sequence[start_pos:end_pos],
                    'Score': round(tract.get('score', 0), 3),
                    'Strand': '+',
                    'Method': 'Curved_DNA_detection',
                    'Pattern_ID': f'CRV_TRACT_{i+1}'
                })
        
        # Remove overlaps within each subclass
        motifs = self._remove_overlaps(motifs)
        
        return motifs

    # -------------------------
    # Top-level scoring API
    # -------------------------
    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
        """
        Returns a combined raw score reflecting:
          - phasing_score for APRs (sum of APR phasing scores)
          - local curvature contribution (sum of local A/T tract scores)
        The sum reflects both number and quality of hits.
        """
        seq = sequence.upper()
        ann = self.annotate_sequence(seq)
        # Sum APR scores
        apr_sum = sum(a['score'] for a in ann.get('aprs', []))
        local_sum = sum(l['score'] for l in ann.get('long_tracts', []))
        return float(apr_sum + local_sum)

    # -------------------------
    # A-tract detection (core)
    # -------------------------
    def find_a_tracts(self, sequence: str, minAT: int = None, max_window: int = None) -> List[Dict[str, Any]]:
        """
        Detect A-tract candidates across the sequence using logic adapted from your C code.

        Returns list of dicts:
          {
            'start', 'end'            : region bounds of the AT-window inspected (0-based, end-exclusive)
            'maxATlen'                : maximal A/AnTn length found inside window
            'maxTlen'                 : maximal T-only run length (not following A)
            'maxATend'                : index (0-based) of end position(where maxATlen ends) relative to full seq (inclusive end index)
            'a_center'                : float center coordinate (1-based in C code; here 0-based float)
            'call'                    : bool whether (maxATlen - maxTlen) >= minAT
            'window_len'              : length of AT window
            'window_seq'              : the AT-window substring
          }

        Implementation detail:
        - We scan for contiguous runs of A/T (AT windows). Within each, we iterate positions and
          compute Alen/Tlen/ATlen per the C algorithm to determine maxATlen and maxTlen.
        - If either forward strand or reverse complement has (maxATlen - maxTlen) >= minAT, we call it an A-tract.
        """
        seq = sequence.upper()
        n = len(seq)
        if minAT is None:
            minAT = self.MIN_AT_TRACT
        if max_window is None:
            max_window = self.MAX_AT_WINDOW  # None allowed

        results: List[Dict[str, Any]] = []

        # Find contiguous A/T windows (length >= minAT)
        for m in re.finditer(r'[AT]{' + str(minAT) + r',}', seq):
            wstart, wend = m.start(), m.end()  # [wstart, wend)
            window_seq = seq[wstart:wend]
            window_len = wend - wstart

            # analyze forward strand window
            maxATlen, maxATend, maxTlen = self._analyze_at_window(window_seq)
            # analyze reverse complement window (to mimic C code check on reverse)
            rc_window = revcomp(window_seq)
            maxATlen_rc, maxATend_rc, maxTlen_rc = self._analyze_at_window(rc_window)

            # compute decisions - apply same logic as C code:
            diff_forward = maxATlen - maxTlen
            diff_rc = maxATlen_rc - maxTlen_rc
            call = False
            chosen_center = None
            chosen_maxATlen = None

            if diff_forward >= minAT or diff_rc >= minAT:
                call = True
                # choose the strand giving larger difference
                if diff_forward >= diff_rc:
                    chosen_maxATlen = maxATlen
                    # compute center coordinate in full sequence (0-based float center)
                    # in C code: a_center = maxATend - ((maxATlen-1)/2) + 1  (1-based)
                    # we'll produce 0-based center = (wstart + maxATend - ((maxATlen-1)/2))
                    chosen_center = (wstart + maxATend) - ((maxATlen - 1) / 2.0)
                else:
                    chosen_maxATlen = maxATlen_rc
                    # maxATend_rc is position in RC sequence; convert to original coords:
                    # RC index i corresponds to original index: wstart + (window_len - 1 - i)
                    # maxATend_rc is index in window_rc (end position index)
                    # In C, they convert similarly; for simplicity compute center via rc mapping
                    i_rc = maxATend_rc
                    # rc_end_original = wstart + (window_len - 1 - i_rc)
                    rc_end_original = wstart + (window_len - 1 - maxATend_rc)
                    chosen_center = rc_end_original - ((chosen_maxATlen - 1) / 2.0)

            results.append({
                'start': wstart,
                'end': wend,
                'window_len': window_len,
                'window_seq': window_seq,
                'maxATlen': int(maxATlen),
                'maxATend': int(wstart + maxATend),
                'maxTlen': int(maxTlen),
                'maxATlen_rc': int(maxATlen_rc),
                'maxATend_rc': int(wstart + (window_len - 1 - maxATend_rc)),
                'maxTlen_rc': int(maxTlen_rc),
                'diff_forward': int(diff_forward),
                'diff_rc': int(diff_rc),
                'call': bool(call),
                'a_center': float(chosen_center) if chosen_center is not None else None,
                'chosen_maxATlen': int(chosen_maxATlen) if chosen_maxATlen is not None else None
            })

        return results

    def _analyze_at_window(self, window_seq: str) -> Tuple[int,int,int]:
        """
        Analyze a contiguous A/T window and return (maxATlen, maxATend_index_in_window, maxTlen)
        Implemented following the logic in your C code:
         - iterate positions; update Alen, Tlen, ATlen, TAlen; track maxATlen, maxTlen and their end positions.
         - maxATend returned as index (0-based) *within the window* of the last position of the max AT run.
        """
        Alen = 0
        Tlen = 0
        ATlen = 0
        TAlen = 0
        maxATlen = 0
        maxTlen = 0
        maxATend = 0
        maxTend = 0
        # we'll iterate from index 0..len(window_seq)-1
        L = len(window_seq)
        # to mimic C code scanning with lookbacks, we iterate straightforwardly
        for i in range(L):
            ch = window_seq[i]
            prev = window_seq[i-1] if i>0 else None
            if ch == 'A':
                Tlen = 0
                TAlen = 0
                # if previous base was T, reset A-run counters per C code
                if prev == 'T':
                    Alen = 1
                    ATlen = 1
                else:
                    Alen += 1
                    ATlen += 1
            elif ch == 'T':
                # if T follows A-run shorter than Alen, it's considered TAlen (T following A)
                if TAlen < Alen:
                    TAlen += 1
                    ATlen += 1
                else:
                    # T is starting a T-only run
                    Tlen += 1
                    TAlen = 0
                    ATlen = 0
                    Alen = 0
            else:
                # non-AT not expected inside this window (we only pass contiguous AT windows)
                Alen = 0
                Tlen = 0
                ATlen = 0
                TAlen = 0
            if ATlen > maxATlen:
                maxATlen = ATlen
                maxATend = i  # end index within window
            if Tlen > maxTlen:
                maxTlen = Tlen
                maxTend = i
        return int(maxATlen), int(maxATend), int(maxTlen)

    # -------------------------
    # APR grouping / phasing
    # -------------------------
    def find_aprs(self, sequence: str, min_tract: int = None, min_apr_tracts: int = None) -> List[Dict[str, Any]]:
        """
        Group a-tract centers into APRs (A-phased repeats).
        Criteria:
          - at least min_apr_tracts centers
          - consecutive center-to-center spacing must be within PHASING_TOL_LOW..PHASING_TOL_HIGH
            (we allow flexible grouping: we slide through centers looking for runs of centers that satisfy spacing)
        Returns list of dicts:
          { 'start_center_idx', 'end_center_idx', 'centers': [...], 'center_positions': [...], 'score': phasing_score, 'n_tracts' }
        """
        if min_tract is None:
            min_tract = self.MIN_AT_TRACT
        if min_apr_tracts is None:
            min_apr_tracts = self.MIN_APR_TRACTS

        # get a-tract calls
        a_calls = [r for r in self.find_a_tracts(sequence, minAT=min_tract) if r['call'] and r['a_center'] is not None]
        centers = [r['a_center'] for r in a_calls]
        centers_sorted = sorted(centers)
        aprs: List[Dict[str, Any]] = []

        if len(centers_sorted) < min_apr_tracts:
            return aprs

        # find runs of centers where consecutive spacing is within tolerance
        i = 0
        while i < len(centers_sorted):
            run = [centers_sorted[i]]
            j = i + 1
            while j < len(centers_sorted):
                spacing = centers_sorted[j] - centers_sorted[j-1]
                if self.PHASING_TOL_LOW <= spacing <= self.PHASING_TOL_HIGH:
                    run.append(centers_sorted[j])
                    j += 1
                else:
                    break
            # if run has enough tracts, call APR
            if len(run) >= min_apr_tracts:
                # score APR by how close spacings are to ideal spacing
                spacings = [run[k+1] - run[k] for k in range(len(run)-1)]
                # closeness = product of gaussian-like terms, but simpler: average deviation
                devs = [abs(sp - self.PHASING_CENTER_SPACING) for sp in spacings]
                # normalized closeness
                mean_dev = sum(devs) / len(devs) if devs else 0.0
                # phasing_score between 0..1: 1 when mean_dev==0, drop linearly with dev up to tolerance
                max_dev_allowed = max(abs(self.PHASING_TOL_HIGH - self.PHASING_CENTER_SPACING),
                                      abs(self.PHASING_CENTER_SPACING - self.PHASING_TOL_LOW))
                phasing_score = max(0.0, 1.0 - (mean_dev / (max_dev_allowed if max_dev_allowed>0 else 1.0)))
                aprs.append({
                    'start_center_idx': i,
                    'end_center_idx': j-1,
                    'center_positions': run,
                    'n_tracts': len(run),
                    'spacings': spacings,
                    'mean_deviation': mean_dev,
                    'score': round(phasing_score, 6)
                })
            i = j

        return aprs

    # -------------------------
    # Local long tract finder
    # -------------------------
    def find_long_tracts(self, sequence: str, min_len: int = None) -> List[Dict[str, Any]]:
        """
        Finds long A-tracts or T-tracts with length >= min_len (default LOCAL_LONG_TRACT).
        Returns list of dicts: {start,end,base,len,score} with score derived from len.
        """
        if min_len is None:
            min_len = self.LOCAL_LONG_TRACT
        seq = sequence.upper()
        results = []
        # A runs
        for m in re.finditer(r'A{' + str(min_len) + r',}', seq):
            ln = m.end() - m.start()
            # simple local score: normalized by (len/(len+6)) to saturate
            score = float(ln) / (ln + 6.0)
            results.append({'start': m.start(), 'end': m.end(), 'base': 'A', 'len': ln, 'score': round(score, 6), 'seq': seq[m.start():m.end()]})
        # T runs
        for m in re.finditer(r'T{' + str(min_len) + r',}', seq):
            ln = m.end() - m.start()
            score = float(ln) / (ln + 6.0)
            results.append({'start': m.start(), 'end': m.end(), 'base': 'T', 'len': ln, 'score': round(score, 6), 'seq': seq[m.start():m.end()]})
        # sort by start
        results.sort(key=lambda x: x['start'])
        return results

    # -------------------------
    # Scoring helpers (interpretability)
    # -------------------------
    def phasing_score(self, apr: Dict[str, Any]) -> float:
        """Return APR phasing score (already stored in apr['score'])."""
        return float(apr.get('score', 0.0))

    def local_curvature_score(self, tract: Dict[str, Any]) -> float:
        """Return local curvature score for a long tract (already stored)."""
        return float(tract.get('score', 0.0))

    # -------------------------
    # Annotate (summary)
    # -------------------------
    def annotate_sequence(self, sequence: str) -> Dict[str, Any]:
        """
        Returns comprehensive annotation:
         - a_tract_windows: raw outputs from find_a_tracts
         - aprs: list of APRs with phasing scores
         - long_tracts: list of local A/T long tracts
         - summary counts and combined score
        """
        seq = sequence.upper()
        a_windows = self.find_a_tracts(seq, minAT=self.MIN_AT_TRACT)
        # filtered called a-tract centers
        a_centers = [w for w in a_windows if w['call'] and w['a_center'] is not None]
        aprs = self.find_aprs(seq, min_tract=self.MIN_AT_TRACT, min_apr_tracts=self.MIN_APR_TRACTS)
        long_tracts = self.find_long_tracts(seq, min_len=self.LOCAL_LONG_TRACT)

        # annotate aprs with constituent windows (optional)
        for apr in aprs:
            apr['constituent_windows'] = []
            for center in apr['center_positions']:
                # find closest a_window with that center
                best = min(a_windows, key=lambda w: abs((w['a_center'] or 1e9) - center))
                apr['constituent_windows'].append(best)

        summary = {
            'n_a_windows': len(a_windows),
            'n_a_centers': len(a_centers),
            'n_aprs': len(aprs),
            'n_long_tracts': len(long_tracts),
            'apr_score_sum': sum(self.phasing_score(a) for a in aprs),
            'long_tract_score_sum': sum(self.local_curvature_score(l) for l in long_tracts),
            'combined_score': sum(self.phasing_score(a) for a in aprs) + sum(self.local_curvature_score(l) for l in long_tracts)
        }

        return {
            'a_tract_windows': a_windows,
            'a_centers': a_centers,
            'aprs': aprs,
            'long_tracts': long_tracts,
            'summary': summary
        }
