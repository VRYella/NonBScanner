"""
Triplex DNA Motif Detector (Mirror repeat strict, content threshold)
===================================================================
Detects three-stranded DNA structures using optimized k-mer seed-and-extend approach.

Key features from Frank-Kamenetskii 1995, Sakamoto 1999, Bacolla 2006:
- Intramolecular triplex: mirror repeats ≥10 bp per arm, loop ≤100 bp, 
  homopurine or homopyrimidine content >90% in arms
- Sticky DNA: pure (GAA)n or (TTC)n (≥12 bp)

PERFORMANCE OPTIMIZATIONS:
- Uses optimized seed-and-extend k-mer index approach from repeat_scanner
- O(n) complexity with k-mer seeding for mirror repeats
- Efficient purine/pyrimidine content filtering
"""

import re
from typing import List, Dict, Any, Tuple

try:
    from .base_detector import BaseMotifDetector
except ImportError:
    class BaseMotifDetector: pass

# Import optimized repeat scanner
try:
    from utils.repeat_scanner import find_mirror_repeats
except ImportError:
    # Fallback if import fails
    find_mirror_repeats = None

class TriplexDetector(BaseMotifDetector):
    """Detector for mirror repeat and sticky DNA triplex motifs using optimized scanner"""

    def get_motif_class_name(self) -> str:
        return "Triplex"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        # Sticky DNA patterns still use regex (simple and efficient)
        # Mirror repeats now use optimized k-mer index from repeat_scanner
        return {
            'triplex_forming_sequences': [
                # GAA sticky DNA - optimized with non-capturing group
                (r'(?:GAA){4,}', 'TRX_5_4', 'GAA repeat', 'Sticky_DNA', 12, 'sticky_dna_score', 0.95, 'Disease-associated repeats', 'Sakamoto 1999'),
                # TTC sticky DNA - optimized with non-capturing group
                (r'(?:TTC){4,}', 'TRX_5_5', 'TTC repeat', 'Sticky_DNA', 12, 'sticky_dna_score', 0.95, 'Disease-associated repeats', 'Sakamoto 1999'),
            ]
        }
    
    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        seq = sequence.upper()
        results = []
        used = [False] * len(seq)
        patterns = self.get_patterns()['triplex_forming_sequences']

        # First, detect mirror repeats using optimized scanner if available
        if find_mirror_repeats is not None:
            from utils.repeat_scanner import find_mirror_repeats as optimized_find
            mirror_results = optimized_find(seq, min_arm=10, max_loop=100, purine_pyrimidine_threshold=0.9)
            
            # Only keep those that pass the triplex threshold (>90% purine or pyrimidine)
            for mr_rec in mirror_results:
                if mr_rec.get('Is_Triplex', False):
                    s = mr_rec['Start'] - 1  # Convert to 0-based
                    e = mr_rec['End']
                    
                    if any(used[s:e]):
                        continue
                    
                    for i in range(s, e):
                        used[i] = True
                    
                    # Determine if purine or pyrimidine
                    pur_frac = mr_rec['Purine_Fraction']
                    pyr_frac = mr_rec['Pyrimidine_Fraction']
                    if pur_frac >= 0.9:
                        subtype = 'Homopurine mirror repeat'
                        pid = 'TRX_MR_PU'
                    else:
                        subtype = 'Homopyrimidine mirror repeat'
                        pid = 'TRX_MR_PY'
                    
                    results.append({
                        'class_name': 'Triplex',
                        'pattern_id': pid,
                        'start': s,
                        'end': e,
                        'length': e - s,
                        'score': self._triplex_potential(seq[s:e]),
                        'matched_seq': seq[s:e],
                        'details': {
                            'type': subtype,
                            'reference': 'Frank-Kamenetskii 1995',
                            'description': 'H-DNA formation',
                            'arm_length': mr_rec['Arm_Length'],
                            'loop_length': mr_rec['Loop'],
                            'purine_fraction': pur_frac,
                            'pyrimidine_fraction': pyr_frac
                        }
                    })
        else:
            # Fallback to regex-based detection for mirror repeats
            # Homopurine mirror repeat
            pat_pu = r'((?:[GA]{1,}){10,})([ATGC]{1,100})((?:[GA]{1,}){10,})'
            for m in re.finditer(pat_pu, seq):
                s, e = m.span()
                if any(used[s:e]):
                    continue
                arm1 = m.group(1)
                arm2 = m.group(3)
                loop = m.group(2)
                if len(arm1) < 10 or len(arm2) < 10 or len(loop) > 100:
                    continue
                pur_ct = sum(1 for b in arm1+arm2 if b in 'AG') / max(1, len(arm1+arm2))
                if pur_ct <= 0.9:
                    continue
                for i in range(s, e):
                    used[i] = True
                results.append({
                    'class_name': 'Triplex',
                    'pattern_id': 'TRX_MR_PU',
                    'start': s,
                    'end': e,
                    'length': e-s,
                    'score': self._triplex_potential(seq[s:e]),
                    'matched_seq': seq[s:e],
                    'details': {
                        'type': 'Homopurine mirror repeat',
                        'reference': 'Frank-Kamenetskii 1995',
                        'description': 'H-DNA formation (homopurine)'
                    }
                })
            
            # Homopyrimidine mirror repeat
            pat_py = r'((?:[CT]{1,}){10,})([ATGC]{1,100})((?:[CT]{1,}){10,})'
            for m in re.finditer(pat_py, seq):
                s, e = m.span()
                if any(used[s:e]):
                    continue
                arm1 = m.group(1)
                arm2 = m.group(3)
                loop = m.group(2)
                if len(arm1) < 10 or len(arm2) < 10 or len(loop) > 100:
                    continue
                pyr_ct = sum(1 for b in arm1+arm2 if b in 'CT') / max(1, len(arm1+arm2))
                if pyr_ct <= 0.9:
                    continue
                for i in range(s, e):
                    used[i] = True
                results.append({
                    'class_name': 'Triplex',
                    'pattern_id': 'TRX_MR_PY',
                    'start': s,
                    'end': e,
                    'length': e-s,
                    'score': self._triplex_potential(seq[s:e]),
                    'matched_seq': seq[s:e],
                    'details': {
                        'type': 'Homopyrimidine mirror repeat',
                        'reference': 'Frank-Kamenetskii 1995',
                        'description': 'H-DNA formation (homopyrimidine)'
                    }
                })

        # Sticky DNA patterns (GAA/TTC) - use regex
        for patinfo in patterns:
            pat, pid, name, cname, minlen, scoretype, cutoff, desc, ref = patinfo
            for m in re.finditer(pat, seq):
                s, e = m.span()
                if any(used[s:e]):
                    continue
                for i in range(s, e):
                    used[i] = True
                results.append({
                    'class_name': cname,
                    'pattern_id': pid,
                    'start': s,
                    'end': e,
                    'length': e-s,
                    'score': self.calculate_score(seq[s:e], patinfo),
                    'matched_seq': seq[s:e],
                    'details': {
                        'type': name,
                        'reference': ref,
                        'description': desc
                    }
                })
        
        results.sort(key=lambda r: r['start'])
        return results

    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        scoring_method = pattern_info[5] if len(pattern_info) > 5 else 'triplex_potential'
        if scoring_method == 'triplex_potential':
            return self._triplex_potential(sequence)
        elif scoring_method == 'sticky_dna_score':
            return self._sticky_dna_score(sequence)
        else:
            return 0.0

    def _triplex_potential(self, sequence: str) -> float:
        """Score: tract length and purine/pyrimidine content (≥90%)"""
        if len(sequence) < 20:
            return 0.0
        pur = sum(1 for b in sequence if b in "AG") / len(sequence)
        pyr = sum(1 for b in sequence if b in "CT") / len(sequence)
        score = (pur if pur > 0.9 else 0) + (pyr if pyr > 0.9 else 0)
        # tract length bonus: scale for very long arms, up to 1.0
        return min(score * len(sequence) / 150, 1.0)
        
    def _sticky_dna_score(self, sequence: str) -> float:
        """Score for sticky DNA: repeat density and length"""
        if len(sequence) < 12:
            return 0.0
        gaa_count = sequence.count("GAA")
        ttc_count = sequence.count("TTC")
        rep_total = gaa_count + ttc_count
        density = (rep_total * 3) / len(sequence)
        extras = sum(len(m.group(0)) for m in re.finditer(r'(?:GAA){2,}', sequence)) + \
                 sum(len(m.group(0)) for m in re.finditer(r'(?:TTC){2,}', sequence))
        cons_bonus = extras / len(sequence) if len(sequence) else 0
        return min(0.7 * density + 0.3 * cons_bonus, 1.0)

    def passes_quality_threshold(self, sequence: str, score: float, pattern_info: Tuple) -> bool:
        """Lower threshold for triplex detection"""
        return score >= 0.2  # Lower threshold for better sensitivity

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Override base method to use sophisticated triplex detection"""
        sequence = sequence.upper().strip()
        motifs = []
        
        # Use the annotate_sequence method which has the sophisticated logic
        # and already includes GAA/TTC detection via patterns
        results = self.annotate_sequence(sequence)
        
        for i, result in enumerate(results):
            motifs.append({
                'ID': f"{sequence_name}_{result['pattern_id']}_{result['start']+1}",
                'Sequence_Name': sequence_name,
                'Class': self.get_motif_class_name(),
                'Subclass': result['details']['type'],
                'Start': result['start'] + 1,  # 1-based coordinates
                'End': result['end'],
                'Length': result['length'],
                'Sequence': result['matched_seq'],
                'Score': round(result['score'], 3),
                'Strand': '+',
                'Method': 'Triplex_detection',
                'Pattern_ID': result['pattern_id']
            })
        
        return motifs
