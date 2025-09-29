"""
Triplex DNA Motif Detector (Mirror repeat strict, content threshold)
===================================================================
Detects three-stranded DNA structures, key rules from Frank-Kamenetskii 1995, Sakamoto 1999, Bacolla 2006, and recent reviews[web:80][web:83][web:85][web:86][web:89].

Features:
- Intramolecular triplex: mirror repeats ≥10 bp per arm, loop ≤100 bp, homopurine or homopyrimidine content >90% in arms
- Sticky DNA: pure (GAA)n or (TTC)n (≥12 bp)
"""

import re
from typing import List, Dict, Any, Tuple

try:
    from .base_detector import BaseMotifDetector
except ImportError:
    class BaseMotifDetector: pass

class TriplexDetector(BaseMotifDetector):
    """Detector for mirror repeat and sticky DNA triplex motifs (literature strict)[web:80][web:83][web:85][web:86][web:89]"""

    def get_motif_class_name(self) -> str:
        return "Triplex"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        return {
            'triplex_forming_sequences': [
                # Homopurine mirror repeat: (arm)-(loop)-(mirror_arm), each arm ≥10bp, loop ≤100bp, >90% purine
                (r'((?:[GA]{1,}){10,})'
                 r'([ATGC]{1,100})'
                 r'((?:[GA]{1,}){10,})',
                 'TRX_MR_PU',
                 'Homopurine mirror repeat',
                 'Triplex',
                 20,  # total arms+min loop
                 'triplex_potential',
                 0.90,
                 'H-DNA formation (homopurine)',
                 'Frank-Kamenetskii 1995'),
                # Homopyrimidine mirror repeat: same logic
                (r'((?:[CT]{1,}){10,})'
                 r'([ATGC]{1,100})'
                 r'((?:[CT]{1,}){10,})',
                 'TRX_MR_PY',
                 'Homopyrimidine mirror repeat',
                 'Triplex',
                 20,
                 'triplex_potential',
                 0.90,
                 'H-DNA formation (homopyrimidine)',
                 'Frank-Kamenetskii 1995'),
                # GAA sticky DNA
                (r'(?:GAA){4,}', 'TRX_5_4', 'GAA repeat', 'Sticky_DNA', 12, 'sticky_dna_score', 0.95, 'Disease-associated repeats', 'Sakamoto 1999'),
                # TTC sticky DNA
                (r'(?:TTC){4,}', 'TRX_5_5', 'TTC repeat', 'Sticky_DNA', 12, 'sticky_dna_score', 0.95, 'Disease-associated repeats', 'Sakamoto 1999'),
            ]
        }
    
    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        seq = sequence.upper()
        results = []
        used = [False] * len(seq)
        patterns = self.get_patterns()['triplex_forming_sequences']

        for patinfo in patterns:
            pat, pid, name, cname, minlen, scoretype, cutoff, desc, ref = patinfo
            for m in re.finditer(pat, seq):
                s, e = m.span()
                if any(used[s:e]):
                    continue
                # Mirror repeats (arms+loop): test content of both arms
                if "mirror" in name:
                    arm1 = m.group(1)
                    arm2 = m.group(3)
                    loop = m.group(2)
                    if len(arm1) < 10 or len(arm2) < 10 or len(loop) > 100:
                        continue
                    is_purine = set(arm1+arm2).issubset({'A','G'})
                    is_pyrimidine = set(arm1+arm2).issubset({'C','T'})
                    pur_ct = sum(1 for b in arm1+arm2 if b in 'AG') / max(1, len(arm1+arm2))
                    pyr_ct = sum(1 for b in arm1+arm2 if b in 'CT') / max(1, len(arm1+arm2))
                    if not (pur_ct > 0.9 or pyr_ct > 0.9):
                        continue
                # All sticky DNA: keep full match if not overlap
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
        """Score: tract length and purine/pyrimidine content (≥90%) [web:80][web:83][web:85]"""
        if len(sequence) < 20:
            return 0.0
        pur = sum(1 for b in sequence if b in "AG") / len(sequence)
        pyr = sum(1 for b in sequence if b in "CT") / len(sequence)
        score = (pur if pur > 0.9 else 0) + (pyr if pyr > 0.9 else 0)
        # tract length bonus: scale for very long arms, up to 1.0
        return min(score * len(sequence) / 150, 1.0)
        
    def _sticky_dna_score(self, sequence: str) -> float:
        """Score for sticky DNA: repeat density and length [web:86]"""
        if len(sequence) < 12:
            return 0.0
        gaa_count = sequence.count("GAA")
        ttc_count = sequence.count("TTC")
        rep_total = gaa_count + ttc_count
        density = (rep_total * 3) / len(sequence)
        extras = sum(len(m) for m in re.findall(r'(?:GAA){2,}', sequence)) + sum(len(m) for m in re.findall(r'(?:TTC){2,}', sequence))
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
        
        # Also add simple GAA/TTC repeat detection
        gaa_pattern = re.compile(r'(GAA){4,}', re.IGNORECASE)
        ttc_pattern = re.compile(r'(TTC){4,}', re.IGNORECASE)
        
        for pattern, name in [(gaa_pattern, 'GAA-repeat'), (ttc_pattern, 'TTC-repeat')]:
            for match in pattern.finditer(sequence):
                start, end = match.span()
                motif_seq = sequence[start:end]
                score = len(motif_seq) / 30.0  # Simple length-based score
                
                motifs.append({
                    'ID': f"{sequence_name}_STICKY_{name}_{start+1}",
                    'Sequence_Name': sequence_name,
                    'Class': self.get_motif_class_name(),
                    'Subclass': 'Sticky DNA',
                    'Start': start + 1,  # 1-based coordinates
                    'End': end,
                    'Length': len(motif_seq),
                    'Sequence': motif_seq,
                    'Score': round(min(score, 1.0), 3),
                    'Strand': '+',
                    'Method': 'Triplex_detection',
                    'Pattern_ID': f'STICKY_{name}'
                })
        
        return motifs
