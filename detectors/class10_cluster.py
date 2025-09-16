"""
Non-B DNA Cluster Regions Detection (Class 10)
Dynamic: any combination of ≥3 motif classes, each occurring 3+ times in 100 nt; no static subclass list
"""

from .base_motif import standardize_motif_output


def find_hotspots(motif_hits, seq_len, window=100, min_count=3) -> list:
    """Find Non-B DNA cluster regions (hotspots)"""
    hotspots = []
    positions = [(hit['Start'], hit['End']) for hit in motif_hits]
    
    for i in range(0, seq_len - window + 1):
        region_start, region_end = i+1, i+window
        count = sum(s <= region_end and e >= region_start for s, e in positions)
        if count >= min_count:
            motifs_in_region = [m for m in motif_hits if m['Start'] <= region_end and m['End'] >= region_start]
            type_div = len({m.get('Subclass', m.get('Subtype', '')) for m in motifs_in_region})
            total_score = sum(float(m.get("Score", m.get("Actual_Score", 0.0))) for m in motifs_in_region)
            hotspots.append({
                "Class": "Non-B_DNA_Cluster",
                "Subclass": "Hotspot",
                "Start": region_start,
                "End": region_end,
                "Length": region_end - region_start + 1,
                "Sequence": "",
                "ScoreMethod": "Hotspot_raw",
                "Score": float(total_score),
                "MotifCount": count,
                "TypeDiversity": type_div,
            })
    
    return merge_hotspots(hotspots)


def merge_hotspots(hotspots) -> list:
    """Merge overlapping hotspots"""
    if not hotspots:
        return []
    
    merged = [hotspots[0]]
    for current in hotspots[1:]:
        last = merged[-1]
        if current['Start'] <= last['End']:
            last['End'] = max(last['End'], current['End'])
            last['Length'] = last['End'] - last['Start'] + 1
            last['MotifCount'] += current['MotifCount']
            last['TypeDiversity'] = max(last['TypeDiversity'], current['TypeDiversity'])
            last['Score'] = float(last['Score']) + float(current['Score'])
        else:
            merged.append(current)
    
    return merged


def find_cluster(motif_hits, seq_len: int, sequence_name: str = "") -> list:
    """Main function to find Non-B DNA cluster regions"""
    hotspots = find_hotspots(motif_hits, seq_len)
    
    # Standardize output format
    standardized_results = []
    for i, motif in enumerate(hotspots, 1):
        standardized_results.append(standardize_motif_output(motif, sequence_name, i))
    
    return standardized_results