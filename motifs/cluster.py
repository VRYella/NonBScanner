"""
Non-B DNA Cluster Regions Detection (Class 10)
Dynamic: any combination of ≥3 motif classes, each occurring 3+ times in 100 nt; no static subclass list
"""

from .base_motif import standardize_motif_output


def find_hotspots(motif_hits, seq_len, window=100, min_count=3) -> list:
    """Find Non-B DNA cluster regions (hotspots) - normalized scoring only"""
    hotspots = []
    positions = [(hit['Start'], hit['End']) for hit in motif_hits]
    
    for i in range(0, seq_len - window + 1):
        region_start, region_end = i+1, i+window
        count = sum(s <= region_end and e >= region_start for s, e in positions)
        if count >= min_count:
            motifs_in_region = [m for m in motif_hits if m['Start'] <= region_end and m['End'] >= region_start]
            type_div = len({m.get('Subclass', m.get('Subtype', '')) for m in motifs_in_region})
            
            # Calculate normalized score based on density and diversity
            density_factor = min(1.0, count / 20.0)  # Normalize by high density (20 motifs in 100bp)
            diversity_factor = min(1.0, type_div / 8.0)  # Normalize by max possible classes (8)
            
            # Get normalized scores from motifs in region
            region_norm_scores = []
            for m in motifs_in_region:
                norm_score = m.get('Normalized_Score', m.get('NormScore', 0.0))
                if norm_score:
                    region_norm_scores.append(float(norm_score))
            
            # Cluster normalized score: average normalized score * density * diversity
            if region_norm_scores:
                avg_norm_score = sum(region_norm_scores) / len(region_norm_scores)
                normalized_score = avg_norm_score * density_factor * diversity_factor
            else:
                normalized_score = density_factor * diversity_factor
            
            # Ensure normalized score is in [0,1] range
            normalized_score = max(0.0, min(1.0, normalized_score))
            
            hotspots.append({
                "Class": "Non-B_DNA_Cluster",
                "Subclass": "Hotspot",
                "Start": region_start,
                "End": region_end,
                "Length": region_end - region_start + 1,
                "Sequence": "",
                "ScoreMethod": "Hotspot_normalized",
                "Normalized_Score": normalized_score,
                "MotifCount": count,
                "TypeDiversity": type_div,
            })
    
    return merge_hotspots(hotspots)


def merge_hotspots(hotspots) -> list:
    """Merge overlapping hotspots - update normalized scores"""
    if not hotspots:
        return []
    
    merged = [hotspots[0]]
    for current in hotspots[1:]:
        last = merged[-1]
        if current['Start'] <= last['End']:
            # Merge regions
            last['End'] = max(last['End'], current['End'])
            last['Length'] = last['End'] - last['Start'] + 1
            last['MotifCount'] += current['MotifCount']
            last['TypeDiversity'] = max(last['TypeDiversity'], current['TypeDiversity'])
            
            # Merge normalized scores (take maximum as it represents peak cluster density)
            last_norm = last.get('Normalized_Score', 0.0)
            current_norm = current.get('Normalized_Score', 0.0)
            last['Normalized_Score'] = max(float(last_norm), float(current_norm))
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