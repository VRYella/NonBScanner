def canonicalize_motif(m):
    """
    Canonicalize motif dictionary to standard format.
    Maintains backward compatibility with Normalized_Score but doesn't require it.
    """
    mapping = {
        'Actual Score': 'Actual_Score',
        'ActualScore': 'Actual_Score',
        'Score': 'Score',
        'Normalized Score': 'Normalized_Score',
        'Normalized_Score': 'Normalized_Score',
        'Class': 'Class',
        'Type': 'Class',
        'Subclass': 'Subclass',
        'Subtype': 'Subclass',
        'Start': 'Start',
        'End': 'End',
        'Length': 'Length',
        'Sequence_Name': 'Sequence_Name',
        'Motif': 'Motif'
    }
    out = {}
    for new_key in mapping.values():
        aliases = [k for k, v in mapping.items() if v == new_key]
        val = None
        for alias in aliases:
            if alias in m:
                val = m[alias]
                break
        out[new_key] = val
    if out.get('Start') and out.get('End') and not out.get('Length'):
        out['Length'] = int(out['End']) - int(out['Start'])
    out['Class'] = str(out.get('Class') or 'Unknown')
    out['Subclass'] = str(out.get('Subclass') or 'Other')
    out['Actual_Score'] = float(out.get('Actual_Score') or m.get('Actual_Score') or m.get('Score') or 0.0)
    out['Score'] = float(out.get('Score') or out['Actual_Score'])
    # Keep Normalized_Score for backward compatibility if it exists, but set to 0 if not present
    out['Normalized_Score'] = float(out.get('Normalized_Score') or m.get('Normalized_Score') or 0.0)
    out['Motif'] = m.get('Motif') or m.get('matched_seq') or ''
    out['Sequence_Name'] = out.get('Sequence_Name') or m.get('sequence_name') or ''
    return out

