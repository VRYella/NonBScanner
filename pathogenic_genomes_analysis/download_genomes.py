#!/usr/bin/env python3
"""
Download small pathogenic genomes from NCBI for NonBScanner analysis.
Selected genomes represent important viral and bacterial pathogens.
"""

import os
import urllib.request
import time

# Directory for genome data
DATA_DIR = "data"
os.makedirs(DATA_DIR, exist_ok=True)

# Selected pathogenic genomes with NCBI accessions
# Selected for small size (<50kb) and biological significance
GENOMES = {
    "SARS-CoV-2": {
        "accession": "NC_045512.2",
        "description": "Severe acute respiratory syndrome coronavirus 2 (COVID-19)",
        "size": "~30kb"
    },
    "Hepatitis_B_Virus": {
        "accession": "NC_003977.2",
        "description": "Hepatitis B virus genotype D",
        "size": "~3.2kb"
    },
    "Human_Papillomavirus_16": {
        "accession": "NC_001526.4",
        "description": "Human papillomavirus type 16 (cervical cancer)",
        "size": "~8kb"
    },
    "Ebola_Virus": {
        "accession": "NC_002549.1",
        "description": "Zaire ebolavirus (Ebola hemorrhagic fever)",
        "size": "~19kb"
    },
    "Influenza_A_H1N1": {
        "accession": "NC_026433.1",
        "description": "Influenza A virus segment 4 (HA gene)",
        "size": "~1.7kb"
    }
}

def download_genome(accession, name):
    """Download genome from NCBI in FASTA format."""
    url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id={accession}"
    output_file = os.path.join(DATA_DIR, f"{name}.fasta")
    
    print(f"Downloading {name} ({accession})...")
    try:
        urllib.request.urlretrieve(url, output_file)
        print(f"  ✓ Saved to {output_file}")
        time.sleep(1)  # Be nice to NCBI servers
        return True
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return False

def main():
    """Download all pathogenic genomes."""
    print("=" * 70)
    print("DOWNLOADING PATHOGENIC GENOMES FOR NONBSCANNER ANALYSIS")
    print("=" * 70)
    print()
    
    success_count = 0
    for name, info in GENOMES.items():
        print(f"Genome: {info['description']}")
        print(f"Size: {info['size']}")
        if download_genome(info['accession'], name):
            success_count += 1
        print()
    
    print("=" * 70)
    print(f"Download complete: {success_count}/{len(GENOMES)} genomes downloaded")
    print("=" * 70)

if __name__ == "__main__":
    main()
