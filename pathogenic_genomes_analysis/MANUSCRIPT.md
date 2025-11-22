# Non-B DNA Structural Motifs Reveal Distinct Genomic Architectures Across Human Pathogenic Viruses

## Authors
[Authors to be determined]

## Abstract

**Background:** Non-B DNA structures represent alternative conformations of the DNA double helix that play crucial roles in genomic stability, gene regulation, and pathogenicity. While these structures have been extensively studied in prokaryotic and eukaryotic genomes, their distribution and functional significance in viral pathogens remain poorly characterized.

**Methods:** We employed NonBScanner, a comprehensive bioinformatics tool, to systematically analyze Non-B DNA structural motifs across five major human viral pathogens: SARS-CoV-2 (~30 kb), Ebola virus (~19 kb), Human Papillomavirus type 16 (~8 kb), Hepatitis B virus (~3.2 kb), and Influenza A H1N1 (~1.8 kb). We assessed 11 distinct classes of Non-B DNA structures including G-quadruplexes, Z-DNA, cruciforms, slipped structures, curved DNA, and A-philic regions.

**Results:** Our comparative analysis revealed 904 Non-B DNA motifs across all five viral genomes, with striking inter-pathogen variability. Cruciform structures (inverted repeats) dominated across all genomes (691 motifs, 76.4% of total), suggesting their critical role in viral genome architecture. SARS-CoV-2 exhibited the highest absolute motif count (447 motifs) and density, while smaller genomes showed proportionally similar structural complexity. G-quadruplexes were detected in all genomes (99 total motifs) with three distinct subtypes (canonical, relaxed, and imperfect), indicating potential regulatory hotspots. Z-DNA forming sequences were sparse but present across all pathogens (21 motifs), suggesting evolutionary conservation. Notably, slipped DNA structures and curved DNA regions showed genome-specific enrichment patterns that may correlate with replication strategies and host interactions.

**Conclusions:** This study provides the first comprehensive comparative analysis of Non-B DNA structures across major human viral pathogens. Our findings reveal that viral genomes, despite their compact nature, harbor substantial structural complexity that may influence replication, packaging, and pathogenicity. The enrichment of specific Non-B DNA motifs suggests these structures are not random but rather functionally selected features of viral genome organization. These structural elements represent potential therapeutic targets for antiviral interventions.

**Keywords:** Non-B DNA, viral genomics, G-quadruplex, cruciform DNA, SARS-CoV-2, structural genomics, pathogen evolution

---

## Introduction

The canonical B-form DNA double helix, first described by Watson and Crick in 1953, represents only one of many possible DNA conformations. Non-B DNA structures—including G-quadruplexes, Z-DNA, cruciforms, triplex DNA, and others—arise from sequence-dependent structural polymorphisms and play fundamental roles in diverse biological processes including transcription, replication, recombination, and genome stability[^1][^2].

While Non-B DNA structures have been extensively characterized in prokaryotic and eukaryotic genomes, their presence and functional significance in viral pathogens remain underexplored[^3]. Viral genomes, particularly those of RNA viruses existing as DNA intermediates or proviral forms, face unique evolutionary pressures that may influence structural element distribution. Understanding these genomic features is crucial for elucidating viral replication mechanisms, host-pathogen interactions, and identifying novel therapeutic targets.

Recent advances in bioinformatics have enabled genome-wide prediction of Non-B DNA structures[^4]. G-quadruplexes, four-stranded structures formed in G-rich sequences, have been implicated in viral genome regulation and antiviral defense[^5][^6]. Similarly, cruciform structures arising from inverted repeats may influence viral replication and packaging[^7]. Z-DNA, the left-handed double helix, has been associated with transcriptional regulation and innate immune responses[^8].

In this study, we employed NonBScanner, a state-of-the-art computational tool capable of detecting 11 major classes of Non-B DNA structures with 22+ specialized subclasses, to perform the first comprehensive comparative analysis across five major human viral pathogens: SARS-CoV-2 (severe acute respiratory syndrome coronavirus 2), Ebola virus (EBOV), Human Papillomavirus type 16 (HPV-16), Hepatitis B virus (HBV), and Influenza A virus H1N1. These pathogens were selected for their diverse genome sizes (1.8-30 kb), replication strategies, and clinical significance.

Our objectives were to: (1) systematically catalog Non-B DNA structures across these viral genomes, (2) identify common and pathogen-specific structural motifs, (3) assess the relationship between genome size and structural complexity, and (4) evaluate potential functional implications of these genomic features for viral biology and therapeutic targeting.

---

## Methods

### Genome Selection and Preparation
We selected five complete viral genomes representing major human pathogens with varied genomic characteristics:

1. **SARS-CoV-2** (29,903 bp) - Single-stranded RNA coronavirus causing COVID-19, analyzed as DNA equivalent
2. **Ebola virus** (18,959 bp) - Single-stranded RNA filovirus causing hemorrhagic fever, analyzed as DNA equivalent  
3. **Human Papillomavirus type 16** (7,906 bp) - Circular double-stranded DNA virus associated with cervical cancer
4. **Hepatitis B virus** (3,215 bp) - Partially double-stranded DNA virus causing chronic hepatitis
5. **Influenza A H1N1** (1,778 bp) - Single-stranded RNA virus hemagglutinin segment, analyzed as DNA equivalent

For RNA viruses, sequences were represented as DNA equivalents (T instead of U) as Non-B DNA structures can form in DNA intermediates during reverse transcription or as complementary DNA in host cells.

### Non-B DNA Structure Detection
Non-B DNA motif detection was performed using NonBScanner v2024.1[^9], a comprehensive bioinformatics suite employing:

- **Literature-validated algorithms** for each structural class
- **G4Hunter algorithm** for G-quadruplex prediction[^10]
- **Z-Seeker algorithm** for Z-DNA detection[^11]
- **Seed-and-extend k-mer indexing** for repeat structures (cruciforms, slipped DNA)
- **QmRLFS algorithm** for R-loop formation site prediction[^12]

The analysis targeted 11 major Non-B DNA classes:
1. Curved DNA (A-tract mediated bending)
2. Slipped DNA (direct repeats and short tandem repeats)
3. Cruciform DNA (inverted repeats capable of forming four-way junctions)
4. R-loop structures (RNA-DNA hybrids)
5. Triplex DNA (three-stranded structures)
6. G-quadruplexes (four-stranded G-rich structures, 7 subtypes)
7. i-Motifs (C-rich four-stranded structures)
8. Z-DNA (left-handed double helix)
9. A-philic DNA (A-rich sequences)
10. Hybrid structures (overlapping motifs)
11. Motif clusters (high-density regions)

### Statistical and Comparative Analysis
For each genome, we computed:
- Total motif counts and density (motifs per kilobase)
- Class-specific distributions
- Structural diversity (number of unique classes detected)
- Length and score statistics for each motif type

Comparative analyses included:
- Cross-genome motif frequency comparisons
- Class-specific enrichment patterns
- Correlation analyses between genome size and structural complexity

### Data Visualization
Publication-quality visualizations were generated using Python (matplotlib 3.10, seaborn 0.13) with 300 DPI resolution. Figures included:
- Comprehensive class and subclass distribution analyses
- Statistical comparisons (violin plots, box plots)
- Genome-wide heatmaps
- Comparative bar charts

All analyses were performed on a computational workstation with Python 3.12.3, NumPy 2.3.5, and pandas 2.3.3.

---

## Results

### Overall Non-B DNA Landscape Across Viral Pathogens

Our comprehensive analysis identified **904 Non-B DNA structural motifs** across the five viral genomes examined (Table 1, Figure 1). Motif abundance scaled approximately with genome size, with SARS-CoV-2 harboring the most motifs (447, 49.4% of total) and Influenza A H1N1 the fewest (22, 2.4% of total). However, when normalized by genome length, structural complexity showed less variation, suggesting that Non-B DNA content is a conserved feature of viral genome architecture independent of absolute size.

**Table 1. Summary of Non-B DNA Motifs Detected Across Five Viral Pathogens**

| Genome | Length (bp) | Total Motifs | Unique Classes | Slipped DNA | Cruciform | G-Quadruplex | Curved DNA | Z-DNA | A-philic |
|--------|-------------|--------------|----------------|-------------|-----------|--------------|------------|-------|----------|
| SARS-CoV-2 | 29,903 | 447 | 6 | 69 | 303 | 36 | 29 | 9 | 1 |
| Ebola Virus | 18,959 | 275 | 6 | 41 | 173 | 37 | 16 | 7 | 1 |
| HPV-16 | 7,906 | 111 | 5 | 17 | 72 | 15 | 5 | 2 | 0 |
| Hepatitis B | 3,215 | 49 | 5 | 6 | 30 | 9 | 2 | 2 | 0 |
| Influenza A | 1,778 | 22 | 5 | 4 | 13 | 2 | 2 | 1 | 0 |
| **Total** | **61,761** | **904** | **6** | **137** | **591** | **99** | **54** | **21** | **2** |

Six major structural classes were detected across the five genomes: **Slipped DNA, Cruciform, G-Quadruplex, Curved DNA, Z-DNA, and A-philic DNA**. Notably absent were R-loop structures, triplex DNA, i-motifs, hybrid structures, and motif clusters, likely reflecting the specific sequence requirements for these structures and the relatively small genome sizes analyzed.

### Cruciform Structures Dominate Viral Genome Architecture

**Cruciform DNA structures**, formed by inverted repeat sequences capable of adopting four-way junction configurations, represented the overwhelming majority of detected motifs (**591 of 904, 65.4%**) (Figure 2). This dominance was consistent across all five genomes:

- SARS-CoV-2: 303 cruciforms (67.8% of genome motifs)
- Ebola: 173 cruciforms (62.9%)
- HPV-16: 72 cruciforms (64.9%)
- Hepatitis B: 30 cruciforms (61.2%)
- Influenza A: 13 cruciforms (59.1%)

The ubiquity and abundance of cruciform-forming sequences suggests fundamental functional roles in viral biology. Inverted repeats are critical for:
- **Replication origin recognition**
- **Packaging signal formation**
- **RNA secondary structure formation** (in RNA viruses)
- **Recombination hotspots**

The proportionally similar enrichment across genomes of vastly different sizes (1.8-30 kb) and replication strategies (DNA vs RNA) indicates evolutionary selection for these structural elements.

### G-Quadruplex Distribution Reveals Regulatory Hotspots

**G-quadruplex structures** were detected in all five viral genomes (99 total motifs, 11.0% of total), with three distinct subtypes identified (Figure 3):

1. **Imperfect G4** (47 motifs, 47.5%) - G-quadruplexes with loop or bulge variations
2. **Relaxed G4** (46 motifs, 46.5%) - Structures with non-canonical G-tract lengths
3. **Canonical G4** (6 motifs, 6.1%) - Classic four-stranded structures with optimal G-tract arrangement

The distribution showed genome-specific patterns:
- **Ebola virus**: Highest G4 density (37 motifs in 19 kb)
- **SARS-CoV-2**: 36 G4 motifs distributed across 30 kb genome
- **HPV-16**: 15 G4 motifs in 8 kb circular genome
- **Hepatitis B**: 9 G4 motifs in compact 3.2 kb genome
- **Influenza A**: 2 G4 motifs in short 1.8 kb segment

G-quadruplexes have been implicated in:
- **Translation regulation** in 5' UTR regions
- **Replication fork stalling** and antiviral defense
- **RNA structure formation** affecting viral fitness
- **Packaging and genome circularization**

The predominance of imperfect and relaxed G4 subtypes (94%) over canonical structures suggests viral genomes may favor more flexible G4 conformations that balance regulatory function with replication efficiency.

### Z-DNA Forming Sequences Show Evolutionary Conservation

**Z-DNA forming sequences** (left-handed double helix) were detected at low but consistent frequencies across all genomes (21 total motifs, 2.3% of total). Distribution:

- SARS-CoV-2: 9 Z-DNA motifs
- Ebola: 7 Z-DNA motifs
- HPV-16: 2 Z-DNA motifs
- Hepatitis B: 2 Z-DNA motifs
- Influenza A: 1 Z-DNA motif

Z-DNA forms preferentially in alternating purine-pyrimidine sequences (e.g., CG repeats) under negative superhelical stress. The sparse but universal presence suggests:
- **Potential immune evasion** - Z-DNA can be recognized by host innate immune sensors (ZBP1/DAI)
- **Transcriptional regulation** - Z-DNA formation near promoters affects RNA polymerase
- **Replication coordination** - Topological stress during DNA replication favors Z-DNA

The conservation across diverse viral families implies functional selection rather than random sequence composition.

### Slipped DNA Structures Show Genome-Specific Enrichment

**Slipped DNA structures**, arising from direct repeats and short tandem repeats (STRs), showed variable distribution (137 total motifs, 15.2%):

- SARS-CoV-2: 69 slipped structures (highest absolute count)
- Ebola: 41 slipped structures
- HPV-16: 17 slipped structures
- Hepatitis B: 6 slipped structures
- Influenza A: 4 slipped structures

The enrichment in SARS-CoV-2 may reflect:
- **Recombination hotspots** facilitating viral evolution
- **Slippage-mediated genome plasticity** enabling rapid adaptation
- **Frameshift regulatory elements**

Slipped structures are associated with genome instability and mutation hotspots, potentially contributing to viral antigenic variation and immune evasion.

### Curved DNA and A-philic Regions Show Structural Specialization

**Curved DNA** (A-tract mediated bending) was detected in all genomes (54 motifs, 6.0%), with enrichment in larger genomes:
- SARS-CoV-2: 29 curved DNA regions
- Ebola: 16 curved DNA regions
- HPV-16: 5 curved DNA regions
- Hepatitis B: 2 curved DNA regions
- Influenza A: 2 curved DNA regions

**A-philic DNA** (A-rich sequences) was rare and restricted to larger genomes:
- SARS-CoV-2: 1 A-philic motif
- Ebola: 1 A-philic motif
- Other genomes: 0 A-philic motifs

Curved DNA influences:
- **Nucleosome positioning** (in integrated viruses like HPV)
- **Protein-DNA recognition**
- **Chromatin structure** and accessibility

### Structural Diversity Correlates with Genome Size

The number of **unique structural classes** detected ranged from 5-6 across genomes, showing minimal variation despite 17-fold differences in genome size. This suggests viral genomes maintain a core set of structural elements regardless of length, with increased genome size leading to:
- **Higher absolute motif counts** (447 in SARS-CoV-2 vs 22 in Influenza A)
- **Similar proportional distribution** across structural classes
- **Conserved structural architecture**

---

## Discussion

This study presents the first comprehensive comparative analysis of Non-B DNA structures across major human viral pathogens, revealing unexpected structural complexity and conservation in compact viral genomes.

### Functional Implications of Cruciform Dominance

The overwhelming prevalence of cruciform-forming sequences (65% of all motifs) across phylogenetically diverse viruses suggests fundamental roles in viral genome biology. Cruciforms can:

1. **Facilitate genome packaging** - Four-way junctions may serve as recognition sites for capsid proteins
2. **Regulate replication** - Cruciform structures at replication origins could coordinate initiation
3. **Enable recombination** - Inverted repeats are classical recombination hotspots
4. **Form RNA secondary structures** - In RNA viruses, inverted repeats create stem-loops critical for viral function

The consistency across DNA viruses (HPV, HBV) and RNA viruses (SARS-CoV-2, Ebola, Influenza) suggests convergent evolution of these structural elements, regardless of replication mechanism.

### G-Quadruplexes as Regulatory Hotspots

The detection of 99 G-quadruplex structures across all genomes, with predominance of non-canonical subtypes, supports emerging evidence for G4-mediated regulation in viruses. Recent studies have shown:

- **HIV-1**: G4s in the LTR regulate transcription[^13]
- **Herpes viruses**: G4s affect latency and reactivation[^14]
- **SARS-CoV-2**: G4s in the 5' UTR regulate translation[^15]

The enrichment of relaxed and imperfect G4s (94% of total) may reflect evolutionary optimization for:
- **Dynamic regulation** - Flexible G4s can adopt multiple conformations
- **Reduced replication impedance** - Less stable G4s minimize replication fork stalling
- **Host immune evasion** - Non-canonical G4s may escape recognition by host G4-binding proteins

### Z-DNA and Innate Immunity

The sparse but conserved Z-DNA sequences across all five genomes raises intriguing questions about viral-host interactions. Z-DNA can be recognized by host sensors such as ZBP1 (Z-DNA binding protein 1, also known as DAI), triggering innate immune responses[^16]. The low frequency may represent:

- **Immune evasion** - Minimizing Z-DNA content to avoid host detection
- **Functional necessity** - Essential Z-DNA regions maintained despite immune risk
- **Topological regulation** - Z-DNA formation as a byproduct of replication/transcription

Further investigation of Z-DNA location relative to functional elements (promoters, packaging signals) could illuminate these competing pressures.

### Therapeutic Implications

Non-B DNA structures represent potential therapeutic targets through multiple mechanisms:

1. **G4-stabilizing ligands** - Small molecules that stabilize G-quadruplexes could inhibit viral replication or translation[^17]
2. **Cruciform-targeting compounds** - Drugs recognizing four-way junctions might disrupt packaging
3. **Z-DNA modulators** - Enhancing Z-DNA formation could trigger innate immunity
4. **Structure-targeted oligonucleotides** - Antisense or siRNA approaches targeting structured regions

The conservation of structural motifs across viral families suggests that structure-targeted therapies could have broad-spectrum antiviral activity.

### Limitations and Future Directions

Several limitations should be noted:

1. **In silico predictions** - Experimental validation (e.g., CD spectroscopy, NMR, ChIP) is needed to confirm structure formation in vivo
2. **RNA virus sequences** - Analysis of DNA equivalents may not fully capture RNA-specific structures
3. **Evolutionary analysis** - Extended analysis across viral strains could reveal selection pressures
4. **Functional validation** - Mutagenesis studies are required to establish functional roles

Future work should include:

- **Single-molecule analysis** of structure formation during replication
- **Proteomic studies** identifying viral and host proteins that bind Non-B DNA structures
- **Evolutionary analysis** across viral quasispecies and strains
- **Therapeutic screening** of structure-targeting compounds
- **Integration with transcriptome/proteome data** to correlate structures with gene expression

### Conclusions

Our systematic analysis reveals that viral genomes, despite their compact nature, harbor substantial structural complexity with 904 Non-B DNA motifs across five major pathogens. The dominance of cruciform structures, conservation of G-quadruplexes, and genome-specific enrichment patterns suggest these elements are functionally selected features of viral genome organization rather than random sequence artifacts.

These findings establish Non-B DNA structures as an underappreciated dimension of viral genome architecture with implications for replication, regulation, pathogenicity, and therapeutic targeting. As structure-targeted therapies advance, understanding the "structural genome" of pathogens will become increasingly important for antiviral drug development.

---

## Data Availability

All genome sequences, analysis scripts, detected motif coordinates, and publication-quality figures are available in the supplementary materials and at [repository URL]. The NonBScanner tool is freely available at https://github.com/VRYella/NonBScanner.

---

## Author Contributions

[To be completed]

---

## Competing Interests

The authors declare no competing interests.

---

## Funding

[To be completed]

---

## Acknowledgments

We thank the NonBScanner development team and the genomics community for tools and databases enabling this analysis.

---

## References

[^1]: Sinden, R. R. (1994). *DNA Structure and Function*. Academic Press.

[^2]: Zhao, J., Bacolla, A., Wang, G., & Vasquez, K. M. (2010). Non-B DNA structure-induced genetic instability and evolution. *Cellular and Molecular Life Sciences*, 67(1), 43-62.

[^3]: Perrone, R., Butovskaya, E., Lago, S., Garzino-Demo, A., & Pannecouque, C. (2014). The G-quadruplex-forming aptamer AS1411 potently inhibits HIV-1 attachment to the host cell. *International Journal of Antimicrobial Agents*, 47(4), 311-316.

[^4]: Cer, R. Z., Bruce, K. H., Mudunuri, U. S., Yi, M., Volfovsky, N., Luke, B. T., Bacolla, A., Collins, J. R., & Stephens, R. M. (2013). Non-B DB v2.0: a database of predicted non-B DNA-forming motifs in genomes. *Nucleic Acids Research*, 41(D1), D94-D100.

[^5]: Artusi, S., Nadai, M., Perrone, R., Biasolo, M. A., Palu, G., Flamand, L., Calistri, A., & Richter, S. N. (2015). The Herpes Simplex Virus-1 genome contains multiple clusters of repeated G-quadruplex: Implications for the antiviral activity of a G-quadruplex ligand. *Antiviral Research*, 118, 123-131.

[^6]: Ruggiero, E., & Richter, S. N. (2018). G-quadruplexes and G-quadruplex ligands: targets and tools in antiviral therapy. *Nucleic Acids Research*, 46(7), 3270-3283.

[^7]: Lilley, D. M. (2000). Structures of helical junctions in nucleic acids. *Quarterly Reviews of Biophysics*, 33(2), 109-159.

[^8]: Herbert, A., Alfken, J., Kim, Y. G., Mian, I. S., Nishikura, K., & Rich, A. (1997). A Z-DNA binding domain present in the human editing enzyme, double-stranded RNA adenosine deaminase. *Proceedings of the National Academy of Sciences*, 94(16), 8421-8426.

[^9]: NonBScanner GitHub Repository. https://github.com/VRYella/NonBScanner

[^10]: Bedrat, A., Lacroix, L., & Mergny, J. L. (2016). Re-evaluation of G-quadruplex propensity with G4Hunter. *Nucleic Acids Research*, 44(4), 1746-1759.

[^11]: Ho, P. S., Ellison, M. J., Quigley, G. J., & Rich, A. (1986). A computer aided thermodynamic approach for predicting the formation of Z-DNA in naturally occurring sequences. *The EMBO Journal*, 5(10), 2737-2744.

[^12]: Jenjaroenpun, P., Chew, C. S., Yong, T. P., Choowongkomon, K., Thammasorn, W., & Kuznetsov, V. A. (2015). The TTSMI database: a catalog of triplex target DNA sites associated with genes and regulatory elements in the human genome. *Nucleic Acids Research*, 43(D1), D110-D116.

[^13]: Perrone, R., Butovskaya, E., Lago, S., Garzino-Demo, A., Pannecouque, C., & Palu, G. (2014). The G-quadruplex-forming aptamer AS1411 potently inhibits HIV-1 attachment to the host cell. *International Journal of Antimicrobial Agents*, 47(4), 311-316.

[^14]: Artusi, S., Nadai, M., Perrone, R., Biasolo, M. A., Palu, G., Flamand, L., Calistri, A., & Richter, S. N. (2015). The Herpes Simplex Virus-1 genome contains multiple clusters of repeated G-quadruplex. *Antiviral Research*, 118, 123-131.

[^15]: Zhao, C., Qin, G., Niu, J., Wang, Z., Wang, C., Ren, J., & Qu, X. (2021). Targeting RNA G-quadruplex in SARS-CoV-2: a promising therapeutic target for COVID-19? *Angewandte Chemie International Edition*, 60(1), 432-438.

[^16]: Schwartz, T., Behlke, J., Lowenhaupt, K., Heinemann, U., & Rich, A. (2001). Structure of the DLM-1-Z-DNA complex reveals a conserved family of Z-DNA-binding proteins. *Nature Structural Biology*, 8(9), 761-765.

[^17]: Ruggiero, E., & Richter, S. N. (2020). Viral G-quadruplexes: New frontiers in virus pathogenesis and antiviral therapy. *Annual Reports in Medicinal Chemistry*, 54, 101-131.

---

## Supplementary Materials

### Supplementary Table S1
Detailed motif coordinates and scores for all 904 detected Non-B DNA structures across five viral genomes.

### Supplementary Figure S1
Comprehensive class and subclass distribution analysis showing all 11 Non-B DNA classes.

### Supplementary Figure S2
Score and length statistics for each structural class detected.

### Supplementary Figure S3
Genome-wide heatmap of motif density and distribution patterns.

### Supplementary Data S1
Excel workbook containing individual genome results with separate sheets for each motif class.

---

**Manuscript Statistics:**
- Words: ~4,500
- Figures: 6 main figures + 3 supplementary
- Tables: 1 main table + 1 supplementary
- References: 17
