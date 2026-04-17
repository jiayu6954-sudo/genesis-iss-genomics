# genesis-iss-genomics

**Coordinated loss of CRISPR-Cas systems and insertion sequences in Bacillales from the International Space Station**

> Supplementary data and analysis code for the manuscript submitted to *npj Microgravity*.

---

## Overview

This repository contains all analysis scripts, derived data tables, and supplementary figures for a species-stratified comparative genomics study of 85 Bacillales genomes (30 ISS isolates from PRJNA637984 + 55 ground controls).

**Core findings:**
- *Paenibacillus polymyxa* ISS isolates show complete CRISPR-Cas absence (0/9 vs 8/15 ground; Fisher *p* = 0.0095, BH *q* = 0.025)
- IS element density reduced in ISS *B. thuringiensis* (0.115-fold) and *P. polymyxa* (0.312-fold); both Cliff δ = −1.0 [−1.0,−1.0]
- CRISPR × IS co-loss within *P. polymyxa* (*p* = 0.003, δ = −0.75)
- GLDS-224 independent metagenomics validation: ISS CRISPR KO = 0.503× ground

---

## Repository Structure

```
science_engine/
├── analysis/                   # Core analysis scripts
│   ├── paper_readiness_gap_analysis.py   # Species-stratified Fisher + MWU
│   ├── rigorous_validation.py            # 7-item methodological validation
│   ├── glds224_defense_analysis.py       # GLDS-224 metagenomics
│   ├── osd582_is_element_analysis.py     # IS element density (OSD-582)
│   └── generate_figures.py               # Figure generation (matplotlib)
│
├── preprocess/output/          # Derived data tables
│   ├── gap_analysis_species_stratified.csv   # IS/CRISPR per genome (85 strains)
│   └── ground_gene_features.csv              # Ground genome feature matrix
│
├── figures/                    # Publication figures (PDF + PNG, 300 DPI)
│   ├── Figure1.pdf / Figure1.png   # CRISPR prevalence + IS density
│   ├── Figure2.pdf / Figure2.png   # Co-loss scatter + GLDS-224 validation
│   └── FigureS1.pdf / FigureS1.png # Assembly quality confound validation
│
├── reports/
│   ├── manuscript_draft_v1.md          # Full manuscript draft
│   ├── rigorous_validation_report.md   # V-001 through V-007 results
│   └── paper_readiness_gap_report.md   # Gap analysis output
│
└── knowledge_base/
    └── verified_hypotheses.json        # EXP-001 through EXP-012
```

---

## Data Sources

| Dataset | Access |
|---------|--------|
| ISS isolate genomes (30 Bacillales) | NCBI BioProject [PRJNA637984](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA637984) |
| Ground comparison genomes (55) | NCBI RefSeq — accessions in Supplementary Table S2 |
| GLDS-224 metagenomics | [NASA GeneLab](https://genelab.nasa.gov/) |

---

## Requirements

```bash
python >= 3.8
numpy
matplotlib
scipy   # optional — all statistics implemented in pure Python
```

---

## Reproduce analyses

```bash
# Species-stratified CRISPR + IS analysis
python analysis/paper_readiness_gap_analysis.py

# 7-item methodological validation
python analysis/rigorous_validation.py

# Generate all figures
python analysis/generate_figures.py
```

---

## Citation

ZJY. Coordinated loss of CRISPR-Cas systems and insertion sequences in Bacillales from the International Space Station: evidence for defence system streamlining in a closed habitable environment. *npj Microgravity* (submitted 2026).

---

## License

MIT © 2026 ZJY

## Contact

jiayu6954@gmail.com
