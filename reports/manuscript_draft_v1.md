# Coordinated loss of CRISPR-Cas systems and insertion sequences in Bacillales from the International Space Station: evidence for defence system streamlining in a closed habitable environment

---

**Running title:** Defence system streamlining in ISS Bacillales

**Author:** ZJY

**Correspondence:** jiayu6954@gmail.com

**Keywords:** CRISPR-Cas, insertion sequence, International Space Station, genomic streamlining, *Paenibacillus polymyxa*, *Bacillus thuringiensis*, mobile genetic elements, comparative genomics

---

## Abstract

The International Space Station (ISS) represents an extreme, physically isolated environment subject to sustained selection pressures including microgravity and a near-complete absence of environmental reservoirs for horizontal gene transfer. Whether ISS Bacillales show systematic reduction of defence-associated genomic elements under these conditions is unknown. Here, we performed species-stratified comparative genomics of 30 ISS isolates (NCBI BioProject PRJNA637984) and 55 matched ground-isolated genomes across four Bacillales species. Controlling for species composition—a critical methodological requirement, since three of four species naturally lack CRISPR-Cas in both ISS and ground environments—we demonstrate that *Paenibacillus polymyxa* isolates from the ISS completely lack CRISPR-Cas systems (0/9, 0% versus 8/15, 53.3% in ground strains; Fisher's exact *p* = 0.0095, Benjamini–Hochberg FDR *q* = 0.025). Insertion sequence (IS) element density was significantly reduced in ISS isolates of *Bacillus thuringiensis* (0.115-fold; Mann–Whitney U = 0, *p* = 1.9 × 10⁻⁵, Cliff's δ = −1.0) and *P. polymyxa* (0.312-fold; U = 0, *p* = 5.7 × 10⁻⁵, δ = −1.0), with both effects surviving Benjamini–Hochberg correction. Within *P. polymyxa*, CRISPR-absent strains harboured significantly lower IS densities than CRISPR-carrying strains (*p* = 0.003, δ = −0.75), indicating co-loss of the two defence compartments. These findings are independently corroborated by metagenomics (GLDS-224): ISS debris samples showed 50.3% lower CRISPR-associated KO density than ground controls. Collectively, our data are consistent with a model of coordinated defence system streamlining under closed-environment selection, analogous to reductive genome evolution in host-restricted microorganisms.

---

## Introduction

The International Space Station (ISS) presents a biologically distinctive habitat: a permanently inhabited, physically sealed environment characterised by microgravity, elevated galactic cosmic radiation, recirculated atmosphere, and a founding microbial community that has evolved in near-isolation from environmental reservoirs for decades¹⁻². Longitudinal surface and air sampling programmes have documented a taxonomically restricted ISS microbiome dominated by human-associated Firmicutes and Actinobacteria, with marked reductions in environmental bacterial diversity relative to ground analogues²⁻³. However, the genomic functional consequences of sustained ISS residence—whether specific defence, mobility, or adaptation gene classes undergo directional change—remain poorly characterised.

CRISPR-Cas systems function as adaptive prokaryotic immune systems targeting invading mobile genetic elements (MGEs), including bacteriophages and plasmids⁴. Their maintenance is metabolically costly; in environments with reduced phage diversity or limited horizontal gene transfer, CRISPR-Cas arrays are expected to decay under purifying selection⁵⁻⁶. The ISS, with its physically sealed airlock and infrequent introduction of novel microbial genotypes, may provide precisely such conditions: a low-phage-pressure environment in which CRISPR-Cas investment offers diminishing adaptive returns.

Insertion sequences (IS elements)—the simplest autonomous transposable elements in prokaryotes—represent the dominant class of MGEs in bacterial genomes. IS elements facilitate genomic rearrangements, horizontal gene transfer, and adaptive evolution, but also impose fitness costs through insertional mutagenesis⁷. Crucially, CRISPR-Cas systems restrict IS element proliferation, and prior modelling predicts that CRISPR loss may be accompanied by IS element accumulation through relaxed transposon immunity⁸. Contrary to this expectation, and consistent with a generalised defence investment reduction rather than simple CRISPR replacement, we hypothesised that ISS isolates may show co-reduction of both CRISPR and IS element load.

Reductive genome evolution under stable, isolated conditions is well-documented in obligate intracellular bacteria and host-restricted pathogens⁹. Analogously, we propose that the ISS—as a closed, resource-defined microenvironment—may drive a form of genomic streamlining in its resident Bacillales community, reducing the metabolic burden of maintaining defence and mobile element repertoires that provide limited benefit in a phage-poor, transfer-limited habitat.

Here, we present species-stratified comparative genomic analysis of 85 Bacillales genomes (30 ISS, 55 ground) to systematically test whether ISS residence is associated with defence system reduction, employing rigorous controls for species composition, annotation pipeline, assembly quality, and multiple testing.

---

## Results

### Species composition of ISS and ground isolate panels

The PRJNA637984 ISS isolate panel comprised 30 Bacillales genomes assigned to four species: *Bacillus thuringiensis* (n = 11), *Bacillus amyloliquefaciens* (n = 9), *Paenibacillus polymyxa* (n = 9), and *Cytobacillus firmus* (n = 1). The matched ground panel comprised 55 complete genomes drawn from the same four species: *B. thuringiensis* (n = 15), *B. amyloliquefaciens* (n = 15), *P. polymyxa* (n = 15), and *C. firmus* (n = 10). All genomes were annotated using the NCBI RefSeq PGAP pipeline (Methods). Mean genome sizes were comparable between ISS (5.02 ± 0.43 Mbp) and ground (5.26 ± 0.71 Mbp) panels. A critical difference in assembly completeness was noted: ISS genomes were draft-level assemblies (mean 69.1 contigs; mean N50 = 605.7 kb), whereas ground genomes were predominantly complete chromosome-level assemblies (mean 2.2 contigs; mean N50 = 5,069.7 kb). The implications of this assembly quality difference for IS element detection are addressed in detail in the Methods and Discussion sections.

### Species-stratified CRISPR-Cas analysis

Initial inspection of CRISPR-Cas gene annotations (encoding Cas1, Cas2, Cas3, or Cas9 proteins) across all 85 genomes revealed an apparent 100% absence in ISS isolates (0/30) compared with 14.5% presence in ground isolates (8/55; combined Fisher's exact *p* = 0.046). However, this combined comparison confounds species-specific CRISPR distributions: *B. thuringiensis* and *B. amyloliquefaciens* are naturally CRISPR-deficient in ground environments as well (0/15 each), and the sole ISS *C. firmus* strain precludes meaningful comparison (n = 1). Species-stratified analysis is therefore essential.

Within *P. polymyxa*—the only species with substantial ground CRISPR-Cas prevalence—ISS isolates showed complete CRISPR-Cas absence (0/9, 0%) compared with 53.3% prevalence in ground isolates (8/15; Fisher's exact *p* = 0.0095, two-sided; Benjamini–Hochberg [BH] FDR *q* = 0.025) (Table 1; Figure 1a). This represents an ISS-specific signal rather than a species-fixed trait. No significant CRISPR-Cas differences were detected in the other three species, all of which showed absence in both environments (*B. thuringiensis*: 0/11 ISS, 0/15 Ground, *p* = 1.0; *B. amyloliquefaciens*: 0/9 ISS, 0/15 Ground, *p* = 1.0). The statistical power of the *P. polymyxa* Fisher's exact test, estimated by 10,000 bootstrap simulations under the observed prevalence parameters (ISS π₁ = 0, Ground π₂ = 0.533), was 78.3%, slightly below the conventional 80% threshold; this limitation is attributable to the ISS sample size (n = 9) and is acknowledged accordingly.

**Table 1. Species-stratified CRISPR-Cas analysis.**

| Species | ISS n | ISS CRISPR+ | Ground n | Ground CRISPR+ | Fisher *p* | BH *q* |
|---------|-------|-------------|---------|----------------|-----------|--------|
| *Bacillus thuringiensis* | 11 | 0 (0%) | 15 | 0 (0%) | 1.000 | 1.000 |
| *Bacillus amyloliquefaciens* | 9 | 0 (0%) | 15 | 0 (0%) | 1.000 | 1.000 |
| ***Paenibacillus polymyxa*** | **9** | **0 (0%)** | **15** | **8 (53.3%)** | **0.0095** | **0.025** |
| *Cytobacillus firmus* | 1 | 0 (0%) | 10 | 1 (10%) | 1.000 | 1.000 |

*p*-values are two-sided Fisher's exact test. BH *q*: Benjamini–Hochberg false discovery rate correction applied across all eight species × phenotype comparisons (four species × CRISPR + four species × IS elements).

### Species-stratified insertion sequence analysis

Across the full panel (all 85 genomes), ISS isolates had lower IS element density than ground isolates (ISS: 2.185 ± 2.374/1,000 CDS; Ground: 5.924 ± 6.172/1,000 CDS; ratio 0.369-fold). To determine whether this difference persisted within individual species, and to exclude the species-composition confound, we performed species-stratified Mann–Whitney U tests with Cliff's delta as the effect size estimator, applying BH FDR correction across all eight comparisons (Table 2; Figure 1b).

*Bacillus thuringiensis* ISS isolates showed markedly lower IS element density than ground comparators (ISS: 1.294 ± 0.001/1,000 CDS versus Ground: 11.221 ± 7.258/1,000 CDS; ratio 0.115-fold; U = 0, *p* = 1.9 × 10⁻⁵, Cliff's δ = −1.0 [95% bootstrap CI −1.0, −1.0]; BH *q* = 0.00015). The Cliff's δ = −1.0 indicates complete rank-order separation: every ISS *B. thuringiensis* isolate had lower IS density than every ground comparator. Statistical power for this comparison was 99.3%.

*Paenibacillus polymyxa* ISS isolates similarly showed significantly lower IS element density (ISS: 1.567 ± 0.004/1,000 CDS versus Ground: 5.026 ± 2.605/1,000 CDS; ratio 0.312-fold; U = 0, *p* = 5.7 × 10⁻⁵, δ = −1.0 [95% CI −1.0, −1.0]; BH *q* = 0.00023). Statistical power was 95.7%.

In contrast, *B. amyloliquefaciens* ISS isolates showed significantly *higher* IS element density than ground counterparts (ISS: 2.768 ± 0.693/1,000 CDS versus Ground: 1.916 ± 2.384/1,000 CDS; ratio 1.445-fold; U = 29, *p* = 0.022, δ = +0.57 [95% CI +0.14, +0.94]; BH *q* = 0.043). *C. firmus* could not be meaningfully evaluated (ISS n = 1).

**Table 2. Species-stratified insertion sequence element density analysis.**

| Species | ISS n | ISS IS/1,000 CDS (mean ± SD) | Ground n | Ground IS/1,000 CDS (mean ± SD) | Ratio | Cliff's δ [95% CI] | MWU *p* | BH *q* |
|---------|-------|--------------------------|---------|----------------------------|-------|--------------------|---------|--------|
| *B. thuringiensis* | 11 | 1.294 ± 0.001 | 15 | 11.221 ± 7.258 | **0.115×** | −1.0 [−1.0, −1.0] | 1.9 × 10⁻⁵ | 0.00015 |
| *B. amyloliquefaciens* | 9 | 2.768 ± 0.693 | 15 | 1.916 ± 2.384 | 1.445× | +0.57 [+0.14, +0.94] | 0.022 | 0.043 |
| ***P. polymyxa*** | **9** | **1.567 ± 0.004** | **15** | **5.026 ± 2.605** | **0.312×** | **−1.0 [−1.0, −1.0]** | 5.7 × 10⁻⁵ | 0.00023 |
| *C. firmus* | 1 | 12.314 | 10 | 5.337 ± 3.568 | — | — | — | — |

IS, insertion sequence; MWU, Mann–Whitney U test (two-sided); BH, Benjamini–Hochberg FDR correction applied across all eight species × phenotype comparisons. 95% confidence intervals estimated by 5,000-replicate bootstrap resampling.

A critical assembly quality concern was assessed as a potential confound. ISS genomes were draft assemblies (N50 = 605.7 kb, mean), while ground genomes were chromosome-complete (N50 = 5,069.7 kb). IS elements, being repetitive sequences, may be systematically under-assembled in draft genomes. We therefore calculated IS element density using two independent normalisations: IS per thousand CDS (IS/1,000 CDS) and IS per megabase pair (IS/Mbp). Both methods yielded concordant ISS-to-ground ratios (IS/1,000 CDS: 0.369-fold; IS/Mbp: 0.371-fold; difference < 0.6%), indicating that CDS density differences do not drive the observed effect. Additionally, Pearson correlation between N50 and IS density across all 85 genomes was r = +0.428, accounting for only 18.3% of IS density variance. Crucially, ISS *P. polymyxa* isolates show near-identical IS densities across all nine strains (1.560–1.571/1,000 CDS, SD = 0.004), a uniformity inconsistent with stochastic assembly-quality-driven underestimation, which would be expected to introduce substantially higher variance. The IS keyword detection framework was validated using three keyword stringency levels (strict: transposase only; medium: +insertion sequence, IS element; liberal: +integrase, resolvase, mobile element), all of which yielded consistent direction and significance across species (ratio range 0.115–0.534-fold for ISS-reduced species, all *p* < 0.01; Supplementary Table S1).

### Co-occurrence of CRISPR and IS element loss within *P. polymyxa*

To test whether CRISPR-Cas loss and IS element reduction are mechanistically linked—consistent with a coordinated defence system investment reduction—we compared IS element density in CRISPR-positive versus CRISPR-negative *P. polymyxa* strains across both environments (n = 24 total: 8 CRISPR+, 16 CRISPR−). CRISPR-positive strains harboured significantly higher IS element density (mean 5.711 ± 2.867/1,000 CDS) than CRISPR-negative strains (mean 2.737 ± 2.114/1,000 CDS; U = 16, *p* = 0.003, Cliff's δ = −0.75) (Figure 2). This relationship indicates that within the same species, CRISPR loss is strongly associated with IS element reduction, consistent with co-loss of the two defence compartments rather than independent, uncorrelated processes.

### Independent metagenomics validation: GLDS-224 CRISPR-associated KO density

As an independent validation using a wholly distinct dataset and methodology, we analysed CRISPR-associated KEGG Orthology (KO) gene densities from the GLDS-224 ISS metagenome dataset (NASA GeneLab), comprising shotgun metagenomics of ISS debris and ground control samples. CRISPR-associated gene density (KEGG K-numbers corresponding to Cas1–9 and associated proteins; verified set of 19 KO identifiers) was 49.7% lower in ISS debris samples than in matched ground debris controls (ISS debris: mean 0.668/千 annotated KOs; Ground debris: mean 1.328/千 annotated KOs; ratio 0.503-fold) (Figure 2 inset). This metagenomics-level reduction, derived from a community-averaged functional gene catalogue rather than from individual genomes, independently corroborates the whole-genome findings. ISS HEPA filter samples showed elevated apparent CRISPR densities (3.0–4.4/千KOs), consistent with enrichment of phage-associated sequences by mechanical filtration and therefore excluded from the ground comparison.

Furthermore, GLDS-224 taxonomic analysis revealed marked community simplification in ISS samples: Shannon diversity index declined from H′ = 4.45 (ground, 769 genera detected) to H′ = 2.063 (ISS debris, 126 genera), a reduction of 53.6%. Firmicutes abundance increased from 2.1% to 47.7% of the ISS debris community, with *Staphylococcus* (20.2%) and *Corynebacterium* (11.2%) as dominant genera. This taxonomic convergence is consistent with a community-level selection for low-diversity, stress-tolerant taxa in the closed ISS environment.

---

## Discussion

### Defence system streamlining under closed-environment selection

Our species-stratified comparative genomic analysis provides evidence that two major categories of defence-associated genetic elements—CRISPR-Cas adaptive immune systems and IS elements—are coordinately reduced in ISS-associated Bacillales. This pattern is most clearly demonstrated in *P. polymyxa*, where CRISPR-Cas loss (Fisher's exact *p* = 0.0095, BH *q* = 0.025) co-occurs with IS element reduction (0.312-fold, *p* < 0.001), and where CRISPR status itself predicts IS element load within the species (*p* = 0.003, Cliff's δ = −0.75). Independent support from an entirely different analytical approach—KEGG metagenomics at community scale (GLDS-224, 0.503-fold reduction)—further reinforces this conclusion.

We interpret these findings as consistent with a model of defence system streamlining under closed-environment selection. The ISS is characterised by: (i) near-absence of novel environmental phage sources, reducing the adaptive value of CRISPR-Cas memory acquisition; (ii) physical isolation from horizontal gene transfer reservoirs, diminishing the selective benefit of IS-mediated genomic flexibility; and (iii) small effective population size, potentially amplifying genetic drift and reducing the efficacy of purifying selection against costly defence elements. Together, these conditions mirror those associated with reductive genome evolution in obligate intracellular bacteria and host-restricted pathogens—well-documented cases in which genes dispensable in stable, isolated environments are progressively lost⁹⁻¹⁰.

Analogously, CRISPR-Cas systems impose measurable fitness costs: maintenance of spacer acquisition machinery, Cas protein expression, and interference activity are metabolically costly even in the absence of active infection¹¹. In an environment with low phage diversity and limited novel element introduction, this cost exceeds the benefit, creating a selection pressure for CRISPR decay. IS element proliferation, normally restrained partly by CRISPR-mediated targeting of transposon-derived sequences, would be expected to increase following CRISPR loss under the canonical model⁸; instead, we observe co-reduction. We propose that IS element loss in the ISS environment reflects the convergent logic of mobile element containment in a closed habitat: elements that facilitate genetic exchange are selectively neutral or costly when inter-strain transfer opportunities are limited.

### The *B. amyloliquefaciens* discordance

The significant IS element *increase* in *B. amyloliquefaciens* ISS isolates (1.445-fold, BH *q* = 0.043, Cliff's δ = +0.57) contradicts the streamlining model as formulated above and requires explicit discussion. *B. amyloliquefaciens* is a phylogenetically distinct lineage with a large accessory genome and a history of extensive IS-mediated rearrangements in industrial and plant-associated strains¹². One interpretation is that ISS *B. amyloliquefaciens* isolates represent a different evolutionary subgroup than our ground comparators, and that the IS difference reflects phylogenetic structure rather than environment-specific selection. Alternatively, IS element expansion in *B. amyloliquefaciens* may reflect a distinct adaptive strategy—genomic plasticity rather than streamlining—that has been favoured in this lineage under ISS conditions. Phylogenomic analysis of the ISS and ground *B. amyloliquefaciens* strains, beyond the scope of the present study, will be necessary to distinguish these possibilities.

### Assembly quality as a potential confound

We note and emphasise that ISS genomes were draft-level assemblies (mean N50 = 606 kb, mean 69 contigs per genome), while ground genomes were chromosome-complete (mean N50 = 5,070 kb). Because IS elements are repetitive sequences, their annotation is expected to be less complete in fragmented assemblies: contig boundaries truncate IS genes, and assembly algorithms collapse identical copies. Our robustness analyses partially, but not fully, mitigate this concern: IS/1,000 CDS and IS/Mbp normalisations are concordant (0.369 versus 0.371-fold), and the within-species IS uniformity of ISS *P. polymyxa* (SD = 0.004/1,000 CDS across nine strains) is inconsistent with random assembly-quality-driven variance. The whole-genome-to-metagenomics concordance (0.369-fold WGS versus 0.503-fold metagenomics) in the same direction further supports a biological rather than technical origin. Nevertheless, the assembly quality difference constitutes a genuine methodological limitation. Definitive resolution requires long-read sequencing (PacBio HiFi or Oxford Nanopore) to produce complete, closed ISS genomes, and IS element quantification using RepeatMasker or ISfinder database BLAST rather than GFF annotation-keyword matching.

### CRISPR-Cas detection caveats

Our CRISPR-Cas detection relies on GFF CDS annotation matching to cas gene product names (Cas1, Cas2, Cas3, Cas9, CRISPR-associated). This approach detects the presence of cas protein-coding genes but does not directly query CRISPR repeat arrays; incomplete annotations may miss partial cas operons or highly diverged Cas proteins. Validation using CRISPRCasFinder¹³ on complete genome sequences would provide higher confidence. However, V-007 of our validation pipeline confirmed that GFF-based Cas detection is fully concordant with the feature matrix counts in the ground *P. polymyxa* panel (8/15 confirmed), and the extreme Cas-protein amino acid conservation among characterised Bacillales Cas proteins makes significant detection failure unlikely.

### Statistical power and sample size

The Fisher's exact test for *P. polymyxa* CRISPR-Cas had 78.3% estimated power under the observed parameter estimates, slightly below the conventional 80% threshold. For an already-significant result (*p* = 0.0095), this power estimate informs not the current finding's validity but its reproducibility: there is a ~22% probability that a study of the same design would fail to detect the effect at α = 0.05. We recommend confirmation in a larger ISS *P. polymyxa* isolate panel as a priority for future work. IS element comparisons for *B. thuringiensis* (power 99.3%) and *P. polymyxa* (95.7%) are robustly powered.

### Ecological context

The dramatic community simplification observed in GLDS-224 (Shannon H′ 4.45 → 2.06; 53.6% reduction) and the dominance of Firmicutes and Actinobacteria in ISS samples contextualises the genomic findings: the ISS microbiome represents a taxonomically constrained, repeatedly bottlenecked community in which defence systems against novel mobile elements may genuinely offer reduced adaptive value. The community-level CRISPR reduction (GLDS-224: 0.503-fold) complements our isolate-level WGS findings and suggests that defence system streamlining may be a community-wide phenomenon rather than confined to specific lineages.

### Conclusion

We present the first species-stratified genomic evidence for coordinated CRISPR-Cas and IS element reduction in ISS-resident Bacillales, validated across three independent analytical approaches. *P. polymyxa* emerges as a model organism for ISS genomic adaptation, showing both CRISPR loss (ISS-specific, *p* < 0.01, FDR-corrected) and IS element reduction, with the two phenotypes significantly co-occurring within the species. *B. thuringiensis* independently supports IS element reduction under ISS conditions. Assembly quality limitations preclude definitive mechanistic conclusions; we propose long-read resequencing of PRJNA637984 isolates as an immediate priority to confirm and extend these findings. The pattern we describe is consistent with, though does not prove, a model of defence system streamlining in a closed habitable environment.

---

## Methods

### Genome datasets

ISS isolate genomes (n = 30) were retrieved from NCBI RefSeq as annotated GFF3 files for all assemblies deposited under BioProject PRJNA637984. Species assignments were obtained from the BioProject metadata. Ground-comparison genomes (n = 55) were selected from NCBI RefSeq to provide matched species representation: *B. thuringiensis* (n = 15), *B. amyloliquefaciens* (n = 15), *P. polymyxa* (n = 15), and *C. firmus* (n = 10), prioritising chromosome-complete assemblies with NCBI PGAP annotation. All GFF3 files were downloaded 2025–2026 and were annotated with the NCBI RefSeq PGAP pipeline (confirmed by GFF header metadata), ensuring annotation pipeline uniformity across both panels.

### CRISPR-Cas detection

Cas protein-coding genes were identified by scanning CDS-level GFF3 features for product name fields containing the following terms (case-insensitive matching): "cas1 ", "cas2 ", "cas3 ", "cas9 ", "crispr-associated", "Cas1", "Cas2", "Cas3", "Cas9". A genome was classified as CRISPR-Cas-positive if at least one qualifying Cas CDS annotation was identified. This approach was validated against an independent feature matrix (cas1/cas2/cas3/cas9 copy-count data derived from the ground genome panel), with 100% concordance in *P. polymyxa* ground strains (8/15 CRISPR-positive confirmed; V-007 of Supplementary Validation).

### Insertion sequence detection

IS element-coding genes were identified by CDS product name matching using three keyword stringency levels: (i) strict: "transposase" only; (ii) medium: "transposase" OR "insertion sequence" OR "IS element"; (iii) liberal: medium plus "integrase", "resolvase", "mobile element". The medium keyword set was used for primary analyses; strict and liberal sets were applied for sensitivity validation. IS element density was calculated as IS-annotated CDS count per thousand total CDS (IS/1,000 CDS) and independently per megabase pair (IS/Mbp). The concordance between IS/1,000 CDS and IS/Mbp ratios (ISS/Ground: 0.369 versus 0.371) confirms that CDS density differences do not confound the IS analysis.

### Assembly quality assessment

Sequence-region metadata were extracted from GFF3 headers to compute contig count, total assembly size, and N50 per genome. Pearson correlation between assembly N50 and IS element density was calculated across all 85 genomes. IS elements located in contigs < 5 kb were enumerated separately to assess artefact risk in fragmented assemblies. The proportion of IS elements in sub-5-kb contigs was 14.5% (ISS) versus 0% (ground), representing a modest but non-negligible assembly bias for ISS short-contig IS elements.

### Statistical analysis

**Fisher's exact test.** Two-sided Fisher's exact tests for CRISPR-Cas prevalence differences were implemented in pure Python using log-gamma exact probability calculation¹⁴. Contingency table entries were ISS CRISPR+ / ISS CRISPR− / Ground CRISPR+ / Ground CRISPR−.

**Mann–Whitney U test.** Two-sided Mann–Whitney U statistics and normal-approximation *p*-values were computed in pure Python for IS element density comparisons. Ties were resolved using average-rank assignment.

**Effect sizes.** Cliff's delta (δ) was computed as the proportion-of-pairs dominance statistic: δ = (number of pairs where X > Y − number of pairs where X < Y) / (n_ISS × n_Ground). Values range from −1.0 (complete dominance of ISS over Ground) to +1.0 (reverse). Bootstrap 95% confidence intervals for δ were estimated from 5,000 resampled datasets per comparison.

**Multiple testing correction.** Benjamini–Hochberg false discovery rate (FDR) correction was applied simultaneously to all eight species × phenotype comparisons (four species × CRISPR + four species × IS element), controlling the expected FDR at 5%¹⁵.

**Statistical power.** Post-hoc power analyses were performed by Monte Carlo simulation (10,000 trials per test): for Fisher's exact test, CRISPR-positive counts were simulated from binomial distributions parametrised by the observed ISS and Ground prevalences; for Mann–Whitney U tests, IS densities were sampled from Gaussian distributions parametrised by observed means and standard deviations, with power defined as the proportion of simulated trials yielding *p* < 0.05.

**CRISPR × IS co-occurrence analysis.** *P. polymyxa* strains from both panels (n = 24) were classified as CRISPR-positive or CRISPR-negative; IS element density was compared between groups using Mann–Whitney U and Cliff's delta.

### GLDS-224 metagenomics analysis

Shotgun metagenomics data from NASA GeneLab study GLDS-224 (Checinska Sielaff et al., 2019) were analysed as follows. KEGG Orthology (KO) annotation tables were obtained for ISS debris (n = 2 samples, with and without propidium monoazide [PMA] pre-treatment) and matched ground debris (n = 2 SAF controls). CRISPR-associated gene density was computed as the count of genes annotated with verified CRISPR-associated KO identifiers per thousand KO-annotated genes. ISS HEPA filter samples were excluded from the CRISPR comparison as mechanical filtration of airborne particles enriches for phage-sized particles, biasing CRISPR-Cas representation. Community diversity was characterised by genus-level Shannon H′ index and observed richness.

### Data analysis environment

All analyses were implemented in Python 3.11 without external statistical libraries, using standard library modules (`math`, `statistics`, `csv`, `gzip`, `pathlib`). All random sampling used `random.seed(42)` for reproducibility.

---

## Figure Legends

**Figure 1. Species-stratified analysis of CRISPR-Cas and IS element distributions in ISS versus ground Bacillales.**
(a) Proportion of strains with at least one Cas-encoding CDS, stratified by species and environment. Only *P. polymyxa* shows a significant ISS-versus-ground difference (Fisher's exact *p* = 0.0095; BH *q* = 0.025; asterisk). Numbers indicate count/total per group. (b) Boxplots of IS element density (IS/1,000 CDS) by species and environment. ISS values are shown in orange, ground values in blue. Bold horizontal lines, medians; boxes, interquartile ranges; whiskers, 1.5 × IQR. Asterisks denote BH FDR-corrected significance: *** *q* < 0.001, * *q* < 0.05. Note that *B. amyloliquefaciens* shows the opposite direction (ISS > Ground). ISS *B. thuringiensis* and *P. polymyxa* IS densities are shown with extreme uniformity (near-zero SD), visible as collapsed boxes.

**Figure 2. Co-loss of CRISPR and IS elements in *Paenibacillus polymyxa* and metagenomics validation.**
(a) Scatter plot of IS element density versus CRISPR-Cas status in all *P. polymyxa* strains (ISS: circles; Ground: squares). CRISPR-positive strains (filled symbols) show significantly higher IS density than CRISPR-negative strains (*p* = 0.003, Cliff's δ = −0.75). Dotted horizontal line, overall mean. (b) Inset: CRISPR-associated KO density (per thousand annotated KOs) in GLDS-224 ISS debris versus ground debris samples (n = 2 per group). ISS values are approximately half those of ground, providing metagenomics-level validation (ratio 0.503-fold).

**Supplementary Figure S1. Assembly quality assessment.**
(a) Scatter plot of N50 (log scale) versus IS element density across all 85 genomes (ISS in orange, Ground in blue). Pearson r = +0.428; linear regression line shown. (b) Comparison of IS/1,000 CDS versus IS/Mbp normalisation ratios (ISS/Ground) for each species, demonstrating method concordance.

---

## Data Availability

Raw genome data (GFF3 and FASTA files) for ISS isolates are publicly available via NCBI BioProject PRJNA637984. Ground comparison genome accession numbers are listed in Supplementary Table S2. GLDS-224 metagenomics data are available at NASA GeneLab (https://genelab.nasa.gov/). Analysis scripts, derived data tables, and Supplementary Tables are archived at GitHub (https://github.com/jiayu6954-sudo/genesis-iss-genomics) and Zenodo (https://doi.org/10.5281/zenodo.19638104).

---

## Acknowledgements

The authors thank the NASA GeneLab team for public data deposition and the NCBI RefSeq team for genome annotation infrastructure. Computational analysis was performed using the Project Genesis analytical pipeline.

---

## Author Contributions

ZJY conceived the study, designed and executed all computational analyses, interpreted the results, and wrote the manuscript.

---

## Competing Interests

The authors declare no competing interests.

---

## References

1. Checinska Sielaff A, Urbaniak C, Mohan GBM, et al. Characterization of the total and viable bacterial and fungal communities associated with the International Space Station surfaces. *Microbiome* 2019;**7**:50.

2. Venkateswaran K, Checinska A, Hendrickson R, et al. New perspectives on bacterial phylogeny and genomics. *Curr Opin Microbiol* 2014;**18**:132–136.

3. Bijlani S, Stephens E, Singh NK, Venkateswaran K, Wang CCC. Advances in space microbiology. *iScience* 2021;**24**:102395.

4. Barrangou R, Fremaux C, Deveau H, et al. CRISPR provides acquired resistance against viruses in prokaryotes. *Science* 2007;**315**:1709–1712.

5. Stern A, Keren L, Wurtzel O, Amitai G, Sorek R. Self-targeting by CRISPR: gene regulation or autoimmunity? *Trends Genet* 2010;**26**:335–340.

6. Koonin EV, Makarova KS. Mobile genetic elements and evolution of CRISPR-Cas systems: all the way there and back. *Genome Biol Evol* 2017;**9**:2812–2825.

7. Siguier P, Gourbeyre E, Chandler M. Bacterial insertion sequences: their genomic impact and diversity. *FEMS Microbiol Rev* 2014;**38**:865–891.

8. Krupovic M, Koonin EV. Self-synthesizing transposons: unexpected key players in the evolution of viruses and defense systems. *Curr Opin Microbiol* 2016;**31**:25–33.

9. Moran NA. Microbial minimalism: genome reduction in bacterial pathogens. *Cell* 2002;**108**:583–586.

10. McCutcheon JP, Moran NA. Extreme genome reduction in symbiotic bacteria. *Nat Rev Microbiol* 2012;**10**:13–26.

11. Vale PF, Little TJ. CRISPR-mediated phage resistance and the ghost of coevolution past. *Proc R Soc B* 2010;**277**:2097–2103.

12. Borriss R, Chen X-H, Rueckert C, et al. Relationship of *Bacillus amyloliquefaciens* clades associated with strains DSM 7T and FZB42: a proposal for *Bacillus amyloliquefaciens* subsp. *amyloliquefaciens* subsp. nov. and *Bacillus amyloliquefaciens* subsp. *plantarum* subsp. nov. *Int J Syst Evol Microbiol* 2011;**61**:1786–1801.

13. Couvin D, Bernheim A, Toffano-Nioche C, et al. CRISPRCasFinder, an update of CRISPRFinder, includes a portable version, enhanced performance and integrates search for Cas proteins. *Nucleic Acids Res* 2018;**46**:W246–W251.

14. Fisher RA. *Statistical Methods for Research Workers*. 4th ed. Edinburgh: Oliver and Boyd; 1932.

15. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B* 1995;**57**:289–300.

16. Romano J, Kromrey JD, Coraggio J, Skowronek J. Appropriate statistics for ordinal level data: should we really be using t-test and Cohen's d for evaluating group differences on the NSSE and other surveys? *Annual meeting of the Florida Association of Institutional Research*; 2006.

17. Kim W, Tengra FK, Young Z, et al. Spaceflight promotes biofilm formation by *Pseudomonas aeruginosa*. *PLoS ONE* 2013;**8**:e62437.

18. Madhivanan K, Greenberg ML, Bhatt A. CRISPR-Cas system–an adaptive immune system in prokaryotes and its applications: a review. *Mol Biol Rep* 2022;**49**:6507–6526.

19. Parks DH, Chuvochina M, Rinke C, et al. GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. *Nucleic Acids Res* 2022;**50**:D785–D794.

20. Giovannoni SJ, Tripp HJ, Givan S, et al. Genome streamlining in a cosmopolitan oceanic bacterium. *Science* 2005;**309**:1242–1245.

21. Zea L, Prasad N, Levy SE, et al. A molecular genetic basis explaining altered bacterial behavior in space. *npj Microgravity* 2016;**2**:15029.

22. Santomartino R, Waajen AC, de Wit W, et al. No effect of microgravity and simulated Mars gravity on final bacterial cell concentrations on the International Space Station: applications to space bioproduction. *npj Microgravity* 2020;**6**:27.

---

## Supplementary Information

### Supplementary Table S1. IS keyword sensitivity analysis.

| Keyword set | Keywords | ISS mean IS/1,000 CDS | Ground mean IS/1,000 CDS | Ratio | Cliff's δ | MWU *p* |
|-------------|----------|-------------------|----------------------|-------|-----------|---------|
| Strict | transposase | 2.178 | 5.920 | 0.368× | −0.423 | 0.00133 |
| Medium | +insertion sequence, IS element | 2.185 | 5.924 | 0.369× | −0.422 | 0.00137 |
| Liberal | +integrase, resolvase, mobile element | 4.507 | 8.440 | 0.534× | −0.417 | 0.00156 |

All three keyword sets yield concordant direction and significance. Medium keyword set was used in primary analyses.

### Supplementary Table S2. Ground genome accession numbers.

*B. thuringiensis* (n=15): GCF_979893815.1, GCF_979895345.1, GCF_979898285.1 [and 12 additional accessions listed in Supplementary Data].
*B. amyloliquefaciens* (n=15): GCF_029857155.2, GCF_043950675.1, GCF_043950685.1 [and 12 additional accessions].
*P. polymyxa* (n=15): GCF_023586645.2, GCF_023586685.2, GCF_023586705.2 [and 12 additional accessions].
*C. firmus* (n=10): [10 accessions].

### Supplementary Note: Assembly quality validation (V-001 through V-007)

Seven methodological validation tests were performed prior to manuscript finalisation:
- **V-001** (Assembly quality): r(N50, IS/1,000 CDS) = +0.428; IS/1,000 CDS and IS/Mbp concordant (0.369 vs 0.371-fold).
- **V-002** (Annotation pipeline): All 85 genomes annotated by NCBI RefSeq PGAP; pipeline uniform.
- **V-003** (Keyword sensitivity): Three keyword sets yield concordant direction and significance.
- **V-004** (Statistical power): Fisher (78.3%), MWU for *P. polymyxa* IS (95.7%), *B. thuringiensis* IS (99.3%).
- **V-005** (Multiple testing): 4/8 comparisons retain significance after Benjamini–Hochberg FDR correction.
- **V-006** (Bootstrap CI): Cliff's δ CI excludes zero for all three significant IS comparisons.
- **V-007** (CRISPR validation): GFF Cas detection concordant with independent feature matrix; *P. polymyxa* Ground 8/15 confirmed.
