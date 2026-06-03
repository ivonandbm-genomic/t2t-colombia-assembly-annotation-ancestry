#  High-quality HiFi de novo assembled and annotated diploid genomes from admixed Colombian individuals

---

## **Overview**

This repository contains the complete bioinformatics workflow for:

- **De novo genome assembly (PacBio HiFi)**
- **Genome annotation**
- **Small variant detection**
- **Global and local ancestry inference**

The analysis was performed on two admixed Colombian individuals from the Valle del Cauca region using long-read sequencing and telomere-to-telomere (T2T) references.

---
## Study Design

To simplify figure labeling and downstream analyses, the following sample aliases were used throughout the repository and manuscript:

| Original sample ID | Alias |
|---|---|
| 017 / 017C | COLM |
| 018 / 018C | COLF |

These aliases are consistently used in figures, population structure analyses, local ancestry plots, and manuscript visualizations.

- **Samples:** 2 individuals (COLM and COLF)
- **Technology:** PacBio HiFi (CCS reads, Q20+)
- **Coverage:** ~30×
- **Reference genomes:**
  - GRCh38
  - T2T-CHM13 v2.0

---

##  **Data Availability**

Due to ethical and privacy constraints, raw genomic data (**FASTQ, BAM, VCF**) are not publicly available.

This repository provides:

- ✔️ Reproducible pipelines  
- ✔️ Scripts used in the analysis  
- ✔️ Summary results and figures  

---

##  **Software Requirements**

- **samtools**
- **Jellyfish v2.3.1**
- **GenomeScope 2.0**
- **hifiasm v0.16.0**
- **gfatools**
- **RagTag v2.0.1**
- **QUAST v5.3.0**
- **BUSCO v5.8.1**
- **Merqury v1.3 + Meryl**
- **SURVIVOR**
-  **pbsv**
- **Sniffles2**
- **NanoVar**
- **cuteSV**
- **SVIM-asm**
- **Sawfish**
- **DeepVariant v1.8.0**
- **bcftools**
- **WhatsHap v2.3**
- **PLINK v1.9**
- **ADMIXTURE v1.3**
- **RFMix v1.5.4**
- **RepeatMasker**
- **RepeatModeler**
- **Liftoff v1.6.3**
  
- 

---
## Repository Structure

```text
scripts/
├── 01_assembly_annotation.sh
├── 02_Ancestry.sh
└── README.md
```

---

## Pipeline Organization

The repository is organized into two main modular pipelines:

### 1. Assembly and Annotation Pipeline

Implemented as a single modular script containing the following analytical sections:

1. BAM preprocessing
2. FASTQ extraction
3. k-mer profiling using Jellyfish and GenomeScope
4. De novo assembly using hifiasm
5. GFA-to-FASTA conversion
6. HiFi read alignment
7. RagTag correction
8. RagTag scaffolding
9. Assembly quality assessment
10. BUSCO completeness analysis
11. Merqury consensus quality evaluation
12. Repeat annotation using RepeatMasker and RepeatModeler
13. Gene annotation transfer using Liftoff
14. Structural variant detection using pbsv
15. Small variant calling using Deepvariant

Each module is internally documented and can be executed independently depending on the analysis stage.

---

### 2. Variant, Population Structure, and Local Ancestry Pipeline

Implemented as a single modular script containing the following analytical sections:

1. Variant filtering and preprocessing
2. Reference panel construction and harmonization
3. PLINK harmonization and quality control
4. Principal component analysis (PCA)
5. Global ancestry inference using ADMIXTURE
6. Variant phasing using WhatsHap
7. Chromosome-wise VCF splitting
8. Preparation of RFMix reference and query files
9. Local ancestry inference using RFMix
10. Population structure visualization
11. Local ancestry plotting

Each section is internally documented and can be executed independently depending on the analysis stage.

## Pipeline Description 


### 1. Read Processing

Input:

* PacBio HiFi BAM files

Steps:

* BAM indexing (pbindex)
* BAM to FASTQ conversion (bam2fastq)

Output:

* FASTQ files

---

### 2. Genome Profiling

Tools:

* Jellyfish
* GenomeScope

Purpose:

* Estimate genome size
* Estimate heterozygosity
* Estimate repeat content

Output:

* K-mer histograms
* GenomeScope reports

---

### 3. De Novo Genome Assembly

Tool:

* Hifiasm

Input:

* PacBio HiFi FASTQ files

Output:

* Primary contigs
* Haplotype-resolved assemblies
* GFA files

---

### 4. Conversion of Assemblies

Tool:

* Gfatools

Purpose:

* Convert GFA assemblies to FASTA format

Output:

* Contig FASTA files

---

### 5. Assembly Refinement

Tool:

* RagTag

Functions:

* Structural correction
* Reference-guided scaffolding

References:

* T2T-CHM13 v2.0
* GRCh38

Output:

* Chromosome-scale scaffold assemblies

---

### 6. Assembly Quality Assessment

Tools:

#### QUAST

Evaluates:

* N50
* Assembly size
* Misassemblies
* Genome coverage

#### BUSCO

Evaluates:

* Conserved ortholog completeness

Database:

* primates_odb12

#### Merqury

Evaluates:

* Consensus quality value (QV)
* Assembly completeness

---

### 7. Repeat Annotation

Tool:

* RepeatMasker

Features:

* LINEs
* SINEs
* LTR elements
* DNA transposons
* Simple repeats

Output:

* GFF annotations
* Repeat statistics

---

### 8. Gene Annotation

Tool:

* Liftoff

Reference Annotation:

* T2T-CHM13 gene models

Output:

* GTF annotations
* Lifted genes
* Unmapped genes

---

### 9. Structural Variant Discovery

Tools:

- pbsv
- Sniffles2
- NanoVar
- cuteSV
- SVIM-asm
- Sawfish

Workflow:

1. Collect SV calls from individual callers.
2. Convert compressed VCFs to standard VCF format.
3. Generate SURVIVOR input files.
4. Merge SV calls across callers.
5. Retain variants supported by at least two methods.
6. Filter final calls according to:


## Output

Merged SV VCFs:

```text
SVs_merge_4/
```

Filtered consensus SV VCFs:

```text
SVs_merge_4/filtered/

### 10. Small Variant Calling

Tool:

* DeepVariant

Model:

* PACBIO

Detected Variants:

* SNPs
* Indels

Output:

* VCF
* gVCF


###  11. Global and Local Ancestry Inference

The workflow integrates:

- Whole-genome sequencing (WGS) variant data
- A harmonized global reference panel (HGDP + 1000 Genomes)
- Population structure analysis (PCA)
- Global ancestry estimation (ADMIXTURE)
- Local ancestry inference (RFMix)
- Integration with structural variation and T2T genome coordinates

#### Input Data

- Phased VCF files (biallelic SNPs, autosomes)
- Reference panel: gnomAD HGDP + 1000 Genomes (v3)
- Genetic maps 
- Population labels:
  - AFR
  - AMR
  - EUR
  - EAS



---

### Variant Filtering

Retained variants:

- Autosomal variants (chr1–22)
- PASS variants only
- Biallelic SNPs
- MAF ≥ 0.05
- Missingness ≤ 5%

---

### Reference Panel Construction

A balanced reference panel was constructed using 1,200 unrelated individuals from HGDP + 1000 Genomes:

- AFR (n = 300)
- AMR (n = 300)
- EUR (n = 300)
- EAS (n = 300)

The AMR panel includes:

#### 1000 Genomes populations

- PEL
- MXL
- CLM
- PUR

#### Indigenous American populations (HGDP)

- Maya
- Pima
- Karitiana
- Suruí
- Colombians

---

### PLINK Processing for PCA and ADMIXTURE

The filtered autosomal SNP datasets from the gnomAD HGDP+1KG reference panel and the  Colombian genomes were converted to PLINK format.

Variant IDs were harmonized, common SNPs were retained, and both datasets were merged.


Quality control was performed using:

- Genotype missingness filtering
- MAF filtering
- LD pruning

#### Data Harmonization

- Intersection of SNPs between target genomes and reference panel
- Strict allele matching
- Chromosome-wise processing

---

### LD Pruning

Performed using PLINK:

```bash
plink --indep-pairwise 50 10 0.2
```

### Principal Component Analysis (PCA)

Performed on the LD-pruned dataset.

Main analyses:

- Population structure visualization
- First 20 PCs retained

---

### Global Ancestry Inference (ADMIXTURE)

ADMIXTURE was run for K = 2–6.

Main analyses:

- Three independent runs per K
- Cross-validation used to determine optimal K
- Estimation of ancestry proportions

#### Output

- Ancestry proportions per individual
- Population structure plots

---

### Local Ancestry Inference (RFMix)

Software:

- RFMix v1.5.4

Mode:

- PopPhased

Phase correction:

- Enabled

#### Parameters

- Generations since admixture: 8
- Trees per window: 100
- Window size: 0.2 cM
- Seed: 12345

#### Input

- Phased VCFs (autosomal SNPs)
- Reference haplotypes (HGDP + 1000 Genomes)

#### Processing

- Chromosome-wise analysis
- Harmonization of SNP positions and alleles

Final dataset:

- 1,199 reference samples
- 1 target genome per run

#### Output

- Local ancestry tracts (cM resolution)
- Per-chromosome ancestry assignments
---

### Output Files

- PCA plots (population structure)
- ADMIXTURE results (global ancestry proportions)
- Local ancestry tracts (RFMix)
- Per-chromosome ancestry summaries (bp and cM)

- 
##  References – Genome Assembly

Cheng, H., Concepcion, G. T., Feng, X., Zhang, H., & Li, H. (2021).  
Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm.  
Nature Methods, 18(2), 170–175.  
https://doi.org/10.1038/s41592-020-01056-5  

Nurk, S., Koren, S., Rhie, A., Rautiainen, M., Bzikadze, A. V., Mikheenko, A.,  
Vollger, M. R., Altemose, N., Uralsky, L., Gershman, A., Aganezov, S.,  
Hoyt, S. J., Diekhans, M., Logsdon, G. A., Alonge, M., Antonarakis, S. E.,  
Borchers, M., Bouffard, G. G., Brooks, S. Y., … Phillippy, A. M. (2022).  
The complete sequence of a human genome.  
Science, 376(6588), 44–53.  
https://doi.org/10.1126/science.abj6987  

Li, H. (2018).  
Minimap2: Pairwise alignment for nucleotide sequences.  
Bioinformatics, 34(18), 3094–3100.  
https://doi.org/10.1093/bioinformatics/bty191  

Rhie, A., Walenz, B. P., Koren, S., & Phillippy, A. M. (2020).  
Merqury: Reference-free quality, completeness, and phasing assessment for genome assemblies.  
Genome Biology, 21, 245.  
https://doi.org/10.1186/s13059-020-02134-9  

Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013).  
QUAST: Quality assessment tool for genome assemblies.  
Bioinformatics, 29(8), 1072–1075.  
https://doi.org/10.1093/bioinformatics/btt086  

Simão, F. A., Waterhouse, R. M., Ioannidis, P., Kriventseva, E. V., & Zdobnov, E. M. (2015).  
BUSCO: Assessing genome assembly and annotation completeness with single-copy orthologs.  
Bioinformatics, 31(19), 3210–3212.  
https://doi.org/10.1093/bioinformatics/btv351  

##  References – Genome Profiling

Marçais, G., & Kingsford, C. (2011).  
A fast, lock-free approach for efficient parallel counting of occurrences of k-mers.  
Bioinformatics, 27(6), 764–770.  
https://doi.org/10.1093/bioinformatics/btr011  

Vurture, G. W., Sedlazeck, F. J., Nattestad, M., Underwood, C. J., Fang, H.,  
Gurtowski, J., & Schatz, M. C. (2017).  
GenomeScope: Fast reference-free genome profiling from short reads.  
Bioinformatics, 33(14), 2202–2204.  
https://doi.org/10.1093/bioinformatics/btx153  

Ranallo-Benavidez, T. R., Jaron, K. S., & Schatz, M. C. (2020).  
GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes.  
Nature Communications, 11, 1432.  
https://doi.org/10.1038/s41467-020-14998-3  


##  References – Genome Annotation

Shumate, A., & Salzberg, S. L. (2021).  
Liftoff: Accurate mapping of gene annotations.  
Bioinformatics, 37(12), 1639–1643.  
https://doi.org/10.1093/bioinformatics/btaa1016  

Flynn, J. M., Hubley, R., Goubert, C., Rosen, J., Clark, A. G., Feschotte, C., & Smit, A. F. A. (2020).  
RepeatModeler2 for automated genomic discovery of transposable element families.  
Proceedings of the National Academy of Sciences, 117(17), 9451–9457.  
https://doi.org/10.1073/pnas.1921046117  

Smit, A. F. A., Hubley, R., & Green, P. (2013–2015).  
RepeatMasker Open-4.0.  
http://www.repeatmasker.org  
## References -Global and local Ancestry Inference

The 1000 Genomes Project Consortium. (2015).
A global reference for human genetic variation. Nature, 526(7571), 68–74.
https://doi.org/10.1038/nature15393

Maples, B. K., Gravel, S., Kenny, E. E., & Bustamante, C. D. (2013).
RFMix: A discriminative modeling approach for rapid and robust local-ancestry inference. American Journal of Human Genetics, 93(2), 278–288.https://doi.org/10.1016/j.ajhg.2013.06.020

Alexander, D. H., Novembre, J., & Lange, K. (2009).
Fast model-based estimation of ancestry in unrelated individuals. Genome Research, 19(9), 1655–1664. https://doi.org/10.1101/gr.094052.109

Shriver, M. D., Smith, M. W., Jin, L., Marcini, A., Akey, J. M., Deka, R., & Ferrell, R. E. (1997).
Ethnic-affiliation estimation by use of population-specific DNA markers. American Journal of Human Genetics, 60(4), 957–964.

Rosenberg, N. A., Pritchard, J. K., Weber, J. L., Cann, H. M., Kidd, K. K., Zhivotovsky, L. A., & Feldman, M. W. (2003).
Genetic structure of human populations. Science, 298(5602), 2381–2385. https://doi.org/10.1126/science.1078311

Kosoy, R., Nassir, R., Tian, C., White, P. A., Butler, L. M., Silva, G., Kittles, R., Alarcon-Riquelme, M. E., Gregersen, P. K., Belmont, J. W., & Seldin, M. F. (2009). Ancestry informative marker sets for determining continental origin and admixture proportions in common populations in America. Human Mutation, 30(1), 69–78. https://doi.org/10.1002/humu.20822

Phillips, C. (2014). Forensic genetic analysis of bio-geographical ancestry. Forensic Science International: Genetics, 12, 49–65.
https://doi.org/10.1016/j.fsigen.2014.05.012

