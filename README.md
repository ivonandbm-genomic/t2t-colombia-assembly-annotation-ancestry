# High-quality diploid chromosome-scale assemblies from admixed Colombian individuals reveal non-uniform structural variation and local ancestry patterns

![Platform](https://img.shields.io/badge/Platform-PacBio%20HiFi-blue)
![Reference](https://img.shields.io/badge/Reference-T2T--CHM13%20v2.0-green)
![Reference](https://img.shields.io/badge/Reference-GRCh38-green)
![License](https://img.shields.io/badge/License-MIT-yellow)

## Official repository of the study

**High-quality diploid chromosome-scale assemblies from admixed Colombian individuals reveal non-uniform structural variation and local ancestry patterns**

---

## Repository Information

| Field       | Description                           |
| ----------- | ------------------------------------- |
| Author      | Ivon Bolaños                          |
| Institution | Universidad del Valle                 |
| Program     | PhD in Biomedical Sciences            |
| Location    | Cali, Colombia                        |
| ORCID       | https://orcid.org/0000-0002-8007-8929 |
| GitHub      | ivonandbm-genomic                     |

---

## Overview

This repository contains the complete bioinformatics workflows used for the assembly, annotation, variant discovery, repeat characterization, methylation profiling, and ancestry inference of two admixed Colombian genomes generated using PacBio HiFi long-read sequencing technology.

The study presents high-quality diploid chromosome-scale assemblies from two individuals from Valle del Cauca, Colombia, providing a genomic resource for investigating structural variation, ancestry composition, repetitive elements, and epigenomic patterns in an underrepresented Latin American population.

---

## Study Design

To simplify figure labeling and downstream analyses, the following aliases were used throughout the repository and manuscript:

| Original Sample ID | Alias |
| ------------------ | ----- |
| 017C               | COLM  |
| 018C               | COLF  |

These aliases are consistently used in assembly analyses, variant discovery, methylation analyses, ancestry inference, figures, and supplementary materials.

### Sequencing Information

| Feature           | Description               |
| ----------------- | ------------------------- |
| Technology        | PacBio HiFi (CCS reads)   |
| Coverage          | ~30×                      |
| Individuals       | 2                         |
| Population        | Admixed Colombian         |
| Geographic Region | Valle del Cauca, Colombia |
| Reference Genomes | GRCh38 and T2T-CHM13 v2.0 |

---

## Scientific Scope

The repository supports analyses related to:

* Genome profiling
* De novo genome assembly
* Chromosome-scale scaffolding
* Assembly quality assessment
* Gene annotation
* Repeat annotation
* Small variant discovery
* Structural variant discovery
* CpG methylation analysis
* Global ancestry inference
* Local ancestry inference
* Comparative analyses against GRCh38 and T2T-CHM13

---

## Workflow Summary

```text
PacBio HiFi Reads
        │
        ▼
Genome Profiling
(Jellyfish + GenomeScope)
        │
        ▼
De Novo Assembly
(Hifiasm)
        │
        ▼
Reference-guided Scaffolding
(RagTag)
        │
        ▼
Assembly Quality Assessment
(QUAST + BUSCO + Merqury)
        │
        ▼
Genome Annotation
(Liftoff)
        │
        ├──────────────┬──────────────┐
        ▼              ▼              ▼
Repeat         Small Variants    Structural Variants
Annotation   DeepVariant/Dipcall  Multi-caller SV
                                      Discovery
        │              │              │
        └──────────────┴──────────────┘
                       │
                       ▼
                Variant Processing
                       │
                       ▼
                Population Structure
               PCA + ADMIXTURE
                       │
                       ▼
                Local Ancestry
                    RFMix
                       
                       
                
```

---

## Repository Structure

```text
t2t-colombia-assembly-annotation-ancestry/

├── README.md
│
├── scripts/
│   ├── 01_assembly_annotation.sh
│   ├── 02_Ancestry.sh
│   └── README.md
│
├── figures/
│
└── supplementary/
```

---

## Main Software

### Genome Assembly

* Hifiasm
* RagTag
* Gfatools

### Assembly Quality Assessment

* QUAST
* BUSCO
* Merqury

### Repeat and Gene Annotation

* RepeatMasker
* RepeatModeler2
* Liftoff

### Small Variant Discovery

* DeepVariant
* Dipcall
* bcftools

### Structural Variant Discovery

* pbsv
* Sniffles2
* NanoVar
* cuteSV
* SVIM-asm
* Sawfish
* SURVIVOR

### Population Genomics

* PLINK
* ADMIXTURE
* WhatsHap
* RFMix

### Genome Profiling

* Jellyfish
* GenomeScope 2.0

---

## Data Availability

Due to ethical, privacy, and consent restrictions associated with human whole-genome sequencing data, raw genomic datasets are not distributed through this GitHub repository.

### Included

* Bioinformatics workflows
* Analysis scripts
* Configuration files
* Documentation
* Summary results
* Figures and supplementary resources

### Not Included

* FASTQ files
* BAM/CRAM files
* Individual-level VCF files
* Individual-level methylation files
* Personally identifiable information
* Clinical metadata

---

## Public Data Deposition

### NCBI BioProject

PRJNA1469122

### Genome Assemblies

| Sample | Assembly Accession |
| ------ | ------------------ |
| COLM   | JBYVEU000000000    |
| COLF   | Pending release    |

---

## Citation

If you use this repository, please cite:

> Bolaños I., Tobar-Tosse F., et al. High-quality diploid chromosome-scale assemblies from admixed Colombian individuals reveal non-uniform structural variation and local ancestry patterns. Submitted manuscript.

---

## License

This project is distributed under the MIT License.

---

## Contact

**Ivon Bolaños**

PhD Candidate in Biomedical Sciences
Universidad del Valle
Cali, Colombia

ORCID: https://orcid.org/0000-0002-8007-8929

