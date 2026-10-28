# EGF splicing kinetics

Are EGF-induced changes in PSI genuine splicing decisions, or kinetic effects of gene induction?

## Data

Kiselev et al. 2015, *Nucleic Acids Research* (GEO: GSE69822).
MCF10A cells stimulated with EGF, bulk RNA-seq time course.

## Repository structure

| Folder     | Content                              | Tracked by git |
|------------|--------------------------------------|----------------|
| `envs/`    | conda environment definitions        | yes            |
| `scripts/` | commands for each pipeline step      | yes            |
| `docs/`    | notes and workflow                   | yes            |
| `data/`    | raw FASTQ, reference genome, GTF     | no             |
| `results/` | BAM files, counts, tables            | partially      |

## Pipeline

0. Repository structure
1. Environment setup (conda)
2. Sample selection (GEO)
3. FASTQ download
4. Quality control (FastQC + MultiQC)
5. Reference genome, HISAT2 index and GTF annotation
6. Alignment (HISAT2 → BAM)
7. Gene-level counts (featureCounts)
8. Differential expression (DESeq2)
9. PSI and ΔPSI (rMATS)
10. PSI × expression integration (phase 1)
11. Exon-intron split analysis, EISA (phase 2)
