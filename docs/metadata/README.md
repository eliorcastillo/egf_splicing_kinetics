# Metadata

- `GSE69822_SraRunTable.csv`: full SRA run table (75 runs), downloaded from SRA Run Selector.
- `samples.tsv`: selected samples for this project (WT, DMSO, EGF time course 0-300 min, n = 18).

## Library characteristics (from GEO sample records)

- Library prep: Illumina TruSeq RNA Sample Prep (poly(A) selection)
- Strandedness: unstranded (to be confirmed from alignments in step 6)
- Layout: single-end, 101 bp
- Replicates: triplicates per time point; biological vs technical not stated in GEO (to check in the paper)

## Quality control (step 4)

- FastQC/MultiQC: good base quality, no adapter contamination (max 0.57%), no overrepresented sequences. No trimming applied.
- High duplication (48-56%) and per-base content bias in the first bases are expected for poly(A) RNA-seq.
- Samples were sequenced on 4 lanes of 3 flowcells (columns `flowcell` and `lane` in `samples.tsv`). Each lane has a technical signature:
  - D1GFLACXX:4 - N calls at cycle 78 (~5% of reads)
  - D1GFLACXX:5 - N calls at cycles 14 and 62 (~7.7% of reads)
  - D1Y2JACXX:3 - GC content peak shifted to 51% (others 46-48%)
  - C1UV8ACXX:1 - no issues
- Lane is not confounded with time point; to be checked as a batch effect in the PCA (step 8).
