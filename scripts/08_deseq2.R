# Step 8 - Differential expression with DESeq2
# Part 1: load counts and sample sheet, QC with PCA.
# Run from the project root (RStudio project).

library(DESeq2)
library(ggplot2)

# --- 1. Count matrix (featureCounts output) -----------------------------
counts_raw <- read.delim("results/counts/gene_counts.txt",
                         comment.char = "#", check.names = FALSE)

gene_info <- counts_raw[, c("Geneid", "gene_name", "gene_type", "Length")]

counts <- as.matrix(counts_raw[, 9:ncol(counts_raw)])
rownames(counts) <- counts_raw$Geneid
colnames(counts) <- sub("\\.bam$", "", basename(colnames(counts)))

# --- 2. Sample sheet -------------------------------------------------------
samples <- read.delim("docs/metadata/samples.tsv")
rownames(samples) <- samples$run
samples$time      <- factor(sprintf("t%03d", samples$time_min))
samples$replicate <- factor(samples$replicate)
samples$lane_id   <- factor(paste(samples$flowcell, samples$lane, sep = ":"))

# Counts and sample sheet must be in the same order
stopifnot(identical(colnames(counts), samples$run))

# --- 3. DESeq2 object and filtering ----------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = samples,
                              design    = ~ time)

keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
dim(dds)

# --- 4. PCA ---------------------------------------------------------------
vsd <- vst(dds, blind = TRUE)

pca <- plotPCA(vsd, intgroup = c("time", "replicate"), returnData = TRUE)
pca$lane_id <- samples[pca$name, "lane_id"]
pct <- round(100 * attr(pca, "percentVar"))

ggplot(pca, aes(PC1, PC2, color = time, shape = replicate)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1 (", pct[1], "%)"),
       y = paste0("PC2 (", pct[2], "%)")) +
  theme_bw()







# PCA notes (2026-10-29):
# - PC1 (64%) + PC2 (19%) trace a time arc: t000/t015 -> t040 -> t090 -> t180 -> t300.
# - Replicate 2 is shifted left (earlier) at every time point, across two lanes:
#   replicate effect, not lane. Use design ~ replicate + time.
# - SRR2061418 (lower assignment rate) is not an outlier.