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



# ===========================================================================
# Part 2: differential expression, design ~ replicate + time
# ===========================================================================

# --- 5. Fit the model -------------------------------------------------------
dds$time <- relevel(dds$time, ref = "t000")
design(dds) <- ~ replicate + time

dds <- DESeq(dds)
resultsNames(dds)


# --- 6. One contrast by hand: t040 vs t000 ----------------------------------
res40 <- results(dds, contrast = c("time", "t040", "t000"), alpha = 0.05)
summary(res40)

res40$gene_name <- gene_info$gene_name[match(rownames(res40), gene_info$Geneid)]
res40 <- res40[order(res40$padj), ]
head(as.data.frame(res40), 10)

# --- 7. Sanity check: known EGF-responsive genes ---------------------------
as.data.frame(res40)[res40$gene_name %in% c("FOS", "EGR1", "MYC", "SRSF3", "GAPDH"),
                     c("gene_name", "baseMean", "log2FoldChange", "padj")]



# --- 8. Time profiles of immediate-early genes ------------------------------
plot_gene <- function(gene_symbol) {
  id <- gene_info$Geneid[gene_info$gene_name == gene_symbol][1]
  d  <- plotCounts(dds, gene = id, intgroup = c("time", "replicate"),
                   returnData = TRUE)
  ggplot(d, aes(time, count, color = replicate, group = replicate)) +
    geom_point(size = 3) + geom_line() +
    scale_y_log10() +
    labs(title = gene_symbol, y = "Normalized counts (log scale)") +
    theme_bw()
}

plot_gene("FOS")
plot_gene("EGR1")


# IEG notes (2026-10-30):
# - FOS and EGR1 are already high at t000 (~10x above t180/t300) in all 3 replicates.
# - Methods: after 18 h starvation, cells were pre-incubated 20 min with DMSO
#   before EGF. t000 = 20 min after handling -> IEG peak from handling, not EGF.
# - t000 is NOT a baseline for IEGs (contrasts vs t000 underestimate or invert them),
#   nor steady state for kinetic modelling of these genes. Likely fine for slower genes.


# --- 9. LRT: genes changing at any time point --------------------------------
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ replicate)
res_lrt <- results(dds_lrt, alpha = 0.05)
summary(res_lrt)

res_lrt$gene_name <- gene_info$gene_name[match(rownames(res_lrt), gene_info$Geneid)]
res_lrt <- res_lrt[order(res_lrt$padj), ]

head(as.data.frame(res_lrt)[, c("gene_name", "baseMean", "padj")], 15)
as.data.frame(res_lrt)[res_lrt$gene_name %in% c("FOS", "EGR1", "MYC", "SRSF3", "GAPDH"),
                       c("gene_name", "baseMean", "padj")]



# --- 10. Effect size: amplitude across the time course ----------------------
norm_log <- log2(counts(dds, normalized = TRUE) + 1)

time_means <- sapply(levels(dds$time), function(t)
  rowMeans(norm_log[, dds$time == t]))

amplitude <- apply(time_means, 1, max) - apply(time_means, 1, min)

res_lrt$amplitude <- amplitude[rownames(res_lrt)]

sig_strong <- res_lrt[which(res_lrt$padj < 0.05 & res_lrt$amplitude >= 1), ]
nrow(sig_strong)

as.data.frame(res_lrt)[res_lrt$gene_name %in% c("FOS", "EGR1", "MYC", "SRSF3", "GAPDH"),
                       c("gene_name", "padj", "amplitude")]



# --- 11. Save results ----------------------------------------------------------
dir.create("results/deseq2", showWarnings = FALSE, recursive = TRUE)

write.csv(as.data.frame(res_lrt), "results/deseq2/lrt_time_results.csv")

time_means_df <- data.frame(gene_name = gene_info$gene_name[match(rownames(time_means), gene_info$Geneid)],
                            time_means, check.names = FALSE)
write.csv(time_means_df, "results/deseq2/time_means_log2.csv")

saveRDS(dds,     "results/deseq2/dds_wald.rds")
saveRDS(dds_lrt, "results/deseq2/dds_lrt.rds")





dds     <- readRDS("results/deseq2/dds_wald.rds")
res_lrt <- read.csv("results/deseq2/lrt_time_results.csv", row.names = 1)
nrow(res_lrt[which(res_lrt$padj < 0.05 & res_lrt$amplitude >= 1), ])
res_lrt[res_lrt$gene_name %in% c("FOS", "EGR1", "MYC", "SRSF3", "GAPDH"), c("gene_name", "padj", "amplitude")]




#-----------------------------------------#

lrt <- read.csv("results/deseq2/lrt_time_results.csv", row.names = 1)
tm  <- read.csv("results/deseq2/time_means_log2.csv",  row.names = 1, check.names = FALSE)

# 1. ¿Cambia y cuánto?
lrt[lrt$gene_name == "SRSF3", c("gene_name", "baseMean", "padj", "amplitude")]

# 2. Su expresión en cada tiempo (log2)
tm[tm$gene_name == "SRSF3", ]


library(DESeq2); library(ggplot2)
dds <- readRDS("results/deseq2/dds_wald.rds")
gene_info <- read.delim("results/counts/gene_counts.txt", comment.char = "#",
                        check.names = FALSE)[, c("Geneid", "gene_name")]

plot_gene <- function(gene_symbol) {
  id <- gene_info$Geneid[gene_info$gene_name == gene_symbol][1]
  d  <- plotCounts(dds, gene = id, intgroup = c("time", "replicate"), returnData = TRUE)
  ggplot(d, aes(time, count, color = replicate, group = replicate)) +
    geom_point(size = 3) + geom_line() + scale_y_log10() +
    labs(title = gene_symbol, y = "Normalized counts (log scale)") + theme_bw()
}
plot_gene("SRSF3")

