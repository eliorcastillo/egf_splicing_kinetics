# 08b_deseq2_checks.R
# Exploratory checks on saved DESeq2 results (run 08_deseq2.R first).
# Reads outputs from results/deseq2/; does not refit any model.

library(DESeq2)
library(ggplot2)

# --- 1. Load saved results ---------------------------------------------------
dds        <- readRDS("results/deseq2/dds_wald.rds")
lrt        <- read.csv("results/deseq2/lrt_time_results.csv", row.names = 1)
time_means <- read.csv("results/deseq2/time_means_log2.csv", row.names = 1,
                       check.names = FALSE)
gene_info  <- read.delim("results/counts/gene_counts.txt", comment.char = "#",
                         check.names = FALSE)[, c("Geneid", "gene_name")]

# --- 2. Genes changing over time ---------------------------------------------
# Significant (padj < 0.05) and with amplitude >= 1 log2 unit
nrow(lrt[which(lrt$padj < 0.05 & lrt$amplitude >= 1), ])

# --- 3. Reference genes ------------------------------------------------------
# Immediate-early genes (FOS, EGR1), MYC, the splicing factor SRSF3,
# and GAPDH as a housekeeping control
ref_genes <- c("FOS", "EGR1", "MYC", "SRSF3", "GAPDH")
lrt[lrt$gene_name %in% ref_genes, c("gene_name", "baseMean", "padj", "amplitude")]

# --- 4. Gene of interest: SRSF3 ----------------------------------------------
# Does it change, and by how much?
lrt[lrt$gene_name == "SRSF3", c("gene_name", "baseMean", "padj", "amplitude")]

# Mean expression per time point (log2)
time_means[time_means$gene_name == "SRSF3", ]

# --- 5. Plot normalized counts over time -------------------------------------
plot_gene <- function(gene_symbol) {
  id <- gene_info$Geneid[gene_info$gene_name == gene_symbol][1]
  d  <- plotCounts(dds, gene = id, intgroup = c("time", "replicate"),
                   returnData = TRUE)
  ggplot(d, aes(time, count, color = replicate, group = replicate)) +
    geom_point(size = 3) + geom_line() + scale_y_log10() +
    labs(title = gene_symbol, y = "Normalized counts (log scale)") +
    theme_bw()
}

plot_gene("SRSF3")
