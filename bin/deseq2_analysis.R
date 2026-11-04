#!/usr/bin/env Rscript
# DESeq2 + fgsea for MCF7: BKM120 vs DMSO (reference level = DMSO)
# usage: Rscript bin/deseq2_analysis.R [samplesheet] [salmon_dir] [tx2gene] [outdir]
suppressPackageStartupMessages({
  library(tximport); library(DESeq2); library(readr)
  library(ggplot2); library(fgsea); library(msigdbr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(i, default) if (length(args) >= i) args[i] else default
samplesheet <- get_arg(1, "assets/samplesheet_mcf7.csv")
salmon_dir  <- get_arg(2, "results_mcf7/salmon")
tx2gene_f   <- get_arg(3, "reference/tx2gene.tsv")
outdir      <- get_arg(4, "analysis_mcf7")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- samples ----
meta <- as.data.frame(read_csv(samplesheet, show_col_types = FALSE))
rownames(meta) <- meta$sample
meta$condition <- factor(meta$condition, levels = c("DMSO", "BKM"))
stopifnot(!any(is.na(meta$condition)))
files <- file.path(salmon_dir, meta$sample, "quant.sf")
names(files) <- meta$sample
stopifnot(all(file.exists(files)))

# ---- gene-level counts ----
tx2gene <- read_tsv(tx2gene_f, col_names = c("tx", "gene", "symbol"), show_col_types = FALSE)
gene_symbol <- tx2gene[!duplicated(tx2gene$gene), c("gene", "symbol")]
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("tx", "gene")], dropInfReps = TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = meta[, "condition", drop = FALSE], design = ~ condition)
keep <- rowSums(counts(dds) >= 10) >= 3   # >=10 counts in at least 3 samples (smallest group size)
dds <- dds[keep, ]
cat("genes kept after filtering:", nrow(dds), "\n")

# ---- differential expression ----
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "BKM", "DMSO"), alpha = 0.05)
print(summary(res))

res_df <- as.data.frame(res)
res_df$gene   <- rownames(res_df)
res_df$symbol <- gene_symbol$symbol[match(res_df$gene, gene_symbol$gene)]
res_df <- res_df[order(res_df$padj), ]
write_csv(res_df, file.path(outdir, "deseq2_BKM_vs_DMSO.csv"))
write_csv(data.frame(gene = rownames(counts(dds, normalized = TRUE)), counts(dds, normalized = TRUE)),
          file.path(outdir, "normalized_counts.csv"))

# ---- PCA ----
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct <- round(100 * attr(pca, "percentVar"))
p_pca <- ggplot(pca, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4) + geom_text(vjust = -1, show.legend = FALSE) +
  xlab(paste0("PC1: ", pct[1], "% variance")) + ylab(paste0("PC2: ", pct[2], "% variance")) +
  theme_bw() + ggtitle("MCF7: PCA of variance-stabilized counts")
ggsave(file.path(outdir, "pca.png"), p_pca, width = 6, height = 5, dpi = 150)

# ---- volcano ----
v <- res_df[!is.na(res_df$padj), ]
v$neglogp <- -log10(pmax(v$padj, 1e-300))
v$sig <- ifelse(v$padj < 0.05 & abs(v$log2FoldChange) >= 1, "padj<0.05 & |log2FC|>=1", "n.s.")
lab <- v[v$sig != "n.s." & !is.na(v$symbol), ]
lab <- head(lab[order(lab$padj), ], 15)
p_vol <- ggplot(v, aes(log2FoldChange, neglogp, color = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_text(data = lab, aes(label = symbol), size = 3, check_overlap = TRUE, color = "black") +
  scale_color_manual(values = c("grey60", "firebrick")) +
  xlab("log2 fold change (BKM / DMSO)") + ylab("-log10 adjusted p") +
  theme_bw() + ggtitle("MCF7: BKM120 vs DMSO")
ggsave(file.path(outdir, "volcano.png"), p_vol, width = 7, height = 6, dpi = 150)

# ---- pathway analysis (fgsea on MSigDB Hallmark, ranked by Wald statistic) ----
rank_df <- res_df[!is.na(res_df$stat) & !is.na(res_df$symbol), ]
rank_df <- rank_df[order(-abs(rank_df$stat)), ]
rank_df <- rank_df[!duplicated(rank_df$symbol), ]
ranks <- setNames(rank_df$stat, rank_df$symbol)

hallmark <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                     error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

set.seed(42)
fg <- fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500)
fg <- fg[order(fg$padj), ]
fg$leadingEdge <- vapply(fg$leadingEdge, paste, character(1), collapse = ";")
write_csv(as.data.frame(fg), file.path(outdir, "fgsea_hallmark.csv"))

top <- head(as.data.frame(fg), 15)
top$pathway <- sub("^HALLMARK_", "", top$pathway)
p_fg <- ggplot(top, aes(reorder(pathway, NES), NES, fill = padj < 0.05)) +
  geom_col() + coord_flip() + xlab(NULL) + ylab("Normalized enrichment score") +
  scale_fill_manual(values = c("grey60", "firebrick"), name = "padj < 0.05") +
  theme_bw() + ggtitle("Top Hallmark pathways (BKM vs DMSO)")
ggsave(file.path(outdir, "fgsea_hallmark.png"), p_fg, width = 8, height = 6, dpi = 150)

cat("\nPI3K/AKT/mTOR Hallmark result:\n")
print(as.data.frame(fg)[fg$pathway == "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
                        c("pathway", "NES", "pval", "padj", "size")])

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("\nDone. Outputs in", outdir, "\n")
