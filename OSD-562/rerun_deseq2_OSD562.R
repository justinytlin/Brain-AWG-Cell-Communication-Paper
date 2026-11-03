#!/usr/bin/env Rscript
# rerun_deseq2_OSD562.R
# ---------------------------------------------------------------
# DESeq2 re-analysis of OSD-562 / GLDS-557 (RRRM-2 hippocampus):
#   ISS-Terminal, 20-week (young) spaceflight (FLT, n = 4)
#   vs ISS-Terminal, 20-week (young) ground control (GC, n = 3)
#
# Called by rerun_deseq2_OSD562.py, which first builds the
# integer counts matrix from the RSEM *.genes.tsv expected_count
# columns (ERCC spike-ins and GENCODE vM37 rRNA / Mt_rRNA /
# rRNA_pseudogene genes removed, zero-sum rows removed).
#
# Settings mirror the GeneLab bulk RNA-seq pipeline (GL-DPPD-7101):
#   design = ~ condition  (GC = reference level)
#   DESeq() defaults: median-of-ratios size factors, parametric
#   dispersion fit, Wald test, independent filtering, BH adjustment;
#   NO lfcShrink.  VST exported with blind = TRUE for PCA.
#
# Usage:
#   Rscript rerun_deseq2_OSD562.R <counts.csv> <coldata.csv> <out_prefix>
# ---------------------------------------------------------------
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript rerun_deseq2_OSD562.R <counts.csv> <coldata.csv> <out_prefix>")
}
counts_file  <- args[1]
coldata_file <- args[2]
out_prefix   <- args[3]

cat("DESeq2 version:", as.character(packageVersion("DESeq2")), "\n")

counts  <- as.matrix(read.csv(counts_file,  row.names = 1, check.names = FALSE))
coldata <- read.csv(coldata_file, row.names = 1, check.names = FALSE)
stopifnot(identical(colnames(counts), rownames(coldata)))
storage.mode(counts) <- "integer"

coldata$condition <- factor(coldata$condition, levels = c("GC", "FLT"))
cat("Samples:\n"); print(coldata)
cat("Genes in counts matrix:", nrow(counts), "\n")

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = coldata,
                              design    = ~ condition)
dds <- DESeq(dds)   # defaults: Wald test, parametric fit, BH

res <- results(dds, contrast = c("condition", "FLT", "GC"),
               alpha = 0.05)         # no lfcShrink (GeneLab default)
res_df <- as.data.frame(res)
res_df$ENSEMBL <- rownames(res_df)

norm_counts <- counts(dds, normalized = TRUE)
norm_df <- as.data.frame(norm_counts)
norm_df$ENSEMBL <- rownames(norm_df)

vsd <- vst(dds, blind = TRUE)
vst_df <- as.data.frame(assay(vsd))

write.csv(res_df,  paste0(out_prefix, "_deseq2_results_raw.csv"), row.names = FALSE)
write.csv(norm_df, paste0(out_prefix, "_normalized_counts.csv"),  row.names = FALSE)
write.csv(vst_df,  paste0(out_prefix, "_vst.csv"),                row.names = TRUE)
write.csv(as.data.frame(sizeFactors(dds)), paste0(out_prefix, "_size_factors.csv"))

cat("\nSize factors:\n"); print(sizeFactors(dds))
cat("\nsummary(results):\n"); summary(res)
cat("Done.\n")
