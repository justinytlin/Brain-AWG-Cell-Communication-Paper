# OSD-562 (GLDS-557, RRRM-2 hippocampus) re-analysis: ISS-T FLT (n = 4) vs ISS-T GC (n = 3)

Date: 2026-11-03

## Why the re-run

The GeneLab differential-expression table
`GLDS-557_rna_seq_differential_expression_rRNArm_GLbulkRNAseq.csv`, contrast
`(12 week & Space Flight & On ISS)v(12 week & Ground Control & On Earth)`, is **not**
the comparison described in the manuscript. Because FY10 is mislabelled "On Earth" in the
ISA, that contrast is 3 FLT (FY3, FY4, FY6) vs 6 GC (ISS-T GY5, GY8, GY10 **plus** LAR GY1,
GY4, GY5). All OSD-562 figures were therefore regenerated from a clean 4 vs 3 contrast:

| Group | Samples |
|---|---|
| FLT (spaceflight, ISS-Terminal, 20 wk) | RRRM2_HPC_FLT_ISS-T_YNG_FY3, FY4, FY6, FY10 |
| GC (ground control, ISS-Terminal, 20 wk) | RRRM2_HPC_GC_ISS-T_YNG_GY5, GY8, GY10 |

## Counts source (primary = GeneLab official counts)

**Primary (current outputs):** GeneLab's official RSEM unnormalized counts table
`GeneLab_counts_download/GLDS-557_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv`
(56,924 rows x 27 samples; GeneLab pipeline NF_RCP 2.0.0, Ensembl 112; rRNA already removed by
GeneLab). This is the quantification that matches the GeneLab DE/VST tables used for OSD-564, so
both datasets are now on the same pipeline version and annotation. Run:
`/opt/anaconda3/bin/python3 rerun_deseq2_OSD562.py` (default `--source genelab`).

**Alternative (kept, suffix `_STARlocal`):** the RSEM `*.genes.tsv` files in
`../OSD-562 Hipathia/STAR/`, with GENCODE vM37 rRNA/Mt_rRNA removed locally
(`--source star`). Outputs from that route were renamed
`OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3_STARlocal.csv`, `OSD562_VST_counts_ISS-T_YNG_STARlocal.csv`,
`OSD562_raw_counts_ISS-T_YNG_STARlocal.csv`, `OSD562_DESeq2_sanity_check_STARlocal.csv`, and the
figures/enrichment built from them live in `_STARlocal_4v3_outputs/`.

The two sources are the same quantification: identical per-sample library totals and identical
per-gene counts for all 28,396 shared genes. The only difference is that GeneLab's "rRNArm" keeps
the two mitochondrial rRNA genes mt-Rnr1 / mt-Rnr2 (ENSMUSG00000064337 / ENSMUSG00000064339,
~32 M counts across the 7 samples) that the vM37 Mt_rRNA filter removed, so the GeneLab-based
matrix has 28,398 genes vs 28,396. log2FC Pearson r between the two 4 v 3 runs = 1.00000
(Spearman 1.00000); size factors agree to 4 decimals; p-values agree to 3-4 significant digits.

## What was done (Task 1: differential expression)

Scripts: `rerun_deseq2_OSD562.py` (driver) + `rerun_deseq2_OSD562.R` (DESeq2).

1. **Counts matrix**: 7 samples from the GeneLab counts table, rounded to integer (as
   `DESeqDataSetFromTximport` does). 92 ERCC spike-in rows dropped; 28,434 zero-sum rows dropped.
   Final matrix: **28,398 genes x 7 samples** (`OSD562_raw_counts_ISS-T_YNG.csv`).
   (`--source star` additionally drops GENCODE vM37 rRNA/Mt_rRNA/rRNA_pseudogene genes, cached in
   `OSD562_rRNA_genes_vM37.tsv`: 354 rRNA + 2 Mt_rRNA, no rRNA_pseudogene biotype in vM37.)
2. **DESeq2 1.52.0 (R 4.6.1)**, installed via BiocManager (Bioconductor 3.23); pydeseq2 fallback
   not needed. `design = ~ condition`, GC = reference level, `DESeq()` defaults (median-of-ratios
   size factors, parametric dispersion fit, Wald test, independent filtering, Benjamini-Hochberg),
   **no lfcShrink** (GeneLab default); `vst(dds, blind = TRUE)` exported for PCA.
   Size factors: FY3 0.587, FY4 0.468, FY6 2.529, FY10 1.759, GY5 0.783, GY8 0.556, GY10 2.052.
   Cook's outlier filtering set p/padj to NA for 36 genes (28,362 genes tested).
3. **Outputs**
   * `OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3.csv` - GeneLab-style: ENSEMBL, SYMBOL, GENENAME,
     ENTREZID (annotation copied from the GeneLab DE table, Ensembl 112), 7 DESeq2-normalized
     count columns, `Group.Mean_(Space Flight)`, `Group.Mean_(Ground Control)`, group SDs,
     `Log2fc_(Space Flight)v(Ground Control)`, lfcSE, Stat,
     `P.value_(Space Flight)v(Ground Control)`, `Adj.p.value_(Space Flight)v(Ground Control)`,
     baseMean. Positive log2FC = higher in spaceflight.
   * `OSD562_VST_counts_ISS-T_YNG.csv` - VST matrix, 7 samples.
   * `OSD562_DESeq2_sanity_check.csv` - GeneLab 3 v 6 vs STARlocal 4 v 3 vs GeneLab 4 v 3.

### Results (GeneLab-counts 4 vs 3 run)

* adj p <= 0.05 & |log2FC| >= 1: **0 up, 0 down**; 0 genes at adj p <= 0.05 at any fold change
  (minimum padj = 0.839).
* nominal p <= 0.05 & |log2FC| >= 1: **197 up, 165 down** (362); nominal p <= 0.05 at any
  fold change: 531 genes.
* PCA (top 2,000 most-variable VST genes, 7 samples): PC1 = 58.5 %, PC2 = 18.8 %, PC3 = 9.8 %.
  **FLT and GC do not separate**: PC1 splits {FY3, FY4, GY5, GY8} (negative) from
  {FY6, FY10, GY10} (positive), i.e. the 7SL-high vs 7SL-low libraries (see below).

### Caveats you should know about

* **The libraries are dominated by 7SL/7SK RNA.** In the GeneLab counts, Rn7s1 + Rn7s2 + Rn7sk
  account for: FY3 78.7 %, FY4 78.7 %, FY6 21.6 %, FY10 36.0 %, GY5 71.4 %, GY8 71.2 %,
  GY10 36.3 % of all counts (library sizes 56-83 M). This is why the DESeq2 size factors span
  0.47-2.53 despite similar read depth, and why PC1 tracks library composition rather than flight
  status. GeneLab's "rRNArm" filter removes rRNA biotypes only, so these genes are retained here
  as in GeneLab.
* **GeneLab's published DE table is from a different quantification than its published counts
  table.** The GeneLab unnormalized counts (and the STAR-folder RSEM files, which are identical)
  reproduce the GeneLab-normalized columns only approximately (log-space r = 0.992-0.998; per-gene
  log2 SD ~ 0.2-0.3; ~1 % of expressed genes differ > 2-fold). Hmga1b is the extreme case: the
  GeneLab DE table has ~936 raw-equivalent reads in FY3, the counts table has 7, so GeneLab's
  3 v 6 top hit (log2FC 7.3, padj 0.0036) is not reproducible from the official counts. The
  GeneLab DE table also has 26,066 genes vs 28,398 here.
* With 0 adjusted-significant genes, ORA steps in the enrichment scripts are skipped
  automatically (< 5 genes); all pathway results are GSEA-prerank.

### Sanity check: GeneLab 3 v 6 vs STARlocal 4 v 3 vs GeneLab-counts 4 v 3 (log2FC / p / padj)

| Gene | GeneLab 3 v 6 | STARlocal 4 v 3 | GeneLab-counts 4 v 3 | same direction (3v6 vs new) |
|---|---|---|---|---|
| Cpt1b | 2.325 / 3.0e-03 / 0.311 | 0.478 / 0.519 / 1.000 | 0.478 / 0.519 / 1.000 | yes |
| Pprc1 | 1.013 / 7.3e-03 / 0.369 | -0.054 / 0.909 / 1.000 | -0.054 / 0.909 / 1.000 | no (~0) |
| Sdc4 | -1.076 / 1.5e-02 / 0.440 | 0.019 / 0.965 / 1.000 | 0.019 / 0.965 / 1.000 | no (~0) |
| Rhof | 0.663 / 3.2e-02 / 0.488 | 0.118 / 0.668 / 1.000 | 0.118 / 0.668 / 1.000 | yes |
| Ghsr | 1.578 / 2.5e-02 / 0.478 | 1.215 / 0.109 / 1.000 | 1.215 / 0.109 / 1.000 | yes |
| Hmga1b | 7.296 / 2.1e-07 / 0.0036 | 0.205 / 0.850 / 1.000 | 0.205 / 0.850 / 1.000 | yes (see caveat) |
| Ttn | 2.878 / 2.1e-06 / 0.012 | 1.277 / 0.0025 / 1.000 | 1.277 / 0.0025 / 1.000 | yes |
| Smarcd1 | 1.368 / 2.0e-06 / 0.012 | 0.199 / 0.526 / 1.000 | 0.199 / 0.526 / 1.000 | yes |
| Vps4a | -0.080 / 0.742 / 0.926 | 0.022 / 0.942 / 1.000 | 0.022 / 0.942 / 1.000 | no |
| Ist1 | -0.258 / 0.177 / 0.640 | -0.186 / 0.507 / 1.000 | -0.187 / 0.507 / 1.000 | yes |
| Chmp1a | -0.239 / 0.296 / 0.724 | 0.040 / 0.893 / 1.000 | 0.040 / 0.893 / 1.000 | no |
| Chmp6 | -0.196 / 0.342 / 0.750 | 0.035 / 0.904 / 1.000 | 0.035 / 0.905 / 1.000 | no |
| Vta1 | not in GeneLab table | | | |
| Rac2 | -0.297 / 0.556 / 0.859 | 0.021 / 0.973 / 1.000 | 0.021 / 0.973 / 1.000 | no |
| Baiap2l2 | 0.480 / 0.660 / 0.898 | -0.234 / 0.760 / 1.000 | -0.234 / 0.760 / 1.000 | no |
| Prkaa1 | -0.048 / 0.742 / 0.926 | 0.005 / 0.982 / 1.000 | 0.005 / 0.982 / 1.000 | no |

Pprc1 and Sdc4 flip sign but are essentially at zero in the new runs (|log2FC| < 0.06).

## Task 2: figure scripts

Outputs built from the GeneLab 3 v 6 table were moved, not deleted, to `_old_GeneLab_3v6_outputs/`; outputs built from the STARlocal 4 v 3 run were moved to `_STARlocal_4v3_outputs/`
(`OSD562_Volcano_plot*.png`, `OSD562_PCA_plot.png`, `OSD562_Significant_Genes.xlsx`,
`KEGG_results_OSD562/`, `GO_BP_results_OSD562/`, `QC_figures/`).
`OSD562_QC_metrics.png` has no generating script in this folder and is unchanged, so it was
left in place.

All scripts now read `OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3.csv` with
`COMP_TAG = '(Space Flight)v(Ground Control)'`:

* `make_figures.py` - PCA from the 7-sample VST (subset before variance ranking), volcano,
  annotated volcano, `OSD562_Significant_Genes.xlsx`, `OSD562_PCA_coordinates.csv`.
* `volcano_custom_genes_OSD562.py` - genes-of-interest volcano
  (`OSD562_Volcano_custom_genes.png` + `OSD562_Volcano_custom_genes_table.csv`); labels all
  adjusted-significant genes (none in this run); OSD-562 genes of interest added; pipe-joined
  symbols matched by splitting on "|"; COX2 row labelled "mt-Co2".
* `KEGG_enrichment_OSD562.py` - KEGG_2019_Mouse, GSEA prerank (sign(log2FC) x -log10 p,
  1000 permutations, seed 42, min 10, max 500) -> `KEGG_results_OSD562/`.
* `GO_BP_enrichment_OSD562.py [BP|CC]` - GO_Biological_Process_2026 -> `GO_BP_results_OSD562/`;
  `CC` argument runs GO_Cellular_Component_2026 -> `GO_CC_results_OSD562/`.
* `filter_GSEA_pathways_OSD562.py` - pulls the harmonized panel terms (and the extra report
  terms) into `OSD562_GSEA_panel_terms_summary.csv`,
  `KEGG_results_OSD562/GSEA_KEGG_pathway_filtered.csv`,
  `GO_BP_results_OSD562/GSEA_GO_BP_mito_filtered.csv`.
* `make_pathway_figures_OSD562.py` - KEGG panels (exact-name match): sEV (7), TNT (3),
  metabolic adaptation (3, supplementary) -> `KEGG_results_OSD562/pathway_figures/`.
* `make_GO_mito_NES_figure_OSD562.py` - GO BP mitochondrial panel (7 terms, matched by GO ID)
  -> `GO_BP_results_OSD562/GSEA_Mito_NES_summary.png`.
* `make_qc_figure_OSD562.py` - unchanged (already restricted to the 7 ISS-T_YNG samples);
  re-run to regenerate `QC_figures/`. RIN (7 samples) = 8.46 +/- 0.67.

Panel figures show every term with its FDR label ("FDR = 0.012", "FDR < 0.001"), asterisk at
FDR <= 0.05, legend "Enriched in spaceflight (NES > 0)" red / "Enriched in ground control
(NES < 0)" blue, titles "<panel>, OSD-562: spaceflight vs ground control (n = 4 vs 3)".

Run order: `rerun_deseq2_OSD562.py` -> `make_figures.py`, `volcano_custom_genes_OSD562.py`,
`KEGG_enrichment_OSD562.py`, `GO_BP_enrichment_OSD562.py BP`, `GO_BP_enrichment_OSD562.py CC`
-> `filter_GSEA_pathways_OSD562.py`, `make_pathway_figures_OSD562.py`,
`make_GO_mito_NES_figure_OSD562.py`, `make_qc_figure_OSD562.py`.

## GSEA results (GeneLab-counts 4 v 3; see OSD562_GSEA_panel_terms_summary.csv)

Ranked list: 19,257 unique symbols. KEGG_2019_Mouse: 295 sets tested, 13 at FDR <= 0.05.
GO_Biological_Process_2026: 3,581 sets, 28 at FDR <= 0.05. GO_Cellular_Component_2026:
320 sets, 14 at FDR <= 0.05. (FDR = 0 means no permutation exceeded the observed score.)
The STARlocal run gave the same calls (NES within ~0.03; see `_STARlocal_4v3_outputs/`).

| Panel | Term | NES | FDR |
|---|---|---|---|
| GO BP mito | Oxidative Phosphorylation (GO:0006119) | -2.398 | < 0.001 |
| GO BP mito | Mitochondrial Respiratory Chain Complex Assembly (GO:0033108) | -2.249 | < 0.001 |
| GO BP mito | Mitochondrial ATP Synthesis Coupled Electron Transport (GO:0042775) | -2.511 | < 0.001 |
| GO BP mito | Proton Motive Force-Driven Mitochondrial ATP Synthesis (GO:0042776) | -2.356 | < 0.001 |
| GO BP mito | Mitochondrial Electron Transport, NADH To Ubiquinone (GO:0006120) | -2.285 | < 0.001 |
| GO BP mito | Mitochondrial Respiratory Chain Complex I Assembly (GO:0032981) | -2.012 | 0.0093 |
| GO BP mito | Mitochondrial Respiratory Chain Complex IV Assembly (GO:0033617) | -1.890 | 0.057 |
| KEGG sEV | Sphingolipid signaling pathway | 0.818 | 1.000 |
| KEGG sEV | Sphingolipid metabolism | 1.077 | 1.000 |
| KEGG sEV | Phospholipase D signaling pathway | 0.814 | 1.000 |
| KEGG sEV | Endocytosis | -0.698 | 0.991 |
| KEGG sEV | Synaptic vesicle cycle | 1.050 | 1.000 |
| KEGG sEV | Phosphatidylinositol signaling system | 0.962 | 1.000 |
| KEGG sEV | Lysosome | -1.227 | 0.424 |
| KEGG TNT | Axon guidance | -1.080 | 0.647 |
| KEGG TNT | Regulation of actin cytoskeleton | -0.894 | 0.927 |
| KEGG TNT | Rap1 signaling pathway | 0.805 | 1.000 |
| KEGG metabolic | AMPK signaling pathway | -0.723 | 1.000 |
| KEGG metabolic | Glycolysis / Gluconeogenesis | 0.789 | 1.000 |
| KEGG metabolic | Oxidative phosphorylation | -2.549 | < 0.001 |
| KEGG (report) | Ribosome | -2.490 | < 0.001 |
| KEGG (report) | Parkinson disease | -2.132 | < 0.001 |
| GO BP (report) | Axon Guidance (GO:0007411) | -1.229 | 0.613 |
| GO BP (report) | Axonogenesis (GO:0007409) | -1.372 | 0.498 |
| GO BP (report) | Regulation Of Calcium Ion-Dependent Exocytosis (GO:0017158) | 1.808 | 0.656 |
| GO BP (report) | Exocytosis (GO:0006887) | not in GO_Biological_Process_2026 | - |
| GO BP (report) | Regulated Exocytosis (GO:0045055) | -0.668 | 1.000 |
| GO BP (report) | Actin Filament Organization (GO:0007015) | 0.981 | 1.000 |
| GO BP (report) | Mitochondrial Gene Expression (GO:0140053) | -1.726 | 0.158 |
| GO CC (report) | Recycling Endosome Membrane (GO:0055038) | 0.732 | 1.000 |
| GO CC (report) | Endocytic Vesicle (GO:0030139) | -0.784 | 1.000 |
| GO CC (report) | Exocytic Vesicle (GO:0070382) | 1.444 | 1.000 |
| GO CC (report) | Mitochondrial Inner Membrane (GO:0005743) | -2.145 | < 0.001 |

Note: `make_pathway_figures_OSD562.py` re-runs prerank from the cached GMT only to draw the
running-score plots; its NES values differ from the results table in the 3rd decimal
(permutation null is collection-dependent). The bar panels and all tables use
`GSEA_KEGG_results.csv` as the source of truth.
