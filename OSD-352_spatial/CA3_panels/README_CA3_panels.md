# CA3_panels – regenerated spatial panels (Fig 2c, 3c, 4c)

Everything in this folder was produced by two scripts, run in this order:

1. `extract_spatial_data.R` (R 4.6.1, base R + `Matrix`; **no Seurat / STutility**) –
   reads `../RR3-brain_ST_final_object.rds` once and writes `spots.csv`, `expr.csv`,
   `structure_report.txt`, `small_pieces.rds` (~0.3 min, ~6 GB RAM).
2. `plot_CA3_panels.py` (`/opt/anaconda3/bin/python3`, pandas/numpy/matplotlib/scipy) –
   builds the figures and the statistics tables from those CSVs.

## Data source and extraction route

* Object: `RR3-brain_ST_final_object.rds` = Seurat v4 object (class `Seurat`, 29 770 spots,
  12 sections) with an STutility `Staffli` object in `@tools$Staffli`.
* Route that worked: `readRDS()` with Seurat absent, then slot access by
  `attr(obj, "<slot>")`; dgCMatrix subsetting needs `library(Matrix)`. No fallback
  installation of Seurat was necessary.
* Expression: assay **SCT**, slot **`data`** (= log1p of SCT-corrected counts; verified
  `max|data - log1p(counts)| = 2e-16`). `expr.csv` also carries the SCT `counts` slot
  (`<gene>_counts`), used only to compute the fraction of spots with >= 1 count.
* Coordinates: `obj@tools$Staffli@meta.data` columns `pixel_x`, `pixel_y` (pixel position on
  the STutility H&E image; identical to `original_x/original_y`, i.e. no image alignment
  transform was applied). Array coordinates `x, y` (= `adj_x, adj_y`) are also exported.
  Maps are drawn with `pixel_y` increasing downward (image orientation), which matches
  the orientation of the author's `ST.FeaturePlot` output in
  `../Pink1_Snph_Itpr1_spatial_combined.png`.
* Metadata columns exported to `spots.csv`: `barcode, section.name, seurat_clusters, annot,
  sample_condition, nCount_SCT, nFeature_SCT, nCount_RNA, nFeature_RNA, sample_name_cond,
  orig.ident, x, y, adj_x, adj_y, pixel_x, pixel_y, sample (Staffli section index 1-12),
  original_x, original_y`.

## Sections, animals, cluster

| Figure label | section.name  | Staffli idx | sample_name_cond (animal) | condition | spots | cluster-8 spots | median genes/spot (nFeature_SCT) |
|---|---|---|---|---|---|---|---|
| F1 | Sample_159_C1 | 7  | CF1 | flight | 2654 | 158 | 5471 |
| F2 | Sample_158_A1 | 1  | CF2 | flight | 2455 | 64  | 953  |
| F3 | Sample_158_C1 | 3  | CF7 | flight | 3170 | 91  | 5651 |
| G1 | Sample_159_B1 | 6  | CG7 | ground | 3330 | 105 | 2002 |
| G2 | Sample_304_C1 | 11 | CG8 | ground | 1301 | 87  | 3324 |
| G3 | Sample_304_A1 | 9  | CG9 | ground | 2298 | 216 | 2843 |

(The other six sections – 158_B1 [CF2], 158_D1 [CF7], 159_A1 [CG7], 159_D1 [CF1], 304_B1 [CG9],
304_D1 [CG8] – are the second section of each animal; they are not drawn but their
cluster-8 statistics are in `CA3_cluster8_per_section_stats.csv`.)

* **CA3 = `seurat_clusters == 8`**, `annot` = "8-Synapses and metabolism in hippocampus CA3
  region" (1512 spots over all 12 sections). Cluster 6 = CA1 pyramidal, 10 = CA1 stratum
  radiatum, 11 = dentate gyrus (not used).
* Labels F1–F3 / G1–G3 follow the order of the original animal IDs (CF1 < CF2 < CF7;
  CG7 < CG8 < CG9); the original ID is printed in parentheses on every panel.

## What is plotted

Each figure = one column per gene (gene name in the header), three rows:

1. **Whole section** – all spots of the section as circles (~0.9 x spot pitch), colour =
   SCT `data` value, **rescaled within each section from 0 (white/faint pink) to that section's
   99th percentile of all spots (dark red)**, values above the 99th percentile saturate.
   Cluster-8 spots have a thin black outline. Each map is labelled with the animal ID and the
   section's median number of detected genes per spot. All maps share the same pixel scale
   (1623 x 1623 px window centred on the section).
2. **CA3 only** – identical layout, identical per-section colour scale; only cluster-8 spots
   are coloured (thin black outline), every other spot is light grey. Label gives the number
   of cluster-8 spots in that section.
3. **Per-animal CDR-adjusted CA3 mean** – one point per animal (flight red, ground blue);
   horizontal bar = mean of the three animals. Computed from **all cluster-8 spots of all 12
   sections (both sections of each animal pooled)**: cellular detection rate
   CDR = nFeature_SCT / 14 630 (fraction of SCT genes detected in the spot), z-scored across the
   1512 cluster-8 spots; ordinary least squares `SCT_data ~ 1 + cdr_z`; adjusted expression per
   spot = residual + intercept; then averaged per animal (CF1 283 spots, CF2 140, CF7 189, CG7
   242, CG8 184, CG9 474). The title gives the MAST coefficient and FDR from Masarapu et al.
   (spot-level MAST hurdle model with CDR covariate and sample as random effect). These
   statistics are NOT recomputed here; the adjusted means are descriptive only, and the
   figure footnote says so.
   The unadjusted means are kept in the supplementary `Fig*_CA3_raw_means.png/.pdf`
   (top row: the single section per animal drawn in the maps; bottom row: both sections pooled).

Shared colour bar: "Low (0)" to "High (section 99th pct.)".
Figure size: 60 mm per gene column + 14 mm margin (Fig2c/Fig3c = 194 x 252 mm, Fig4c =
134 x 270 mm); fonts 5.5–9 pt; PNG and PDF at 300 dpi; scatter layers are rasterised in the
PDF, text stays vector (fonttype 42).

## Outputs

* `Fig2c_Pink1_Snph_Itpr1_CA3.png / .pdf`
* `Fig3c_Rhob_Sec14l1_Ubl3_CA3.png / .pdf`
* `Fig4c_Arpc2_Arap2_CA3.png / .pdf`
* `Fig2c_Pink1_Snph_Itpr1_CA3_raw_means.png / .pdf`, `Fig3c_Rhob_Sec14l1_Ubl3_CA3_raw_means.png / .pdf`,
  `Fig4c_Arpc2_Arap2_CA3_raw_means.png / .pdf` – supplementary, unadjusted per-animal means.
* `CA3_cluster8_per_section_stats.csv` – for every gene x all 12 sections: n spots, n cluster-8
  spots, median nFeature_SCT, cluster-8 mean/median SCT `data`, fraction of cluster-8 spots with
  >= 1 SCT count, mean SCT count, section 99th percentile used for the colour scale.
* `CA3_cluster8_per_animal_adjusted.csv` – per gene x animal (both sections pooled): n cluster-8
  spots, raw mean, CDR-adjusted mean, CDR+log-depth-adjusted mean, fraction expressing, mean CDR.
* `CA3_cluster8_per_animal_summary.csv` – per gene: the six raw per-animal cluster-8 means
  (figure section), flight/ground means and difference, sign agreement with MAST, fraction of
  expressing spots, the pooled-both-sections raw comparison, and the CDR-adjusted columns
  (`adj_F1(CF1)` ... `adj_G3(CG9)`, `adj_flight_mean`, `adj_ground_mean`, `adj_diff_F_minus_G`,
  `adj_agrees_with_MAST`, `adj2_cdr_depth_diff`, `adj2_agrees_with_MAST`, `cdr_slope_per_SD`,
  `r_cdr_vs_expr`, `frac_expr_*_pooled`).
* `spots.csv`, `expr.csv`, `structure_report.txt`, `extract_log.txt`, `plot_log.txt`,
  `small_pieces.rds` (meta.data + coordinates + expression subset, so the 1.7 GB object
  need not be reloaded).

## Raw per-animal cluster-8 means (SCT `data`), figure sections (supplementary figure, top row)

| gene | MAST coef / FDR | F1 | F2 | F3 | G1 | G2 | G3 | flight mean | ground mean | F - G | same sign as MAST? |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Pink1   | -0.104 / 0.043 | 1.532 | 0.234 | 1.541 | 0.504 | 1.059 | 0.880 | 1.102 | 0.815 | +0.288 | **no** |
| Snph    | -0.170 / 0.049 | 1.524 | 0.214 | 1.340 | 0.745 | 1.156 | 0.936 | 1.026 | 0.946 | +0.080 | **no** (yes when both sections per animal are pooled: -0.057) |
| Itpr1   | +0.246 / 0.040 | 0.714 | 0.104 | 0.823 | 0.026 | 0.390 | 0.392 | 0.547 | 0.270 | +0.277 | yes |
| Rhob    | +0.488 / 0.027 | 1.522 | 0.428 | 1.731 | 0.616 | 0.506 | 0.480 | 1.227 | 0.534 | +0.693 | yes |
| Sec14l1 | +0.298 / 0.043 | 0.983 | 0.127 | 1.051 | 0.212 | 0.404 | 0.326 | 0.720 | 0.314 | +0.406 | yes |
| Ubl3    | +0.345 / 0.048 | 1.323 | 0.257 | 1.550 | 0.501 | 0.620 | 0.537 | 1.043 | 0.553 | +0.491 | yes |
| Arpc2   | +0.461 / 0.039 | 1.576 | 0.334 | 2.011 | 0.619 | 0.637 | 0.602 | 1.307 | 0.619 | +0.688 | yes |
| Arap2   | +0.190 / 0.048 | 0.392 | 0.043 | 0.586 | 0.079 | 0.172 | 0.132 | 0.340 | 0.128 | +0.212 | yes |

**Why the raw means are not in the main figure:** they track sequencing depth. The two deepest
sections (CF1 159_C1 and CF7 158_C1, ~5500 genes/spot) are both flight, the shallowest (CF2
158_A1, 953 genes/spot) is flight, and the ground sections sit in between (2000–3300
genes/spot). SCT `data` is not depth-independent at this spread (correlation of expression with
CDR across the 1512 cluster-8 spots r = 0.36–0.67 for the 8 genes), so for every gene,
including Pink1 and Snph which MAST calls down in flight, F1 and F3 are the highest raw points
and F2 the lowest. The spatial maps avoid this by per-section scaling; the quantitative row
uses the CDR adjustment below.

## CDR-adjusted per-animal cluster-8 means (main figure, bottom row)

All cluster-8 spots, both sections per animal. CDR = nFeature_SCT / 14 630, z-scored over the
1512 cluster-8 spots; per-animal mean CDR: CF1 0.341, CF2 0.113, CF7 0.374 | CG7 0.261,
CG8 0.213, CG9 0.218. Model: `SCT_data ~ 1 + cdr_z` (OLS); adjusted = residual + intercept.

| gene | MAST coef / FDR | CDR slope (per SD) | F1 (CF1) | F2 (CF2) | F3 (CF7) | G1 (CG7) | G2 (CG8) | G3 (CG9) | flight mean | ground mean | F - G | same sign as MAST? | F - G with log(nCount_SCT) added | frac. expressing F / G |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Pink1   | -0.104 / 0.043 | 0.353 | 1.100 | 0.840 | 1.068 | 1.051 | 1.220 | 1.032 | 1.002 | 1.101 | -0.099 | yes | -0.022 (yes) | 0.77 / 0.84 |
| Snph    | -0.170 / 0.049 | 0.289 | 0.986 | 0.740 | 0.950 | 1.078 | 1.170 | 1.022 | 0.892 | 1.090 | -0.198 | yes | -0.141 (yes) | 0.75 / 0.89 |
| Itpr1   | +0.246 / 0.040 | 0.192 | 0.464 | 0.440 | 0.579 | 0.185 | 0.358 | 0.413 | 0.494 | 0.319 | +0.176 | yes | +0.201 (yes) | 0.51 / 0.31 |
| Rhob    | +0.488 / 0.027 | 0.430 | 1.134 | 1.238 | 1.110 | 1.024 | 0.639 | 0.630 | 1.161 | 0.764 | +0.396 | yes | +0.432 (yes) | 0.88 / 0.63 |
| Sec14l1 | +0.298 / 0.043 | 0.299 | 0.591 | 0.682 | 0.690 | 0.507 | 0.479 | 0.440 | 0.654 | 0.475 | +0.179 | yes | +0.193 (yes) | 0.66 / 0.47 |
| Ubl3    | +0.345 / 0.048 | 0.396 | 1.069 | 0.983 | 1.027 | 0.894 | 0.749 | 0.699 | 1.027 | 0.781 | +0.246 | yes | +0.298 (yes) | 0.80 / 0.67 |
| Arpc2   | +0.461 / 0.039 | 0.472 | 1.118 | 1.264 | 1.350 | 1.038 | 0.735 | 0.793 | 1.244 | 0.855 | +0.389 | yes | +0.440 (yes) | 0.86 / 0.70 |
| Arap2   | +0.190 / 0.048 | 0.150 | 0.297 | 0.314 | 0.412 | 0.195 | 0.212 | 0.223 | 0.341 | 0.210 | +0.131 | yes | +0.146 (yes) | 0.41 / 0.23 |

After CDR adjustment all eight genes agree in direction with the MAST coefficient (the raw
means disagreed for Pink1 and Snph), and the adjusted flight-minus-ground differences are of
the same order as the MAST coefficients. Adding log(nCount_SCT) as a second regressor
(`expression ~ cdr_z + log depth`, centred) changes no conclusion; cdr_z and log depth are
strongly collinear, so the individual slopes of that model are not interpretable (the CDR slope
flips sign for Pink1/Snph while log depth absorbs the effect). The fraction of expressing
cluster-8 spots (>= 1 SCT count; mean of the 3 animals per group, both sections pooled) is
unadjusted and depth-sensitive, but it also goes the same way as MAST for all eight genes
(lower in flight for Pink1 and Snph, higher for the other six).

Suggested legend wording for the bottom row: "Points, per-animal mean of CDR-adjusted SCT
log-expression over CA3 (cluster 8) spots, both sections per animal; bar, group mean.
CDR-adjusted = residuals from expression ~ cellular detection rate fitted across all CA3 spots;
descriptive only. Statistics (coef, FDR) are from the spot-level MAST mixed model of Masarapu
et al."

## Main-text maps in the original layout (`plot_whole_section_maps.py`)
`Fig2c_*_maps.png/.pdf`, `Fig3c_*_maps.png/.pdf`, `Fig4c_*_maps.png/.pdf` reproduce the layout of the original `<gene>_spatial_combined` figures (gene columns; Spaceflight | Ground control; three sections each; one Low/High colour bar) with two corrections: colour is rescaled within each section (0 to the section's 99th percentile of SCT log-normalised expression) and the CA3 boundary (cluster-8 spots adjacent to a non-cluster-8 spot) is outlined in black. Sections are one per animal in the fixed order F1 (CF1, 159_C1), F2 (CF2, 158_A1), F3 (CF7, 158_C1) and G1 (CG7, 159_B1), G2 (CG8, 304_C1), G3 (CG9, 304_A1). Use these for the main text; the three-row `*_CA3.png` panels and `*_raw_means.png` go to the supplement.
