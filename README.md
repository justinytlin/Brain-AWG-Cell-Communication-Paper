# Brain-AWG-Cell-Communication-Paper

Analysis code for *Short-Duration Spaceflight Alters Cellular Transport Pathways in the Hippocampus* (Lin et al.), which asks whether spaceflight alters mitochondrial, tunneling-nanotube (TNT) and small-extracellular-vesicle (sEV) transcriptional programs in the mouse hippocampus using three NASA Open Science Data Repository datasets.

| Folder | Dataset | What is here |
|---|---|---|
| `OSD-564/` | GLDS-569, Rodent Research-10, right hippocampus, FLT n=5 vs GC n=5 | volcano/PCA, GSEApy prerank (KEGG_2019_Mouse, GO 2026), g:Profiler ORA figure scripts, QC figure; GSEA result tables and g:Profiler exports |
| `OSD-562/` | GLDS-557, RRRM-2, right hippocampus, ISS-terminal 20-week cohort, FLT n=4 vs GC n=3 | DESeq2 re-derivation of the 4-vs-3 contrast (`rerun_deseq2_OSD562.py` + `.R`, see `README_OSD562_rerun.md` for why), the resulting DE table, GSEA scripts and results, figure scripts |
| `OSD-352_bulk/` | GLDS-352, Rodent Research-3, pooled whole-brain sections, FLT n=3 vs GC n=3 (supplementary) | volcano/PCA, GSEA scripts and results, RSeQC QC figure |
| `OSD-352_spatial/` | Masarapu et al. 2024 Visium spatial transcriptomics of RR-3 brain | `CA3_panels/extract_spatial_data.R` (base-R extraction of expression, cluster and coordinates from the published Seurat object), `CA3_panels/plot_CA3_panels.py` (per-section-scaled maps, CA3-only maps, CDR-adjusted per-animal means), exported spot tables, `README_CA3_panels.md` |

## Inputs (not stored here)
- GeneLab processed differential-expression, VST and QC tables and ISA metadata: https://osdr.nasa.gov/bio/repo/data/studies/OSD-352, OSD-562, OSD-564
- OSD-562 counts: `GLDS-557_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv` from OSDR (place in `OSD-562/GeneLab_counts_download/`); `rerun_deseq2_OSD562.py --source star` alternatively rebuilds the matrix from the per-sample RSEM files plus the GENCODE vM37 GTF, giving identical counts
- Spatial object `RR3-brain_ST_final_object.rds`: Mendeley Data doi 10.17632/fjxrcbh672.1 (Masarapu et al. 2024, Nat Commun 15:4778); DE statistics from its Supplementary Data 8 (`41467_2024_48916_MOESM11_ESM.xlsx`, CC BY 4.0)
- Per-sample QC: `GLDS-569`, `GLDS-557` and `GLDS-352_rna_seq_qc_metrics_GLbulkRNAseq.csv` from OSDR

## Environment
- Python 3 with pandas, numpy, matplotlib, scikit-learn, gseapy 1.3.0, openpyxl, adjustText
- R 4.6 with DESeq2 1.52 (OSD-562 re-run) and Matrix (spatial extraction); Seurat is not required

## Conventions
- log2FC is spaceflight relative to ground control; positive = up in flight.
- Adjusted-significant: adj p <= 0.05 and |log2FC| >= 1; nominal: p <= 0.05 and |log2FC| >= 1.
- GSEA: gseapy prerank on sign(log2FC) x -log10(p), 1000 permutations, seed 42, gene-set size 10-500. NES > 0 = enriched in spaceflight.
- Identical curated panels are drawn for every dataset: 7 GO BP mitochondrial terms, 7 KEGG sEV pathways, 3 KEGG TNT pathways.

## Running
Each script is run from inside its folder with the input files placed alongside (paths are relative). See `OSD-562/README_OSD562_rerun.md` and `OSD-352_spatial/CA3_panels/README_CA3_panels.md` for the two multi-step pipelines.

## Journal-ready figures
`srep_style.py` (shared typeface/size/line-width rules for Scientific Reports), `srep_nes_panel.py` (all NES bar panels) and `srep_panels.py` (volcano, PCA and QC panels) are copied into each dataset folder so scripts run from their own directory. `assemble_figures.py` places the per-panel PNGs 1:1 onto 183 mm-wide canvases with bold lower-case panel letters and writes `Figures_SciRep/<figure>.png/.tif/.pdf` at 300 dpi. The assembled PNGs are committed in `Figures_SciRep/`.

## Figures (`Figures_SciRep/`)
Numbered as cited in the manuscript.

| File | Content |
|---|---|
| `Figure_1` | OSD-564 volcano (curated genes) and PCA |
| `Figure_2` | Mitochondria: OSD-564 GO BP GSEA (a), OSD-564 g:Profiler ORA, down (b), OSD-562 GO BP GSEA (c), OSD-352 CA3 maps of Pink1, Snph, Itpr1 (d) |
| `Figure_3` | sEVs: OSD-564 KEGG GSEA (a), ORA (b), OSD-352 CA3 maps of Rhob, Sec14l1, Ubl3 (c) |
| `Figure_4` | TNTs: OSD-564 KEGG GSEA (a), ORA (b), OSD-352 CA3 maps of Arpc2, Arap2 (c) |
| `Supplementary_Figure_1` | Per-sample sequencing QC: OSD-564 (a), OSD-562 (b), OSD-352 read distribution (c) |
| `Supplementary_Figure_2` | OSD-562 volcano (a) and PCA (b); OSD-352 whole-brain volcano (c) and PCA (d) |
| `Supplementary_Figure_3` | OSD-352 whole-brain GSEA: KEGG TNT (a), KEGG sEV (b), GO BP mitochondrial (c) |
| `Supplementary_Figure_4` | Additional OSD-564 GSEA: KEGG metabolic adaptation (a), GO BP protrusion/exocytosis (b), GO CC compartments (c) |
| `Supplementary_Figure_5` | OSD-562 KEGG GSEA: TNT (a), sEV (b), metabolic adaptation (c) |
| `Supplementary_Figure_6`-`8` | OSD-352 CA3: whole-section and CA3-only maps with per-animal CDR-adjusted means, for Pink1/Snph/Itpr1 (6), Rhob/Sec14l1/Ubl3 (7), Arpc2/Arap2 (8) |
| `Supplementary_Figure_9`-`11` | Unadjusted per-animal CA3 means for the same three gene groups |

## Tables (`Tables/`)
Built by `build_supplementary_tables.py` (run from the project root, like `assemble_figures.py`). Every workbook has a README sheet describing sources and methods.

| File | Content |
|---|---|
| `Main_Tables_1-6.xlsx` | Main-text Tables 1-6 with legends, exported from the manuscript |
| `Supplementary_Table_S1_OSD-562_DE.xlsx` | OSD-562 re-derived DESeq2 (4 vs 3): 45 curated genes and all 28,398 genes |
| `Supplementary_Table_S2_OSD-352_bulk_DE.xlsx` | OSD-352 whole-brain GeneLab DESeq2 (3 vs 3): 24 curated genes and all genes |
| `Supplementary_Table_S3_GSEA.xlsx` | Complete GSEApy prerank results, one sheet per dataset and library, with the figure/table each displayed term appears in |
| `Supplementary_Table_S4` | **Not yet included**: g:Profiler ORA for OSD-564, pending re-run on the adjusted-p < 0.05 gene lists |
| `Supplementary_Table_S5_OSD-352_spatial.xlsx` | 25 hippocampal spatial DEGs from Masarapu et al. 2024 Supplementary Data 8, and per-animal CA3 means |
| `Supplementary_Table_S6_curated_lists_DRAFT.xlsx` | Curated genes and pathways with categories and placement; rationale and references still to be written |
| `Supplementary_Table_S7_sample_QC.xlsx` | Per-sample RIN, read depth, mapping, rRNA and genes detected; OSD-562 7SL/7SK fraction |
