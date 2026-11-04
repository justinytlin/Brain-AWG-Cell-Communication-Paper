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
- Spatial object `RR3-brain_ST_final_object.rds`: Mendeley Data doi 10.17632/fjxrcbh672.1 (Masarapu et al. 2024, Nat Commun 15:4778); DE statistics from its Supplementary Data 8

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
