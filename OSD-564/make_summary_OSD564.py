"""
KEGG GSEA NES panels for OSD-564 (harmonized term sets shared with OSD-562).

Panels
  sEV        (7): Sphingolipid signaling pathway; Sphingolipid metabolism;
                  Phospholipase D signaling pathway; Endocytosis; Synaptic vesicle
                  cycle; Phosphatidylinositol signaling system; Lysosome
  TNT        (3): Axon guidance; Regulation of actin cytoskeleton; Rap1 signaling pathway
  Metabolic  (3): AMPK signaling pathway; Glycolysis / Gluconeogenesis;
                  Oxidative phosphorylation   (supplementary)

The former "Abeta / Tau" panel (mis-titled OSD-562) was dropped: it belonged to an
abandoned framing of the manuscript.

Also writes KEGG_results/pathway_figures/KEGG_leading_edge_genes.xlsx for the
13 harmonized pathways.
"""
import os
import warnings
import pandas as pd

from osd564_nes_panel import plot_nes_panel

warnings.filterwarnings('ignore')

# ──────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────
GROUPS = {
    'sEV': {
        'pathways': [
            "Sphingolipid signaling pathway",
            "Sphingolipid metabolism",
            "Phospholipase D signaling pathway",
            "Endocytosis",
            "Synaptic vesicle cycle",
            "Phosphatidylinositol signaling system",
            "Lysosome",
        ],
        'panel':    'KEGG GSEA, sEV-related pathways',
        'filename': 'GSEA_sEV_pathways_NES.png',
    },
    'TNT': {
        'pathways': [
            "Axon guidance",
            "Regulation of actin cytoskeleton",
            "Rap1 signaling pathway",
        ],
        'panel':    'KEGG GSEA, TNT-related pathways',
        'filename': 'GSEA_TNT_pathways_NES.png',
    },
    'Metabolic_adaptation': {
        'pathways': [
            "AMPK signaling pathway",
            "Glycolysis / Gluconeogenesis",
            "Oxidative phosphorylation",
        ],
        'panel':    'KEGG GSEA, metabolic adaptation pathways (supplementary)',
        'filename': 'GSEA_metabolic_adaptation_pathways_NES.png',
    },
}

TARGET_PATHWAYS = []
for cfg in GROUPS.values():
    for p in cfg['pathways']:
        if p not in TARGET_PATHWAYS:
            TARGET_PATHWAYS.append(p)

GSEA_CSV = 'KEGG_results/GSEA_KEGG_results.csv'
OUT_DIR  = 'KEGG_results/pathway_figures'
os.makedirs(OUT_DIR, exist_ok=True)

# ──────────────────────────────────────────────────────────────────
# 1. Load GSEA results (exact KEGG name match)
# ──────────────────────────────────────────────────────────────────
print("Loading GSEA results ...")
gsea_df = pd.read_csv(GSEA_CSV, index_col=0)
gsea_df['Term'] = gsea_df['Term'].str.strip()
gsea_df = gsea_df.set_index('Term')

records = []
for term in TARGET_PATHWAYS:
    if term not in gsea_df.index:
        print(f"  [MISSING] {term}")
        continue
    row = gsea_df.loc[term]
    genes = [g.strip() for g in str(row['Lead_genes']).split(';')
             if g.strip()] if pd.notna(row['Lead_genes']) else []
    records.append({
        'label':    term,
        'NES':      float(row['NES']),
        'FDR':      float(row['FDR q-val']),
        'NOM_pval': float(row['NOM p-val']),
        'Tag_pct':  row['Tag %'],
        'n_lead':   len(genes),
        'genes':    genes,
    })
    print(f"  {term:<42s} NES={row['NES']:+.3f}  FDR={row['FDR q-val']:.4f}")

plot_df = (pd.DataFrame(records)
           .sort_values('NES', ascending=True)
           .reset_index(drop=True))
print(f"  {len(plot_df)} pathways loaded.")

# ──────────────────────────────────────────────────────────────────
# 2. Grouped NES panels
# ──────────────────────────────────────────────────────────────────
print("\nCreating grouped NES figures ...")
for group_key, cfg in GROUPS.items():
    sub = plot_df[plot_df['label'].isin(cfg['pathways'])]
    recs = sub[['label', 'NES', 'FDR']].to_dict('records')
    out_path = f"{OUT_DIR}/{cfg['filename']}"
    plot_nes_panel(recs, cfg['panel'], out_path)

# ──────────────────────────────────────────────────────────────────
# 3. Leading edge genes — Excel file
# ──────────────────────────────────────────────────────────────────
print("\nCreating leading edge genes Excel file ...")

xl_path = f'{OUT_DIR}/KEGG_leading_edge_genes.xlsx'
with pd.ExcelWriter(xl_path, engine='openpyxl') as writer:

    # ── Sheet 1: Summary (one row per pathway) ──────────────────
    summary = plot_df[['label', 'NES', 'FDR', 'NOM_pval', 'Tag_pct', 'n_lead']].copy()
    summary.columns = ['Pathway', 'NES', 'FDR q-val', 'NOM p-val', 'Tag %', 'N_leading_genes']
    panel_of = {p: k for k, cfg in GROUPS.items() for p in cfg['pathways']}
    summary.insert(1, 'Panel', summary['Pathway'].map(panel_of))
    summary['Leading_genes'] = plot_df['genes'].apply(lambda g: '; '.join(g))
    summary.to_excel(writer, sheet_name='Summary', index=False)

    # ── Sheet 2: Long format (one row per gene) ─────────────────
    long_rows = []
    for _, row in plot_df.iterrows():
        for gene in row['genes']:
            long_rows.append({
                'Pathway':   row['label'],
                'Panel':     panel_of.get(row['label'], ''),
                'NES':       round(row['NES'], 4),
                'FDR q-val': round(row['FDR'], 4),
                'Gene':      gene,
            })
    pd.DataFrame(long_rows).to_excel(writer, sheet_name='Genes_long_format', index=False)

    # ── Sheets 3+: One sheet per pathway ────────────────────────
    for _, row in plot_df.sort_values('label').iterrows():
        sheet_name = row['label'].replace('/', '-')[:31]   # Excel max 31 chars, no '/'
        gene_df = pd.DataFrame({
            'Gene':    row['genes'],
            'Pathway': row['label'],
            'NES':     round(row['NES'], 4),
            'FDR':     round(row['FDR'], 4),
        })
        gene_df.to_excel(writer, sheet_name=sheet_name, index=False)

print(f"  Saved -> {xl_path}")
print("\nAll done.")
