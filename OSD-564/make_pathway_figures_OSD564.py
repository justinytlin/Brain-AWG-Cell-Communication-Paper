"""
KEGG GSEA enrichment plots + combined NES summary for the harmonized OSD-564 panels.

Harmonized KEGG sets (shared with OSD-562; exact KEGG names):
  sEV       : Sphingolipid signaling pathway; Sphingolipid metabolism; Phospholipase D
              signaling pathway; Endocytosis; Synaptic vesicle cycle;
              Phosphatidylinositol signaling system; Lysosome
  TNT       : Axon guidance; Regulation of actin cytoskeleton; Rap1 signaling pathway
  Metabolic : AMPK signaling pathway; Glycolysis / Gluconeogenesis; Oxidative phosphorylation

The former "Abeta / Tau" framing (Alzheimer disease plot mis-titled OSD-562) was
dropped from the panels.  "Alzheimer disease" is still drawn as a plain enrichment
plot so that the old mis-titled file is overwritten with a correctly titled one.

GSEA prerank is re-run (seed 42) only to obtain the running-ES profile for the
enrichment plots; NES / FDR / p shown in every figure are taken from the saved
KEGG_results/GSEA_KEGG_results.csv so the figures match the tables.
"""
import os
import warnings
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gseapy as gp
from gseapy import gseaplot

from osd564_nes_panel import plot_nes_panel, STUDY_TAG

warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'Arial',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300,
})

# ──────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────
PANELS = {
    'sEV': [
        "Sphingolipid signaling pathway",
        "Sphingolipid metabolism",
        "Phospholipase D signaling pathway",
        "Endocytosis",
        "Synaptic vesicle cycle",
        "Phosphatidylinositol signaling system",
        "Lysosome",
    ],
    'TNT': [
        "Axon guidance",
        "Regulation of actin cytoskeleton",
        "Rap1 signaling pathway",
    ],
    'Metabolic adaptation': [
        "AMPK signaling pathway",
        "Glycolysis / Gluconeogenesis",
        "Oxidative phosphorylation",
    ],
}
TARGET_PATHWAYS = []
for lst in PANELS.values():
    for p in lst:
        if p not in TARGET_PATHWAYS:
            TARGET_PATHWAYS.append(p)

# Extra enrichment plots that are not part of any NES panel
EXTRA_ENRICHMENT_PLOTS = ["Alzheimer disease"]

OUT_DIR    = 'KEGG_results/pathway_figures'
GMT_FILE   = 'KEGG_results/GSEA_prerank/gene_sets.gmt'
RANKED_CSV = 'KEGG_results/GSEA_ranked_gene_list.csv'
RESULT_CSV = 'KEGG_results/GSEA_KEGG_results.csv'
PERM_NUM   = 1000

os.makedirs(OUT_DIR, exist_ok=True)

# ──────────────────────────────────────────────────────────────────
# 1. Saved GSEA statistics (source of truth for NES / FDR)
# ──────────────────────────────────────────────────────────────────
print("Loading saved GSEA KEGG results ...")
res = pd.read_csv(RESULT_CSV, index_col=0)
res['Term'] = res['Term'].str.strip()
res = res.set_index('Term')
print(f"  {len(res)} pathways in {RESULT_CSV}")

# ──────────────────────────────────────────────────────────────────
# 2. Load ranked gene list and re-run prerank (for the ES profiles)
# ──────────────────────────────────────────────────────────────────
print("Loading ranked gene list ...")
rnk = pd.read_csv(RANKED_CSV, index_col=0).squeeze('columns')
rnk.name = 'rank_metric'
rnk = rnk.sort_values(ascending=False)
rnk = rnk[~rnk.index.duplicated(keep='first')]
print(f"  {len(rnk):,} genes  |  range [{rnk.min():.2f}, {rnk.max():.2f}]")

print(f"\nRunning GSEA prerank (permutations = {PERM_NUM}) ...")
gsea_res = gp.prerank(
    rnk=rnk,
    gene_sets=GMT_FILE,
    outdir=None,
    seed=42,
    permutation_num=PERM_NUM,
    min_size=10,
    max_size=500,
    threads=4,
    verbose=False,
)
print(f"  Done. {len(gsea_res.results)} pathways computed.")

# ──────────────────────────────────────────────────────────────────
# 3. Enrichment plots (exact KEGG name match)
# ──────────────────────────────────────────────────────────────────
print("\nGenerating enrichment plots ...")
found, missing = 0, 0

for term in TARGET_PATHWAYS + EXTRA_ENRICHMENT_PLOTS:
    if term not in gsea_res.results or term not in res.index:
        print(f"  [MISSING]  '{term}'")
        missing += 1
        continue

    safe_name = (term.replace('/', '-').replace(' ', '_')
                 .replace(':', '').replace(',', '')
                 .replace('(', '').replace(')', ''))
    ofname = f"{OUT_DIR}/GSEA_{safe_name}.png"

    saved = res.loc[term]
    plot_kwargs = dict(gsea_res.results[term])
    # keep the statistics identical to the saved table
    plot_kwargs['nes']  = float(saved['NES'])
    plot_kwargs['pval'] = float(saved['NOM p-val'])
    plot_kwargs['fdr']  = float(saved['FDR q-val'])

    try:
        gseaplot(
            rank_metric=gsea_res.ranking,
            term=f'KEGG: {term}\n{STUDY_TAG}',
            ofname=ofname,
            **plot_kwargs,
        )
        plt.close('all')
        print(f"  [OK]  {term}  ->  {ofname}")
        found += 1
    except Exception as e:
        print(f"  [ERROR]  '{term}': {e}")
        missing += 1

print(f"\nFinished.  {found} plots saved to:  {OUT_DIR}/")
if missing:
    print(f"           {missing} pathway(s) not found or failed.")

# ──────────────────────────────────────────────────────────────────
# 4. Combined summary figure — NES bar plot for the 13 harmonized pathways
# ──────────────────────────────────────────────────────────────────
print("\nCreating combined pathway summary figure ...")
records = []
for panel, lst in PANELS.items():
    for term in lst:
        if term not in res.index:
            print(f"  [MISSING] {term}")
            continue
        row = res.loc[term]
        records.append({
            'label': f'{term}  [{panel}]',
            'NES':   float(row['NES']),
            'FDR':   float(row['FDR q-val']),
        })

summary_path = f'{OUT_DIR}/GSEA_selected_pathways_NES_summary.png'
plot_nes_panel(records, 'KEGG GSEA, harmonized sEV / TNT / metabolic pathways',
               summary_path, fig_w=11.0)
print("Done.")
