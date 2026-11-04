"""
make_pathway_figures_OSD562.py
------------------------------
KEGG GSEA panel figures for OSD-562 (re-run DESeq2, ISS-T FLT n = 4 vs
ISS-T GC n = 3), using the HARMONIZED panel sets shared with OSD-564 / OSD-352
(identical term order in every dataset):

  * KEGG sEV (7), KEGG TNT (3), KEGG metabolic adaptation (3, supplementary)

Terms are matched by EXACT KEGG name (case-insensitive) in
KEGG_results_OSD562/GSEA_KEGG_results.csv (the source of truth for NES/FDR);
every term is drawn with its FDR label whether or not significant (asterisk
at FDR <= 0.05, pale bar above 0.05).  The panels are drawn by the shared
srep_nes_panel.nes_panel (Scientific Reports style, 90 x 60 mm, 300 dpi PNG +
vector PDF); xlim (-2.2, 2.2) for sEV / TNT and (-3.2, 3.2) for metabolic.

The individual running-score (enrichment) plots require re-running
gseapy.prerank from the cached GMT; that is NOT done by default.  Pass
``--enrichment-plots`` to regenerate them (same seed, identical NES/FDR).
"""
import os
import sys
import warnings
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import srep_style as srep
srep.apply()
import matplotlib.pyplot as plt

from srep_nes_panel import nes_panel, XLIM

warnings.filterwarnings('ignore')

# ──────────────────────────────────────────────────────────────────
# Configuration (harmonized with OSD-564 / OSD-352; fixed term order)
# ──────────────────────────────────────────────────────────────────
CATEGORIES = {
    'sEV': {
        'title': 'KEGG sEV pathways, OSD-562 (n = 4 vs 3)',
        'xlim':  XLIM['sEV'],
        'pathways': ['Sphingolipid signaling pathway', 'Sphingolipid metabolism',
                     'Phospholipase D signaling pathway', 'Endocytosis',
                     'Synaptic vesicle cycle', 'Phosphatidylinositol signaling system',
                     'Lysosome'],
    },
    'TNT': {
        'title': 'KEGG TNT pathways, OSD-562 (n = 4 vs 3)',
        'xlim':  XLIM['TNT'],
        'pathways': ['Axon guidance', 'Regulation of actin cytoskeleton',
                     'Rap1 signaling pathway'],
    },
    'Metabolic': {
        'title': 'KEGG metabolic adaptation, OSD-562 (n = 4 vs 3)',
        'xlim':  XLIM['metabolic'],
        'pathways': ['AMPK signaling pathway', 'Glycolysis / Gluconeogenesis',
                     'Oxidative phosphorylation'],
    },
}

RESULTS_CSV = 'KEGG_results_OSD562/GSEA_KEGG_results.csv'
OUT_DIR     = 'KEGG_results_OSD562/pathway_figures'
GMT_FILE    = 'KEGG_results_OSD562/GSEA_prerank/gene_sets.gmt'
RANKED_CSV  = 'KEGG_results_OSD562/GSEA_ranked_gene_list.csv'
PERM_NUM    = 1000
os.makedirs(OUT_DIR, exist_ok=True)


def exact_term(name, keys):
    for k in keys:
        if k.strip().lower() == name.lower():
            return k
    return None


# ──────────────────────────────────────────────────────────────────
# 1. Results table (source of truth for NES / FDR)
# ──────────────────────────────────────────────────────────────────
res = pd.read_csv(RESULTS_CSV, index_col=0)
res['Term'] = res['Term'].str.strip()
res = res.set_index('Term')

all_pathways = list(dict.fromkeys(p for c in CATEGORIES.values() for p in c['pathways']))

# ──────────────────────────────────────────────────────────────────
# 2. (optional) running-score plots; requires a prerank re-run from the cached GMT
# ──────────────────────────────────────────────────────────────────
if '--enrichment-plots' in sys.argv:
    import gseapy as gp
    from gseapy import gseaplot
    print("Loading ranked gene list ...")
    rnk = pd.read_csv(RANKED_CSV, index_col=0).squeeze('columns')
    rnk = rnk.sort_values(ascending=False)
    rnk = rnk[~rnk.index.duplicated(keep='first')]
    print(f"  {srep.fmt_int(len(rnk))} genes")
    print(f"Re-running GSEA prerank from cached GMT (permutations = {PERM_NUM}) ...")
    gsea_res = gp.prerank(rnk=rnk, gene_sets=GMT_FILE, outdir=None, seed=42,
                          permutation_num=PERM_NUM, min_size=10, max_size=500,
                          threads=4, verbose=False)
    rerun = gsea_res.res2d.set_index('Term')
    print("\nIndividual enrichment plots ...")
    for term in all_pathways:
        matched = exact_term(term, gsea_res.results.keys())
        if matched is None:
            print(f"  [MISSING] '{term}'")
            continue
        d_nes = abs(float(rerun.at[matched, 'NES']) - float(res.at[matched, 'NES']))
        if d_nes > 1e-6:
            print(f"  [WARN] NES differs between re-run and results table for {term}: {d_nes:.2e}")
        safe = term.replace('/', '-').replace(' ', '_').replace(':', '').replace(',', '')
        ofname = f"{OUT_DIR}/GSEA_{safe}.png"
        try:
            gseaplot(rank_metric=gsea_res.ranking, term=matched, ofname=ofname,
                     **gsea_res.results[matched])
            plt.close('all')
            print(f"  [OK] {term}")
        except Exception as e:
            print(f"  [ERROR] '{term}': {e}")
else:
    print("Skipping running-score plots (no GSEA re-run); pass --enrichment-plots to regenerate them.")

# ──────────────────────────────────────────────────────────────────
# 3. Per-panel NES bar figures (harmonized order, first term at the top)
# ──────────────────────────────────────────────────────────────────
print("\nPanel NES figures ...")
for cat_key, cat in CATEGORIES.items():
    records = []
    for term in cat['pathways']:
        matched = exact_term(term, res.index.tolist())
        if matched is None:
            print(f"  [MISSING] {term}")
            continue
        row = res.loc[matched]
        records.append({'label': term, 'NES': float(row['NES']), 'FDR': float(row['FDR q-val'])})
    if not records:
        continue
    pd.DataFrame(records).to_csv(f'{OUT_DIR}/GSEA_{cat_key}_NES_summary.csv', index=False)

    terms = [{'label': r['label'], 'nes': r['NES'], 'fdr': r['FDR']} for r in records]
    nes_panel(terms, cat['title'], f'{OUT_DIR}/GSEA_{cat_key}_NES_summary', xlim=cat['xlim'])
    print(pd.DataFrame(records).to_string(index=False))

print(f"\nDone. All outputs in: {OUT_DIR}/")
