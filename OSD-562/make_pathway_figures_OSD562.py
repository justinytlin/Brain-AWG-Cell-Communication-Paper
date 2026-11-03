"""
make_pathway_figures_OSD562.py
------------------------------
KEGG GSEA panel figures for OSD-562 (re-run DESeq2, ISS-T FLT n = 4 vs
ISS-T GC n = 3), using the HARMONIZED panel sets shared with OSD-564:

  * KEGG sEV (7), KEGG TNT (3), KEGG metabolic adaptation (3, supplementary)

Terms are matched by EXACT KEGG name (case-insensitive) in
KEGG_results_OSD562/GSEA_KEGG_results.csv; every term is drawn with its FDR
label whether or not significant (asterisk at FDR <= 0.05).  Individual
enrichment (running-score) plots are regenerated from the cached GMT with
the same seed, so NES/FDR are identical to the results table.
"""
import os
import warnings
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import gseapy as gp
from gseapy import gseaplot

warnings.filterwarnings('ignore')
plt.rcParams.update({
    'font.family': 'Arial',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300,
})

# ──────────────────────────────────────────────────────────────────
# Configuration (harmonized with OSD-564)
# ──────────────────────────────────────────────────────────────────
SUBTITLE = 'OSD-562: spaceflight vs ground control (n = 4 vs 3)'
CATEGORIES = {
    'sEV': {
        'panel': 'KEGG sEV pathways',
        'pathways': ['Sphingolipid signaling pathway', 'Sphingolipid metabolism',
                     'Phospholipase D signaling pathway', 'Endocytosis',
                     'Synaptic vesicle cycle', 'Phosphatidylinositol signaling system',
                     'Lysosome'],
    },
    'TNT': {
        'panel': 'KEGG TNT pathways',
        'pathways': ['Axon guidance', 'Regulation of actin cytoskeleton',
                     'Rap1 signaling pathway'],
    },
    'Metabolic': {
        'panel': 'KEGG metabolic adaptation pathways',
        'pathways': ['AMPK signaling pathway', 'Glycolysis / Gluconeogenesis',
                     'Oxidative phosphorylation'],
    },
}

RESULTS_CSV = 'KEGG_results_OSD562/GSEA_KEGG_results.csv'
OUT_DIR     = 'KEGG_results_OSD562/pathway_figures'
GMT_FILE    = 'KEGG_results_OSD562/GSEA_prerank/gene_sets.gmt'
RANKED_CSV  = 'KEGG_results_OSD562/GSEA_ranked_gene_list.csv'
PERM_NUM    = 1000
RED, BLUE   = '#C62828', '#1565C0'
os.makedirs(OUT_DIR, exist_ok=True)


def fdr_label(q):
    return 'FDR < 0.001' if q < 0.001 else f'FDR = {q:.3f}'


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

# ──────────────────────────────────────────────────────────────────
# 2. Re-run prerank from cached GMT for the running-score plots
# ──────────────────────────────────────────────────────────────────
print("Loading ranked gene list ...")
rnk = pd.read_csv(RANKED_CSV, index_col=0).squeeze('columns')
rnk = rnk.sort_values(ascending=False)
rnk = rnk[~rnk.index.duplicated(keep='first')]
print(f"  {len(rnk):,} genes")
print(f"Re-running GSEA prerank from cached GMT (permutations = {PERM_NUM}) ...")
gsea_res = gp.prerank(rnk=rnk, gene_sets=GMT_FILE, outdir=None, seed=42,
                      permutation_num=PERM_NUM, min_size=10, max_size=500,
                      threads=4, verbose=False)
rerun = gsea_res.res2d.set_index('Term')

all_pathways = list(dict.fromkeys(p for c in CATEGORIES.values() for p in c['pathways']))
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

# ──────────────────────────────────────────────────────────────────
# 3. Per-panel NES bar figures
# ──────────────────────────────────────────────────────────────────
legend_handles = [
    mpatches.Patch(color=RED,  alpha=1.0,  label='Enriched in spaceflight (NES > 0)'),
    mpatches.Patch(color=BLUE, alpha=1.0,  label='Enriched in ground control (NES < 0)'),
    mpatches.Patch(color='grey', alpha=1.0,  label='FDR <= 0.05 (*)'),
    mpatches.Patch(color='grey', alpha=0.45, label='FDR > 0.05'),
]

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
    plot_df = pd.DataFrame(records).sort_values('NES').reset_index(drop=True)
    pd.DataFrame(records).to_csv(f'{OUT_DIR}/GSEA_{cat_key}_NES_summary.csv', index=False)

    n = len(plot_df)
    fig, ax = plt.subplots(figsize=(10, max(3.2, n * 0.6 + 2.0)))
    for i, row in plot_df.iterrows():
        sig = row['FDR'] <= 0.05
        c = RED if row['NES'] > 0 else BLUE
        ax.barh(i, row['NES'], color=c, alpha=1.0 if sig else 0.45, edgecolor='white', height=0.68)
        lab = fdr_label(row['FDR']) + (' *' if sig else '')
        off = 0.05 if row['NES'] >= 0 else -0.05
        ax.text(row['NES'] + off, i, lab, va='center',
                ha='left' if row['NES'] >= 0 else 'right', fontsize=9, color='black')
    ax.set_yticks(range(n))
    ax.set_yticklabels(plot_df['label'], fontsize=11)
    ax.axvline(0, color='#333333', lw=1.0)
    lim = max(1.0, plot_df['NES'].abs().max()) * 1.6
    ax.set_xlim(-lim, lim)
    ax.xaxis.grid(True, color='#e8e8e8', lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.set_xlabel('Normalized enrichment score (NES)', fontsize=12)
    ax.set_title(f"{cat['panel']}, {SUBTITLE}", fontsize=12, fontweight='bold', pad=10)
    ax.legend(handles=legend_handles, fontsize=9, frameon=True, framealpha=0.92,
              edgecolor='#cccccc', loc='lower right', labelspacing=0.5)
    plt.tight_layout()
    out = f'{OUT_DIR}/GSEA_{cat_key}_NES_summary.png'
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved → {out}")
    print(plot_df.to_string(index=False))

print(f"\nDone. All outputs in: {OUT_DIR}/")
