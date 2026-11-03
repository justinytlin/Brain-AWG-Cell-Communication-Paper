"""
KEGG_enrichment_OSD562.py
-------------------------
KEGG (KEGG_2019_Mouse) ORA + GSEA-prerank for OSD-562, re-run DESeq2
contrast ISS-T FLT 20 wk (n = 4) vs ISS-T GC 20 wk (n = 3)
(OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3.csv).

GSEA prerank: ranking metric sign(log2FC) * -log10(nominal p), 1000
permutations, seed 42, min size 10, max size 500.
"""
import os
import time
import warnings
import numpy as np
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
# Configuration
# ──────────────────────────────────────────────────────────────────
DE_FILE   = 'OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3.csv'
COMP_TAG  = '(Space Flight)v(Ground Control)'
ADJ_THR   = 0.05
FC_THR    = 1.0
KEGG_DB   = 'KEGG_2019_Mouse'
N_TOP     = 15
PERM_NUM  = 1000
OUT_DIR   = 'KEGG_results_OSD562'
SUBTITLE  = 'OSD-562: spaceflight vs ground control (n = 4 vs 3)'

os.makedirs(OUT_DIR, exist_ok=True)

# ──────────────────────────────────────────────────────────────────
# 1. Load & clean DE data
# ──────────────────────────────────────────────────────────────────
print("Loading DE data (re-run 4 vs 3) …")
de = pd.read_csv(DE_FILE, index_col=0, low_memory=False)

df = de[['SYMBOL',
         f'Log2fc_{COMP_TAG}',
         f'P.value_{COMP_TAG}',
         f'Adj.p.value_{COMP_TAG}']].copy()
df.columns = ['symbol', 'log2fc', 'pval', 'adjpval']
for col in ['log2fc', 'pval', 'adjpval']:
    df[col] = pd.to_numeric(df[col], errors='coerce')
df = df.dropna(subset=['symbol', 'log2fc', 'pval', 'adjpval'])
df = df[df['symbol'].str.strip() != '']

# Deduplicate: keep lowest adj p (then nominal p) per symbol
df = df.sort_values(['adjpval', 'pval']).drop_duplicates(subset='symbol', keep='first')
print(f"Total genes after deduplication: {len(df):,}")

# ──────────────────────────────────────────────────────────────────
# 2. ORA gene lists (adj.p ≤ 0.05, |log2FC| ≥ 1)
# ──────────────────────────────────────────────────────────────────
sig_up   = df[(df['adjpval'] <= ADJ_THR) & (df['log2fc'] >=  FC_THR)]['symbol'].tolist()
sig_down = df[(df['adjpval'] <= ADJ_THR) & (df['log2fc'] <= -FC_THR)]['symbol'].tolist()
sig_all  = sig_up + sig_down
print(f"\nORA gene lists — UP: {len(sig_up)}, DOWN: {len(sig_down)}, ALL: {len(sig_all)}")

# ──────────────────────────────────────────────────────────────────
# 3. ORA with Enrichr / KEGG  (skipped automatically if < 5 genes)
# ──────────────────────────────────────────────────────────────────
print("\nRunning ORA (Enrichr KEGG) …")


def run_ora(gene_list, label):
    if len(gene_list) < 5:
        print(f"  {label}: only {len(gene_list)} genes — skipping.")
        return None
    for attempt in range(4):
        try:
            time.sleep(3 * (attempt + 1))
            enr = gp.enrichr(gene_list=gene_list, gene_sets=[KEGG_DB],
                             organism='mouse', outdir=None, verbose=False)
            break
        except Exception as e:
            if attempt == 3:
                print(f"  {label}: Enrichr failed after 4 attempts — {e}")
                return None
            print(f"  {label}: attempt {attempt+1} failed ({e}), retrying …")
    res = enr.results.copy()
    res = res[res['Gene_set'] == KEGG_DB].copy()
    res['neg_log10_adjp'] = -np.log10(res['Adjusted P-value'].clip(lower=1e-300))
    res['Gene_count'] = res['Overlap'].str.split('/').str[0].astype(int)
    res['Pathway_size'] = res['Overlap'].str.split('/').str[1].astype(int)
    res['Gene_ratio'] = res['Gene_count'] / res['Pathway_size']
    res = res.sort_values('Adjusted P-value')
    out_path = f'{OUT_DIR}/ORA_{label}_KEGG.csv'
    res.to_csv(out_path, index=False)
    n_sig = (res['Adjusted P-value'] <= 0.05).sum()
    print(f"  {label}: {n_sig} sig. pathways (adj.p ≤ 0.05) / {len(res)} total → {out_path}")
    return res


ora_up   = run_ora(sig_up,   'Upregulated')
ora_down = run_ora(sig_down, 'Downregulated')
ora_all  = run_ora(sig_all,  'All_significant')

# ──────────────────────────────────────────────────────────────────
# 4. ORA bar plots
# ──────────────────────────────────────────────────────────────────
def plot_ora_bar(res, label, color, title):
    if res is None or len(res) == 0:
        return
    top = res.head(N_TOP).copy().sort_values('neg_log10_adjp', ascending=True)
    term_labels = top['Term'].str.replace(r'\s*KEGG.*$', '', regex=True).str.strip()
    fig, ax = plt.subplots(figsize=(10, max(4, len(top) * 0.50 + 1.2)))
    bars = ax.barh(term_labels, top['neg_log10_adjp'],
                   color=color, edgecolor='white', height=0.72, alpha=0.88)
    for bar, (_, row) in zip(bars, top.iterrows()):
        ax.text(bar.get_width() + 0.04, bar.get_y() + bar.get_height() / 2,
                f"n={row['Gene_count']}", va='center', ha='left', fontsize=8.5)
    ax.axvline(-np.log10(0.05), color='#9E9E9E', lw=1.2, ls='--', alpha=0.8, label='adj.p = 0.05')
    ax.set_xlabel('-log10 (adjusted p-value)', fontsize=11)
    ax.set_title(f'KEGG pathway ORA, {title}\n{SUBTITLE}', fontsize=11.5, fontweight='bold', pad=10)
    ax.legend(fontsize=9, frameon=False, loc='lower right')
    plt.tight_layout()
    path = f'{OUT_DIR}/ORA_{label}_barplot.png'
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved → {path}")


print("\nPlotting ORA bar plots …")
plot_ora_bar(ora_up,   'Upregulated',     '#C62828', 'up in spaceflight (adj p<=0.05, log2FC>=1)')
plot_ora_bar(ora_down, 'Downregulated',   '#1565C0', 'down in spaceflight (adj p<=0.05, log2FC<=-1)')
plot_ora_bar(ora_all,  'All_significant', '#6A1B9A', 'all significant genes (adj p<=0.05, |log2FC|>=1)')

# ──────────────────────────────────────────────────────────────────
# 5. Ranked gene list for GSEA: sign(log2FC) × −log10(nominal p)
# ──────────────────────────────────────────────────────────────────
print("\nPreparing ranked gene list for GSEA …")
df['rank_metric'] = np.sign(df['log2fc']) * (-np.log10(df['pval'].clip(lower=1e-300)))
ranked = (df[['symbol', 'rank_metric']].dropna()
          .sort_values('rank_metric', ascending=False)
          .set_index('symbol')['rank_metric'])
print(f"Ranked list: {len(ranked):,} genes  |  range [{ranked.min():.2f}, {ranked.max():.2f}]")
ranked.to_csv(f'{OUT_DIR}/GSEA_ranked_gene_list.csv', header=True)

# ──────────────────────────────────────────────────────────────────
# 6. GSEA prerank
# ──────────────────────────────────────────────────────────────────
print(f"\nRunning GSEA prerank (permutations = {PERM_NUM}) …")
gsea_res = gp.prerank(
    rnk=ranked, gene_sets=KEGG_DB, outdir=f'{OUT_DIR}/GSEA_prerank',
    seed=42, permutation_num=PERM_NUM, min_size=10, max_size=500,
    threads=4, verbose=False,
)
gsea_df = gsea_res.res2d.copy()
gsea_df.to_csv(f'{OUT_DIR}/GSEA_KEGG_results.csv')
n_fdr25 = (gsea_df['FDR q-val'] <= 0.25).sum()
n_fdr05 = (gsea_df['FDR q-val'] <= 0.05).sum()
print(f"GSEA done.  {len(gsea_df)} gene sets | FDR ≤ 0.05: {n_fdr05}  |  FDR ≤ 0.25: {n_fdr25}")

# ──────────────────────────────────────────────────────────────────
# 7. GSEA NES summary bar plot (top |NES| each direction)
# ──────────────────────────────────────────────────────────────────
print("\nPlotting GSEA NES summary …")


def fdr_label(q):
    return 'FDR < 0.001' if q < 0.001 else f'FDR = {q:.3f}'


n_each = N_TOP // 2
top_pos = gsea_df[gsea_df['NES'] > 0].nlargest(n_each, 'NES')
top_neg = gsea_df[gsea_df['NES'] < 0].nsmallest(n_each, 'NES')
top_gsea = pd.concat([top_neg, top_pos]).sort_values('NES', ascending=True).reset_index(drop=True)

if len(top_gsea) > 0:
    fig, ax = plt.subplots(figsize=(11, max(5, len(top_gsea) * 0.52 + 1.8)))
    for i, row in top_gsea.iterrows():
        q = float(row['FDR q-val'])
        c = '#C62828' if row['NES'] > 0 else '#1565C0'
        ax.barh(i, row['NES'], color=c, alpha=1.0 if q <= 0.05 else 0.45,
                edgecolor='white', height=0.72)
        lab = fdr_label(q) + (' *' if q <= 0.05 else '')
        ax.text(row['NES'] + (0.06 if row['NES'] > 0 else -0.06), i, lab,
                va='center', ha='left' if row['NES'] > 0 else 'right', fontsize=8)
    ax.set_yticks(range(len(top_gsea)))
    ax.set_yticklabels(top_gsea['Term'].str.strip().tolist(), fontsize=9)
    ax.axvline(0, color='black', lw=0.9)
    lim = top_gsea['NES'].abs().max() * 1.45
    ax.set_xlim(-lim, lim)
    ax.set_xlabel('Normalized enrichment score (NES)', fontsize=11)
    ax.set_title(f'KEGG pathway GSEA, {SUBTITLE}\n(ranked by sign(log2FC) x -log10 p)',
                 fontsize=11.5, fontweight='bold', pad=10)
    handles = [mpatches.Patch(color='#C62828', alpha=1.0,  label='Enriched in spaceflight (NES > 0)'),
               mpatches.Patch(color='#1565C0', alpha=1.0,  label='Enriched in ground control (NES < 0)'),
               mpatches.Patch(color='grey',    alpha=1.0,  label='FDR <= 0.05 (*)'),
               mpatches.Patch(color='grey',    alpha=0.45, label='FDR > 0.05')]
    ax.legend(handles=handles, fontsize=9, frameon=True, framealpha=0.9,
              edgecolor='#cccccc', loc='lower right')
    plt.tight_layout()
    path = f'{OUT_DIR}/GSEA_NES_barplot.png'
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved → {path}")

# ──────────────────────────────────────────────────────────────────
# 8. Individual GSEA enrichment plots for top hits
# ──────────────────────────────────────────────────────────────────
print("\nGenerating individual GSEA enrichment plots …")
top_terms = gsea_df.nsmallest(6, 'FDR q-val')['Term'].tolist()
for term in top_terms:
    try:
        safe_name = term[:70].replace('/', '-').replace(' ', '_').replace(':', '')
        ofname = f"{OUT_DIR}/GSEA_enrich_{safe_name}.png"
        gseaplot(rank_metric=gsea_res.ranking, term=term, ofname=ofname,
                 **gsea_res.results[term])
        plt.close('all')
        print(f"  Saved enrichment plot: {term[:70]}")
    except Exception as e:
        print(f"  Could not plot '{term[:60]}': {e}")

# ──────────────────────────────────────────────────────────────────
# 9. Save combined Excel summary
# ──────────────────────────────────────────────────────────────────
print("\nWriting Excel summary …")
xl_path = f'{OUT_DIR}/KEGG_enrichment_summary.xlsx'
with pd.ExcelWriter(xl_path, engine='openpyxl') as writer:
    for name, r in [('ORA_Upregulated', ora_up), ('ORA_Downregulated', ora_down),
                    ('ORA_All_significant', ora_all)]:
        if r is not None:
            r.to_excel(writer, sheet_name=name, index=False)
    gsea_df.sort_values('FDR q-val').to_excel(writer, sheet_name='GSEA_results')
print(f"  Saved → {xl_path}")
print("\nDone. All outputs in:", OUT_DIR)
