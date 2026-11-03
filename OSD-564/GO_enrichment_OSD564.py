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
COMP_TAG  = '(Space Flight)v(Ground Control)'   # log2(FLT / GC)
ADJ_THR   = 0.05
FC_THR    = 1.0
GO_DBS    = {
    'BP': 'GO_Biological_Process_2026',
    'CC': 'GO_Cellular_Component_2026',
    'MF': 'GO_Molecular_Function_2026',
}
GO_COLORS = {
    'BP': '#2E7D32',   # green
    'CC': '#F57F17',   # amber
    'MF': '#6A1B9A',   # purple
}
N_TOP     = 15
PERM_NUM  = 1000
OUT_DIR   = 'GO_results'

os.makedirs(OUT_DIR, exist_ok=True)

# ──────────────────────────────────────────────────────────────────
# 1. Load & clean DE data
# ──────────────────────────────────────────────────────────────────
print("Loading DE data …")
de = pd.read_csv(
    'GLDS-569_rna_seq_differential_expression_rRNArm_GLbulkRNAseq.csv',
    index_col=0, low_memory=False,
)

df = de[['SYMBOL',
         f'Log2fc_{COMP_TAG}',
         f'P.value_{COMP_TAG}',
         f'Adj.p.value_{COMP_TAG}']].copy()
df.columns = ['symbol', 'log2fc', 'pval', 'adjpval']

for col in ['log2fc', 'pval', 'adjpval']:
    df[col] = pd.to_numeric(df[col], errors='coerce')

df = df.dropna(subset=['symbol', 'log2fc', 'pval', 'adjpval'])
df = df[df['symbol'].str.strip() != '']
df = df.sort_values('adjpval').drop_duplicates(subset='symbol', keep='first')
print(f"Total genes after deduplication: {len(df):,}")

# ──────────────────────────────────────────────────────────────────
# 2. Gene lists for ORA (adj.p ≤ 0.05, |log2FC| ≥ 1)
# ──────────────────────────────────────────────────────────────────
sig_up   = df[(df['adjpval'] <= ADJ_THR) & (df['log2fc'] >=  FC_THR)]['symbol'].tolist()
sig_down = df[(df['adjpval'] <= ADJ_THR) & (df['log2fc'] <= -FC_THR)]['symbol'].tolist()
sig_all  = sig_up + sig_down
print(f"ORA gene lists — UP: {len(sig_up)}, DOWN: {len(sig_down)}, ALL: {len(sig_all)}")

# ──────────────────────────────────────────────────────────────────
# 3. ORA helper
# ──────────────────────────────────────────────────────────────────
def run_ora(gene_list, label, databases):
    """Run Enrichr ORA for a list of GO databases, return dict keyed by DB abbreviation."""
    if len(gene_list) < 5:
        print(f"  {label}: only {len(gene_list)} genes — skipping.")
        return {}
    results = {}
    for abbr, db_name in databases.items():
        print(f"  {label} / {abbr} ({db_name}) …", end=' ', flush=True)
        for attempt in range(4):
            try:
                time.sleep(4 * (attempt + 1))
                enr = gp.enrichr(
                    gene_list=gene_list,
                    gene_sets=[db_name],
                    organism='mouse',
                    outdir=None,
                    verbose=False,
                )
                break
            except Exception as e:
                if attempt == 3:
                    print(f"FAILED ({e})")
                    enr = None
                    break
                print(f"retry {attempt+1} …", end=' ', flush=True)
        if enr is None:
            continue
        res = enr.results.copy()
        res = res[res['Gene_set'] == db_name].copy()
        res['neg_log10_adjp'] = -np.log10(res['Adjusted P-value'].clip(lower=1e-300))
        res['Gene_count']   = res['Overlap'].str.split('/').str[0].astype(int)
        res['Pathway_size'] = res['Overlap'].str.split('/').str[1].astype(int)
        res['Gene_ratio']   = res['Gene_count'] / res['Pathway_size']
        res = res.sort_values('Adjusted P-value')
        out_path = f'{OUT_DIR}/ORA_{label}_{abbr}.csv'
        res.to_csv(out_path, index=False)
        n_sig = (res['Adjusted P-value'] <= 0.05).sum()
        print(f"{n_sig} sig. / {len(res)} total")
        results[abbr] = res
    return results

# ──────────────────────────────────────────────────────────────────
# 4. Run ORA for each direction
# ──────────────────────────────────────────────────────────────────
print("\n── ORA: Upregulated ──")
ora_up   = run_ora(sig_up,   'Upregulated',    GO_DBS)
print("\n── ORA: Downregulated ──")
ora_down = run_ora(sig_down, 'Downregulated',  GO_DBS)
print("\n── ORA: All significant ──")
ora_all  = run_ora(sig_all,  'All_significant', GO_DBS)

# ──────────────────────────────────────────────────────────────────
# 5. ORA bar plots — one per GO category, showing up & down together
# ──────────────────────────────────────────────────────────────────
def clean_term(term_str):
    """Strip GO ID suffix like (GO:0006412) and trim whitespace."""
    import re
    return re.sub(r'\s*\(GO:\d+\)\s*$', '', term_str).strip()

def plot_ora_bar(res, label, color, title, abbr):
    if res is None or len(res) == 0:
        return
    top = res.head(N_TOP).copy().sort_values('neg_log10_adjp', ascending=True)
    term_labels = [clean_term(t) for t in top['Term']]

    fig, ax = plt.subplots(figsize=(11, max(4, len(top) * 0.50 + 1.2)))
    bars = ax.barh(term_labels, top['neg_log10_adjp'],
                   color=color, edgecolor='white', height=0.72, alpha=0.88)

    for bar, (_, row) in zip(bars, top.iterrows()):
        ax.text(bar.get_width() + 0.04,
                bar.get_y() + bar.get_height() / 2,
                f"n={row['Gene_count']}",
                va='center', ha='left', fontsize=8.5)

    ax.axvline(-np.log10(0.05), color='#9E9E9E', lw=1.2, ls='--', alpha=0.8,
               label='adj.p = 0.05')
    ax.set_xlabel('-log₁₀ (Adjusted P-value)', fontsize=11)
    ax.set_title(f'GO {abbr} ORA — {title}\n(OSD-564: Spaceflight vs Ground Control)',
                 fontsize=11.5, fontweight='bold', pad=10)
    ax.legend(fontsize=9, frameon=False, loc='lower right')
    plt.tight_layout()

    path = f'{OUT_DIR}/ORA_{label}_{abbr}_barplot.png'
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved → {path}")

print("\nPlotting ORA bar plots …")
direction_map = {
    'Upregulated':    (ora_up,   '#C62828', 'Upregulated in Spaceflight'),
    'Downregulated':  (ora_down, '#1565C0', 'Downregulated in Spaceflight'),
    'All_significant':(ora_all,  '#37474F', 'All Significant Genes'),
}
for label, (ora_dict, color, title) in direction_map.items():
    for abbr, go_color in GO_COLORS.items():
        res = ora_dict.get(abbr)
        plot_ora_bar(res, label, go_color, title, abbr)

# ──────────────────────────────────────────────────────────────────
# 6. Combined ORA dot plot — top 10 per GO category (up + down)
# ──────────────────────────────────────────────────────────────────
def build_combined_dot_plot(ora_up, ora_down, abbr, go_color, title_abbr):
    frames = []
    if ora_up.get(abbr) is not None and len(ora_up[abbr]):
        top_u = ora_up[abbr].head(10).copy()
        top_u['Direction'] = 'Upregulated in Spaceflight'
        frames.append(top_u)
    if ora_down.get(abbr) is not None and len(ora_down[abbr]):
        top_d = ora_down[abbr].head(10).copy()
        top_d['Direction'] = 'Downregulated in Spaceflight'
        frames.append(top_d)
    if not frames:
        return

    comb = pd.concat(frames, ignore_index=True)
    comb['Term_short'] = comb['Term'].apply(clean_term)
    palette = {
        'Upregulated in Spaceflight':   '#C62828',
        'Downregulated in Spaceflight': '#1565C0',
    }
    colors = comb['Direction'].map(palette)

    fig, ax = plt.subplots(figsize=(11, max(5, len(comb) * 0.44 + 1.5)))
    ax.scatter(
        comb['Gene_ratio'],
        range(len(comb)),
        c=colors,
        s=comb['Gene_count'] * 14,
        alpha=0.82,
        edgecolors='white',
        linewidths=0.5,
        zorder=3,
    )
    ax.set_yticks(range(len(comb)))
    ax.set_yticklabels(comb['Term_short'], fontsize=8.5)
    ax.set_xlabel('Gene Ratio (overlap / term size)', fontsize=11)
    ax.set_title(f'GO {title_abbr} ORA — Top Enriched Terms\n(OSD-564: Spaceflight vs Ground Control)',
                 fontsize=11.5, fontweight='bold', pad=10)

    up_patch   = mpatches.Patch(color='#C62828', alpha=0.88, label='Upregulated in Spaceflight')
    down_patch = mpatches.Patch(color='#1565C0', alpha=0.88, label='Downregulated in Spaceflight')
    for sz in [5, 15, 30]:
        ax.scatter([], [], s=sz * 14, c='grey', alpha=0.6, label=f'n = {sz}')
    ax.legend(handles=[up_patch, down_patch], fontsize=8.5,
              frameon=True, framealpha=0.9, edgecolor='#cccccc', loc='lower right')
    ax.grid(axis='x', color='#e0e0e0', lw=0.8, zorder=0)
    plt.tight_layout()

    path = f'{OUT_DIR}/ORA_combined_{abbr}_dotplot.png'
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved → {path}")

print("\nPlotting combined ORA dot plots …")
for abbr, go_color in GO_COLORS.items():
    build_combined_dot_plot(ora_up, ora_down, abbr, go_color, abbr)

# ──────────────────────────────────────────────────────────────────
# 7. Ranked gene list for GSEA
#    Metric: sign(log2FC) × −log10(nominal p-value)
# ──────────────────────────────────────────────────────────────────
print("\nPreparing ranked gene list for GSEA …")
df['rank_metric'] = np.sign(df['log2fc']) * (-np.log10(df['pval'].clip(lower=1e-300)))
ranked = (df[['symbol', 'rank_metric']]
          .dropna()
          .sort_values('rank_metric', ascending=False)
          .set_index('symbol')['rank_metric'])
print(f"Ranked list: {len(ranked):,} genes  |  range [{ranked.min():.2f}, {ranked.max():.2f}]")
ranked.to_csv(f'{OUT_DIR}/GSEA_ranked_gene_list.csv', header=True)

# ──────────────────────────────────────────────────────────────────
# 8. GSEA prerank — one run per GO category
# ──────────────────────────────────────────────────────────────────
gsea_results = {}

for abbr, db_name in GO_DBS.items():
    print(f"\nRunning GSEA prerank — {abbr} ({db_name}) …")
    gsea_res = gp.prerank(
        rnk=ranked,
        gene_sets=db_name,
        outdir=f'{OUT_DIR}/GSEA_{abbr}_prerank',
        seed=42,
        permutation_num=PERM_NUM,
        min_size=10,
        max_size=500,
        threads=4,
        verbose=False,
    )
    gsea_df = gsea_res.res2d.copy()
    gsea_df.to_csv(f'{OUT_DIR}/GSEA_{abbr}_results.csv')
    n_fdr05 = (gsea_df['FDR q-val'] <= 0.05).sum()
    n_fdr25 = (gsea_df['FDR q-val'] <= 0.25).sum()
    print(f"  FDR ≤ 0.05: {n_fdr05}  |  FDR ≤ 0.25: {n_fdr25}")
    gsea_results[abbr] = (gsea_res, gsea_df)

# ──────────────────────────────────────────────────────────────────
# 9. GSEA NES summary bar plots — one per GO category
# ──────────────────────────────────────────────────────────────────
print("\nPlotting GSEA NES summaries …")

for abbr, (gsea_res, gsea_df) in gsea_results.items():
    go_color_pos = '#C62828'
    go_color_neg = '#1565C0'

    gsea_df['-log10FDR'] = -np.log10(gsea_df['FDR q-val'].clip(lower=1e-10))
    n_each  = N_TOP // 2
    top_pos = gsea_df[gsea_df['NES'] > 0].nlargest(n_each, 'NES')
    top_neg = gsea_df[gsea_df['NES'] < 0].nsmallest(n_each, 'NES')
    top_g   = (pd.concat([top_neg, top_pos])
               .sort_values('NES', ascending=True)
               .reset_index(drop=True))

    if len(top_g) == 0:
        print(f"  {abbr}: no results to plot.")
        continue

    term_labels = [clean_term(t) for t in top_g['Term']]

    fig, ax = plt.subplots(figsize=(12, max(5, len(top_g) * 0.52 + 1.8)))

    for i, row in top_g.iterrows():
        is_sig = row['FDR q-val'] <= 0.25
        c = go_color_pos if row['NES'] > 0 else go_color_neg
        ax.barh(i, row['NES'], color=c, alpha=1.0 if is_sig else 0.35,
                edgecolor='white', height=0.72)
        fdr_str = f"FDR={row['FDR q-val']:.3f}" + (' *' if row['FDR q-val'] <= 0.05 else '')
        x_off = 0.06 if row['NES'] > 0 else -0.06
        ha    = 'left' if row['NES'] > 0 else 'right'
        ax.text(row['NES'] + x_off, i, fdr_str, va='center', ha=ha, fontsize=8)

    ax.set_yticks(range(len(top_g)))
    ax.set_yticklabels(term_labels, fontsize=9)
    ax.axvline(0, color='black', lw=0.9)
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=11)
    ax.set_title(
        f'GO {abbr} GSEA — OSD-564\nSpaceflight vs Ground Control  (ranked by sign × −log₁₀p)',
        fontsize=11.5, fontweight='bold', pad=10)

    up_patch    = mpatches.Patch(color=go_color_pos, alpha=0.9,  label='Enriched in Spaceflight')
    down_patch  = mpatches.Patch(color=go_color_neg, alpha=0.9,  label='Enriched in Ground Control')
    faded_patch = mpatches.Patch(color='grey',       alpha=0.35, label='FDR > 0.25 (not sig.)')
    ax.legend(handles=[up_patch, down_patch, faded_patch],
              fontsize=9, frameon=True, framealpha=0.9, edgecolor='#cccccc', loc='lower right')

    plt.tight_layout()
    path = f'{OUT_DIR}/GSEA_{abbr}_NES_barplot.png'
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved → {path}")

# ──────────────────────────────────────────────────────────────────
# 10. Individual GSEA enrichment plots for top 6 per GO category
# ──────────────────────────────────────────────────────────────────
print("\nGenerating individual GSEA enrichment plots …")

for abbr, (gsea_res, gsea_df) in gsea_results.items():
    top_terms = gsea_df.nsmallest(6, 'FDR q-val')['Term'].tolist()
    if not top_terms:
        top_terms = gsea_df.nlargest(6, 'NES')['Term'].tolist()
    for term in top_terms:
        try:
            safe = term[:65].replace('/', '-').replace(' ', '_').replace(':', '').replace('(', '').replace(')', '')
            ofname = f"{OUT_DIR}/GSEA_{abbr}_enrich_{safe}.png"
            gseaplot(
                rank_metric=gsea_res.ranking,
                term=term,
                ofname=ofname,
                **gsea_res.results[term],
            )
            plt.close('all')
            print(f"  [{abbr}] {term[:65]}")
        except Exception as e:
            print(f"  [{abbr}] Could not plot '{term[:55]}': {e}")

# ──────────────────────────────────────────────────────────────────
# 11. Combined Excel summary
# ──────────────────────────────────────────────────────────────────
print("\nWriting Excel summary …")

xl_path = f'{OUT_DIR}/GO_enrichment_summary.xlsx'
with pd.ExcelWriter(xl_path, engine='openpyxl') as writer:
    for label, ora_dict in [('UP', ora_up), ('DOWN', ora_down), ('ALL', ora_all)]:
        for abbr, res in ora_dict.items():
            sheet = f'ORA_{label}_{abbr}'[:31]
            res.to_excel(writer, sheet_name=sheet, index=False)
    for abbr, (_, gsea_df) in gsea_results.items():
        sheet = f'GSEA_{abbr}'[:31]
        gsea_df.sort_values('FDR q-val').to_excel(writer, sheet_name=sheet)

print(f"  Saved → {xl_path}")
print("\nDone. All outputs in:", OUT_DIR)
