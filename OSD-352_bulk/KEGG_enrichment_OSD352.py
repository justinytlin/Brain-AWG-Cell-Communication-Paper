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
# NOTE: use the (Space Flight)v(Ground Control) contrast so that positive log2FC = up in
#       spaceflight (matches OSD-562/OSD-564). The reversed tag '(Ground Control)v(Space Flight)'
#       flips every NES sign.
COMP_TAG  = '(Space Flight)v(Ground Control)'   # log2(FLT / GC)
ADJ_THR   = 0.05
FC_THR    = 1.0
KEGG_DB   = 'KEGG_2019_Mouse'
N_TOP     = 15
PERM_NUM  = 1000
OUT_DIR   = 'KEGG_results_OSD352'

os.makedirs(OUT_DIR, exist_ok=True)

# ──────────────────────────────────────────────────────────────────
# 1. Load & clean DE data
# ──────────────────────────────────────────────────────────────────
print("Loading DE data …")
de = pd.read_csv(
    'GLDS-352_rna_seq_differential_expression.csv',
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

# Deduplicate: keep lowest adj p-value per symbol
df = df.sort_values('adjpval').drop_duplicates(subset='symbol', keep='first')
print(f"Total genes after deduplication: {len(df):,}")

# ──────────────────────────────────────────────────────────────────
# 2. ORA gene lists (adj.p ≤ 0.05, |log2FC| ≥ 1)
# ──────────────────────────────────────────────────────────────────
sig_up   = df[(df['adjpval'] <= ADJ_THR) & (df['log2fc'] >=  FC_THR)]['symbol'].tolist()
sig_down = df[(df['adjpval'] <= ADJ_THR) & (df['log2fc'] <= -FC_THR)]['symbol'].tolist()
sig_all  = sig_up + sig_down

print(f"\nORA gene lists — UP: {len(sig_up)}, DOWN: {len(sig_down)}, ALL: {len(sig_all)}")

# ──────────────────────────────────────────────────────────────────
# 3. ORA with Enrichr / KEGG
# ──────────────────────────────────────────────────────────────────
print("\nRunning ORA (Enrichr KEGG) …")

def run_ora(gene_list, label):
    if len(gene_list) < 5:
        print(f"  {label}: only {len(gene_list)} genes — skipping.")
        return None
    for attempt in range(4):
        try:
            time.sleep(3 * (attempt + 1))   # back-off: 3 s, 6 s, 9 s, 12 s
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=[KEGG_DB],
                organism='mouse',
                outdir=None,
                verbose=False,
            )
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
        ax.text(bar.get_width() + 0.04,
                bar.get_y() + bar.get_height() / 2,
                f"n={row['Gene_count']}",
                va='center', ha='left', fontsize=8.5)

    threshold_x = -np.log10(0.05)
    ax.axvline(threshold_x, color='#9E9E9E', lw=1.2, ls='--', alpha=0.8,
               label='adj.p = 0.05')
    ax.set_xlabel('-log₁₀ (Adjusted P-value)', fontsize=11)
    ax.set_title(f'KEGG Pathway ORA — {title}\n(OSD-352: Spaceflight vs Ground Control)',
                 fontsize=11.5, fontweight='bold', pad=10)
    ax.legend(fontsize=9, frameon=False, loc='lower right')
    plt.tight_layout()

    path = f'{OUT_DIR}/ORA_{label}_barplot.png'
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved → {path}")

print("\nPlotting ORA bar plots …")
plot_ora_bar(ora_up,   'Upregulated',    '#C62828', 'Upregulated in Spaceflight (adj.p ≤ 0.05, log2FC ≥ 1)')
plot_ora_bar(ora_down, 'Downregulated',  '#1565C0', 'Downregulated in Spaceflight (adj.p ≤ 0.05, log2FC ≤ -1)')
plot_ora_bar(ora_all,  'All_significant','#6A1B9A', 'All Significant Genes (adj.p ≤ 0.05, |log2FC| ≥ 1)')

# ──────────────────────────────────────────────────────────────────
# 5. ORA dot plot (combined up & down, top 10 each)
# ──────────────────────────────────────────────────────────────────
def build_combined_dot_plot(ora_up, ora_down):
    frames = []
    if ora_up is not None and len(ora_up):
        top_u = ora_up.head(10).copy()
        top_u['Direction'] = 'Upregulated in Spaceflight'
        frames.append(top_u)
    if ora_down is not None and len(ora_down):
        top_d = ora_down.head(10).copy()
        top_d['Direction'] = 'Downregulated in Spaceflight'
        frames.append(top_d)
    if not frames:
        return

    comb = pd.concat(frames, ignore_index=True)
    comb['Term_short'] = comb['Term'].str.replace(r'\s*KEGG.*$', '', regex=True).str.strip()

    palette = {'Upregulated in Spaceflight': '#C62828', 'Downregulated in Spaceflight': '#1565C0'}
    colors  = comb['Direction'].map(palette)

    fig, ax = plt.subplots(figsize=(11, max(5, len(comb) * 0.42 + 1.5)))
    sc = ax.scatter(
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
    ax.set_yticklabels(comb['Term_short'], fontsize=9)
    ax.set_xlabel('Gene Ratio (overlap / pathway size)', fontsize=11)
    ax.set_title('KEGG Pathway ORA — Top Enriched Pathways\n(OSD-352: Spaceflight vs Ground Control)',
                 fontsize=11.5, fontweight='bold', pad=10)

    up_patch   = mpatches.Patch(color='#C62828', alpha=0.88, label='Upregulated in Spaceflight')
    down_patch = mpatches.Patch(color='#1565C0', alpha=0.88, label='Downregulated in Spaceflight')

    # Size legend
    for sz in [5, 15, 30]:
        ax.scatter([], [], s=sz * 14, c='grey', alpha=0.6, label=f'n = {sz}')

    handles_extra = [
        mpatches.Patch(color='grey', alpha=0, label=''),
        up_patch, down_patch,
    ]
    ax.legend(fontsize=8.5, frameon=True, framealpha=0.9, edgecolor='#cccccc',
              loc='lower right', labelspacing=0.8)
    ax.grid(axis='x', color='#e0e0e0', lw=0.8, zorder=0)
    plt.tight_layout()

    path = f'{OUT_DIR}/ORA_combined_dotplot.png'
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved → {path}")

print("\nPlotting combined ORA dot plot …")
build_combined_dot_plot(ora_up, ora_down)

# ──────────────────────────────────────────────────────────────────
# 6. Ranked gene list for GSEA
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
# 7. GSEA prerank
# ──────────────────────────────────────────────────────────────────
print(f"\nRunning GSEA prerank (permutations = {PERM_NUM}) …")

gsea_res = gp.prerank(
    rnk=ranked,
    gene_sets=KEGG_DB,
    outdir=f'{OUT_DIR}/GSEA_prerank',
    seed=42,
    permutation_num=PERM_NUM,
    min_size=10,
    max_size=500,
    threads=4,
    verbose=False,
)

gsea_df = gsea_res.res2d.copy()
gsea_df.to_csv(f'{OUT_DIR}/GSEA_KEGG_results.csv')

n_fdr25 = (gsea_df['FDR q-val'] <= 0.25).sum()
n_fdr05 = (gsea_df['FDR q-val'] <= 0.05).sum()
print(f"GSEA done.  FDR ≤ 0.05: {n_fdr05}  |  FDR ≤ 0.25: {n_fdr25}")

# ──────────────────────────────────────────────────────────────────
# 8. GSEA NES summary bar plot
# ──────────────────────────────────────────────────────────────────
print("\nPlotting GSEA NES summary …")

gsea_df['-log10FDR'] = -np.log10(gsea_df['FDR q-val'].clip(lower=1e-10))

n_each = N_TOP // 2
top_pos = gsea_df[gsea_df['NES'] > 0].nlargest(n_each, 'NES')
top_neg = gsea_df[gsea_df['NES'] < 0].nsmallest(n_each, 'NES')
top_gsea = pd.concat([top_neg, top_pos]).sort_values('NES', ascending=True).reset_index(drop=True)

if len(top_gsea) > 0:
    term_labels = top_gsea['Term'].str.strip().tolist()

    fig, ax = plt.subplots(figsize=(11, max(5, len(top_gsea) * 0.52 + 1.8)))

    for i, row in top_gsea.iterrows():
        is_sig = row['FDR q-val'] <= 0.25
        c = '#C62828' if row['NES'] > 0 else '#1565C0'
        ax.barh(i, row['NES'], color=c, alpha=1.0 if is_sig else 0.35,
                edgecolor='white', height=0.72)
        fdr_label = f"FDR={row['FDR q-val']:.3f}" + (' *' if row['FDR q-val'] <= 0.05 else '')
        x_off = 0.06 if row['NES'] > 0 else -0.06
        ha    = 'left'   if row['NES'] > 0 else 'right'
        ax.text(row['NES'] + x_off, i, fdr_label,
                va='center', ha=ha, fontsize=8)

    ax.set_yticks(range(len(top_gsea)))
    ax.set_yticklabels(term_labels, fontsize=9)
    ax.axvline(0, color='black', lw=0.9)
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=11)
    ax.set_title('KEGG Pathway GSEA, OSD-352 whole-brain bulk RNA-seq\nSpaceflight vs ground control (n = 3 vs 3; ranked by sign × −log₁₀p)',
                 fontsize=11.5, fontweight='bold', pad=10)

    up_patch    = mpatches.Patch(color='#C62828', alpha=0.9,  label='Enriched in spaceflight (NES > 0)')
    down_patch  = mpatches.Patch(color='#1565C0', alpha=0.9,  label='Enriched in ground control (NES < 0)')
    faded_patch = mpatches.Patch(color='grey',    alpha=0.35, label='FDR > 0.25 (not significant)')
    ax.legend(handles=[up_patch, down_patch, faded_patch],
              fontsize=9, frameon=True, framealpha=0.9, edgecolor='#cccccc',
              loc='lower right')

    plt.tight_layout()
    path = f'{OUT_DIR}/GSEA_NES_barplot.png'
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved → {path}")

# ──────────────────────────────────────────────────────────────────
# 9. Individual GSEA enrichment plots for top hits
# ──────────────────────────────────────────────────────────────────
print("\nGenerating individual GSEA enrichment plots …")

top_terms = gsea_df.nsmallest(6, 'FDR q-val')['Term'].tolist()
if not top_terms:
    top_terms = gsea_df.nlargest(6, 'NES')['Term'].tolist()

for term in top_terms:
    try:
        safe_name = term[:70].replace('/', '-').replace(' ', '_').replace(':', '')
        ofname = f"{OUT_DIR}/GSEA_enrich_{safe_name}.png"
        gseaplot(
            rank_metric=gsea_res.ranking,
            term=term,
            ofname=ofname,
            **gsea_res.results[term],
        )
        plt.close('all')
        print(f"  Saved enrichment plot: {term[:70]}")
    except Exception as e:
        print(f"  Could not plot '{term[:60]}': {e}")

# ──────────────────────────────────────────────────────────────────
# 10. Save combined Excel summary
# ──────────────────────────────────────────────────────────────────
print("\nWriting Excel summary …")

xl_path = f'{OUT_DIR}/KEGG_enrichment_summary.xlsx'
with pd.ExcelWriter(xl_path, engine='openpyxl') as writer:
    if ora_up is not None:
        ora_up.to_excel(writer, sheet_name='ORA_Upregulated', index=False)
    if ora_down is not None:
        ora_down.to_excel(writer, sheet_name='ORA_Downregulated', index=False)
    if ora_all is not None:
        ora_all.to_excel(writer, sheet_name='ORA_All_significant', index=False)
    gsea_df.sort_values('FDR q-val').to_excel(writer, sheet_name='GSEA_results')

print(f"  Saved → {xl_path}")
print("\nDone. All outputs in:", OUT_DIR)
