"""
OSD-564 RNA-seq QC summary (12 metrics, 10 libraries).

Laid out 2 columns x 6 rows with large fonts so metric names and sample IDs stay
legible when the figure is reduced to single-column width.  Dashed lines mark the
FLT and GC group means and are labeled with the mean value.
"""
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'Arial',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300,
})

# ──────────────────────────────────────────────────────────────────
# Load & annotate
# ──────────────────────────────────────────────────────────────────
df = pd.read_csv('GLDS-569_rna_seq_qc_metrics_GLbulkRNAseq.csv')

df['condition'] = df['sample'].apply(lambda s: 'FLT' if '_FLT_' in s else 'GC')
# FLT libraries first, then GC, each in sample order
df = (df.assign(_o=df['condition'].map({'FLT': 0, 'GC': 1}))
        .sort_values(['_o', 'sample']).drop(columns='_o').reset_index(drop=True))
# short sample IDs: F1 ... F9 (FLT), G1 ... G9 (GC)
df['label'] = df['sample'].apply(lambda s: s.split('_')[-1])

COLORS = {'FLT': '#E07B39', 'GC': '#3A7DC9'}

# ──────────────────────────────────────────────────────────────────
# Panel definitions: (column, y-label, title, optional scale fn, value fmt)
# ──────────────────────────────────────────────────────────────────
panels = [
    ('read_depth',                'Read depth (M)',       'Read depth',                    lambda x: x / 1e6, '{:.0f}'),
    ('rin',                       'RIN',                  'RNA integrity number (RIN)',    None, '{:.1f}'),
    ('rrna_contamination',        'rRNA (%)',             'rRNA fraction',                 None, '{:.1f}'),
    ('raw_percent_duplicates_f',  'Duplicates (%)',       'Duplicate reads (raw, R1)',     None, '{:.1f}'),
    ('gene_detected_gt10',        'Genes detected',       'Genes detected (>10 counts)',   None, '{:,.0f}'),
    ('uniquely_mapped_percent',   'Uniquely mapped (%)',  'Unique mapping rate (STAR)',    None, '{:.1f}'),
    ('multimapped_percent',       'Multi-mapped (%)',     'Multi-mapping rate (STAR)',     None, '{:.1f}'),
    ('unmapped_tooshort_percent', 'Unmapped (%)',         'Unmapped, too short (STAR)',    None, '{:.1f}'),
    ('mean_genebody_cov_5_20',    "5' coverage (%)",      "Gene-body coverage, 5' end",    None, '{:.1f}'),
    ('mean_genebody_cov_80_95',   "3' coverage (%)",      "Gene-body coverage, 3' end",    None, '{:.1f}'),
    ('ratio_genebody_cov_3_to_5', "3'/5' ratio",          "Gene-body 3'/5' bias",          None, '{:.2f}'),
    ('pct_sense',                 'Sense reads (%)',      'Strandedness (% sense)',        None, '{:.1f}'),
]

N_COLS = 2
N_ROWS = int(np.ceil(len(panels) / N_COLS))   # 6 rows

FS_TITLE, FS_LABEL, FS_TICK, FS_MEAN = 14, 12.5, 12, 10.5

fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(9.5, 3.1 * N_ROWS))
axes = axes.flatten()

x = np.arange(len(df))
bar_w = 0.72

for i, (col, ylabel, title, scale_fn, fmt) in enumerate(panels):
    ax = axes[i]
    vals = df[col].astype(float).copy()
    if scale_fn:
        vals = vals.apply(scale_fn)

    bar_colors = [COLORS[c] for c in df['condition']]
    ax.bar(x, vals, color=bar_colors, edgecolor='white', linewidth=0.4,
           width=bar_w, zorder=2)

    # Group-mean dashed lines, labeled with the mean value
    y_lo, y_hi = vals.min(), vals.max()
    for cond in ['FLT', 'GC']:
        idx = df.index[df['condition'] == cond].tolist()
        mean_val = vals.iloc[idx].mean()
        x_min = min(idx) - bar_w / 2
        x_max = max(idx) + bar_w / 2
        ax.hlines(mean_val, x_min, x_max, colors=COLORS[cond],
                  linewidths=1.8, linestyles='--', zorder=4)
        ax.text((x_min + x_max) / 2, mean_val,
                f'{cond} mean {fmt.format(mean_val)}',
                ha='center', va='bottom', fontsize=FS_MEAN, color=COLORS[cond],
                fontweight='bold', zorder=5,
                bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.8))

    # y-range: start at zero only where that keeps the bars readable
    span = y_hi - y_lo
    if y_lo > 0 and span < 0.35 * y_hi:
        ax.set_ylim(y_lo - 0.6 * span, y_hi + 0.75 * span)
    else:
        ax.set_ylim(0, y_hi * 1.22)

    ax.set_xticks(x)
    ax.set_xticklabels(df['label'], fontsize=FS_TICK)
    for lab, cond in zip(ax.get_xticklabels(), df['condition']):
        lab.set_color(COLORS[cond])
        lab.set_fontweight('bold')
    ax.set_ylabel(ylabel, fontsize=FS_LABEL)
    ax.set_title(title, fontsize=FS_TITLE, fontweight='bold', pad=6)
    ax.tick_params(axis='y', labelsize=FS_TICK)
    ax.set_xlim(-0.6, len(df) - 0.4)
    ax.grid(axis='y', linewidth=0.5, color='#dddddd', zorder=1)
    ax.set_axisbelow(True)

for j in range(len(panels), len(axes)):
    axes[j].axis('off')

# ──────────────────────────────────────────────────────────────────
# Legend, titles & layout
# ──────────────────────────────────────────────────────────────────
legend_handles = [
    mpatches.Patch(color=COLORS['FLT'], label='Spaceflight (FLT; F1-F9)'),
    mpatches.Patch(color=COLORS['GC'],  label='Ground control (GC; G1-G9)'),
    Line2D([0], [0], color='#555555', lw=1.8, ls='--', label='group mean (dashed, labeled)'),
]
fig.legend(handles=legend_handles, loc='lower center', ncol=3,
           fontsize=FS_LABEL, frameon=False, bbox_to_anchor=(0.5, -0.012))

fig.suptitle('OSD-564 RNA-seq quality-control metrics, RR-10 FLT (n=5) vs GC (n=5)\n'
             'GC libraries show lower unique-mapping and higher rRNA fraction than FLT',
             fontsize=FS_TITLE + 1, fontweight='bold', y=1.005)

plt.tight_layout(h_pad=2.2, w_pad=2.0, rect=(0, 0.012, 1, 0.985))

os.makedirs('QC_figures', exist_ok=True)
for ext in ('pdf', 'png'):
    out = f'QC_figures/OSD564_QC_summary.{ext}'
    plt.savefig(out, bbox_inches='tight', dpi=300)
    print(f'Saved: {out}')
