import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'Arial',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300,
})

# ──────────────────────────────────────────────────────────────────
# Load & annotate  — keep only 12-week ISS-T (young) samples
# ──────────────────────────────────────────────────────────────────
df = pd.read_csv('GLDS-557_rna_seq_qc_metrics_GLbulkRNAseq.csv')

df = df[df['sample'].str.contains('ISS-T_YNG', na=False)].copy()
df = df.reset_index(drop=True)

df['condition'] = df['sample'].apply(lambda s: 'FLT' if '_FLT_' in s else 'GC')
df['label'] = df['sample'].apply(
    lambda s: ('FLT-' if '_FLT_' in s else 'GC-') + s.split('_')[-1]
)

COLORS = {'FLT': '#E07B39', 'GC': '#3A7DC9'}

# ──────────────────────────────────────────────────────────────────
# Panel definitions: (column, y-label, title, optional scale fn)
# ──────────────────────────────────────────────────────────────────
panels = [
    ('read_depth',               'Read Depth (M)',         'Read Depth',                   lambda x: x / 1e6),
    ('rin',                      'RIN Score',              'RNA Integrity Number (RIN)',    None),
    ('rrna_contamination',       'rRNA (%)',               'rRNA Contamination',            None),
    ('raw_percent_duplicates_f', 'Duplicates (%)',         '% Duplicate Reads (Raw)',       None),
    ('gene_detected_gt10',       'Genes Detected',         'Genes Detected (>10 counts)',   None),
    ('uniquely_mapped_percent',  'Uniquely Mapped (%)',    'Unique Mapping Rate (STAR)',     None),
    ('multimapped_percent',      'Multi-mapped (%)',       'Multi-mapping Rate (STAR)',      None),
    ('unmapped_tooshort_percent','Unmapped (%)',           'Unmapped – Too Short (STAR)',    None),
    ('mean_genebody_cov_5_20',   "5' Coverage (%)",        "Gene Body Coverage – 5' End",   None),
    ('mean_genebody_cov_80_95',  "3' Coverage (%)",        "Gene Body Coverage – 3' End",   None),
    ('ratio_genebody_cov_3_to_5',"3'/5' Ratio",            "Gene Body 3'/5' Bias",          None),
    ('pct_sense',                'Sense Reads (%)',        'Strandedness (% Sense)',        None),
]

N_COLS = 4
N_ROWS = len(panels) // N_COLS  # exactly 3 rows for 12 panels

fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(16, N_ROWS * 3.5))
axes = axes.flatten()

x = np.arange(len(df))
bar_w = 0.72

for i, (col, ylabel, title, scale_fn) in enumerate(panels):
    ax = axes[i]
    vals = df[col].copy()
    if scale_fn:
        vals = vals.apply(scale_fn)

    bar_colors = [COLORS[c] for c in df['condition']]
    ax.bar(x, vals, color=bar_colors, edgecolor='white', linewidth=0.4, width=bar_w, zorder=2)

    # Group mean lines
    for cond, grp in df.groupby('condition'):
        idx = grp.index.tolist()
        mean_val = vals.iloc[idx].mean()
        x_min = min(idx) - bar_w / 2
        x_max = max(idx) + bar_w / 2
        ax.hlines(mean_val, x_min, x_max,
                  colors=COLORS[cond], linewidths=1.5, linestyles='--', zorder=3)

    ax.set_xticks(x)
    ax.set_xticklabels(df['label'], rotation=45, ha='right', fontsize=7)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_title(title, fontsize=8.5, fontweight='bold', pad=4)
    ax.tick_params(axis='both', labelsize=7)
    ax.set_xlim(-0.6, len(df) - 0.4)
    ax.grid(axis='y', linewidth=0.4, color='#dddddd', zorder=1)
    ax.set_axisbelow(True)

# ──────────────────────────────────────────────────────────────────
# Legend & layout
# ──────────────────────────────────────────────────────────────────
legend_patches = [
    mpatches.Patch(color=COLORS['FLT'], label='Spaceflight (FLT)'),
    mpatches.Patch(color=COLORS['GC'],  label='Ground Control (GC)'),
]
fig.legend(handles=legend_patches, loc='lower center', ncol=2,
           fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.01))

fig.suptitle('OSD-562 RNA-seq Quality Control Metrics\n12-Week Spaceflight vs Ground Control (ISS-T, Young)',
             fontsize=13, fontweight='bold', y=1.01)

plt.tight_layout(h_pad=3.5, w_pad=2.0)

os.makedirs('QC_figures', exist_ok=True)
for fmt in ('pdf', 'png'):
    out = f'QC_figures/OSD562_QC_summary.{fmt}'
    plt.savefig(out, bbox_inches='tight', dpi=300)
    print(f'Saved: {out}')
