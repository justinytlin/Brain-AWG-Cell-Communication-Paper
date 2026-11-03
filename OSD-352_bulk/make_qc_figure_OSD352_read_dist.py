"""
OSD-352 (RR-3 whole-brain bulk RNA-seq) RSeQC read-distribution QC summary,
FLT (n = 3) vs GC (n = 3).  Fonts / figure size enlarged so that metric names
and sample IDs are legible at manuscript size.
"""
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
# Load & annotate
# ──────────────────────────────────────────────────────────────────
df = pd.read_csv(
    'read_dist_multiqc_report/read_dist_multiqc_data/multiqc_rseqc_read_distribution.txt',
    sep='\t'
)
df.rename(columns={'Sample': 'sample'}, inplace=True)

df['condition'] = df['sample'].apply(lambda s: 'FLT' if '_FLT_' in s else 'GC')
# e.g. RR3_BRN_FLT_F1_... → FLT-F1
df['label'] = df['sample'].apply(
    lambda s: ('FLT-' if '_FLT_' in s else 'GC-') + s.split('_')[-3]
)

# Sort: FLT first, then GC (alphabetically within each group)
df['condition'] = pd.Categorical(df['condition'], categories=['FLT', 'GC'], ordered=True)
df = df.sort_values(['condition', 'label']).reset_index(drop=True)
df['condition'] = df['condition'].astype(str)
print("Samples:", df['label'].tolist())

COLORS = {'FLT': '#E07B39', 'GC': '#3A7DC9'}

# ──────────────────────────────────────────────────────────────────
# Panel definitions: (column, y-label, title)
# ──────────────────────────────────────────────────────────────────
panels = [
    ('cds_exons_tag_pct',        'Tags (%)', 'CDS exons'),
    ('5_utr_exons_tag_pct',      'Tags (%)', "5' UTR exons"),
    ('3_utr_exons_tag_pct',      'Tags (%)', "3' UTR exons"),
    ('introns_tag_pct',          'Tags (%)', 'Introns'),
    ('tss_up_1kb_tag_pct',       'Tags (%)', 'TSS upstream 1 kb'),
    ('tss_up_5kb_tag_pct',       'Tags (%)', 'TSS upstream 5 kb'),
    ('tss_up_10kb_tag_pct',      'Tags (%)', 'TSS upstream 10 kb'),
    ('tes_down_1kb_tag_pct',     'Tags (%)', 'TES downstream 1 kb'),
    ('tes_down_5kb_tag_pct',     'Tags (%)', 'TES downstream 5 kb'),
    ('tes_down_10kb_tag_pct',    'Tags (%)', 'TES downstream 10 kb'),
    ('other_intergenic_tag_pct', 'Tags (%)', 'Other intergenic'),
]

N_COLS = 4
N_ROWS = 3  # 12 slots for 11 panels

FS_TITLE  = 13   # panel titles (metric names)
FS_TICK   = 11   # sample IDs / y ticks
FS_YLABEL = 11.5
FS_LEGEND = 13
FS_SUP    = 17

fig, axes = plt.subplots(N_ROWS, N_COLS, figsize=(22, N_ROWS * 5.2))
axes = axes.flatten()

x = np.arange(len(df))
bar_w = 0.72

for i, (col, ylabel, title) in enumerate(panels):
    ax = axes[i]
    vals = df[col].copy()

    bar_colors = [COLORS[c] for c in df['condition']]
    ax.bar(x, vals, color=bar_colors, edgecolor='white', linewidth=0.5, width=bar_w, zorder=2)

    # Group mean lines
    for cond, grp in df.groupby('condition'):
        idx = grp.index.tolist()
        mean_val = vals.iloc[idx].mean()
        ax.hlines(mean_val, min(idx) - bar_w / 2, max(idx) + bar_w / 2,
                  colors=COLORS[cond], linewidths=2.0, linestyles='--', zorder=3)

    ax.set_xticks(x)
    ax.set_xticklabels(df['label'], rotation=45, ha='right', fontsize=FS_TICK)
    ax.set_ylabel(ylabel, fontsize=FS_YLABEL)
    ax.set_title(title, fontsize=FS_TITLE, fontweight='bold', pad=8)
    ax.tick_params(axis='y', labelsize=FS_TICK)
    ax.set_xlim(-0.6, len(df) - 0.4)
    ax.grid(axis='y', linewidth=0.5, color='#dddddd', zorder=1)
    ax.set_axisbelow(True)

# Hide unused 12th panel
axes[-1].set_visible(False)

# ──────────────────────────────────────────────────────────────────
# Legend & layout
# ──────────────────────────────────────────────────────────────────
legend_patches = [
    mpatches.Patch(color=COLORS['FLT'], label='Spaceflight (FLT, n = 3)'),
    mpatches.Patch(color=COLORS['GC'],  label='Ground control (GC, n = 3)'),
    plt.Line2D([0], [0], color='#555555', lw=2.0, ls='--', label='Group mean'),
]
fig.legend(handles=legend_patches, loc='lower center', ncol=3,
           fontsize=FS_LEGEND, frameon=False, bbox_to_anchor=(0.5, -0.005))

fig.suptitle('OSD-352 bulk RNA-seq read distribution (RSeQC), FLT (n=3) vs GC (n=3)',
             fontsize=FS_SUP, fontweight='bold', y=1.005)

plt.tight_layout(h_pad=3.5, w_pad=2.5, rect=(0, 0.02, 1, 1))

os.makedirs('QC_figures', exist_ok=True)
for fmt in ('pdf', 'png'):
    out = f'QC_figures/OSD352_read_dist_QC_summary.{fmt}'
    plt.savefig(out, bbox_inches='tight', dpi=300)
    print(f'Saved: {out}')
