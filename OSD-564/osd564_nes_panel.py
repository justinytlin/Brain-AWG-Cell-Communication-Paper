"""
Shared styling for the OSD-564 GSEA NES bar panels.

Used by: make_summary_OSD564.py, make_pathway_figures_OSD564.py,
         make_GO_mito_NES_figure.py, make_GO_MF_themed_figures.py

Harmonized conventions (shared with the OSD-562 figure scripts):
  * every requested term is drawn whether or not it is significant
  * FDR label on every bar: "FDR = 0.012", or "FDR < 0.001" when FDR rounds to 0
  * asterisk appended at FDR <= 0.05; non-significant bars drawn lighter
  * red  = Enriched in spaceflight (NES > 0)
  * blue = Enriched in ground control (NES < 0)
  * title = "<panel>, OSD-564: spaceflight vs ground control (n = 5 vs 5)"
"""
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors

plt.rcParams.update({
    'font.family': 'Arial',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300,
})

RED  = '#C62828'
BLUE = '#1565C0'
STUDY_TAG = 'OSD-564: spaceflight vs ground control (n = 5 vs 5)'


def panel_title(panel_name):
    return f'{panel_name},\n{STUDY_TAG}'


def fdr_label(fdr):
    """'FDR < 0.001' when the permutation FDR is 0 (or rounds to 0.000)."""
    if fdr < 0.0005:
        return 'FDR < 0.001'
    return f'FDR = {fdr:.3f}'


def bar_rgba(nes, fdr):
    base  = BLUE if nes < 0 else RED
    alpha = 1.00 if fdr <= 0.05 else 0.45
    r, g, b = mcolors.to_rgb(base)
    return (r, g, b, alpha)


def plot_nes_panel(records, panel_name, out_path, subtitle=None,
                   fig_w=10.0, label_fontsize=11):
    """
    records : list of dicts with keys 'label', 'NES', 'FDR'
    Draws a horizontal NES bar chart sorted by NES (most negative at bottom).
    """
    df = (pd.DataFrame(records)
          .sort_values('NES', ascending=True)
          .reset_index(drop=True))
    n = len(df)
    if n == 0:
        print(f'  [SKIP] no terms for {panel_name}')
        return None

    fig_h = max(3.4, n * 0.58 + 2.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for i, row in df.iterrows():
        nes, fdr = float(row['NES']), float(row['FDR'])
        ax.barh(i, nes, color=bar_rgba(nes, fdr), edgecolor='white',
                height=0.68, zorder=2)
        txt = fdr_label(fdr) + (' *' if fdr <= 0.05 else '')
        pad = 0.06 if nes >= 0 else -0.06
        ax.text(nes + pad, i, txt, va='center',
                ha='left' if nes >= 0 else 'right',
                fontsize=9.5, color='#222222', zorder=4)

    ax.set_yticks(range(n))
    ax.set_yticklabels(df['label'], fontsize=label_fontsize)
    ax.axvline(0, color='#222222', lw=1.1, zorder=3)
    ax.xaxis.grid(True, color='#ebebeb', lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)

    x_min = min(0.0, df['NES'].min())
    x_max = max(0.0, df['NES'].max())
    # room for the FDR labels outside the bar tips
    ax.set_xlim(x_min - (1.05 if x_min < 0 else 0.35),
                x_max + (1.05 if x_max > 0 else 0.35))
    ax.set_ylim(-0.7, n - 0.3)

    ax.set_xlabel('Normalized enrichment score (NES)', fontsize=11)
    title = panel_title(panel_name)
    if subtitle:
        title += f'\n{subtitle}'
    ax.set_title(title, fontsize=11.5, fontweight='bold', pad=12)

    handles = [
        mpatches.Patch(color=RED,  alpha=1.00, label='Enriched in spaceflight (NES > 0)'),
        mpatches.Patch(color=BLUE, alpha=1.00, label='Enriched in ground control (NES < 0)'),
        mpatches.Patch(color='grey', alpha=0.45,
                       label='lighter bar: FDR > 0.05;  * = FDR ≤ 0.05'),
    ]
    fig.legend(handles=handles, fontsize=9, frameon=False, ncol=1,
               loc='lower center', bbox_to_anchor=(0.5, 0.0))

    # reserve ~0.85 in at the bottom for the figure-level legend
    plt.tight_layout(rect=(0, 0.85 / fig_h, 1, 1))
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  Saved -> {out_path}')
    return df
