"""
Legacy styling helpers for the OSD-564 bar panels (Scientific Reports style).

NOTE: the GSEA NES panels are no longer drawn here.  Every NES panel in the
project (OSD-564, OSD-562, OSD-352) is drawn by the shared srep_nes_panel.py
(nes_panel), which fixes the term order, the axis limits and the layout across
datasets.  This module is kept only for make_gprofiler_ORA_figures_OSD564.py
(label_case, wrap_lines, fit_xlim_to_texts, ...); plot_nes_panel is unused.

Conventions (shared with the OSD-562 figure scripts):
  * every requested term is drawn whether or not it is significant
  * FDR label on every bar via srep.fdr_label ("FDR = 0.012" / "FDR < 0.001"),
    asterisk appended at FDR <= 0.05; non-significant bars drawn in the pale colours
  * red  = Enriched in spaceflight (NES > 0)
  * blue = Enriched in ground control (NES < 0)
  * y-axis term names in sentence case (GO IDs dropped; they are in the tables),
    wrapped at ~38 characters
  * one 8 pt title line: "OSD-564, spaceflight vs ground control (n = 5 vs 5)"
  * panel size srep.PANEL_MM["nes"] (90 x 60 mm) for up to 7 terms; taller panels
    (+7 mm per extra term) for the 10- and 13-term supplementary panels
  * written as 300 dpi PNG + vector PDF by srep.save
"""
import os
import re
import textwrap
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import srep_style as srep; srep.apply()
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator

RED      = srep.COLORS['flight']
BLUE     = srep.COLORS['ground']
RED_NS   = srep.COLORS['nom_up']
BLUE_NS  = srep.COLORS['nom_down']
STUDY_TAG = 'OSD-564, spaceflight vs ground control (n = 5 vs 5)'

_GO_ID = re.compile(r'\s*\((GO|KEGG|REAC|WP|HP):?\d+\)\s*$')


def label_case(label):
    """srep.sentence_case plus: keep single upper-case letters that stand alone
    ("Phospholipase D", "Myosin V")."""
    words = label.split(' ')
    out = srep.sentence_case(label).split(' ')
    for i, (orig, new) in enumerate(zip(words, out)):
        core = orig.strip('(),;:[]')
        if len(core) == 1 and core.isupper():
            out[i] = orig
    return ' '.join(out)


def wrap_lines(text, min_width, max_lines=2):
    """Wrap at the smallest width >= min_width that needs at most `max_lines`
    lines (plain textwrap at a fixed width can split a 75-character label into
    three uneven lines)."""
    kw = dict(break_on_hyphens=False, break_long_words=False)
    for w in range(min_width, max(min_width, len(text)) + 1):
        lines = textwrap.wrap(text, w, **kw)
        if len(lines) <= max_lines:
            return '\n'.join(lines)
    return '\n'.join(textwrap.wrap(text, min_width, **kw))


def term_label(term, width=30, keep_id=False):
    """Sentence-case term name (GO ID stripped unless keep_id) wrapped for the
    y-axis: at most two lines, at least `width` characters per line."""
    name = term.strip() if keep_id else _GO_ID.sub('', term.strip())
    return wrap_lines(label_case(name), width)


def fdr_text(fdr):
    return srep.fdr_label(fdr) + (' *' if fdr <= 0.05 else '')


def bar_color(nes, fdr):
    if nes < 0:
        return BLUE if fdr <= 0.05 else BLUE_NS
    return RED if fdr <= 0.05 else RED_NS


def fit_xlim_to_texts(fig, ax, anchored, x_lo_data, x_hi_data, pad_mm=1.5,
                      bar_pad_frac=0.03, passes=2, lo_fixed=None):
    """Set ax xlim so that every anchored text lies inside the axes.

    anchored : list of (text, x_anchor_data, side) with side 'right' when the
               text extends to the right of its anchor and 'left' otherwise.
               Texts must have in_layout=False so the layout engine does not
               shrink the axes to make room for them.
    x_lo_data, x_hi_data : extremes of the bars that must stay inside.
    lo_fixed : pin the lower limit (e.g. 0 for a -log10 p axis).

    The required range R = hi - lo satisfies
        lo = min(x_lo - b R, min_i(a_i - f_i R)),  hi = max(x_hi + b R, max_j(a_j + f_j R))
    with f = (text width + pad) / axes width; solved by bisection on R.
    """
    pad_px = pad_mm * fig.dpi / 25.4
    for _ in range(passes):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        W = ax.get_window_extent(r).width
        left, right = [], []
        for t, a, side in anchored:
            f = (t.get_window_extent(r).width + pad_px) / W
            (right if side == 'right' else left).append((a, f))

        def lo_of(R):
            if lo_fixed is not None:
                return lo_fixed
            return min([x_lo_data - bar_pad_frac * R] + [a - f * R for a, f in left])

        def hi_of(R):
            return max([x_hi_data + bar_pad_frac * R] + [a + f * R for a, f in right])

        def h(R):
            return hi_of(R) - lo_of(R) - R

        R_a, R_b = 1e-6, 1e4
        if h(R_b) > 0:
            print('  [WARN] FDR labels cannot fit in the axes width; panel too narrow')
            R = R_b
        else:
            for _ in range(60):
                R_m = 0.5 * (R_a + R_b)
                if h(R_m) > 0:
                    R_a = R_m
                else:
                    R_b = R_m
            R = R_b
        ax.set_xlim(lo_of(R), hi_of(R))


def xlabel_no_clip(fig, ax, label):
    """Centred x-axis label; right-aligned to the axes when a centred label would
    run past the panel edge (narrow axes next to a wide y-label column)."""
    ax.set_xlabel(label)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fw = fig.get_window_extent(r).width
    if ax.xaxis.label.get_window_extent(r).x1 > fw - 2:
        ax.set_xlabel(label, loc='right')


def panel_size_mm(n_terms):
    w, h = srep.PANEL_MM['nes']
    if n_terms > 7:
        h = h + 7 * (n_terms - 7)
    return (w, h)


def plot_nes_panel(records, panel_name, out_path, subtitle=None, size_mm=None,
                   wrap=30, title=STUDY_TAG):
    """
    records : list of dicts with keys 'label', 'NES', 'FDR'
    Draws a horizontal NES bar chart sorted by NES (most negative at bottom)
    and writes <out_path without extension>.png and .pdf.
    panel_name / subtitle are only echoed to the console (the figure legend
    carries them); the single in-panel line is `title`.
    """
    df = (pd.DataFrame(records)
          .sort_values('NES', ascending=True)
          .reset_index(drop=True))
    n = len(df)
    if n == 0:
        print(f'  [SKIP] no terms for {panel_name}')
        return None

    w_mm, h_mm = size_mm or panel_size_mm(n)
    fig, ax = plt.subplots(figsize=srep.figsize_mm(w_mm, h_mm), constrained_layout=True)

    # FDR text at the bar tip when all bars point the same way; in a mixed-sign
    # panel the text sits on the empty side of the zero line instead (tip labels
    # on both sides would leave almost no room for the bars at 90 mm width)
    mixed = bool((df['NES'] > 0).any() and (df['NES'] < 0).any())
    anchored = []
    for i, row in df.iterrows():
        nes, fdr = float(row['NES']), float(row['FDR'])
        ax.barh(i, nes, color=bar_color(nes, fdr), edgecolor='none', linewidth=0,
                height=0.68, zorder=2)
        if mixed:
            x_anchor, side = 0.0, ('left' if nes >= 0 else 'right')
        else:
            x_anchor, side = nes, ('right' if nes >= 0 else 'left')
        dx = 3 if side == 'right' else -3
        t = ax.annotate(fdr_text(fdr), xy=(x_anchor, i), xytext=(dx, 0),
                        textcoords='offset points', va='center',
                        ha='left' if side == 'right' else 'right',
                        color='black', zorder=4, annotation_clip=False)
        t.set_in_layout(False)       # xlim is widened for it below instead
        anchored.append((t, x_anchor, side))

    ax.set_yticks(range(n))
    ax.set_yticklabels([term_label(l, wrap) for l in df['label']], linespacing=1.0)
    ax.axvline(0, color='black', lw=srep.LINE_PT, zorder=3)
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)

    x_min = float(min(0.0, df['NES'].min()))
    x_max = float(max(0.0, df['NES'].max()))
    ax.set_xlim(x_min - 0.1, x_max + 0.1)
    ax.set_ylim(-0.6, n - 0.4)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10]))
    if title:
        fig.suptitle(title, x=0.01, ha='left')

    handles = []
    if (df['NES'] > 0).any():
        handles.append(mpatches.Patch(color=RED, label='Enriched in spaceflight (NES > 0)'))
    if (df['NES'] < 0).any():
        handles.append(mpatches.Patch(color=BLUE, label='Enriched in ground control (NES < 0)'))
    fig.legend(handles=handles, loc='outside lower left', ncol=1, frameon=False)

    fit_xlim_to_texts(fig, ax, anchored, x_min, x_max)
    xlabel_no_clip(fig, ax, 'Normalized enrichment score (NES)')
    fit_xlim_to_texts(fig, ax, anchored, x_min, x_max)

    base = os.path.splitext(out_path)[0]
    srep.save(fig, base)
    plt.close(fig)
    print(f'  Saved -> {base}.png / .pdf   ({panel_name}; {w_mm} x {h_mm} mm)')
    return df
