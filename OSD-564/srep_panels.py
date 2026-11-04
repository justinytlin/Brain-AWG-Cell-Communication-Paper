"""
Shared panel drawers for the three bulk RNA-seq datasets (OSD-564, OSD-562, OSD-352)
so that the PCA, volcano and QC panels are identical in convention and can sit side
by side in one journal figure.  An identical copy of this file sits in every figure
folder next to srep_style.py; the dataset scripts only prepare data and call

    pca_panel(coords, labels, groups, pc_var, title, out_base)        60 x 60 mm
    volcano_panel(df, size_key, title, label_genes, highlight_genes,
                  out_base, counts, ...)                             118 x 100 / 90 x 85 mm
    qc_grid(df, panels, out_base)                                     183 x 120 mm, 4 x 3

Conventions (all lettering Arial 8 pt via srep_style.apply()):

PCA      filled circles, 1 pt black edge, red = spaceflight, blue = ground control;
         short sample IDs in 8 pt grey next to every point; no zero lines; axis labels
         "PC1 (42.0% variance)"; frameless legend "Spaceflight" / "Ground control"
         placed where the axes are empty; one 8 pt title line (dataset only).

Volcano  grey "Not significant"; pale red/blue circles for nominal p <= 0.05 and
         |log2FC| >= 1; red/blue stars for adjusted p <= 0.05 and |log2FC| >= 1;
         orange filled circles with 1 pt black edge for curated genes (optional orange
         open diamonds for an "extra" set); dashed grey 1 pt threshold lines at
         nominal p = 0.05 and |log2FC| = 1; frameless legend below the axes with
         identical wording everywhere (adjusted entries with n = 0 are omitted and
         replaced by a "No gene at adjusted p <= 0.05" text line); italic 8 pt gene
         labels placed by adjustText with 1 pt leader lines and a white backing only
         where the label would sit on data points; symmetric x-limits.

         NOTE: one legend entry ("Nominal p <= 0.05, |log2FC| >= 1, down in spaceflight
         (n = 244)") is ~76 mm wide at 8 pt, so two legend columns do not fit in either
         panel width; the legend is one column, left-aligned, in all three datasets.

QC       4 x 3 grid of bar panels, every y-axis from 0, 1 pt dashed black group-mean
         lines, sample IDs as x tick labels, metric titles in sentence case with
         units, one shared frameless legend below the grid, no title line.
"""
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator
from adjustText import adjust_text
import srep_style as srep

# italic gene symbols inside legend strings ($\it{Vps4a}$) must also be Arial
plt.rcParams.update({
    'mathtext.fontset': 'custom',
    'mathtext.rm': 'Arial',
    'mathtext.it': 'Arial:italic',
    'mathtext.bf': 'Arial:bold',
})

MINUS = '−'                      # the same Unicode minus in every axis label
THR_GREY = '#9e9e9e'
LEADER_GREY = '#555555'
LABEL_GREY = '#555555'

FC_THR, NOM_THR, ADJ_THR = 1.0, 0.05, 0.05
SEED = 0.10                           # initial label offset (axes fraction) before adjustText

XLAB = 'log2 fold change (spaceflight/ground control)'
YLAB = f'{MINUS}log10(nominal p)'

GROUP_NAMES = {'flight': 'Spaceflight', 'ground': 'Ground control'}


# ──────────────────────────────────────────────────────────────────
# helpers
# ──────────────────────────────────────────────────────────────────
def _group_key(g):
    """Map 'Spaceflight' / 'FLT' / 'flight' / 'Ground control' / 'GC' / 'ground' to a COLORS key."""
    s = str(g).lower()
    if s in ('flight', 'flt', 'spaceflight', 'space flight', 'f'):
        return 'flight'
    if s in ('ground', 'gc', 'ground control', 'g'):
        return 'ground'
    raise ValueError(f'unknown group label: {g!r}')


def _dot(color, label, marker='o', ms=4.0, mec='none', mfc=None, mew=0):
    return Line2D([], [], marker=marker, ls='', markersize=ms, color=color,
                  markerfacecolor=color if mfc is None else mfc, markeredgecolor=mec,
                  markeredgewidth=mew, label=label)


def _text_only(label):
    """Legend entry without a symbol (used for the 'No gene at adjusted p' line)."""
    return Line2D([], [], ls='', marker='', label=label)


def _title(fig, title):
    """One 8 pt title line at the top-left of the figure (constrained layout reserves it)."""
    fig.suptitle(title, x=0.01, ha='left')


# ──────────────────────────────────────────────────────────────────
# A. PCA panel (60 x 60 mm)
# ──────────────────────────────────────────────────────────────────
def pca_panel(coords, labels, groups, pc_var, title, out_base):
    """
    coords : (n, 2) array-like or DataFrame with the PC1 / PC2 scores (first two columns)
    labels : short sample IDs (F1 ... G9), one per row of coords
    groups : group per sample ('Spaceflight'/'flight'/'FLT' or 'Ground control'/'ground'/'GC')
    pc_var : (PC1 %, PC2 %) explained variance
    title  : dataset line, e.g. 'OSD-564, spaceflight vs ground control'
    """
    xy = np.asarray(pd.DataFrame(coords).iloc[:, :2], dtype=float)
    labels = list(labels)
    keys = [_group_key(g) for g in groups]
    assert len(labels) == len(keys) == len(xy)

    fig, ax = plt.subplots(figsize=srep.figsize_mm('pca'), constrained_layout=True)
    for key in ('flight', 'ground'):
        idx = [i for i, k in enumerate(keys) if k == key]
        if not idx:
            continue
        ax.scatter(xy[idx, 0], xy[idx, 1], s=28, color=srep.COLORS[key], marker='o',
                   edgecolors='black', linewidths=srep.LINE_PT, zorder=5,
                   label=GROUP_NAMES[key])
    ax.set_xlabel(f'PC1 ({pc_var[0]:.1f}% variance)')
    ax.set_ylabel(f'PC2 ({pc_var[1]:.1f}% variance)')
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    # head-room for the sample IDs (upper right of each point) and the legend
    x0, x1 = xy[:, 0].min(), xy[:, 0].max()
    y0, y1 = xy[:, 1].min(), xy[:, 1].max()
    dx, dy = (x1 - x0) or 1.0, (y1 - y0) or 1.0
    ax.set_xlim(x0 - 0.18 * dx, x1 + 0.30 * dx)
    ax.set_ylim(y0 - 0.18 * dy, y1 + 0.30 * dy)
    _title(fig, title)
    ax.legend(loc='best', frameon=False, handletextpad=0.4, labelspacing=0.3)

    # sample IDs (8 pt grey, no leader lines): each ID is placed next to its point at the
    # first of eight candidate offsets (upper right first) whose box neither covers a
    # point nor an already placed ID nor the legend, nor leaves the axes
    _place_pca_ids(fig, ax, xy, labels)
    srep.save(fig, out_base)
    plt.close(fig)
    print(f'Saved -> {out_base}.png / .pdf')


def _place_pca_ids(fig, ax, xy, labels, gap_pt=2.5, marker_r_pt=3.5):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    ax_bb = ax.get_window_extent(r)
    leg = ax.get_legend()
    leg_bb = leg.get_window_extent(r) if leg is not None else None
    px = fig.dpi / 72.0
    pts = ax.transData.transform(xy)
    mr = marker_r_pt * px                       # marker radius in pixels (+ gap)
    gap = gap_pt * px
    # (dx, dy, ha, va) in pixels: upper right, upper left, lower right, lower left, right, left, top, bottom
    cands = [(mr + gap, mr, 'left', 'bottom'), (-(mr + gap), mr, 'right', 'bottom'),
             (mr + gap, -mr, 'left', 'top'), (-(mr + gap), -mr, 'right', 'top'),
             (mr + gap, 0, 'left', 'center'), (-(mr + gap), 0, 'right', 'center'),
             (0, mr + gap, 'center', 'bottom'), (0, -(mr + gap), 'center', 'top')]
    placed = []
    for (cx, cy), lab in zip(pts, labels):
        best, best_score = None, None
        for dx, dy, ha, va in cands:
            t = ax.text(0, 0, lab, color=LABEL_GREY, ha=ha, va=va, zorder=6)
            t.set_in_layout(False)
            t.set_position(ax.transData.inverted().transform((cx + dx, cy + dy)))
            bb = t.get_window_extent(r)
            score = 0.0
            if bb.x0 < ax_bb.x0 or bb.x1 > ax_bb.x1 or bb.y0 < ax_bb.y0 or bb.y1 > ax_bb.y1:
                score += 100.0
            for (qx, qy) in pts:                   # any marker inside the (grown) box
                if bb.x0 - mr < qx < bb.x1 + mr and bb.y0 - mr < qy < bb.y1 + mr:
                    score += 10.0
            for ob in placed:
                if bb.overlaps(ob):
                    score += 10.0
            if leg_bb is not None and bb.overlaps(leg_bb):
                score += 10.0
            if best is None or score < best_score:
                if best is not None:
                    best.remove()
                best, best_score, best_bb = t, score, bb
            else:
                t.remove()
            if score == 0:
                break
        placed.append(best_bb)


# ──────────────────────────────────────────────────────────────────
# B. Volcano panel (118 x 100 mm or 90 x 85 mm)
# ──────────────────────────────────────────────────────────────────
def volcano_masks(df):
    """Significance classes over all tested genes (columns log2fc, pval, adjpval)."""
    adj_up   = (df['log2fc'] >=  FC_THR) & (df['adjpval'] <= ADJ_THR)
    adj_down = (df['log2fc'] <= -FC_THR) & (df['adjpval'] <= ADJ_THR)
    nom_up   = (df['log2fc'] >=  FC_THR) & (df['pval'] <= NOM_THR) & ~adj_up
    nom_down = (df['log2fc'] <= -FC_THR) & (df['pval'] <= NOM_THR) & ~adj_down
    not_sig  = ~(adj_up | adj_down | nom_up | nom_down)
    return adj_up, adj_down, nom_up, nom_down, not_sig


def volcano_counts(df):
    adj_up, adj_down, nom_up, nom_down, _ = volcano_masks(df)
    return {'n_nom_up': int(nom_up.sum()), 'n_nom_down': int(nom_down.sum()),
            'n_adj_up': int(adj_up.sum()), 'n_adj_down': int(adj_down.sum())}


def legend_labels(counts):
    f = srep.fmt_int
    return {
        'ns':       'Not significant',
        'nom_up':   f'Nominal p ≤ 0.05, |log2FC| ≥ 1, up in spaceflight (n = {f(counts["n_nom_up"])})',
        'nom_down': f'Nominal p ≤ 0.05, |log2FC| ≥ 1, down in spaceflight (n = {f(counts["n_nom_down"])})',
        'adj_up':   f'Adjusted p ≤ 0.05, |log2FC| ≥ 1, up in spaceflight (n = {f(counts["n_adj_up"])})',
        'adj_down': f'Adjusted p ≤ 0.05, |log2FC| ≥ 1, down in spaceflight (n = {f(counts["n_adj_down"])})',
        'none':     'No gene at adjusted p ≤ 0.05',
        'none_up':  'No gene at adjusted p ≤ 0.05, |log2FC| ≥ 1, up in spaceflight',
        'none_down': 'No gene at adjusted p ≤ 0.05, |log2FC| ≥ 1, down in spaceflight',
        'thr':      'Thresholds (p = 0.05, |log2FC| = 1)',
    }


def _ylog(df):
    return -np.log10(df['pval'].clip(lower=1e-300))


def _draw_layers(ax, df, counts):
    adj_up, adj_down, nom_up, nom_down, not_sig = volcano_masks(df)
    y = _ylog(df)
    L = legend_labels(counts)
    ax.scatter(df.loc[not_sig, 'log2fc'], y[not_sig], s=3, color=srep.COLORS['ns'],
               alpha=0.5, linewidths=0, rasterized=True)
    ax.scatter(df.loc[nom_up, 'log2fc'], y[nom_up], s=7, color=srep.COLORS['nom_up'],
               linewidths=0, rasterized=True)
    ax.scatter(df.loc[nom_down, 'log2fc'], y[nom_down], s=7, color=srep.COLORS['nom_down'],
               linewidths=0, rasterized=True)
    ax.scatter(df.loc[adj_up, 'log2fc'], y[adj_up], s=34, color=srep.COLORS['up'],
               linewidths=0, marker='*', zorder=5)
    ax.scatter(df.loc[adj_down, 'log2fc'], y[adj_down], s=34, color=srep.COLORS['down'],
               linewidths=0, marker='*', zorder=5)
    ax.axhline(-np.log10(NOM_THR), color=THR_GREY, lw=srep.LINE_PT, ls='--', zorder=1)
    ax.axvline( FC_THR, color=THR_GREY, lw=srep.LINE_PT, ls='--', zorder=1)
    ax.axvline(-FC_THR, color=THR_GREY, lw=srep.LINE_PT, ls='--', zorder=1)
    ax.set_xlabel(XLAB)
    ax.set_ylabel(YLAB)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=8, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))

    handles = [_dot(srep.COLORS['ns'], L['ns']),
               _dot(srep.COLORS['nom_up'], L['nom_up']),
               _dot(srep.COLORS['nom_down'], L['nom_down'])]
    if counts['n_adj_up']:
        handles.append(_dot(srep.COLORS['up'], L['adj_up'], marker='*', ms=7))
    if counts['n_adj_down']:
        handles.append(_dot(srep.COLORS['down'], L['adj_down'], marker='*', ms=7))
    if not counts['n_adj_up'] and not counts['n_adj_down']:
        handles.append(_text_only(L['none']))
    elif not counts['n_adj_up']:
        handles.append(_text_only(L['none_up']))
    elif not counts['n_adj_down']:
        handles.append(_text_only(L['none_down']))
    handles.append(Line2D([], [], color=THR_GREY, lw=srep.LINE_PT, ls='--', label=L['thr']))
    return handles


def _limits(ax, df, pad_x, pad_top):
    """Symmetric x-limits (quantile based) and a y-range that always includes the top gene."""
    y = _ylog(df)
    half = max(abs(df['log2fc'].quantile(0.001)), abs(df['log2fc'].quantile(0.999))) + pad_x
    ax.set_xlim(-half, half)
    ax.set_ylim(-0.2, max(y.quantile(0.9998) * 1.15, y.max() * 1.08) * pad_top)


def _separate_labels(fig, ax, texts, pad_pt=1.0, iters=80):
    """Deterministic clean-up after adjustText: while two label boxes overlap, move the
    upper one up and the lower one down by half the overlap plus a pad, keeping every
    label inside the axes.  Works in display pixels."""
    pad = pad_pt * fig.dpi / 72.0
    for _ in range(iters):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        ax_bb = ax.get_window_extent(r)
        bbs = [t.get_window_extent(r).expanded(1.02, 1.02) for t in texts]
        moved = False
        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                a, b = bbs[i], bbs[j]
                if a.overlaps(b):
                    dy = min(a.y1, b.y1) - max(a.y0, b.y0) + pad
                    up, dn = (i, j) if a.y0 >= b.y0 else (j, i)
                    for k, sgn in ((up, +1), (dn, -1)):
                        x, yy = ax.transData.transform(texts[k].get_position())
                        if bbs[k].y1 + sgn * dy / 2 > ax_bb.y1 or bbs[k].y0 + sgn * dy / 2 < ax_bb.y0:
                            continue
                        texts[k].set_position(ax.transData.inverted().transform((x, yy + sgn * dy / 2)))
                        moved = True
                    break
            if moved:
                break
        if not moved:
            return


def _push_inside(fig, ax, texts, pad_pt=2.0):
    """Shift every label so that its box lies inside the axes (plus a small pad)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    ax_bb = ax.get_window_extent(r)
    pad = pad_pt * fig.dpi / 72.0
    for t in texts:
        bb = t.get_window_extent(r)
        dx = dy = 0.0
        if bb.x0 < ax_bb.x0 + pad:
            dx = ax_bb.x0 + pad - bb.x0
        elif bb.x1 > ax_bb.x1 - pad:
            dx = ax_bb.x1 - pad - bb.x1
        if bb.y0 < ax_bb.y0 + pad:
            dy = ax_bb.y0 + pad - bb.y0
        elif bb.y1 > ax_bb.y1 - pad:
            dy = ax_bb.y1 - pad - bb.y1
        if dx or dy:
            x, y = ax.transData.transform(t.get_position())
            t.set_position(ax.transData.inverted().transform((x + dx, y + dy)))


def _place_labels(fig, ax, df, label_df, avoid_df):
    """Italic labels placed by adjustText, 1 pt leader lines, white backing only where
    the label would sit on data points."""
    if label_df is None or len(label_df) == 0:
        return
    ly = _ylog(label_df)
    # seed every label a short distance from its point, radially away from the centroid
    # of the labelled points (in axes fractions) so that the leader lines fan outwards
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    frac = np.column_stack([(label_df['log2fc'].values - x0) / (x1 - x0), (ly.values - y0) / (y1 - y0)])
    cen = frac.mean(axis=0)
    texts, anchors = [], []
    for (px, py), (fx, fy), name in zip(zip(label_df['log2fc'], ly), frac, label_df['display']):
        v = np.array([fx, fy]) - cen
        v = v / np.linalg.norm(v) if np.linalg.norm(v) > 1e-6 else np.array([0.0, 1.0])
        v[1] = max(v[1], 0.2)                       # never seed straight below a point
        v = v / np.linalg.norm(v)
        sx, sy = x0 + (fx + SEED * v[0]) * (x1 - x0), y0 + (fy + SEED * v[1]) * (y1 - y0)
        txt = ax.text(sx, sy, name, ha='center', va='center', fontstyle='italic',
                      color='black', zorder=8, clip_on=False)
        txt.set_in_layout(False)
        texts.append(txt)
        anchors.append((px, py))
    fig.canvas.draw()                      # final axes geometry before the repulsion
    # static points the labels are repelled from: the labelled / highlighted genes
    adjust_text(texts, ax=ax, x=avoid_df['log2fc'].values, y=_ylog(avoid_df).values,
                force_text=(0.5, 0.9), force_static=(0.3, 0.6),
                force_pull=(0.01, 0.02), force_explode=(0.3, 0.6),
                expand=(1.3, 1.6), max_move=(300, 300), iter_lim=2000,
                ensure_inside_axes=True, min_arrow_len=99999)
    _separate_labels(fig, ax, texts)
    _push_inside(fig, ax, texts)

    # leader lines from each point to its (moved) label
    for txt, (px, py) in zip(texts, anchors):
        lx, ly_ = txt.get_position()
        if np.hypot(lx - px, ly_ - py) > 0.08:
            ax.annotate('', xy=(px, py), xytext=(lx, ly_), xycoords='data', textcoords='data',
                        zorder=6, arrowprops=dict(arrowstyle='-', color=LEADER_GREY,
                                                  lw=srep.LINE_PT, shrinkA=2, shrinkB=2))

    # white backing only where the label box covers data points (box grown by ~3 pt so
    # that the larger star / circle markers count as covered)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    pts = ax.transData.transform(np.column_stack([df['log2fc'].values, _ylog(df).values]))
    grow = 3.0 * fig.dpi / 72.0
    for txt in texts:
        bb = txt.get_window_extent(r)
        bb = matplotlib.transforms.Bbox([[bb.x0 - grow, bb.y0 - grow], [bb.x1 + grow, bb.y1 + grow]])
        inside = ((pts[:, 0] >= bb.x0) & (pts[:, 0] <= bb.x1) &
                  (pts[:, 1] >= bb.y0) & (pts[:, 1] <= bb.y1))
        if inside.any():
            txt.set_bbox(dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.85))


def volcano_panel(df, size_key, title, label_genes, highlight_genes, out_base, counts=None,
                  extra_genes=None, extra_label=None, highlight_label='Curated genes'):
    """
    df              : ALL tested genes; columns log2fc, pval, adjpval (legend counts use this)
    size_key        : 'volcano' (118 x 100 mm) or 'volcano_supp' (90 x 85 mm)
    title           : one 8 pt title line (dataset, contrast, n)
    label_genes     : DataFrame (log2fc, pval, display) of genes to label, or None
    highlight_genes : DataFrame (log2fc, pval) drawn as orange circles, 1 pt black edge, or None
    counts          : dict n_nom_up / n_nom_down / n_adj_up / n_adj_down (default: from df)
    extra_genes     : optional DataFrame drawn as orange open diamonds (legend text extra_label)
    """
    counts = counts or volcano_counts(df)
    fig, ax = plt.subplots(figsize=srep.figsize_mm(size_key), constrained_layout=True)
    handles = _draw_layers(ax, df, counts)

    cur = srep.COLORS['curated']
    has_labels = label_genes is not None and len(label_genes) > 0
    if highlight_genes is not None and len(highlight_genes):
        ax.scatter(highlight_genes['log2fc'], _ylog(highlight_genes), s=16, color=cur,
                   linewidths=srep.LINE_PT, edgecolors='black', zorder=6)
        handles.append(_dot(cur, f'{highlight_label} (n = {srep.fmt_int(len(highlight_genes))})',
                            ms=4.5, mec='black', mew=srep.LINE_PT))
    if extra_genes is not None and len(extra_genes):
        ax.scatter(extra_genes['log2fc'], _ylog(extra_genes), s=24, color='white',
                   linewidths=srep.LINE_PT, marker='D', edgecolors=cur, zorder=7)
        handles.append(_dot(cur, extra_label, marker='D', ms=5, mfc='white', mec=cur,
                            mew=srep.LINE_PT))

    # labelled panels get wider flanks and more head-room for the labels
    _limits(ax, df, pad_x=1.6 if has_labels else 0.4, pad_top=1.2 if has_labels else 1.0)
    _title(fig, title)
    fig.legend(handles=handles, loc='outside lower left', ncol=1, frameon=False,
               handletextpad=0.5, labelspacing=0.25, columnspacing=1.0)

    if has_labels:
        # labels are repelled from the labelled / highlighted genes and from the stars
        adj_up, adj_down, *_ = volcano_masks(df)
        avoid = [label_genes, df.loc[adj_up | adj_down]]
        if highlight_genes is not None and len(highlight_genes):
            avoid.append(highlight_genes)
        if extra_genes is not None and len(extra_genes):
            avoid.append(extra_genes)
        avoid_df = pd.concat([a[['log2fc', 'pval']] for a in avoid]).drop_duplicates()
        _place_labels(fig, ax, df, label_genes, avoid_df)

    srep.save(fig, out_base)
    plt.close(fig)
    print(f'Saved -> {out_base}.png / .pdf')


# ──────────────────────────────────────────────────────────────────
# C. QC grid (183 x 120 mm, 4 x 3)
# ──────────────────────────────────────────────────────────────────
def qc_grid(df, panels, out_base, n_cols=4, n_rows=3):
    """
    df     : one row per library, columns 'label' (short sample ID), 'condition'
             ('FLT'/'GC' or 'flight'/'ground'), plus the metric columns; rows in plot order
    panels : list of (column, title with unit, scale_fn or None, integer_axis)
    """
    keys = [_group_key(c) for c in df['condition']]
    n_f, n_g = keys.count('flight'), keys.count('ground')
    fig, axes = plt.subplots(n_rows, n_cols, figsize=srep.figsize_mm('qc'), constrained_layout=True)
    axes = axes.flatten()
    x = np.arange(len(df))
    bar_w = 0.72
    for i, (col, title, scale_fn, integer) in enumerate(panels):
        ax = axes[i]
        vals = df[col].astype(float).copy()
        if scale_fn:
            vals = vals.apply(scale_fn)
        ax.bar(x, vals, color=[srep.COLORS[k] for k in keys], edgecolor='none', linewidth=0,
               width=bar_w, zorder=2)
        for key in ('flight', 'ground'):            # 1 pt dashed group-mean lines
            idx = [j for j, k in enumerate(keys) if k == key]
            if idx:
                ax.hlines(vals.iloc[idx].mean(), min(idx) - bar_w / 2, max(idx) + bar_w / 2,
                          colors='black', linewidths=srep.LINE_PT, linestyles='--', zorder=4)
        ax.set_ylim(0, vals.max() * 1.12)           # every bar axis starts at zero
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))
        if integer:
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: srep.fmt_int(v)))
        ax.set_xticks(x)
        ax.set_xticklabels(df['label'], rotation=90)
        ax.set_xlim(-0.6, len(df) - 0.4)
        ax.tick_params(axis='x', length=0)
        ax.set_title(srep.sentence_case(title))
    for j in range(len(panels), len(axes)):
        axes[j].axis('off')
    handles = [
        mpatches.Patch(color=srep.COLORS['flight'], label=f'Spaceflight (n = {n_f})'),
        mpatches.Patch(color=srep.COLORS['ground'], label=f'Ground control (n = {n_g})'),
        Line2D([0], [0], color='black', lw=srep.LINE_PT, ls='--', label='Group mean'),
    ]
    fig.legend(handles=handles, loc='outside lower center', ncol=3, frameon=False,
               handlelength=1.2, handletextpad=0.4, columnspacing=1.5)
    srep.save(fig, out_base)
    plt.close(fig)
    print(f'Saved -> {out_base}.png / .pdf')
