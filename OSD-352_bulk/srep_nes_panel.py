"""
srep_nes_panel.py -- the ONE GSEA NES bar panel used by every dataset folder
(OSD-564, OSD-562, OSD-352).  An identical copy sits in each figure folder so
that the scripts can be run from their own directory; keep the copies in sync.

    from srep_nes_panel import nes_panel
    nes_panel(terms, title, out_base, xlim=None)

terms    : ORDERED list of dicts {label, nes, fdr}; the first term is drawn at
           the top of the panel.  The order is never changed here (the
           harmonized term order is fixed by the calling script so that panels
           from different datasets can sit side by side).
title    : one 8 pt line, left-aligned, sentence case, e.g.
           "GO BP mitochondrial terms, OSD-564 (n = 5 vs 5)"
out_base : output path without extension; <out_base>.png (300 dpi) and
           <out_base>.pdf (vector, Type 42 fonts) are written by srep.save.
           A PNG/PDF that is about to be overwritten is first copied, once, to
           the folder's "_pre_unify/" subfolder (same relative path).
xlim     : symmetric x-limits shared by all panels of the same kind:
           XLIM["mito"] / XLIM["metabolic"] = (-3.2, 3.2),
           XLIM["sEV"]  / XLIM["TNT"]       = (-2.2, 2.2);  None = (-3.2, 3.2).

Panel rules (Scientific Reports style, srep_style.py: Arial 8 pt, 1 pt lines,
white background):
  * 90 x 60 mm (srep.PANEL_MM["nes"]); +6 mm per term beyond seven
  * the axes rectangle is identical in every panel (left edge at 48 mm so the
    wrapped term names fit; right edge 1.5 mm from the panel edge)
  * term names: srep.sentence_case (single capital tokens such as
    "Phospholipase D" and "NADH" kept, see FIXUPS), set IDs such as
    "(GO:0006119)" dropped, wrapped at 34 characters to at most two lines
  * bars: srep.COLORS["up"] for NES > 0, srep.COLORS["down"] for NES < 0;
    full colour at FDR <= 0.05, 40 % alpha above; 1 pt black line at NES = 0
  * FDR text: srep.fdr_label(fdr) + " *" at FDR <= 0.05, just beyond the bar
    tip on the bar's own side; if it would run past the x-limit it is placed
    on the other side of the zero line instead (so nothing is ever clipped)
  * x-axis: "Normalized enrichment score (NES)", ticks every 1; no y tick
    marks; left spine removed (the zero line serves); bottom spine 1 pt
  * legend below the axes, no frame, 8 pt, two rows (one row does not fit at
    90 mm): red patch "Enriched in spaceflight (NES > 0)", blue patch
    "Enriched in ground control (NES < 0)", and the note
    "* FDR <= 0.05 (pale bars: FDR > 0.05)"
"""
import os
import re
import math
import shutil
import textwrap
import matplotlib
matplotlib.use('Agg')
import srep_style as srep
srep.apply()
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import MultipleLocator

# ── shared constants ─────────────────────────────────────────────────────────
XLIM_DEFAULT = (-3.2, 3.2)
XLIM = {
    'mito':      (-3.2, 3.2),
    'metabolic': (-3.2, 3.2),
    'sEV':       (-2.2, 2.2),
    'TNT':       (-2.2, 2.2),
}
WRAP_CHARS   = 34      # wrap width for term names
MAX_LINES    = 2
BASE_TERMS   = 7       # up to this many terms: 90 x 60 mm
EXTRA_MM     = 6.0     # extra panel height per term beyond BASE_TERMS
AXES_LEFT_MM = 48.0    # axes left edge (term-name column to the left of it)
RIGHT_MM     = 1.5     # panel edge -> axes right edge
EDGE_MM      = 1.0     # title / legend inset from the panel edge
GAP_MM       = 1.0     # gap between title and axes, and between x-axis block and legend
ALPHA_NS     = 0.4     # bar alpha at FDR > 0.05
BAR_HEIGHT   = 0.68
FDR_PAD_PT   = 2.0     # gap between bar tip (or zero line) and the FDR text
SIG          = 0.05

COL_UP   = srep.COLORS['up']
COL_DOWN = srep.COLORS['down']
XLABEL   = 'Normalized enrichment score (NES)'
LEGEND_UP   = 'Enriched in spaceflight (NES > 0)'
LEGEND_DOWN = 'Enriched in ground control (NES < 0)'
LEGEND_NOTE = ('* FDR ≤ 0.05', '(pale bars: FDR > 0.05)')
BACKUP_DIR  = '_pre_unify'

# srep.sentence_case lower-cases stand-alone capital letters ("Phospholipase D");
# single capitals are restored generically in label_case, and this map catches
# any remaining known case (keys are the lower-cased label).
FIXUPS = {
    'phospholipase d signaling pathway': 'Phospholipase D signaling pathway',
    'mitochondrial electron transport, nadh to ubiquinone':
        'Mitochondrial electron transport, NADH to ubiquinone',
}
_SET_ID = re.compile(r'\s*\((GO|KEGG|REAC|WP|HP):?\d+\)\s*$')


# ── text helpers ─────────────────────────────────────────────────────────────
def label_case(name):
    """srep.sentence_case, keeping single upper-case letters that stand alone in
    the source ("Phospholipase D", "Myosin V") and applying FIXUPS."""
    name = name.strip()
    src = name.split(' ')
    out = srep.sentence_case(name).split(' ')
    for i, (a, b) in enumerate(zip(src, out)):
        core = a.strip('(),;:[]')
        if len(core) == 1 and core.isupper():
            out[i] = a
    s = ' '.join(out)
    return FIXUPS.get(s.lower(), s)


def wrap_label(text, width=WRAP_CHARS, max_lines=MAX_LINES):
    """Wrap at `width` characters; if that needs more than `max_lines` lines,
    widen until it does not (never more than two lines for a y label)."""
    kw = dict(break_long_words=False, break_on_hyphens=False)
    w = width
    while True:
        lines = textwrap.wrap(text, w, **kw)
        if len(lines) <= max_lines or w >= len(text):
            return '\n'.join(lines)
        w += 1


def term_label(name):
    """Y-axis label for a term: set ID dropped, sentence case, wrapped."""
    return wrap_label(label_case(_SET_ID.sub('', str(name))))


def fdr_text(fdr):
    if fdr is None or (isinstance(fdr, float) and math.isnan(fdr)):
        return 'not tested'
    fdr = float(fdr)
    return srep.fdr_label(fdr) + (' *' if fdr <= SIG else '')


def panel_size_mm(n_terms):
    w, h = srep.PANEL_MM['nes']
    if n_terms > BASE_TERMS:
        h = h + EXTRA_MM * (n_terms - BASE_TERMS)
    return (w, h)


def _get(d, *keys):
    for k in keys:
        if k in d:
            return d[k]
    raise KeyError(keys[0])


def _backup(out_base):
    """Copy an existing <out_base>.png/.pdf to _pre_unify/ (once; never overwritten)."""
    for ext in ('.png', '.pdf'):
        src = os.path.abspath(out_base + ext)
        if not os.path.exists(src):
            continue
        rel = os.path.relpath(src, os.getcwd())
        if rel.startswith('..'):
            dst = os.path.join(os.path.dirname(src), BACKUP_DIR, os.path.basename(src))
        else:
            dst = os.path.join(os.getcwd(), BACKUP_DIR, rel)
        if os.path.exists(dst):
            continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)
        print(f'  backed up {rel} -> {os.path.relpath(dst, os.getcwd())}')


# ── the panel ────────────────────────────────────────────────────────────────
def nes_panel(terms, title, out_base, xlim=None):
    """Draw the harmonized NES bar panel and write <out_base>.png / .pdf.
    Returns the (width, height) in mm of the panel written, or None if empty."""
    recs = []
    for t in terms:
        nes = _get(t, 'nes', 'NES')
        fdr = _get(t, 'fdr', 'FDR')
        nes = None if nes is None or (isinstance(nes, float) and math.isnan(nes)) else float(nes)
        fdr = None if fdr is None or (isinstance(fdr, float) and math.isnan(fdr)) else float(fdr)
        recs.append((str(_get(t, 'label', 'Term', 'term')), nes, fdr))
    n = len(recs)
    if n == 0:
        print(f'  [SKIP] no terms for {out_base}')
        return None

    lo, hi = tuple(xlim) if xlim is not None else XLIM_DEFAULT
    amax = max([abs(v) for _, v, _ in recs if v is not None], default=0.0)
    if amax > min(abs(lo), abs(hi)):
        new = math.ceil(amax * 10) / 10 + 0.2
        print(f'  [WARN] |NES| = {amax:.2f} exceeds xlim {(lo, hi)}; widened to +/-{new:.1f}')
        lo, hi = -new, new

    W, H = panel_size_mm(n)
    fig = plt.figure(figsize=srep.figsize_mm(W, H))
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    px_mm = fig.dpi / 25.4
    fp = FontProperties()                       # rc font: Arial 8 pt

    # title: one line, top-left
    ttl = fig.text(EDGE_MM / W, 1 - EDGE_MM / H, srep.sentence_case(title),
                   ha='left', va='top')
    fig.canvas.draw()
    tb = ttl.get_window_extent(r)
    if tb.width / px_mm > W - 2 * EDGE_MM:
        print(f'  [WARN] title wider than the panel: "{title}"')
    top_mm = (fig.bbox.height - tb.y0) / px_mm + GAP_MM

    # legend: two rows at the lower left (the two patch labels alone are wider
    # than 90 mm in one row), the FDR note aligned to the same rows at the right
    handles = [mpatches.Patch(color=COL_UP, label=LEGEND_UP),
               mpatches.Patch(color=COL_DOWN, label=LEGEND_DOWN)]
    leg = fig.legend(handles=handles, loc='lower left',
                     bbox_to_anchor=(EDGE_MM / W, EDGE_MM / H), ncol=1, frameon=False,
                     handlelength=1.0, handleheight=0.7, handletextpad=0.4,
                     labelspacing=0.3, borderpad=0, borderaxespad=0)
    fig.canvas.draw()
    for lt, note in zip(leg.get_texts(), LEGEND_NOTE):
        bb = lt.get_window_extent(r)
        _, _, desc = r.get_text_width_height_descent(lt.get_text(), lt.get_fontproperties(), False)
        fig.text(1 - RIGHT_MM / W, (bb.y0 + desc) / fig.bbox.height, note,
                 ha='right', va='baseline')
    leg_top_mm = leg.get_window_extent(r).y1 / px_mm

    # axes: provisional rectangle, then measured x-axis block -> final rectangle
    ax = fig.add_axes([AXES_LEFT_MM / W, 0.4, 1 - (AXES_LEFT_MM + RIGHT_MM) / W, 0.4])
    ax.set_xlim(lo, hi)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.set_xlabel(XLABEL, loc='right')
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_linewidth(srep.LINE_PT)
    ax.tick_params(axis='y', length=0)
    ax.set_ylim(-0.5, n - 0.5)
    ys = [n - 1 - i for i in range(n)]               # first term at the top
    ax.set_yticks(ys)
    labels = [term_label(lab) for lab, _, _ in recs]
    ax.set_yticklabels(labels, linespacing=1.0)
    fig.canvas.draw()
    # x-axis block = tick marks + tick labels + axis label (measured from the
    # axes bottom; the tight bbox alone excludes the tick marks)
    x_block_mm = (ax.get_window_extent(r).y0 - ax.xaxis.get_tightbbox(r).y0) / px_mm
    lab_w_mm = max(t.get_window_extent(r).width for t in ax.get_yticklabels()) / px_mm
    left_mm = AXES_LEFT_MM
    avail = AXES_LEFT_MM - 3.5 / 72 * 25.4 - 0.5    # minus the y tick pad
    if lab_w_mm > avail:
        left_mm = AXES_LEFT_MM + (lab_w_mm - avail)
        print(f'  [WARN] term names wider than the label column ({lab_w_mm:.1f} mm); '
              f'axes shifted right by {left_mm - AXES_LEFT_MM:.1f} mm')
    bottom_mm = leg_top_mm + GAP_MM + x_block_mm
    ax.set_position([left_mm / W, bottom_mm / H,
                     1 - (left_mm + RIGHT_MM) / W, 1 - (bottom_mm + top_mm) / H])
    fig.canvas.draw()

    # bars, zero line
    for (lab, nes, fdr), y in zip(recs, ys):
        if nes is None:
            continue
        sig = fdr is not None and fdr <= SIG
        ax.barh(y, nes, height=BAR_HEIGHT, color=COL_UP if nes > 0 else COL_DOWN,
                alpha=1.0 if sig else ALPHA_NS, edgecolor='none', linewidth=0, zorder=2)
    ax.axvline(0, color='black', lw=srep.LINE_PT, zorder=3)

    # FDR text just beyond the bar tip on the bar's own side; flipped to the other
    # side of the zero line when it would run past the x-limit
    ax_w_px = ax.get_window_extent(r).width
    units_per_px = (hi - lo) / ax_w_px
    pad = FDR_PAD_PT / 72 * fig.dpi * units_per_px
    # decide the text side once per sign so that the FDR column is uniform within a panel
    items = []
    for (lab, nes, fdr), y in zip(recs, ys):
        s = fdr_text(fdr)
        w = r.get_text_width_height_descent(s, fp, False)[0] * units_per_px
        items.append((0.0 if nes is None else nes, y, s, w))
    pos_flip = any(v + pad + w > hi for v, _, _, w in items if v >= 0)
    neg_flip = any(v - pad - w < lo for v, _, _, w in items if v < 0)
    for v, y, s, w in items:
        if v >= 0:
            x, ha = (-pad, 'right') if pos_flip else (v + pad, 'left')
        else:
            x, ha = (pad, 'left') if neg_flip else (v - pad, 'right')
        ax.text(x, y, s, ha=ha, va='center', color='black', zorder=4, clip_on=False)

    _backup(out_base)
    os.makedirs(os.path.dirname(os.path.abspath(out_base)), exist_ok=True)
    srep.save(fig, out_base)
    plt.close(fig)
    print(f'  Saved -> {out_base}.png / .pdf  ({W:.0f} x {H:.0f} mm, {n} terms, '
          f'xlim = ({lo:g}, {hi:g}))')
    return (W, H)
