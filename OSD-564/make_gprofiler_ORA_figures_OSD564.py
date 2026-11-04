"""
g:Profiler ORA bar panels for OSD-564 (spaceflight vs ground control), Scientific
Reports style (90 x 65 mm, 8 pt Arial, no truncated bars).

Inputs: the two g:Profiler exports in gprofiler/ (up- and down-regulated DEG queries).
Terms are matched by term_id (GO:xxxxxxx / KEGG:xxxxx) so that e.g. "Ribosome" picks
the KEGG pathway and not the GO:CC term of the same name.  Labels are the real term
names from the export, formatted "Term name (ID) [source]" in sentence case.

Panels
  TNTs : up  = actin cytoskeleton organization; plasma membrane bounded cell projection
               organization; small GTPase binding; intracellular transport; cell-cell signaling
         down= translation at synapse; cytosolic ribosome (GO:CC); Ribosome (KEGG)
  sEVs : up  = vesicle-mediated transport; regulated exocytosis; vesicle organization;
               vesicle localization; extracellular organelle; protein localization to
               plasma membrane
         down= (none: no vesicle-specific term is present in the downregulated export)
  mito : down= oxidative phosphorylation; aerobic respiration; mitochondrial respiratory
               chain complex assembly; mitochondrial ATP synthesis coupled electron
               transport; proton motive force-driven mitochondrial ATP synthesis;
               mitochondrial respiratory chain complex I assembly; mitochondrial
               cytochrome c oxidase assembly; mitochondrial inner membrane (GO:CC);
               mitochondrial envelope (GO:CC)

Generic terms enriched in both directions ("organelle organization", "synapse",
"intracellular membrane-bounded organelle") are intentionally excluded.

Figure-legend text (not drawn in the panel): input = adjusted p < 0.05 DEGs,
410 up / 496 down; g:SCS-adjusted p; dashed vertical line = adjusted p = 0.05;
dashed horizontal line separates the up- and down-regulated queries.
"""
import os
import warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import srep_style as srep; srep.apply()
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator

from osd564_nes_panel import (label_case, wrap_lines, fit_xlim_to_texts,
                              xlabel_no_clip, STUDY_TAG)

warnings.filterwarnings('ignore')

# ── Paths ──────────────────────────────────────────────────────────
GPROFILER_DIR = 'gprofiler'
UP_CSV = os.path.join(
    GPROFILER_DIR,
    'Supplemental DEG Tables (OSD-562, OSD-564, OSD-352)'
    ' - OSD-564 g_Profiler Upregulated Pathways.csv',
)
DOWN_CSV = os.path.join(
    GPROFILER_DIR,
    'Supplemental DEG Tables (OSD-562, OSD-564, OSD-352)'
    ' - OSD-564 g_Profiler Downregulated Pathways.csv',
)

# ── Load ORA results ───────────────────────────────────────────────
print("Loading g:Profiler ORA results ...")
up_df   = pd.read_csv(UP_CSV,   skiprows=2)
down_df = pd.read_csv(DOWN_CSV, skiprows=2)
print(f"  up: {len(up_df)} terms   down: {len(down_df)} terms")


def make_lookup(df):
    """term_id -> row (first occurrence)."""
    lu = {}
    for _, row in df.iterrows():
        key = str(row['term_id']).strip()
        if key not in lu:
            lu[key] = row
    return lu


up_lu   = make_lookup(up_df)
down_lu = make_lookup(down_df)


def neglog10_padj(row):
    """Prefer the exported column; fall back to computing it from adjusted_p_value."""
    v = row.get('negative_log10_of_adjusted_p_value', np.nan)
    if pd.isna(v):
        v = -np.log10(float(row['adjusted_p_value']))
    return float(v)


# ── Figure definitions (term_id, expected term_name for sanity check) ──
FIGURES = [
    {
        'key':   'TNTs',
        'up_terms': [
            ('GO:0030036', 'actin cytoskeleton organization'),
            ('GO:0120036', 'plasma membrane bounded cell projection organization'),
            ('GO:0031267', 'small GTPase binding'),
            ('GO:0046907', 'intracellular transport'),
            ('GO:0007267', 'cell-cell signaling'),
        ],
        'dn_terms': [
            ('GO:0140241', 'translation at synapse'),
            ('GO:0022626', 'cytosolic ribosome'),
            ('KEGG:03010', 'Ribosome'),
        ],
        'fname': os.path.join(GPROFILER_DIR, 'gprofiler_ORA_TNTs_pathways'),
    },
    {
        'key':   'sEVs',
        'up_terms': [
            ('GO:0016192', 'vesicle-mediated transport'),
            ('GO:0045055', 'regulated exocytosis'),
            ('GO:0016050', 'vesicle organization'),
            ('GO:0051648', 'vesicle localization'),
            ('GO:0043230', 'extracellular organelle'),
            ('GO:0072659', 'protein localization to plasma membrane'),
        ],
        # No vesicle-specific term (synaptic vesicle / exocytosis / secretory vesicle)
        # is present in the downregulated export, so the down side is empty.
        'dn_terms': [],
        'fname': os.path.join(GPROFILER_DIR, 'gprofiler_ORA_sEVs_pathways'),
    },
    {
        'key':   'mito',
        'up_terms': [],
        'dn_terms': [
            ('GO:0006119', 'oxidative phosphorylation'),
            ('GO:0009060', 'aerobic respiration'),
            ('GO:0033108', 'mitochondrial respiratory chain complex assembly'),
            ('GO:0042775', 'mitochondrial ATP synthesis coupled electron transport'),
            ('GO:0042776', 'proton motive force-driven mitochondrial ATP synthesis'),
            ('GO:0032981', 'mitochondrial respiratory chain complex I assembly'),
            ('GO:0033617', 'mitochondrial cytochrome c oxidase assembly'),
            ('GO:0005743', 'mitochondrial inner membrane'),
            ('GO:0005740', 'mitochondrial envelope'),
        ],
        'fname': os.path.join(GPROFILER_DIR, 'gprofiler_ORA_mito_downregulated_pathways'),
    },
]

# ── Styling constants ──────────────────────────────────────────────
UP_COLOR    = srep.COLORS['up']        # over-represented among genes upregulated in spaceflight
DN_COLOR    = srep.COLORS['down']      # over-represented among genes downregulated in spaceflight
UP_COLOR_NS = srep.COLORS['nom_up']
DN_COLOR_NS = srep.COLORS['nom_down']
GREY        = '#777777'
SIG_THR     = -np.log10(0.05)          # ~1.30
WRAP        = 38                       # characters per label line


def wrap_label(term_name, term_id, source, min_w=WRAP):
    """'Term name (GO:xxxxxxx) [GO:BP]' in sentence case, wrapped to at most two
    lines of at least min_w characters."""
    return wrap_lines(label_case(f'{term_name} ({term_id}) [{source}]'), min_w)


def collect(term_list, lookup, direction):
    recs = []
    for term_id, expected_name in term_list:
        row = lookup.get(term_id)
        if row is None:
            print(f"  [MISSING {direction.upper()}] {term_id} ({expected_name})")
            continue
        name = str(row['term_name']).strip()
        if name.lower() != expected_name.lower():
            print(f"  [NAME CHECK] {term_id}: export says '{name}', script expected '{expected_name}'")
        recs.append({
            'name':      name,
            'id':        term_id,
            'source':    str(row['source']).strip(),
            'val':       neglog10_padj(row),
            'padj':      float(row['adjusted_p_value']),
            'n':         int(row['intersection_size']),
            'direction': direction,
        })
    return recs


# ── Generate one figure per category ──────────────────────────────
for fig_def in FIGURES:
    print(f"\nGenerating: {fig_def['key']}")

    up_recs = collect(fig_def['up_terms'], up_lu,   'up')
    dn_recs = collect(fig_def['dn_terms'], down_lu, 'dn')

    # Most significant at the top within each direction
    up_recs.sort(key=lambda x: x['val'], reverse=True)
    dn_recs.sort(key=lambda x: x['val'], reverse=True)
    combined = up_recs + dn_recs
    n = len(combined)
    if n == 0:
        print("  No data - skipping.")
        continue

    for r in combined:
        print(f"    {r['direction']:>2s}  {r['val']:6.2f}  {r['name']} ({r['id']}) [{r['source']}]  n={r['n']}")

    # top row = first record  ->  reverse the y positions
    y_pos = {i: n - 1 - i for i in range(n)}
    y_labels = [wrap_label(r['name'], r['id'], r['source']) for r in combined]

    # 90 x 65 mm holds up to 7 rows of two-line 8 pt labels; the 8- and 9-term
    # panels need ~6 mm per row plus ~19 mm for title, x-axis and a one-row
    # legend (+3.5 mm for a second legend row)
    w_mm, h_mm = srep.PANEL_MM['ora']
    n_legend_rows = int(bool(up_recs)) + int(bool(dn_recs))
    h_mm = max(h_mm, 19 + 6 * n + 3.5 * (n_legend_rows - 1))
    fig, ax = plt.subplots(figsize=srep.figsize_mm(w_mm, h_mm), constrained_layout=True)

    # ── Bars (full length, never truncated) + printed value ─────────
    anchored = []
    for i, rec in enumerate(combined):
        y   = y_pos[i]
        sig = rec['padj'] <= 0.05
        if rec['direction'] == 'up':
            color = UP_COLOR if sig else UP_COLOR_NS
        else:
            color = DN_COLOR if sig else DN_COLOR_NS
        val = rec['val']
        ax.barh(y, val, color=color, edgecolor='none', linewidth=0, height=0.70, zorder=2)
        t = ax.annotate(f'{val:.2f}', xy=(val, y), xytext=(3, 0), textcoords='offset points',
                        va='center', ha='left', color='black', zorder=4, annotation_clip=False)
        t.set_in_layout(False)
        anchored.append((t, val, 'right'))

    # ── Group separator (up vs down queries) ───────────────────────
    n_up, n_dn = len(up_recs), len(dn_recs)
    if n_up > 0 and n_dn > 0:
        sep_y = y_pos[n_up - 1] - 0.5
        ax.axhline(sep_y, color=GREY, lw=srep.LINE_PT, linestyle='--', zorder=3)

    # ── Significance threshold (adjusted p = 0.05) ─────────────────
    ax.axvline(SIG_THR, color=GREY, lw=srep.LINE_PT, linestyle='--', zorder=1)

    # ── Axis formatting ────────────────────────────────────────────
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_yticks([y_pos[i] for i in range(n)])
    ax.set_yticklabels(y_labels, linespacing=1.0)
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
    fig.suptitle(STUDY_TAG, x=0.01, ha='left')

    # ── Legend (only the directions present) ──────────────────────
    handles = []
    if n_up:
        handles.append(mpatches.Patch(color=UP_COLOR,
                       label='Over-represented among genes upregulated in spaceflight'))
    if n_dn:
        handles.append(mpatches.Patch(color=DN_COLOR,
                       label='Over-represented among genes downregulated in spaceflight'))
    fig.legend(handles=handles, loc='outside lower left', ncol=1, frameon=False)

    x_max = max(r['val'] for r in combined)
    fit_xlim_to_texts(fig, ax, anchored, 0.0, x_max, lo_fixed=0.0)
    xlabel_no_clip(fig, ax, '-log10(adjusted p)')
    fit_xlim_to_texts(fig, ax, anchored, 0.0, x_max, lo_fixed=0.0)

    srep.save(fig, fig_def['fname'])
    plt.close(fig)
    print(f"  Saved -> {fig_def['fname']}.png / .pdf   ({w_mm} x {h_mm} mm)")

print("\nDone.")
