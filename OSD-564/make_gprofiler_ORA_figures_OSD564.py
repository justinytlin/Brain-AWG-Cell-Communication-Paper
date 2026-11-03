"""
g:Profiler ORA bar panels for OSD-564 (spaceflight vs ground control).

Inputs: the two g:Profiler exports in gprofiler/ (up- and down-regulated DEG queries).
Terms are matched by term_id (GO:xxxxxxx / KEGG:xxxxx) so that e.g. "Ribosome" picks
the KEGG pathway and not the GO:CC term of the same name.  Labels are the real term
names from the export, formatted "term name (ID) [source]".

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
"""
import os
import textwrap
import warnings
import pandas as pd
import numpy as np
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
SUBTITLE = 'input: adjusted p < 0.05 DEGs, 410 up / 496 down; g:SCS-adjusted p'

FIGURES = [
    {
        'key':   'TNTs',
        'title': 'g:Profiler ORA, TNT-related terms, OSD-564: spaceflight vs ground control',
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
        'fname': os.path.join(GPROFILER_DIR, 'gprofiler_ORA_TNTs_pathways.png'),
    },
    {
        'key':   'sEVs',
        'title': 'g:Profiler ORA, sEV-related terms, OSD-564: spaceflight vs ground control',
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
        'fname': os.path.join(GPROFILER_DIR, 'gprofiler_ORA_sEVs_pathways.png'),
    },
    {
        'key':   'mito',
        'title': 'g:Profiler ORA, mitochondrial terms, OSD-564: spaceflight vs ground control',
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
        'fname': os.path.join(GPROFILER_DIR, 'gprofiler_ORA_mito_downregulated_pathways.png'),
    },
]

# ── Styling constants ──────────────────────────────────────────────
UP_COLOR  = '#C62828'   # red  - over-represented among genes upregulated in spaceflight
DN_COLOR  = '#1565C0'   # blue - over-represented among genes downregulated in spaceflight
X_CAP     = 20.0        # bars truncated here; the true value is printed as the label
SIG_THR   = -np.log10(0.05)  # ~1.30
FOOTNOTE  = 'bars truncated at 20; label shows -log10 adjusted p'


def wrap_label(term_name, term_id, source, max_w=46):
    """'term name (GO:xxxxxxx) [GO:BP]' wrapped to max_w characters per line."""
    return '\n'.join(textwrap.wrap(f'{term_name} ({term_id}) [{source}]', max_w))


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

    fig_h = max(4.2, n * 0.62 + 3.0)
    fig, ax = plt.subplots(figsize=(11.5, fig_h))

    # ── Draw bars ──────────────────────────────────────────────────
    for i, rec in enumerate(combined):
        y     = y_pos[i]
        color = UP_COLOR if rec['direction'] == 'up' else DN_COLOR
        alpha = 1.00 if rec['padj'] <= 0.05 else 0.45
        val   = rec['val']
        bw    = min(val, X_CAP)

        ax.barh(y, bw, color=color, alpha=alpha,
                edgecolor='white', height=0.70, zorder=2)

        if val > X_CAP:
            # visible truncation mark: two white slanted lines across the bar end
            for dx in (-0.55, -0.15):
                ax.plot([X_CAP + dx - 0.25, X_CAP + dx + 0.25], [y - 0.42, y + 0.42],
                        color='white', lw=2.6, solid_capstyle='butt', zorder=3)
            ax.text(X_CAP + 0.35, y, f'{val:.1f}  (truncated)',
                    va='center', ha='left', fontsize=10, color=color,
                    fontweight='bold', zorder=4)
        else:
            sig_str = ' *' if rec['padj'] <= 0.05 else ''
            ax.text(bw + 0.25, y, f'{val:.2f}{sig_str}',
                    va='center', ha='left', fontsize=10, color='#333333', zorder=4)

    # ── Group separator ────────────────────────────────────────────
    n_up, n_dn = len(up_recs), len(dn_recs)
    if n_up > 0 and n_dn > 0:
        sep_y = y_pos[n_up - 1] - 0.5
        ax.axhline(sep_y, color='#777777', lw=1.2, linestyle='--', alpha=0.8, zorder=3)

    # ── Significance threshold line ────────────────────────────────
    ax.axvline(SIG_THR, color='#666666', lw=0.9, linestyle=':', alpha=0.85, zorder=1)
    ax.text(SIG_THR + 0.12, n - 0.35, 'adj. p = 0.05', va='top', ha='left',
            fontsize=9, color='#666666', style='italic')

    # ── Axis formatting ────────────────────────────────────────────
    ax.set_xlim(0, X_CAP + 5.5)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_yticks([y_pos[i] for i in range(n)])
    ax.set_yticklabels(y_labels, fontsize=10.5)

    xticks = list(range(0, int(X_CAP) + 1, 5))
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(t) for t in xticks], fontsize=10)
    # visible cap marker on the axis
    ax.axvline(X_CAP, color='#999999', lw=0.9, linestyle='--', alpha=0.9, zorder=1)
    ax.text(X_CAP, -0.62, 'cap', ha='center', va='bottom', fontsize=8.5, color='#777777')

    ax.set_xlabel('$-\\log_{10}$(adjusted p-value), g:Profiler ORA', fontsize=11)
    ax.set_title(f"{fig_def['title']}\n{SUBTITLE}", fontsize=11.5, fontweight='bold', pad=10)

    ax.xaxis.grid(True, color='#eeeeee', lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0, pad=5)

    # ── Legend + footnote ─────────────────────────────────────────
    legend_handles = [
        mpatches.Patch(color=UP_COLOR, label='Over-represented among genes upregulated in spaceflight'),
        mpatches.Patch(color=DN_COLOR, label='Over-represented among genes downregulated in spaceflight'),
    ]
    fig.legend(handles=legend_handles, fontsize=9.5, frameon=False,
               loc='lower center', bbox_to_anchor=(0.5, 0.25 / fig_h), ncol=1)
    fig.text(0.5, 0.04 / fig_h, FOOTNOTE + ';  * = adjusted p ≤ 0.05',
             ha='center', va='bottom', fontsize=9, color='#555555', style='italic')

    # reserve ~0.95 in at the bottom for the legend + footnote
    plt.tight_layout(rect=(0, 0.95 / fig_h, 1, 1))
    plt.savefig(fig_def['fname'], dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved -> {fig_def['fname']}")

print("\nDone.")
