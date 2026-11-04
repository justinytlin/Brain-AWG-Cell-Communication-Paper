"""
KEGG GSEA enrichment plots + combined NES summary for the harmonized OSD-564 panels.

Harmonized KEGG sets (shared with OSD-562 / OSD-352; exact KEGG names, fixed order):
  sEV       : Sphingolipid signaling pathway; Sphingolipid metabolism; Phospholipase D
              signaling pathway; Endocytosis; Synaptic vesicle cycle;
              Phosphatidylinositol signaling system; Lysosome
  TNT       : Axon guidance; Regulation of actin cytoskeleton; Rap1 signaling pathway
  Metabolic : AMPK signaling pathway; Glycolysis / Gluconeogenesis; Oxidative phosphorylation

The former "Abeta / Tau" framing (Alzheimer disease plot mis-titled OSD-562) was
dropped from the panels.  "Alzheimer disease" is still drawn as a plain enrichment
plot so that the old mis-titled file is overwritten with a correctly titled one.

GSEA prerank is re-run (seed 42) only to obtain the running-ES profile for the
enrichment plots; NES / FDR / p shown in every figure are taken from the saved
KEGG_results/GSEA_KEGG_results.csv so the figures match the tables.
Pass ``--nes-only`` to skip the prerank re-run and the enrichment plots and only
redraw the combined NES summary panel (shared srep_nes_panel.nes_panel).
"""
import os
import sys
import warnings
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import srep_style as srep; srep.apply()
import matplotlib.pyplot as plt
import matplotlib.text
import matplotlib.lines
import matplotlib.collections
from matplotlib.ticker import FuncFormatter

from srep_nes_panel import nes_panel, label_case, XLIM

warnings.filterwarnings('ignore')

NES_ONLY  = '--nes-only' in sys.argv[1:]
STUDY_TAG = 'OSD-564, spaceflight vs ground control (n = 5 vs 5)'   # enrichment-plot title line
ENRICH_SIZE_MM = (90, 80)   # running-ES enrichment plots (not assembled into a figure)


def restyle_gseaplot(fig, nes, pval, fdr):
    """gseapy hard-codes 14-16 pt bold lettering, 0.5 pt lines, 'Pval: 0.000e+00'
    statistics and 8 % margins; bring the figure to the journal style (8 pt regular,
    >= 1 pt lines, sentence case, journal number formats, nothing clipped)."""
    stat_text = {
        'NES:':  f'NES = {nes:.2f}',
        'Pval:': ('p < 0.001' if pval < 0.001 else f'p = {pval:.3f}'),
        'FDR:':  srep.fdr_label(fdr),
    }
    for t in fig.findobj(matplotlib.text.Text):
        t.set_fontsize(srep.FONT_SIZE)
        t.set_fontweight('normal')
        for key, new in stat_text.items():
            if t.get_text().startswith(key):
                t.set_text(new)
        if t.get_text().startswith('Zero score at'):
            t.set_text('Zero score at ' + srep.fmt_int(float(t.get_text().split()[-1])))
    for l in fig.findobj(matplotlib.lines.Line2D):
        if l.get_linewidth() < srep.LINE_PT:
            l.set_linewidth(srep.LINE_PT)
    for c in fig.findobj(matplotlib.collections.LineCollection):
        c.set_linewidth(srep.LINE_PT)
    for ax in fig.axes:
        ax.set_xlabel(label_case(ax.get_xlabel()))
        ax.set_ylabel(label_case(ax.get_ylabel()))
        if ax.get_xlabel() == 'Gene rank':
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: srep.fmt_int(v)))
        for sp in ax.spines.values():
            sp.set_linewidth(srep.LINE_PT)
        # gseapy stacks the axes in [0.1, 0.9] of the figure; remap so that the
        # x-axis label and the two-line title are not clipped at 90 x 80 mm
        pos = ax.get_position()
        y0 = 0.15 + (pos.y0 - 0.1) * (0.88 - 0.15) / 0.8
        h  = pos.height * (0.88 - 0.15) / 0.8
        ax.set_position([0.15, y0, 0.82, h])
    if fig._suptitle is not None:
        fig._suptitle.set_fontsize(srep.FONT_SIZE)
        fig._suptitle.set_fontweight('normal')
        fig._suptitle.set_y(0.985)

# ──────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────
PANELS = {
    'sEV': [
        "Sphingolipid signaling pathway",
        "Sphingolipid metabolism",
        "Phospholipase D signaling pathway",
        "Endocytosis",
        "Synaptic vesicle cycle",
        "Phosphatidylinositol signaling system",
        "Lysosome",
    ],
    'TNT': [
        "Axon guidance",
        "Regulation of actin cytoskeleton",
        "Rap1 signaling pathway",
    ],
    'Metabolic adaptation': [
        "AMPK signaling pathway",
        "Glycolysis / Gluconeogenesis",
        "Oxidative phosphorylation",
    ],
}
TARGET_PATHWAYS = []
for lst in PANELS.values():
    for p in lst:
        if p not in TARGET_PATHWAYS:
            TARGET_PATHWAYS.append(p)

# Extra enrichment plots that are not part of any NES panel
EXTRA_ENRICHMENT_PLOTS = ["Alzheimer disease"]

OUT_DIR    = 'KEGG_results/pathway_figures'
GMT_FILE   = 'KEGG_results/GSEA_prerank/gene_sets.gmt'
RANKED_CSV = 'KEGG_results/GSEA_ranked_gene_list.csv'
RESULT_CSV = 'KEGG_results/GSEA_KEGG_results.csv'
PERM_NUM   = 1000

os.makedirs(OUT_DIR, exist_ok=True)

# ──────────────────────────────────────────────────────────────────
# 1. Saved GSEA statistics (source of truth for NES / FDR)
# ──────────────────────────────────────────────────────────────────
print("Loading saved GSEA KEGG results ...")
res = pd.read_csv(RESULT_CSV, index_col=0)
res['Term'] = res['Term'].str.strip()
res = res.set_index('Term')
print(f"  {len(res)} pathways in {RESULT_CSV}")


def enrichment_plots():
    """Re-run prerank (seed 42) for the running-ES profiles and draw one
    enrichment plot per harmonized pathway (NES / p / FDR from the saved table)."""
    import gseapy as gp
    from gseapy import gseaplot

    # ── 2. Load ranked gene list and re-run prerank (for the ES profiles) ──
    print("Loading ranked gene list ...")
    rnk = pd.read_csv(RANKED_CSV, index_col=0).squeeze('columns')
    rnk.name = 'rank_metric'
    rnk = rnk.sort_values(ascending=False)
    rnk = rnk[~rnk.index.duplicated(keep='first')]
    print(f"  {len(rnk):,} genes  |  range [{rnk.min():.2f}, {rnk.max():.2f}]")

    print(f"\nRunning GSEA prerank (permutations = {PERM_NUM}) ...")
    gsea_res = gp.prerank(
        rnk=rnk,
        gene_sets=GMT_FILE,
        outdir=None,
        seed=42,
        permutation_num=PERM_NUM,
        min_size=10,
        max_size=500,
        threads=4,
        verbose=False,
    )
    print(f"  Done. {len(gsea_res.results)} pathways computed.")

    # ── 3. Enrichment plots (exact KEGG name match) ──
    print("\nGenerating enrichment plots ...")
    found, missing = 0, 0

    for term in TARGET_PATHWAYS + EXTRA_ENRICHMENT_PLOTS:
        if term not in gsea_res.results or term not in res.index:
            print(f"  [MISSING]  '{term}'")
            missing += 1
            continue

        safe_name = (term.replace('/', '-').replace(' ', '_')
                     .replace(':', '').replace(',', '')
                     .replace('(', '').replace(')', ''))
        obase = f"{OUT_DIR}/GSEA_{safe_name}"

        saved = res.loc[term]
        plot_kwargs = dict(gsea_res.results[term])
        # keep the statistics identical to the saved table
        plot_kwargs['nes']  = float(saved['NES'])
        plot_kwargs['pval'] = float(saved['NOM p-val'])
        plot_kwargs['fdr']  = float(saved['FDR q-val'])

        try:
            axes = gseaplot(
                rank_metric=gsea_res.ranking,
                term=f'KEGG: {label_case(term)}\n{STUDY_TAG}',
                ofname=None,
                figsize=srep.figsize_mm(*ENRICH_SIZE_MM),
                **plot_kwargs,
            )
            fig = axes[0].figure
            restyle_gseaplot(fig, plot_kwargs['nes'], plot_kwargs['pval'], plot_kwargs['fdr'])
            srep.save(fig, obase)
            plt.close('all')
            print(f"  [OK]  {term}  ->  {obase}.png / .pdf")
            found += 1
        except Exception as e:
            print(f"  [ERROR]  '{term}': {e}")
            missing += 1

    print(f"\nFinished.  {found} plots saved to:  {OUT_DIR}/")
    if missing:
        print(f"           {missing} pathway(s) not found or failed.")


if NES_ONLY:
    print("--nes-only: skipping the prerank re-run and the enrichment plots.")
else:
    enrichment_plots()

# ──────────────────────────────────────────────────────────────────
# 4. Combined summary figure — NES bar panel for the 13 harmonized pathways
#    (harmonized order, first term at the top; 13 terms -> 90 x 96 mm)
# ──────────────────────────────────────────────────────────────────
print("\nCreating combined pathway summary figure ...")
terms = []
for panel, lst in PANELS.items():
    for term in lst:
        if term not in res.index:
            print(f"  [MISSING] {term}")
            continue
        row = res.loc[term]
        terms.append({'label': f'{term} [{panel}]',
                      'nes': float(row['NES']), 'fdr': float(row['FDR q-val'])})

nes_panel(terms, 'KEGG sEV, TNT and metabolic pathways, OSD-564 (n = 5 vs 5)',
          f'{OUT_DIR}/GSEA_selected_pathways_NES_summary', xlim=XLIM['metabolic'])
print("Done.")
