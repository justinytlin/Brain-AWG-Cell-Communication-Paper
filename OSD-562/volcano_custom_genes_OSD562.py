"""
volcano_custom_genes_OSD562.py
------------------------------
Volcano plot for OSD-562 (ISS-T spaceflight vs ground control, re-run
DESeq2 4 vs 3; OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3.csv) with a user-specified
set of genes highlighted, a curated subset of them labelled (italic gene
symbols, 1 pt leader lines), plus labels on ALL adjusted-significant genes
(adj p <= 0.05, |log2FC| >= 1).

Pipe-separated entries (e.g. "Ndufb4|Ndufb4c|Ndufb4b") are resolved by
trying each alias in turn against the DE table SYMBOL column, which is
itself split on "|" so pipe-joined symbols in the table also match.

Style: Scientific Reports (srep_style.py) -- Arial 8 pt everywhere, sentence
case, no line < 1 pt; 90 x 85 mm panel saved as 300 dpi PNG + vector PDF.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import srep_style as srep
srep.apply()
import srep_panels                      # shared PCA / volcano / QC panel drawers (identical in all datasets)
import warnings
warnings.filterwarnings('ignore')

# ──────────────────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────────────────
DE_FILE  = 'OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3.csv'
COMP_TAG = '(Space Flight)v(Ground Control)'
FC_THR   = 1.0    # |log2FC| threshold for dotted lines
ADJ_THR  = 0.05
NOM_THR  = 0.05
DATASET_LINE = 'OSD-562, spaceflight vs ground control (n = 4 vs 3)'

# Genes to highlight -- pipe-separated entries = aliases (first match wins)
ANNOTATION_ENTRIES = [
    "Adam8", "Apoc1", "Aqp1", "Arpc2", "Arap2", "Baiap2l2", "Bcat1",
    "C4a", "Casr", "Ect2", "Epcam", "Folr1", "Gas2", "Gas2l3", "Ghsr",
    "Itga2", "Lpl", "Myh7b", "Rac2", "Rhof", "S100a4", "S100a8",
    "Smtnl2", "Timp1", "Tmsb4x", "Tnni2", "Wdr1", "Wnt7a", "Zfas1",
    "Mrps25", "Chchd1", "Rmrp",
    "Ndufb4|Ndufb4c|Ndufb4b",
    "Ndufaf6",
    "Atp5mk|Gm35857",
    "Pet100",
    "COX2",
    # OSD-562-specific genes of interest
    "Cpt1b", "Pprc1", "Sdc4", "Vps4a", "Ist1",
    "Ttn", "Hmga1b", "Smarcd1",
]
# Curated mitochondrial genes (always labelled)
CURATED_MITO = ["Mrps25", "Chchd1", "Rmrp", "Ndufb4|Ndufb4c|Ndufb4b", "Ndufaf6",
                "Atp5mk|Gm35857", "Pet100", "COX2"]
# Labelled subset: 90 x 85 mm is too small to label every highlighted gene at 8 pt
# without overlaps, so only these (plus all adjusted-significant genes) carry a label.
LABEL_ENTRIES = ["Cpt1b", "Pprc1", "Sdc4", "Rhof", "Ghsr", "Vps4a", "Ist1", "Rac2",
                 "Baiap2l2", "Ttn", "Hmga1b", "Smarcd1"] + CURATED_MITO
# Display-name overrides (table symbol -> label on the figure)
DISPLAY_OVERRIDES = {'COX2': 'mt-Co2', 'COX1': 'mt-Co1', 'COX3': 'mt-Co3'}


# ──────────────────────────────────────────────────────────────────────────────
# 1. Load DE data
# ──────────────────────────────────────────────────────────────────────────────
print("Loading DE data (re-run 4 vs 3) ...")
de = pd.read_csv(DE_FILE, index_col=0, low_memory=False)

df_v = de[['SYMBOL',
           f'Log2fc_{COMP_TAG}',
           f'P.value_{COMP_TAG}',
           f'Adj.p.value_{COMP_TAG}']].copy()
df_v.columns = ['symbol', 'log2fc', 'pval', 'adjpval']
for col in ['log2fc', 'pval', 'adjpval']:
    df_v[col] = pd.to_numeric(df_v[col], errors='coerce')
df_v = df_v.dropna(subset=['log2fc', 'pval', 'adjpval'])
df_v['-log10p'] = -np.log10(df_v['pval'].clip(lower=1e-300))
print(f"DE table loaded: {srep.fmt_int(len(df_v))} genes")

# ──────────────────────────────────────────────────────────────────────────────
# 2. Significance categories
# ──────────────────────────────────────────────────────────────────────────────
adj_sig_up   = (df_v['log2fc'] >=  FC_THR) & (df_v['adjpval'] <= ADJ_THR)
adj_sig_down = (df_v['log2fc'] <= -FC_THR) & (df_v['adjpval'] <= ADJ_THR)
nom_sig_up   = (df_v['log2fc'] >=  FC_THR) & (df_v['pval'] <= NOM_THR) & ~adj_sig_up
nom_sig_down = (df_v['log2fc'] <= -FC_THR) & (df_v['pval'] <= NOM_THR) & ~adj_sig_down
not_sig      = ~(adj_sig_up | adj_sig_down | nom_sig_up | nom_sig_down)

n_adj_up, n_adj_down = int(adj_sig_up.sum()), int(adj_sig_down.sum())
n_nom_up, n_nom_down = int(nom_sig_up.sum()), int(nom_sig_down.sum())
print(f"Adj.p<=0.05 & |FC|>=1 -- UP: {n_adj_up}  DOWN: {n_adj_down}")
print(f"Nom.p<=0.05 & |FC|>=1 -- UP: {n_nom_up}  DOWN: {n_nom_down}")

# ──────────────────────────────────────────────────────────────────────────────
# 3. Resolve annotation genes (pipe-separated aliases on both sides)
# ──────────────────────────────────────────────────────────────────────────────
symbol_sets = df_v['symbol'].fillna('').apply(
    lambda s: {a.strip().lower() for a in str(s).split('|') if a.strip()})

found_rows, not_found = [], []
for entry in ANNOTATION_ENTRIES:
    aliases = [a.strip() for a in entry.split('|')]
    matched, display_label = None, aliases[0]
    for alias in aliases:
        hits = df_v[symbol_sets.apply(lambda ss, a=alias.lower(): a in ss)]
        if not hits.empty:
            matched = hits.loc[hits['pval'].idxmin()]   # lowest p if several rows
            display_label = DISPLAY_OVERRIDES.get(alias, alias)
            break
    if matched is not None:
        row_copy = matched.copy()
        row_copy['display_label'] = display_label
        row_copy['ensembl'] = matched.name
        row_copy['labelled'] = entry in LABEL_ENTRIES
        found_rows.append(row_copy)
    else:
        not_found.append(entry)

df_annot = pd.DataFrame(found_rows).reset_index(drop=True)
print(f"\nAnnotation genes found: {len(df_annot)} / {len(ANNOTATION_ENTRIES)}")
if not_found:
    print(f"Not found in DE table: {not_found}")

# All adjusted-significant genes get a label too
df_adj = df_v[adj_sig_up | adj_sig_down].copy()
df_adj['display_label'] = [
    DISPLAY_OVERRIDES.get(str(s).split('|')[0], str(s).split('|')[0]) if isinstance(s, str) else idx
    for s, idx in zip(df_adj['symbol'], df_adj.index)]
df_adj['ensembl'] = df_adj.index
df_adj['labelled'] = True
df_adj = df_adj[~df_adj['ensembl'].isin(df_annot['ensembl'] if len(df_annot) else [])]
print(f"Adjusted-significant genes to label: {len(df_adj)}")

# ──────────────────────────────────────────────────────────────────────────────
# 4. Draw volcano plot (90 x 85 mm)
# ──────────────────────────────────────────────────────────────────────────────
y_col = '-log10p'
COUNTS = {'n_nom_up': n_nom_up, 'n_nom_down': n_nom_down,
          'n_adj_up': n_adj_up, 'n_adj_down': n_adj_down}
assert COUNTS == srep_panels.volcano_counts(df_v)

# ──────────────────────────────────────────────────────────────────────────────
# 5. Gene labels: the LABEL_ENTRIES subset (+ every adjusted-significant gene),
#    italic 8 pt, placed by adjustText with 1 pt leader lines (shared panel drawer)
# ──────────────────────────────────────────────────────────────────────────────
all_df = pd.concat([df_annot, df_adj], ignore_index=True) if len(df_adj) else df_annot
label_df = all_df[all_df['labelled']].reset_index(drop=True)
label_df = label_df.rename(columns={'display_label': 'display'})

# ──────────────────────────────────────────────────────────────────────────────
# 6. Draw + save (shared 90 x 85 mm volcano panel; curated genes = orange circles)
# ──────────────────────────────────────────────────────────────────────────────
srep_panels.volcano_panel(df_v, 'volcano_supp', DATASET_LINE, label_genes=label_df,
                          highlight_genes=df_annot, out_base='OSD562_Volcano_custom_genes',
                          counts=COUNTS)

print("\nAnnotation gene summary:")
summary = all_df[['display_label', 'ensembl', 'log2fc', 'pval', 'adjpval', 'labelled']].copy()
summary.columns = ['Gene', 'ENSEMBL', 'log2FC', 'p-value', 'adj. p-value', 'labelled_on_figure']
summary = summary.sort_values('p-value')
summary.to_csv('OSD562_Volcano_custom_genes_table.csv', index=False)
print(summary.to_string(index=False))
