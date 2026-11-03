"""
volcano_custom_genes_OSD562.py
──────────────────────────────
Volcano plot for OSD-562 (ISS-T spaceflight vs ground control, re-run
DESeq2 4 vs 3; OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3.csv) with a user-specified
set of genes highlighted and labelled, plus labels on ALL adjusted-
significant genes (adj p <= 0.05, |log2FC| >= 1).

Pipe-separated entries (e.g. "Ndufb4|Ndufb4c|Ndufb4b") are resolved by
trying each alias in turn against the DE table SYMBOL column, which is
itself split on "|" so pipe-joined symbols in the table also match.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'Arial',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300,
})

# ──────────────────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────────────────
DE_FILE  = 'OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3.csv'
COMP_TAG = '(Space Flight)v(Ground Control)'
FC_THR   = 1.0    # |log2FC| threshold for dashed lines
ADJ_THR  = 0.05
NOM_THR  = 0.05

# Genes to highlight — pipe-separated entries = aliases (first match wins)
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
]
# Display-name overrides (table symbol -> label on the figure)
DISPLAY_OVERRIDES = {'COX2': 'mt-Co2', 'COX1': 'mt-Co1', 'COX3': 'mt-Co3'}

# ──────────────────────────────────────────────────────────────────────────────
# 1. Load DE data
# ──────────────────────────────────────────────────────────────────────────────
print("Loading DE data (re-run 4 vs 3) …")
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
print(f"DE table loaded: {len(df_v):,} genes")

# ──────────────────────────────────────────────────────────────────────────────
# 2. Significance categories
# ──────────────────────────────────────────────────────────────────────────────
adj_sig_up   = (df_v['log2fc'] >=  FC_THR) & (df_v['adjpval'] <= ADJ_THR)
adj_sig_down = (df_v['log2fc'] <= -FC_THR) & (df_v['adjpval'] <= ADJ_THR)
nom_sig_up   = (df_v['log2fc'] >=  FC_THR) & (df_v['pval'] <= NOM_THR) & ~adj_sig_up
nom_sig_down = (df_v['log2fc'] <= -FC_THR) & (df_v['pval'] <= NOM_THR) & ~adj_sig_down
not_sig      = ~(adj_sig_up | adj_sig_down | nom_sig_up | nom_sig_down)

n_adj_up, n_adj_down = adj_sig_up.sum(), adj_sig_down.sum()
n_nom_up, n_nom_down = nom_sig_up.sum(), nom_sig_down.sum()
print(f"Adj.p<=0.05 & |FC|>=1 — UP: {n_adj_up}  DOWN: {n_adj_down}")
print(f"Nom.p<=0.05 & |FC|>=1 — UP: {n_nom_up}  DOWN: {n_nom_down}")

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
df_adj = df_adj[~df_adj['ensembl'].isin(df_annot['ensembl'] if len(df_annot) else [])]
print(f"Adjusted-significant genes to label: {len(df_adj)}")

# ──────────────────────────────────────────────────────────────────────────────
# 4. Draw volcano plot
# ──────────────────────────────────────────────────────────────────────────────
y_col = '-log10p'
fig, ax = plt.subplots(figsize=(14, 11))
fig.subplots_adjust(top=0.88)

ax.scatter(df_v.loc[not_sig, 'log2fc'], df_v.loc[not_sig, y_col],
           s=7, color='#BDBDBD', alpha=0.35, linewidths=0, rasterized=True,
           label='Not significant')
ax.scatter(df_v.loc[nom_sig_up, 'log2fc'], df_v.loc[nom_sig_up, y_col],
           s=18, color='#EF9A9A', alpha=0.75, linewidths=0.3,
           edgecolors='#C62828', rasterized=True,
           label=f'Nom. sig. up in spaceflight (p<=0.05, |log2FC|>=1; n={n_nom_up})')
ax.scatter(df_v.loc[nom_sig_down, 'log2fc'], df_v.loc[nom_sig_down, y_col],
           s=18, color='#90CAF9', alpha=0.75, linewidths=0.3,
           edgecolors='#1565C0', rasterized=True,
           label=f'Nom. sig. down in spaceflight (p<=0.05, |log2FC|>=1; n={n_nom_down})')
ax.scatter(df_v.loc[adj_sig_up, 'log2fc'], df_v.loc[adj_sig_up, y_col],
           s=90, color='#C62828', alpha=1.0, linewidths=0.5,
           edgecolors='#8B0000', marker='*', zorder=5,
           label=f'Adj. sig. up (adj p<=0.05, |log2FC|>=1; n={n_adj_up})')
ax.scatter(df_v.loc[adj_sig_down, 'log2fc'], df_v.loc[adj_sig_down, y_col],
           s=90, color='#1565C0', alpha=1.0, linewidths=0.5,
           edgecolors='#0D47A1', marker='*', zorder=5,
           label=f'Adj. sig. down (adj p<=0.05, |log2FC|>=1; n={n_adj_down})')
if len(df_annot):
    ax.scatter(df_annot['log2fc'], df_annot[y_col],
               s=70, color='#FF6F00', alpha=1.0, linewidths=0.9,
               edgecolors='#E65100', zorder=7,
               label=f'Genes of interest (n={len(df_annot)})')

ymax_plot  = df_v[y_col].quantile(0.9998)
y_axis_top = ymax_plot * 1.15
ax.axhline(-np.log10(NOM_THR), color='#9E9E9E', lw=1.0, ls='--', alpha=0.7,
           label=f'Nominal p = {NOM_THR}')
ax.axvline( FC_THR, color='#9E9E9E', lw=1.0, ls=':', alpha=0.7)
ax.axvline(-FC_THR, color='#9E9E9E', lw=1.0, ls=':', alpha=0.7)

ax.set_xlabel('log2 fold change (spaceflight / ground control)', fontsize=12)
ax.set_ylabel('-log10 (nominal p-value)', fontsize=12)
ax.set_title('OSD-562 Volcano Plot, ISS-T FLT 20 Weeks (n=4) vs ISS-T GC 20 Weeks (n=3)',
             fontsize=11.5, fontweight='bold', pad=14)
ax.set_xlim(df_v['log2fc'].quantile(0.001) - 2.0, df_v['log2fc'].quantile(0.999) + 1.5)
ax.set_ylim(-0.2, y_axis_top * 1.12)
ax.legend(fontsize=9.5, frameon=True, framealpha=0.92, edgecolor='#cccccc',
          loc='upper left', markerscale=1.3, labelspacing=0.6,
          handlelength=1.8, borderpad=0.8)

# ──────────────────────────────────────────────────────────────────────────────
# 5. Gene labels with adjustText
# ──────────────────────────────────────────────────────────────────────────────
label_df = pd.concat([df_annot, df_adj], ignore_index=True) if len(df_adj) else df_annot
texts, dot_coords = [], []
for _, row in label_df.iterrows():
    txt = ax.text(row['log2fc'], row[y_col], row['display_label'],
                  fontsize=8, ha='center', va='center',
                  color='#111111', fontweight='bold', zorder=9, clip_on=False,
                  bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.85))
    texts.append(txt)
    dot_coords.append((row['log2fc'], row[y_col]))

if texts:
    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax,
                    x=label_df['log2fc'].values, y=label_df[y_col].values,
                    force_text=(0.6, 1.2), force_static=(0.2, 0.4),
                    force_pull=(0.01, 0.02), force_explode=(3.0, 5.0),
                    expand=(1.5, 2.2), max_move=(400, 400),
                    iter_lim=1000, min_arrow_len=99999)
        print("adjustText applied successfully.")
    except ImportError:
        print("adjustText not installed — labels placed at point positions.")
        for txt in texts:
            x0, y0 = txt.get_position()
            txt.set_position((x0, y0 + 0.3))
    for txt, (px, py) in zip(texts, dot_coords):
        lx, ly = txt.get_position()
        if np.hypot(lx - px, ly - py) > 0.08:
            ax.annotate('', xy=(px, py), xytext=(lx, ly),
                        xycoords='data', textcoords='data', zorder=6,
                        arrowprops=dict(arrowstyle='-', color='#555555',
                                        lw=0.8, shrinkA=4, shrinkB=4))

# ──────────────────────────────────────────────────────────────────────────────
# 6. Save
# ──────────────────────────────────────────────────────────────────────────────
out_path = 'OSD562_Volcano_custom_genes.png'
plt.savefig(out_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"\nSaved → {out_path}")

print("\nAnnotation gene summary:")
summary = label_df[['display_label', 'ensembl', 'log2fc', 'pval', 'adjpval']].copy()
summary.columns = ['Gene', 'ENSEMBL', 'log2FC', 'p-value', 'adj. p-value']
summary = summary.sort_values('p-value')
summary.to_csv('OSD562_Volcano_custom_genes_table.csv', index=False)
print(summary.to_string(index=False))
