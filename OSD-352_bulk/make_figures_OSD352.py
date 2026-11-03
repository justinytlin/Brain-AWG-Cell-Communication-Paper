"""
OSD-352 (GLDS-352, RR-3) whole-brain bulk RNA-seq: PCA, volcano plots and the
significant-gene spreadsheet.  FLT (n = 3: F1, F2, F7) vs GC (n = 3: G7, G8, G9).

Conventions (harmonized with OSD-562 / OSD-564):
  * contrast = (Space Flight)v(Ground Control): positive log2FC = up in spaceflight
  * volcano legend counts are computed over ALL genes with non-missing statistics
    (22,264 Ensembl IDs), i.e. before any symbol-based filtering
  * when one gene symbol maps to several Ensembl IDs, the row with the highest
    mean normalized count is kept (applied before annotation / export)
"""
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'Arial',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300,
})

# ──────────────────────────────────────────────────────────────────
# 1.  Load data
# ──────────────────────────────────────────────────────────────────
print("Loading normalized counts …")
counts_raw = pd.read_csv('GLDS-352_rna_seq_Normalized_Counts.csv', index_col=0)

# log2(x + 1) transform as a substitute for VST
counts = np.log2(counts_raw + 1)

print("Loading differential expression data …")
de = pd.read_csv('GLDS-352_rna_seq_differential_expression.csv',
                 index_col=0, low_memory=False)
de.index.name = 'ENSEMBL'

print(f"Counts matrix: {counts.shape[0]:,} genes × {counts.shape[1]} samples")
print(f"DE table:      {de.shape[0]:,} genes × {de.shape[1]} columns")

# ──────────────────────────────────────────────────────────────────
# 2.  Group definitions
# ──────────────────────────────────────────────────────────────────
GROUP_INFO = {
    'Spaceflight (FLT)': {
        'samples': ['RR3_BRN_FLT_F1', 'RR3_BRN_FLT_F2', 'RR3_BRN_FLT_F7'],
        'color':   '#D32F2F',
    },
    'Ground control (GC)': {
        'samples': ['RR3_BRN_GC_G7', 'RR3_BRN_GC_G8', 'RR3_BRN_GC_G9'],
        'color':   '#1565C0',
    },
}

# ──────────────────────────────────────────────────────────────────
# 3.  PCA  ── top 2,000 most-variable genes, z-scored
# ──────────────────────────────────────────────────────────────────
print("\nRunning PCA …")
TOP_N = 2000
gene_var = counts.var(axis=1)
top_genes = gene_var.nlargest(TOP_N).index
X = counts.loc[top_genes].T          # samples × genes

X_scaled = StandardScaler().fit_transform(X)

pca = PCA(n_components=6, random_state=42)
pca_coords = pca.fit_transform(X_scaled)
pca_df = pd.DataFrame(pca_coords, index=X.index,
                      columns=[f'PC{i+1}' for i in range(6)])

var_exp = pca.explained_variance_ratio_ * 100
print(f"PC1: {var_exp[0]:.1f}%   PC2: {var_exp[1]:.1f}%   PC3: {var_exp[2]:.1f}%")

fig, ax = plt.subplots(figsize=(7, 5.8))

for gname, ginfo in GROUP_INFO.items():
    samps = [s for s in ginfo['samples'] if s in pca_df.index]
    ax.scatter(pca_df.loc[samps, 'PC1'], pca_df.loc[samps, 'PC2'],
               s=160, color=ginfo['color'], marker='o',
               edgecolors='black', linewidths=1.4, alpha=1.0, zorder=5,
               label=gname)
    for s in samps:
        ax.annotate(s.split('_')[-1], (pca_df.loc[s, 'PC1'], pca_df.loc[s, 'PC2']),
                    textcoords='offset points', xytext=(8, 6), fontsize=8.5,
                    color='#333333', zorder=6)

ax.set_xlabel(f'PC1 ({var_exp[0]:.1f}% variance)', fontsize=12)
ax.set_ylabel(f'PC2 ({var_exp[1]:.1f}% variance)', fontsize=12)
fig.suptitle('OSD-352 PCA, whole-brain bulk, FLT (n=3) vs GC (n=3)',
             fontsize=12, fontweight='bold', y=0.955)
ax.set_title('log2(normalized counts + 1), top 2,000 variable genes, z-scored',
             fontsize=10, color='#444444', pad=8)

ax.legend(fontsize=10, frameon=True, framealpha=0.9,
          edgecolor='#cccccc', loc='best', markerscale=1.1,
          labelspacing=1.0, handlelength=1.5, borderpad=0.8)

ax.axhline(0, color='#cccccc', lw=0.8, ls='--', zorder=0)
ax.axvline(0, color='#cccccc', lw=0.8, ls='--', zorder=0)

plt.tight_layout(rect=(0, 0, 1, 0.93))
pca_path = 'OSD352_PCA_plot.png'
plt.savefig(pca_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved → {pca_path}")

# ──────────────────────────────────────────────────────────────────
# 4.  DE table: all tested genes (for legend counts) and the
#     symbol-deduplicated table (for annotation / export)
# ──────────────────────────────────────────────────────────────────
print("\nPreparing DE tables …")

# Convention: Log2fc_(Space Flight)v(Ground Control) = log2(Spaceflight / Ground Control)
# Positive = up in spaceflight, negative = down in spaceflight
COMP_TAG = '(Space Flight)v(Ground Control)'

log2fc_col = f'Log2fc_{COMP_TAG}'
pval_col   = f'P.value_{COMP_TAG}'
adjp_col   = f'Adj.p.value_{COMP_TAG}'

for c in [log2fc_col, pval_col, adjp_col, 'SYMBOL']:
    if c not in de.columns:
        raise KeyError(f"column not found in DE table: {c}")

df_all = de[['SYMBOL', 'GENENAME', 'ENTREZID', log2fc_col, pval_col, adjp_col]].copy()
df_all.columns = ['symbol', 'genename', 'entrezid', 'log2fc', 'pval', 'adjpval']
for c in ['log2fc', 'pval', 'adjpval']:
    df_all[c] = pd.to_numeric(df_all[c], errors='coerce')

# (a) ALL TESTED GENES: non-missing statistics, symbol may be missing
df_all = df_all.dropna(subset=['log2fc', 'pval', 'adjpval'])
df_all['-log10p']    = -np.log10(df_all['pval'].clip(lower=1e-300))
df_all['-log10adjp'] = -np.log10(df_all['adjpval'].clip(lower=1e-300))
n_tested = len(df_all)
n_no_symbol = int(df_all['symbol'].isna().sum())
print(f"Tested genes (non-missing statistics): {n_tested:,}  "
      f"({n_no_symbol:,} without a gene symbol)")

# (b) SYMBOL-DEDUPLICATED table: one row per symbol, keeping the Ensembl ID with
#     the highest mean normalized count (e.g. Gcat: ENSMUSG00000006378 is kept,
#     ENSMUSG00000116378 — counts of 1 in five of six samples — is dropped)
mean_counts = counts_raw.mean(axis=1)
df_sym = df_all.dropna(subset=['symbol']).copy()
df_sym = df_sym[df_sym['symbol'].astype(str).str.strip() != '']
df_sym['mean_norm_count'] = mean_counts.reindex(df_sym.index).values

dup_mask = df_sym['symbol'].duplicated(keep=False)
dup_symbols = sorted(df_sym.loc[dup_mask, 'symbol'].unique())
n_rows_before = len(df_sym)
df_sym = (df_sym.sort_values(['symbol', 'mean_norm_count'], ascending=[True, False])
                .drop_duplicates(subset='symbol', keep='first'))
df_sym = df_sym.loc[df_all.index.intersection(df_sym.index)]   # restore original order
print(f"Symbol deduplication: {len(dup_symbols)} symbols mapped to >1 Ensembl ID; "
      f"{n_rows_before - len(df_sym)} rows dropped (highest mean normalized count kept).")
print("  Affected symbols:", ', '.join(dup_symbols))
if 'Gcat' in dup_symbols:
    kept = df_sym.index[df_sym['symbol'] == 'Gcat'].tolist()
    print(f"  Gcat → kept {kept}")

# Significance thresholds
FC_THR  = 1.0   # |log2FC| >= 1
ADJ_THR = 0.05  # adjusted p-value
NOM_THR = 0.05  # nominal p-value

def classify_masks(d):
    adj_up   = (d['log2fc'] >=  FC_THR) & (d['adjpval'] <= ADJ_THR)
    adj_down = (d['log2fc'] <= -FC_THR) & (d['adjpval'] <= ADJ_THR)
    nom_up   = (d['log2fc'] >=  FC_THR) & (d['pval'] <= NOM_THR) & ~adj_up
    nom_down = (d['log2fc'] <= -FC_THR) & (d['pval'] <= NOM_THR) & ~adj_down
    not_sig  = ~(adj_up | adj_down | nom_up | nom_down)
    return adj_up, adj_down, nom_up, nom_down, not_sig

# Legend counts: over ALL tested genes (not after dropping symbol-less rows)
adj_sig_up, adj_sig_down, nom_sig_up, nom_sig_down, not_sig = classify_masks(df_all)
n_adj_up, n_adj_down = int(adj_sig_up.sum()), int(adj_sig_down.sum())
n_nom_up, n_nom_down = int(nom_sig_up.sum()), int(nom_sig_down.sum())
print(f"\nLegend counts (n = {n_tested:,} tested genes):")
print(f"  Adj.p<=0.05 & |log2FC|>=1 — UP: {n_adj_up}  DOWN: {n_adj_down}")
print(f"  Nom.p<=0.05 & |log2FC|>=1 (not adj.) — UP: {n_nom_up}  DOWN: {n_nom_down}")

# Adjusted-significant genes (any |log2FC|) to be labelled by name on every volcano
adj_any = df_sym[df_sym['adjpval'] <= ADJ_THR].copy()
print("  Adjusted-significant genes (adj.p<=0.05): "
      + ', '.join(f"{r.symbol} ({r.log2fc:+.2f}, adj p {r.adjpval:.3f})" for r in adj_any.itertuples()))

# ──────────────────────────────────────────────────────────────────
# 5.  Volcano helpers
# ──────────────────────────────────────────────────────────────────
Y_COL = '-log10p'   # y-axis uses the nominal p-value (more informative at n = 3 vs 3)
X_LABEL = 'log2 fold change (spaceflight / ground control)'
Y_LABEL = '-log10 (nominal p-value)'
TITLE_MAIN = 'OSD-352 Volcano Plot, RR-3 whole-brain bulk RNA-seq, FLT (n=3) vs GC (n=3), ~18 weeks'
LEGEND_NOTE = f'n counted over all tested genes ({n_tested:,} with non-missing statistics)'

def draw_volcano_layers(ax):
    """Background + significance layers; legend counts are over all tested genes."""
    ax.scatter(df_all.loc[not_sig, 'log2fc'], df_all.loc[not_sig, Y_COL],
               s=7, color='#BDBDBD', alpha=0.40, linewidths=0, rasterized=True,
               label='Not significant')
    ax.scatter(df_all.loc[nom_sig_up, 'log2fc'], df_all.loc[nom_sig_up, Y_COL],
               s=18, color='#EF9A9A', alpha=0.80, linewidths=0.3,
               edgecolors='#C62828', rasterized=True,
               label=f'Nom. sig. up in spaceflight (p<=0.05, |log2FC|>=1; n={n_nom_up})')
    ax.scatter(df_all.loc[nom_sig_down, 'log2fc'], df_all.loc[nom_sig_down, Y_COL],
               s=18, color='#90CAF9', alpha=0.80, linewidths=0.3,
               edgecolors='#1565C0', rasterized=True,
               label=f'Nom. sig. down in spaceflight (p<=0.05, |log2FC|>=1; n={n_nom_down})')
    ax.scatter(df_all.loc[adj_sig_up, 'log2fc'], df_all.loc[adj_sig_up, Y_COL],
               s=110, color='#C62828', alpha=1.0, linewidths=0.5,
               edgecolors='#8B0000', marker='*', zorder=5,
               label=f'Adj. sig. up in spaceflight (adj p<=0.05, |log2FC|>=1; n={n_adj_up})')
    ax.scatter(df_all.loc[adj_sig_down, 'log2fc'], df_all.loc[adj_sig_down, Y_COL],
               s=110, color='#1565C0', alpha=1.0, linewidths=0.5,
               edgecolors='#0D47A1', marker='*', zorder=5,
               label=f'Adj. sig. down in spaceflight (adj p<=0.05, |log2FC|>=1; n={n_adj_down})')

    ax.axhline(-np.log10(NOM_THR), color='#9E9E9E', lw=1.0, ls='--', alpha=0.7,
               label=f'Nominal p = {NOM_THR}')
    ax.axvline( FC_THR, color='#9E9E9E', lw=1.0, ls=':', alpha=0.7)
    ax.axvline(-FC_THR, color='#9E9E9E', lw=1.0, ls=':', alpha=0.7)
    ax.set_xlabel(X_LABEL, fontsize=12)
    ax.set_ylabel(Y_LABEL, fontsize=12)

def label_points(ax, d, fontsize=8, connector=True):
    """Place gene-name labels (adjustText if available) with connector lines."""
    if len(d) == 0:
        return
    texts, coords = [], []
    for _, row in d.iterrows():
        txt = ax.text(row['log2fc'], row[Y_COL], row['label'],
                      fontsize=fontsize, ha='center', va='center',
                      color='#111111', fontweight='bold', zorder=8, clip_on=False,
                      bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.85))
        texts.append(txt)
        coords.append((row['log2fc'], row[Y_COL]))
    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax,
                    x=d['log2fc'].values, y=d[Y_COL].values,
                    force_text=(0.5, 1.0), force_static=(0.15, 0.3),
                    force_pull=(0.01, 0.02), force_explode=(3.0, 5.0),
                    expand=(1.4, 2.0), max_move=(300, 300), iter_lim=1000,
                    min_arrow_len=99999)
    except ImportError:
        for txt in texts:
            x0, y0 = txt.get_position()
            txt.set_position((x0, y0 + 0.25))
    if connector:
        for txt, (px, py) in zip(texts, coords):
            lx, ly = txt.get_position()
            if np.hypot(lx - px, ly - py) > 0.08:
                ax.annotate('', xy=(px, py), xytext=(lx, ly),
                            xycoords='data', textcoords='data', zorder=6,
                            arrowprops=dict(arrowstyle='-', color='#555555',
                                            lw=0.8, shrinkA=4, shrinkB=4))

def select_annotation(requested):
    """
    Match a curated list against the symbol-deduplicated table.  Requested entries
    may be pipe-joined aliases ('Ndufb4|Ndufb4c'); DE symbols may be pipe-joined too
    ('Atp5mk|Gm35857').  Both sides are split on '|' and matched case-insensitively.
    The display label is the first alias of the requested entry.
    """
    alias_to_label = {}
    for entry in requested:
        parts = [p.strip() for p in entry.split('|') if p.strip()]
        for p in parts:
            alias_to_label.setdefault(p.lower(), parts[0])

    rows = []
    for ens, row in df_sym.iterrows():
        for part in str(row['symbol']).split('|'):
            lab = alias_to_label.get(part.strip().lower())
            if lab is not None:
                r = row.copy()
                r['label'] = lab
                r['ENSEMBL'] = ens
                rows.append(r)
                break
    out = pd.DataFrame(rows)
    if len(out):
        out = out.sort_values('pval').drop_duplicates(subset='label', keep='first')
    wanted = {e.split('|')[0] for e in requested}
    missing = sorted(wanted - set(out['label'])) if len(out) else sorted(wanted)
    return out, missing

def add_adj_sig_labels(d):
    """Make sure every adjusted-significant gene (e.g. Tcf7l1) is labelled."""
    extra = adj_any.copy()
    extra['label'] = extra['symbol']
    extra['ENSEMBL'] = extra.index
    if len(d) == 0:
        return extra.reset_index(drop=True)
    extra = extra[~extra['label'].isin(d['label'])]
    return pd.concat([d, extra], ignore_index=True)

def y_top(pad_top=1.15):
    """Upper y-limit: quantile-based, but always including the most significant gene."""
    return max(df_all[Y_COL].quantile(0.9998) * pad_top, df_all[Y_COL].max() * 1.10)

def xy_limits(ax, pad_x=0.4, pad_top=1.15):
    ax.set_xlim(df_all['log2fc'].quantile(0.001) - pad_x,
                df_all['log2fc'].quantile(0.999) + pad_x)
    ax.set_ylim(-0.2, y_top(pad_top))

def finish(ax, title_second_line=None, legend_loc='upper right', legend_fs=7.5):
    title = TITLE_MAIN if title_second_line is None else f'{TITLE_MAIN}\n{title_second_line}'
    ax.set_title(title, fontsize=11, fontweight='bold', pad=14)
    leg = ax.legend(fontsize=legend_fs, frameon=True, framealpha=0.92,
                    edgecolor='#cccccc', loc=legend_loc, ncol=1,
                    title=LEGEND_NOTE, title_fontsize=legend_fs - 0.5)
    leg._legend_box.align = 'left'

# ──────────────────────────────────────────────────────────────────
# 6a. Volcano plot (adjusted-significant gene labelled)
# ──────────────────────────────────────────────────────────────────
print("\nPreparing volcano plot …")
fig, ax = plt.subplots(figsize=(9, 7.5))
draw_volcano_layers(ax)
xy_limits(ax)
label_points(ax, add_adj_sig_labels(pd.DataFrame()), fontsize=9)
finish(ax)
volcano_path = 'OSD352_Volcano_plot.png'
plt.savefig(volcano_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved → {volcano_path}")

# ──────────────────────────────────────────────────────────────────
# 6b. Volcano plot with the OSD-352-specific curated gene annotations
# ──────────────────────────────────────────────────────────────────
print("\nPreparing annotated volcano plot (curated OSD-352 genes) …")

annotation_genes = [
    "Rac2", "Myo5c", "Psd4", "Epha2", "Itga2", "Rab43", "Ap1m2",
    "Cd209a", "Cd209g", "Cd300lf", "Klk8", "Mapk12", "Adora2a", "Nos2", "Agt",
    "Cpt1b", "Gcat", "Hpdl", "Bbox1", "Osgin1", "Pask", "Aqp9", "Sdc4",
]

df_annot, missing = select_annotation(annotation_genes)
df_annot = add_adj_sig_labels(df_annot)
print(f"Found {len(df_annot)} annotation genes (incl. adjusted-significant) in DE data; "
      f"missing: {missing if missing else 'none'}")
for r in df_annot.sort_values('pval').itertuples():
    print(f"  {r.label:<8s} {r.ENSEMBL}  log2FC={r.log2fc:+.3f}  p={r.pval:.4g}  adj p={r.adjpval:.3g}")

fig, ax = plt.subplots(figsize=(14, 11))
draw_volcano_layers(ax)
ax.scatter(df_annot['log2fc'], df_annot[Y_COL],
           s=60, color='#FF6F00', alpha=1.0, linewidths=0.8,
           edgecolors='#E65100', zorder=6,
           label=f'Annotated genes ({len(df_annot)})')
ax.set_xlim(df_all['log2fc'].quantile(0.001) - 2.5, df_all['log2fc'].quantile(0.999) + 1.5)
ax.set_ylim(-0.2, y_top(1.15 * 1.12))
label_points(ax, df_annot, fontsize=8)
finish(ax, title_second_line='Curated OSD-352 gene annotations',
       legend_loc='upper left', legend_fs=10)
annot_volcano_path = 'OSD352_Volcano_plot_annotated.png'
plt.savefig(annot_volcano_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved → {annot_volcano_path}")

# ──────────────────────────────────────────────────────────────────
# 6c. Volcano plot with mitochondrial / cytoskeletal gene annotations
# ──────────────────────────────────────────────────────────────────
print("\nPreparing mitochondrial/cytoskeletal annotated volcano plot …")

annotation_genes_v2 = [
    "Adam8", "Apoc1", "Aqp1", "Arpc2", "Arap2", "Baiap2l2", "Bcat1",
    "C4a", "Casr", "Ect2", "Epcam", "Folr1", "Gas2", "Gas2l3", "Ghsr",
    "Itga2", "Lpl", "Myh7b", "Rac2", "Rhof", "S100a4", "S100a8",
    "Smtnl2", "Timp1", "Tmsb4x", "Tnni2", "Wdr1", "Wnt7a", "Zfas1",
    "Mrps25", "Chchd1", "Rmrp", "Ndufb4|Ndufb4c|Ndufb4b", "Ndufaf6",
    "Atp5mk|Gm35857", "Pet100", "COX2|mt-Co2",
]

df_annot_v2, missing_v2 = select_annotation(annotation_genes_v2)
df_annot_v2 = add_adj_sig_labels(df_annot_v2)
print(f"Found {len(df_annot_v2)} annotation genes (incl. adjusted-significant) in DE data; "
      f"missing: {missing_v2 if missing_v2 else 'none'}")

fig, ax = plt.subplots(figsize=(14, 11))
draw_volcano_layers(ax)
ax.scatter(df_annot_v2['log2fc'], df_annot_v2[Y_COL],
           s=60, color='#FF6F00', alpha=1.0, linewidths=0.8,
           edgecolors='#E65100', zorder=6,
           label=f'Annotated genes ({len(df_annot_v2)})')
ax.set_xlim(df_all['log2fc'].quantile(0.001) - 2.5, df_all['log2fc'].quantile(0.999) + 1.5)
ax.set_ylim(-0.2, y_top(1.15 * 1.12))
label_points(ax, df_annot_v2, fontsize=8)
finish(ax, title_second_line='Mitochondrial / cytoskeletal gene annotations',
       legend_loc='upper left', legend_fs=10)
annot_volcano_v2_path = 'OSD352_Volcano_plot_annotated_mito_cyto.png'
plt.savefig(annot_volcano_v2_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved → {annot_volcano_v2_path}")

# ──────────────────────────────────────────────────────────────────
# 7.  Spreadsheet of significant genes (symbol-deduplicated table)
# ──────────────────────────────────────────────────────────────────
print("\nExporting significant genes spreadsheet …")

df_export = df_sym.copy()
df_export = df_export.rename(columns={
    'symbol': 'SYMBOL', 'genename': 'GENENAME', 'entrezid': 'ENTREZID',
    'log2fc': 'log2FC_SpaceflightVsGC', 'pval': 'P_value', 'adjpval': 'Adj_P_value',
    'mean_norm_count': 'Mean_normalized_count',
})

def classify(row):
    fc, p, adjp = row['log2FC_SpaceflightVsGC'], row['P_value'], row['Adj_P_value']
    if   adjp <= ADJ_THR and fc >=  FC_THR: return 'Adj.sig UP'
    elif adjp <= ADJ_THR and fc <= -FC_THR: return 'Adj.sig DOWN'
    elif p    <= NOM_THR and fc >=  FC_THR: return 'Nom.sig UP'
    elif p    <= NOM_THR and fc <= -FC_THR: return 'Nom.sig DOWN'
    else: return 'Not significant'

df_export['Significance'] = df_export.apply(classify, axis=1)
df_export['Direction']    = df_export['log2FC_SpaceflightVsGC'].apply(
    lambda x: 'UP in spaceflight' if x > 0 else 'DOWN in spaceflight')

EXPORT_COLS = ['ENSEMBL', 'SYMBOL', 'GENENAME', 'ENTREZID', 'Mean_normalized_count',
               'log2FC_SpaceflightVsGC', 'P_value', 'Adj_P_value', 'Significance', 'Direction']

adj_sig = df_export[df_export['Significance'].str.startswith('Adj.sig')].sort_values('Adj_P_value')
nom_sig = df_export[df_export['Significance'].str.startswith('Nom.sig')].sort_values('P_value')
all_sig = pd.concat([adj_sig, nom_sig])

print(f"  Adj.sig:  {len(adj_sig)} genes  ({(adj_sig['Direction']=='UP in spaceflight').sum()} UP, "
      f"{(adj_sig['Direction']=='DOWN in spaceflight').sum()} DOWN)")
print(f"  Nom.sig:  {len(nom_sig)} genes  ({(nom_sig['Direction']=='UP in spaceflight').sum()} UP, "
      f"{(nom_sig['Direction']=='DOWN in spaceflight').sum()} DOWN)")

curated = df_export[df_export['SYMBOL'].isin(
    ['Rac2', 'Itga2', 'Cpt1b', 'Hpdl', 'Bbox1', 'Pask', 'Sdc4', 'Tcf7l1', 'Gcat'])].reset_index()

dedup_report = pd.DataFrame({
    'Symbol': dup_symbols,
    'Ensembl IDs (all)': [', '.join(df_all.index[df_all['symbol'] == s]) for s in dup_symbols],
    'Ensembl ID kept (highest mean normalized count)':
        [df_sym.index[df_sym['symbol'] == s][0] for s in dup_symbols],
})

readme = pd.DataFrame({'Note': [
    'OSD-352 (GLDS-352, RR-3) whole-brain bulk RNA-seq, FLT (n = 3) vs GC (n = 3).',
    f'Contrast: {COMP_TAG}; positive log2FC = up in spaceflight.',
    f'Tested genes with non-missing statistics: {n_tested:,} (volcano legend counts use this set).',
    f'Rows without a gene symbol ({n_no_symbol:,}) are excluded from the gene sheets.',
    f'Symbols mapping to several Ensembl IDs ({len(dup_symbols)}): the Ensembl ID with the highest '
    'mean normalized count is kept (see sheet Deduplicated symbols).',
    'Significance: Adj.sig = adj p <= 0.05 and |log2FC| >= 1; Nom.sig = nominal p <= 0.05 and |log2FC| >= 1 (not Adj.sig).',
]})

xl_path = 'OSD352_Significant_Genes.xlsx'
with pd.ExcelWriter(xl_path, engine='openpyxl') as writer:
    readme.to_excel(writer, sheet_name='README', index=False)
    all_sig.reset_index()[EXPORT_COLS].to_excel(writer, sheet_name='All Significant', index=False)
    adj_sig.reset_index()[EXPORT_COLS].to_excel(writer, sheet_name='Adj.p Significant (0.05)', index=False)
    nom_sig.reset_index()[EXPORT_COLS].to_excel(writer, sheet_name='Nominal p Significant (0.05)', index=False)
    curated[EXPORT_COLS].to_excel(writer, sheet_name='Curated genes', index=False)
    dedup_report.to_excel(writer, sheet_name='Deduplicated symbols', index=False)
    for ws in writer.sheets.values():
        for col_cells in ws.columns:
            width = max(len(str(c.value)) if c.value is not None else 0 for c in col_cells)
            ws.column_dimensions[col_cells[0].column_letter].width = min(width + 2, 70)

print(f"Saved → {xl_path}")
print("\nCurated genes after deduplication:")
print(curated[['ENSEMBL', 'SYMBOL', 'log2FC_SpaceflightVsGC', 'P_value', 'Adj_P_value', 'Direction']]
      .to_string(index=False))
print("\nDone.")
