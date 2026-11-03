"""
make_figures.py  (OSD-562, re-run DESeq2 4 vs 3)
------------------------------------------------
PCA (7-sample VST from rerun_deseq2_OSD562.py), volcano plots and the
significant-gene spreadsheet for
    ISS-T FLT 20 weeks (n = 4)  vs  ISS-T GC 20 weeks (n = 3)
using OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3.csv (not the GeneLab 3 v 6 table).
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'Arial',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300,
})

DE_FILE  = 'OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3.csv'
VST_FILE = 'OSD562_VST_counts_ISS-T_YNG.csv'
COMP_TAG = '(Space Flight)v(Ground Control)'
TITLE_N  = 'ISS-T FLT 20 Weeks (n=4) vs ISS-T GC 20 Weeks (n=3)'

# ──────────────────────────────────────────────────────────────────
# 1.  Load data
# ──────────────────────────────────────────────────────────────────
print("Loading VST counts (7-sample re-run) …")
vst = pd.read_csv(VST_FILE, index_col=0)

print("Loading differential expression data (re-run) …")
de = pd.read_csv(DE_FILE, index_col=0, low_memory=False)

print(f"VST matrix: {vst.shape[0]:,} genes × {vst.shape[1]} samples")
print(f"DE table:   {de.shape[0]:,} genes × {de.shape[1]} columns")

# ──────────────────────────────────────────────────────────────────
# 2.  Group definitions
# ──────────────────────────────────────────────────────────────────
GROUP_INFO = {
    'ISS-T FLT': {
        'samples': ['RRRM2_HPC_FLT_ISS-T_YNG_FY3',
                    'RRRM2_HPC_FLT_ISS-T_YNG_FY4',
                    'RRRM2_HPC_FLT_ISS-T_YNG_FY6',
                    'RRRM2_HPC_FLT_ISS-T_YNG_FY10'],
        'color':   '#D32F2F',
    },
    'ISS-T GC': {
        'samples': ['RRRM2_HPC_GC_ISS-T_YNG_GY5',
                    'RRRM2_HPC_GC_ISS-T_YNG_GY8',
                    'RRRM2_HPC_GC_ISS-T_YNG_GY10'],
        'color':   '#1565C0',
    },
}
SAMPLES_7 = GROUP_INFO['ISS-T FLT']['samples'] + GROUP_INFO['ISS-T GC']['samples']
missing = [s for s in SAMPLES_7 if s not in vst.columns]
assert not missing, f"VST is missing samples: {missing}"
vst = vst[SAMPLES_7]            # subset to the 7 samples BEFORE variance ranking

# ──────────────────────────────────────────────────────────────────
# 3.  PCA  ── top 2,000 most-variable genes (7 samples only)
# ──────────────────────────────────────────────────────────────────
print("\nRunning PCA …")
TOP_N = 2000
gene_var = vst.var(axis=1)
top_genes = gene_var.nlargest(TOP_N).index
X = vst.loc[top_genes].T          # samples × genes

X_scaled = StandardScaler().fit_transform(X)
n_comp = min(10, X.shape[0] - 1)
pca = PCA(n_components=n_comp, random_state=42)
pca_coords = pca.fit_transform(X_scaled)
pca_df = pd.DataFrame(pca_coords, index=X.index,
                      columns=[f'PC{i+1}' for i in range(n_comp)])

var_exp = pca.explained_variance_ratio_ * 100
print(f"PC1: {var_exp[0]:.1f}%   PC2: {var_exp[1]:.1f}%   PC3: {var_exp[2]:.1f}%")
pca_df.assign(group=['ISS-T FLT'] * 4 + ['ISS-T GC'] * 3).to_csv('OSD562_PCA_coordinates.csv')

fig, ax = plt.subplots(figsize=(7, 5.5))
for gname, ginfo in GROUP_INFO.items():
    samps = [s for s in ginfo['samples'] if s in pca_df.index]
    ax.scatter(pca_df.loc[samps, 'PC1'], pca_df.loc[samps, 'PC2'],
               s=160, color=ginfo['color'], marker='o',
               edgecolors='black', linewidths=1.4, alpha=1.0,
               zorder=5, label=gname)
    for s in samps:
        ax.annotate(s.split('_')[-1], (pca_df.at[s, 'PC1'], pca_df.at[s, 'PC2']),
                    xytext=(6, 6), textcoords='offset points', fontsize=8, color='#444444')

ax.set_xlabel(f'PC1 ({var_exp[0]:.1f}% variance)', fontsize=12)
ax.set_ylabel(f'PC2 ({var_exp[1]:.1f}% variance)', fontsize=12)
ax.set_title(f'OSD-562 PCA Plot, {TITLE_N}',
             fontsize=11.5, fontweight='bold', pad=12)
ax.legend(fontsize=10, frameon=True, framealpha=0.9,
          edgecolor='#cccccc', loc='best', markerscale=1.3,
          labelspacing=1.2, handlelength=1.5, borderpad=0.8)
ax.axhline(0, color='#cccccc', lw=0.8, ls='--', zorder=0)
ax.axvline(0, color='#cccccc', lw=0.8, ls='--', zorder=0)
plt.tight_layout()
pca_path = 'OSD562_PCA_plot.png'
plt.savefig(pca_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved → {pca_path}")

# ──────────────────────────────────────────────────────────────────
# 4.  Volcano plot  ── ISS-T Young FLT vs GC
# ──────────────────────────────────────────────────────────────────
print("\nPreparing volcano plot …")
log2fc_col = f'Log2fc_{COMP_TAG}'
pval_col   = f'P.value_{COMP_TAG}'
adjp_col   = f'Adj.p.value_{COMP_TAG}'
symbol_col = 'SYMBOL'
for c in [log2fc_col, pval_col, adjp_col, symbol_col]:
    assert c in de.columns, f"column not found: {c}"

df_v = de[[symbol_col, log2fc_col, pval_col, adjp_col]].copy()
df_v.columns = ['symbol', 'log2fc', 'pval', 'adjpval']
for c in ['log2fc', 'pval', 'adjpval']:
    df_v[c] = pd.to_numeric(df_v[c], errors='coerce')
df_v = df_v.dropna(subset=['log2fc', 'pval', 'adjpval'])
df_v['-log10p'] = -np.log10(df_v['pval'].clip(lower=1e-300))

FC_THR, ADJ_THR, NOM_THR = 1.0, 0.05, 0.05
adj_sig_up   = (df_v['log2fc'] >=  FC_THR) & (df_v['adjpval'] <= ADJ_THR)
adj_sig_down = (df_v['log2fc'] <= -FC_THR) & (df_v['adjpval'] <= ADJ_THR)
nom_sig_up   = (df_v['log2fc'] >=  FC_THR) & (df_v['pval'] <= NOM_THR) & ~adj_sig_up
nom_sig_down = (df_v['log2fc'] <= -FC_THR) & (df_v['pval'] <= NOM_THR) & ~adj_sig_down
not_sig      = ~(adj_sig_up | adj_sig_down | nom_sig_up | nom_sig_down)

n_adj_up, n_adj_down = adj_sig_up.sum(), adj_sig_down.sum()
n_nom_up, n_nom_down = nom_sig_up.sum(), nom_sig_down.sum()
print(f"Adj.p≤0.05 & |FC|≥1 — UP: {n_adj_up}  DOWN: {n_adj_down}")
print(f"Nom.p≤0.05 & |FC|≥1 — UP: {n_nom_up}  DOWN: {n_nom_down}")

LEG = {
    'nom_up':   f'Nom. sig. up in spaceflight (p<=0.05, |log2FC|>=1; n={n_nom_up})',
    'nom_down': f'Nom. sig. down in spaceflight (p<=0.05, |log2FC|>=1; n={n_nom_down})',
    'adj_up':   f'Adj. sig. up (adj p<=0.05, |log2FC|>=1; n={n_adj_up})',
    'adj_down': f'Adj. sig. down (adj p<=0.05, |log2FC|>=1; n={n_adj_down})',
}
XLAB = 'log2 fold change (spaceflight / ground control)'
YLAB = '-log10 (nominal p-value)'


def draw_volcano_base(ax, y_col='-log10p'):
    ax.scatter(df_v.loc[not_sig, 'log2fc'], df_v.loc[not_sig, y_col],
               s=7, color='#BDBDBD', alpha=0.40, linewidths=0, rasterized=True,
               label='Not significant')
    ax.scatter(df_v.loc[nom_sig_up, 'log2fc'], df_v.loc[nom_sig_up, y_col],
               s=18, color='#EF9A9A', alpha=0.80, linewidths=0.3,
               edgecolors='#C62828', rasterized=True, label=LEG['nom_up'])
    ax.scatter(df_v.loc[nom_sig_down, 'log2fc'], df_v.loc[nom_sig_down, y_col],
               s=18, color='#90CAF9', alpha=0.80, linewidths=0.3,
               edgecolors='#1565C0', rasterized=True, label=LEG['nom_down'])
    ax.scatter(df_v.loc[adj_sig_up, 'log2fc'], df_v.loc[adj_sig_up, y_col],
               s=90, color='#C62828', alpha=1.0, linewidths=0.5,
               edgecolors='#8B0000', marker='*', zorder=5, label=LEG['adj_up'])
    ax.scatter(df_v.loc[adj_sig_down, 'log2fc'], df_v.loc[adj_sig_down, y_col],
               s=90, color='#1565C0', alpha=1.0, linewidths=0.5,
               edgecolors='#0D47A1', marker='*', zorder=5, label=LEG['adj_down'])
    ax.axhline(-np.log10(NOM_THR), color='#9E9E9E', lw=1.0, ls='--', alpha=0.7,
               label=f'Nominal p = {NOM_THR}')
    ax.axvline( FC_THR, color='#9E9E9E', lw=1.0, ls=':', alpha=0.7)
    ax.axvline(-FC_THR, color='#9E9E9E', lw=1.0, ls=':', alpha=0.7)
    ax.set_xlabel(XLAB, fontsize=12)
    ax.set_ylabel(YLAB, fontsize=12)


y_col = '-log10p'
y_axis_top = df_v[y_col].quantile(0.9998) * 1.15

fig, ax = plt.subplots(figsize=(9, 7.5))
draw_volcano_base(ax)
ax.set_title(f'OSD-562 Volcano Plot, {TITLE_N}', fontsize=11.5, fontweight='bold', pad=14)
ax.set_xlim(df_v['log2fc'].quantile(0.001) - 0.4, df_v['log2fc'].quantile(0.999) + 0.4)
ax.set_ylim(-0.2, y_axis_top * 1.25)
ax.legend(fontsize=7.5, frameon=True, framealpha=0.92, edgecolor='#cccccc',
          loc='upper left', ncol=1)
volcano_path = 'OSD562_Volcano_plot.png'
plt.savefig(volcano_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved → {volcano_path}")

# ──────────────────────────────────────────────────────────────────
# 4b. Volcano plot with custom gene annotations
#     (all adjusted-significant genes are labelled as well)
# ──────────────────────────────────────────────────────────────────
print("\nPreparing annotated volcano plot …")
annotation_genes = [
    "Sh3bp1", "Stard13", "Ehbp1l1", "Filip1", "Myo1f", "Arhgap28",
    "Sdc4", "Peg10", "Gga3", "Arfgef1", "Doc2b", "Fam222a",
    "Adnp", "Nos1ap", "Abca2", "Pycard", "Ripk1", "Pla2g4a", "Hpgds",
    "Dhtkd1", "Mrs2", "Slc25a21", "Slc25a34", "Cpt1b", "Alas2", "Pmaip1", "Pprc1",
    "Rhof", "Ghsr", "Vps4a", "Ist1", "Rac2", "Baiap2l2",
]
annot_mask = df_v['symbol'].str.lower().isin([g.lower() for g in annotation_genes])
df_annot = df_v[annot_mask].copy()
df_adj   = df_v[adj_sig_up | adj_sig_down].copy()
print(f"Found {len(df_annot)} of {len(annotation_genes)} annotation genes; "
      f"{len(df_adj)} adjusted-significant genes to label")

fig, ax = plt.subplots(figsize=(14, 11))
draw_volcano_base(ax)
ax.scatter(df_annot['log2fc'], df_annot[y_col],
           s=60, color='#FF6F00', alpha=1.0, linewidths=0.8,
           edgecolors='#E65100', zorder=6,
           label=f'Annotated genes (n={len(df_annot)})')
ax.set_title(f'OSD-562 Volcano Plot – Custom Gene Annotations\n{TITLE_N}',
             fontsize=11.5, fontweight='bold', pad=14)
ax.set_xlim(df_v['log2fc'].quantile(0.001) - 2.5, df_v['log2fc'].quantile(0.999) + 1.5)
ax.set_ylim(-0.2, y_axis_top * 1.3)

label_df = pd.concat([df_annot, df_adj]).drop_duplicates()
texts, annot_coords = [], []
for _, row in label_df.iterrows():
    lab = row['symbol'] if isinstance(row['symbol'], str) else row.name
    txt = ax.text(row['log2fc'], row[y_col], lab,
                  fontsize=8, ha='center', va='center',
                  color='#111111', fontweight='bold', zorder=8, clip_on=False,
                  bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.85))
    texts.append(txt)
    annot_coords.append((row['log2fc'], row[y_col]))

if texts:
    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax,
                    x=label_df['log2fc'].values, y=label_df[y_col].values,
                    force_text=(0.5, 1.0), force_static=(0.15, 0.3),
                    force_pull=(0.01, 0.02), force_explode=(3.0, 5.0),
                    expand=(1.4, 2.0), max_move=(300, 300),
                    iter_lim=1000, min_arrow_len=99999)
    except ImportError:
        for txt in texts:
            x0, y0 = txt.get_position()
            txt.set_position((x0, y0 + 0.25))
    for txt, (px, py) in zip(texts, annot_coords):
        lx, ly = txt.get_position()
        if np.hypot(lx - px, ly - py) > 0.08:
            ax.annotate('', xy=(px, py), xytext=(lx, ly), xycoords='data',
                        textcoords='data', zorder=6,
                        arrowprops=dict(arrowstyle='-', color='#555555',
                                        lw=0.8, shrinkA=4, shrinkB=4))

ax.legend(fontsize=9.5, frameon=True, framealpha=0.92, edgecolor='#cccccc',
          loc='upper left', markerscale=1.4, labelspacing=0.6,
          handlelength=1.8, borderpad=0.8)
annot_volcano_path = 'OSD562_Volcano_plot_annotated.png'
plt.savefig(annot_volcano_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved → {annot_volcano_path}")

# ──────────────────────────────────────────────────────────────────
# 5.  Spreadsheet of significant genes
# ──────────────────────────────────────────────────────────────────
print("\nExporting significant genes spreadsheet …")
annot_cols = [c for c in ['SYMBOL', 'GENENAME', 'ENTREZID'] if c in de.columns]
df_export = de[annot_cols].copy()
df_export.index.name = 'ENSEMBL'
df_export['log2FC_SpaceflightVsGC'] = pd.to_numeric(de[log2fc_col], errors='coerce')
df_export['P_value']                 = pd.to_numeric(de[pval_col],   errors='coerce')
df_export['Adj_P_value']             = pd.to_numeric(de[adjp_col],   errors='coerce')
df_export['Group_Mean_Spaceflight']  = de['Group.Mean_(Space Flight)']
df_export['Group_Mean_GroundControl'] = de['Group.Mean_(Ground Control)']
df_export = df_export.dropna(subset=['log2FC_SpaceflightVsGC', 'P_value', 'Adj_P_value'])


def classify(row):
    fc, p, adjp = row['log2FC_SpaceflightVsGC'], row['P_value'], row['Adj_P_value']
    if   adjp <= 0.05 and fc >=  1.0: return 'Adj.sig UP'
    elif adjp <= 0.05 and fc <= -1.0: return 'Adj.sig DOWN'
    elif p    <= 0.05 and fc >=  1.0: return 'Nom.sig UP'
    elif p    <= 0.05 and fc <= -1.0: return 'Nom.sig DOWN'
    return 'Not significant'


df_export['Significance'] = df_export.apply(classify, axis=1)
df_export['Direction']    = df_export['log2FC_SpaceflightVsGC'].apply(
    lambda x: 'UP in Spaceflight' if x > 0 else 'DOWN in Spaceflight')

adj_sig = df_export[df_export['Significance'].str.startswith('Adj.sig')].sort_values('Adj_P_value')
nom_sig = df_export[df_export['Significance'].str.startswith('Nom.sig')].sort_values('P_value')
all_sig = pd.concat([adj_sig, nom_sig]).reset_index()
# nominal p <= 0.05 regardless of fold change (for reference)
nom_anyfc = df_export[df_export['P_value'] <= 0.05].sort_values('P_value').reset_index()

print(f"  Adj.sig:  {len(adj_sig)} genes  ({(adj_sig['Direction']=='UP in Spaceflight').sum()} UP, "
      f"{(adj_sig['Direction']=='DOWN in Spaceflight').sum()} DOWN)")
print(f"  Nom.sig:  {len(nom_sig)} genes  ({(nom_sig['Direction']=='UP in Spaceflight').sum()} UP, "
      f"{(nom_sig['Direction']=='DOWN in Spaceflight').sum()} DOWN)")
print(f"  Nom.p<=0.05 any FC: {len(nom_anyfc)} genes")

xl_path = 'OSD562_Significant_Genes.xlsx'
with pd.ExcelWriter(xl_path, engine='openpyxl') as writer:
    all_sig.to_excel(writer, sheet_name='All Significant', index=False)
    adj_sig.reset_index().to_excel(writer, sheet_name='Adj.p Significant (0.05)', index=False)
    nom_sig.reset_index().to_excel(writer, sheet_name='Nominal p Significant (0.05)', index=False)
    nom_anyfc.to_excel(writer, sheet_name='Nominal p (0.05) any FC', index=False)
    pd.DataFrame({
        'Metric': ['Contrast', 'DE table', 'Adj.sig UP', 'Adj.sig DOWN',
                   'Nom.sig UP (|log2FC|>=1)', 'Nom.sig DOWN (|log2FC|>=1)', 'Nom.p<=0.05 any FC'],
        'Value': ['ISS-T FLT 20wk (n=4) vs ISS-T GC 20wk (n=3), DESeq2 Wald, BH', DE_FILE,
                  int((adj_sig['Direction']=='UP in Spaceflight').sum()),
                  int((adj_sig['Direction']=='DOWN in Spaceflight').sum()),
                  int((nom_sig['Direction']=='UP in Spaceflight').sum()),
                  int((nom_sig['Direction']=='DOWN in Spaceflight').sum()),
                  len(nom_anyfc)],
    }).to_excel(writer, sheet_name='Summary', index=False)
print(f"Saved → {xl_path}")
print("\nDone.")
