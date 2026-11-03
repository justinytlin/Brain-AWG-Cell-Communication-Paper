"""
OSD-564 (GLDS-569, RR-10) volcano plots, PCA and significant-gene export.

Contrast: (Space Flight) v (Ground Control)  ->  log2FC > 0 = up in spaceflight.

Legend counts are computed over ALL genes with non-missing statistics (18,517 tested
genes), including rows whose SYMBOL is missing.  Curated-gene matching splits
pipe-joined SYMBOL entries ("Ndufb4|Ndufb4c|Ndufb4b") and matches the first token
first (then any token); the display label is the first token.  "COX2"
(ENSMUSG00000064354) is the mitochondrially encoded mt-Co2 and is displayed as such.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
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

# ──────────────────────────────────────────────────────────────────
# 1.  Load data
# ──────────────────────────────────────────────────────────────────
print("Loading VST counts ...")
vst = pd.read_csv('GLDS-569_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv', index_col=0)

print("Loading differential expression data ...")
de = pd.read_csv('GLDS-569_rna_seq_differential_expression_rRNArm_GLbulkRNAseq.csv',
                 index_col=0, low_memory=False)
de.index.name = 'ENSEMBL'

print(f"VST matrix: {vst.shape[0]:,} genes x {vst.shape[1]} samples")
print(f"DE table:   {de.shape[0]:,} genes x {de.shape[1]} columns")

# ──────────────────────────────────────────────────────────────────
# 2.  Group definitions
# ──────────────────────────────────────────────────────────────────
GROUP_INFO = {
    'RR-10 FLT': {
        'samples': ['RR10_HPC_FLT_WT_F1', 'RR10_HPC_FLT_WT_F3', 'RR10_HPC_FLT_WT_F5',
                    'RR10_HPC_FLT_WT_F7', 'RR10_HPC_FLT_WT_F9'],
        'color':   '#D32F2F',
    },
    'RR-10 GC': {
        'samples': ['RR10_HPC_GC_WT_G1', 'RR10_HPC_GC_WT_G3', 'RR10_HPC_GC_WT_G5',
                    'RR10_HPC_GC_WT_G7', 'RR10_HPC_GC_WT_G9'],
        'color':   '#1565C0',
    },
}

# ──────────────────────────────────────────────────────────────────
# 3.  PCA  -- top 2,000 most-variable genes, z-scored
# ──────────────────────────────────────────────────────────────────
print("\nRunning PCA ...")
TOP_N = 2000
gene_var = vst.var(axis=1)
top_genes = gene_var.nlargest(TOP_N).index
X = vst.loc[top_genes].T          # samples x genes

X_scaled = StandardScaler().fit_transform(X)

pca = PCA(n_components=10, random_state=42)
pca_coords = pca.fit_transform(X_scaled)
pca_df = pd.DataFrame(pca_coords, index=X.index,
                      columns=[f'PC{i+1}' for i in range(10)])

var_exp = pca.explained_variance_ratio_ * 100
print(f"PC1: {var_exp[0]:.1f}%   PC2: {var_exp[1]:.1f}%   PC3: {var_exp[2]:.1f}%")

fig, ax = plt.subplots(figsize=(7, 5.5))

for gname, ginfo in GROUP_INFO.items():
    samps = [s for s in ginfo['samples'] if s in pca_df.index]
    ax.scatter(pca_df.loc[samps, 'PC1'], pca_df.loc[samps, 'PC2'],
               s=160, color=ginfo['color'], marker='o',
               edgecolors='black', linewidths=1.4, zorder=5,
               label=f'{gname} (n={len(samps)})')

ax.set_xlabel(f'PC1 ({var_exp[0]:.1f}% variance)', fontsize=12)
ax.set_ylabel(f'PC2 ({var_exp[1]:.1f}% variance)', fontsize=12)
ax.set_title('OSD-564 PCA, RR-10 FLT (n=5) vs GC (n=5)\n'
             f'VST counts, top {TOP_N:,} variable genes, z-scored',
             fontsize=11.5, fontweight='bold', pad=12)

ax.legend(fontsize=10, frameon=True, framealpha=0.9,
          edgecolor='#cccccc', loc='best', markerscale=1.1,
          labelspacing=1.0, handlelength=1.5, borderpad=0.8)

ax.axhline(0, color='#cccccc', lw=0.8, ls='--', zorder=0)
ax.axvline(0, color='#cccccc', lw=0.8, ls='--', zorder=0)

plt.tight_layout()
pca_path = 'OSD564_PCA_plot.png'
plt.savefig(pca_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved -> {pca_path}")

# ──────────────────────────────────────────────────────────────────
# 4.  Volcano data  -- RR-10 FLT vs GC
# ──────────────────────────────────────────────────────────────────
print("\nPreparing volcano plot ...")

COMP_TAG   = '(Space Flight)v(Ground Control)'   # log2(FLT / GC): positive = up in flight
log2fc_col = f'Log2fc_{COMP_TAG}'
pval_col   = f'P.value_{COMP_TAG}'
adjp_col   = f'Adj.p.value_{COMP_TAG}'
symbol_col = 'SYMBOL'

for c in [log2fc_col, pval_col, adjp_col, symbol_col]:
    if c not in de.columns:
        raise KeyError(f"column not found in DE table: {c}")

df_v = de[[symbol_col, log2fc_col, pval_col, adjp_col]].copy()
df_v.columns = ['symbol', 'log2fc', 'pval', 'adjpval']
df_v['symbol']  = df_v['symbol'].fillna('').astype(str)
for c in ['log2fc', 'pval', 'adjpval']:
    df_v[c] = pd.to_numeric(df_v[c], errors='coerce')
# keep every gene with statistics -- missing SYMBOL does NOT drop a gene
df_v = df_v.dropna(subset=['log2fc', 'pval', 'adjpval'])
df_v['-log10p'] = -np.log10(df_v['pval'].clip(lower=1e-300))
df_v['first_symbol'] = df_v['symbol'].str.split('|').str[0]
N_TESTED = len(df_v)
print(f"Genes with statistics (tested): {N_TESTED:,}  "
      f"(of which {int((df_v['symbol'] == '').sum()):,} without SYMBOL)")

# Significance thresholds
FC_THR  = 1.0   # |log2FC| >= 1
ADJ_THR = 0.05  # adjusted p-value
NOM_THR = 0.05  # nominal p-value

adj_sig_up   = (df_v['log2fc'] >=  FC_THR) & (df_v['adjpval'] <= ADJ_THR)
adj_sig_down = (df_v['log2fc'] <= -FC_THR) & (df_v['adjpval'] <= ADJ_THR)
nom_sig_up   = (df_v['log2fc'] >=  FC_THR) & (df_v['pval'] <= NOM_THR) & ~adj_sig_up
nom_sig_down = (df_v['log2fc'] <= -FC_THR) & (df_v['pval'] <= NOM_THR) & ~adj_sig_down
not_sig      = ~(adj_sig_up | adj_sig_down | nom_sig_up | nom_sig_down)

n_adj_up, n_adj_down = int(adj_sig_up.sum()), int(adj_sig_down.sum())
n_nom_up, n_nom_down = int(nom_sig_up.sum()), int(nom_sig_down.sum())
print(f"Adj.p<=0.05 & |FC|>=1  -- UP: {n_adj_up}  DOWN: {n_adj_down}")
print(f"Nom.p<=0.05 & |FC|>=1 (not adj.) -- UP: {n_nom_up}  DOWN: {n_nom_down}")

XLAB   = 'log2 fold change (spaceflight / ground control)'
YLAB   = '-log10(nominal p-value)'
TITLE  = 'OSD-564 Volcano Plot, RR-10 FLT (n=5) vs RR-10 GC (n=5), 18-19 weeks'
LEGEND_TITLE = f'n counted over all tested genes (N = {N_TESTED:,})'

LBL_NOM_UP   = f'Nom. sig. up in spaceflight (p<=0.05, |log2FC|>=1, not adj. sig.), n = {n_nom_up}'
LBL_NOM_DOWN = f'Nom. sig. down in spaceflight (p<=0.05, |log2FC|>=1, not adj. sig.), n = {n_nom_down}'
LBL_ADJ_UP   = f'Adj. sig. up in spaceflight (adj p<=0.05, |log2FC|>=1), n = {n_adj_up}'
LBL_ADJ_DOWN = f'Adj. sig. down in spaceflight (adj p<=0.05, |log2FC|>=1), n = {n_adj_down}'

y_col = '-log10p'


def draw_volcano_layers(ax, star_size=90):
    """Background, nominal and adjusted layers shared by both volcano figures."""
    ax.scatter(df_v.loc[not_sig, 'log2fc'], df_v.loc[not_sig, y_col],
               s=7, color='#BDBDBD', alpha=0.40, linewidths=0, rasterized=True,
               label='Not significant')
    ax.scatter(df_v.loc[nom_sig_up, 'log2fc'], df_v.loc[nom_sig_up, y_col],
               s=18, color='#EF9A9A', alpha=0.80, linewidths=0.3,
               edgecolors='#C62828', rasterized=True, label=LBL_NOM_UP)
    ax.scatter(df_v.loc[nom_sig_down, 'log2fc'], df_v.loc[nom_sig_down, y_col],
               s=18, color='#90CAF9', alpha=0.80, linewidths=0.3,
               edgecolors='#1565C0', rasterized=True, label=LBL_NOM_DOWN)
    ax.scatter(df_v.loc[adj_sig_up, 'log2fc'], df_v.loc[adj_sig_up, y_col],
               s=star_size, color='#C62828', alpha=1.0, linewidths=0.5,
               edgecolors='#8B0000', marker='*', zorder=5, label=LBL_ADJ_UP)
    ax.scatter(df_v.loc[adj_sig_down, 'log2fc'], df_v.loc[adj_sig_down, y_col],
               s=star_size, color='#1565C0', alpha=1.0, linewidths=0.5,
               edgecolors='#0D47A1', marker='*', zorder=5, label=LBL_ADJ_DOWN)
    ax.axhline(-np.log10(NOM_THR), color='#9E9E9E', lw=1.0, ls='--', alpha=0.7,
               label=f'Nominal p = {NOM_THR}')
    ax.axvline( FC_THR, color='#9E9E9E', lw=1.0, ls=':', alpha=0.7)
    ax.axvline(-FC_THR, color='#9E9E9E', lw=1.0, ls=':', alpha=0.7)
    ax.set_xlabel(XLAB, fontsize=12)
    ax.set_ylabel(YLAB, fontsize=12)


# ── 4a. plain volcano ─────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 7.5))
draw_volcano_layers(ax)
ax.set_title(TITLE, fontsize=11.5, fontweight='bold', pad=14)

y_axis_top = df_v[y_col].quantile(0.9998) * 1.15
ax.set_xlim(df_v['log2fc'].quantile(0.001) - 0.4,
            df_v['log2fc'].quantile(0.999) + 0.4)
ax.set_ylim(-0.2, y_axis_top)

# legend below the axes so it never hides data points
ax.legend(fontsize=8, frameon=True, framealpha=0.92, edgecolor='#cccccc',
          loc='upper center', ncol=2, bbox_to_anchor=(0.5, -0.09),
          title=LEGEND_TITLE, title_fontsize=8)

volcano_path = 'OSD564_Volcano_plot.png'
plt.savefig(volcano_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved -> {volcano_path}")

# ──────────────────────────────────────────────────────────────────
# 4b. Volcano plot with curated gene annotations
# ──────────────────────────────────────────────────────────────────
print("\nPreparing annotated volcano plot ...")

# Curated genes: canonical symbols (first token of a pipe-joined SYMBOL)
CURATED_GENES = [
    "Adam8", "Apoc1", "Aqp1", "Arpc2", "Arap2", "Baiap2l2", "Bcat1", "C4a", "Casr",
    "Ect2", "Epcam", "Folr1", "Gas2", "Gas2l3", "Ghsr", "Itga2",
    "Lpl", "Myh7b", "Rac2", "Rhof", "S100a4", "S100a8", "Smtnl2",
    "Timp1", "Tmsb4x", "Tnni2", "Wdr1", "Wnt7a", "Zfas1",
    "Mrps25", "Chchd1", "Rmrp",
    "Ndufb4",    # SYMBOL "Ndufb4|Ndufb4c|Ndufb4b", ENSMUSG00000022820 (not the Ndufb4c paralog)
    "Ndufaf6",
    "Atp5mk",    # SYMBOL "Atp5mk|Gm35857", ENSMUSG00000071528
    "Pet100",
    "COX2",      # ENSMUSG00000064354, mitochondrially encoded mt-Co2
]
# Extra genes of interest shown with a distinct marker (adjusted-significant, |log2FC| < 1)
EXTRA_GENES = ["Vps4a", "Ist1", "Prkaa1"]

DISPLAY_NAME = {'COX2': 'mt-Co2'}   # SYMBOL -> label shown on the figure


def match_genes(gene_list, df):
    """
    For each curated gene, find DE rows whose SYMBOL (split on '|') contains it:
    first-token matches are preferred; any-token matches are used only when no
    first-token match exists.  If several rows match, keep the smallest p-value.
    Returns a DataFrame indexed by ENSEMBL with a 'display' column.
    """
    toks = df['symbol'].str.split('|')
    first = toks.str[0].str.lower()
    out = []
    for g in gene_list:
        gl = g.lower()
        m = first == gl
        if not m.any():
            m = toks.apply(lambda t: gl in [x.lower() for x in t])
        if not m.any():
            print(f"  [NOT FOUND] {g}")
            continue
        hit = df[m].sort_values('pval').iloc[[0]].copy()
        hit['curated'] = g
        out.append(hit)
    res = pd.concat(out)
    res['display'] = res['first_symbol'].map(lambda s: DISPLAY_NAME.get(s, s))
    return res


df_annot = match_genes(CURATED_GENES, df_v)
df_extra = match_genes(EXTRA_GENES, df_v)
print(f"Found {len(df_annot)} of {len(CURATED_GENES)} curated genes; "
      f"{len(df_extra)} of {len(EXTRA_GENES)} extra genes")
for ens, r in pd.concat([df_annot, df_extra]).iterrows():
    print(f"   {r['display']:<10s} {ens}  SYMBOL='{r['symbol']}'  "
          f"log2FC={r['log2fc']:+.2f}  p={r['pval']:.2e}  adj p={r['adjpval']:.4f}")

fig, ax = plt.subplots(figsize=(14, 11))
draw_volcano_layers(ax)

ax.scatter(df_annot['log2fc'], df_annot[y_col],
           s=60, color='#FF6F00', alpha=1.0, linewidths=0.8,
           edgecolors='#E65100', zorder=6,
           label=f'Curated genes (n={len(df_annot)})')
ax.scatter(df_extra['log2fc'], df_extra[y_col],
           s=110, color='#2E7D32', alpha=1.0, linewidths=0.9, marker='D',
           edgecolors='#1B5E20', zorder=7,
           label=f'Vps4a, Ist1, Prkaa1 (adj p<=0.05, |log2FC|<1) (n={len(df_extra)})')

ax.set_title(TITLE + '\ncurated gene annotations',
             fontsize=11.5, fontweight='bold', pad=14)

y_axis_top_a = df_v[y_col].quantile(0.9998) * 1.15
ax.set_xlim(df_v['log2fc'].quantile(0.001) - 2.5,
            df_v['log2fc'].quantile(0.999) + 1.5)
ax.set_ylim(-0.2, y_axis_top_a * 1.12)

df_lab = pd.concat([df_annot, df_extra])
texts, annot_coords = [], []
for _, row in df_lab.iterrows():
    txt = ax.text(row['log2fc'], row[y_col], row['display'],
                  fontsize=8, ha='center', va='center',
                  color='#111111', fontweight='bold', zorder=8, clip_on=False,
                  bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.85))
    texts.append(txt)
    annot_coords.append((row['log2fc'], row[y_col]))

try:
    from adjustText import adjust_text
    adjust_text(texts, ax=ax,
                x=df_lab['log2fc'].values, y=df_lab[y_col].values,
                force_text=(0.5, 1.0), force_static=(0.15, 0.3),
                force_pull=(0.01, 0.02), force_explode=(3.0, 5.0),
                expand=(1.4, 2.0), max_move=(300, 300), iter_lim=1000,
                min_arrow_len=99999)
except ImportError:
    for txt in texts:
        x0, y0 = txt.get_position()
        txt.set_position((x0, y0 + 0.25))

# Connector lines from each dot to its (repositioned) label
for txt, (px, py) in zip(texts, annot_coords):
    lx, ly = txt.get_position()
    if np.hypot(lx - px, ly - py) > 0.08:
        ax.annotate('', xy=(px, py), xytext=(lx, ly), xycoords='data', textcoords='data',
                    zorder=6, arrowprops=dict(arrowstyle='-', color='#555555',
                                              lw=0.8, shrinkA=4, shrinkB=4))

# legend below the axes so it never hides data points or labels
ax.legend(fontsize=9.5, frameon=True, framealpha=0.92, edgecolor='#cccccc',
          loc='upper center', ncol=2, bbox_to_anchor=(0.5, -0.06),
          markerscale=1.2, labelspacing=0.6, handlelength=1.8, borderpad=0.8,
          title=LEGEND_TITLE, title_fontsize=9.5)

annot_volcano_path = 'OSD564_Volcano_plot_annotated.png'
plt.savefig(annot_volcano_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved -> {annot_volcano_path}")

# ──────────────────────────────────────────────────────────────────
# 5.  Spreadsheet of significant genes (same gene set / thresholds as the legend)
# ──────────────────────────────────────────────────────────────────
print("\nExporting significant genes spreadsheet ...")

df_export = de[['SYMBOL', 'GENENAME', 'ENTREZID']].copy()
df_export['log2FC_SpaceflightVsGC'] = pd.to_numeric(de[log2fc_col], errors='coerce')
df_export['P_value']                 = pd.to_numeric(de[pval_col],   errors='coerce')
df_export['Adj_P_value']             = pd.to_numeric(de[adjp_col],   errors='coerce')
df_export = df_export.dropna(subset=['log2FC_SpaceflightVsGC', 'P_value', 'Adj_P_value'])
assert len(df_export) == N_TESTED
df_export.insert(1, 'Display_symbol',
                 df_export['SYMBOL'].fillna('').astype(str).str.split('|').str[0]
                 .map(lambda s: DISPLAY_NAME.get(s, s)))


def classify(row):
    fc, p, adjp = row['log2FC_SpaceflightVsGC'], row['P_value'], row['Adj_P_value']
    if   adjp <= ADJ_THR and fc >=  FC_THR: return 'Adj.sig UP'
    elif adjp <= ADJ_THR and fc <= -FC_THR: return 'Adj.sig DOWN'
    elif p    <= NOM_THR and fc >=  FC_THR: return 'Nom.sig UP'
    elif p    <= NOM_THR and fc <= -FC_THR: return 'Nom.sig DOWN'
    else: return 'Not significant'


df_export['Significance'] = df_export.apply(classify, axis=1)
df_export['Direction']    = df_export['log2FC_SpaceflightVsGC'].apply(
    lambda x: 'UP in Spaceflight' if x > 0 else 'DOWN in Spaceflight')

adj_sig = df_export[df_export['Significance'].str.startswith('Adj.sig')].sort_values('Adj_P_value')
nom_sig = df_export[df_export['Significance'].str.startswith('Nom.sig')].sort_values('P_value')

xl_adj_up   = int((adj_sig['Significance'] == 'Adj.sig UP').sum())
xl_adj_down = int((adj_sig['Significance'] == 'Adj.sig DOWN').sum())
xl_nom_up   = int((nom_sig['Significance'] == 'Nom.sig UP').sum())
xl_nom_down = int((nom_sig['Significance'] == 'Nom.sig DOWN').sum())
print(f"  Adj.sig:  {len(adj_sig)} genes  ({xl_adj_up} UP, {xl_adj_down} DOWN)")
print(f"  Nom.sig (not adj.):  {len(nom_sig)} genes  ({xl_nom_up} UP, {xl_nom_down} DOWN)")
assert (xl_adj_up, xl_adj_down, xl_nom_up, xl_nom_down) == (n_adj_up, n_adj_down, n_nom_up, n_nom_down), \
    "xlsx counts differ from volcano legend counts"

col_order = ['ENSEMBL', 'SYMBOL', 'Display_symbol', 'GENENAME', 'ENTREZID',
             'log2FC_SpaceflightVsGC', 'P_value', 'Adj_P_value', 'Significance', 'Direction']
all_sig = pd.concat([adj_sig, nom_sig]).reset_index()[col_order]

counts = pd.DataFrame([
    ['Genes with statistics (tested)', N_TESTED, ''],
    ['Adj.sig UP   (adj p<=0.05, log2FC>=1)',  n_adj_up,   'Adj.p Significant (0.05)'],
    ['Adj.sig DOWN (adj p<=0.05, log2FC<=-1)', n_adj_down, 'Adj.p Significant (0.05)'],
    ['Nom.sig UP   (p<=0.05, log2FC>=1, not adj.sig)',  n_nom_up,   'Nominal p Significant (0.05)'],
    ['Nom.sig DOWN (p<=0.05, log2FC<=-1, not adj.sig)', n_nom_down, 'Nominal p Significant (0.05)'],
], columns=['Category', 'n', 'Sheet'])
counts['Note'] = 'counts match the volcano legend; computed over all tested genes incl. missing SYMBOL'

xl_path = 'OSD564_Significant_Genes.xlsx'
with pd.ExcelWriter(xl_path, engine='openpyxl') as writer:
    counts.to_excel(writer, sheet_name='Counts', index=False)
    all_sig.to_excel(writer, sheet_name='All Significant', index=False)
    adj_sig.reset_index()[col_order].to_excel(writer, sheet_name='Adj.p Significant (0.05)', index=False)
    nom_sig.reset_index()[col_order].to_excel(writer, sheet_name='Nominal p Significant (0.05)', index=False)

print(f"Saved -> {xl_path}")
print("\nDone.")
