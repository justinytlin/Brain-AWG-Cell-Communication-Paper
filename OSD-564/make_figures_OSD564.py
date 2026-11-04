"""
OSD-564 (GLDS-569, RR-10) volcano plots, PCA and significant-gene export,
Scientific Reports style (srep_style: Arial 8 pt, >= 1 pt lines, white background).

Contrast: (Space Flight) v (Ground Control)  ->  log2FC > 0 = up in spaceflight.

Legend counts are computed over ALL genes with non-missing statistics (18,517 tested
genes), including rows whose SYMBOL is missing.  Curated-gene matching splits
pipe-joined SYMBOL entries ("Ndufb4|Ndufb4c|Ndufb4b") and matches the first token
first (then any token); the display label is the first token.  "COX2"
(ENSMUSG00000064354) is the mitochondrially encoded mt-Co2 and is displayed as such.

Panels: OSD564_Volcano_plot (118 x 100 mm), OSD564_Volcano_plot_annotated (118 x 100 mm),
OSD564_PCA_plot (60 x 60 mm); each written as 300 dpi PNG + vector PDF.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import srep_style as srep; srep.apply()
import srep_panels                      # shared PCA / volcano / QC panel drawers (identical in all datasets)
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

STUDY_TAG = 'OSD-564, spaceflight vs ground control (n = 5 vs 5)'
PCA_TAG   = 'OSD-564, spaceflight vs ground control'   # 60 mm panel: the n goes in the figure legend

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
    'Spaceflight': {
        'samples': ['RR10_HPC_FLT_WT_F1', 'RR10_HPC_FLT_WT_F3', 'RR10_HPC_FLT_WT_F5',
                    'RR10_HPC_FLT_WT_F7', 'RR10_HPC_FLT_WT_F9'],
        'color':   srep.COLORS['flight'],
    },
    'Ground control': {
        'samples': ['RR10_HPC_GC_WT_G1', 'RR10_HPC_GC_WT_G3', 'RR10_HPC_GC_WT_G5',
                    'RR10_HPC_GC_WT_G7', 'RR10_HPC_GC_WT_G9'],
        'color':   srep.COLORS['ground'],
    },
}

# ──────────────────────────────────────────────────────────────────
# 3.  PCA  -- top 2000 most-variable genes, z-scored
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

# shared 60 x 60 mm PCA panel: short sample IDs (F1 ... G9) next to every point
pca_samples = [s for g in GROUP_INFO.values() for s in g['samples'] if s in pca_df.index]
pca_groups  = [g for g, info in GROUP_INFO.items() for s in info['samples'] if s in pca_df.index]
srep_panels.pca_panel(pca_df.loc[pca_samples, ['PC1', 'PC2']],
                      labels=[s.split('_')[-1] for s in pca_samples],
                      groups=pca_groups, pc_var=(var_exp[0], var_exp[1]),
                      title=PCA_TAG, out_base='OSD564_PCA_plot')

# ──────────────────────────────────────────────────────────────────
# 4.  Volcano data  -- spaceflight vs ground control
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
print(f"Genes with statistics (tested): {srep.fmt_int(N_TESTED)}  "
      f"(of which {srep.fmt_int((df_v['symbol'] == '').sum())} without SYMBOL)")

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

# legend counts (over all tested genes) handed to the shared panel drawer
COUNTS = {'n_nom_up': n_nom_up, 'n_nom_down': n_nom_down,
          'n_adj_up': n_adj_up, 'n_adj_down': n_adj_down}
assert COUNTS == srep_panels.volcano_counts(df_v)

# ── 4a. plain volcano (118 x 100 mm, no gene labels) ──────────────
srep_panels.volcano_panel(df_v, 'volcano', STUDY_TAG, label_genes=None, highlight_genes=None,
                          out_base='OSD564_Volcano_plot', counts=COUNTS)

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

# shared 118 x 100 mm volcano: all curated genes (orange circles) + Vps4a / Ist1 / Prkaa1
# (orange open diamonds) are labelled; adjustText places the labels with leader lines
LBL_EXTRA = 'Vps4a, Ist1, Prkaa1 (adjusted p ≤ 0.05, |log2FC| < 1)'
df_lab = pd.concat([df_annot, df_extra])
srep_panels.volcano_panel(df_v, 'volcano', STUDY_TAG, label_genes=df_lab, highlight_genes=df_annot,
                          out_base='OSD564_Volcano_plot_annotated', counts=COUNTS,
                          extra_genes=df_extra, extra_label=LBL_EXTRA)

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
