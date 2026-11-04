"""
make_figures.py  (OSD-562, re-run DESeq2 4 vs 3)
------------------------------------------------
PCA (7-sample VST from rerun_deseq2_OSD562.py), volcano plots and the
significant-gene spreadsheet for
    ISS-T FLT 20 weeks (n = 4)  vs  ISS-T GC 20 weeks (n = 3)
using OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3.csv (not the GeneLab 3 v 6 table).

Style: Scientific Reports (srep_style.py) -- Arial 8 pt everywhere, sentence
case, no line < 1 pt; PCA 60 x 60 mm, volcano panels 90 x 85 mm, each saved
as 300 dpi PNG + vector PDF (srep.save).
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import srep_style as srep
srep.apply()
import srep_panels                      # shared PCA / volcano / QC panel drawers (identical in all datasets)
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

DE_FILE  = 'OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3.csv'
VST_FILE = 'OSD562_VST_counts_ISS-T_YNG.csv'
COMP_TAG = '(Space Flight)v(Ground Control)'
DATASET_LINE = 'OSD-562, spaceflight vs ground control (n = 4 vs 3)'
PCA_LINE     = 'OSD-562, spaceflight vs ground control'   # 60 mm panel: the n goes in the figure legend


# ──────────────────────────────────────────────────────────────────
# 1.  Load data
# ──────────────────────────────────────────────────────────────────
print("Loading VST counts (7-sample re-run) ...")
vst = pd.read_csv(VST_FILE, index_col=0)

print("Loading differential expression data (re-run) ...")
de = pd.read_csv(DE_FILE, index_col=0, low_memory=False)

print(f"VST matrix: {srep.fmt_int(vst.shape[0])} genes x {vst.shape[1]} samples")
print(f"DE table:   {srep.fmt_int(de.shape[0])} genes x {de.shape[1]} columns")

# ──────────────────────────────────────────────────────────────────
# 2.  Group definitions
# ──────────────────────────────────────────────────────────────────
GROUP_INFO = {
    'Spaceflight': {
        'samples': ['RRRM2_HPC_FLT_ISS-T_YNG_FY3',
                    'RRRM2_HPC_FLT_ISS-T_YNG_FY4',
                    'RRRM2_HPC_FLT_ISS-T_YNG_FY6',
                    'RRRM2_HPC_FLT_ISS-T_YNG_FY10'],
        'color':   srep.COLORS['flight'],
    },
    'Ground control': {
        'samples': ['RRRM2_HPC_GC_ISS-T_YNG_GY5',
                    'RRRM2_HPC_GC_ISS-T_YNG_GY8',
                    'RRRM2_HPC_GC_ISS-T_YNG_GY10'],
        'color':   srep.COLORS['ground'],
    },
}
SAMPLES_7 = GROUP_INFO['Spaceflight']['samples'] + GROUP_INFO['Ground control']['samples']
missing = [s for s in SAMPLES_7 if s not in vst.columns]
assert not missing, f"VST is missing samples: {missing}"
vst = vst[SAMPLES_7]            # subset to the 7 samples BEFORE variance ranking

# ──────────────────────────────────────────────────────────────────
# 3.  PCA  -- top 2000 most-variable genes (7 samples only)
# ──────────────────────────────────────────────────────────────────
print("\nRunning PCA ...")
TOP_N = 2000
gene_var = vst.var(axis=1)
top_genes = gene_var.nlargest(TOP_N).index
X = vst.loc[top_genes].T          # samples x genes

X_scaled = StandardScaler().fit_transform(X)
n_comp = min(10, X.shape[0] - 1)
pca = PCA(n_components=n_comp, random_state=42)
pca_coords = pca.fit_transform(X_scaled)
pca_df = pd.DataFrame(pca_coords, index=X.index,
                      columns=[f'PC{i+1}' for i in range(n_comp)])

var_exp = pca.explained_variance_ratio_ * 100
print(f"PC1: {var_exp[0]:.1f}%   PC2: {var_exp[1]:.1f}%   PC3: {var_exp[2]:.1f}%")
pca_df.assign(group=['ISS-T FLT'] * 4 + ['ISS-T GC'] * 3).to_csv('OSD562_PCA_coordinates.csv')

# shared 60 x 60 mm PCA panel: short sample IDs (FY3 ... GY10) next to every point
pca_samples = [s for g in GROUP_INFO.values() for s in g['samples'] if s in pca_df.index]
pca_groups  = [g for g, info in GROUP_INFO.items() for s in info['samples'] if s in pca_df.index]
srep_panels.pca_panel(pca_df.loc[pca_samples, ['PC1', 'PC2']],
                      labels=[s.split('_')[-1] for s in pca_samples],
                      groups=pca_groups, pc_var=(var_exp[0], var_exp[1]),
                      title=PCA_LINE, out_base='OSD562_PCA_plot')

# ──────────────────────────────────────────────────────────────────
# 4.  Volcano plot  -- ISS-T Young FLT vs GC
# ──────────────────────────────────────────────────────────────────
print("\nPreparing volcano plot ...")
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

n_adj_up, n_adj_down = int(adj_sig_up.sum()), int(adj_sig_down.sum())
n_nom_up, n_nom_down = int(nom_sig_up.sum()), int(nom_sig_down.sum())
print(f"Adj.p<=0.05 & |FC|>=1 -- UP: {n_adj_up}  DOWN: {n_adj_down}")
print(f"Nom.p<=0.05 & |FC|>=1 -- UP: {n_nom_up}  DOWN: {n_nom_down}")

# legend counts (over all tested genes) handed to the shared panel drawer
COUNTS = {'n_nom_up': n_nom_up, 'n_nom_down': n_nom_down,
          'n_adj_up': n_adj_up, 'n_adj_down': n_adj_down}
assert COUNTS == srep_panels.volcano_counts(df_v)
y_col = '-log10p'

# plain 90 x 85 mm volcano (no gene labels); adjusted entries with n = 0 become a text line
srep_panels.volcano_panel(df_v, 'volcano', DATASET_LINE, label_genes=None, highlight_genes=None,
                          out_base='OSD562_Volcano_plot', counts=COUNTS)

# ──────────────────────────────────────────────────────────────────
# 4b. Volcano plot with custom gene annotations
#     (all adjusted-significant genes are labelled as well)
# ──────────────────────────────────────────────────────────────────
print("\nPreparing annotated volcano plot ...")
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

# shared 90 x 85 mm volcano: the annotation genes (orange circles) and every
# adjusted-significant gene are labelled (adjustText, 1 pt leader lines)
label_df = pd.concat([df_annot, df_adj]).drop_duplicates()
label_df['display'] = [r['symbol'] if isinstance(r['symbol'], str) else idx
                       for idx, r in label_df.iterrows()]
srep_panels.volcano_panel(df_v, 'volcano', DATASET_LINE, label_genes=label_df,
                          highlight_genes=df_annot, out_base='OSD562_Volcano_plot_annotated',
                          counts=COUNTS)

# ──────────────────────────────────────────────────────────────────
# 5.  Spreadsheet of significant genes
# ──────────────────────────────────────────────────────────────────
print("\nExporting significant genes spreadsheet ...")
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
print(f"Saved -> {xl_path}")
print("\nDone.")
