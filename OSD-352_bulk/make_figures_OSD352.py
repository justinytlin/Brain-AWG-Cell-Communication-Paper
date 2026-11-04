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
import srep_style as srep
srep.apply()                      # Scientific Reports style: Arial 8 pt, >= 1 pt lines, white, no legend frames
import srep_panels                # shared PCA / volcano / QC panel drawers (identical in all datasets)
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
warnings.filterwarnings('ignore')

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
    'Spaceflight': {
        'samples': ['RR3_BRN_FLT_F1', 'RR3_BRN_FLT_F2', 'RR3_BRN_FLT_F7'],
        'color':   srep.COLORS['flight'],
    },
    'Ground control': {
        'samples': ['RR3_BRN_GC_G7', 'RR3_BRN_GC_G8', 'RR3_BRN_GC_G9'],
        'color':   srep.COLORS['ground'],
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

# shared 60 x 60 mm PCA panel: short sample IDs (F1, F2, F7, G7, G8, G9) next to every point;
# the "top 2,000 variable genes" setting goes in the figure legend
pca_samples = [s for g in GROUP_INFO.values() for s in g['samples'] if s in pca_df.index]
pca_groups  = [g for g, info in GROUP_INFO.items() for s in info['samples'] if s in pca_df.index]
srep_panels.pca_panel(pca_df.loc[pca_samples, ['PC1', 'PC2']],
                      labels=[s.split('_')[-1] for s in pca_samples],
                      groups=pca_groups, pc_var=(var_exp[0], var_exp[1]),
                      title='OSD-352, spaceflight vs ground control',
                      out_base='OSD352_PCA_plot')

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
# one short in-panel line; the full description belongs to the figure legend.  The full
# "OSD-352 whole-brain bulk, spaceflight vs ground control (n = 3 vs 3)" is 86 mm wide at
# 8 pt and does not fit the 90 mm panel; this 80 mm variant keeps the parallel wording.
TITLE_MAIN = 'OSD-352 whole brain, spaceflight vs ground control (n = 3 vs 3)'

# legend counts (over all tested genes) handed to the shared panel drawer
COUNTS = {'n_nom_up': n_nom_up, 'n_nom_down': n_nom_down,
          'n_adj_up': n_adj_up, 'n_adj_down': n_adj_down}
assert COUNTS == srep_panels.volcano_counts(df_all)

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

def draw_annotated(df_points, df_labels, base):
    """Shared 90 x 85 mm volcano: annotated genes as orange circles (1 pt black edge),
    the df_labels subset labelled in italics by adjustText with 1 pt leader lines."""
    srep_panels.volcano_panel(df_all, 'volcano_supp', TITLE_MAIN,
                              label_genes=df_labels.rename(columns={'label': 'display'}),
                              highlight_genes=df_points, out_base=base, counts=COUNTS)

# ──────────────────────────────────────────────────────────────────
# 6a. Volcano plot (no gene labels)
# ──────────────────────────────────────────────────────────────────
print("\nPreparing volcano plot …")
srep_panels.volcano_panel(df_all, 'volcano_supp', TITLE_MAIN, label_genes=None,
                          highlight_genes=None, out_base='OSD352_Volcano_plot', counts=COUNTS)

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

# At 90 x 85 mm not every annotated gene can carry a label.  All annotated genes are
# still drawn (orange); labels are drawn for the genes discussed in the text plus the
# well-separated remainder, so that no label overlaps another.
MUST_LABEL = ['Tcf7l1', 'Rac2', 'Itga2', 'Cpt1b', 'Hpdl', 'Bbox1', 'Pask', 'Sdc4']
EXTRA_LABEL = ['Psd4', 'Osgin1', 'Myo5c', 'Adora2a', 'Agt', 'Gcat']
df_lab = df_annot[df_annot['label'].isin(MUST_LABEL + EXTRA_LABEL)]
missing_must = sorted(set(MUST_LABEL) - set(df_lab['label']))
if missing_must:
    print(f"  [WARN] required labels not found in the DE table: {missing_must}")

draw_annotated(df_annot, df_lab, 'OSD352_Volcano_plot_annotated')

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
    "Atp5mk|Gm35857", "Pet100", "mt-Co2|COX2",   # displayed as mt-Co2 (as in OSD-562 / OSD-564)
]

df_annot_v2, missing_v2 = select_annotation(annotation_genes_v2)
df_annot_v2 = add_adj_sig_labels(df_annot_v2)
print(f"Found {len(df_annot_v2)} annotation genes (incl. adjusted-significant) in DE data; "
      f"missing: {missing_v2 if missing_v2 else 'none'}")
for r in df_annot_v2.sort_values('pval').itertuples():
    print(f"  {r.label:<8s} {r.ENSEMBL}  log2FC={r.log2fc:+.3f}  p={r.pval:.4g}  adj p={r.adjpval:.3g}")

# Labels: the adjusted-significant gene, every mitochondrial gene of the list, and the
# cytoskeletal / annotated genes that pass p <= 0.05 and |log2FC| >= 1 (the rest stay
# as unlabelled orange points to avoid overlapping text at 90 x 85 mm).
MITO_LABELS = ['Mrps25', 'Chchd1', 'Rmrp', 'Ndufb4', 'Ndufaf6', 'Atp5mk', 'Pet100', 'mt-Co2']
sig_v2 = (df_annot_v2['pval'] <= NOM_THR) & (df_annot_v2['log2fc'].abs() >= FC_THR)
df_lab_v2 = df_annot_v2[sig_v2 | df_annot_v2['label'].isin(MITO_LABELS + ['Tcf7l1'])]
print(f"  Labelled: {', '.join(df_lab_v2['label'])}")

draw_annotated(df_annot_v2, df_lab_v2, 'OSD352_Volcano_plot_annotated_mito_cyto')

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
