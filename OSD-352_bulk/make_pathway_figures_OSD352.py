"""
Harmonized GSEA NES summary panels for OSD-352 (RR-3 whole-brain bulk RNA-seq,
FLT n = 3 vs GC n = 3).

Panel sets are IDENTICAL to the OSD-564 / OSD-562 figures:
  * GO BP mitochondrial (7 terms, matched by GO ID)
  * KEGG sEV pathways   (7 terms, matched by exact KEGG name)
  * KEGG TNT pathways   (3 terms, matched by exact KEGG name)

Inputs (produced by KEGG_enrichment_OSD352.py / GO_BP_enrichment_OSD352.py, both
run on the (Space Flight)v(Ground Control) contrast so NES > 0 = enriched in
spaceflight):
  KEGG_results_OSD352/GSEA_KEGG_results.csv
  GO_BP_results_OSD352/GSEA_GO_BP_results.csv
  KEGG_results_OSD352/GSEA_ranked_gene_list.csv   (+ GSEA_prerank/gene_sets.gmt)
  GO_BP_results_OSD352/GSEA_ranked_gene_list.csv  (+ GSEA_prerank/gene_sets.gmt)

Every term in a panel is drawn with its FDR label whether or not it is
significant; an asterisk marks FDR <= 0.05.  The NES panels are drawn by the
shared srep_nes_panel.nes_panel (identical layout and term order in OSD-564,
OSD-562 and OSD-352; xlim (-3.2, 3.2) for the mitochondrial panel, (-2.2, 2.2)
for the sEV and TNT panels).
"""
import os
import re
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import srep_style as srep
srep.apply()                      # Scientific Reports style: Arial 8 pt, >= 1 pt lines, white, no legend frames
import matplotlib.pyplot as plt
import gseapy as gp
from gseapy import gseaplot

from srep_nes_panel import nes_panel, XLIM

warnings.filterwarnings('ignore')

# `--nes-only`: redraw the NES summary panels from the saved GSEA tables without
# re-running prerank for the running-ES curves or rewriting the Excel summary.
NES_ONLY = '--nes-only' in sys.argv[1:]

# ──────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────
DATASET_TAG = 'OSD-352 whole-brain bulk RNA-seq: spaceflight vs ground control (n = 3 vs 3)'
PERM_NUM    = 1000
SEED        = 42

KEGG_DIR = 'KEGG_results_OSD352'
GO_DIR   = 'GO_BP_results_OSD352'

LIBS = {
    'KEGG': {
        'results_csv': f'{KEGG_DIR}/GSEA_KEGG_results.csv',
        'ranked_csv':  f'{KEGG_DIR}/GSEA_ranked_gene_list.csv',
        'gmt':         f'{KEGG_DIR}/GSEA_prerank/gene_sets.gmt',
        'plot_dir':    f'{KEGG_DIR}/pathway_figures',
    },
    'GO_BP': {
        'results_csv': f'{GO_DIR}/GSEA_GO_BP_results.csv',
        'ranked_csv':  f'{GO_DIR}/GSEA_ranked_gene_list.csv',
        'gmt':         f'{GO_DIR}/GSEA_prerank/gene_sets.gmt',
        'plot_dir':    f'{GO_DIR}/pathway_figures',
    },
}

# Harmonized panels. For GO BP the key is the GO ID (matched against the
# "(GO:xxxxxxx)" suffix of the Enrichr term); for KEGG the key is the exact
# KEGG_2019_Mouse term name (case-insensitive exact match).
PANELS = {
    'Mito': {
        'lib':      'GO_BP',
        'panel':    'GO BP mitochondrial pathways',
        'title':    'GO BP mitochondrial terms, OSD-352 whole-brain bulk (n = 3 vs 3)',
        'xlim':     XLIM['mito'],
        'out_png':  f'{GO_DIR}/GSEA_Mito_NES_summary.png',
        'terms': [
            ('GO:0006119', 'Oxidative Phosphorylation'),
            ('GO:0033108', 'Mitochondrial Respiratory Chain Complex Assembly'),
            ('GO:0042775', 'Mitochondrial ATP Synthesis Coupled Electron Transport'),
            ('GO:0042776', 'Proton Motive Force-Driven Mitochondrial ATP Synthesis'),
            ('GO:0006120', 'Mitochondrial Electron Transport, NADH to Ubiquinone'),
            ('GO:0032981', 'Mitochondrial Respiratory Chain Complex I Assembly'),
            ('GO:0033617', 'Mitochondrial Respiratory Chain Complex IV Assembly'),
        ],
    },
    'sEVs_Docking': {
        'lib':      'KEGG',
        'panel':    'KEGG sEV pathways',
        'title':    'KEGG sEV pathways, OSD-352 whole-brain bulk (n = 3 vs 3)',
        'xlim':     XLIM['sEV'],
        'out_png':  f'{KEGG_DIR}/pathway_figures/GSEA_sEVs_Docking_NES_summary.png',
        'terms': [
            ('Sphingolipid signaling pathway',        'Sphingolipid signaling pathway'),
            ('Sphingolipid metabolism',               'Sphingolipid metabolism'),
            ('Phospholipase D signaling pathway',     'Phospholipase D signaling pathway'),
            ('Endocytosis',                           'Endocytosis'),
            ('Synaptic vesicle cycle',                'Synaptic vesicle cycle'),
            ('Phosphatidylinositol signaling system', 'Phosphatidylinositol signaling system'),
            ('Lysosome',                              'Lysosome'),
        ],
    },
    'TNTs_Factin': {
        'lib':      'KEGG',
        'panel':    'KEGG TNT pathways',
        'title':    'KEGG TNT pathways, OSD-352 whole-brain bulk (n = 3 vs 3)',
        'xlim':     XLIM['TNT'],
        'out_png':  f'{KEGG_DIR}/pathway_figures/GSEA_TNTs_Factin_NES_summary.png',
        'terms': [
            ('Axon guidance',                    'Axon guidance'),
            ('Regulation of actin cytoskeleton', 'Regulation of actin cytoskeleton'),
            ('Rap1 signaling pathway',           'Rap1 signaling pathway'),
        ],
    },
}

for lib in LIBS.values():
    os.makedirs(lib['plot_dir'], exist_ok=True)

# ──────────────────────────────────────────────────────────────────
# 1. Load GSEA result tables
# ──────────────────────────────────────────────────────────────────
def go_id_of(term):
    m = re.search(r'\((GO:\d{7})\)\s*$', str(term))
    return m.group(1) if m else None

results = {}
for lib_name, lib in LIBS.items():
    if not os.path.exists(lib['results_csv']):
        raise FileNotFoundError(
            f"{lib['results_csv']} not found: run the enrichment script for {lib_name} first.")
    df = pd.read_csv(lib['results_csv'], index_col=0)
    df['Term'] = df['Term'].astype(str).str.strip()
    if lib_name == 'GO_BP':
        df['match_key'] = df['Term'].map(go_id_of)
    else:
        df['match_key'] = df['Term'].str.lower()
    results[lib_name] = df
    print(f"{lib_name}: {len(df):,} gene sets loaded from {lib['results_csv']}")

def lookup(lib_name, key):
    """Return the result row for a harmonized key (GO ID or exact KEGG name)."""
    df = results[lib_name]
    k = key if lib_name == 'GO_BP' else key.lower()
    hit = df[df['match_key'] == k]
    if len(hit) == 0:
        return None
    if len(hit) > 1:
        print(f"  [WARN] {key}: {len(hit)} matches, using the first ({hit.iloc[0]['Term']})")
    return hit.iloc[0]

# ──────────────────────────────────────────────────────────────────
# 2. Collect the panel records
# ──────────────────────────────────────────────────────────────────
def fdr_label(fdr):
    if pd.isna(fdr):
        return 'not tested'
    return srep.fdr_label(float(fdr))   # 'FDR < 0.001' (incl. FDR == 0) or 'FDR = 0.88'

panel_tables = {}
print("\nHarmonized panel values:")
for pkey, pinfo in PANELS.items():
    recs = []
    print(f"  [{pinfo['panel']}]")
    for key, label in pinfo['terms']:
        row = lookup(pinfo['lib'], key)
        if row is None:
            print(f"    [MISSING] {label} ({key}) — shown as 'not tested'")
            recs.append({'key': key, 'label': label, 'Term': None,
                         'NES': np.nan, 'NOM p-val': np.nan, 'FDR': np.nan,
                         'FWER p-val': np.nan, 'Size': np.nan,
                         'Leading edge size': 0, 'Lead_genes': ''})
            continue
        n_lead = len(str(row['Lead_genes']).split(';')) if pd.notna(row['Lead_genes']) else 0
        size = np.nan
        tag_str = str(row.get('Tag %', ''))
        if '/' in tag_str:
            size = int(tag_str.split('/')[1])
        recs.append({'key': key, 'label': label, 'Term': row['Term'],
                     'NES': float(row['NES']), 'NOM p-val': float(row['NOM p-val']),
                     'FDR': float(row['FDR q-val']), 'FWER p-val': float(row['FWER p-val']),
                     'Size': size, 'Leading edge size': n_lead,
                     'Lead_genes': str(row['Lead_genes']) if pd.notna(row['Lead_genes']) else ''})
        print(f"    {label:<58s} NES={float(row['NES']):+.3f}  p={float(row['NOM p-val']):.4f}  "
              f"FDR={float(row['FDR q-val']):.4f}")
    panel_tables[pkey] = pd.DataFrame(recs)

# ──────────────────────────────────────────────────────────────────
# 3. NES summary bar plot per panel (shared srep_nes_panel.nes_panel)
#    All terms are drawn in the listed (harmonized) order, first term at the
#    top, so panels are directly comparable across OSD-564 / OSD-562 / OSD-352.
#    A term missing from the result table is drawn without a bar ("not tested").
# ──────────────────────────────────────────────────────────────────
print("\nCreating NES summary panels ...")
for pkey, pinfo in PANELS.items():
    tbl = panel_tables[pkey]
    terms = [{'label': row['label'], 'nes': row['NES'], 'fdr': row['FDR']}
             for _, row in tbl.iterrows()]
    nes_panel(terms, pinfo['title'], os.path.splitext(pinfo['out_png'])[0], xlim=pinfo['xlim'])

if NES_ONLY:
    print("\n--nes-only: skipping the running-ES plots (no prerank re-run) and the Excel summary.")
    sys.exit(0)

# ──────────────────────────────────────────────────────────────────
# 4. Individual enrichment (running-ES) plots for every harmonized term
#    The curves need the prerank object; re-run prerank restricted to the
#    harmonized gene sets (ES/NES per set are unchanged, FDR is not, so the
#    NES / p / FDR shown on the plot are taken from the full-library CSV).
# ──────────────────────────────────────────────────────────────────
def read_gmt(path):
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                sets[parts[0].strip()] = [g for g in parts[2:] if g]
    return sets

def safe_name(term):
    return (term.replace('/', '-').replace(' ', '_').replace(':', '')
                .replace(',', '').replace("'", '').replace('(', '').replace(')', ''))

print("\nGenerating individual enrichment plots ...")
for lib_name, lib in LIBS.items():
    wanted = {}
    for pkey, pinfo in PANELS.items():
        if pinfo['lib'] != lib_name:
            continue
        for _, row in panel_tables[pkey].iterrows():
            if row['Term'] is not None and pd.notna(row['NES']):
                wanted[row['Term']] = row
    if not wanted:
        continue
    if not (os.path.exists(lib['gmt']) and os.path.exists(lib['ranked_csv'])):
        print(f"  [SKIP] {lib_name}: GMT or ranked list missing")
        continue

    gmt_all = read_gmt(lib['gmt'])
    sub = {t: gmt_all[t] for t in wanted if t in gmt_all}
    for t in [t for t in wanted if t not in gmt_all]:
        print(f"  [WARN] {lib_name}: '{t}' not in {lib['gmt']}")
    if not sub:
        continue

    rnk = pd.read_csv(lib['ranked_csv'], index_col=0).squeeze('columns')
    rnk = rnk.sort_values(ascending=False)
    rnk = rnk[~rnk.index.duplicated(keep='first')]

    gres = gp.prerank(rnk=rnk, gene_sets=sub, outdir=None, seed=SEED,
                      permutation_num=PERM_NUM, min_size=10, max_size=500,
                      threads=4, verbose=False)

    for term, row in wanted.items():
        if term not in gres.results:
            print(f"  [MISSING] {lib_name}: {term}")
            continue
        r = dict(gres.results[term])
        # Use the full-library statistics for the annotation text
        r['nes']  = float(row['NES'])
        r['pval'] = float(row['NOM p-val'])
        r['fdr']  = float(row['FDR'])
        ofname = f"{lib['plot_dir']}/GSEA_{safe_name(term)}.png"
        try:
            gseaplot(rank_metric=gres.ranking, term=term, ofname=ofname, **r)
            plt.close('all')
            print(f"  [OK]  {term}  →  {ofname}")
        except Exception as e:
            print(f"  [ERROR] {term}: {e}")

# ──────────────────────────────────────────────────────────────────
# 5. Excel summary of the harmonized panels
# ──────────────────────────────────────────────────────────────────
xl_path = f'{KEGG_DIR}/OSD352_pathway_GSEA_summary.xlsx'
print(f"\nWriting Excel summary → {xl_path}")

def autosize(ws):
    for col_cells in ws.columns:
        width = max(len(str(c.value)) if c.value is not None else 0 for c in col_cells)
        ws.column_dimensions[col_cells[0].column_letter].width = min(width + 2, 60)

with pd.ExcelWriter(xl_path, engine='openpyxl') as writer:
    all_rows = []
    for pkey, pinfo in PANELS.items():
        t = panel_tables[pkey].copy()
        t.insert(0, 'Panel', pinfo['panel'])
        t.insert(1, 'Library', 'GO_Biological_Process_2026' if pinfo['lib'] == 'GO_BP' else 'KEGG_2019_Mouse')
        t = t.rename(columns={'key': 'Match key', 'label': 'Pathway', 'Term': 'Matched term',
                              'FDR': 'FDR q-val', 'Size': 'Gene set size',
                              'Lead_genes': 'Leading edge genes'})
        t['Leading edge genes'] = t['Leading edge genes'].str.replace(';', ', ')
        t['FDR label'] = t['FDR q-val'].map(fdr_label)
        cols = ['Panel', 'Library', 'Match key', 'Pathway', 'Matched term', 'NES', 'NOM p-val',
                'FDR q-val', 'FDR label', 'FWER p-val', 'Gene set size', 'Leading edge size',
                'Leading edge genes']
        t = t[[c for c in cols if c in t.columns]]
        sheet = pinfo['panel'][:31]
        t.to_excel(writer, sheet_name=sheet, index=False)
        autosize(writer.sheets[sheet])
        all_rows.append(t)
    pd.concat(all_rows, ignore_index=True).to_excel(writer, sheet_name='All panels', index=False)
    autosize(writer.sheets['All panels'])

print("\nDone.")
