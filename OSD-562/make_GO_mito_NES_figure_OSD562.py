"""
make_GO_mito_NES_figure_OSD562.py
---------------------------------
GO BP mitochondrial panel (7 terms, HARMONIZED with OSD-564 / OSD-352, matched
by GO ID, identical term order) from GO_BP_results_OSD562/GSEA_GO_BP_results.csv
(GO_Biological_Process_2026, re-run DESeq2 ISS-T FLT n = 4 vs ISS-T GC n = 3).
NES / FDR are taken from that table (no GSEA re-run).

The panel itself is drawn by the shared srep_nes_panel.nes_panel (Scientific
Reports style: Arial 8 pt, 90 x 60 mm, 300 dpi PNG + vector PDF).
"""
import os
import re
import warnings
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import srep_style as srep
srep.apply()

from srep_nes_panel import nes_panel, XLIM

warnings.filterwarnings('ignore')

TARGET_TERMS = [  # (GO ID, display label) -- fixed order, first term at the top
    ('GO:0006119', 'Oxidative Phosphorylation'),
    ('GO:0033108', 'Mitochondrial Respiratory Chain Complex Assembly'),
    ('GO:0042775', 'Mitochondrial ATP Synthesis Coupled Electron Transport'),
    ('GO:0042776', 'Proton Motive Force-Driven Mitochondrial ATP Synthesis'),
    ('GO:0006120', 'Mitochondrial Electron Transport, NADH To Ubiquinone'),
    ('GO:0032981', 'Mitochondrial Respiratory Chain Complex I Assembly'),
    ('GO:0033617', 'Mitochondrial Respiratory Chain Complex IV Assembly'),
]
TITLE    = 'GO BP mitochondrial terms, OSD-562 (n = 4 vs 3)'
GSEA_CSV = 'GO_BP_results_OSD562/GSEA_GO_BP_results.csv'
OUT_DIR  = 'GO_BP_results_OSD562'
os.makedirs(OUT_DIR, exist_ok=True)

print("Loading GSEA GO BP results ...")
gsea_df = pd.read_csv(GSEA_CSV, index_col=0)
gsea_df['Term'] = gsea_df['Term'].str.strip()

records = []
for go_id, label in TARGET_TERMS:
    hit = gsea_df[gsea_df['Term'].str.contains(re.escape(f'({go_id})'), regex=True)]
    if hit.empty:
        print(f"  [MISSING] {go_id} {label}")
        continue
    row = hit.iloc[0]
    records.append({'GO_ID': go_id, 'label': label, 'Term_in_library': row['Term'],
                    'NES': float(row['NES']), 'FDR': float(row['FDR q-val']),
                    'NOM_p': float(row['NOM p-val'])})
    print(f"  [OK] {go_id} {label}  ->  NES={float(row['NES']):.3f}  FDR={float(row['FDR q-val']):.4f}")

pd.DataFrame(records).to_csv(f'{OUT_DIR}/GSEA_Mito_NES_summary.csv', index=False)

terms = [{'label': r['label'], 'nes': r['NES'], 'fdr': r['FDR']} for r in records]
nes_panel(terms, TITLE, f'{OUT_DIR}/GSEA_Mito_NES_summary', xlim=XLIM['mito'])
print("Done.")
