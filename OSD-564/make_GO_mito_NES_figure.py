"""
GO BP GSEA, mitochondrial terms -- harmonized 7-term panel (identical term order
in OSD-564, OSD-562 and OSD-352).  Terms are matched by GO ID in
GO_results/GSEA_BP_results.csv; NES / FDR are taken from that table (no GSEA re-run).
Output: GO_results/GSEA_mito_pathways_NES.png / .pdf (90 x 60 mm, Scientific
Reports style, drawn by the shared srep_nes_panel.nes_panel).
"""
import os
import warnings
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import srep_style as srep; srep.apply()

from srep_nes_panel import nes_panel, XLIM

warnings.filterwarnings('ignore')

# ──────────────────────────────────────────────────────────────────
# Configuration: harmonized GO BP mitochondrial panel (7 terms, fixed order;
# the first term is drawn at the top of the panel)
# ──────────────────────────────────────────────────────────────────
TARGET_GO_IDS = [
    'GO:0006119',   # Oxidative phosphorylation
    'GO:0033108',   # Mitochondrial respiratory chain complex assembly
    'GO:0042775',   # Mitochondrial ATP synthesis coupled electron transport
    'GO:0042776',   # Proton motive force-driven mitochondrial ATP synthesis
    'GO:0006120',   # Mitochondrial electron transport, NADH to ubiquinone
    'GO:0032981',   # Mitochondrial respiratory chain complex I assembly
    'GO:0033617',   # Mitochondrial respiratory chain complex IV assembly
]
TITLE    = 'GO BP mitochondrial terms, OSD-564 (n = 5 vs 5)'
GSEA_CSV = 'GO_results/GSEA_BP_results.csv'
OUT_DIR  = 'GO_results'
os.makedirs(OUT_DIR, exist_ok=True)

# ──────────────────────────────────────────────────────────────────
# 1. Load GSEA results and match by GO ID
# ──────────────────────────────────────────────────────────────────
print("Loading GSEA BP results ...")
gsea_df = pd.read_csv(GSEA_CSV, index_col=0)
gsea_df['Term'] = gsea_df['Term'].str.strip()

terms = []
for go_id in TARGET_GO_IDS:
    hits = gsea_df[gsea_df['Term'].str.contains(f'({go_id})', regex=False)]
    if hits.empty:
        print(f"  [MISSING] {go_id}")
        continue
    row = hits.iloc[0]
    terms.append({'label': row['Term'], 'nes': float(row['NES']), 'fdr': float(row['FDR q-val'])})
    print(f"  {row['Term']:<70s} NES={row['NES']:+.3f}  FDR={row['FDR q-val']:.4f}")

# ──────────────────────────────────────────────────────────────────
# 2. Plot (shared panel; xlim (-3.2, 3.2) for every mitochondrial panel)
# ──────────────────────────────────────────────────────────────────
nes_panel(terms, TITLE, f'{OUT_DIR}/GSEA_mito_pathways_NES', xlim=XLIM['mito'])
print("Done.")
