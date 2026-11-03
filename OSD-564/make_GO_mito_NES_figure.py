"""
GO BP GSEA, mitochondrial terms -- harmonized 7-term panel (shared with OSD-562).
Terms are matched by GO ID in GO_results/GSEA_BP_results.csv.
"""
import os
import warnings
import pandas as pd

from osd564_nes_panel import plot_nes_panel

warnings.filterwarnings('ignore')

# ──────────────────────────────────────────────────────────────────
# Configuration: harmonized GO BP mitochondrial panel (7 terms)
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

GSEA_CSV = 'GO_results/GSEA_BP_results.csv'
OUT_DIR  = 'GO_results'
os.makedirs(OUT_DIR, exist_ok=True)

# ──────────────────────────────────────────────────────────────────
# 1. Load GSEA results and match by GO ID
# ──────────────────────────────────────────────────────────────────
print("Loading GSEA BP results ...")
gsea_df = pd.read_csv(GSEA_CSV, index_col=0)
gsea_df['Term'] = gsea_df['Term'].str.strip()

records = []
for go_id in TARGET_GO_IDS:
    hits = gsea_df[gsea_df['Term'].str.contains(f'({go_id})', regex=False)]
    if hits.empty:
        print(f"  [MISSING] {go_id}")
        continue
    row = hits.iloc[0]
    records.append({
        'label': row['Term'],
        'NES':   float(row['NES']),
        'FDR':   float(row['FDR q-val']),
    })
    print(f"  {row['Term']:<70s} NES={row['NES']:+.3f}  FDR={row['FDR q-val']:.4f}")

# ──────────────────────────────────────────────────────────────────
# 2. Plot
# ──────────────────────────────────────────────────────────────────
out_path = f'{OUT_DIR}/GSEA_mito_pathways_NES.png'
plot_nes_panel(records, 'GO BP GSEA, mitochondrial terms', out_path,
               fig_w=11.0, label_fontsize=10.5)
print("Done.")
