"""
make_GO_mito_NES_figure_OSD562.py
---------------------------------
GO BP mitochondrial panel (7 terms, HARMONIZED with OSD-564, matched by GO
ID) from GO_BP_results_OSD562/GSEA_GO_BP_results.csv (GO_Biological_Process_2026,
re-run DESeq2 ISS-T FLT n = 4 vs ISS-T GC n = 3).  Every term is drawn with
its FDR label whether or not significant (asterisk at FDR <= 0.05).
"""
import os
import re
import warnings
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings('ignore')
plt.rcParams.update({
    'font.family': 'Arial',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300,
})

TARGET_TERMS = [  # (GO ID, display label)
    ('GO:0006119', 'Oxidative Phosphorylation'),
    ('GO:0033108', 'Mitochondrial Respiratory Chain Complex Assembly'),
    ('GO:0042775', 'Mitochondrial ATP Synthesis Coupled Electron Transport'),
    ('GO:0042776', 'Proton Motive Force-Driven Mitochondrial ATP Synthesis'),
    ('GO:0006120', 'Mitochondrial Electron Transport, NADH To Ubiquinone'),
    ('GO:0032981', 'Mitochondrial Respiratory Chain Complex I Assembly'),
    ('GO:0033617', 'Mitochondrial Respiratory Chain Complex IV Assembly'),
]
SUBTITLE = 'OSD-562: spaceflight vs ground control (n = 4 vs 3)'
PANEL    = 'GO BP mitochondrial pathways'
GSEA_CSV = 'GO_BP_results_OSD562/GSEA_GO_BP_results.csv'
OUT_DIR  = 'GO_BP_results_OSD562'
RED, BLUE = '#C62828', '#1565C0'
os.makedirs(OUT_DIR, exist_ok=True)


def fdr_label(q):
    return 'FDR < 0.001' if q < 0.001 else f'FDR = {q:.3f}'


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
    print(f"  [OK] {go_id} {label}  →  NES={float(row['NES']):.3f}  FDR={float(row['FDR q-val']):.4f}")

pd.DataFrame(records).to_csv(f'{OUT_DIR}/GSEA_Mito_NES_summary.csv', index=False)
plot_df = pd.DataFrame(records).sort_values('NES').reset_index(drop=True)
n = len(plot_df)

legend_handles = [
    mpatches.Patch(color=RED,  alpha=1.0,  label='Enriched in spaceflight (NES > 0)'),
    mpatches.Patch(color=BLUE, alpha=1.0,  label='Enriched in ground control (NES < 0)'),
    mpatches.Patch(color='grey', alpha=1.0,  label='FDR <= 0.05 (*)'),
    mpatches.Patch(color='grey', alpha=0.45, label='FDR > 0.05'),
]

fig, ax = plt.subplots(figsize=(11, max(3.2, n * 0.6 + 2.0)))
for i, row in plot_df.iterrows():
    sig = row['FDR'] <= 0.05
    c = RED if row['NES'] > 0 else BLUE
    ax.barh(i, row['NES'], color=c, alpha=1.0 if sig else 0.45, edgecolor='white', height=0.68)
    lab = fdr_label(row['FDR']) + (' *' if sig else '')
    off = 0.05 if row['NES'] >= 0 else -0.05
    ax.text(row['NES'] + off, i, lab, va='center',
            ha='left' if row['NES'] >= 0 else 'right', fontsize=9, color='black')

ax.set_yticks(range(n))
ax.set_yticklabels(plot_df['label'], fontsize=10.5)
ax.axvline(0, color='#222222', lw=1.2)
lim = max(1.0, plot_df['NES'].abs().max()) * 1.6
ax.set_xlim(-lim, lim)
ax.xaxis.grid(True, color='#ebebeb', lw=0.8, zorder=0)
ax.set_axisbelow(True)
ax.set_xlabel('Normalized enrichment score (NES)', fontsize=11)
ax.set_title(f'{PANEL}, {SUBTITLE}', fontsize=12, fontweight='bold', pad=12)
ax.legend(handles=legend_handles, fontsize=9, frameon=True, framealpha=0.92,
          edgecolor='#cccccc', loc='lower right', labelspacing=0.5)
plt.tight_layout()
out_path = f'{OUT_DIR}/GSEA_Mito_NES_summary.png'
plt.savefig(out_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"\nSaved → {out_path}")
