"""
filter_GSEA_pathways_OSD562.py
------------------------------
Pull the HARMONIZED panel sets (identical to the OSD-564 panels) out of the
OSD-562 GSEA result tables (re-run DESeq2, 4 vs 3):

  * KEGG sEV (7), KEGG TNT (3), KEGG metabolic adaptation (3)  -> exact KEGG names
  * GO BP mitochondrial (7)                                     -> matched by GO ID
  * extra reporting terms (GO BP / GO CC / KEGG)                -> for the report

Writes:
  KEGG_results_OSD562/GSEA_KEGG_pathway_filtered.csv
  GO_BP_results_OSD562/GSEA_GO_BP_mito_filtered.csv
  OSD562_GSEA_panel_terms_summary.csv   (everything, one table)
"""
import re
import pandas as pd

KEGG_CATEGORIES = {
    'sEV': {
        'title': 'KEGG sEV pathways',
        'pathways': ['Sphingolipid signaling pathway', 'Sphingolipid metabolism',
                     'Phospholipase D signaling pathway', 'Endocytosis',
                     'Synaptic vesicle cycle', 'Phosphatidylinositol signaling system',
                     'Lysosome'],
    },
    'TNT': {
        'title': 'KEGG TNT pathways',
        'pathways': ['Axon guidance', 'Regulation of actin cytoskeleton',
                     'Rap1 signaling pathway'],
    },
    'Metabolic': {
        'title': 'KEGG metabolic adaptation pathways (supplementary)',
        'pathways': ['AMPK signaling pathway', 'Glycolysis / Gluconeogenesis',
                     'Oxidative phosphorylation'],
    },
}

GO_BP_MITO = [  # (GO ID, display name)
    ('GO:0006119', 'Oxidative Phosphorylation'),
    ('GO:0033108', 'Mitochondrial Respiratory Chain Complex Assembly'),
    ('GO:0042775', 'Mitochondrial ATP Synthesis Coupled Electron Transport'),
    ('GO:0042776', 'Proton Motive Force-Driven Mitochondrial ATP Synthesis'),
    ('GO:0006120', 'Mitochondrial Electron Transport, NADH To Ubiquinone'),
    ('GO:0032981', 'Mitochondrial Respiratory Chain Complex I Assembly'),
    ('GO:0033617', 'Mitochondrial Respiratory Chain Complex IV Assembly'),
]

EXTRA = {
    'KEGG': ['Ribosome', 'Parkinson disease', 'Oxidative phosphorylation'],
    'GO_BP': [('GO:0007411', 'Axon Guidance'), ('GO:0007409', 'Axonogenesis'),
              ('GO:0017158', 'Regulation Of Calcium Ion-Dependent Exocytosis'),
              ('GO:0006887', 'Exocytosis'), ('GO:0045055', 'Regulated Exocytosis'),
              ('GO:0007015', 'Actin Filament Organization'),
              ('GO:0140053', 'Mitochondrial Gene Expression')],
    'GO_CC': [('GO:0055038', 'Recycling Endosome Membrane'), ('GO:0030139', 'Endocytic Vesicle'),
              ('GO:0070382', 'Exocytic Vesicle'), ('GO:0005743', 'Mitochondrial Inner Membrane')],
}

FILES = {'KEGG': 'KEGG_results_OSD562/GSEA_KEGG_results.csv',
         'GO_BP': 'GO_BP_results_OSD562/GSEA_GO_BP_results.csv',
         'GO_CC': 'GO_CC_results_OSD562/GSEA_GO_CC_results.csv'}


def load(path):
    g = pd.read_csv(path, index_col=0)
    g['Term'] = g['Term'].str.strip()
    return g


def kegg_exact(g, name):
    hit = g[g['Term'].str.lower() == name.lower()]
    return hit.iloc[0] if len(hit) else None


def go_by_id(g, go_id):
    hit = g[g['Term'].str.contains(re.escape(f'({go_id})'), regex=True)]
    return hit.iloc[0] if len(hit) else None


def rec(panel, query, row, library):
    if row is None:
        return {'Panel': panel, 'Query': query, 'Library': library, 'Term': 'NOT IN RESULTS',
                'NES': None, 'FDR q-val': None, 'NOM p-val': None, 'Lead_genes': None}
    return {'Panel': panel, 'Query': query, 'Library': library, 'Term': row['Term'],
            'NES': float(row['NES']), 'FDR q-val': float(row['FDR q-val']),
            'NOM p-val': float(row['NOM p-val']), 'Lead_genes': row['Lead_genes']}


tables = {k: load(v) for k, v in FILES.items()}
rows = []
for key, cat in KEGG_CATEGORIES.items():
    for p in cat['pathways']:
        rows.append(rec(cat['title'], p, kegg_exact(tables['KEGG'], p), 'KEGG_2019_Mouse'))
for go_id, name in GO_BP_MITO:
    rows.append(rec('GO BP mitochondrial', f'{name} ({go_id})', go_by_id(tables['GO_BP'], go_id),
                    'GO_Biological_Process_2026'))
for p in EXTRA['KEGG']:
    rows.append(rec('Report: KEGG', p, kegg_exact(tables['KEGG'], p), 'KEGG_2019_Mouse'))
for go_id, name in EXTRA['GO_BP']:
    rows.append(rec('Report: GO BP', f'{name} ({go_id})', go_by_id(tables['GO_BP'], go_id),
                    'GO_Biological_Process_2026'))
for go_id, name in EXTRA['GO_CC']:
    rows.append(rec('Report: GO CC', f'{name} ({go_id})', go_by_id(tables['GO_CC'], go_id),
                    'GO_Cellular_Component_2026'))

summary = pd.DataFrame(rows)
summary.to_csv('OSD562_GSEA_panel_terms_summary.csv', index=False)
summary[summary['Library'] == 'KEGG_2019_Mouse'].to_csv(
    'KEGG_results_OSD562/GSEA_KEGG_pathway_filtered.csv', index=False)
summary[summary['Panel'] == 'GO BP mitochondrial'].to_csv(
    'GO_BP_results_OSD562/GSEA_GO_BP_mito_filtered.csv', index=False)

pd.set_option('display.width', 250)
print(summary[['Panel', 'Query', 'NES', 'FDR q-val', 'NOM p-val']].to_string(index=False))
print('\nWrote OSD562_GSEA_panel_terms_summary.csv, KEGG_results_OSD562/GSEA_KEGG_pathway_filtered.csv, '
      'GO_BP_results_OSD562/GSEA_GO_BP_mito_filtered.csv')
