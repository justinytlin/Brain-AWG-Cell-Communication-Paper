"""
OSD-564 RNA-seq QC summary (12 metrics, 10 libraries), Scientific Reports style.

183 x 88 mm, 4 columns x 3 rows (panel a of Supplementary Figure 1; 10 sample IDs need the wider columns) of bar panels drawn by the shared srep_panels.qc_grid
(identical conventions in OSD-564 / OSD-562 / OSD-352): every y-axis starts at 0,
8 pt Arial throughout, sample IDs (F1 ... F9, G1 ... G9) as tick labels, 1 pt dashed
black group-mean lines, one shared legend below the grid, no title line.
Written as 300 dpi PNG + vector PDF.

Figure-legend text (not drawn): GC libraries show lower unique mapping and a
higher rRNA fraction than FLT; dashed lines = group means.
"""
import os
import warnings
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import srep_style as srep; srep.apply()
import srep_panels

warnings.filterwarnings('ignore')

# ──────────────────────────────────────────────────────────────────
# Load & annotate
# ──────────────────────────────────────────────────────────────────
df = pd.read_csv('GLDS-569_rna_seq_qc_metrics_GLbulkRNAseq.csv')

df['condition'] = df['sample'].apply(lambda s: 'FLT' if '_FLT_' in s else 'GC')
# short sample IDs: F1 ... F9 (FLT), G1 ... G9 (GC)
df['label'] = df['sample'].apply(lambda s: s.split('_')[-1])
# FLT libraries first, then GC, each in numeric sample order
df = (df.assign(_o=df['condition'].map({'FLT': 0, 'GC': 1}),
                _n=df['label'].str.extract(r'(\d+)')[0].astype(int))
        .sort_values(['_o', '_n']).drop(columns=['_o', '_n']).reset_index(drop=True))

# ──────────────────────────────────────────────────────────────────
# Panel definitions: (column, metric name with unit, optional scale fn, integer axis)
# -- the same twelve GeneLab metrics and titles as make_qc_figure_OSD562.py
# ──────────────────────────────────────────────────────────────────
panels = [
    ('read_depth',                'Total reads (millions)',          lambda x: x / 1e6, False),
    ('rin',                       'RNA integrity number (RIN)',      None, False),
    ('rrna_contamination',        'rRNA contamination (%)',          None, False),
    ('raw_percent_duplicates_f',  'Duplicate reads, R1 (%)',         None, False),
    ('gene_detected_gt10',        'Genes detected (> 10 counts)',    None, True),
    ('uniquely_mapped_percent',   'Uniquely mapped reads (%)',       None, False),
    ('multimapped_percent',       'Multi-mapped reads (%)',          None, False),
    ('unmapped_tooshort_percent', 'Unmapped, too short (%)',         None, False),
    ('mean_genebody_cov_5_20',    "5′ gene-body coverage (%)",       None, False),
    ('mean_genebody_cov_80_95',   "3′ gene-body coverage (%)",       None, False),
    ('ratio_genebody_cov_3_to_5', "Gene-body 3′/5′ ratio",           None, False),
    ('pct_sense',                 'Sense-strand reads (%)',          None, False),
]

os.makedirs('QC_figures', exist_ok=True)
srep_panels.qc_grid(df, panels, 'QC_figures/OSD564_QC_summary',
                    n_cols=4, n_rows=3, size_key='qc_4x3')
