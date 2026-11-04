"""
make_qc_figure_OSD562.py
------------------------
RNA-seq QC metric grid (4 x 3) for the OSD-562 ISS-terminal 20-week samples:
spaceflight (n = 4) vs ground control (n = 3), from the GeneLab QC table
GLDS-557_rna_seq_qc_metrics_GLbulkRNAseq.csv.

Drawn by the shared srep_panels.qc_grid (identical conventions in OSD-564 / OSD-562 /
OSD-352): 183 x 120 mm, every y-axis starts at 0, Arial 8 pt, sample IDs (FY3 ... GY10)
as tick labels, 1 pt dashed black group-mean lines, one shared legend below the grid,
no title line.  Saved as 300 dpi PNG + vector PDF.
"""
import os
import warnings
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import srep_style as srep
srep.apply()
import srep_panels

warnings.filterwarnings('ignore')

# ──────────────────────────────────────────────────────────────────
# Load & annotate  -- keep only ISS-terminal (young, 20-week) samples
# ──────────────────────────────────────────────────────────────────
df = pd.read_csv('GLDS-557_rna_seq_qc_metrics_GLbulkRNAseq.csv')

df = df[df['sample'].str.contains('ISS-T_YNG', na=False)].copy()
df['condition'] = df['sample'].apply(lambda s: 'FLT' if '_FLT_' in s else 'GC')
df['label'] = df['sample'].apply(lambda s: s.split('_')[-1])   # sample ID (FY3, GY5, ...)
# FLT libraries first, then GC, each in numeric sample order (FY3, FY4, FY6, FY10, ...)
df = (df.assign(_o=df['condition'].map({'FLT': 0, 'GC': 1}),
                _n=df['label'].str.extract(r'(\d+)')[0].astype(int))
        .sort_values(['_o', '_n']).drop(columns=['_o', '_n']).reset_index(drop=True))

# ──────────────────────────────────────────────────────────────────
# Panel definitions: (column, metric name with unit, optional scale fn, integer axis)
# -- the same twelve GeneLab metrics and titles as make_qc_figure_OSD564.py
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
srep_panels.qc_grid(df, panels, 'QC_figures/OSD562_QC_summary')
