"""
OSD-352 (RR-3 whole-brain bulk RNA-seq) RSeQC read-distribution QC summary,
FLT (n = 3) vs GC (n = 3).

Drawn by the shared srep_panels.qc_grid (identical conventions in OSD-564 / OSD-562 /
OSD-352): 183 x 120 mm, 4 x 3 grid (11 metric panels, 12th slot empty), every y-axis
starts at 0, Arial 8 pt, sample IDs (F1, F2, F7, G7, G8, G9) as tick labels, 1 pt dashed
black group-mean lines, one shared legend below the grid, no title line.
Saved as 300 dpi PNG + vector PDF.
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
# Load & annotate
# ──────────────────────────────────────────────────────────────────
df = pd.read_csv(
    'read_dist_multiqc_report/read_dist_multiqc_data/multiqc_rseqc_read_distribution.txt',
    sep='\t'
)
df.rename(columns={'Sample': 'sample'}, inplace=True)

df['condition'] = df['sample'].apply(lambda s: 'FLT' if '_FLT_' in s else 'GC')
# e.g. RR3_BRN_FLT_F1_read_dist -> F1 (short sample ID, as in the PCA panel)
df['label'] = df['sample'].apply(lambda s: s.split('_')[-3])

# FLT first, then GC, each in numeric sample order
df = (df.assign(_o=df['condition'].map({'FLT': 0, 'GC': 1}),
                _n=df['label'].str.extract(r'(\d+)')[0].astype(int))
        .sort_values(['_o', '_n']).drop(columns=['_o', '_n']).reset_index(drop=True))
print("Samples:", df['label'].tolist())

# ──────────────────────────────────────────────────────────────────
# Panel definitions: (column, metric name with unit, optional scale fn, integer axis)
# ──────────────────────────────────────────────────────────────────
panels = [
    ('cds_exons_tag_pct',        'CDS exons (%)',             None, False),
    ('5_utr_exons_tag_pct',      "5′ UTR exons (%)",          None, False),
    ('3_utr_exons_tag_pct',      "3′ UTR exons (%)",          None, False),
    ('introns_tag_pct',          'Introns (%)',               None, False),
    ('tss_up_1kb_tag_pct',       'TSS upstream 1 kb (%)',     None, False),
    ('tss_up_5kb_tag_pct',       'TSS upstream 5 kb (%)',     None, False),
    ('tss_up_10kb_tag_pct',      'TSS upstream 10 kb (%)',    None, False),
    ('tes_down_1kb_tag_pct',     'TES downstream 1 kb (%)',   None, False),
    ('tes_down_5kb_tag_pct',     'TES downstream 5 kb (%)',   None, False),
    ('tes_down_10kb_tag_pct',    'TES downstream 10 kb (%)',  None, False),
    ('other_intergenic_tag_pct', 'Other intergenic (%)',      None, False),
]

os.makedirs('QC_figures', exist_ok=True)
srep_panels.qc_grid(df, panels, 'QC_figures/OSD352_read_dist_QC_summary')
