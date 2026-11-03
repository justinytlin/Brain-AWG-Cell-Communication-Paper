#!/opt/anaconda3/bin/python3
"""
rerun_deseq2_OSD562.py
----------------------
Re-run the OSD-562 / GLDS-557 (RRRM-2 hippocampus) differential
expression analysis as a clean 4 vs 3 contrast:

    FLT = RRRM2_HPC_FLT_ISS-T_YNG_{FY3, FY4, FY6, FY10}   (spaceflight,  ISS-Terminal, 20 wk)
    GC  = RRRM2_HPC_GC_ISS-T_YNG_{GY5, GY8, GY10}         (ground ctrl,  ISS-Terminal, 20 wk)

Why: the GeneLab contrast
  "(12 week & Space Flight & On ISS)v(12 week & Ground Control & On Earth)"
is 3 FLT (FY3, FY4, FY6; FY10 is mislabelled "On Earth" in the ISA) vs
6 GC (3 ISS-T + 3 LAR), so it is not the ISS-T FLT vs ISS-T GC comparison
the manuscript describes.

Steps
  1. Read the RSEM gene-level files (STAR folder), take expected_count,
     round to integer (as DESeqDataSetFromTximport does).
  2. Drop ERCC spike-ins and GENCODE vM37 rRNA / Mt_rRNA / rRNA_pseudogene
     genes (to mimic GeneLab's "rRNArm" tables), drop zero-sum rows.
  3. Run DESeq2 (R, rerun_deseq2_OSD562.R) with design ~ condition,
     GC as reference, Wald test, BH, no lfcShrink; export vst(blind=TRUE).
  4. Format a GeneLab-style DE table, annotated with SYMBOL / GENENAME /
     ENTREZID from the existing GeneLab table (same Ensembl 112 annotation).
  5. Sanity-check selected genes against the GeneLab contrast.

Outputs (this folder):
  OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3.csv
  OSD562_VST_counts_ISS-T_YNG.csv
  OSD562_raw_counts_ISS-T_YNG.csv          (integer input matrix, for the record)
  OSD562_rRNA_genes_vM37.tsv               (cached rRNA list from the GTF)
  OSD562_DESeq2_sanity_check.csv
"""
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

HERE     = os.path.dirname(os.path.abspath(__file__))
STAR_DIR = os.path.join(os.path.dirname(HERE), 'OSD-562 Hipathia', 'STAR')
GTF      = os.path.join(os.path.dirname(HERE), 'OSD-562 Hipathia', 'gencode.vM37.annotationcopy.gtf')
RSCRIPT  = '/Library/Frameworks/R.framework/Resources/bin/Rscript'
R_SCRIPT = os.path.join(HERE, 'rerun_deseq2_OSD562.R')
GENELAB_DE = os.path.join(HERE, 'GLDS-557_rna_seq_differential_expression_rRNArm_GLbulkRNAseq.csv')

FLT = ['RRRM2_HPC_FLT_ISS-T_YNG_FY3', 'RRRM2_HPC_FLT_ISS-T_YNG_FY4',
       'RRRM2_HPC_FLT_ISS-T_YNG_FY6', 'RRRM2_HPC_FLT_ISS-T_YNG_FY10']
GC  = ['RRRM2_HPC_GC_ISS-T_YNG_GY5', 'RRRM2_HPC_GC_ISS-T_YNG_GY8',
       'RRRM2_HPC_GC_ISS-T_YNG_GY10']
SAMPLES = FLT + GC

RRNA_TYPES = {'rRNA', 'Mt_rRNA', 'rRNA_pseudogene'}
RRNA_CACHE = os.path.join(HERE, 'OSD562_rRNA_genes_vM37.tsv')

OUT_DE  = os.path.join(HERE, 'OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3.csv')
OUT_VST = os.path.join(HERE, 'OSD562_VST_counts_ISS-T_YNG.csv')
OUT_RAW = os.path.join(HERE, 'OSD562_raw_counts_ISS-T_YNG.csv')
OUT_SAN = os.path.join(HERE, 'OSD562_DESeq2_sanity_check.csv')

OLD_TAG = '(12 week & Space Flight & On ISS)v(12 week & Ground Control & On Earth)'
NEW_TAG = '(Space Flight)v(Ground Control)'

os.chdir(HERE)

# ------------------------------------------------------------------
# 1. Counts matrix from RSEM expected_count
# ------------------------------------------------------------------
print('Building counts matrix from RSEM expected_count ...')
cols = {}
for s in SAMPLES:
    f = os.path.join(STAR_DIR, f'GLDS-557_rna_seq_{s}.genes.tsv')
    t = pd.read_csv(f, sep='\t', usecols=['gene_id', 'expected_count'], index_col='gene_id')
    cols[s] = t['expected_count']
counts = pd.DataFrame(cols)
counts = counts.round(0).astype(int)
print(f'  {counts.shape[0]:,} gene rows x {counts.shape[1]} samples')

# ------------------------------------------------------------------
# 2. Filters: ERCC, rRNA (GENCODE vM37), zero-sum rows
# ------------------------------------------------------------------
n0 = len(counts)
is_ercc = counts.index.str.startswith('ERCC-')
counts = counts[~is_ercc]
print(f'  dropped {is_ercc.sum()} ERCC spike-in rows')

if not os.path.exists(RRNA_CACHE):
    print(f'  parsing GTF for rRNA genes (one-off, cached to {os.path.basename(RRNA_CACHE)}) ...')
    rows = []
    pat_id   = re.compile(r'gene_id "([^"]+)"')
    pat_type = re.compile(r'gene_type "([^"]+)"')
    pat_name = re.compile(r'gene_name "([^"]+)"')
    with open(GTF) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.split('\t', 9)
            if parts[2] != 'gene':
                continue
            gtype = pat_type.search(parts[8])
            if gtype and gtype.group(1) in RRNA_TYPES:
                gid  = pat_id.search(parts[8]).group(1)
                gnm  = pat_name.search(parts[8])
                rows.append((gid, gtype.group(1), gnm.group(1) if gnm else ''))
    pd.DataFrame(rows, columns=['gene_id', 'gene_type', 'gene_name']).to_csv(
        RRNA_CACHE, sep='\t', index=False)
rrna = pd.read_csv(RRNA_CACHE, sep='\t')
rrna_ids = set(rrna['gene_id'].str.replace(r'\.\d+$', '', regex=True))
print(f'  rRNA-type genes in GTF: {len(rrna_ids)} '
      f'({rrna["gene_type"].value_counts().to_dict()})')

stripped = counts.index.str.replace(r'\.\d+$', '', regex=True)
is_rrna = stripped.isin(rrna_ids)
print(f'  dropped {is_rrna.sum()} rRNA rows present in RSEM output')
counts = counts[~is_rrna]

nonzero = counts.sum(axis=1) > 0
print(f'  dropped {(~nonzero).sum()} zero-sum rows')
counts = counts[nonzero]
print(f'  final matrix: {counts.shape[0]:,} genes x {counts.shape[1]} samples')
counts.index.name = 'ENSEMBL'
counts.to_csv(OUT_RAW)

coldata = pd.DataFrame({'condition': ['FLT'] * len(FLT) + ['GC'] * len(GC)},
                       index=SAMPLES)
coldata.index.name = 'sample'
coldata_file = os.path.join(HERE, '_deseq2_tmp_coldata.csv')
coldata.to_csv(coldata_file)

# ------------------------------------------------------------------
# 3. DESeq2 in R
# ------------------------------------------------------------------
prefix = os.path.join(HERE, '_deseq2_tmp')
print('\nRunning DESeq2 (R) ...')
cmd = [RSCRIPT, R_SCRIPT, OUT_RAW, coldata_file, prefix]
r = subprocess.run(cmd, capture_output=True, text=True)
print(r.stdout)
if r.returncode != 0:
    print(r.stderr)
    sys.exit('DESeq2 failed')

res  = pd.read_csv(prefix + '_deseq2_results_raw.csv', index_col='ENSEMBL')
norm = pd.read_csv(prefix + '_normalized_counts.csv', index_col='ENSEMBL')
vst  = pd.read_csv(prefix + '_vst.csv', index_col=0)

# ------------------------------------------------------------------
# 4. GeneLab-style table
# ------------------------------------------------------------------
print('Formatting GeneLab-style DE table ...')
gl = pd.read_csv(GENELAB_DE, index_col=0, low_memory=False)
annot = gl[['SYMBOL', 'GENENAME', 'ENTREZID']]

out = pd.DataFrame(index=counts.index)
out.index.name = 'ENSEMBL'
out = out.join(annot, how='left')
out = out.join(norm[SAMPLES], how='left')
out[f'Group.Mean_(Space Flight)']   = norm[FLT].mean(axis=1)
out[f'Group.Mean_(Ground Control)'] = norm[GC].mean(axis=1)
out[f'Group.Stdev_(Space Flight)']   = norm[FLT].std(axis=1, ddof=1)
out[f'Group.Stdev_(Ground Control)'] = norm[GC].std(axis=1, ddof=1)
out[f'Log2fc_{NEW_TAG}']      = res['log2FoldChange']
out[f'lfcSE_{NEW_TAG}']       = res['lfcSE']
out[f'Stat_{NEW_TAG}']        = res['stat']
out[f'P.value_{NEW_TAG}']     = res['pvalue']
out[f'Adj.p.value_{NEW_TAG}'] = res['padj']
out['baseMean']               = res['baseMean']
out.to_csv(OUT_DE)
print(f'  wrote {OUT_DE}  ({len(out):,} genes; '
      f'{out["SYMBOL"].notna().sum():,} with SYMBOL)')

vst.index.name = ''
vst[SAMPLES].to_csv(OUT_VST)
print(f'  wrote {OUT_VST}')

padj = out[f'Adj.p.value_{NEW_TAG}']; lfc = out[f'Log2fc_{NEW_TAG}']; p = out[f'P.value_{NEW_TAG}']
print(f'\n  adj.p<=0.05 & |log2FC|>=1 : up={((padj<=.05)&(lfc>=1)).sum()}  down={((padj<=.05)&(lfc<=-1)).sum()}')
print(f'  adj.p<=0.05 (any FC)      : {(padj<=.05).sum()}')
print(f'  nom.p<=0.05 & |log2FC|>=1 : up={((p<=.05)&(lfc>=1)).sum()}  down={((p<=.05)&(lfc<=-1)).sum()}')
print(f'  nom.p<=0.05 (any FC)      : {(p<=.05).sum()}')

# ------------------------------------------------------------------
# 5. Sanity check vs GeneLab contrast
# ------------------------------------------------------------------
check = ['Cpt1b', 'Pprc1', 'Sdc4', 'Rhof', 'Ghsr', 'Hmga1b', 'Ttn', 'Smarcd1',
         'Vps4a', 'Ist1', 'Chmp1a', 'Chmp6', 'Vta1', 'Rac2', 'Baiap2l2', 'Prkaa1']
rows = []
for g in check:
    ids = gl.index[gl['SYMBOL'] == g]
    for gid in ids:
        rows.append({
            'SYMBOL': g, 'ENSEMBL': gid,
            'old_log2FC': gl.at[gid, f'Log2fc_{OLD_TAG}'],
            'old_p':      gl.at[gid, f'P.value_{OLD_TAG}'],
            'old_padj':   gl.at[gid, f'Adj.p.value_{OLD_TAG}'],
            'new_log2FC': out[f'Log2fc_{NEW_TAG}'].get(gid, np.nan),
            'new_p':      out[f'P.value_{NEW_TAG}'].get(gid, np.nan),
            'new_padj':   out[f'Adj.p.value_{NEW_TAG}'].get(gid, np.nan),
        })
san = pd.DataFrame(rows)
san['same_direction'] = np.sign(san['old_log2FC']) == np.sign(san['new_log2FC'])
san.to_csv(OUT_SAN, index=False)
pd.set_option('display.width', 200)
print('\nSanity check (GeneLab 3v6 contrast vs new 4v3):')
print(san.to_string(index=False))

for f in [coldata_file, prefix + '_deseq2_results_raw.csv', prefix + '_normalized_counts.csv',
          prefix + '_vst.csv', prefix + '_size_factors.csv']:
    if os.path.exists(f):
        os.remove(f)
print('\nDone.')
