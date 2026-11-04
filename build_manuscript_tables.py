#!/usr/bin/env python3
"""Build manuscript Tables 1-6 (Word + Markdown) from the regenerated result files."""
import json, math, sys
import numpy as np
from docx import Document
from docx.shared import Pt, Mm
from docx.enum.text import WD_ALIGN_PARAGRAPH

P = './table_inputs/'      # JSON extracts of the DE, GSEA and g:Profiler values (table2/3/4.json)
OUT = './Manuscript_Tables'  # writes OUT.docx and OUT.md
t2 = json.load(open(P + 'table2.json'))
t3 = json.load(open(P + 'table3.json'))
t4 = json.load(open(P + 'table4.json'))

SUP = str.maketrans('0123456789-', '⁰¹²³⁴⁵⁶⁷⁸⁹⁻')


def fmt_p(p):
    if p < 1e-3:
        m, e = f'{p:.1e}'.split('e')
        return f'{m} × 10{str(int(e)).translate(SUP)}'
    return f'{p:.3f}'


def fmt_fdr(q):
    return '< 0.001' if q < 1e-3 else f'{q:.3f}'


def fmt_mean(m, s):
    d = 0 if m >= 100 else 1                      # one precision per cell, chosen from the mean
    def f(v):
        return f'{v:,.{d}f}' if v >= 10000 else f'{v:.{d}f}'
    return f'{f(m)} ± {f(s)}'


def fmt_fc(v, nd=2):
    return f'{v:+.{nd}f}'.replace('-', '−')


# ───────────────────────────── Table 1 ─────────────────────────────
T1_HEAD = ['Dataset (mission)', 'Assay, tissue', 'Strain, sex', 'GC n', 'FLT n', 'Days in orbit', 'Age at euthanasia', 'Dose']
T1 = [
    ['OSD-564 (RR-10)', 'Bulk RNA-seq, right hippocampus', 'B6129SF2/J, female', '5', '5', '28', '18 to 19 weeks', '9.67 mGy'],
    ['OSD-562 (RRRM-2, ISS-terminal 20-week cohort)', 'Bulk RNA-seq, right hippocampus', 'C57BL/6NTac, female', '3', '4', '53 to 56', '20 weeks', 'about 17 mGy'],
    ['OSD-352 (RR-3)', 'Visium spatial transcriptomics, coronal hemisphere sections through hippocampus (Masarapu et al. 2024)', 'BALB/c, female', '3', '3', '39 to 42', 'about 18 weeks', '8.97 mGy'],
    ['OSD-352 (RR-3), supplementary', 'Bulk RNA-seq, pooled whole-brain sections', 'BALB/c, female', '3', '3', '39 to 42', 'about 18 weeks', '8.97 mGy'],
]
T1_LEG = ('Table 1. Datasets analysed. GC, ground control; FLT, spaceflight; n, animals per group included in the '
          'analyses reported here; Days in orbit and Age at euthanasia are the ranges across animals recorded in the '
          'OSDR metadata; Dose, cumulative absorbed radiation dose recorded for the mission. OSD-352 spatial statistics '
          'are taken from Masarapu et al. 2024. No statistical test is reported in this table.')

# ───────────────────────────── Table 2 ─────────────────────────────
T2_HEAD = ['Category', 'Gene', 'log2FC', 'GC mean ± SD', 'FLT mean ± SD', 'p', 'Adjusted p', 'Significance']
order = {'Mitochondria': 0, 'sEVs': 1, 'TNTs': 2}
t2s = sorted(t2, key=lambda r: (order[r['category']], abs(r['log2fc']) < 1, -r['log2fc']))
T2 = []
T2_ITALIC_COL = 1
for r in t2s:
    sig = 'Adjusted' if r['padj'] <= 0.05 else 'Nominal'
    gene = r['gene'] + ('ᵃ' if abs(r['log2fc']) < 1 else '')
    T2.append([r['category'] if r['category'] != 'Mitochondria' else 'Mitochondria', gene, fmt_fc(r['log2fc']),
               fmt_mean(r['gc_mean'], r['gc_sd']), fmt_mean(r['flt_mean'], r['flt_sd']), fmt_p(r['p']), fmt_p(r['padj']), sig])
T2_LEG = ('Table 2. Curated mitochondria-, TNT- and sEV-associated genes changed in OSD-564. log2FC, log2 fold change of '
          'spaceflight (FLT) over ground control (GC) from the GeneLab DESeq2 pipeline (Wald test; n = 5 FLT and 5 GC '
          'hippocampi); mean ± SD, DESeq2 median-of-ratios normalized counts per group (n = 5 animals); adjusted p, '
          'Benjamini-Hochberg. Rows marked Nominal have p ≤ 0.05 but adjusted p > 0.05 and are reported as '
          'hypothesis-generating only; all other rows have adjusted p ≤ 0.05 and |log2FC| ≥ 1. ᵃAdjusted p ≤ 0.05 with '
          '|log2FC| < 1, included because of their roles in ESCRT-dependent vesicle budding (Vps4a, Ist1) or energy '
          'sensing (Prkaa1). mt-Co2 is the mitochondrially encoded cytochrome c oxidase subunit II. Category assignment '
          'follows Supplementary Table S6.')

# ───────────────────────────── Table 3 ─────────────────────────────
GO_ROWS = ['Oxidative Phosphorylation', 'Mitochondrial Respiratory Chain Complex Assembly',
           'Mitochondrial ATP Synthesis Coupled Electron Transport', 'Proton Motive Force-Driven Mitochondrial ATP Synthesis',
           'Mitochondrial Electron Transport, NADH to Ubiquinone', 'Mitochondrial Respiratory Chain Complex IV Assembly',
           'Mitochondrial Respiratory Chain Complex I Assembly']
GO_NAME = {  # sentence case as in the ontology
    'Oxidative Phosphorylation': 'Oxidative phosphorylation (GO:0006119)',
    'Mitochondrial Respiratory Chain Complex Assembly': 'Mitochondrial respiratory chain complex assembly (GO:0033108)',
    'Mitochondrial ATP Synthesis Coupled Electron Transport': 'Mitochondrial ATP synthesis coupled electron transport (GO:0042775)',
    'Proton Motive Force-Driven Mitochondrial ATP Synthesis': 'Proton motive force-driven mitochondrial ATP synthesis (GO:0042776)',
    'Mitochondrial Electron Transport, NADH to Ubiquinone': 'Mitochondrial electron transport, NADH to ubiquinone (GO:0006120)',
    'Mitochondrial Respiratory Chain Complex IV Assembly': 'Mitochondrial respiratory chain complex IV assembly (GO:0033617)',
    'Mitochondrial Respiratory Chain Complex I Assembly': 'Mitochondrial respiratory chain complex I assembly (GO:0032981)',
}
for k, v in GO_NAME.items():   # keep the GO ID that the result file carries, not a typed one
    gid = t3['go564'][k]['term'].split('(')[-1].rstrip(')')
    GO_NAME[k] = v.split(' (GO:')[0] + f' ({gid})'


def direction(nes):
    return 'Spaceflight' if nes > 0 else 'Ground control'


T3A_HEAD = ['GO Biological Process term', 'Genes in set (564 / 562)', 'OSD-564 NES', 'OSD-564 FDR', 'OSD-562 NES', 'OSD-562 FDR', 'Enriched in']
T3A = []
for k in GO_ROWS:
    a, b = t3['go564'][k], t3['go562'][k]
    T3A.append([GO_NAME[k], f"{a['size']} / {b['size']}", fmt_fc(a['nes']), fmt_fdr(a['fdr']), fmt_fc(b['nes']), fmt_fdr(b['fdr']), direction(a['nes'])])

KEGG_ROWS = [('sEV', 'Sphingolipid signaling pathway'), ('sEV', 'Sphingolipid metabolism'), ('sEV', 'Lysosome'),
             ('sEV', 'Phospholipase D signaling pathway'), ('sEV', 'Endocytosis'), ('sEV', 'Synaptic vesicle cycle'),
             ('sEV', 'Phosphatidylinositol signaling system'),
             ('TNT', 'Axon guidance'), ('TNT', 'Regulation of actin cytoskeleton'), ('TNT', 'Rap1 signaling pathway'),
             ('Metabolic adaptation', 'AMPK signaling pathway'), ('Metabolic adaptation', 'Glycolysis / Gluconeogenesis'),
             ('Mitochondria', 'Oxidative phosphorylation'), ('Mitochondria', 'Ribosome'), ('Mitochondria', 'Parkinson disease')]
T3B_HEAD = ['Category', 'KEGG pathway', 'Genes in set (564 / 562)', 'OSD-564 NES', 'OSD-564 FDR', 'OSD-562 NES', 'OSD-562 FDR', 'Enriched in (OSD-564)']
T3B = []
for cat, k in KEGG_ROWS:
    a, b = t3['kegg564'][k], t3['kegg562'][k]
    T3B.append([cat, k.replace(' / ', '/'), f"{a['size']} / {b['size']}", fmt_fc(a['nes']), fmt_fdr(a['fdr']), fmt_fc(b['nes']), fmt_fdr(b['fdr']), direction(a['nes'])])
T3_LEG = ('Table 3. Gene set enrichment of curated (a) GO Biological Process mitochondrial terms and (b) KEGG pathways in '
          'OSD-564 and OSD-562. NES, normalized enrichment score from GSEApy 1.3.0 prerank, genes ranked by '
          'sign(log2FC) × −log10 nominal p, 1,000 gene-set permutations, gene-set size 10 to 500, seed 42; NES > 0, '
          'enriched in spaceflight; NES < 0, enriched in ground control (Enriched in column). FDR, permutation-based Benjamini-Hochberg q '
          'value; "< 0.001" denotes q = 0 after 1,000 permutations. NES carries no standard error; its uncertainty is '
          'given by the permutation null distribution from which FDR is derived. Genes in set, members of the term '
          'present in each ranked list. OSD-564, n = 5 vs 5 (GeneLab DESeq2 statistics; GO_Biological_Process_2026, '
          'KEGG_2019_Mouse); OSD-562, n = 4 vs 3 (re-derived DESeq2 contrast, Methods). In OSD-562 no sEV, TNT or '
          'metabolic pathway reached FDR ≤ 0.05. Terms were selected before analysis (Supplementary Table S6); complete '
          'results are in Supplementary Table S3.')

# ───────────────────────────── Table 4 ─────────────────────────────
T4_HEAD = ['Input list', 'Source', 'Term', '−log10 adj. p', 'Term size', 'Intersection', 'Query size']
pick_down = ['mitochondrial envelope', 'mitochondrial inner membrane', 'oxidative phosphorylation',
             'mitochondrial respiratory chain complex assembly', 'aerobic respiration',
             'mitochondrial ATP synthesis coupled electron transport', 'mitochondrial cytochrome c oxidase assembly',
             'proton motive force-driven mitochondrial ATP synthesis', 'mitochondrial respiratory chain complex I assembly']
pick_up = ['vesicle-mediated transport', 'cell-cell signaling', 'intracellular transport', 'regulated exocytosis',
           'small GTPase binding', 'plasma membrane bounded cell projection organization', 'actin cytoskeleton organization',
           'vesicle organization', 'extracellular organelle', 'vesicle localization', 'protein localization to plasma membrane']
pick_ribo = ['cytosolic ribosome', 'Ribosome', 'translation at synapse']
allrows = t4['down_mito'] + t4['up'] + t4['down_ribo']
allrows.append(dict(term='mitochondrial cytochrome c oxidase assembly', source='GO:BP', id='GO:0033617', nlp=3.326181, size=24, inter=7, query=496))


def find4(name, source=None):
    for r in allrows:
        if r['term'].lower() == name.lower() and (source is None or r['source'] == source):
            return r
    raise KeyError(name)


T4 = []
for n in pick_down:
    r = find4(n, 'GO:BP' if n == 'oxidative phosphorylation' else None)
    T4.append(['Down in spaceflight', r['source'], f"{r['term']} ({r['id']})", f"{r['nlp']:.2f}", str(r['size']), str(r['inter']), str(r['query'])])
for n in pick_up:
    r = find4(n)
    T4.append(['Up in spaceflight', r['source'], f"{r['term']} ({r['id']})", f"{r['nlp']:.2f}", str(r['size']), str(r['inter']), str(r['query'])])
for n in pick_ribo:
    r = find4(n, 'KEGG' if n == 'Ribosome' else None)
    T4.append(['Down in spaceflight', r['source'], f"{r['term']} ({r['id']})", f"{r['nlp']:.2f}", str(r['size']), str(r['inter']), str(r['query'])])
T4_LEG = ('Table 4. Over-representation of GO and KEGG terms among OSD-564 genes upregulated or downregulated in '
          'spaceflight (GeneLab DESeq2 adjusted p < 0.05, any fold change). Values are −log10 of the g:SCS-adjusted p from '
          'g:Profiler g:GOSt (Mus musculus, annotated-genome background; GO:BP, GO:CC, GO:MF, KEGG, Reactome, '
          'WikiPathways, CORUM). Term size, number of annotated genes in the term; Intersection, number of input genes in '
          'the term; Query size, number of input genes with an annotation in that source [these are the sizes in the '
          'current export, 391 to 502; update after re-running g:Profiler from the exact lists, Methods]. Input list '
          'indicates which gene set produced the term; the three ribosomal terms share the same cytosolic '
          'ribosomal-protein genes. Full exports are in Supplementary Table S4.')

# ───────────────────────────── Table 5 ─────────────────────────────
import pandas as pd
ca3 = pd.read_csv('OSD-352_spatial/CA3_panels/CA3_cluster8_per_animal_summary.csv')
cat5 = {'Pink1': 'Mitochondria', 'Snph': 'Mitochondria', 'Itpr1': 'Mitochondria', 'Rhob': 'sEVs', 'Sec14l1': 'sEVs', 'Ubl3': 'sEVs', 'Arpc2': 'TNTs', 'Arap2': 'TNTs'}
T5_HEAD = ['Category', 'Gene', 'Cluster', 'MAST coefficient', 'FDR', 'CDR-adjusted mean ± SD, FLT', 'CDR-adjusted mean ± SD, GC', 'Difference']
T5 = []
for g in ['Pink1', 'Snph', 'Itpr1', 'Rhob', 'Sec14l1', 'Ubl3', 'Arpc2', 'Arap2']:
    r = ca3[ca3['gene'] == g].iloc[0]
    f = np.array([r['adj_F1(CF1)'], r['adj_F2(CF2)'], r['adj_F3(CF7)']]); c = np.array([r['adj_G1(CG7)'], r['adj_G2(CG8)'], r['adj_G3(CG9)']])
    T5.append([cat5[g], g, '8 (CA3)', fmt_fc(r['MAST_coef'], 4), f"{r['MAST_FDR']:.3f}", f"{f.mean():.2f} ± {f.std(ddof=1):.2f}", f"{c.mean():.2f} ± {c.std(ddof=1):.2f}", fmt_fc(f.mean() - c.mean())])
T5_LEG = ('Table 5. Spatially resolved differential expression in hippocampal CA3 (spatial cluster 8) of OSD-352. MAST '
          'coefficient, flight-versus-ground coefficient from a MAST hurdle mixed model with sample as random effect, and '
          'FDR, Benjamini-Hochberg adjusted p, both from Masarapu et al. 2024 Supplementary Data 8 (three flight and three '
          'ground animals, two coronal sections each; spot-level model; positive, higher in spaceflight). CDR-adjusted '
          'mean ± SD, per-animal mean SCTransform log-normalized expression in cluster-8 spots (both sections per animal) '
          'after linear adjustment for cellular detection rate, mean and standard deviation over n = 3 animals per group; '
          'these values are descriptive and were not used for inference. Category assignment follows Supplementary Table S6.')

# ───────────────────────────── Table 6 ─────────────────────────────
T6_HEAD = ['Hypothesis', 'OSD-564 (28 d, hippocampus, n = 5/5)', 'OSD-562 (53 to 56 d, hippocampus, n = 4/3)', 'OSD-352 CA3 spatial (39 to 42 d, n = 3/3)', 'OSD-352 whole-brain bulk (supp., n = 3/3)']
T6 = [
    ['Mitochondrial respiration suppressed', 'Yes: 7 of 7 GO terms, KEGG oxidative phosphorylation and ORA, FDR ≤ 0.05', 'Yes: 6 of 7 GO terms and KEGG oxidative phosphorylation, FDR ≤ 0.01', 'Pink1 and Snph down, Itpr1 up (FDR < 0.05)', 'No: all terms trend toward spaceflight, n.s.'],
    ['Actin and protrusion (TNT) programs up', 'Yes: 3 KEGG pathways, GO BP terms, 8 genes at adjusted p ≤ 0.05', 'No: NES −1.08 to +0.81, n.s.', 'Arpc2 and Arap2 up (FDR < 0.05)', 'No: n.s.; Rac2 and Itga2 nominally down'],
    ['Vesicle and sEV programs up', 'Yes: 7 KEGG pathways, GO CC terms and ORA', 'No: NES −1.23 to +1.08, n.s.', 'Rhob, Sec14l1 and Ubl3 up (FDR < 0.05)', 'No: n.s.'],
]
T6_LEG = ('Table 6. Support for each hypothesis by dataset. Yes, the curated pathway panel reached FDR ≤ 0.05 (GSEA) or '
          'the named genes reached adjusted p ≤ 0.05; No, no term or gene in the panel reached that threshold; n.s., not '
          'significant. Entries summarize Tables 2 to 5 and Supplementary Tables S1 to S5; no new statistics are computed '
          'here. d, days in orbit; n, flight/ground animals.')

TEXT_MM = 176                      # portrait US Letter with 20 mm margins
W1 = [34, 44, 24, 10, 10, 14, 20, 20]
W2 = [22, 18, 14, 25, 25, 20, 24, 28]
W3A = [62, 20, 16, 16, 16, 16, 30]
W3B = [22, 40, 20, 15, 15, 15, 15, 34]
W4 = [22, 14, 76, 16, 16, 20, 12]
W5 = [20, 16, 14, 22, 14, 34, 34, 22]
W6 = [30, 38, 36, 36, 36]
TABLES = [
    ('Table 1', T1_HEAD, T1, T1_LEG, None, W1),
    ('Table 2', T2_HEAD, T2, T2_LEG, 1, W2),
    ('Table 3a', T3A_HEAD, T3A, None, None, W3A),
    ('Table 3b', T3B_HEAD, T3B, T3_LEG, None, W3B),
    ('Table 4', T4_HEAD, T4, T4_LEG, None, W4),
    ('Table 5', T5_HEAD, T5, T5_LEG, 1, W5),
    ('Table 6', T6_HEAD, T6, T6_LEG, None, W6),
]
for _, h, _, _, _, w in TABLES:
    assert len(w) == len(h) and abs(sum(w) - TEXT_MM) <= 1, h

# ───────────────────────────── Word ─────────────────────────────
doc = Document()
st = doc.styles['Normal']; st.font.name = 'Arial'; st.font.size = Pt(9)
for s in doc.sections:                                  # portrait US Letter, 20 mm margins -> 176 mm text width
    s.page_width, s.page_height = Mm(215.9), Mm(279.4)
    s.left_margin = s.right_margin = Mm(20); s.top_margin = s.bottom_margin = Mm(20)
p = doc.add_paragraph(); r = p.add_run('Manuscript tables (copy each table and its legend into the main document)'); r.bold = True
doc.add_paragraph('Values are taken directly from the regenerated result files (see STEP_BY_STEP_REVISION_GUIDE.md, Step 3). '
                  'Tables are plain Word tables without colour; gene symbols are italic. Square-bracketed text in a legend is a note to resolve before submission.')


def add_table(title, head, rows, italic_col, widths):
    p = doc.add_paragraph(); p.add_run(title).bold = True
    tb = doc.add_table(rows=1, cols=len(head)); tb.style = 'Table Grid'; tb.autofit = False
    for i, h in enumerate(head):
        c = tb.rows[0].cells[i]; c.text = ''; c.paragraphs[0].add_run(h).bold = True
    for row in rows:
        cells = tb.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ''
            if italic_col is not None and i == italic_col:
                base, mark = (v[:-1], v[-1]) if v.endswith('ᵃ') else (v, '')
                cells[i].paragraphs[0].add_run(base).italic = True
                if mark:
                    cells[i].paragraphs[0].add_run('a').font.superscript = True
            else:
                cells[i].paragraphs[0].add_run(v)
    for row in tb.rows:
        for c, w in zip(row.cells, widths):
            c.width = Mm(w)
            for para in c.paragraphs:
                para.paragraph_format.space_after = Pt(0)
                for run in para.runs:
                    run.font.size = Pt(8)


for title, head, rows, leg, ital, widths in TABLES:
    add_table(title, head, rows, ital, widths)
    if leg:
        doc.add_paragraph(leg)
    doc.add_paragraph()
doc.save(OUT + '.docx')

# ───────────────────────────── Markdown ─────────────────────────────
md = ['# Manuscript tables', '', 'Built from the regenerated result files (guide Step 3). Copy from the Word version (`Manuscript_Tables.docx`) to keep table structure; this file is for reading and diffing.', '']
for title, head, rows, leg, ital, widths in TABLES:
    md += [f'## {title}', '', '| ' + ' | '.join(head) + ' |', '|' + '---|' * len(head)]
    for row in rows:
        cells = [(f'*{v[:-1]}*ᵃ' if v.endswith('ᵃ') else f'*{v}*') if (ital is not None and i == ital) else v for i, v in enumerate(row)]
        md.append('| ' + ' | '.join(cells) + ' |')
    md.append('')
    if leg:
        md += [leg, '']
open(OUT + '.md', 'w').write('\n'.join(md))
print('wrote', OUT + '.docx / .md')
for title, head, rows, *_ in TABLES:
    print(f'{title}: {len(rows)} rows x {len(head)} cols')
