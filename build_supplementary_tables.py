#!/usr/bin/env python3
"""
Build the manuscript's main Tables 1-6 (Excel export of the tables in the revised manuscript)
and Supplementary Tables S1-S5, S6 (draft) and S7 from the result files.

S4 (g:Profiler ORA) uses the existing exports in the OSD-564 gprofiler/ folder. Their titles say "Raw DEGs
p<0.05", but every intersection gene has adjusted p < 0.05 (any fold change), so the input matches
the Methods; the title text is a mislabel.

Run from the project root:  /opt/anaconda3/bin/python3 build_supplementary_tables.py
Outputs go to Tables/.
"""
import os, re
import numpy as np
import pandas as pd
import docx
from openpyxl import load_workbook
from openpyxl.styles import Font, Alignment

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "Tables")
os.makedirs(OUT, exist_ok=True)
P564 = os.path.join(ROOT, "OSD-564 Figures")
P562 = os.path.join(ROOT, "OSD-562 Figures")
P352 = os.path.join(ROOT, "OSD-352", "NEW FIGURES")
PST = os.path.join(ROOT, "OSD-352", "NASA_RR3_Brain-main", "manuscript_plots", "CA3_panels")
MANUSCRIPT_NAME = "Paper REVISED v4 (clean).docx"
MANUSCRIPT = os.path.join(ROOT, MANUSCRIPT_NAME)
MASARAPU_SD8 = os.path.expanduser("~/Downloads/41467_2024_48916_MOESM11_ESM.xlsx")

FC_THR, P_THR = 1.0, 0.05

# ──────────────────────────────────────────────────────────────────
# Category of each curated gene, as stated in the manuscript (Tables 2 and 5 and Results text).
# Genes not listed here have no category in the manuscript ("Unassigned"; see S6 draft).
# ──────────────────────────────────────────────────────────────────
CATEGORY = {}
for g in ["Mrps25", "Pet100", "Chchd1", "Ndufb4", "Atp5mk", "Ndufaf6", "mt-Co2", "Rmrp", "Prkaa1",
          "Pink1", "Snph", "Itpr1", "Cpt1b", "Hpdl", "Bbox1", "Pask"]:
    CATEGORY[g] = "Mitochondria"
for g in ["Epcam", "Casr", "Ghsr", "Wnt7a", "Folr1", "Lpl", "Bcat1", "Aqp1", "Apoc1", "Timp1", "Zfas1",
          "C4a", "S100a8", "Vps4a", "Ist1", "Rhob", "Sec14l1", "Ubl3"]:
    CATEGORY[g] = "sEV"
for g in ["Baiap2l2", "Rac2", "Adam8", "Myh7b", "Gas2l3", "Itga2", "Ect2", "Smtnl2", "Rhof", "Gas2",
          "S100a4", "Tnni2", "Wdr1", "Tmsb4x", "Arpc2", "Arap2", "Rapgef5"]:
    CATEGORY[g] = "TNT"
cat = lambda g: CATEGORY.get(g, "Unassigned")


def sig_call(fc, p, q):
    if pd.notna(q) and q <= P_THR and abs(fc) >= FC_THR:
        return "Adjusted"
    if pd.notna(p) and p <= P_THR and abs(fc) >= FC_THR:
        return "Nominal"
    return "n.s."


def readme_df(rows):
    return pd.DataFrame(rows, columns=["Item", "Description"])


def write_book(path, sheets, wide_cols=None):
    """sheets: list of (name, DataFrame). README sheets get wrapped text."""
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        for name, df in sheets:
            df.to_excel(xw, sheet_name=name, index=False)
    wb = load_workbook(path)
    for ws in wb.worksheets:
        ws.freeze_panes = "A2"
        for c in ws[1]:
            c.font = Font(bold=True)
        for col in ws.columns:
            letter = col[0].column_letter
            n = max(len(str(c.value)) if c.value is not None else 0 for c in col[:200])
            ws.column_dimensions[letter].width = min(max(10, n + 2), 60)
        if ws.title.lower().startswith("readme") or ws.title.lower().startswith("legend"):
            ws.column_dimensions["A"].width = 28
            ws.column_dimensions["B"].width = 120
            for row in ws.iter_rows(min_row=2):
                for c in row:
                    c.alignment = Alignment(wrap_text=True, vertical="top")
    wb.save(path)
    print(f"wrote {os.path.relpath(path, ROOT)}: " + ", ".join(f"{n} ({len(d)})" for n, d in sheets))


# ══════════════════════════════════════════════════════════════════
# Main Tables 1-6, exported from the revised manuscript
# ══════════════════════════════════════════════════════════════════
def cell_text(cell):
    out = []
    for p in cell.paragraphs:
        s = ""
        for r in p.runs:
            t = r.text
            if r.font.superscript and t.isalpha():
                t = "".join({"a": "ᵃ", "b": "ᵇ", "c": "ᶜ"}.get(ch, ch) for ch in t)
            s += t
        out.append(s)
    return "\n".join(out).strip()


def main_tables():
    d = docx.Document(MANUSCRIPT)
    legends = [p.text for p in d.paragraphs if re.match(r"Table \d\. ", p.text)]
    names = ["Table 1", "Table 2", "Table 3a", "Table 3b", "Table 4", "Table 5", "Table 6"]
    assert len(d.tables) == len(names) and len(legends) == 6, (len(d.tables), len(legends))
    sheets = []
    for name, t in zip(names, d.tables):
        rows = [[cell_text(c) for c in r.cells] for r in t.rows]
        sheets.append((name, pd.DataFrame(rows[1:], columns=rows[0])))
    leg_rows = [(f"Table {i + 1}", leg) for i, leg in enumerate(legends)]
    leg_rows.append(("Note", f"Exported from '{MANUSCRIPT_NAME}'."))
    write_book(os.path.join(OUT, "Main_Tables_1-6.xlsx"), [("Legends", readme_df(leg_rows))] + sheets)


# ══════════════════════════════════════════════════════════════════
# S1  OSD-562 differential expression (re-derived 4 FLT vs 3 GC)
# ══════════════════════════════════════════════════════════════════
def table_s1():
    tag = "(Space Flight)v(Ground Control)"
    de = pd.read_csv(os.path.join(P562, "OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3.csv"), low_memory=False)
    cur = pd.read_csv(os.path.join(P562, "OSD562_Volcano_custom_genes_table.csv"))
    m = cur[["Gene", "ENSEMBL"]].merge(de, on="ENSEMBL", how="left")
    s = pd.DataFrame({
        "Gene": m["Gene"], "Ensembl ID": m["ENSEMBL"], "Category": m["Gene"].map(cat),
        "GC mean": m["Group.Mean_(Ground Control)"], "GC SD": m["Group.Stdev_(Ground Control)"],
        "FLT mean": m["Group.Mean_(Space Flight)"], "FLT SD": m["Group.Stdev_(Space Flight)"],
        "log2FC (FLT/GC)": m[f"Log2fc_{tag}"], "lfcSE": m[f"lfcSE_{tag}"], "Wald statistic": m[f"Stat_{tag}"],
        "p": m[f"P.value_{tag}"], "Adjusted p": m[f"Adj.p.value_{tag}"],
    })
    s["Significance"] = [sig_call(a, b, c) for a, b, c in zip(s["log2FC (FLT/GC)"], s["p"], s["Adjusted p"])]
    s = s.sort_values(["Category", "Gene"]).reset_index(drop=True)

    allg = de.rename(columns=lambda c: c.replace(f"_{tag}", "").replace("_(Space Flight)", " FLT")
                     .replace("_(Ground Control)", " GC"))
    readme = readme_df([
        ("Dataset", "OSD-562 (GLDS-557), Rodent Research Reference Mission-2, right hippocampus, ISS-terminal "
                    "young (20-week) cohort; C57BL/6NTac females."),
        ("Samples", "Spaceflight (FLT): FY3, FY4, FY6, FY10. Ground control (GC): GY5, GY8, GY10."),
        ("Input", "GeneLab RSEM unnormalized counts with rRNA removed "
                  "(GLDS-557_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv), rounded to integers, "
                  "ERCC spike-ins and zero-sum rows removed."),
        ("Model", "DESeq2 1.52.0 (R 4.6.1), design ~ condition, GC as reference, Wald test, Benjamini-Hochberg "
                  "adjustment, no log2FC shrinkage. Script: OSD-562/rerun_deseq2_OSD562.py (+ .R)."),
        ("Why re-derived", "The repository's processed contrast compared 3 FLT with 6 GC animals (omitting FY10 and "
                           "pooling live-animal-return with ISS-terminal controls). See OSD-562/README_OSD562_rerun.md."),
        ("Means and SD", "DESeq2 size-factor-normalized counts, per group."),
        ("Significance", "Adjusted: adjusted p <= 0.05 and |log2FC| >= 1. Nominal: p <= 0.05 and |log2FC| >= 1. "
                         "No gene reached adjusted significance."),
        ("Curated genes", "The 45 genes marked on the OSD-562 volcano plot (Supplementary Figure 2a). Category as "
                          "assigned in the manuscript; 'Unassigned' genes are OSD-562 genes of interest without a "
                          "category (see Supplementary Table S6)."),
        ("All genes", f"Complete DESeq2 result, {len(de):,} rows, including per-sample normalized counts."),
    ])
    write_book(os.path.join(OUT, "Supplementary_Table_S1_OSD-562_DE.xlsx"),
               [("README", readme), ("Curated genes", s), ("All genes", allg)])
    return s


# ══════════════════════════════════════════════════════════════════
# S2  OSD-352 whole-brain bulk RNA-seq (GeneLab DESeq2, 3 FLT vs 3 GC)
# ══════════════════════════════════════════════════════════════════
OSD352_CURATED = ["Rac2", "Myo5c", "Psd4", "Epha2", "Itga2", "Rab43", "Ap1m2", "Cd209a", "Cd209g", "Cd300lf",
                  "Klk8", "Mapk12", "Adora2a", "Nos2", "Agt", "Cpt1b", "Gcat", "Hpdl", "Bbox1", "Osgin1", "Pask",
                  "Aqp9", "Sdc4", "Tcf7l1"]


def table_s2():
    tag = "(Space Flight)v(Ground Control)"
    de = pd.read_csv(os.path.join(P352, "GLDS-352_rna_seq_differential_expression.csv"), low_memory=False)
    nc = pd.read_csv(os.path.join(P352, "GLDS-352_rna_seq_Normalized_Counts.csv"), index_col=0)
    de["mean_norm"] = nc.mean(axis=1).reindex(de["ENSEMBL"]).values
    keep = de.dropna(subset=[f"P.value_{tag}"])
    dedup = (keep.dropna(subset=["SYMBOL"]).sort_values(["SYMBOL", "mean_norm"], ascending=[True, False])
             .drop_duplicates("SYMBOL"))
    m = dedup.set_index("SYMBOL").reindex(OSD352_CURATED)
    missing = [g for g in OSD352_CURATED if pd.isna(m.loc[g, "ENSEMBL"])]
    assert not missing, missing
    s = pd.DataFrame({
        "Gene": OSD352_CURATED, "Ensembl ID": m["ENSEMBL"].values, "Category": [cat(g) for g in OSD352_CURATED],
        "GC mean": m["Group.Mean_(Ground Control)"].values, "GC SD": m["Group.Stdev_(Ground Control)"].values,
        "FLT mean": m["Group.Mean_(Space Flight)"].values, "FLT SD": m["Group.Stdev_(Space Flight)"].values,
        "log2FC (FLT/GC)": m[f"Log2fc_{tag}"].values, "Wald statistic": m[f"Stat_{tag}"].values,
        "p": m[f"P.value_{tag}"].values, "Adjusted p": m[f"Adj.p.value_{tag}"].values,
    })
    s["Significance"] = [sig_call(a, b, c) for a, b, c in zip(s["log2FC (FLT/GC)"], s["p"], s["Adjusted p"])]
    s.loc[(s["Adjusted p"] <= P_THR) & (s["Significance"] != "Adjusted"), "Significance"] = "Adjusted (|log2FC| < 1)"
    s = s.sort_values(["Category", "Gene"]).reset_index(drop=True)

    cols = ["ENSEMBL", "SYMBOL", "GENENAME", "ENTREZID"] + [c for c in de.columns if c.startswith("RR3_BRN_")] + \
           [f"Log2fc_{tag}", f"Stat_{tag}", f"P.value_{tag}", f"Adj.p.value_{tag}",
            "Group.Mean_(Ground Control)", "Group.Mean_(Space Flight)",
            "Group.Stdev_(Ground Control)", "Group.Stdev_(Space Flight)"]
    allg = de[cols].rename(columns=lambda c: c.replace(f"_{tag}", "").replace("_(Space Flight)", " FLT")
                           .replace("_(Ground Control)", " GC"))
    n_dup = keep.dropna(subset=["SYMBOL"])["SYMBOL"].duplicated(keep=False)
    readme = readme_df([
        ("Dataset", "OSD-352 (GLDS-352), Rodent Research-3, BALB/c females; bulk RNA-seq of pooled whole-brain "
                    "sections (not dissected hippocampus). Reported in the supplement only."),
        ("Samples", "FLT: F1, F2, F7. GC: G7, G8, G9."),
        ("Source", "GeneLab processed differential expression (GLDS-352_rna_seq_differential_expression.csv) and "
                   "normalized counts, GeneLab bulk RNA-seq pipeline (DESeq2). Signs are spaceflight over ground "
                   "control."),
        ("Duplicate symbols", f"{keep.dropna(subset=['SYMBOL']).loc[n_dup, 'SYMBOL'].nunique()} gene symbols map to "
                              "more than one Ensembl ID; for the curated sheet the ID with the highest mean "
                              "normalized count is kept. The All genes sheet keeps every row."),
        ("Significance", "Adjusted: adjusted p <= 0.05 and |log2FC| >= 1. Nominal: p <= 0.05 and |log2FC| >= 1. "
                         "Tcf7l1 is the only gene at adjusted p <= 0.05."),
        ("Curated genes", "The 24 genes marked on the OSD-352 volcano plot (Supplementary Figure 2c)."),
        ("All genes", f"Complete GeneLab result, {len(de):,} rows ({len(keep):,} with statistics)."),
    ])
    write_book(os.path.join(OUT, "Supplementary_Table_S2_OSD-352_bulk_DE.xlsx"),
               [("README", readme), ("Curated genes", s), ("All genes", allg)])
    return s


# ══════════════════════════════════════════════════════════════════
# S3  Complete GSEA results
# ══════════════════════════════════════════════════════════════════
GO_MITO = ["GO:0006119", "GO:0033108", "GO:0042775", "GO:0042776", "GO:0006120", "GO:0033617", "GO:0032981"]
KEGG_SEV = ["Sphingolipid signaling pathway", "Sphingolipid metabolism", "Lysosome",
            "Phospholipase D signaling pathway", "Endocytosis", "Synaptic vesicle cycle",
            "Phosphatidylinositol signaling system"]
KEGG_TNT = ["Axon guidance", "Regulation of actin cytoskeleton", "Rap1 signaling pathway"]
KEGG_MET = ["AMPK signaling pathway", "Glycolysis / Gluconeogenesis", "Oxidative phosphorylation"]
KEGG_T3_OTHER = ["Ribosome", "Parkinson disease"]
GO_BP_PROTR = ["Actin Filament Organization", "Axonogenesis", "Axon Guidance", "Neuron Projection Morphogenesis",
               "Regulated Exocytosis", "Synaptic Vesicle Recycling"]
GO_CC_COMP = ["Recycling Endosome Membrane", "Endocytic Vesicle", "Clathrin-Coated Vesicle", "Exocytic Vesicle",
              "Mitochondrial Inner Membrane", "Ribosome"]

norm = lambda s: re.sub(r"\s+", " ", re.sub(r"\s*/\s*", " / ", str(s))).strip().lower()


def shown_in(dataset, library, term):
    gid = re.search(r"\(GO:\d+\)", term)
    gid = gid.group(0)[1:-1] if gid else None
    name = norm(re.sub(r"\s*\(GO:\d+\)\s*$", "", term))
    tags = []
    if library == "GO BP":
        if gid in GO_MITO:
            tags += {"OSD-564": ["Fig. 2a", "Table 3a"], "OSD-562": ["Fig. 2c", "Table 3a"],
                     "OSD-352 bulk": ["Supp. Fig. 3c"]}[dataset]
        if dataset == "OSD-564" and name in map(norm, GO_BP_PROTR):
            tags.append("Supp. Fig. 4b")
    elif library == "GO CC":
        if dataset == "OSD-564" and name in map(norm, GO_CC_COMP):
            tags.append("Supp. Fig. 4c")
    elif library == "KEGG":
        t3 = dataset in ("OSD-564", "OSD-562")
        if name in map(norm, KEGG_SEV):
            tags += {"OSD-564": ["Fig. 3a"], "OSD-562": ["Supp. Fig. 5b"], "OSD-352 bulk": ["Supp. Fig. 3b"]}[dataset]
        if name in map(norm, KEGG_TNT):
            tags += {"OSD-564": ["Fig. 4a"], "OSD-562": ["Supp. Fig. 5a"], "OSD-352 bulk": ["Supp. Fig. 3a"]}[dataset]
        if name in map(norm, KEGG_MET):
            tags += {"OSD-564": ["Supp. Fig. 4a"], "OSD-562": ["Supp. Fig. 5c"], "OSD-352 bulk": []}[dataset]
        if t3 and name in map(norm, KEGG_SEV + KEGG_TNT + KEGG_MET + KEGG_T3_OTHER):
            tags.append("Table 3b")
    return "; ".join(tags)


def table_s3():
    files = [
        ("OSD-564", "KEGG", f"{P564}/KEGG_results/GSEA_KEGG_results.csv"),
        ("OSD-564", "GO BP", f"{P564}/GO_results/GSEA_BP_results.csv"),
        ("OSD-564", "GO CC", f"{P564}/GO_results/GSEA_CC_results.csv"),
        ("OSD-562", "KEGG", f"{P562}/KEGG_results_OSD562/GSEA_KEGG_results.csv"),
        ("OSD-562", "GO BP", f"{P562}/GO_BP_results_OSD562/GSEA_GO_BP_results.csv"),
        ("OSD-562", "GO CC", f"{P562}/GO_CC_results_OSD562/GSEA_GO_CC_results.csv"),
        ("OSD-352 bulk", "KEGG", f"{P352}/KEGG_results_OSD352/GSEA_KEGG_results.csv"),
        ("OSD-352 bulk", "GO BP", f"{P352}/GO_BP_results_OSD352/GSEA_GO_BP_results.csv"),
    ]
    sheets, counts = [], []
    for ds, lib, f in files:
        r = pd.read_csv(f)
        out = pd.DataFrame({
            "Term": r["Term"].str.replace(r"\s*\(GO:\d+\)\s*$", "", regex=True),
            "Term ID": r["Term"].str.extract(r"\((GO:\d+)\)")[0],
            "Genes in set (ranked list)": r["Tag %"].str.split("/").str[1].astype(int),
            "ES": r["ES"], "NES": r["NES"], "Nominal p": r["NOM p-val"], "FDR q": r["FDR q-val"],
            "FWER p": r["FWER p-val"], "Tag %": r["Tag %"], "Gene %": r["Gene %"],
            "Leading-edge genes": r["Lead_genes"],
            "Shown in": [shown_in(ds, lib, t) for t in r["Term"]],
        })
        if lib == "KEGG":
            out = out.drop(columns="Term ID")
        out = out.sort_values(["FDR q", "NES"], ascending=[True, False]).reset_index(drop=True)
        name = f"{ds.replace(' bulk', '')} {lib}"
        sheets.append((name, out))
        counts.append(f"{name}: {len(out):,} terms, {(out['FDR q'] <= 0.05).sum()} at FDR <= 0.05, "
                      f"{(out['Shown in'] != '').sum()} shown in figures/tables")
    readme = readme_df([
        ("Method", "GSEApy 1.3.0 prerank. Genes ranked by sign(log2FC) x -log10(nominal p) from each dataset's "
                   "DESeq2 result; 1,000 gene-set permutations; gene-set size 10 to 500; seed 42."),
        ("Libraries", "Enrichr KEGG_2019_Mouse, GO_Biological_Process_2026 and GO_Cellular_Component_2026. "
                      "GO CC was run for OSD-564 and OSD-562 only."),
        ("DE input", "OSD-564: GeneLab DESeq2 (GLDS-569, 5 vs 5). OSD-562: re-derived DESeq2 (4 vs 3; Supplementary "
                     "Table S1). OSD-352 bulk: GeneLab DESeq2 (GLDS-352, 3 vs 3; Supplementary Table S2)."),
        ("Sign", "NES > 0, enriched in spaceflight; NES < 0, enriched in ground control."),
        ("Significance", "FDR q <= 0.05. FDR q = 0 means q < 0.001 after 1,000 permutations."),
        ("Genes in set", "Members of the term present in the ranked list (denominator of Tag %)."),
        ("Shown in", "Main or supplementary figure panel and main table in which the term is displayed. Curated "
                     "panels were chosen before analysis (Supplementary Table S6)."),
        ("Sheets", "\n".join(counts)),
    ])
    write_book(os.path.join(OUT, "Supplementary_Table_S3_GSEA.xlsx"), [("README", readme)] + sheets)
    return sheets


# ══════════════════════════════════════════════════════════════════
# S4  g:Profiler over-representation analysis, OSD-564
# ══════════════════════════════════════════════════════════════════
GPROF = os.path.join(P564, "gprofiler", "Supplemental DEG Tables (OSD-562, OSD-564, OSD-352) - OSD-564 g_Profiler {}.csv")
GPROF_SOURCES = ["GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC", "WP", "TF", "MIRNA", "HPA", "CORUM", "HP"]
# Terms drawn in the ORA figure panels (make_gprofiler_ORA_figures_OSD564.py); all of them are in Table 4.
ORA_SHOWN = {
    "down": {**{t: "Fig. 2b" for t in ["GO:0005740", "GO:0005743", "GO:0006119", "GO:0033108", "GO:0009060",
                                       "GO:0042775", "GO:0033617", "GO:0042776", "GO:0032981"]},
             **{t: "Fig. 4b" for t in ["GO:0022626", "KEGG:03010", "GO:0140241"]}},
    "up": {**{t: "Fig. 3b" for t in ["GO:0016192", "GO:0045055", "GO:0016050", "GO:0043230", "GO:0051648",
                                     "GO:0072659"]},
           **{t: "Fig. 4b" for t in ["GO:0007267", "GO:0046907", "GO:0031267", "GO:0120036", "GO:0030036"]}},
}


def read_gprofiler(which):
    raw = pd.read_csv(GPROF.format(which), header=None, dtype=str)
    hdr = next(i for i in range(10) if "term_id" in raw.iloc[i].tolist())
    return pd.read_csv(GPROF.format(which), header=hdr)


def table_s4():
    de = pd.read_csv(os.path.join(P564, "GLDS-569_rna_seq_differential_expression_rRNArm_GLbulkRNAseq.csv"))
    fc = "Log2fc_(Ground Control)v(Space Flight)"            # GeneLab contrast is GC/FLT: flip the sign
    p, q = fc.replace("Log2fc_", "P.value_"), fc.replace("Log2fc_", "Adj.p.value_")
    de["log2FC"] = -de[fc]
    sel = de[de[q] < P_THR].copy()
    sel["Direction"] = np.where(sel["log2FC"] > 0, "Up in spaceflight", "Down in spaceflight")
    genes = (sel.rename(columns={"ENSEMBL": "Ensembl ID", "SYMBOL": "Gene", p: "p", q: "Adjusted p",
                                 "log2FC": "log2FC (FLT/GC)"})
             [["Direction", "Ensembl ID", "Gene", "log2FC (FLT/GC)", "p", "Adjusted p"]]
             .sort_values(["Direction", "Adjusted p"], ascending=[False, True]))
    n_up, n_dn = (genes.Direction == "Up in spaceflight").sum(), (genes.Direction == "Down in spaceflight").sum()
    # g:Profiler reports human-style upper-case symbols; restore mouse capitalisation from the DE table
    case = dict(zip(de["SYMBOL"].dropna().str.upper(), de["SYMBOL"].dropna()))
    sheets, counts = [], []
    for which, key, label in [("Upregulated Pathways", "up", "Up in spaceflight"),
                              ("Downregulated Pathways", "down", "Down in spaceflight")]:
        g = read_gprofiler(which)
        padj = g["adjusted_p_value"]
        assert np.allclose(10 ** -g["negative_log10_of_adjusted_p_value"], padj, rtol=5e-3), which
        assert (padj < 0.05).all()
        shown = g["term_id"].map(ORA_SHOWN[key])
        out = pd.DataFrame({
            "Source": g["source"], "Term ID": g["term_id"], "Term name": g["term_name"],
            "Adjusted p (g:SCS)": padj, "−log10 adjusted p": g["negative_log10_of_adjusted_p_value"].round(2),
            "Term size": g["term_size"], "Query size": g["query_size"], "Intersection size": g["intersection_size"],
            "Effective domain size": g["effective_domain_size"],
            "Driver term": np.where(g["highlighted"].astype(str).str.lower() == "true", "Yes", ""),
            "Shown in": np.where(shown.notna(), shown.fillna("") + "; Table 4", ""),
            "Intersection genes": g["intersections"].fillna("").map(
                lambda s: ", ".join(case.get(x, x) for x in s.split(",") if x)),
        })
        out["_src"] = out["Source"].map({s: i for i, s in enumerate(GPROF_SOURCES)}).fillna(99)
        out = out.sort_values(["_src", "Adjusted p (g:SCS)"]).drop(columns="_src")
        assert (out["Shown in"] != "").sum() == len(ORA_SHOWN[key]), "a figure term is missing from the export"
        sheets.append((label, out))
        counts.append(f"{label}: {len(out)} terms (" +
                      ", ".join(f"{s} {n}" for s, n in out["Source"].value_counts().reindex(GPROF_SOURCES).dropna().astype(int).items()) + ")")
    readme = readme_df([
        ("Contents", "Complete g:Profiler g:GOSt over-representation results for OSD-564 hippocampus (spaceflight vs "
                     "ground control): every term with g:SCS-adjusted p < 0.05, run separately for genes up- and "
                     "downregulated in spaceflight. Table 4 and Figures 2b, 3b and 4b show curated subsets."),
        ("Input genes", f"GeneLab DESeq2 results for GLDS-569 (rRNA-removed; n = 5 FLT vs 5 GC): genes with "
                        f"adjusted p < 0.05, any fold change, split by direction in spaceflight ({n_up} up, {n_dn} down; "
                        f"sheet 'Input genes'). Genes without a g:Profiler annotation in a source do not count "
                        f"towards that source's query size."),
        ("Settings", "Organism Mus musculus; statistical domain scope, annotated genes only; significance threshold, "
                     "g:SCS 0.05; data sources as listed in the Methods. Sources with significant terms: "
                     "GO:BP, GO:CC and GO:MF (Gene Ontology), KEGG, REAC (Reactome), WP (WikiPathways), "
                     "TF (TRANSFAC motifs), CORUM (protein complexes) and HP (Human Phenotype Ontology)."),
        ("Adjusted p (g:SCS)", "g:Profiler's multiple-testing-corrected p, as exported. −log10 adjusted p is the "
                               "value plotted in Figures 2b, 3b and 4b and reported in Table 4."),
        ("Sizes", "Term size, genes annotated to the term; Query size, input genes annotated in that source; "
                  "Intersection size, input genes in the term; Effective domain size, annotated genes in the "
                  "source (the statistical background)."),
        ("Driver term", "Yes, term highlighted by g:Profiler's driver-term filtering (GO only)."),
        ("Shown in", "Figure panel and main table in which the term is displayed."),
        ("Intersection genes", "Input genes annotated to the term. Mouse gene-symbol capitalisation restored from "
                               "the DE table (g:Profiler returns upper-case symbols)."),
        ("Sheets", "\n".join(counts)),
    ])
    write_book(os.path.join(OUT, "Supplementary_Table_S4_OSD-564_gProfiler_ORA.xlsx"),
               [("README", readme)] + sheets + [("Input genes", genes)])


# ══════════════════════════════════════════════════════════════════
# S5  OSD-352 spatial DE (Masarapu et al. 2024, Supplementary Data 8) + per-animal CA3
# ══════════════════════════════════════════════════════════════════
REGION = {6: "CA1 pyramidal layer", 8: "CA3", 10: "CA1 stratum radiatum", 11: "Dentate gyrus"}
SHOWN_SPATIAL = {"Pink1": "Fig. 2d; Table 5; Supp. Figs 6, 9", "Snph": "Fig. 2d; Table 5; Supp. Figs 6, 9",
                 "Itpr1": "Fig. 2d; Table 5; Supp. Figs 6, 9", "Rhob": "Fig. 3c; Table 5; Supp. Figs 7, 10",
                 "Sec14l1": "Fig. 3c; Table 5; Supp. Figs 7, 10", "Ubl3": "Fig. 3c; Table 5; Supp. Figs 7, 10",
                 "Arpc2": "Fig. 4c; Table 5; Supp. Figs 8, 11", "Arap2": "Fig. 4c; Table 5; Supp. Figs 8, 11",
                 "Rapgef5": "Results text"}


def table_s5():
    m = pd.read_excel(MASARAPU_SD8, header=4).dropna(subset=["primerid"])
    m["cluster"] = m["cluster"].astype(int)
    h = m[m["cluster"].isin(REGION) & (m["fdr"] < 0.05)].copy()
    deg = pd.DataFrame({
        "Gene": h["primerid"], "Cluster": h["cluster"], "Region": h["cluster"].map(REGION),
        "Coefficient (FLT vs GC)": h["coef"], "p": h["p_val"], "FDR": h["fdr"],
        "Category": [CATEGORY.get(g, "") for g in h["primerid"]],
        "Shown in": [SHOWN_SPATIAL.get(g, "") for g in h["primerid"]],
    }).sort_values(["Cluster", "FDR"]).reset_index(drop=True)
    pa = pd.read_csv(os.path.join(PST, "CA3_cluster8_per_animal_adjusted.csv"))
    pa = pa.rename(columns={"animal_id": "Animal", "label": "Figure label", "condition": "Condition",
                            "gene": "Gene", "n_sections": "Sections", "n_ca3_spots": "CA3 spots",
                            "raw_mean": "Unadjusted mean", "adj_cdr_mean": "CDR-adjusted mean",
                            "adj_cdr_depth_mean": "CDR + log-depth adjusted mean",
                            "frac_expressing": "Fraction of spots expressing", "mean_cdr": "Mean CDR"})
    pa = pa[["Gene", "Animal", "Figure label", "Condition", "Sections", "CA3 spots", "Unadjusted mean",
             "CDR-adjusted mean", "CDR + log-depth adjusted mean", "Fraction of spots expressing", "Mean CDR"]]
    summ = (pa.groupby(["Gene", "Condition"])["CDR-adjusted mean"].agg(["mean", "std"]).unstack())
    summ.columns = [f"{c[1].capitalize()} {c[0]}" for c in summ.columns]
    summ = summ.reset_index()
    readme = readme_df([
        ("Source", "Spatial DEG statistics are reproduced from Masarapu, Y. et al. Nat. Commun. 15, 4778 (2024) "
                   "(ref. 21), Supplementary Data 8, under a CC BY 4.0 licence. They were not recomputed."),
        ("Statistics", "Coefficient, p and FDR as reported in Supplementary Data 8 for spaceflight versus ground "
                       "control within each spatial cluster (3 FLT and 3 GC animals, two coronal sections each; "
                       "spot-level model as described in ref. 21). Note: the caption of Supplementary Data 8 "
                       "describes a two-sided Wilcoxon rank-sum test, while the coefficient columns and ref. 21's "
                       "Methods describe the MAST model used in this manuscript; check before submission."),
        ("Clusters", "The four hippocampal clusters were examined: 6 (CA1 pyramidal layer), 8 (CA3), "
                     "10 (CA1 stratum radiatum), 11 (dentate gyrus). Clusters 10 and 11 contain no gene at "
                     f"FDR < 0.05. Cluster 6: {(deg.Cluster == 6).sum()} genes; cluster 8: {(deg.Cluster == 8).sum()} genes."),
        ("Per-animal CA3", "Descriptive only. Expression = SCT log-normalized (Seurat SCT 'data' slot) in all cluster-8 "
                           "spots of both sections per animal. Cellular detection rate (CDR) = detected genes per "
                           "spot / 14,630, z-scored across the 1,512 cluster-8 spots; adjusted expression = "
                           "residual + intercept of OLS expression ~ 1 + CDR_z, averaged per animal. Covers the eight "
                           "genes in Table 5. Script: OSD-352_spatial/CA3_panels/plot_CA3_panels.py."),
    ])
    write_book(os.path.join(OUT, "Supplementary_Table_S5_OSD-352_spatial.xlsx"),
               [("README", readme), ("Spatial DEGs", deg), ("Per-animal CA3", pa), ("CA3 group summary", summ)])
    return deg


# ══════════════════════════════════════════════════════════════════
# S6  Curated gene and pathway lists (DRAFT: rationale and references to be completed by the authors)
# ══════════════════════════════════════════════════════════════════
OSD564_CURATED = ["Adam8", "Apoc1", "Aqp1", "Baiap2l2", "Bcat1", "C4a", "Casr", "Ect2", "Epcam", "Folr1", "Gas2",
                  "Gas2l3", "Ghsr", "Itga2", "Lpl", "Myh7b", "Rac2", "Rhof", "S100a4", "S100a8", "Smtnl2", "Timp1",
                  "Tmsb4x", "Tnni2", "Wdr1", "Wnt7a", "Zfas1", "Mrps25", "Chchd1", "Rmrp", "Ndufb4", "Ndufaf6",
                  "Atp5mk", "Pet100", "mt-Co2", "Arpc2", "Arap2", "Vps4a", "Ist1", "Prkaa1"]
TODO = "TO BE COMPLETED BY AUTHORS"


def table_s6(s1, s2, deg, gsea_sheets):
    t2 = set(OSD564_CURATED) - {"Arpc2", "Arap2"}
    genes = sorted(set(OSD564_CURATED) | set(s1["Gene"]) | set(s2["Gene"]) | set(SHOWN_SPATIAL))
    rows = []
    for g in genes:
        shown = []
        if g in OSD564_CURATED:
            shown.append("Fig. 1a (OSD-564)")
        if g in t2:
            shown.append("Table 2")
        if g in set(s1["Gene"]):
            shown.append("Supp. Fig. 2a (OSD-562)")
        if g in set(s2["Gene"]):
            shown.append("Supp. Fig. 2c (OSD-352 bulk)")
        if g in SHOWN_SPATIAL:
            shown.append(SHOWN_SPATIAL[g] + " (OSD-352 spatial)")
        rows.append((g, cat(g), "; ".join(shown), "", "", ""))
    gdf = pd.DataFrame(rows, columns=["Gene", "Category", "Shown in", "Function (one line)",
                                      "Rationale for category", "Reference(s)"])
    prows = []
    for name, df in gsea_sheets:
        ds, lib = name.split(" ", 1)
        for r in df[df["Shown in"] != ""].itertuples():
            term = r.Term
            if lib == "GO BP" and getattr(r, "_2", None) in GO_MITO or (lib == "KEGG" and norm(term) == "oxidative phosphorylation"):
                c = "Mitochondria"
            elif lib == "KEGG" and norm(term) in map(norm, KEGG_SEV):
                c = "sEV"
            elif lib == "KEGG" and norm(term) in map(norm, KEGG_TNT):
                c = "TNT"
            elif lib == "KEGG" and norm(term) in map(norm, KEGG_MET):
                c = "Metabolic adaptation"
            elif lib == "KEGG" and norm(term) == "ribosome":
                c = "Ribosome/translation"
            elif lib == "KEGG" and norm(term) == "parkinson disease":
                c = "Mitochondria"
            elif norm(term) in map(norm, ["Actin Filament Organization", "Axonogenesis", "Axon Guidance",
                                          "Neuron Projection Morphogenesis"]):
                c = "TNT"
            elif norm(term) in map(norm, ["Regulated Exocytosis", "Synaptic Vesicle Recycling",
                                          "Recycling Endosome Membrane", "Endocytic Vesicle",
                                          "Clathrin-Coated Vesicle", "Exocytic Vesicle"]):
                c = "sEV"
            elif norm(term) == "mitochondrial inner membrane":
                c = "Mitochondria"
            elif norm(term) == "ribosome":
                c = "Ribosome/translation"
            else:
                c = ""
            prows.append((term, getattr(r, "_2", "") if lib != "KEGG" else "", lib, c))
    pdf = (pd.DataFrame(prows, columns=["Term", "Term ID", "Library", "Category"])
           .drop_duplicates(["Term", "Library"]).sort_values(["Library", "Category", "Term"]))
    pdf["Rationale"] = ""
    pdf["Reference(s)"] = ""
    gdf.loc[:, ["Function (one line)", "Rationale for category", "Reference(s)"]] = TODO
    pdf.loc[:, ["Rationale", "Reference(s)"]] = TODO
    readme = readme_df([
        ("Status", "DRAFT. Gene lists, categories (as stated in the manuscript) and figure/table placement are "
                   "filled in from the analysis scripts. Function, rationale and references must be written by the "
                   "authors before submission."),
        ("Unassigned genes", "Genes labelled on a volcano plot that the manuscript does not place in the mitochondria, "
                             "sEV or TNT category. Assign a category with a rationale, or describe them as "
                             "dataset-specific genes of interest."),
        ("Panel sizes", f"The volcano label sets differ by dataset: OSD-564 {len(OSD564_CURATED)} genes "
                        f"(37 curated + Vps4a, Ist1, Prkaa1), OSD-562 {len(s1)}, OSD-352 bulk {len(s2)}. The text "
                        "should explain why."),
        ("References", "References cited only in this table go in the Supplementary Information's own reference "
                       "list, not the main reference list."),
    ])
    write_book(os.path.join(OUT, "Supplementary_Table_S6_curated_lists_DRAFT.xlsx"),
               [("README", readme), ("Genes", gdf), ("Pathways and terms", pdf)])


# ══════════════════════════════════════════════════════════════════
# S7  Per-sample RNA and library QC
# ══════════════════════════════════════════════════════════════════
def table_s7():
    qc = []
    for ds, f, pat in [("OSD-564", f"{P564}/GLDS-569_rna_seq_qc_metrics_GLbulkRNAseq.csv", r"RR10_HPC_"),
                       ("OSD-562", f"{P562}/GLDS-557_rna_seq_qc_metrics_GLbulkRNAseq.csv", r"_ISS-T_YNG_"),
                       ("OSD-352 bulk", f"{P352}/GLDS-352_rna_seq_qc_metrics_GLbulkRNAseq.csv", r"RR3_BRN_")]:
        q = pd.read_csv(f)
        q = q[q["sample"].str.contains(pat)].copy()
        q["Dataset"] = ds
        qc.append(q)
    q = pd.concat(qc, ignore_index=True)
    out = pd.DataFrame({
        "Dataset": q["Dataset"], "Sample": q["sample"], "Animal": q["sample"].str.split("_").str[-1],
        "Group": np.where(q["sample"].str.contains("_FLT_"), "Spaceflight", "Ground control"),
        "RIN": q["rin"], "Read depth": q["read_depth"], "rRNA contamination (%)": q["rrna_contamination"],
        "Uniquely mapped (%)": q["uniquely_mapped_percent"], "Multi-mapped (%)": q["multimapped_percent"],
        "Genes detected (> 10 counts)": q["gene_detected_gt10"],
    })
    # 7SL/7SK fraction for OSD-562 (same computation as rerun_deseq2_OSD562.py)
    de = pd.read_csv(os.path.join(P562, "OSD562_DESeq2_ISS-T_YNG_FLT4_vs_GC3.csv"), usecols=["ENSEMBL", "SYMBOL"])
    raw = pd.read_csv(os.path.join(P562, "OSD562_raw_counts_ISS-T_YNG.csv"), index_col=0)
    ids = de.loc[de["SYMBOL"].isin(["Rn7s1", "Rn7s2", "Rn7sk"]), "ENSEMBL"]
    frac = raw.loc[raw.index.intersection(ids)].sum() / raw.sum() * 100
    out["7SL/7SK RNA (% of counts)"] = out["Sample"].map(frac).round(1)
    readme_rows = [
        ("Source", "GeneLab per-sample QC metrics tables (GLDS-569, GLDS-557, GLDS-352 "
                   "*_rna_seq_qc_metrics_GLbulkRNAseq.csv; NASA OSDR)."),
        ("Samples", "Only the animals analysed in the manuscript: OSD-564 all 10; OSD-562 the ISS-terminal young "
                    "cohort (FY3, FY4, FY6, FY10; GY5, GY8, GY10); OSD-352 bulk all 6."),
        ("7SL/7SK RNA", "OSD-562 only: Rn7s1 + Rn7s2 + Rn7sk counts as a percentage of all counts in the GeneLab "
                        "RSEM rRNA-removed count table used for the re-derived DESeq2 contrast."),
    ]
    for ds, g in out.groupby("Dataset", sort=False):
        readme_rows.append((f"RIN, {ds}", f"mean {g.RIN.mean():.2f} ± {g.RIN.std():.2f} (SD), range "
                                          f"{g.RIN.min():.1f}-{g.RIN.max():.1f}"))
    write_book(os.path.join(OUT, "Supplementary_Table_S7_sample_QC.xlsx"),
               [("README", readme_df(readme_rows)), ("Per-sample QC", out)])
    return out


if __name__ == "__main__":
    main_tables()
    s1 = table_s1()
    s2 = table_s2()
    gsea = table_s3()
    table_s4()
    deg = table_s5()
    table_s6(s1, s2, deg, gsea)
    s7 = table_s7()
