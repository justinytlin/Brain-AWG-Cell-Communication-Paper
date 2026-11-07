#!/usr/bin/env python3
"""
Assemble the final Scientific Reports figures: each figure is ONE image file
(PNG 300 dpi, LZW TIFF, PDF) built from the per-panel PNGs that the folder
scripts export at fixed physical sizes (see srep_style.PANEL_MM).  Panels are
placed 1:1 (no rescaling) on a white canvas 183 mm wide (double column) and
each panel gets a bold lower-case letter in 8 pt Arial, the same size as all
other lettering in the figures.

Run from the project root:  /opt/anaconda3/bin/python3 assemble_figures.py
Outputs go to Figures_SciRep/.
"""
import os, sys
from PIL import Image, ImageDraw, ImageFont

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "Figures_SciRep")
DPI = 300
PX = DPI / 25.4                      # pixels per mm
FONT_PT = 8
FONT_PX = int(round(FONT_PT / 72 * DPI))
BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
WIDTH_MM = 183

P564 = os.path.join(ROOT, "OSD-564 Figures")
P562 = os.path.join(ROOT, "OSD-562 Figures")
P352 = os.path.join(ROOT, "OSD-352", "NEW FIGURES")
PST = os.path.join(ROOT, "OSD-352", "NASA_RR3_Brain-main", "manuscript_plots", "CA3_panels")

# (letter, file, x_mm, y_mm, expected_w_mm, expected_h_mm)
FIGURES = {
    "Figure_1": [
        ("a", f"{P564}/OSD564_Volcano_plot_annotated.png", 0, 0, 118, 100),
        ("b", f"{P564}/OSD564_PCA_plot.png", 123, 0, 60, 60),
    ],
    "Figure_2": [
        ("a", f"{P564}/GO_results/GSEA_mito_pathways_NES.png", 0, 0, 90, 60),
        ("b", f"{P564}/gprofiler/gprofiler_ORA_mito_downregulated_pathways.png", 93, 0, 90, 65),
        ("c", f"{P562}/GO_BP_results_OSD562/GSEA_Mito_NES_summary.png", 0, 66, 90, 60),
        ("d", f"{PST}/Fig2c_Pink1_Snph_Itpr1_maps.png", 0, 130, 183, 112),
    ],
    "Figure_3": [
        ("a", f"{P564}/KEGG_results/pathway_figures/GSEA_sEV_pathways_NES.png", 0, 0, 90, 60),
        ("b", f"{P564}/gprofiler/gprofiler_ORA_sEVs_pathways.png", 93, 0, 90, 65),
        ("c", f"{PST}/Fig3c_Rhob_Sec14l1_Ubl3_maps.png", 0, 69, 183, 112),
    ],
    "Figure_4": [
        ("a", f"{P564}/KEGG_results/pathway_figures/GSEA_TNT_pathways_NES.png", 0, 0, 90, 60),
        ("b", f"{P564}/gprofiler/gprofiler_ORA_TNTs_pathways.png", 93, 0, 90, 65),
        ("c", f"{PST}/Fig4c_Arpc2_Arap2_maps.png", 0, 74, 122, 105),
    ],
    "Supplementary_Figure_1": [          # read-quality metrics, one row of panels per dataset
        ("a", f"{P564}/QC_figures/OSD564_QC_summary.png", 0, 0, 183, 88),
        ("b", f"{P562}/QC_figures/OSD562_QC_summary.png", 0, 92, 183, 70),
        ("c", f"{P352}/QC_figures/OSD352_read_dist_QC_summary.png", 0, 166, 183, 70),
    ],
    "Supplementary_Figure_2": [          # same layout as Figure 1, one row per dataset
        ("a", f"{P562}/OSD562_Volcano_custom_genes.png", 0, 0, 118, 100),
        ("b", f"{P562}/OSD562_PCA_plot.png", 123, 0, 60, 60),
        ("c", f"{P352}/OSD352_Volcano_plot_annotated.png", 0, 104, 118, 100),
        ("d", f"{P352}/OSD352_PCA_plot.png", 123, 104, 60, 60),
    ],
    # Supplementary Figures are numbered in order of first citation in the manuscript (v3):
    # 3 = OSD-352 bulk GSEA, 4 = extra OSD-564 GSEA, 5 = OSD-562 KEGG GSEA,
    # 6-8 = CDR-adjusted per-animal CA3 means, 9-11 = unadjusted per-animal CA3 means.
    "Supplementary_Figure_3": [
        ("a", f"{P352}/KEGG_results_OSD352/pathway_figures/GSEA_TNTs_Factin_NES_summary.png", 0, 0, 90, 60),
        ("b", f"{P352}/KEGG_results_OSD352/pathway_figures/GSEA_sEVs_Docking_NES_summary.png", 93, 0, 90, 60),
        ("c", f"{P352}/GO_BP_results_OSD352/GSEA_Mito_NES_summary.png", 0, 63, 90, 60),
    ],
    "Supplementary_Figure_4": [
        ("a", f"{P564}/KEGG_results/pathway_figures/GSEA_metabolic_adaptation_pathways_NES.png", 0, 0, 90, 60),
        ("b", f"{P564}/GO_results/GSEA_BP_protrusion_exocytosis_NES_barplot.png", 93, 0, 90, 60),
        ("c", f"{P564}/GO_results/GSEA_CC_compartments_NES_barplot.png", 0, 63, 90, 60),
    ],
    "Supplementary_Figure_5": [
        ("a", f"{P562}/KEGG_results_OSD562/pathway_figures/GSEA_TNT_NES_summary.png", 0, 0, 90, 60),
        ("b", f"{P562}/KEGG_results_OSD562/pathway_figures/GSEA_sEV_NES_summary.png", 93, 0, 90, 60),
        ("c", f"{P562}/KEGG_results_OSD562/pathway_figures/GSEA_Metabolic_NES_summary.png", 0, 63, 90, 60),
    ],
    "Supplementary_Figure_6": [("", f"{PST}/Fig2c_Pink1_Snph_Itpr1_CA3.png", 0, 0, 183, 235)],
    "Supplementary_Figure_7": [("", f"{PST}/Fig3c_Rhob_Sec14l1_Ubl3_CA3.png", 0, 0, 183, 235)],
    "Supplementary_Figure_8": [("", f"{PST}/Fig4c_Arpc2_Arap2_CA3.png", 0, 0, 122, 235)],
    "Supplementary_Figure_9": [("", f"{PST}/Fig2c_Pink1_Snph_Itpr1_CA3_raw_means.png", 0, 0, 183, 90)],
    "Supplementary_Figure_10": [("", f"{PST}/Fig3c_Rhob_Sec14l1_Ubl3_CA3_raw_means.png", 0, 0, 183, 90)],
    "Supplementary_Figure_11": [("", f"{PST}/Fig4c_Arpc2_Arap2_CA3_raw_means.png", 0, 0, 122, 90)],
}

LETTER_GUTTER_MM = 3.0   # each lettered panel is pushed down by this much; the letter sits in the gap above it


def mm2px(v):
    return int(round(v * PX))


def load_panel(path, w_mm, h_mm):
    im = Image.open(path).convert("RGB")
    exp = (mm2px(w_mm), mm2px(h_mm))
    tol = 0.03
    if abs(im.width / exp[0] - 1) > tol or abs(im.height / exp[1] - 1) > tol:
        # never rescale (that would change the 8 pt lettering); place at the panel's true size and warn
        print(f"  note: {os.path.basename(path)} is {im.width / PX:.1f} x {im.height / PX:.1f} mm "
              f"(layout assumed {w_mm} x {h_mm}); placed 1:1")
    return im


def assemble(name, panels):
    font = ImageFont.truetype(BOLD, FONT_PX)
    ims = []
    for l, f, x, y, w, h in panels:
        im = load_panel(f, w, h)
        xp, yp = mm2px(x), mm2px(y)
        # push the panel down if it would overlap a panel already placed above it
        for _, im0, x0, y0 in ims:
            if xp < x0 + im0.width and x0 < xp + im.width and yp < y0 + im0.height and y0 < yp:
                yp = y0 + im0.height + mm2px(1)
        if l:
            yp += mm2px(LETTER_GUTTER_MM)
        ims.append((l, im, xp, yp))
    W = mm2px(WIDTH_MM) if len(panels) > 1 or panels[0][4] >= WIDTH_MM else ims[0][1].width
    H = max(y + im.height for _, im, _, y in ims)
    canvas = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(canvas)
    for letter, im, x, y in ims:
        canvas.paste(im, (x, y))
        if letter:
            draw.text((x + mm2px(0.5), y - mm2px(LETTER_GUTTER_MM) + mm2px(0.3)), letter, fill="black", font=font)
    os.makedirs(OUT, exist_ok=True)
    base = os.path.join(OUT, name)
    canvas.save(base + ".png", dpi=(DPI, DPI))
    canvas.save(base + ".tif", compression="tiff_lzw", dpi=(DPI, DPI))
    canvas.save(base + ".pdf", resolution=DPI)
    print(f"{name}: {W / PX:.0f} x {H / PX:.0f} mm, {len(panels)} panel(s) -> {base}.png/.tif/.pdf")


if __name__ == "__main__":
    wanted = sys.argv[1:] or list(FIGURES)
    for n in wanted:
        missing = [f for _, f, *_ in FIGURES[n] if not os.path.exists(f)]
        if missing:
            print(f"{n}: SKIPPED, missing {missing}")
            continue
        assemble(n, FIGURES[n])
