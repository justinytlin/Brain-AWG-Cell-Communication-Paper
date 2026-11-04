#!/usr/bin/env python3
"""
Whole-section spatial expression maps in the layout of the original
"<gene>_spatial_combined" figures (gene columns; Spaceflight | Ground control;
three sections each; one Low/High colour bar), but scientifically accurate:

  * colour = SCT log-normalized expression (log1p of SCT-corrected counts),
    rescaled WITHIN EACH SECTION from 0 to that section's 99th percentile, so the
    3- to 5-fold difference in sequencing depth between sections does not make
    deep sections look uniformly "high" for every gene;
  * hippocampal CA3 (ST cluster 8, the only region the Table statistics refer to)
    is outlined with a 1 pt black boundary line (union of the cluster-8 spot
    footprints, see ca3_outline);
  * each section carries its animal label (F1..F3, G1..G3; original ID in
    parentheses);
  * one section per animal, animals in a fixed order in every figure;
  * every section is drawn at the same scale; a 1 mm scale bar (10 Visium spot
    pitches of 100 um) with its length written on it sits below the bottom-left
    section of every gene block.

Style: srep_style (Scientific Reports: Arial 8 pt everywhere, no bold, no line
thinner than 1 pt, sentence case, white background).  Panel sizes come from
srep.figsize_mm("maps3") = 183 x 112 mm (three genes) and "maps2" = 122 x 105 mm
(two genes); the assembly script places them 1:1 and adds the panel letters.

Inputs: spots.csv and expr.csv written by extract_spatial_data.R (same folder).
Outputs: Fig2c_*_maps.png/.pdf, Fig3c_*_maps.png/.pdf, Fig4c_*_maps.png/.pdf
Run:  /opt/anaconda3/bin/python3 plot_whole_section_maps.py
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srep_style as srep
srep.apply()
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree
from shapely.geometry import Point
from shapely.ops import unary_union

HERE = os.path.dirname(os.path.abspath(__file__))
PCT_HIGH = 99
CA3_CLUSTER = 8
SPOT_PITCH_UM = 100.0                 # Visium spot centre-to-centre distance
SCALE_BAR_MM = 1.0                    # length of the scale bar
CMAP = LinearSegmentedColormap.from_list("greyred", ["#d9d9d9", "#f4a582", "#d6604d", "#b2182b", "#8b0000"])
# (section, label, original animal id); one section per animal, fixed order
FLIGHT = [("Sample_159_C1", "F1", "CF1"), ("Sample_158_A1", "F2", "CF2"), ("Sample_158_C1", "F3", "CF7")]
GROUND = [("Sample_159_B1", "G1", "CG7"), ("Sample_304_C1", "G2", "CG8"), ("Sample_304_A1", "G3", "CG9")]
PANELS = {"Fig2c": ["Pink1", "Snph", "Itpr1"],
          "Fig3c": ["Rhob", "Sec14l1", "Ubl3"],
          "Fig4c": ["Arpc2", "Arap2"]}
MM = srep.MM
FS = srep.FONT_SIZE
LW = srep.LINE_PT
BAR_STRIP = 0.14                      # extra window height (fraction of rng) below row-3 maps for the scale bar

spots = pd.read_csv(os.path.join(HERE, "spots.csv"))
expr = pd.read_csv(os.path.join(HERE, "expr.csv"))
df = spots.merge(expr, on="barcode")
sub = {s: df[df["section.name"] == s].copy() for s, _, _ in FLIGHT + GROUND}
# common data range so every section is drawn at the same physical scale
rng = max(max(d.pixel_x.max() - d.pixel_x.min(), d.pixel_y.max() - d.pixel_y.min()) for d in sub.values()) * 1.04
# median nearest-neighbour spot spacing in image pixels = 100 um (Visium pitch)
nn_px = np.median([np.median(cKDTree(d[["pixel_x", "pixel_y"]].values)
                             .query(d[["pixel_x", "pixel_y"]].values, k=2)[0][:, 1]) for d in sub.values()])
PX_PER_MM = nn_px * 1000.0 / SPOT_PITCH_UM          # 1 mm = 10 spot pitches


def marker_area(ax, fig, scale=0.92):
    """scatter 's' so that spots just touch (diameter = nearest-neighbour spacing).
    fig.canvas.draw() / ax.apply_aspect() must have been called so the axes bbox is final."""
    bbox = ax.get_window_extent(renderer=fig.canvas.get_renderer())
    pts_per_unit = bbox.width / rng * 72.0 / fig.dpi
    return (nn_px * pts_per_unit * scale) ** 2


def ca3_outline(xy):
    """Boundary polygon(s) of the CA3 region: union of discs (radius 0.58 x spot pitch,
    so that adjacent spots merge) around the cluster-8 spots, lightly smoothed."""
    geom = unary_union([Point(x, y).buffer(0.58 * nn_px, 12) for x, y in xy])
    geom = geom.buffer(0.3 * nn_px).buffer(-0.3 * nn_px)
    return [geom] if geom.geom_type == "Polygon" else list(geom.geoms)


def draw_scale_bar(ax, x_right, y_bar):
    """1 mm bar (1.5 pt) with its length written directly above it; data coordinates."""
    x0 = x_right - SCALE_BAR_MM * PX_PER_MM
    ax.add_line(Line2D([x0, x_right], [y_bar, y_bar], color="black", linewidth=1.5,
                       solid_capstyle="butt", zorder=5))
    ax.text(0.5 * (x0 + x_right), y_bar - 0.035 * rng, f"{SCALE_BAR_MM:g} mm", ha="center", va="bottom",
            fontsize=FS, zorder=5)


def draw_section(ax, fig, d, gene, label, scale_bar=False, bottom_strip=0.0):
    v = d[f"{gene}_data"].values
    vmax = max(np.percentile(v, PCT_HIGH), 1e-6)
    norm = Normalize(0, vmax, clip=True)
    ca3 = (d["seurat_clusters"] == CA3_CLUSTER).values
    xy = d[["pixel_x", "pixel_y"]].values
    cx, cy = 0.5 * (d.pixel_x.min() + d.pixel_x.max()), 0.5 * (d.pixel_y.min() + d.pixel_y.max())
    ax.set_xlim(cx - rng / 2, cx + rng / 2)
    ax.set_ylim(cy + rng / 2 + bottom_strip * rng, cy - rng / 2)   # image orientation (y grows downward)
    ax.set_aspect("equal")
    ax.apply_aspect()                                            # final bbox before sizing markers
    s = marker_area(ax, fig)
    ax.scatter(xy[:, 0], xy[:, 1], c=v, cmap=CMAP, norm=norm, s=s, linewidths=0, rasterized=True, zorder=1)
    for poly in ca3_outline(xy[ca3]):
        for ring in [poly.exterior, *poly.interiors]:
            rx, ry = ring.xy
            ax.plot(rx, ry, color="black", linewidth=LW, solid_joinstyle="round", zorder=3)
    if scale_bar:
        draw_scale_bar(ax, cx + rng / 2 - 0.03 * rng, cy + rng / 2 + (bottom_strip - 0.03) * rng)
    ax.axis("off")
    ax.set_title(label, fontsize=FS, pad=2, loc="left")


def build(fig_id, genes):
    ng = len(genes)
    fig = plt.figure(figsize=srep.figsize_mm("maps3" if ng == 3 else "maps2"))
    fw_mm, fh_mm = fig.get_size_inches() / MM
    top_mm = 4.0 + 4.0 + 3.5           # gene title | Spaceflight/Ground control | animal label
    bot_mm = 12.0                      # colour bar + its label
    gs = fig.add_gridspec(nrows=3, ncols=ng * 3 - 1,
                          width_ratios=sum([[1, 1, 0.3]] * ng, [])[:-1],
                          height_ratios=[1, 1, 1 + BAR_STRIP],
                          left=1.0 / fw_mm, right=1 - 1.0 / fw_mm,
                          top=1 - top_mm / fh_mm, bottom=bot_mm / fh_mm, wspace=0.04, hspace=0.16)
    fig.canvas.draw()
    blocks = []
    for gi, gene in enumerate(genes):
        c0 = gi * 3
        axes = {}
        for r in range(3):
            sF, labF, origF = FLIGHT[r]
            sG, labG, origG = GROUND[r]
            strip = BAR_STRIP if r == 2 else 0.0
            axF = fig.add_subplot(gs[r, c0])
            draw_section(axF, fig, sub[sF], gene, f"{labF} ({origF})", scale_bar=(r == 2), bottom_strip=strip)
            axG = fig.add_subplot(gs[r, c0 + 1])
            draw_section(axG, fig, sub[sG], gene, f"{labG} ({origG})", bottom_strip=strip)
            axes[(r, 0)], axes[(r, 1)] = axF, axG
        pF, pG = axes[(0, 0)].get_position(), axes[(0, 1)].get_position()
        y_head = pF.y1 + 3.5 / fh_mm                              # above the animal label
        fig.text(0.5 * (pF.x0 + pF.x1), y_head, srep.sentence_case("Spaceflight"), ha="center", va="bottom", fontsize=FS)
        fig.text(0.5 * (pG.x0 + pG.x1), y_head, srep.sentence_case("Ground control"), ha="center", va="bottom", fontsize=FS)
        fig.text(0.5 * (pF.x0 + pG.x1), y_head + 4.0 / fh_mm, gene, ha="center", va="bottom", fontsize=FS,
                 fontstyle="italic")
        blocks.append((pF.x0, pG.x1, axes[(2, 0)].get_position().y0, y_head + 4.0 / fh_mm + 3.0 / fh_mm))
    # vertical separators between gene blocks (1 pt)
    for (x0a, x1a, y0, y1), (x0b, _, _, _) in zip(blocks[:-1], blocks[1:]):
        xs = 0.5 * (x1a + x0b)
        fig.add_artist(Line2D([xs, xs], [y0, y1], transform=fig.transFigure, color="#9a9a9a", linewidth=LW))
    # colour bar (one label line; ticks Low / High)
    cb_w, cb_h = 30.0 / fw_mm, 3.0 / fh_mm
    cax = fig.add_axes([0.5 - cb_w / 2, 3.5 / fh_mm, cb_w, cb_h])
    cb = fig.colorbar(ScalarMappable(norm=Normalize(0, 1), cmap=CMAP), cax=cax, orientation="horizontal")
    cb.set_ticks([0, 1]); cb.set_ticklabels(["Low", "High"]); cb.ax.tick_params(labelsize=FS, length=0, pad=1.5)
    cb.outline.set_linewidth(LW)
    fig.text(0.5, 3.5 / fh_mm + cb_h + 1.2 / fh_mm, srep.sentence_case("Expression (scaled within each section)"),
             ha="center", va="bottom", fontsize=FS)
    base = os.path.join(HERE, f"{fig_id}_{'_'.join(genes)}_maps")
    srep.save(fig, base)
    plt.close(fig)
    print("wrote", base + ".png / .pdf")


if __name__ == "__main__":
    print(f"spot pitch (median nearest-neighbour distance) = {nn_px:.2f} px = {SPOT_PITCH_UM:g} um; "
          f"1 mm = {PX_PER_MM:.1f} px; map window = {rng:.0f} px")
    for fid, genes in PANELS.items():
        build(fid, genes)
