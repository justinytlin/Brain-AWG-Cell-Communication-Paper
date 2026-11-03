#!/usr/bin/env python3
"""
Whole-section spatial expression maps in the layout of the original
"<gene>_spatial_combined" figures (gene columns; Spaceflight | Ground control;
three sections each; one Low/High colour bar), but scientifically accurate:

  * colour = SCT log-normalised expression (log1p of SCT-corrected counts),
    rescaled WITHIN EACH SECTION from 0 to that section's 99th percentile, so the
    3- to 5-fold difference in sequencing depth between sections does not make
    deep sections look uniformly "high" for every gene;
  * hippocampal CA3 spots (ST cluster 8, the only region the Table statistics
    refer to) are drawn with a thin black outline;
  * each section carries its animal label (F1..F3, G1..G3; original ID in
    parentheses);
  * one section per animal, animals in a fixed order in every figure.

Inputs: spots.csv and expr.csv written by extract_spatial_data.R (same folder).
Outputs: Fig2c_*_maps.png/.pdf, Fig3c_*_maps.png/.pdf, Fig4c_*_maps.png/.pdf
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
PCT_HIGH = 99
CA3_CLUSTER = 8
CMAP = LinearSegmentedColormap.from_list("greyred", ["#d9d9d9", "#f4a582", "#d6604d", "#b2182b", "#8b0000"])
# (section, label, original animal id); one section per animal, fixed order
FLIGHT = [("Sample_159_C1", "F1", "CF1"), ("Sample_158_A1", "F2", "CF2"), ("Sample_158_C1", "F3", "CF7")]
GROUND = [("Sample_159_B1", "G1", "CG7"), ("Sample_304_C1", "G2", "CG8"), ("Sample_304_A1", "G3", "CG9")]
PANELS = {"Fig2c": ["Pink1", "Snph", "Itpr1"],
          "Fig3c": ["Rhob", "Sec14l1", "Ubl3"],
          "Fig4c": ["Arpc2", "Arap2"]}

spots = pd.read_csv(os.path.join(HERE, "spots.csv"))
expr = pd.read_csv(os.path.join(HERE, "expr.csv"))
df = spots.merge(expr, on="barcode")
sub = {s: df[df["section.name"] == s].copy() for s, _, _ in FLIGHT + GROUND}
# common data range so every section is drawn at the same physical scale
rng = max(max(d.pixel_x.max() - d.pixel_x.min(), d.pixel_y.max() - d.pixel_y.min()) for d in sub.values()) * 1.04
nn_px = np.median([np.median(cKDTree(d[["pixel_x", "pixel_y"]].values)
                             .query(d[["pixel_x", "pixel_y"]].values, k=2)[0][:, 1]) for d in sub.values()])


def marker_area(ax, fig, scale=0.92):
    """scatter 's' so that spots just touch (diameter = nearest-neighbour spacing)."""
    bbox = ax.get_window_extent(renderer=fig.canvas.get_renderer())
    pts_per_unit = bbox.width / rng * 72.0 / fig.dpi
    return (nn_px * pts_per_unit * scale) ** 2


def draw_section(ax, fig, d, gene, label):
    v = d[f"{gene}_data"].values
    vmax = max(np.percentile(v, PCT_HIGH), 1e-6)
    norm = Normalize(0, vmax, clip=True)
    ca3 = (d["seurat_clusters"] == CA3_CLUSTER).values
    xy = d[["pixel_x", "pixel_y"]].values
    # CA3 boundary = CA3 spots with at least one non-CA3 neighbour (or a missing neighbour)
    tree = cKDTree(xy)
    nbrs = tree.query_ball_point(xy[ca3], r=1.25 * nn_px)
    boundary_idx = np.where(ca3)[0][[any((not ca3[j]) for j in nb) or len(nb) < 6 for nb in nbrs]]
    bnd = np.zeros(len(d), bool); bnd[boundary_idx] = True
    s = marker_area(ax, fig)
    ax.scatter(xy[~bnd, 0], xy[~bnd, 1], c=v[~bnd], cmap=CMAP, norm=norm, s=s, linewidths=0, rasterized=True)
    ax.scatter(xy[bnd, 0], xy[bnd, 1], c=v[bnd], cmap=CMAP, norm=norm, s=s,
               edgecolors="black", linewidths=0.6, rasterized=True)
    cx, cy = 0.5 * (d.pixel_x.min() + d.pixel_x.max()), 0.5 * (d.pixel_y.min() + d.pixel_y.max())
    ax.set_xlim(cx - rng / 2, cx + rng / 2)
    ax.set_ylim(cy + rng / 2, cy - rng / 2)          # image orientation (y grows downward)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(label, fontsize=7, color="#333333", pad=2, loc="left")


def build(fig_id, genes):
    ng = len(genes)
    col_w, row_h = 1.55, 1.55                       # inches per section panel
    fig_w = ng * (2 * col_w + 0.35) + 0.2
    fig_h = 3 * row_h + 1.75
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=300)
    gs = fig.add_gridspec(nrows=3, ncols=ng * 3 - 1,
                          width_ratios=sum([[1, 1, 0.22]] * ng, [])[:-1],
                          left=0.01, right=0.99, top=0.86, bottom=0.19, wspace=0.05, hspace=0.16)
    fig.canvas.draw()
    for gi, gene in enumerate(genes):
        c0 = gi * 3
        for r in range(3):
            sF, labF, origF = FLIGHT[r]
            sG, labG, origG = GROUND[r]
            axF = fig.add_subplot(gs[r, c0]); draw_section(axF, fig, sub[sF], gene, f"{labF} ({origF})")
            axG = fig.add_subplot(gs[r, c0 + 1]); draw_section(axG, fig, sub[sG], gene, f"{labG} ({origG})")
            if r == 0:
                axF.text(0.5, 1.12, "Spaceflight", transform=axF.transAxes, ha="center", va="bottom", fontsize=9)
                axG.text(0.5, 1.12, "Ground control", transform=axG.transAxes, ha="center", va="bottom", fontsize=9)
        # gene title centred over its two sub-columns
        xl = fig.add_subplot(gs[0, c0]).get_position(); xr = fig.add_subplot(gs[0, c0 + 1]).get_position()
        for a in fig.axes[-2:]:
            a.remove()
        fig.text(0.5 * (xl.x0 + xr.x1), 0.955, gene, ha="center", va="center", fontsize=13, fontweight="bold", fontstyle="italic")
        if gi < ng - 1:
            xsep = xr.x1 + 0.5 * (0.22 / (2 + 0.22)) * (xr.x1 - xl.x0) / 2 * 1.0
            fig.add_artist(plt.Line2D([xr.x1 + 0.012, xr.x1 + 0.012], [0.19, 0.97], transform=fig.transFigure,
                                      color="#9a9a9a", linewidth=0.8))
    # colour bar
    cax = fig.add_axes([0.5 - 0.09, 0.075, 0.18, 0.018])
    cb = fig.colorbar(ScalarMappable(norm=Normalize(0, 1), cmap=CMAP), cax=cax, orientation="horizontal")
    cb.set_ticks([0, 1]); cb.set_ticklabels(["Low", "High"]); cb.ax.tick_params(labelsize=8, length=0)
    cb.outline.set_linewidth(0.4)
    fig.text(0.5, 0.108, "Expression (scaled within each section)", ha="center", va="bottom", fontsize=8.5)
    fig.text(0.5, 0.012,
             f"SCT log-normalised expression per Visium spot, rescaled within each section from 0 to its {PCT_HIGH}th percentile; "
             "one coronal section per animal. Black outline = hippocampal CA3 (ST cluster 8), the region tested in Table 5. "
             "Data: Masarapu et al. 2024 (OSD-352, RR-3).",
             ha="center", va="bottom", fontsize=6.2, color="#444444", wrap=True)
    base = os.path.join(HERE, f"{fig_id}_{'_'.join(genes)}_maps")
    fig.savefig(base + ".png", dpi=300); fig.savefig(base + ".pdf", dpi=300)
    plt.close(fig)
    print("wrote", base + ".png")


if __name__ == "__main__":
    for fid, genes in PANELS.items():
        build(fid, genes)
