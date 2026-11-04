#!/usr/bin/env python3
"""
plot_CA3_panels.py

Builds the spatial-transcriptomics figure panels (Fig2c, Fig3c, Fig4c) for the
hippocampal CA3 (ST cluster 8) story from the CSVs exported by
extract_spatial_data.R (spots.csv, expr.csv in the same folder).

For every gene, one column with three rows:
  row 1  whole-section maps  : all spots of the section, colour = SCT log1p
                               expression scaled PER SECTION (0 -> section's
                               99th percentile); the CA3 (cluster 8) region is
                               outlined with a 1 pt black boundary line.
  row 2  CA3-only maps       : only cluster-8 spots coloured (same per-section
                               scaling), every other spot in light grey, same
                               1 pt boundary line.
  row 3  quantitative panel  : per-animal mean of detection-rate (CDR)-adjusted
                               SCT log-expression over ALL cluster-8 spots of
                               both sections of the animal (3 flight animals,
                               red vs 3 ground animals, blue); mean +/- SD across
                               the three animals as an error bar; MAST coef / FDR
                               from Masarapu et al. annotated; y axis starts at 0.
                               CDR adjustment: expression ~ 1 + z(nFeature_SCT /
                               14630) fitted by OLS across all cluster-8 spots;
                               adjusted = residual + intercept.

Map layout inside each block: 3 rows (animals 1..3) x 2 columns (Spaceflight |
Ground control).  All maps share one data window, i.e. one scale; a 1 mm scale
bar (10 Visium spot pitches of 100 um) with its length written on it sits below
the bottom-left section of every block.

Supplementary Fig*_CA3_raw_means.png/.pdf: the unadjusted per-animal means
(figure section only, and both sections pooled), same error bars.

Style: srep_style (Scientific Reports: Arial 8 pt for all lettering, no bold,
no line thinner than 1 pt, sentence case, white background).  Panel sizes:
srep.figsize_mm("ca3_3") = 183 x 235 mm, "ca3_2" = 122 x 235 mm,
"raw_means3" = 183 x 90 mm, "raw_means2" = 122 x 90 mm.  Long explanatory
footnotes are not drawn (they belong in the figure legend); panel letters are
added by the assembly script.

Outputs (PNG 300 dpi + PDF) and CSVs with per-section / per-animal statistics.
Run:  /opt/anaconda3/bin/python3 plot_CA3_panels.py
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
from matplotlib import gridspec, colors
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from scipy.spatial import cKDTree
from shapely.geometry import Point
from shapely.ops import unary_union

HERE = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------------
# Design decisions / constants
# ----------------------------------------------------------------------------
CA3_CLUSTER = 8                      # seurat_clusters value for hippocampal CA3
PCT_HIGH = 99                        # per-section upper limit of colour scale
EXPR_SLOT = "data"                   # SCT "data" = log1p(SCT-corrected counts)
SPOT_PITCH_UM = 100.0                # Visium spot centre-to-centre distance
SCALE_BAR_MM = 1.0                   # scale-bar length

# section -> (animal label used in the figure, original animal id from sample_name_cond)
FLIGHT = [("Sample_159_C1", "F1", "CF1"), ("Sample_158_A1", "F2", "CF2"), ("Sample_158_C1", "F3", "CF7")]
GROUND = [("Sample_159_B1", "G1", "CG7"), ("Sample_304_C1", "G2", "CG8"), ("Sample_304_A1", "G3", "CG9")]

PANELS = {
    "Fig2c": ["Pink1", "Snph", "Itpr1"],
    "Fig3c": ["Rhob", "Sec14l1", "Ubl3"],
    "Fig4c": ["Arpc2", "Arap2"],
}
# Masarapu et al. MAST (spot level, sample = random effect) coefficient / FDR
MAST = {
    "Pink1": (-0.1039, 0.043), "Snph": (-0.1696, 0.049), "Itpr1": (+0.2455, 0.040),
    "Rhob": (+0.4878, 0.027), "Sec14l1": (+0.2976, 0.043), "Ubl3": (+0.3445, 0.048),
    "Arpc2": (+0.4608, 0.039), "Arap2": (+0.1898, 0.048),
}

COL_FLIGHT = srep.COLORS["flight"]
COL_GROUND = srep.COLORS["ground"]
COL_CONTEXT = "#dcdcdc"              # non-CA3 spots on the CA3-only maps
CMAP = colors.LinearSegmentedColormap.from_list(
    "greyred", ["#d9d9d9", "#f4a582", "#d6604d", "#b2182b", "#8b0000"])   # grey -> red, as in the whole-section maps

MM = srep.MM
FS = srep.FONT_SIZE                  # the ONE font size (8 pt) for every piece of lettering
LW = srep.LINE_PT                    # minimum line width (1 pt)
LABEL_STRIP = 0.27                   # extra window height (fraction of WIN) above each map for its two-line label
BAR_STRIP = 0.14                     # extra window height below row-3 maps for the scale bar


def sc(label):
    """sentence case, applied line by line (srep.sentence_case splits on spaces only)."""
    return "\n".join(srep.sentence_case(l) for l in label.split("\n"))


# ----------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------
spots = pd.read_csv(os.path.join(HERE, "spots.csv"))
expr = pd.read_csv(os.path.join(HERE, "expr.csv"))
assert (spots.barcode == expr.barcode).all()
df = pd.concat([spots, expr.drop(columns="barcode")], axis=1)
df["is_ca3"] = df.seurat_clusters == CA3_CLUSTER

SECTIONS = FLIGHT + GROUND
sec_names = [s for s, _, _ in SECTIONS]
sub = {s: df[df["section.name"] == s].copy() for s in sec_names}

# common spatial window so all maps share the same scale / spot size (pixels of
# the STutility image; pixel_y grows downward -> y axis is inverted when drawn)
rng_x = max(d.pixel_x.max() - d.pixel_x.min() for d in sub.values())
rng_y = max(d.pixel_y.max() - d.pixel_y.min() for d in sub.values())
WIN = 1.04 * max(rng_x, rng_y)

# nearest-neighbour spot distance (pixels) = 100 um Visium pitch -> marker diameter and scale bar
nn_px = np.median([np.median(cKDTree(d[["pixel_x", "pixel_y"]].values)
                             .query(d[["pixel_x", "pixel_y"]].values, k=2)[0][:, 1]) for d in sub.values()])
PX_PER_MM = nn_px * 1000.0 / SPOT_PITCH_UM          # 1 mm = 10 spot pitches

n_spots = df.groupby("section.name").size()
n_ca3 = df[df.is_ca3].groupby("section.name").size()
med_genes = df.groupby("section.name").nFeature_SCT.median()


def section_stats(gene):
    """per-section descriptive statistics in cluster-8 spots (all 12 sections)."""
    rows = []
    for s, d in df.groupby("section.name"):
        c = d[d.is_ca3]
        rows.append(dict(
            gene=gene, section=s, animal_id=d.sample_name_cond.iloc[0].split("_")[0],
            condition=d.sample_condition.iloc[0], n_spots=len(d), n_ca3_spots=len(c),
            median_nFeature_SCT=d.nFeature_SCT.median(),
            ca3_mean_SCT_data=c[f"{gene}_data"].mean(),
            ca3_median_SCT_data=c[f"{gene}_data"].median(),
            ca3_frac_expressing=(c[f"{gene}_counts"] > 0).mean(),
            ca3_mean_SCT_counts=c[f"{gene}_counts"].mean(),
            section_p99_SCT_data=np.percentile(d[f"{gene}_data"], PCT_HIGH),
            section_mean_SCT_data=d[f"{gene}_data"].mean(),
            in_figure=s in sec_names,
        ))
    return pd.DataFrame(rows)


N_GENES_SCT = 14630                  # genes in the SCT assay (rows of the SCT matrix)
ANIMAL_LABEL = {"CF1": "F1", "CF2": "F2", "CF7": "F3", "CG7": "G1", "CG8": "G2", "CG9": "G3"}
ANIMALS_F = ["CF1", "CF2", "CF7"]
ANIMALS_G = ["CG7", "CG8", "CG9"]


def _ols_adjust(y, X):
    """OLS y ~ 1 + X ; returns (adjusted = residual + intercept, coefficients)."""
    A = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ beta
    return resid + beta[0], beta


def animal_adjusted(gene):
    """Detection-rate-adjusted per-animal statistics.

    Uses ALL cluster-8 spots of all 12 sections (both sections of each animal).
    CDR = nFeature_SCT / N_GENES_SCT (fraction of SCT genes detected in the spot),
    z-scored across the CA3 spots.  Model 1: SCT_data ~ 1 + cdr_z.
    Model 2: SCT_data ~ 1 + cdr_z + log(nCount_SCT) (centred).
    Adjusted expression = residual + intercept; then averaged per animal.
    Returns (per-animal DataFrame, dict of model coefficients).
    """
    c = df[df.is_ca3].copy()
    c["animal_id"] = c.sample_name_cond.str.split("_").str[0]
    y = c[f"{gene}_data"].values.astype(float)
    cdr = c.nFeature_SCT.values / N_GENES_SCT
    cdr_z = (cdr - cdr.mean()) / cdr.std(ddof=1)
    logd = np.log(c.nCount_SCT.values.astype(float))
    logd_c = logd - logd.mean()
    c["adj_cdr"], b1 = _ols_adjust(y, cdr_z[:, None])
    c["adj_cdr_depth"], b2 = _ols_adjust(y, np.column_stack([cdr_z, logd_c]))
    c["expressing"] = c[f"{gene}_counts"] > 0
    per = (c.groupby(["animal_id", "sample_condition"])
             .agg(n_ca3_spots=("expressing", "size"),
                  raw_mean=(f"{gene}_data", "mean"),
                  adj_cdr_mean=("adj_cdr", "mean"),
                  adj_cdr_depth_mean=("adj_cdr_depth", "mean"),
                  frac_expressing=("expressing", "mean"),
                  mean_cdr=("nFeature_SCT", lambda v: (v / N_GENES_SCT).mean()))
             .reset_index().rename(columns={"sample_condition": "condition"}))
    per["label"] = per.animal_id.map(ANIMAL_LABEL)
    per["gene"] = gene
    coefs = dict(cdr_intercept=b1[0], cdr_slope=b1[1],
                 cdrdepth_intercept=b2[0], cdrdepth_slope_cdr=b2[1], cdrdepth_slope_logdepth=b2[2],
                 r_cdr_expr=np.corrcoef(cdr_z, y)[0, 1])
    return per, coefs


# ----------------------------------------------------------------------------
# Drawing helpers
# ----------------------------------------------------------------------------
def marker_area(ax, fig, scale=0.9):
    """scatter `s` (points^2) so that a marker spans ~`scale` x the spot pitch.
    fig.canvas.draw() / ax.apply_aspect() must have been called so the axes bbox is final."""
    bbox = ax.get_window_extent(renderer=fig.canvas.get_renderer())
    ax_w_pt = bbox.width * 72.0 / fig.dpi
    diam_pt = nn_px * scale * ax_w_pt / WIN
    return diam_pt ** 2


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
    ax.text(0.5 * (x0 + x_right), y_bar - 0.035 * WIN, f"{SCALE_BAR_MM:g} mm", ha="center", va="bottom",
            fontsize=FS, zorder=5)


def draw_map(ax, fig, d, gene, vmax, mode, label, scale_bar=False, bottom_strip=0.0):
    """mode = 'whole' or 'ca3'."""
    v = d[f"{gene}_{EXPR_SLOT}"].values
    ca3 = d.is_ca3.values
    xy = d[["pixel_x", "pixel_y"]].values
    cx = 0.5 * (d.pixel_x.min() + d.pixel_x.max())
    cy = 0.5 * (d.pixel_y.min() + d.pixel_y.max())
    ax.set_xlim(cx - WIN / 2, cx + WIN / 2)
    # image coordinates: y down; extra room above for the label, below (row 3) for the scale bar
    ax.set_ylim(cy + WIN / 2 + bottom_strip * WIN, cy - WIN / 2 - LABEL_STRIP * WIN)
    ax.set_aspect("equal")
    ax.apply_aspect()                                            # final bbox before sizing markers
    s = marker_area(ax, fig)
    norm = colors.Normalize(vmin=0, vmax=vmax)
    if mode == "whole":
        ax.scatter(xy[:, 0], xy[:, 1], c=v, cmap=CMAP, norm=norm, s=s, linewidths=0, rasterized=True, zorder=1)
    else:
        ax.scatter(xy[~ca3, 0], xy[~ca3, 1], color=COL_CONTEXT, s=s, linewidths=0, rasterized=True, zorder=1)
        ax.scatter(xy[ca3, 0], xy[ca3, 1], c=v[ca3], cmap=CMAP, norm=norm, s=s, linewidths=0,
                   rasterized=True, zorder=2)
    for poly in ca3_outline(xy[ca3]):
        for ring in [poly.exterior, *poly.interiors]:
            rx, ry = ring.xy
            ax.plot(rx, ry, color="black", linewidth=LW, solid_joinstyle="round", zorder=3)
    if scale_bar:
        draw_scale_bar(ax, cx + WIN / 2 - 0.03 * WIN, cy + WIN / 2 + (bottom_strip - 0.03) * WIN)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.text(0.0, 1.0, label, transform=ax.transAxes, ha="left", va="top", fontsize=FS, linespacing=1.15)


def draw_points(ax, gene, flight_vals, ground_vals, labels_f, labels_g, ylabel):
    """3 flight (red) vs 3 ground (blue) per-animal values; mean +/- SD (n = 3 animals)
    as a 1 pt error bar with 1 pt caps next to the points; y axis starts at 0."""
    xs = {"flight": 0.0, "ground": 1.5}
    jitter = np.array([-0.45, -0.15, 0.15])          # the three animals
    x_err = 0.5                                      # error bar to the right of the points
    for cond, vals, labs, col in (("flight", flight_vals, labels_f, COL_FLIGHT),
                                  ("ground", ground_vals, labels_g, COL_GROUND)):
        vals = np.asarray(vals, float)
        ax.scatter(xs[cond] + jitter, vals, s=16, color=col, linewidths=0, zorder=3)
        for lab, x, y in zip(labs, xs[cond] + jitter, vals):
            ax.annotate(lab, (x, y), xytext=(0, 3), textcoords="offset points",
                        ha="center", va="bottom", fontsize=FS, color="#333333")
        ax.errorbar(xs[cond] + x_err, vals.mean(), yerr=vals.std(ddof=1), fmt="none", ecolor=col,
                    elinewidth=LW, capsize=2.5, capthick=LW, zorder=2)
        ax.hlines(vals.mean(), xs[cond] + x_err - 0.12, xs[cond] + x_err + 0.12, color=col, linewidth=LW, zorder=2)
    coef, fdr = MAST[gene]
    ax.set_title(f"MAST coef {coef:+.2f}, FDR = {fdr:.3f}".replace("-", "−"),   # true minus sign
                 fontsize=FS, pad=3, loc="center")
    ax.set_xticks([xs["flight"], xs["ground"]])
    ax.set_xticklabels([sc("Spaceflight"), sc("Ground control")], fontsize=FS)
    ax.set_xlim(-0.7, xs["ground"] + 0.75)
    allv = np.concatenate([flight_vals, ground_vals]).astype(float)
    ax.set_ylim(0, allv.max() * 1.3)                 # never truncated: starts at 0
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.tick_params(axis="y", labelsize=FS, length=3, width=LW, pad=1.5)
    ax.tick_params(axis="x", length=0, pad=2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(LW)
    ax.set_ylabel(ylabel, fontsize=FS, labelpad=2)


def draw_quant(ax, gene, per):
    """main-figure quantitative panel: CDR-adjusted per-animal means (both sections pooled)."""
    p = per.set_index("animal_id")
    draw_points(ax, gene,
                [p.loc[a, "adj_cdr_mean"] for a in ANIMALS_F], [p.loc[a, "adj_cdr_mean"] for a in ANIMALS_G],
                [ANIMAL_LABEL[a] for a in ANIMALS_F], [ANIMAL_LABEL[a] for a in ANIMALS_G],
                "CDR-adjusted mean\nSCT log-expression")          # CA3 spots, both sections: legend


def draw_quant_raw(ax, gene, stats, per, which):
    """supplementary raw panels: which = 'section' (one section per animal, as drawn in the
    maps) or 'pooled' (all cluster-8 spots of both sections of each animal)."""
    if which == "section":
        st = stats.set_index("section")
        fv = [st.loc[s, "ca3_mean_SCT_data"] for s, _, _ in FLIGHT]
        gv = [st.loc[s, "ca3_mean_SCT_data"] for s, _, _ in GROUND]
        ylabel = "Mean SCT log-expression\n(CA3 spots, figure section)"
    else:
        p = per.set_index("animal_id")
        fv = [p.loc[a, "raw_mean"] for a in ANIMALS_F]; gv = [p.loc[a, "raw_mean"] for a in ANIMALS_G]
        ylabel = "Mean SCT log-expression\n(CA3 spots, both sections)"
    draw_points(ax, gene, fv, gv, [l for _, l, _ in FLIGHT], [l for _, l, _ in GROUND], ylabel)


# ----------------------------------------------------------------------------
# Figure builder
# ----------------------------------------------------------------------------
def build_figure(fig_id, genes, all_stats, all_per):
    n = len(genes)
    fig = plt.figure(figsize=srep.figsize_mm("ca3_3" if n == 3 else "ca3_2"))
    fw_mm, fh_mm = fig.get_size_inches() / MM
    left_mm = 9.0                                        # rotated row labels
    # rows: gene header | whole-section block | CA3 block | colour bar | quant
    # row heights in mm: gene header | whole-section block | CA3 block | colour bar | quant row
    quant_mm, cbar_mm, head_mm = 38.0, 11.0, 3.0
    block_mm = 0.5 * (fh_mm - 2.5 - head_mm - cbar_mm - quant_mm - 4.0)
    outer = gridspec.GridSpec(5, n, figure=fig,
                              height_ratios=[head_mm, block_mm, block_mm, cbar_mm, quant_mm],
                              left=left_mm / fw_mm, right=1 - 1.0 / fw_mm,
                              bottom=1.5 / fh_mm, top=1 - 1.0 / fh_mm,
                              hspace=0.02, wspace=0.14)
    fig.canvas.draw()
    for j, gene in enumerate(genes):
        vmax = {s: max(np.percentile(sub[s][f"{gene}_{EXPR_SLOT}"], PCT_HIGH), 1e-6) for s in sec_names}
        # gene header (italic, same 8 pt as everything else, no bold)
        axh = fig.add_subplot(outer[0, j]); axh.axis("off")
        axh.text(0.5, 0.5, gene, ha="center", va="center", fontsize=FS, fontstyle="italic")
        for r, mode in ((1, "whole"), (2, "ca3")):
            inner = gridspec.GridSpecFromSubplotSpec(4, 2, subplot_spec=outer[r, j],
                                                     height_ratios=[0.12, 1 + LABEL_STRIP, 1 + LABEL_STRIP,
                                                                    1 + LABEL_STRIP + BAR_STRIP],
                                                     hspace=0.02, wspace=0.03)
            for k, head in enumerate(("Spaceflight", "Ground control")):
                axc = fig.add_subplot(inner[0, k]); axc.axis("off")
                axc.text(0.5, 0.0, sc(head), ha="center", va="bottom", fontsize=FS)
            for i in range(3):
                for k, lst in enumerate((FLIGHT, GROUND)):
                    sec, lab, orig = lst[i]
                    ax = fig.add_subplot(inner[i + 1, k])
                    if mode == "whole":
                        label = f"{lab} ({orig})\n{srep.fmt_int(int(med_genes[sec]))} genes/spot"
                    else:
                        label = f"{lab} ({orig})\nCA3 n = {srep.fmt_int(n_ca3[sec])}"
                    draw_map(ax, fig, sub[sec], gene, vmax[sec], mode, label,
                             scale_bar=(i == 2 and k == 0), bottom_strip=BAR_STRIP if i == 2 else 0.0)
        bq = outer[4, j].get_position(fig)
        axq = fig.add_axes([bq.x0 + 13.0 / fw_mm, bq.y0 + 4.0 / fh_mm, bq.width - 13.0 / fw_mm,
                            bq.height - 4.0 / fh_mm - 4.0 / fh_mm])   # room for x tick labels / title
        draw_quant(axq, gene, all_per[gene])

    # row labels (left margin)
    for r, txt in ((1, "Whole section\n(all spots; CA3 = black outline)"),
                   (2, "CA3 only\n(cluster 8 colored; other spots gray)")):
        bb = outer[r, 0].get_position(fig)
        fig.text(3.5 / fw_mm, 0.5 * (bb.y0 + bb.y1), sc(txt), rotation=90, ha="center", va="center", fontsize=FS)

    # shared colour bar (one label line)
    bb0 = outer[3, 0].get_position(fig); bb1 = outer[3, n - 1].get_position(fig)
    cw = 40.0 / fw_mm
    cax = fig.add_axes([0.5 * (bb0.x0 + bb1.x1) - cw / 2, bb0.y0 + 4.5 / fh_mm, cw, 2.5 / fh_mm])
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=CMAP, norm=colors.Normalize(0, 1),
                                          orientation="horizontal", ticks=[0, 1])
    cb.ax.set_xticklabels(["Low", "High"], fontsize=FS)   # 0 and the section's 99th percentile (legend)
    cb.ax.tick_params(length=0, pad=1.5)
    cb.outline.set_linewidth(LW)
    cax.set_title(sc("SCT log-normalized expression (scaled within each section)"), fontsize=FS, pad=2, loc="center")

    base = os.path.join(HERE, f"{fig_id}_{'_'.join(genes)}_CA3")
    srep.save(fig, base)
    plt.close(fig)
    return base


def build_raw_figure(fig_id, genes, all_stats, all_per):
    """Supplementary: unadjusted per-animal cluster-8 means (figure section only, and both
    sections pooled) so the raw numbers behind the CDR-adjusted panel stay visible."""
    n = len(genes)
    fig = plt.figure(figsize=srep.figsize_mm("raw_means3" if n == 3 else "raw_means2"))
    fw_mm, fh_mm = fig.get_size_inches() / MM
    gs = gridspec.GridSpec(3, n, figure=fig, height_ratios=[1.0, 30.0, 30.0],
                           left=14.0 / fw_mm, right=1 - 1.0 / fw_mm, bottom=5.0 / fh_mm, top=1 - 1.5 / fh_mm,
                           hspace=0.42, wspace=0.5)
    fig.canvas.draw()
    for j, gene in enumerate(genes):
        axh = fig.add_subplot(gs[0, j]); axh.axis("off")
        axh.text(0.5, 0.5, gene, ha="center", va="center", fontsize=FS, fontstyle="italic")
        ax1 = fig.add_subplot(gs[1, j]); draw_quant_raw(ax1, gene, all_stats[gene], all_per[gene], "section")
        ax2 = fig.add_subplot(gs[2, j]); draw_quant_raw(ax2, gene, all_stats[gene], all_per[gene], "pooled")
        ax2.set_title("", loc="center")                    # MAST annotation once per gene
    base = os.path.join(HERE, f"{fig_id}_{'_'.join(genes)}_CA3_raw_means")
    srep.save(fig, base)
    plt.close(fig)
    return base


# ----------------------------------------------------------------------------
if __name__ == "__main__":
    genes_all = [g for gl in PANELS.values() for g in gl]
    all_stats = {g: section_stats(g) for g in genes_all}
    stats_df = pd.concat(all_stats.values(), ignore_index=True)
    stats_df.to_csv(os.path.join(HERE, "CA3_cluster8_per_section_stats.csv"), index=False)

    print("Spots per section / cluster-8 spots per section / median genes per spot:")
    print(pd.DataFrame({"n_spots": n_spots, "n_cluster8": n_ca3, "median_nFeature_SCT": med_genes}))
    print(f"\nSpot pitch (median nearest-neighbour distance) = {nn_px:.2f} px = {SPOT_PITCH_UM:g} um; "
          f"1 mm = {PX_PER_MM:.1f} px; map window = {WIN:.0f} px")

    print("\nPer-animal cluster-8 mean SCT log-expression (figure sections only) vs MAST direction:")
    rows = []
    for g in genes_all:
        st = all_stats[g]
        f = st[st.in_figure & (st.condition == "flight")]
        gr = st[st.in_figure & (st.condition == "ground")]
        fm, gm = f.ca3_mean_SCT_data.mean(), gr.ca3_mean_SCT_data.mean()
        coef, fdr = MAST[g]
        # pooled version using both sections of each animal
        st2 = st.groupby(["animal_id", "condition"])[["ca3_mean_SCT_data", "n_ca3_spots"]].apply(
            lambda d: np.average(d.ca3_mean_SCT_data, weights=d.n_ca3_spots)).reset_index(name="m")
        fm2 = st2[st2.condition == "flight"].m.mean(); gm2 = st2[st2.condition == "ground"].m.mean()
        rows.append(dict(
            gene=g, MAST_coef=coef, MAST_FDR=fdr,
            **{f"{lab}({orig})": round(st.set_index('section').loc[s, 'ca3_mean_SCT_data'], 3)
               for s, lab, orig in SECTIONS},
            flight_mean=round(fm, 3), ground_mean=round(gm, 3), diff_F_minus_G=round(fm - gm, 3),
            agrees_with_MAST=np.sign(fm - gm) == np.sign(coef),
            frac_expr_flight=round(f.ca3_frac_expressing.mean(), 3),
            frac_expr_ground=round(gr.ca3_frac_expressing.mean(), 3),
            diff_both_sections_per_animal=round(fm2 - gm2, 3),
            agrees_both_sections=np.sign(fm2 - gm2) == np.sign(coef),
        ))
    summ = pd.DataFrame(rows)
    pd.set_option("display.width", 300); pd.set_option("display.max_columns", 80)
    print(summ.to_string(index=False))

    print("\nFraction of cluster-8 spots with >=1 SCT count, per figure section:")
    print(stats_df[stats_df.in_figure].pivot(index="gene", columns="section", values="ca3_frac_expressing").round(2))

    # ---- detection-rate-adjusted per-animal values (all 12 sections, pooled per animal) ----
    all_per, all_coefs, adj_rows = {}, {}, []
    for g in genes_all:
        per, coefs = animal_adjusted(g)
        all_per[g], all_coefs[g] = per, coefs
        p = per.set_index("animal_id")
        coef, fdr = MAST[g]
        rec = dict(gene=g, MAST_coef=coef, MAST_FDR=fdr,
                   cdr_slope_per_SD=round(coefs["cdr_slope"], 4), r_cdr_vs_expr=round(coefs["r_cdr_expr"], 3))
        for a in ANIMALS_F + ANIMALS_G:
            rec[f"adj_{ANIMAL_LABEL[a]}({a})"] = round(p.loc[a, "adj_cdr_mean"], 3)
        fa = np.mean([p.loc[a, "adj_cdr_mean"] for a in ANIMALS_F]); ga = np.mean([p.loc[a, "adj_cdr_mean"] for a in ANIMALS_G])
        fa2 = np.mean([p.loc[a, "adj_cdr_depth_mean"] for a in ANIMALS_F]); ga2 = np.mean([p.loc[a, "adj_cdr_depth_mean"] for a in ANIMALS_G])
        fr = np.mean([p.loc[a, "raw_mean"] for a in ANIMALS_F]); gr_ = np.mean([p.loc[a, "raw_mean"] for a in ANIMALS_G])
        rec.update(adj_flight_mean=round(fa, 3), adj_ground_mean=round(ga, 3), adj_diff_F_minus_G=round(fa - ga, 3),
                   adj_agrees_with_MAST=np.sign(fa - ga) == np.sign(coef),
                   adj2_cdr_depth_diff=round(fa2 - ga2, 3), adj2_agrees_with_MAST=np.sign(fa2 - ga2) == np.sign(coef),
                   raw_pooled_flight_mean=round(fr, 3), raw_pooled_ground_mean=round(gr_, 3),
                   frac_expr_flight_pooled=round(np.mean([p.loc[a, "frac_expressing"] for a in ANIMALS_F]), 3),
                   frac_expr_ground_pooled=round(np.mean([p.loc[a, "frac_expressing"] for a in ANIMALS_G]), 3))
        adj_rows.append(rec)
    adj = pd.DataFrame(adj_rows)
    print("\nCDR-adjusted per-animal cluster-8 means (all cluster-8 spots, both sections per animal pooled):")
    print(adj.to_string(index=False))
    print("\nModel coefficients (expression ~ 1 + cdr_z [+ log nCount_SCT]), all cluster-8 spots:")
    print(pd.DataFrame(all_coefs).T.round(4).to_string())
    per_df = pd.concat(all_per.values(), ignore_index=True)
    per_df["n_sections"] = 2
    per_df.to_csv(os.path.join(HERE, "CA3_cluster8_per_animal_adjusted.csv"), index=False)
    print("\nPer-animal CDR (fraction of SCT genes detected, mean over cluster-8 spots):")
    print(per_df[per_df.gene == genes_all[0]][["animal_id", "label", "condition", "n_ca3_spots", "mean_cdr"]].to_string(index=False))

    # mean +/- SD (n = 3 animals per group) of the values drawn in the quantitative panels
    sd_rows = []
    for g in genes_all:
        p = all_per[g].set_index("animal_id")
        for col, name in (("adj_cdr_mean", "CDR_adjusted"), ("raw_mean", "raw_pooled")):
            fv = np.array([p.loc[a, col] for a in ANIMALS_F]); gv = np.array([p.loc[a, col] for a in ANIMALS_G])
            sd_rows.append(dict(gene=g, measure=name, flight_mean=fv.mean(), flight_sd=fv.std(ddof=1),
                                ground_mean=gv.mean(), ground_sd=gv.std(ddof=1), n_per_group=3))
    print("\nMean +/- SD across the three animals per group (as drawn by the error bars):")
    print(pd.DataFrame(sd_rows).round(3).to_string(index=False))

    summ = summ.merge(adj.drop(columns=["MAST_coef", "MAST_FDR"]), on="gene")
    summ.to_csv(os.path.join(HERE, "CA3_cluster8_per_animal_summary.csv"), index=False)

    for fig_id, genes in PANELS.items():
        base = build_figure(fig_id, genes, all_stats, all_per)
        print("wrote", base + ".png / .pdf")
        base = build_raw_figure(fig_id, genes, all_stats, all_per)
        print("wrote", base + ".png / .pdf")
