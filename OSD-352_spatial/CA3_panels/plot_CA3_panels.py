#!/usr/bin/env python3
"""
plot_CA3_panels.py

Builds the spatial-transcriptomics figure panels (Fig2c, Fig3c, Fig4c) for the
hippocampal CA3 (ST cluster 8) story from the CSVs exported by
extract_spatial_data.R (spots.csv, expr.csv in the same folder).

For every gene, one column with three rows:
  row 1  whole-section maps  : all spots of the section, colour = SCT log1p
                               expression scaled PER SECTION (0 -> section's
                               99th percentile); cluster-8 (CA3) spots outlined
                               in black.
  row 2  CA3-only maps       : only cluster-8 spots coloured (same per-section
                               scaling), every other spot in light grey.
  row 3  quantitative panel  : per-animal mean of detection-rate (CDR)-adjusted
                               SCT log-expression over ALL cluster-8 spots of
                               both sections of the animal (3 flight animals,
                               red vs 3 ground animals, blue), mean bar, MAST
                               coef / FDR from Masarapu et al. annotated.
                               CDR adjustment: expression ~ 1 + z(nFeature_SCT /
                               14630) fitted by OLS across all cluster-8 spots;
                               adjusted = residual + intercept.

Map layout inside each block: 3 rows (animals 1..3) x 2 columns (Flight | Ground).

Supplementary Fig*_CA3_raw_means.png/.pdf: the unadjusted per-animal means
(figure section only, and both sections pooled).

Outputs (PNG + PDF, 300 dpi) and CSVs with per-section / per-animal statistics.
Run:  /opt/anaconda3/bin/python3 plot_CA3_panels.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec, colors
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------------
# Design decisions / constants
# ----------------------------------------------------------------------------
CA3_CLUSTER = 8                      # seurat_clusters value for hippocampal CA3
PCT_HIGH = 99                        # per-section upper limit of colour scale
EXPR_SLOT = "data"                   # SCT "data" = log1p(SCT-corrected counts)

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

COL_FLIGHT = "#c0392b"
COL_GROUND = "#2b6cb0"
COL_CONTEXT = "#dcdcdc"              # non-CA3 spots on the CA3-only maps
CMAP = colors.LinearSegmentedColormap.from_list(
    "expr", plt.get_cmap("Reds")(np.linspace(0.07, 1.0, 256)))   # faint pink -> dark red

MM = 1 / 25.4
COL_W_MM = 60.0                      # width of one gene column (3 genes -> 180 mm)
FS_GENE, FS_HEAD, FS_LAB, FS_TICK, FS_FOOT = 9, 7, 6, 6, 5.5

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": FS_TICK,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300,
})

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

# nearest-neighbour spot distance (pixels) -> marker diameter
nn_px = np.median([np.median(cKDTree(d[["pixel_x", "pixel_y"]].values)
                             .query(d[["pixel_x", "pixel_y"]].values, k=2)[0][:, 1]) for d in sub.values()])

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
    """scatter `s` (points^2) so that a marker spans ~`scale` x the spot pitch."""
    bbox = ax.get_window_extent(renderer=fig.canvas.get_renderer())
    ax_w_pt = bbox.width * 72.0 / fig.dpi
    diam_pt = nn_px * scale * ax_w_pt / WIN
    return diam_pt ** 2


def draw_map(ax, fig, d, gene, vmax, mode, label):
    """mode = 'whole' or 'ca3'."""
    v = d[f"{gene}_{EXPR_SLOT}"].values
    ca3 = d.is_ca3.values
    s = marker_area(ax, fig)
    norm = colors.Normalize(vmin=0, vmax=vmax)
    if mode == "whole":
        ax.scatter(d.pixel_x[~ca3], d.pixel_y[~ca3], c=v[~ca3], cmap=CMAP, norm=norm,
                   s=s, linewidths=0, rasterized=True, zorder=1)
        ax.scatter(d.pixel_x[ca3], d.pixel_y[ca3], c=v[ca3], cmap=CMAP, norm=norm,
                   s=s, edgecolors="black", linewidths=0.13, rasterized=True, zorder=2)
    else:
        ax.scatter(d.pixel_x[~ca3], d.pixel_y[~ca3], color=COL_CONTEXT, s=s, linewidths=0,
                   rasterized=True, zorder=1)
        ax.scatter(d.pixel_x[ca3], d.pixel_y[ca3], c=v[ca3], cmap=CMAP, norm=norm,
                   s=s, edgecolors="black", linewidths=0.12, rasterized=True, zorder=2)
    cx = 0.5 * (d.pixel_x.min() + d.pixel_x.max())
    cy = 0.5 * (d.pixel_y.min() + d.pixel_y.max())
    ax.set_xlim(cx - WIN / 2, cx + WIN / 2)
    ax.set_ylim(cy + WIN / 2, cy - WIN / 2 - 0.16 * WIN)   # image coordinates: y down; extra room for label
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.text(0.01, 0.99, label, transform=ax.transAxes, ha="left", va="top",
            fontsize=FS_LAB, linespacing=1.1)


def draw_points(ax, gene, flight_vals, ground_vals, labels_f, labels_g, ylabel, clip_zero=True):
    """3 flight (red) vs 3 ground (blue) per-animal values with a mean bar."""
    xs = {"flight": 0, "ground": 1}
    jitter = np.array([-0.18, 0.0, 0.18])
    for cond, vals, labs, col in (("flight", flight_vals, labels_f, COL_FLIGHT),
                                  ("ground", ground_vals, labels_g, COL_GROUND)):
        vals = np.asarray(vals, float)
        ax.scatter(xs[cond] + jitter, vals, s=22, color=col, edgecolors="white", linewidths=0.5, zorder=3)
        for lab, x, y in zip(labs, xs[cond] + jitter, vals):
            ax.annotate(lab, (x, y), xytext=(0, 4), textcoords="offset points",
                        ha="center", va="bottom", fontsize=FS_LAB - 1, color="#333333")
        ax.hlines(vals.mean(), xs[cond] - 0.32, xs[cond] + 0.32, color=col, linewidth=1.4, zorder=2)
    coef, fdr = MAST[gene]
    ax.set_title(f"MAST coef {coef:+.2f}, FDR {fdr:.3f}", fontsize=FS_LAB, color="#333333", pad=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Flight", "Ground"], fontsize=FS_HEAD)
    ax.set_xlim(-0.6, 1.6)
    allv = np.concatenate([flight_vals, ground_vals]).astype(float)
    lo, hi = allv.min(), allv.max()
    pad = 0.35 * max(hi - lo, 0.1)
    ax.set_ylim(max(0, lo - pad) if clip_zero else lo - pad, hi + pad)
    ax.tick_params(axis="y", labelsize=FS_TICK, length=2, pad=1.5)
    ax.tick_params(axis="x", length=0, pad=2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_ylabel(ylabel, fontsize=FS_LAB, labelpad=2)


def draw_quant(ax, gene, per):
    """main-figure quantitative panel: CDR-adjusted per-animal means (both sections pooled)."""
    p = per.set_index("animal_id")
    draw_points(ax, gene,
                [p.loc[a, "adj_cdr_mean"] for a in ANIMALS_F], [p.loc[a, "adj_cdr_mean"] for a in ANIMALS_G],
                [ANIMAL_LABEL[a] for a in ANIMALS_F], [ANIMAL_LABEL[a] for a in ANIMALS_G],
                "CDR-adjusted mean\nSCT log-expr. (CA3 spots)", clip_zero=False)


def draw_quant_raw(ax, gene, stats, per, which):
    """supplementary raw panels: which = 'section' (one section per animal, as drawn in the
    maps) or 'pooled' (all cluster-8 spots of both sections of each animal)."""
    if which == "section":
        st = stats.set_index("section")
        fv = [st.loc[s, "ca3_mean_SCT_data"] for s, _, _ in FLIGHT]
        gv = [st.loc[s, "ca3_mean_SCT_data"] for s, _, _ in GROUND]
        ylabel = "Raw mean SCT log-expr.\n(CA3 spots, figure section)"
    else:
        p = per.set_index("animal_id")
        fv = [p.loc[a, "raw_mean"] for a in ANIMALS_F]; gv = [p.loc[a, "raw_mean"] for a in ANIMALS_G]
        ylabel = "Raw mean SCT log-expr.\n(CA3 spots, both sections)"
    draw_points(ax, gene, fv, gv, [l for _, l, _ in FLIGHT], [l for _, l, _ in GROUND], ylabel)


# ----------------------------------------------------------------------------
# Figure builder
# ----------------------------------------------------------------------------
def build_figure(fig_id, genes, all_stats, all_per):
    n = len(genes)
    fig_w = COL_W_MM * n + 14                            # + left margin for row labels
    fig_h = 252 + (3 - n) * 18
    fig = plt.figure(figsize=(fig_w * MM, fig_h * MM))
    # rows: gene header | whole-section block | CA3 block | colourbar | quant | footnote
    outer = gridspec.GridSpec(6, n, figure=fig,
                              height_ratios=[0.06, 1.0, 1.0, 0.12, 0.52, 0.17 * 3 / n],
                              left=14 / fig_w, right=1 - 3 / fig_w, bottom=0.005, top=0.985,
                              hspace=0.08, wspace=0.10)
    for j, gene in enumerate(genes):
        stats = all_stats[gene]
        vmax = {s: max(np.percentile(sub[s][f"{gene}_{EXPR_SLOT}"], PCT_HIGH), 1e-6) for s in sec_names}
        # gene header
        axh = fig.add_subplot(outer[0, j]); axh.axis("off")
        axh.text(0.5, 0.5, gene, ha="center", va="center", fontsize=FS_GENE, fontweight="bold",
                 fontstyle="italic")
        for r, mode in ((1, "whole"), (2, "ca3")):
            inner = gridspec.GridSpecFromSubplotSpec(4, 2, subplot_spec=outer[r, j],
                                                     height_ratios=[0.09, 1, 1, 1],
                                                     hspace=0.02, wspace=0.02)
            for k, head in enumerate(("Flight", "Ground")):
                axc = fig.add_subplot(inner[0, k]); axc.axis("off")
                axc.text(0.5, 0.1, head, ha="center", va="bottom", fontsize=FS_HEAD,
                         color=COL_FLIGHT if k == 0 else COL_GROUND, fontweight="bold")
            for i in range(3):
                for k, lst in enumerate((FLIGHT, GROUND)):
                    sec, lab, orig = lst[i]
                    ax = fig.add_subplot(inner[i + 1, k])
                    label = f"{lab} ({orig})\nmedian {int(med_genes[sec])} genes/spot"
                    if mode == "ca3":
                        label = f"{lab} ({orig})\nCA3 n = {int(n_ca3[sec])}"
                    draw_map(ax, fig, sub[sec], gene, vmax[sec], mode, label)
        axq = fig.add_subplot(outer[4, j])
        draw_quant(axq, gene, all_per[gene])

    # row labels (left margin)
    for r, txt in ((1, "Whole section\n(all spots; CA3 = black outline)"),
                   (2, "CA3 only\n(cluster 8 coloured; other spots grey)")):
        bb = outer[r, 0].get_position(fig)
        fig.text(3.0 / fig_w, 0.5 * (bb.y0 + bb.y1), txt, rotation=90, ha="center", va="center",
                 fontsize=FS_HEAD)

    # shared colourbar
    bb0 = outer[3, 0].get_position(fig); bb1 = outer[3, n - 1].get_position(fig)
    cw = 0.45 * (bb1.x1 - bb0.x0)
    cax = fig.add_axes([0.5 * (bb0.x0 + bb1.x1) - cw / 2, bb0.y0 + 0.62 * (bb0.y1 - bb0.y0),
                        cw, 0.17 * (bb0.y1 - bb0.y0)])
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=CMAP, norm=colors.Normalize(0, 1),
                                          orientation="horizontal", ticks=[0, 1])
    cb.ax.set_xticklabels(["Low (0)", f"High (section {PCT_HIGH}th pct.)"], fontsize=FS_TICK)
    cb.ax.tick_params(length=2, pad=1.5)
    cb.outline.set_linewidth(0.4)
    cax.set_title("SCT log-normalised expression, scaled per section", fontsize=FS_LAB, pad=2)

    # footnote
    foot = (f"Colour: Seurat SCT assay, 'data' slot = log1p(SCT-corrected counts), rescaled within each section "
            f"from 0 to that section's {PCT_HIGH}th percentile so sequencing-depth differences between sections "
            f"do not dominate. CA3 = ST cluster {CA3_CLUSTER}\n('Synapses and metabolism in hippocampus CA3 region'). "
            "Maps show one section per animal (original animal ID in parentheses). "
            "Bottom row: points = per-animal mean of CDR-adjusted expression, bar = mean of the 3 animals. "
            "CDR-adjusted = residuals (plus intercept) from expression ~ cellular detection rate "
            "(fraction of SCT genes detected per spot, z-scored) fitted across all CA3 spots, pooled over both "
            "sections per animal; descriptive, the reported statistics (MAST coef / FDR) are from the spot-level "
            "MAST mixed model of Masarapu et al. (sample as random effect) and are not derived from these means.")
    import textwrap
    foot = "\n".join(textwrap.wrap(foot.replace("\n", " "), width=int(0.86 * (fig_w - 17))))
    bb = outer[5, 0].get_position(fig)
    fig.text(14 / fig_w, bb.y1 - 0.002, foot, ha="left", va="top", fontsize=FS_FOOT, linespacing=1.25)

    base = os.path.join(HERE, f"{fig_id}_{'_'.join(genes)}_CA3")
    fig.savefig(base + ".png", dpi=300, facecolor="white")
    fig.savefig(base + ".pdf", dpi=300, facecolor="white")
    plt.close(fig)
    return base


def build_raw_figure(fig_id, genes, all_stats, all_per):
    """Supplementary: unadjusted per-animal cluster-8 means (figure section only, and both
    sections pooled) so the raw numbers behind the CDR-adjusted panel stay visible."""
    n = len(genes)
    fig_w = COL_W_MM * n + 14
    fig_h = 95
    fig = plt.figure(figsize=(fig_w * MM, fig_h * MM))
    gs = gridspec.GridSpec(3, n, figure=fig, height_ratios=[1, 1, 0.3],
                           left=14 / fig_w, right=1 - 3 / fig_w, bottom=0.02, top=0.93,
                           hspace=0.55, wspace=0.10)
    for j, gene in enumerate(genes):
        ax1 = fig.add_subplot(gs[0, j]); draw_quant_raw(ax1, gene, all_stats[gene], all_per[gene], "section")
        ax1.set_title(f"{gene}   " + ax1.get_title(), fontsize=FS_LAB, color="#333333", pad=3)
        ax2 = fig.add_subplot(gs[1, j]); draw_quant_raw(ax2, gene, all_stats[gene], all_per[gene], "pooled")
    foot = ("Unadjusted per-animal means of SCT log-expression over cluster-8 (CA3) spots. Top: the single "
            "section per animal drawn in the maps; bottom: all cluster-8 spots of both sections of each animal. "
            "These raw means track sequencing depth (median genes/spot: CF1 5471/4159, CF2 953/1740, CF7 5651/5367; "
            "CG7 2002/5086, CG8 3324/2582, CG9 2843/2685) and are shown for transparency only; the main figure uses "
            "CDR-adjusted values.")
    import textwrap
    foot = "\n".join(textwrap.wrap(foot, width=int(0.86 * (fig_w - 17))))
    fig.text(14 / fig_w, gs[2, 0].get_position(fig).y1, foot, ha="left", va="top", fontsize=FS_FOOT, linespacing=1.25)
    base = os.path.join(HERE, f"{fig_id}_{'_'.join(genes)}_CA3_raw_means")
    fig.savefig(base + ".png", dpi=300, facecolor="white")
    fig.savefig(base + ".pdf", dpi=300, facecolor="white")
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
    print(f"\nSpot pitch (median nearest-neighbour distance) = {nn_px:.1f} px; map window = {WIN:.0f} px")

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

    summ = summ.merge(adj.drop(columns=["MAST_coef", "MAST_FDR"]), on="gene")
    summ.to_csv(os.path.join(HERE, "CA3_cluster8_per_animal_summary.csv"), index=False)

    for fig_id, genes in PANELS.items():
        base = build_figure(fig_id, genes, all_stats, all_per)
        print("wrote", base + ".png / .pdf")
        base = build_raw_figure(fig_id, genes, all_stats, all_per)
        print("wrote", base + ".png / .pdf")
