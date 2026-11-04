"""
Shared figure style for Scientific Reports submissions.

Every figure script in this project imports this module (a copy sits in each
figure folder so scripts can be run from their own directory) and calls
``apply()`` before creating figures.  It enforces:

* one sans-serif typeface (Arial, falling back to Helvetica) at ONE size (8 pt)
  for every piece of lettering in every figure;
* lettering in sentence case (``sentence_case``);
* no line thinner than 1 pt (axes, ticks, bars, outlines, grid);
* white background, no top/right spines, no legend frames, no decorative effects;
* panel sizes in millimetres so that panels can be assembled 1:1 into a single
  figure file without rescaling (``PANEL_MM``, ``figsize_mm``);
* thousands separators only for numbers with five or more digits (``fmt_int``);
* 300 dpi PNG + vector PDF export with embedded (Type 42) fonts (``save``).
"""
import re
import matplotlib as mpl

FONT_SIZE = 8          # pt, used for ALL lettering in ALL figures
LINE_PT = 1.0          # minimum line width (pt)
MM = 1.0 / 25.4        # inches per mm
DOUBLE_COL_MM = 183.0  # Scientific Reports double-column width
SINGLE_COL_MM = 88.0

# physical panel sizes (width, height) in mm; the assembly script places panels 1:1
PANEL_MM = {
    "volcano":      (118, 100),
    "pca":          (60, 60),
    "nes":          (90, 60),    # harmonized NES bar panels (7, 7 or 3 terms)
    "ora":          (90, 65),    # g:Profiler bar panels
    "maps3":        (183, 112),  # whole-section spatial maps, three genes
    "maps2":        (122, 105),  # whole-section spatial maps, two genes
    "ca3_3":        (183, 235),  # three-row CA3 panel, three genes (supplementary, own figure)
    "ca3_2":        (122, 235),
    "raw_means3":   (183, 90),
    "raw_means2":   (122, 90),
    "qc":           (183, 120),  # QC metric grid
    "volcano_supp": (90, 85),
}

COLORS = {
    "flight": "#c0392b", "ground": "#2e6fbd",
    "up": "#c0392b", "down": "#2e6fbd",
    "ns": "#bdbdbd", "nom_up": "#f1a9a0", "nom_down": "#a9c8ec",
    "curated": "#f39c12", "context": "#dedede",
}


def apply():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Helvetica Neue", "Liberation Sans", "DejaVu Sans"],
        "font.size": FONT_SIZE, "axes.titlesize": FONT_SIZE, "axes.labelsize": FONT_SIZE,
        "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE, "legend.fontsize": FONT_SIZE,
        "legend.title_fontsize": FONT_SIZE, "figure.titlesize": FONT_SIZE,
        "axes.titleweight": "normal", "axes.titlelocation": "left", "axes.titlepad": 4,
        "axes.linewidth": LINE_PT, "lines.linewidth": LINE_PT, "patch.linewidth": LINE_PT,
        "grid.linewidth": LINE_PT, "hatch.linewidth": LINE_PT,
        "xtick.major.width": LINE_PT, "ytick.major.width": LINE_PT,
        "xtick.minor.width": LINE_PT, "ytick.minor.width": LINE_PT,
        "xtick.major.size": 3, "ytick.major.size": 3, "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
        "legend.frameon": False, "legend.handlelength": 1.2, "legend.borderaxespad": 0.3,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "savefig.dpi": 300, "figure.dpi": 100,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "mathtext.default": "regular", "axes.unicode_minus": True,
        "scatter.edgecolors": "none",
    })


def figsize_mm(key_or_w, h=None):
    """figsize in inches from a PANEL_MM key or an explicit (w, h) in mm."""
    if h is None:
        w, h = PANEL_MM[key_or_w]
    else:
        w = key_or_w
    return (w * MM, h * MM)


# tokens that keep their capitalisation inside a sentence-case label
_PROTECTED = {
    "ATP", "NADH", "NAD", "FAD", "GTP", "GTPase", "GO", "BP", "CC", "MF", "KEGG", "NES", "FDR", "DNA", "RNA", "mRNA",
    "rRNA", "ER", "SNARE", "ESCRT", "AMPK", "PI3K", "Akt", "MAPK", "Rap1", "Ras", "Rho", "Arp2/3", "Wnt", "Notch",
    "CA1", "CA3", "ST", "SCT", "CDR", "MAST", "PCA", "PC1", "PC2", "PC3", "RR-3", "RR-10", "RRRM-2", "ISS", "ISS-T",
    "FLT", "GC", "OSD-352", "OSD-562", "OSD-564", "GLDS-352", "GLDS-557", "GLDS-569", "RIN", "QC", "RSeQC",
    "DESeq2", "GSEA", "GSEApy", "ORA", "Visium", "Parkinson", "Alzheimer", "Huntington", "Golgi", "Rab", "Arf",
    "Ca2+", "ATPase", "NAFLD", "H+", "TNT", "TNTs", "sEV", "sEVs", "EV", "EVs", "RNA-seq", "I", "II", "III", "IV", "V",
    "VST", "UMI", "g:Profiler", "g:SCS", "NASA", "OSDR", "GeneLab", "Enrichr", "FoxO", "HIF-1", "mTOR", "TGF-beta",
    "NF-kappa", "cAMP", "cGMP", "PKG", "Jak-STAT", "ErbB", "VEGF", "Fc", "PPAR", "Hippo", "Hedgehog", "Toll-like",
}
_ROMAN = re.compile(r"^[IVX]+$")


def sentence_case(label: str) -> str:
    """Lower-case every word except the first and protected tokens (acronyms,
    gene symbols are passed through unchanged by callers, roman numerals,
    tokens containing digits or a colon such as GO:0006119)."""
    if not label:
        return label
    words = label.split(" ")
    out = []
    for i, w in enumerate(words):
        core = w.strip("(),;:[]")
        keep = (core in _PROTECTED or _ROMAN.match(core) or any(ch.isdigit() for ch in core)
                or ":" in core or (core.isupper() and len(core) > 1))
        if i == 0:
            out.append(w if keep else (w[:1].upper() + w[1:].lower()))
        else:
            out.append(w if keep else w.lower())
    return " ".join(out)


def fmt_int(n) -> str:
    """Thousands separator only for five or more digits (journal style)."""
    n = int(round(n))
    return f"{n:,}" if abs(n) >= 10000 else str(n)


def fdr_label(q) -> str:
    return "FDR < 0.001" if q is None or q < 0.001 else f"FDR = {q:.3f}"


def panel_letter(fig, ax, letter, dx_mm=-6.0, dy_mm=2.0):
    """Bold lower-case panel letter at the top-left of an axes, same size as all other lettering."""
    bbox = ax.get_position()
    fw, fh = fig.get_size_inches()
    x = bbox.x0 + dx_mm * MM / fw
    y = bbox.y1 + dy_mm * MM / fh
    fig.text(x, y, letter, fontsize=FONT_SIZE, fontweight="bold", ha="left", va="bottom")


def save(fig, base, tiff=False):
    """300 dpi PNG + vector PDF (+ optional LZW TIFF) with a white background."""
    fig.savefig(base + ".png", dpi=300, facecolor="white")
    fig.savefig(base + ".pdf", dpi=300, facecolor="white")
    if tiff:
        from PIL import Image
        Image.open(base + ".png").save(base + ".tif", compression="tiff_lzw", dpi=(300, 300))
