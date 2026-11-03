options(repos = c(CRAN = "https://cloud.r-project.org"))

library(Seurat)
library(SeuratObject)
library(STutility)
library(ggplot2)
library(patchwork)
library(cowplot)

setwd("/Users/justinlin/Desktop/Research /TNT sEV Paper/OSD-352/NASA_RR3_Brain-main/manuscript_plots")

st_brain <- readRDS("RR3-brain_ST_final_object.rds")

genes_of_interest <- c("Arpc2", "Rhob", "Sec14l1", "Arap2", "Ubl3", "Epha4", "Pink1", "Ly6h", "Snph", "Itpr1")
DefaultAssay(st_brain) <- "SCT"
st_brain@meta.data$no_label <- ""  # empty label to suppress auto section index

# Sections in exact top-to-bottom display order (section.name uses "Sample_" prefix)
flight_sections <- c("Sample_158_A1", "Sample_158_C1", "Sample_159_C1")
ground_sections <- c("Sample_159_B1", "Sample_304_A1", "Sample_304_C1")

# Pre-subset one object per section to guarantee plot order
flight_objs <- lapply(flight_sections, function(s)
  SubsetSTData(st_brain, expression = section.name == s))
ground_objs <- lapply(ground_sections, function(s)
  SubsetSTData(st_brain, expression = section.name == s))

# Collect gene panels for combined figure
all_gene_panels <- list()

for (gene in genes_of_interest) {
  shared_max <- max(GetAssayData(st_brain, assay = "SCT", slot = "data")[gene, ])

  # Build a single-section ggplot: white spot border, shared fill scale
  make_plot <- function(obj, show_legend = TRUE) {
    p  <- ST.FeaturePlot(obj, features = gene, pt.size = 1.8, ncol = 1,
                         show.sb = FALSE, label.by = "no_label")
    gg <- p[[1]]
    gg$layers[[1]]$aes_params$colour <- "white"  # white border instead of black
    gg <- gg +
      scale_fill_gradientn(colours = c("lightgrey", "red"),
                           limits = c(0, shared_max), na.value = "lightgrey",
                           name = "Expression") +
      theme(plot.title = element_blank())
    if (!show_legend) gg <- gg + theme(legend.position = "none")
    gg
  }

  # ── Individual PDF ──────────────────────────────────────────────────────────
  p_f <- wrap_plots(lapply(flight_objs, make_plot), ncol = 1)
  p_g <- wrap_plots(lapply(ground_objs, make_plot), ncol = 1)

  p_combined <- (p_f | p_g) +
    plot_annotation(
      title = gene,
      caption = "Left: Spaceflight    Right: Ground control",
      theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
                    plot.caption = element_text(size = 12, hjust = 0.5))
    )
  ggsave(paste0(gene, "_spatial_flight_vs_ground.pdf"), p_combined,
         dpi = 300, width = 10, height = 12)
  message("Saved individual: ", gene)

  # ── Panel for combined figure ───────────────────────────────────────────────
  p_f_nl <- wrap_plots(lapply(flight_objs, function(o) make_plot(o, show_legend = FALSE)), ncol = 1)
  p_g_nl <- wrap_plots(lapply(ground_objs, function(o) make_plot(o, show_legend = FALSE)), ncol = 1)

  # Condition sub-labels
  lbl_f  <- ggdraw() + draw_label("Spaceflight",    fontface = "plain", size = 22, hjust = 0.5)
  lbl_g  <- ggdraw() + draw_label("Ground control", fontface = "plain", size = 22, hjust = 0.5)
  lbl_gn <- ggdraw() + draw_label(gene,             fontface = "bold",  size = 32, hjust = 0.5)

  col_f <- plot_grid(lbl_f, p_f_nl, ncol = 1, rel_heights = c(0.05, 0.95))
  col_g <- plot_grid(lbl_g, p_g_nl, ncol = 1, rel_heights = c(0.05, 0.95))

  gene_panel <- plot_grid(
    lbl_gn,
    plot_grid(col_f, col_g, nrow = 1),
    ncol = 1, rel_heights = c(0.04, 0.96)
  )
  all_gene_panels[[gene]] <- gene_panel
}

# ── Helper: build and save a combined figure for a subset of genes ───────────
build_combined <- function(gene_subset, filename) {
  panels <- all_gene_panels[gene_subset]
  n <- length(panels)

  legend_src <- ggplot(data.frame(x = 0, y = 0, z = 0.5), aes(x, y, fill = z)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = c("lightgrey", "red"),
      name    = "Expression",
      limits  = c(0, 1),
      breaks  = c(0, 1),
      labels  = c("Low", "High")
    ) +
    guides(fill = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = unit(10, "cm"), barheight = unit(0.8, "cm")
    )) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 18),
          legend.text  = element_text(size = 16))
  shared_legend <- get_legend(legend_src)

  divider <- ggplot() +
    geom_vline(xintercept = 0.5, color = "gray50", linewidth = 0.6) +
    theme_void()

  panels_with_dividers <- list()
  for (i in seq_along(panels)) {
    panels_with_dividers[[length(panels_with_dividers) + 1]] <- panels[[i]]
    if (i < n)
      panels_with_dividers[[length(panels_with_dividers) + 1]] <- divider
  }
  divider_widths <- c(rep(c(1, 0.08), n - 1), 1)

  combined_panels <- plot_grid(plotlist = panels_with_dividers, nrow = 1,
                               rel_widths = divider_widths)
  final_fig <- plot_grid(combined_panels, shared_legend,
                         ncol = 1, rel_heights = c(0.88, 0.12))

  fig_width <- 14 * n / 2   # scale width proportionally (14 cm per 2 genes)
  ggsave(filename, final_fig, dpi = 300, width = fig_width, height = 14)
  message("Saved combined figure: ", filename)
}

# ── Figure 1: Arpc2 + Arap2 ──────────────────────────────────────────────────
build_combined(c("Arpc2", "Arap2"), "Arpc2_Arap2_spatial_combined.pdf")

# ── Figure 2: Rhob + Sec14l1 + Ubl3 ─────────────────────────────────────────
build_combined(c("Rhob", "Sec14l1", "Ubl3"), "Rhob_Sec14l1_Ubl3_spatial_combined.pdf")

# ── Figure 3: Epha4 + Pink1 + Ly6h ──────────────────────────────────────────
build_combined(c("Epha4", "Pink1", "Ly6h"), "Epha4_Pink1_Ly6h_spatial_combined.pdf")

# ── Figure 4: Pink1 + Snph + Itpr1 ───────────────────────────────────────────
build_combined(c("Pink1", "Snph", "Itpr1"), "Pink1_Snph_Itpr1_spatial_combined.pdf")
