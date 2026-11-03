#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# extract_spatial_data.R
#
# Export per-spot metadata, pixel coordinates and SCT expression for a set of
# genes from the Masarapu et al. RR3 brain Seurat/STutility object WITHOUT
# loading Seurat, SeuratObject or STutility (none are installed).
#
# Route: base-R readRDS() + attr(obj, "<slot>") on the S4 object.  R restores
# S4 objects even if the defining package is not available (the class is
# simply "unregistered"); slots are plain attributes.  Only `Matrix` is needed
# so that the dgCMatrix expression matrices can be subset by row name.
#
# Outputs (written next to this script, in CA3_panels/):
#   spots.csv        one row per spot: barcode, section.name, seurat_clusters,
#                    annot, sample_condition, nCount_SCT, nFeature_SCT,
#                    coordinates from the Staffli meta.data (x, y, adj_x, adj_y,
#                    pixel_x, pixel_y, sample) -- see structure_report.txt for
#                    the exact column names that were present.
#   expr.csv         one row per spot: SCT "data" (log1p corrected counts) for
#                    the genes below as <gene>_data and SCT "counts" as
#                    <gene>_counts (if present).
#   structure_report.txt   what was found in the object (slots, dims, columns).
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(Matrix))

out_dir <- "/Users/justinlin/Desktop/Research /TNT sEV Paper/OSD-352/NASA_RR3_Brain-main/manuscript_plots/CA3_panels"
rds     <- "/Users/justinlin/Desktop/Research /TNT sEV Paper/OSD-352/NASA_RR3_Brain-main/manuscript_plots/RR3-brain_ST_final_object.rds"

genes <- c("Pink1", "Snph", "Itpr1", "Rhob", "Sec14l1", "Ubl3", "Arpc2", "Arap2",
           "Rapgef5", "Nnat", "Adcy1")

report <- file.path(out_dir, "structure_report.txt")
logf <- function(...) { msg <- paste0(..., collapse = ""); cat(msg, "\n"); cat(msg, "\n", file = report, append = TRUE) }
cat("", file = report)

t0 <- Sys.time()
logf("Reading ", rds, " ...")
obj <- readRDS(rds)
logf("readRDS done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
logf("class(obj): ", paste(class(obj), collapse = ","), " ; isS4: ", isS4(obj))
logf("top-level attributes: ", paste(names(attributes(obj)), collapse = ", "))

slot_ <- function(x, name) attr(x, name, exact = TRUE)

# ---- meta.data ------------------------------------------------------------
md <- slot_(obj, "meta.data")
logf("meta.data: ", nrow(md), " spots x ", ncol(md), " columns")
logf("meta.data columns: ", paste(colnames(md), collapse = ", "))
for (cc in c("section.name", "seurat_clusters", "sample_condition", "annot")) {
  if (cc %in% colnames(md)) {
    tb <- table(md[[cc]], useNA = "ifany")
    logf("  table(", cc, "): ", paste(names(tb), tb, sep = "=", collapse = "; "))
  } else logf("  !! column ", cc, " not found")
}

# ---- assays ---------------------------------------------------------------
assays <- slot_(obj, "assays")
logf("assays: ", paste(names(assays), collapse = ", "))
sct <- assays[["SCT"]]
logf("SCT assay attributes: ", paste(names(attributes(sct)), collapse = ", "))
sct_data   <- slot_(sct, "data")
sct_counts <- slot_(sct, "counts")
logf("SCT data: class ", paste(class(sct_data), collapse = ","), " dim ", paste(dim(sct_data), collapse = " x "))
logf("SCT counts: class ", paste(class(sct_counts), collapse = ","), " dim ", paste(dim(sct_counts), collapse = " x "))
logf("SCT data range (first 1000 cols): ", paste(range(sct_data[, 1:min(1000, ncol(sct_data))]), collapse = " - "))
# Check that data really is log1p(counts): compare for a handful of spots
present <- genes[genes %in% rownames(sct_data)]
missing <- setdiff(genes, rownames(sct_data))
logf("genes present in SCT: ", paste(present, collapse = ", "))
if (length(missing)) logf("!! genes missing from SCT: ", paste(missing, collapse = ", "))
chk <- max(abs(as.numeric(sct_data[present, 1:200]) - log1p(as.numeric(sct_counts[present, 1:200]))))
logf("max |data - log1p(counts)| over first 200 spots of selected genes: ", signif(chk, 4),
     "  (0 means data == log1p(SCT counts))")

# ---- Staffli (STutility) coordinates --------------------------------------
tools <- slot_(obj, "tools")
logf("tools: ", paste(names(tools), collapse = ", "))
coords <- NULL
coord_source <- NA
if ("Staffli" %in% names(tools)) {
  st <- tools[["Staffli"]]
  logf("Staffli attributes: ", paste(names(attributes(st)), collapse = ", "))
  smd <- slot_(st, "meta.data")
  logf("Staffli meta.data: ", nrow(smd), " x ", ncol(smd), " ; columns: ", paste(colnames(smd), collapse = ", "))
  logf("Staffli samplenames: ", paste(slot_(st, "samplenames"), collapse = ", "))
  dims <- slot_(st, "dims")
  if (!is.null(dims)) logf("Staffli dims (per sample image size): ", paste(capture.output(str(dims)), collapse = " | "))
  logf("Staffli meta.data head: ", paste(capture.output(print(head(smd, 3))), collapse = " | "))
  logf("Staffli rownames match meta.data rownames: ", identical(rownames(smd), rownames(md)),
       " ; all in: ", all(rownames(md) %in% rownames(smd)))
  coords <- smd[rownames(md), , drop = FALSE]
  coord_source <- "tools$Staffli@meta.data"
}
imgs <- slot_(obj, "images")
logf("images slot names: ", paste(names(imgs), collapse = ", "), " (length ", length(imgs), ")")
if (is.null(coords) && length(imgs)) {
  # Fallback: Seurat VisiumV1 image coordinates
  cl <- lapply(names(imgs), function(n) { cc <- slot_(imgs[[n]], "coordinates"); cc$image <- n; cc })
  coords <- do.call(rbind, cl)
  coords <- coords[rownames(md), , drop = FALSE]
  coord_source <- "images[[*]]@coordinates"
}
logf("coordinate source used: ", coord_source)

# ---- spots.csv ------------------------------------------------------------
keep <- intersect(c("section.name", "seurat_clusters", "annot", "sample_condition",
                    "nCount_SCT", "nFeature_SCT", "nCount_RNA", "nFeature_RNA",
                    "sample_name_cond", "clusters", "orig.ident"), colnames(md))
spots <- data.frame(barcode = rownames(md), md[, keep, drop = FALSE], check.names = FALSE)
if (!is.null(coords)) {
  colnames(coords) <- paste0(colnames(coords))
  spots <- cbind(spots, coords)
}
write.csv(spots, file.path(out_dir, "spots.csv"), row.names = FALSE)
logf("wrote spots.csv: ", nrow(spots), " x ", ncol(spots), " ; columns: ", paste(colnames(spots), collapse = ", "))

# ---- expr.csv -------------------------------------------------------------
ed <- as.matrix(t(sct_data[present, , drop = FALSE]))
colnames(ed) <- paste0(present, "_data")
ec <- as.matrix(t(sct_counts[present, , drop = FALSE]))
colnames(ec) <- paste0(present, "_counts")
stopifnot(identical(rownames(ed), rownames(md)))
expr <- data.frame(barcode = rownames(ed), ed, ec, check.names = FALSE)
write.csv(expr, file.path(out_dir, "expr.csv"), row.names = FALSE)
logf("wrote expr.csv: ", nrow(expr), " x ", ncol(expr))

# ---- also keep a compact RDS of the small pieces so the 1.7 GB object never
#      has to be reloaded for follow-up questions -----------------------------
saveRDS(list(meta.data = md, coords = coords, expr = expr, coord_source = coord_source,
             assay_names = names(assays), sct_genes = rownames(sct_data)),
        file.path(out_dir, "small_pieces.rds"))
logf("Total time: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
