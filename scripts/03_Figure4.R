###############################################################################
# Script:      03_Figure4.R
# Project:     Single-cell RNA-seq analysis of mouse kidney (Sham / ADN /
#              ADN + Stattic)
# Description: Subset analysis of proximal straight tubule (PST) cells.
#              Re-embedding within PST, between-condition DEG analysis,
#              curated heat-map of condition-defining genes, and per-gene
#              violin plots with Kruskal-Wallis + Dunn's post-hoc test.
#
# Generates the following manuscript panels:
#   Figure 4A  -- PST UMAP coloured by condition (and PST highlight inset)
#   Figure 4B  -- Heat-map of condition-defining gene clusters in PST
#   Figure 4C  -- Per-gene violin plots with Kruskal-Wallis + Dunn's test
#                 (tubular integrity: Acsm2, Lrp2, Slc34a1;
#                  pro-ferroptosis : Spp1, Vdac2, Vdac3)
#
# Input  : Stattic_annotated.rds                       (output of script 02)
#          Gene-set text files derived from `Table S1_Gene List.xlsx`
#          (place exported .txt files inside <gene_list_dir>)
# Output : PST subset Seurat object, DEG tables, and the figure panels above.
###############################################################################


# ---------------------------------------------------------------------------
# 1. Libraries
# ---------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(FSA)            # dunnTest
library(ggsignif)       # geom_signif

set.seed(1234)


# ---------------------------------------------------------------------------
# 2. Configuration
# ---------------------------------------------------------------------------
out_dir       <- "path/to/results/RDS"
fig_dir       <- "path/to/results/Figures"
deg_dir       <- "path/to/results/Documents"
gene_list_dir <- "path/to/results/GeneSets"

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(deg_dir, showWarnings = FALSE, recursive = TRUE)

celltype_levels <- c("PCT", "PST", "LOH", "CNT", "IC", "Podo", "Mesangial",
                     "Endo", "GEC", "Fibroblast", "Myeloid",
                     "T & NK", "B", "Proliferating")

condition_cols  <- c(sham        = "#62C15B",
                     ADN         = "#E66342",
                     ADN_Stattic = "#3F6EFF")


# ---------------------------------------------------------------------------
# 3. Load the annotated object from script 02
# ---------------------------------------------------------------------------
Stattic <- readRDS(file.path(out_dir, "Stattic_annotated.rds"))


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>  Figure 4A (inset)  <<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Highlight PST on the full-kidney UMAP. All other cell types are greyed
# out so that the spatial location of PST is unambiguous.

pst_highlight_cols           <- setNames(rep("grey", length(celltype_levels)),
                                         celltype_levels)
pst_highlight_cols["PST"]    <- "#D78D31"

fig4A_inset <- DimPlot(Stattic, reduction = "umap",
                       cols = pst_highlight_cols)

ggsave(file.path(fig_dir, "Fig4A_inset_PST_highlight.pdf"),
       fig4A_inset, width = 6, height = 5)
# ===========================================================================


# ---------------------------------------------------------------------------
# 4. PST subset + re-embedding
# ---------------------------------------------------------------------------
# Re-scale and re-embed using only PST cells. We reuse the integrated Harmony
# embedding rather than re-running PCA from scratch, so that batch effects
# remain corrected within the subset. A very low clustering resolution is
# used because PST is treated here as a single population; downstream
# differences are tested between conditions, not between sub-clusters.

PST <- subset(Stattic, idents = "PST")

PST <- ScaleData(PST, verbose = FALSE)
PST <- RunPCA(PST, npcs = 30, verbose = FALSE)
PST <- RunUMAP(PST, reduction = "harmony", dims = 1:30)
PST <- FindNeighbors(PST, reduction = "harmony", dims = 1:30)

set.seed(1234)
PST <- FindClusters(PST, resolution = 0.01)   # yields a single cluster ("0")

PST <- RenameIdents(PST, `0` = "PST")
PST$PST_CellType <- Idents(PST)

# Lock the condition factor order so that downstream split.by / group.by
# plots and AverageExpression columns always render in
# sham -> ADN -> ADN_Stattic order.
PST$condition <- factor(as.character(PST$condition),
                        levels = c("sham", "ADN", "ADN_Stattic"))

saveRDS(PST, file.path(out_dir, "PST_subset.rds"))


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure 4A  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# PST UMAP coloured by condition.

Idents(PST) <- "condition"

fig4A <- DimPlot(PST, reduction = "umap", pt.size = 0.1,
                 label = FALSE, group.by = "condition",
                 cols = condition_cols)

ggsave(file.path(fig_dir, "Fig4A_PST_UMAP_by_condition.pdf"),
       fig4A, width = 5, height = 4)
# ===========================================================================


# ---------------------------------------------------------------------------
# 5. Between-condition DEG analysis in PST
# ---------------------------------------------------------------------------
# Output tables are used externally to curate the gene clusters shown in
# Figure 4B. Genes are reported with |log2FC| > 1 and detection in at least
# 25% of cells in either group.

run_DEG <- function(obj, ident.1, ident.2, out_path) {
  Idents(obj) <- "condition"
  degs <- FindMarkers(obj,
                      ident.1         = ident.1,
                      ident.2         = ident.2,
                      only.pos        = FALSE,
                      min.pct         = 0.25,
                      logfc.threshold = 1)
  write.table(degs, file = out_path, sep = "\t",
              quote = FALSE, col.names = NA)
  invisible(degs)
}

run_DEG(PST, "sham", "ADN",
        file.path(deg_dir, "DEGs_sham_vs_ADN_PST_log2FC1.txt"))

run_DEG(PST, "sham", "ADN_Stattic",
        file.path(deg_dir, "DEGs_sham_vs_ADNStattic_PST_log2FC1.txt"))

run_DEG(PST, "ADN",  "ADN_Stattic",
        file.path(deg_dir, "DEGs_ADN_vs_ADNStattic_PST_log2FC1.txt"))


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure 4B  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Row-scaled average-expression heat-map across the three conditions for the
# four DEG-derived gene clusters (35 + 198 + 16 + 16 genes). The
# ADN_Stattic-specific up-regulated cluster was not retained because it
# contained too few genes to support a robust signature. The final, curated
# gene list is provided as a sheet of `Table S1_Gene List.xlsx`; export it
# to <gene_list_dir>/final_PST_Heatmap.txt before running the next block.

draw_avg_heatmap <- function(obj, gene_file, out_pdf,
                             group_by = "condition",
                             width = 4, height = 10) {
  genes <- read.table(file.path(gene_list_dir, gene_file),
                      header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)[, 1]

  # group.by is passed explicitly so the column order is independent of
  # whatever Idents() the caller happens to have set.
  avg_mat <- AverageExpression(obj, features = genes,
                               group.by = group_by,
                               assays   = "RNA")$RNA
  avg_mat <- as.matrix(avg_mat)

  # Keep only genes that exist in the assay, and preserve the input order.
  genes_present <- genes[genes %in% rownames(avg_mat)]
  avg_mat       <- avg_mat[genes_present, , drop = FALSE]

  pdf(out_pdf, width = width, height = height)
  pheatmap(avg_mat,
           scale         = "row",
           show_rownames = TRUE,
           cluster_rows  = FALSE,
           cluster_cols  = FALSE,
           color         = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
  dev.off()
}

draw_avg_heatmap(PST,
                 gene_file = "final_PST_Heatmap.txt",
                 group_by  = "condition",
                 out_pdf   = file.path(fig_dir, "Fig4B_PST_heatmap.pdf"))
# ===========================================================================


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure 4C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Per-gene violin plots grouped by condition with Kruskal-Wallis and
# Dunn's post-hoc test (Bonferroni-adjusted).

# --- helpers --------------------------------------------------------------
signif_label <- function(p) {
  if (is.na(p))         "ns"
  else if (p < 0.0001)  "****"
  else if (p < 0.001)   "***"
  else if (p < 0.01)    "**"
  else if (p < 0.05)    "*"
  else                  "ns"
}

dunn_annotations <- function(dunn_res, comparisons) {
  sapply(comparisons, function(pair) {
    p1 <- paste(pair[1], "-", pair[2])
    p2 <- paste(pair[2], "-", pair[1])
    row <- dunn_res$res[dunn_res$res$Comparison %in% c(p1, p2), ]
    signif_label(row$P.adj)
  })
}

violin_with_dunn <- function(obj, gene,
                             comparisons,
                             y_positions,
                             ylim_max,
                             cols = condition_cols) {
  expr_df <- data.frame(
    expression = FetchData(obj, vars = gene)[, 1],
    condition  = obj$condition
  )

  kw   <- kruskal.test(expression ~ condition, data = expr_df)
  message(sprintf("[%s] Kruskal-Wallis p = %.3g", gene, kw$p.value))

  dunn <- dunnTest(expression ~ condition, data = expr_df,
                   method = "bonferroni")
  print(dunn$res)

  ann  <- dunn_annotations(dunn, comparisons)

  VlnPlot(obj, features = gene,
          cols     = cols,
          group.by = "condition",
          pt.size  = 0) +
    geom_signif(comparisons = comparisons,
                annotations = ann,
                y_position  = y_positions) +
    ylim(0, ylim_max)
}

# --- per-gene plotting ----------------------------------------------------
my_comparisons <- list(c("sham", "ADN"),
                       c("ADN",  "ADN_Stattic"),
                       c("sham", "ADN_Stattic"))

# Per-gene y-axis ranges and bracket positions used in the manuscript figure.
# Edit these if your data fall outside the printed ranges.
genes_4C <- list(
  Acsm2   = list(ylim = 7, y_pos = c(5.5, 5.5, 6.5)),
  Lrp2    = list(ylim = 5, y_pos = c(4.0, 4.0, 4.5)),
  Slc34a1 = list(ylim = 5, y_pos = c(4.0, 4.0, 4.5)),
  Spp1    = list(ylim = 9, y_pos = c(7.0, 7.0, 8.0)),
  Vdac2   = list(ylim = 5, y_pos = c(4.0, 4.0, 4.5)),
  Vdac3   = list(ylim = 5, y_pos = c(4.0, 4.0, 4.5))
)

for (g in names(genes_4C)) {
  cfg <- genes_4C[[g]]
  p   <- violin_with_dunn(PST, g,
                          comparisons = my_comparisons,
                          y_positions = cfg$y_pos,
                          ylim_max    = cfg$ylim)
  ggsave(file.path(fig_dir, paste0("Fig4C_", g, "_violin.pdf")),
         p, width = 4, height = 4)
}
# ===========================================================================


# ---------------------------------------------------------------------------
# 6. Save sessionInfo
# ---------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           con = file.path(out_dir, "sessionInfo_03_PST.txt"))

###############################################################################
# End of script 03.
###############################################################################
