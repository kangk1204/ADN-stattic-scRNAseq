###############################################################################
# Script:      04_FigureS8.R
# Project:     Single-cell RNA-seq analysis of mouse kidney (Sham / ADN /
#              ADN + Stattic)
# Description: Subset analysis of glomerular endothelial cells (GEC).
#              Re-embedding within GEC, identification of GEC0/GEC1
#              sub-populations, between-condition proportion testing, and
#              characterisation of the two sub-populations by DEG and
#              marker expression.
#
# Generates the following manuscript panels:
#   Figure S8A -- GEC subset UMAP (GEC0 vs GEC1)
#   Figure S8B -- scProportionTest (sham vs ADN, ADN vs ADN_Stattic)
#   Figure S8C -- Volcano plot of GEC0 vs GEC1 DEGs
#                 (DEG TABLE only -- the published volcano was drawn in
#                 GraphPad Prism from the exported table; an optional
#                 EnhancedVolcano block is provided as an alternative.)
#   Figure S8D -- Dot plot of GEC0 / GEC1 marker genes
#
# Input  : Stattic_annotated.rds  (output of script 02)
# Output : GEC subset Seurat object, GEC0-vs-GEC1 DEG table, and the figure
#          panels above.
###############################################################################


# ---------------------------------------------------------------------------
# 1. Libraries
# ---------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(scProportionTest)
library(EnhancedVolcano)    # only required if running the optional
                            # R-based volcano in the Fig S8C section

set.seed(1234)


# ---------------------------------------------------------------------------
# 2. Configuration
# ---------------------------------------------------------------------------
out_dir  <- "path/to/results/RDS"
fig_dir  <- "path/to/results/Figures"
deg_dir  <- "path/to/results/Documents"

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(deg_dir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# 3. Load the annotated object from script 02
# ---------------------------------------------------------------------------
Stattic <- readRDS(file.path(out_dir, "Stattic_annotated.rds"))


# ---------------------------------------------------------------------------
# 4. GEC subset + re-embedding
# ---------------------------------------------------------------------------
# As with the PST subset, re-scaling and re-running PCA refresh the within-
# subset variable structure, while UMAP / FindNeighbors use the parent
# Harmony embedding so that batch correction is preserved.

GEC <- subset(Stattic, idents = "GEC")

GEC <- ScaleData(GEC, verbose = FALSE)
GEC <- RunPCA(GEC, npcs = 30, verbose = FALSE)
GEC <- RunUMAP(GEC, reduction = "harmony", dims = 1:30)
GEC <- FindNeighbors(GEC, reduction = "harmony", dims = 1:30)

set.seed(1234)
GEC <- FindClusters(GEC, resolution = 0.1)   # yields two clusters: 0, 1

GEC               <- RenameIdents(GEC, `0` = "GEC0", `1` = "GEC1")
GEC$GEC_CellType  <- Idents(GEC)
Idents(GEC)       <- "GEC_CellType"
levels(GEC)       <- c("GEC0", "GEC1")

# Lock the condition factor order so that downstream split.by / group.by
# plots and proportion tests always render in sham -> ADN -> ADN_Stattic.
GEC$condition <- factor(as.character(GEC$condition),
                        levels = c("sham", "ADN", "ADN_Stattic"))

saveRDS(GEC, file.path(out_dir, "GEC_subset.rds"))


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure S8A  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
figS8A <- DimPlot(GEC, reduction = "umap", pt.size = 1,
                  label = TRUE, label.size = 5)

ggsave(file.path(fig_dir, "FigS8A_GEC_UMAP.pdf"),
       figS8A, width = 5, height = 4)
# ===========================================================================


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure S8B  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Cell-type proportion permutation test for GEC0 / GEC1.
# Left  : sham vs ADN          -- effect of disease
# Right : ADN  vs ADN_Stattic  -- effect of STAT3 inhibition

run_proptest <- function(obj, s1, s2) {
  pt <- sc_utils(obj)
  pt <- permutation_test(pt,
                         cluster_identity = "GEC_CellType",
                         sample_1         = s1,
                         sample_2         = s2,
                         sample_identity  = "condition")
  permutation_plot(pt)
}

figS8B_left  <- run_proptest(GEC, "sham", "ADN")
figS8B_right <- run_proptest(GEC, "ADN",  "ADN_Stattic")

ggsave(file.path(fig_dir, "FigS8B_left_sham_vs_ADN.pdf"),
       figS8B_left,  width = 4, height = 3)
ggsave(file.path(fig_dir, "FigS8B_right_ADN_vs_ADNStattic.pdf"),
       figS8B_right, width = 4, height = 3)
# ===========================================================================


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure S8C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Volcano plot of GEC0 vs GEC1 DEGs.
#   Positive log2FC -> up in GEC0
#   Negative log2FC -> up in GEC1
#
# WORKFLOW NOTE
# -------------
# In the published manuscript, Figure S8C was produced in GraphPad Prism.
# The DEG table written below was used as the source data: the columns
# `avg_log2FC` and `p_val_adj` were exported and replotted in Prism with
# manual styling (point colours, gene labels, fonts, axis ticks).
#
# This R script therefore generates only the SOURCE DATA TABLE used by
# the figure. An optional `EnhancedVolcano` block is provided at the end
# of this section for users who prefer a fully R-based volcano. The two
# routes share the same point coordinates (avg_log2FC vs -log10(p_adj))
# but the styling will not match the manuscript figure exactly.

DEGs_GEC0_GEC1 <- FindMarkers(GEC,
                              ident.1   = "GEC0",
                              ident.2   = "GEC1",
                              only.pos  = FALSE,
                              min.pct   = 0.25)

write.table(DEGs_GEC0_GEC1,
            file = file.path(deg_dir, "DEGs_GEC0_vs_GEC1.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

# Cutoffs reported in the manuscript: |log2FC| > 1 and adjusted p < 0.05.
# Genes labelled in the manuscript volcano:
#   GEC0-specific (positive log2FC): Spp1, Fabp4, Tmsb4x, Rplp0, Rps2, Rpl3
#   GEC1-specific (negative log2FC): Dlc1, Kdr, Adgrf5, Cdh5, Tie1, Epas1


# --- Optional: R-based volcano (NOT used in the manuscript) ----------------
# Switch the `if (FALSE)` to `if (TRUE)` to render an R volcano. The point
# positions reproduce the Prism figure exactly; cosmetic details (label
# placement via drawConnectors, fonts, tick spacing, colour boundaries)
# will differ.
if (FALSE) {

  fc_cut  <- 1
  p_cut   <- 0.05

  keyvals <- ifelse(
    DEGs_GEC0_GEC1$avg_log2FC >  fc_cut & DEGs_GEC0_GEC1$p_val_adj < p_cut, "red3",
    ifelse(DEGs_GEC0_GEC1$avg_log2FC < -fc_cut & DEGs_GEC0_GEC1$p_val_adj < p_cut,
           "turquoise", "gray")
  )
  names(keyvals)[keyvals == "red3"]      <- "GEC0-specific"
  names(keyvals)[keyvals == "turquoise"] <- "GEC1-specific"
  names(keyvals)[keyvals == "gray"]      <- "ns"

  labelled_genes <- c(
    "Spp1", "Fabp4", "Tmsb4x", "Rplp0", "Rps2", "Rpl3",
    "Dlc1", "Kdr", "Adgrf5", "Cdh5", "Tie1", "Epas1"
  )

  figS8C_alt <- EnhancedVolcano(
    DEGs_GEC0_GEC1,
    lab            = rownames(DEGs_GEC0_GEC1),
    x              = "avg_log2FC",
    y              = "p_val_adj",
    selectLab      = labelled_genes,
    colCustom      = keyvals,
    pCutoff        = p_cut,
    FCcutoff       = fc_cut,
    pointSize      = 1,
    labSize        = 4.0,
    labCol         = "black",
    drawConnectors = TRUE,
    legendLabSize  = 10
  ) + ylim(0, 100) + xlim(-8, 8)

  ggsave(file.path(fig_dir, "FigS8C_alternative_R_volcano.pdf"),
         figS8C_alt, width = 7, height = 6)

}
# ===========================================================================


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure S8D  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Dot plot of canonical GEC0 / GEC1 markers across the two sub-populations.

markers_S8D <- c(
  # GEC0
  "Spp1", "Fabp4", "Tmsb4x", "Rplp0", "Rps2", "Rpl3",
  # GEC1
  "Dlc1", "Kdr", "Adgrf5", "Tie1", "Epas1"
)

# Make sure the y-axis order is GEC0 -> GEC1
levels(GEC) <- c("GEC0", "GEC1")

figS8D <- DotPlot(GEC, features = markers_S8D)

ggsave(file.path(fig_dir, "FigS8D_GEC_markers_dotplot.pdf"),
       figS8D, width = 7, height = 3)
# ===========================================================================


# ---------------------------------------------------------------------------
# 5. Save sessionInfo
# ---------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           con = file.path(out_dir, "sessionInfo_04_GEC.txt"))

###############################################################################
# End of script 04.
###############################################################################
