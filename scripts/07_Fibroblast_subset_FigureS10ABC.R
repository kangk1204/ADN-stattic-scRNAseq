###############################################################################
# Script:      07_Fibroblast_subset_FigureS10ABC.R
# Project:     Single-cell RNA-seq analysis of mouse kidney (Sham / ADN /
#              ADN + Stattic)
# Description: Subset analysis of fibroblasts. Re-embedding within the
#              Fibroblast identity yields four sub-populations:
#              Quiescent Fibroblast / Proto-Myofibroblast / Myofibroblast /
#              Pericyte.
#
# Generates the following manuscript panels:
#   Figure S10A -- Fibroblast subset UMAP (4 sub-populations)
#   Figure S10B -- Dot plot of fibroblast / pericyte markers
#   Figure S10C -- scProportionTest (sham vs ADN, ADN vs ADN_Stattic)
#
# Input  : Stattic_annotated.rds  (output of script 02)
# Output : Fibro_subset.rds + the figure panels above.
###############################################################################


# ---------------------------------------------------------------------------
# 1. Libraries
# ---------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(scCustomize)        # DotPlot_scCustom
library(scProportionTest)

set.seed(1234)


# ---------------------------------------------------------------------------
# 2. Configuration
# ---------------------------------------------------------------------------
out_dir  <- "path/to/results/RDS"
fig_dir  <- "path/to/results/Figures"

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

fibro_levels <- c("Quiescent Fibroblast",
                  "Proto-Myofibroblast",
                  "Myofibroblast",
                  "Pericyte")


# ---------------------------------------------------------------------------
# 3. Load the annotated object from script 02
# ---------------------------------------------------------------------------
Stattic <- readRDS(file.path(out_dir, "Stattic_annotated.rds"))


# ---------------------------------------------------------------------------
# 4. Fibroblast subset + re-embedding
# ---------------------------------------------------------------------------
Fibro <- subset(Stattic, idents = "Fibroblast")

Fibro <- ScaleData(Fibro, verbose = FALSE)
Fibro <- RunPCA(Fibro, npcs = 30, verbose = FALSE)
Fibro <- RunUMAP(Fibro, reduction = "harmony", dims = 1:30)
Fibro <- FindNeighbors(Fibro, reduction = "harmony", dims = 1:30)

set.seed(1234)
Fibro <- FindClusters(Fibro, resolution = 0.1)   # yields four clusters: 0-3

Fibro <- RenameIdents(
  Fibro,
  `0` = "Myofibroblast",
  `1` = "Proto-Myofibroblast",
  `2` = "Pericyte",
  `3` = "Quiescent Fibroblast"
)

Fibro$Fibro_CellType <- Idents(Fibro)
Idents(Fibro)         <- "Fibro_CellType"
levels(Fibro)         <- fibro_levels

# Lock the condition factor order so that downstream split.by / group.by
# plots and proportion tests always render in sham -> ADN -> ADN_Stattic.
Fibro$condition <- factor(as.character(Fibro$condition),
                          levels = c("sham", "ADN", "ADN_Stattic"))

saveRDS(Fibro, file.path(out_dir, "Fibro_subset.rds"))


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure S10A  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
figS10A <- DimPlot(Fibro, reduction = "umap", pt.size = 1,
                   label = TRUE, label.size = 5)

ggsave(file.path(fig_dir, "FigS10A_Fibro_UMAP.pdf"),
       figS10A, width = 5, height = 4)
# ===========================================================================


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure S10B  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Dot plot of canonical fibroblast / pericyte markers across the four
# sub-populations.

markers_S10B <- c(
  # Pan-fibroblast / mesenchymal
  "Cd34", "Vim", "Gsn",
  # Quiescent / activation-associated
  "Clu", "Snai1", "Hmgb2",
  # Pro-myofibroblast / TGF-beta axis
  "Pdgfrb", "Tgfbr1", "Tnc", "Thy1",
  # Pericyte / contractile
  "Myh11", "Kcna5", "Cspg4", "Acta2", "Rgs5"
)

figS10B <- DotPlot_scCustom(Fibro, features = markers_S10B,
                            x_lab_rotate = TRUE)

ggsave(file.path(fig_dir, "FigS10B_Fibro_markers.pdf"),
       figS10B, width = 8, height = 3)
# ===========================================================================


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure S10C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Cell-type proportion permutation test for the four fibroblast sub-
# populations.
#   Left  : sham vs ADN          -- effect of disease
#   Right : ADN  vs ADN_Stattic  -- effect of STAT3 inhibition

run_proptest <- function(obj, s1, s2) {
  pt <- sc_utils(obj)
  pt <- permutation_test(pt,
                         cluster_identity = "Fibro_CellType",
                         sample_1         = s1,
                         sample_2         = s2,
                         sample_identity  = "condition")
  permutation_plot(pt)
}

figS10C_left  <- run_proptest(Fibro, "sham", "ADN")
figS10C_right <- run_proptest(Fibro, "ADN",  "ADN_Stattic")

ggsave(file.path(fig_dir, "FigS10C_left_sham_vs_ADN.pdf"),
       figS10C_left,  width = 4, height = 3)
ggsave(file.path(fig_dir, "FigS10C_right_ADN_vs_ADNStattic.pdf"),
       figS10C_right, width = 4, height = 3)
# ===========================================================================


# ---------------------------------------------------------------------------
# 5. Save sessionInfo
# ---------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           con = file.path(out_dir, "sessionInfo_07_Fibroblast.txt"))

###############################################################################
# End of script 07.
# Next: 08_MacroFibro_CellChat_Figure5F_S10DE.R
###############################################################################
