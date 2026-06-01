###############################################################################
# Script:      02_Figure3_FigureS3.R
# Project:     Single-cell RNA-seq analysis of mouse kidney (Sham / ADN /
#              ADN + Stattic)
# Description: Higher-resolution clustering, manual cell-type annotation,
#              Stat3 expression visualisation, gene-set module scores,
#              cell-type proportion testing.
#
# Generates the following manuscript panels:
#   Figure 3A  -- UMAP coloured by cell type (all 25,160 cells)
#   Figure 3B  -- UMAP split by condition
#   Figure 3C  -- Stat3 violin grouped by condition (3 violins)
#   Figure 3D  -- Stat3 violin per cell type, split by condition
#   Figure 3E  -- Module score UMAPs: lipid oxidation / epigenetic regulation /
#                 inflammatory response / ECM organisation
#   Figure 3F  -- Ferroptosis module score UMAP + per cell-type violin
#   Figure S3B -- Average-expression heat-map of canonical cell-type markers
#   Figure S3C -- scProportionTest (sham vs ADN, ADN vs ADN_Stattic)
#   Figure S3D -- Group 1-4 module score UMAPs
#
#   Figure S3A is produced by 01_QC_and_Integration.R and is not regenerated
#   here.
#
# Input  : Kidney_HarmonyIntegrated_dim30_res1.rds  (output of script 01)
#          Gene-set text files derived from `Table S1_Gene List.xlsx`
#          (provided alongside this repository). Each sheet of the workbook
#          must be exported as a tab-delimited .txt file -- one gene symbol
#          per row, with a header line -- and placed in <gene_list_dir>.
#          See the "Helper" section below for the expected file names.
# Output : Cell-type annotated Seurat object + all figure panels listed above
###############################################################################


# ---------------------------------------------------------------------------
# 1. Libraries
# ---------------------------------------------------------------------------
# scRNA-seq analysis.

library(Seurat)
library(scCustomize)        # FeaturePlot_scCustom
library(scProportionTest)   # cell-type proportion permutation test
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(1234)


# ---------------------------------------------------------------------------
# 2. Configuration
# ---------------------------------------------------------------------------
# Paths follow the same convention as script 01.
out_dir       <- "path/to/results/RDS"
fig_dir       <- "path/to/results/Figures"
deg_dir       <- "path/to/results/Documents"
gene_list_dir <- "path/to/results/GeneSets"   # holds the *.txt gene lists

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(deg_dir, showWarnings = FALSE, recursive = TRUE)

# Cell-type colour palette is intentionally not hard-coded here -- Seurat's
# default ggplot palette already matches the manuscript; override only if
# you need exact colour reproduction.

# Final cell-type order used throughout Figure 3 / Figure S3
celltype_levels <- c("PCT", "PST", "LOH", "CNT", "IC", "Podo", "Mesangial",
                     "Endo", "GEC", "Fibroblast", "Myeloid",
                     "T & NK", "B", "Proliferating")

# Condition colours used in Figure 3C
condition_cols  <- c(sham        = "#62C15B",
                     ADN         = "#E66342",
                     ADN_Stattic = "#3F6EFF")


# ---------------------------------------------------------------------------
# 3. Load the integrated object from script 01
# ---------------------------------------------------------------------------
Stattic <- readRDS(file.path(out_dir, "Kidney_HarmonyIntegrated_dim30_res1.rds"))


# ---------------------------------------------------------------------------
# 4. Higher-resolution clustering for finer cell-type discrimination
# ---------------------------------------------------------------------------
# Script 01 saved clusters at resolution = 1.0 (~30 clusters). Annotation of
# closely related kidney cell types (e.g. PCT vs PST, LOH subsets, myeloid
# subsets) requires more granularity, so we re-cluster at resolution = 1.8
# without re-running PCA / UMAP. This yields ~33 clusters (0-32) under our
# environment; small differences in cluster count are expected when the
# Seurat / Harmony / igraph versions differ.

set.seed(1234)
Stattic <- FindClusters(Stattic, resolution = 1.8)

# Sanity check
DimPlot(Stattic, reduction = "umap", label = TRUE, pt.size = 0.01)


# ---------------------------------------------------------------------------
# 5. Cluster marker discovery (used only to build the annotation)
# ---------------------------------------------------------------------------
# Cluster -> cell-type assignment below was made by inspecting the top
# differentially expressed genes per cluster against the canonical kidney
# marker panel (next section). Re-running this on a different machine may
# return slightly different cluster numbers -- always verify the markers
# before reusing the RenameIdents() block downstream.

DefaultAssay(Stattic) <- "RNA"

Stattic.markers <- FindAllMarkers(
  Stattic,
  min.pct          = 0.25,
  only.pos         = TRUE,
  logfc.threshold  = 1
)

write.table(
  Stattic.markers,
  file = file.path(deg_dir, "Stattic_res1.8_clusterMarkers_log2FC1.txt"),
  sep  = "\t", quote = FALSE
)


# ---------------------------------------------------------------------------
# 6. Canonical kidney cell-type marker panels
# ---------------------------------------------------------------------------
# Kept as plain lists for transparency. Used by VlnPlot / DotPlot to verify

markers_PT          <- c("Slc27a2", "Slc34a1", "Slc5a2", "Akr1c21")
markers_InjuredPT   <- c("Lrp2", "Spp2", "Slc22a8")
markers_LOH         <- c("Umod", "Cryab", "Slc12a1")
markers_CNT         <- c("Trpv5", "Calb1", "Kcnj1", "Slc8a1", "Aqp2")
markers_Podo        <- c("Nphs2", "Nphs1", "Clic3")
markers_GEC         <- c("Ehd3")
markers_IC          <- c("Atp6v1g3", "Slc4a1", "Slc26a4")
markers_Endo        <- c("Cd34", "Pecam1", "Kdr")
markers_Fibroblast  <- c("Col1a1", "Col6a3", "Pdgfra", "Dcn")
markers_Mesangial   <- c("Gata3", "Pdgfrb")
markers_Macrophage  <- c("Cd74", "Cd83", "C1qb", "Cd68", "Ccl12", "Clec4a1")
markers_DC          <- c("Cd209a", "Ccl17", "Clec10a", "Flt3", "Xcr1")
markers_Granulocyte <- c("S100a8", "S100a9", "Cxcr2", "Csf3r", "Mmp9")
markers_NKT         <- c("Icos", "Cd4", "Cd28", "Cd3g", "Cd3e", "Themis")
markers_B           <- c("Cd19", "Cd79a", "Fcer2a", "Cd79b")


# ---------------------------------------------------------------------------
# 7. Cell-type annotation
# ---------------------------------------------------------------------------
# Note: macrophage, DC and granulocyte clusters were merged into a single
# "Myeloid" identity for clarity, since further myeloid subsetting is
# performed in a separate downstream script (see Figure 5).

Stattic <- RenameIdents(
  Stattic,
  `0`  = "PST",          `1`  = "Myeloid",      `2`  = "Endo",
  `3`  = "LOH",          `4`  = "Myeloid",      `5`  = "PCT",
  `6`  = "CNT",          `7`  = "Myeloid",      `8`  = "T & NK",
  `9`  = "CNT",          `10` = "T & NK",       `11` = "Myeloid",
  `12` = "PST",          `13` = "GEC",          `14` = "Myeloid",
  `15` = "T & NK",       `16` = "IC",           `17` = "Myeloid",
  `18` = "Endo",         `19` = "GEC",          `20` = "Fibroblast",
  `21` = "Proliferating",`22` = "B",            `23` = "Myeloid",
  `24` = "LOH",          `25` = "LOH",          `26` = "Podo",
  `27` = "IC",           `28` = "Mesangial",    `29` = "Myeloid",
  `30` = "Fibroblast",   `31` = "Fibroblast",   `32` = "Myeloid"
)

# Persist the annotation as a metadata column and lock the display order.
# (The original draft contained two consecutive levels() calls with mismatched
# vectors; the first one would silently fail because it referenced identities
# -- "Mono & Macro", "DC", "Granulocyte" -- that were already merged into
# "Myeloid" by RenameIdents. Only the consolidated order is kept here.)
Stattic$CellType <- Idents(Stattic)
Idents(Stattic)  <- "CellType"
levels(Stattic)  <- celltype_levels

# Lock the condition factor order so that downstream split.by / group.by
# plots always render in sham -> ADN -> ADN_Stattic order, even after
# the object is reloaded by other scripts.
Stattic$condition <- factor(as.character(Stattic$condition),
                            levels = c("sham", "ADN", "ADN_Stattic"))

saveRDS(Stattic, file.path(out_dir, "Stattic_annotated.rds"))


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure 3A & 3B  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# 3A: UMAP coloured by cell type (all 25,160 cells)
fig3A <- DimPlot(Stattic, reduction = "umap", pt.size = 0.01,
                 label = TRUE, label.size = 5)
ggsave(file.path(fig_dir, "Fig3A_UMAP_celltype.pdf"),
       fig3A, width = 7, height = 6)

# 3B: same UMAP split by condition (Sham 9,505 / ADN 9,369 / ADN_Stattic 6,286)
fig3B <- DimPlot(Stattic, reduction = "umap", pt.size = 0.01,
                 label = FALSE, split.by = "condition")
ggsave(file.path(fig_dir, "Fig3B_UMAP_split_condition.pdf"),
       fig3B, width = 15, height = 5)
# ===========================================================================


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure 3C & 3D  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# IMPORTANT: in the original draft the two violin plots were labelled the
# wrong way round (the `split.by` plot was tagged as 3C and the `group.by`
# plot as 3D). Cross-checking with the published figure shows it should be
# the opposite: 3C is the three-condition summary, 3D is the per-cell-type
# split. The labels have been corrected here.

# 3C: Stat3 expression as three condition-level violins (single panel)
fig3C <- VlnPlot(Stattic, features = "Stat3",
                 group.by = "condition", pt.size = 0,
                 cols = condition_cols)
ggsave(file.path(fig_dir, "Fig3C_Stat3_byCondition.pdf"),
       fig3C, width = 4, height = 4)

# 3D: Stat3 expression per cell type, split by condition
fig3D <- VlnPlot(Stattic, features = "Stat3",
                 split.by = "condition", pt.size = 0,
                 cols = condition_cols)
ggsave(file.path(fig_dir, "Fig3D_Stat3_perCellType_splitCondition.pdf"),
       fig3D, width = 10, height = 4)
# ===========================================================================


# ---------------------------------------------------------------------------
# Helper: add a module score from a one-column gene-list text file
# ---------------------------------------------------------------------------
# All gene sets used in this study are provided as a single Excel workbook,
# `Table S1_Gene List.xlsx`, deposited alongside this repository. Each
# sheet of the workbook corresponds to one gene set and contains one
# official mouse gene symbol per row (single column with a header row).
# Before running this script, export each sheet to a tab-delimited .txt
# file and place the .txt files inside <gene_list_dir>, using the file
# names listed below.
#
# Files referenced below:
#   LipidOxidation.txt        (n =  46)  -- Figure 3E
#   EpigeneticRegulation.txt  (n =  48)  -- Figure 3E
#   InflammatoryResponse.txt  (n = 143)  -- Figure 3E
#   ECMOrganization.txt       (n =  76)  -- Figure 3E
#   Ferroptosis.txt           (n =  48)  -- Figure 3F
#   group1.txt ... group4.txt            -- Figure S3D

add_module_score <- function(obj, gene_file, score_name) {
  genes  <- read.table(file.path(gene_list_dir, gene_file),
                       header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
  AddModuleScore(obj, features = list(genes[, 1]),
                 name = paste0(score_name, "_score"))
}


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure 3E  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Four module-score UMAPs.  Note: AddModuleScore() appends a "1" suffix to
# the user-supplied name, so the column to plot is "<name>_score1".

module_sets_3E <- c(
  LipidOxidation       = "LipidOxidation.txt",
  EpigeneticRegulation = "EpigeneticRegulation.txt",
  InflammatoryResponse = "InflammatoryResponse.txt",
  ECMOrganization      = "ECMOrganization.txt"
)

for (nm in names(module_sets_3E)) {
  obj <- add_module_score(Stattic, module_sets_3E[nm], nm)
  p   <- FeaturePlot_scCustom(obj, reduction = "umap",
                              features  = paste0(nm, "_score1"), pt.size = 0.01)
  ggsave(file.path(fig_dir, paste0("Fig3E_", nm, ".pdf")),
         p, width = 5, height = 4)
}
# ===========================================================================


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure 3F  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Left panel : ferroptosis module score on UMAP (all cells)
# Right panel: ferroptosis module score per cell type, restricted to the
#              epithelial and stromal compartments. Immune populations
#              (Myeloid, T & NK, B, Proliferating) are excluded because
#              ferroptosis is a feature of parenchymal injury and including
#              the immune cells would visually compress the y-axis range
#              of the cell types of interest. PST and GEC are highlighted
#              as the ferroptosis-prone populations (dashed outlines in 3F).

Ferroptosis.Module <- add_module_score(Stattic, "Ferroptosis.txt", "Ferroptosis")

# --- 3F left: UMAP ---------------------------------------------------------
fig3F_left <- FeaturePlot_scCustom(
  Ferroptosis.Module, reduction = "umap",
  features  = "Ferroptosis_score1", pt.size = 0.01
)
ggsave(file.path(fig_dir, "Fig3F_Ferroptosis_UMAP.pdf"),
       fig3F_left, width = 5, height = 4)

# --- 3F right: violin (epithelial + stromal only) --------------------------
# Subset to non-immune cell types and lock the display order.
ferro_celltypes_3F <- c("PCT", "PST", "LOH", "CNT", "IC", "Podo",          # epithelial
                        "Mesangial", "Endo", "GEC", "Fibroblast")          # stromal

Ferroptosis.Module.parenchymal <- subset(Ferroptosis.Module,
                                         idents = ferro_celltypes_3F)
Idents(Ferroptosis.Module.parenchymal) <- factor(
  as.character(Idents(Ferroptosis.Module.parenchymal)),
  levels = ferro_celltypes_3F
)

fig3F_right <- VlnPlot(Ferroptosis.Module.parenchymal,
                       features = "Ferroptosis_score1",
                       pt.size  = 0) +
  ylim(-0.4, 1.2) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.1)

ggsave(file.path(fig_dir, "Fig3F_Ferroptosis_VlnPlot.pdf"),
       fig3F_right, width = 7, height = 4)
# ===========================================================================


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure S3B  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Row-scaled average-expression heat-map of canonical kidney markers across
# all annotated cell types.

heatmap_genes <- c(
  # PT markers (PCT and PST share several)
  "Slc9a3", "Slc4a4", "Ppara", "Acox1", "Egfr", "Slc34a1",
  "Slc22a2", "Slc22a8", "Ggt1", "Slc22a6", "Slc7a13",
  # LOH
  "Slc12a1", "Umod",
  # CNT
  "Aqp2", "Slc12a3", "Calb1",
  # IC
  "Atp6v1g3", "Slc26a4", "Atp6v0d2",
  # Podo
  "Clic3", "Nphs2", "Nphs1",
  # Mesangial
  "Sfrp2", "Itga8", "Gata3",
  # Endo
  "Cd93", "Igfbp3", "Robo4",
  # GEC
  "Ehd3", "Tbx3",
  # Fibroblast
  "Col6a3", "Dcn", "Col1a1",
  "Acta2", "Ppp1r14a",
  # Myeloid
  "Ccl12", "Clec4a1", "Cd68",
  "S100a9", "S100a8", "Cxcr2", "Mmp9",
  # T & NK
  "Nkg7", "Icos", "Cd3e", "Cd3g", "Themis", "Il7r",
  # B
  "Cd79a", "Iglc2", "Iglc3", "Pou2af1",
  # Proliferating
  "Mki67", "Hist1h3c"
)

# AverageExpression returns a list keyed by assay; pull the RNA matrix and
# enforce the manual gene order. Genes that are absent from the assay are
# dropped (rather than producing NA rows that would render as grey stripes
# after row-wise scaling in pheatmap).
avg_mat <- AverageExpression(Stattic, features = heatmap_genes,
                             assays = "RNA")$RNA
avg_mat <- as.matrix(avg_mat)

genes_present <- heatmap_genes[heatmap_genes %in% rownames(avg_mat)]
missing_genes <- setdiff(heatmap_genes, genes_present)
if (length(missing_genes) > 0) {
  message("FigS3B: ", length(missing_genes),
          " gene(s) not found and dropped: ",
          paste(missing_genes, collapse = ", "))
}
avg_mat <- avg_mat[genes_present, , drop = FALSE]

pdf(file.path(fig_dir, "FigS3B_marker_heatmap.pdf"),
    width = 6, height = 12)
pheatmap(avg_mat,
         scale          = "row",
         show_rownames  = TRUE,
         cluster_rows   = FALSE,
         cluster_cols   = FALSE,
         color          = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()
# ===========================================================================


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure S3C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Cell-type proportion permutation test.
# Left  : sham vs ADN          -- effect of disease
# Right : ADN  vs ADN_Stattic  -- effect of STAT3 inhibition

run_proptest <- function(obj, s1, s2) {
  pt <- sc_utils(obj)
  pt <- permutation_test(pt,
                         cluster_identity = "CellType",
                         sample_1         = s1,
                         sample_2         = s2,
                         sample_identity  = "condition")
  permutation_plot(pt)
}

figS3C_left  <- run_proptest(Stattic, "sham", "ADN")
figS3C_right <- run_proptest(Stattic, "ADN",  "ADN_Stattic")

ggsave(file.path(fig_dir, "FigS3C_left_sham_vs_ADN.pdf"),
       figS3C_left,  width = 5, height = 5)
ggsave(file.path(fig_dir, "FigS3C_right_ADN_vs_ADNStattic.pdf"),
       figS3C_right, width = 5, height = 5)
# ===========================================================================


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure S3D  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Module-score UMAPs for the four gene groups defined in Figure 2E
# (group1-4 correspond to the four functional gene clusters identified
# from the bulk RNA-seq analysis: glucocorticoid biosynthesis / fatty acid
# metabolism / epigenetic regulation / inflammatory & ECM).

module_sets_S3D <- paste0("group", 1:4, ".txt")
names(module_sets_S3D) <- paste0("group", 1:4)

for (nm in names(module_sets_S3D)) {
  obj <- add_module_score(Stattic, module_sets_S3D[nm], nm)
  p   <- FeaturePlot_scCustom(obj, reduction = "umap",
                              features  = paste0(nm, "_score1"), pt.size = 0.01)
  ggsave(file.path(fig_dir, paste0("FigS3D_", nm, "_score.pdf")),
         p, width = 5, height = 4)
}
# ===========================================================================


# ---------------------------------------------------------------------------
# 8. Save sessionInfo for full reproducibility
# ---------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           con = file.path(out_dir, "sessionInfo_02_Annotation.txt"))

###############################################################################
# End of script 02.
# Next: 03_PST_subset.R  (Figure 4 -- proximal straight tubule subclustering)
###############################################################################
