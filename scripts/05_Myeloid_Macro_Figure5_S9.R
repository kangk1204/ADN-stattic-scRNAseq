###############################################################################
# Script:      05_Myeloid_Macro_Figure5_S9.R
# Project:     Single-cell RNA-seq analysis of mouse kidney (Sham / ADN /
#              ADN + Stattic)
# Description: Subset analysis of myeloid cells (Monocyte / Macrophage / DC /
#              Granulocyte) followed by finer sub-clustering of the
#              Monocyte-Macrophage compartment into six populations.
#              Group-4 module score localisation, marker characterisation,
#              and per-condition module-score testing.
#
# Generates the following manuscript panels:
#   Figure 5A  -- Group-4 module score in myeloid cells
#                 (UMAP + violin per myeloid subtype)
#   Figure 5B  -- UMAP of the six Mono/Macrophage sub-populations
#   Figure 5C  -- Average-expression heat-map of subtype-defining markers
#   Figure 5D  -- Module scores (Inflammatory response / ECM organization /
#                 IL-6/JAK/STAT3 signalling) per Mono/Macrophage subtype,
#                 split by condition, with Kruskal-Wallis + Dunn's test
#   Figure S9A -- Stacked violin of Monocyte / MDM / KRM markers across the
#                 six Mono/Macrophage sub-populations
#   Figure S9B -- Group-4 (n = 742) module score per Mono/Macrophage subtype
#                 with the same statistical layout as Figure 5D
#
# Note   : Figure 5E, S9C are produced by Script 06 (CellChat).
#          Figure 5F (THBS1, includes fibroblast subtypes) requires a
#          separate Macro + Fibroblast CellChat run and is not generated
#          by either script.
#
# Input  : Stattic_annotated.rds  (output of script 02)
#          Module-score gene-set .txt files in <gene_list_dir>:
#            group4.txt
#            group4_GeneList.txt              (n = 742)
#            group4_InflammatoryResponse.txt  (n = 143)
#            group4_ECM.txt                   (n =  76)
#            HALLMARK_IL6_JAK_STAT3_SIGNALING.txt (n = 85)
#          (sheets of `Table S1_Gene List.xlsx`, exported as tab-delimited
#          .txt with one gene symbol per row and a header line.)
# Output : Myeloid + Macro subset Seurat objects, and the figure panels above.
###############################################################################


# ---------------------------------------------------------------------------
# 1. Libraries
# ---------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scCustomize)    # FeaturePlot_scCustom, Stacked_VlnPlot
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
gene_list_dir <- "path/to/results/GeneSets"

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

condition_cols <- c(sham        = "#62C15B",
                    ADN         = "#E66342",
                    ADN_Stattic = "#3F6EFF")

# Mono/Macrophage sub-population display order (used in Fig 5B/C/D, S9A/B)
macro_levels <- c("Monocyte",
                  "Arg1+ Macrophage",
                  "Ccl4+ Macrophage",
                  "Dlc+ Macrophage",
                  "C1qb+ Macrophage",
                  "Cxcl10+ Macrophage")


# ---------------------------------------------------------------------------
# 3. Load the annotated object from script 02
# ---------------------------------------------------------------------------
Stattic <- readRDS(file.path(out_dir, "Stattic_annotated.rds"))


# ---------------------------------------------------------------------------
# 4. Myeloid subset -- first pass clustering and QC trimming
# ---------------------------------------------------------------------------
# An initial sub-clustering of the Myeloid identity is run at resolution 0.6
# only to identify low-quality clusters whose cells show elevated
# percent.mt and reduced nFeature_RNA. These cells (clusters 3 and 7 in
# the first pass) are removed, and the remaining cells are re-embedded.

Myeloid <- subset(Stattic, idents = "Myeloid")

Myeloid <- ScaleData(Myeloid, verbose = FALSE)
Myeloid <- RunPCA(Myeloid, npcs = 30, verbose = FALSE)
Myeloid <- RunUMAP(Myeloid, reduction = "harmony", dims = 1:30)
Myeloid <- FindNeighbors(Myeloid, reduction = "harmony", dims = 1:30)

set.seed(1234)
Myeloid <- FindClusters(Myeloid, resolution = 0.6)

# QC inspection (used to choose the clusters to drop)
VlnPlot(Myeloid, features = c("percent.mt", "nFeature_RNA", "Ptprc"),
        pt.size = 0)

# Drop the low-quality clusters and re-embed
keep_clusters <- as.character(c(0, 1, 2, 4, 5, 6, 8, 9, 10))
Myeloid <- subset(Myeloid, idents = keep_clusters)

Myeloid <- ScaleData(Myeloid, verbose = FALSE)
Myeloid <- RunPCA(Myeloid, npcs = 30, verbose = FALSE)
Myeloid <- RunUMAP(Myeloid, reduction = "harmony", dims = 1:30)
Myeloid <- FindNeighbors(Myeloid, reduction = "harmony", dims = 1:30)

set.seed(1234)
Myeloid <- FindClusters(Myeloid, resolution = 0.6)


# ---------------------------------------------------------------------------
# 5. Myeloid annotation (4 broad classes)
# ---------------------------------------------------------------------------
# Cluster -> identity assignments below were chosen by inspecting per-
# cluster marker genes against the canonical myeloid panel
# (e.g. Vcan/Ly6c2 for Monocyte; Cd68/Cd74/C1q for Macrophage; Vdr, Xcr1,
# Ccr7 for DC subsets; S100a8/9, Csf3r for Granulocyte).

Myeloid <- RenameIdents(
  Myeloid,
  `0` = "Macrophage",   `1` = "Macrophage",   `2` = "Macrophage",
  `3` = "DC",           `4` = "Granulocyte",  `5` = "Monocyte",
  `6` = "DC",           `7` = "Macrophage",   `8` = "DC",
  `9` = "Macrophage",   `10` = "DC"
)

Myeloid$Myeloid_CellType <- Idents(Myeloid)
Idents(Myeloid)           <- "Myeloid_CellType"
levels(Myeloid)           <- c("Monocyte", "Macrophage", "DC", "Granulocyte")

saveRDS(Myeloid, file.path(out_dir, "Myeloid_subset.rds"))


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure 5A  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Group-4 module score localised on the myeloid UMAP, plus a per-subtype
# violin to show that the score peaks in Macrophages.

group4 <- read.table(file.path(gene_list_dir, "group4.txt"),
                     header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)

Myeloid <- AddModuleScore(Myeloid,
                          features = list(group4[, 1]),
                          name     = "group4_score")

fig5A_umap <- FeaturePlot_scCustom(Myeloid, reduction = "umap",
                                   features  = "group4_score1",
                                   na_cutoff = NA, pt.size = 0.01)

fig5A_vln  <- VlnPlot(Myeloid, features = "group4_score1", pt.size = 0)

ggsave(file.path(fig_dir, "Fig5A_group4_UMAP.pdf"),
       fig5A_umap, width = 5, height = 4)
ggsave(file.path(fig_dir, "Fig5A_group4_violin.pdf"),
       fig5A_vln,  width = 5, height = 3)
# ===========================================================================


# ---------------------------------------------------------------------------
# 6. Macro subset (Monocyte + Macrophage) -- finer sub-clustering
# ---------------------------------------------------------------------------
# DC and Granulocyte populations are excluded. Remaining cells are re-
# embedded and clustered at resolution 0.4, yielding six sub-populations.

Macro <- subset(Myeloid, idents = c("Monocyte", "Macrophage"))

Macro <- ScaleData(Macro, verbose = FALSE)
Macro <- RunPCA(Macro, npcs = 30, verbose = FALSE)
Macro <- RunUMAP(Macro, reduction = "harmony", dims = 1:30)
Macro <- FindNeighbors(Macro, reduction = "harmony", dims = 1:30)

set.seed(1234)
Macro <- FindClusters(Macro, resolution = 0.4)

Macro <- RenameIdents(
  Macro,
  `0` = "Ccl4+ Macrophage",
  `1` = "Dlc+ Macrophage",
  `2` = "C1qb+ Macrophage",
  `3` = "Monocyte",
  `4` = "Arg1+ Macrophage",
  `5` = "Cxcl10+ Macrophage"
)

Macro$Macro_CellType <- Idents(Macro)
Idents(Macro)         <- "Macro_CellType"
levels(Macro)         <- macro_levels

# Lock the condition factor order so that downstream split.by / group.by
# plots and faceted module-score panels always render in
# sham -> ADN -> ADN_Stattic.
Macro$condition <- factor(as.character(Macro$condition),
                          levels = c("sham", "ADN", "ADN_Stattic"))

saveRDS(Macro, file.path(out_dir, "Macro_subset.rds"))


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure 5B  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
fig5B <- DimPlot(Macro, reduction = "umap", pt.size = 0.01,
                 label = TRUE, repel = TRUE, label.size = 4)

ggsave(file.path(fig_dir, "Fig5B_Macro_UMAP.pdf"),
       fig5B, width = 5, height = 4)
# ===========================================================================


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure 5C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Row-scaled average-expression heat-map of subtype-defining markers across
# the six Mono/Macrophage sub-populations (3 markers x 6 subtypes).

heatmap_genes_5C <- c(
  # Monocyte
  "Chil3",  "Vcan",   "Ly6c2",
  # Arg1+ Macrophage
  "Arg1",   "Vegfa",  "Fn1",
  # Ccl4+ Macrophage
  "Ccl4",   "Cd83",   "Atf3",
  # Dlc+ Macrophage
  "Dlc1",   "Nrp1",   "Psd3",
  # C1qb+ Macrophage
  "C1qb",   "C1qa",   "C1qc",
  # Cxcl10+ Macrophage
  "Cxcl10", "Ccl5",   "Ifit3b"
)

avg_mat_5C <- AverageExpression(Macro, features = heatmap_genes_5C,
                                assays = "RNA")$RNA
avg_mat_5C <- as.matrix(avg_mat_5C)

# Preserve manual order (only genes actually present in the assay are kept)
genes_present <- heatmap_genes_5C[heatmap_genes_5C %in% rownames(avg_mat_5C)]
avg_mat_5C    <- avg_mat_5C[genes_present, , drop = FALSE]

pdf(file.path(fig_dir, "Fig5C_Macro_marker_heatmap.pdf"),
    width = 5, height = 7)
pheatmap(avg_mat_5C,
         scale         = "row",
         show_rownames = TRUE,
         cluster_rows  = FALSE,
         cluster_cols  = FALSE,
         color         = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
dev.off()
# ===========================================================================


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure S9A  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Stacked violin of Monocyte / Monocyte-derived macrophage (MDM) /
# Kidney-resident macrophage (KRM) markers across the six sub-populations.

markers_S9A <- c(
  # Monocyte markers
  "Chil3", "Vcan",
  # MDM markers
  "Ccr2", "Cd14", "Itgax", "Clec4e",
  # KRM markers
  "Adgre1", "Fcgr1", "Cx3cr1", "Cd81"
)

figS9A <- Stacked_VlnPlot(Macro, features = markers_S9A, pt.size = 0)

ggsave(file.path(fig_dir, "FigS9A_stacked_violin.pdf"),
       figS9A, width = 5, height = 8)
# ===========================================================================


# ---------------------------------------------------------------------------
# 7. Helpers for module-score violin plots (Figure 5D, S9B)
# ---------------------------------------------------------------------------
# Per-cell-type Kruskal-Wallis with Dunn's post-hoc test (Bonferroni-
# adjusted), faceted by sub-population. Returns a single ggplot.

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

plot_module_violin_facet <- function(obj, score_col,
                                     comparisons, y_pos,
                                     ylim_range, title) {
  plot_data <- data.frame(
    expr      = obj@meta.data[[score_col]],
    condition = factor(obj$condition,         levels = names(condition_cols)),
    cell_type = factor(obj$Macro_CellType,    levels = macro_levels)
  )

  # Per cell-type Kruskal-Wallis + Dunn's
  ann_df <- do.call(rbind, lapply(macro_levels, function(ct) {
    sub_data <- plot_data[plot_data$cell_type == ct, ]
    if (nrow(sub_data) < 3 ||
        length(unique(sub_data$condition)) < 2) return(NULL)

    kw <- kruskal.test(expr ~ condition, data = sub_data)
    message(sprintf("[%s | %s] Kruskal-Wallis p = %.3g",
                    title, ct, kw$p.value))

    dunn <- dunnTest(expr ~ condition, data = sub_data,
                     method = "bonferroni")
    ann  <- dunn_annotations(dunn, comparisons)

    data.frame(
      cell_type  = factor(ct, levels = macro_levels),
      group1     = sapply(comparisons, `[`, 1),
      group2     = sapply(comparisons, `[`, 2),
      annotation = ann,
      y_position = y_pos,
      stringsAsFactors = FALSE
    )
  }))

  ggplot(plot_data, aes(x = condition, y = expr, fill = condition)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.1) +
    facet_wrap(~ cell_type, nrow = 1, strip.position = "bottom") +
    scale_fill_manual(values = condition_cols) +
    geom_signif(
      data = ann_df,
      aes(xmin        = group1,
          xmax        = group2,
          annotations = annotation,
          y_position  = y_position),
      manual      = TRUE,
      inherit.aes = FALSE
    ) +
    coord_cartesian(ylim = ylim_range) +
    labs(title = title, y = "module score", x = "") +
    theme_classic() +
    theme(
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      strip.background = element_blank(),
      strip.placement  = "outside",
      strip.text       = element_text(size = 10),
      legend.title     = element_blank()
    )
}


# ---------------------------------------------------------------------------
# 8. Add the four module scores to the Macro object
# ---------------------------------------------------------------------------
add_module <- function(obj, gene_file, score_name) {
  genes <- read.table(file.path(gene_list_dir, gene_file),
                      header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)[, 1]
  AddModuleScore(obj, features = list(genes),
                 name = paste0(score_name, "_score"))
}

Macro <- add_module(Macro, "group4_InflammatoryResponse.txt",  "Inflammatory")
Macro <- add_module(Macro, "group4_ECM.txt",                   "ECM")
Macro <- add_module(Macro, "HALLMARK_IL6_JAK_STAT3_SIGNALING.txt",
                                                               "IL6JAKSTAT3")
Macro <- add_module(Macro, "group4_GeneList.txt",              "Group4")


my_comparisons <- list(c("sham", "ADN"),
                       c("ADN",  "ADN_Stattic"),
                       c("sham", "ADN_Stattic"))


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure 5D  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Three module-score panels per Mono/Macrophage subtype.

# y-axis ranges and bracket positions follow the manuscript figure; adjust
# if your data fall outside the printed ranges.
panels_5D <- list(
  list(score = "Inflammatory_score1",
       title = "Inflammatory response (n = 143)",
       ylim  = c(-0.2, 0.85),
       y_pos = c(0.65, 0.65, 0.78)),
  list(score = "ECM_score1",
       title = "ECM organization (n = 76)",
       ylim  = c(-0.15, 0.45),
       y_pos = c(0.32, 0.32, 0.40)),
  list(score = "IL6JAKSTAT3_score1",
       title = "IL-6/JAK/STAT3 signaling (n = 85)",
       ylim  = c(-0.20, 0.50),
       y_pos = c(0.42, 0.42, 0.47))
)

for (p in panels_5D) {
  g <- plot_module_violin_facet(
    obj         = Macro,
    score_col   = p$score,
    comparisons = my_comparisons,
    y_pos       = p$y_pos,
    ylim_range  = p$ylim,
    title       = p$title
  )
  fname <- gsub("_score1$", "", p$score)
  ggsave(file.path(fig_dir, paste0("Fig5D_", fname, ".pdf")),
         g, width = 12, height = 3.2)
}
# ===========================================================================


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure S9B  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Same layout as Figure 5D, applied to the full group-4 gene list.

figS9B <- plot_module_violin_facet(
  obj         = Macro,
  score_col   = "Group4_score1",
  comparisons = my_comparisons,
  y_pos       = c(0.40, 0.40, 0.46),
  ylim_range  = c(-0.05, 0.50),
  title       = "group 4 (n = 742)"
)

ggsave(file.path(fig_dir, "FigS9B_group4.pdf"),
       figS9B, width = 12, height = 3.2)
# ===========================================================================


# ---------------------------------------------------------------------------
# 9. Save sessionInfo
# ---------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           con = file.path(out_dir, "sessionInfo_05_Myeloid_Macro.txt"))

###############################################################################
# End of script 05.
# Next: 06_CellChat_Figure5E_S9C.R
###############################################################################
