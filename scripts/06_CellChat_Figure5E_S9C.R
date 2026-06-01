###############################################################################
# Script:      06_CellChat_Figure5E_S9C.R
# Project:     Single-cell RNA-seq analysis of mouse kidney (Sham / ADN /
#              ADN + Stattic)
# Description: Cell-cell communication inference with CellChat across the
#              ferroptosis-relevant compartments. PT (PCT/PST), GEC
#              (GEC0/GEC1) and Mono/Macrophage subtypes are merged into
#              a single object, re-integrated, and CellChat is run
#              independently for each condition (sham / ADN / ADN_Stattic).
#
# Generates the following manuscript panels:
#   Figure 5E  -- SPP1 source-restricted circle plots
#                 (PST -> Mac and GEC0 -> Mac, in three conditions)
#   Figure S9C -- Outgoing vs Incoming interaction-strength scatter,
#                 in three conditions
#
# Note   : Figure 5F (THBS1, includes fibroblast subtypes) is not produced
#          by this script. It requires a separate CellChat run on the
#          Macro + Fibroblast subset.
#
# Input  : PST_subset.rds, GEC_subset.rds, Macro_subset.rds
#          (outputs of scripts 03, 04, 05)
# Output : Per-condition CellChat objects + the figure panels above.
###############################################################################


# ---------------------------------------------------------------------------
# 1. Libraries
# ---------------------------------------------------------------------------
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(circlize)
library(dplyr)

options(stringsAsFactors = FALSE)
set.seed(1234)


# ---------------------------------------------------------------------------
# 2. Configuration
# ---------------------------------------------------------------------------
out_dir  <- "path/to/results/RDS"
fig_dir  <- "path/to/results/Figures"

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Final cell-type display order (used as factor levels throughout the script)
final_levels <- c("PCT", "PST",
                  "GEC0", "GEC1",
                  "Monocyte", "Arg1+ Macrophage", "Ccl4+ Macrophage",
                  "Dlc+ Macrophage", "C1qb+ Macrophage", "Cxcl10+ Macrophage")

conditions   <- c("sham", "ADN", "ADN_Stattic")


# ---------------------------------------------------------------------------
# 3. Load the three subset Seurat objects
# ---------------------------------------------------------------------------
PT    <- readRDS(file.path(out_dir, "PST_subset.rds"))
GEC   <- readRDS(file.path(out_dir, "GEC_subset.rds"))
Macro <- readRDS(file.path(out_dir, "Macro_subset.rds"))


# ---------------------------------------------------------------------------
# 4. Merge + re-integrate
# ---------------------------------------------------------------------------
# After subsetting and re-scaling each compartment separately, the three
# objects no longer share a common normalisation / embedding space. We
# merge them, run the standard Seurat workflow, and re-integrate with
# Harmony so that CellChat operates on a comparable expression space.

ferro_merged <- merge(
  x            = PT,
  y            = list(GEC, Macro),
  add.cell.ids = c("PT", "GEC", "Macro"),
  project      = "Ferroptotic_Macro_CellChat"
)

ferro_merged <- NormalizeData(ferro_merged)
ferro_merged <- FindVariableFeatures(ferro_merged)
ferro_merged <- ScaleData(ferro_merged)
ferro_merged <- RunPCA(ferro_merged)

ferro_merged <- IntegrateLayers(
  object         = ferro_merged,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  verbose        = FALSE
)
ferro_merged[["RNA"]] <- JoinLayers(ferro_merged[["RNA"]])

set.seed(1234)
ferro_merged <- FindNeighbors(ferro_merged, reduction = "harmony", dims = 1:30)
ferro_merged <- FindClusters(ferro_merged, resolution = 1)
ferro_merged <- RunUMAP(ferro_merged, reduction = "harmony", dims = 1:30)


# ---------------------------------------------------------------------------
# 5. Build a single `final_CellType` column from the three subset labels
# ---------------------------------------------------------------------------
# Each cell carries a non-NA value in exactly one of PST_CellType,
# GEC_CellType, Macro_CellType (the others are NA after merge). We coalesce
# them into a single column.

PT_ct  <- as.character(ferro_merged$PST_CellType)
GEC_ct <- as.character(ferro_merged$GEC_CellType)
Mac_ct <- as.character(ferro_merged$Macro_CellType)

ferro_merged$final_CellType <- ifelse(!is.na(PT_ct),  PT_ct,
                              ifelse(!is.na(GEC_ct), GEC_ct, Mac_ct))

Idents(ferro_merged) <- "final_CellType"
levels(ferro_merged) <- final_levels

ferro_merged$condition <- factor(ferro_merged$condition, levels = conditions)

saveRDS(ferro_merged,
        file.path(out_dir, "Ferroptotic_Macro_combined.rds"))


# ---------------------------------------------------------------------------
# 6. Split by condition + create CellChat objects
# ---------------------------------------------------------------------------
seu_by_cond <- lapply(conditions, function(c) {
  obj <- subset(ferro_merged, condition == c)
  Idents(obj) <- "final_CellType"
  levels(obj) <- final_levels
  obj$final_CellType <- Idents(obj)
  obj
})
names(seu_by_cond) <- conditions

cc_by_cond <- lapply(seu_by_cond, function(obj) {
  createCellChat(object = obj,
                 meta     = obj@meta.data,
                 group.by = "final_CellType",
                 assay    = "RNA")
})


# ---------------------------------------------------------------------------
# 7. Run the CellChat inference pipeline per condition
# ---------------------------------------------------------------------------
# computeCommunProb is the slow step (typically 10-30 min per condition
# with 4 workers); intermediate checkpoints are saved.

CellChatDB <- CellChatDB.mouse

run_cellchat <- function(cc, ckpt_name) {
  cc@DB <- CellChatDB
  cc <- subsetData(cc)

  future::plan("multisession", workers = 4)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- smoothData(cc, adj = PPI.mouse)

  saveRDS(cc, file.path(out_dir,
                        paste0("cellchat_", ckpt_name, "_pre_communProb.rds")))

  cc <- computeCommunProb(cc, type = "triMean")

  saveRDS(cc, file.path(out_dir,
                        paste0("cellchat_", ckpt_name, "_post_communProb.rds")))

  cc <- filterCommunication(cc, min.cells = 1)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")

  saveRDS(cc, file.path(out_dir,
                        paste0("cellchat_", ckpt_name, "_final.rds")))

  cc
}

cc_by_cond <- lapply(conditions, function(c) {
  run_cellchat(cc_by_cond[[c]], ckpt_name = c)
})
names(cc_by_cond) <- conditions

# To re-load already-computed objects without re-running the slow step:
#   cc_by_cond <- lapply(conditions, function(c)
#     readRDS(file.path(out_dir, paste0("cellchat_", c, "_final.rds"))))
#   names(cc_by_cond) <- conditions


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure 5E  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# SPP1 source-restricted circle plots.
#   Sources : PST and GEC0
#   Targets : all populations
#   Layout  : 2 sources x 3 conditions = 6 panels
#
# IMPORTANT: source indices are looked up DYNAMICALLY per condition.
# CellChat re-orders / drops cell types whose count is zero in a given
# condition, so a hard-coded `sources.use = 2` would silently point at
# the wrong cell type (e.g. GEC0 rather than PST) when one of the early
# levels is missing. `match()` against `levels(cc@idents)` is robust to
# this.
#
# Two versions are saved per panel:
#   - "labelled"  : edge weights printed on the circle (label.edge = TRUE)
#   - "unlabelled": clean version for figure assembly (label.edge = FALSE)

source_names <- c("PST", "GEC0")

plot_spp1_source <- function(cc, src_name, label_edge = TRUE) {
  idx <- match(src_name, levels(cc@idents))
  if (is.na(idx)) {
    warning(sprintf("'%s' not present in this CellChat object - skipped",
                    src_name))
    return(invisible(NULL))
  }
  netVisual_aggregate(
    cc,
    signaling        = "SPP1",
    sources.use      = idx,
    layout           = "circle",
    vertex.label.cex = if (label_edge) 0.5 else 0.1,
    edge.label.cex   = 0.8,
    label.edge       = label_edge,
    arrow.width      = 2,
    arrow.size       = 1
  )
}

for (cnd in conditions) {
  for (src_name in source_names) {
    for (lab in c(TRUE, FALSE)) {
      tag <- if (lab) "labelled" else "unlabelled"
      pdf(file.path(fig_dir,
                    paste0("Fig5E_SPP1_", src_name, "_", cnd, "_", tag, ".pdf")),
          width = 6, height = 6)
      plot_spp1_source(cc_by_cond[[cnd]],
                       src_name   = src_name,
                       label_edge = lab)
      dev.off()
    }
  }
}
# ===========================================================================


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure S9C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Outgoing vs incoming signalling strength scatter plot, per condition.

plot_signaling_scatter <- function(cc, title) {
  netAnalysis_signalingRole_scatter(cc, label.size = 5) +
    xlim(0, 12.5) + ylim(0, 13) +
    labs(title = title) +
    theme(legend.position = "right")
}

figS9C <- wrap_plots(
  plot_signaling_scatter(cc_by_cond$sham,        "sham"),
  plot_signaling_scatter(cc_by_cond$ADN,         "ADN"),
  plot_signaling_scatter(cc_by_cond$ADN_Stattic, "ADN_Stattic"),
  nrow = 1
)

ggsave(file.path(fig_dir, "FigS9C_signaling_role_scatter.pdf"),
       figS9C, width = 15, height = 5)
# ===========================================================================


# ===========================================================================
# OPTIONAL: exploratory CellChat plots NOT included in the published figures
# ---------------------------------------------------------------------------
# The following blocks were used during analysis to inspect the global
# communication network and to generate diagnostic plots. They do not
# correspond to any published panel and can be skipped when reproducing
# Figure 5 / Figure S9. Kept here for completeness / sanity-checking.
# ===========================================================================

if (FALSE) {

  ## (a) Aggregated interaction networks per condition (counts and weights)
  for (cnd in conditions) {
    cc <- cc_by_cond[[cnd]]
    grp <- as.numeric(table(cc@idents))
    par(mfrow = c(1, 2), xpd = TRUE)
    netVisual_circle(cc@net$count,  vertex.weight = grp, weight.scale = TRUE,
                     label.edge = FALSE,
                     title.name = paste0(cnd, " - number of interactions"))
    netVisual_circle(cc@net$weight, vertex.weight = grp, weight.scale = TRUE,
                     label.edge = FALSE,
                     title.name = paste0(cnd, " - interaction strength"))
  }

  ## (b) Heat-map of all-vs-all interaction strengths
  for (cnd in conditions) {
    print(netVisual_heatmap(cc_by_cond[[cnd]], measure = "weight"))
  }

  ## (c) SPP1 pathway heat-map per condition
  for (cnd in conditions) {
    print(netVisual_heatmap(cc_by_cond[[cnd]], signaling = "SPP1",
                            color.heatmap = "Greys"))
  }

  ## (d) SPP1 sender/receiver/mediator/influencer profile per condition
  for (cnd in conditions) {
    netAnalysis_signalingRole_network(
      cc_by_cond[[cnd]],
      signaling = "SPP1",
      width     = 8, height = 2.5,
      font.size = 10
    )
  }

}


# ---------------------------------------------------------------------------
# 8. Save sessionInfo
# ---------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           con = file.path(out_dir, "sessionInfo_06_CellChat.txt"))

###############################################################################
# End of script 06.
###############################################################################
