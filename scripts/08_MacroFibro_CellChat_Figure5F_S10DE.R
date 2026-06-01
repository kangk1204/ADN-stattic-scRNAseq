###############################################################################
# Script:      08_MacroFibro_CellChat_Figure5F_S10DE.R
# Project:     Single-cell RNA-seq analysis of mouse kidney (Sham / ADN /
#              ADN + Stattic)
# Description: Two complementary analyses centred on Macrophage <-> Fibroblast
#              communication:
#
#  Part 1 -- Macrophage + Fibroblast CellChat (per condition)
#            > Figure 5F  : THBS signalling role network
#            > Figure S10D: THBS source-restricted circle plot
#                           (Arg1+ Macrophage as the source)
#
#  Part 2 -- Combined PT + GEC + Macro + Fibro Seurat object for cross-
#            compartment expression visualisation
#            > Figure S10E: Spp1 / Cd44 / Thbs1 violins, split by condition
#
# Input  : Macro_subset.rds        (output of script 05)
#          Fibro_subset.rds        (output of script 07)
#          PST_subset.rds          (output of script 03)
#          GEC_subset.rds          (output of script 04)
# Output : Per-condition Macro+Fibro CellChat objects, combined four-subset
#          Seurat object, and the figure panels above.
###############################################################################


# ---------------------------------------------------------------------------
# 1. Libraries
# ---------------------------------------------------------------------------
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(scCustomize)
library(dplyr)

options(stringsAsFactors = FALSE)
set.seed(1234)


# ---------------------------------------------------------------------------
# 2. Configuration
# ---------------------------------------------------------------------------
out_dir  <- "path/to/results/RDS"
fig_dir  <- "path/to/results/Figures"

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

conditions     <- c("sham", "ADN", "ADN_Stattic")

condition_cols <- c(sham        = "#62C15B",
                    ADN         = "#E66342",
                    ADN_Stattic = "#3F6EFF")

# Cell-type display order for the Macro + Fibro CellChat objects
mf_levels <- c("Monocyte",
               "Arg1+ Macrophage",  "Ccl4+ Macrophage",
               "Dlc+ Macrophage",   "C1qb+ Macrophage", "Cxcl10+ Macrophage",
               "Quiescent Fibroblast", "Proto-Myofibroblast",
               "Myofibroblast",        "Pericyte")

# Cell-type display order for the four-subset combined object (Part 2)
all_levels <- c("PCT", "PST",
                "GEC0", "GEC1",
                "Monocyte",
                "Arg1+ Macrophage",  "Ccl4+ Macrophage",
                "Dlc+ Macrophage",   "C1qb+ Macrophage", "Cxcl10+ Macrophage",
                "Quiescent Fibroblast", "Proto-Myofibroblast",
                "Myofibroblast",        "Pericyte")


# ===========================================================================
# PART 1 - Macro + Fibro CellChat (Figure 5F, Figure S10D)
# ===========================================================================

# ---------------------------------------------------------------------------
# 3. Load Macro and Fibro objects + merge + re-integrate
# ---------------------------------------------------------------------------
Macro <- readRDS(file.path(out_dir, "Macro_subset.rds"))
Fibro <- readRDS(file.path(out_dir, "Fibro_subset.rds"))

mf_merged <- merge(
  x            = Macro,
  y            = list(Fibro),
  add.cell.ids = c("Macro", "Fibro"),
  project      = "Macro_Fibro_CellChat"
)

mf_merged <- NormalizeData(mf_merged)
mf_merged <- FindVariableFeatures(mf_merged)
mf_merged <- ScaleData(mf_merged)
mf_merged <- RunPCA(mf_merged)

mf_merged <- IntegrateLayers(
  object         = mf_merged,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  verbose        = FALSE
)
mf_merged[["RNA"]] <- JoinLayers(mf_merged[["RNA"]])

set.seed(1234)
mf_merged <- FindNeighbors(mf_merged, reduction = "harmony", dims = 1:30)
mf_merged <- FindClusters(mf_merged, resolution = 1)
mf_merged <- RunUMAP(mf_merged, reduction = "harmony", dims = 1:30)


# ---------------------------------------------------------------------------
# 4. Build a single `final_CellType` column from the two subset labels
# ---------------------------------------------------------------------------
Mac_ct <- as.character(mf_merged$Macro_CellType)
Fib_ct <- as.character(mf_merged$Fibro_CellType)
mf_merged$final_CellType <- ifelse(!is.na(Mac_ct), Mac_ct, Fib_ct)

Idents(mf_merged) <- "final_CellType"
levels(mf_merged) <- mf_levels
mf_merged$condition <- factor(mf_merged$condition, levels = conditions)

saveRDS(mf_merged, file.path(out_dir, "Macro_Fibro_combined.rds"))


# ---------------------------------------------------------------------------
# 5. Split by condition + create CellChat objects
# ---------------------------------------------------------------------------
seu_by_cond <- lapply(conditions, function(c) {
  obj <- subset(mf_merged, condition == c)
  Idents(obj)         <- "final_CellType"
  levels(obj)         <- mf_levels
  obj$final_CellType  <- Idents(obj)
  obj
})
names(seu_by_cond) <- conditions

cc_by_cond <- lapply(seu_by_cond, function(obj) {
  createCellChat(object   = obj,
                 meta     = obj@meta.data,
                 group.by = "final_CellType",
                 assay    = "RNA")
})


# ---------------------------------------------------------------------------
# 6. Run CellChat per condition (computeCommunProb is the slow step)
# ---------------------------------------------------------------------------
CellChatDB <- CellChatDB.mouse

run_cellchat <- function(cc, ckpt) {
  cc@DB <- CellChatDB
  cc <- subsetData(cc)

  future::plan("multisession", workers = 4)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- smoothData(cc, adj = PPI.mouse)

  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = 1)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")

  saveRDS(cc, file.path(out_dir,
                        paste0("cellchat_MF_", ckpt, "_final.rds")))
  cc
}

cc_by_cond <- lapply(conditions, function(c) run_cellchat(cc_by_cond[[c]], c))
names(cc_by_cond) <- conditions

# To re-load already-computed objects without re-running computeCommunProb:
#   cc_by_cond <- lapply(conditions, function(c)
#     readRDS(file.path(out_dir, paste0("cellchat_MF_", c, "_final.rds"))))
#   names(cc_by_cond) <- conditions


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure 5F  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# THBS pathway signalling-role network (Sender / Receiver / Mediator /
# Influencer x cell type), one panel per condition. CellChat refers to
# the THBS1-containing pathway as "THBS".

for (cnd in conditions) {
  pdf(file.path(fig_dir, paste0("Fig5F_THBS_roleNetwork_", cnd, ".pdf")),
      width = 8, height = 2.5)
  netAnalysis_signalingRole_network(
    cc_by_cond[[cnd]],
    signaling = "THBS",
    width     = 8,
    height    = 2.5,
    font.size = 10
  )
  dev.off()
}
# ===========================================================================


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure S10D  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# THBS signalling restricted to Arg1+ Macrophage as the source.
# The index is looked up DYNAMICALLY per condition (CellChat may drop
# cell types whose count is zero, which would shift hard-coded indices).

for (cnd in conditions) {
  cc  <- cc_by_cond[[cnd]]
  idx <- match("Arg1+ Macrophage", levels(cc@idents))
  if (is.na(idx)) {
    warning(sprintf("'Arg1+ Macrophage' not present in %s - skipped", cnd))
    next
  }

  pdf(file.path(fig_dir, paste0("FigS10D_THBS_Arg1plus_", cnd, ".pdf")),
      width = 6, height = 6)
  netVisual_aggregate(
    cc,
    signaling        = "THBS",
    sources.use      = idx,
    layout           = "circle",
    vertex.label.cex = 0.8,
    edge.label.cex   = 0.8,
    label.edge       = FALSE,
    arrow.width      = 2,
    arrow.size       = 1
  )
  dev.off()
}
# ===========================================================================


# ===========================================================================
# PART 2 - Combined PT + GEC + Macro + Fibro object (Figure S10E)
# ===========================================================================

# ---------------------------------------------------------------------------
# 7. Load the four subset objects + merge + re-integrate
# ---------------------------------------------------------------------------
PT  <- readRDS(file.path(out_dir, "PST_subset.rds"))
GEC <- readRDS(file.path(out_dir, "GEC_subset.rds"))
# Macro and Fibro are already loaded above

ferro_all <- merge(
  x            = PT,
  y            = list(GEC, Macro, Fibro),
  add.cell.ids = c("PT", "GEC", "Macro", "Fibro"),
  project      = "Ferroptotic_Macro_Fibro"
)

ferro_all <- NormalizeData(ferro_all)
ferro_all <- FindVariableFeatures(ferro_all)
ferro_all <- ScaleData(ferro_all)
ferro_all <- RunPCA(ferro_all)

ferro_all <- IntegrateLayers(
  object         = ferro_all,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  verbose        = FALSE
)
ferro_all[["RNA"]] <- JoinLayers(ferro_all[["RNA"]])


# ---------------------------------------------------------------------------
# 8. Build a single `final_CellType` column from the four subset labels
# ---------------------------------------------------------------------------
PT_ct  <- as.character(ferro_all$PST_CellType)
GEC_ct <- as.character(ferro_all$GEC_CellType)
Mac_ct <- as.character(ferro_all$Macro_CellType)
Fib_ct <- as.character(ferro_all$Fibro_CellType)

ferro_all$final_CellType <-
  ifelse(!is.na(PT_ct),  PT_ct,
  ifelse(!is.na(GEC_ct), GEC_ct,
  ifelse(!is.na(Mac_ct), Mac_ct, Fib_ct)))

Idents(ferro_all)         <- "final_CellType"
levels(ferro_all)         <- all_levels
ferro_all$condition       <- factor(ferro_all$condition, levels = conditions)

saveRDS(ferro_all,
        file.path(out_dir, "Ferroptotic_Macro_Fibro.rds"))


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure S10E  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# Three-feature violin plot of the SPP1 / CD44 / THBS1 axis across all 14
# parenchymal + stromal + myeloid populations, with cells split by
# condition.

figS10E <- VlnPlot(
  ferro_all,
  features = c("Spp1", "Cd44", "Thbs1"),
  cols     = condition_cols,
  split.by = "condition",
  pt.size  = 0
)

ggsave(file.path(fig_dir, "FigS10E_Spp1_Cd44_Thbs1.pdf"),
       figS10E, width = 24, height = 4)
# ===========================================================================


# ---------------------------------------------------------------------------
# 9. Save sessionInfo
# ---------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           con = file.path(out_dir, "sessionInfo_08_MacroFibro_CellChat.txt"))

###############################################################################
# End of script 08.
###############################################################################
