###############################################################################
# Script:      01_QC_and_Integration.R
# Project:     Single-cell RNA-seq analysis of mouse kidney (Sham / ADN /
#              ADN + Stattic)
# Description: Quality control, ambient RNA removal (SoupX), doublet detection
#              (scDblFinder), filtering, merging, and Harmony integration of
#              three 10x Genomics scRNA-seq samples.
#
# Sample design:
#   - sham         : control kidney
#   - ADN          : adenine-induced nephropathy
#   - ADN_Stattic  : adenine-induced nephropathy + Stattic (STAT3 inhibitor)
#
# Input  : Cell Ranger v7.2.0 outputs (filtered_feature_bc_matrix and the full
#          sample folder containing both raw_feature_bc_matrix and
#          filtered_feature_bc_matrix, required by SoupX::load10X)
# Output : An integrated Seurat v5 object ready for cell-type annotation.
#
# Tested with:
#   R 4.4.2, Seurat 5.3.0, SoupX 1.6.2, scDblFinder 1.20.2, harmony 1.2.3
###############################################################################


# ---------------------------------------------------------------------------
# 1. Libraries
# ---------------------------------------------------------------------------
library(Seurat)        # v5.3.0  -- core single-cell analysis
library(SoupX)         # ambient RNA correction
library(scDblFinder)   # doublet detection
library(dplyr)
library(ggplot2)
library(patchwork)
library(scCustomize)   # convenience plotting helpers (optional)

set.seed(1234)         # global seed; additional seeds are set before each
                       # stochastic step (UMAP, clustering, scDblFinder).


# ---------------------------------------------------------------------------
# 2. User configuration
# ---------------------------------------------------------------------------
# Edit only this block to reproduce the analysis on a different machine.
# All downstream code references these variables, so paths are not hard-coded
# anywhere else.

# Root directory containing the three Cell Ranger sample folders
data_root   <- "path/to/cellranger/HN00216204_10X_RawData_Outs"

# Output directory for intermediate and final RDS files
out_dir     <- "path/to/results/RDS"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Sample names must match folder names under data_root
samples     <- c("sham", "ADN", "ADN_Stattic")

# Per-sample paths
sample_dirs <- file.path(data_root, samples, "outs")
filt_dirs   <- file.path(sample_dirs, "filtered_feature_bc_matrix")
names(sample_dirs) <- samples
names(filt_dirs)   <- samples


# ---------------------------------------------------------------------------
# 3. Ambient RNA removal with SoupX
# ---------------------------------------------------------------------------
# SoupX::load10X requires the full Cell Ranger "outs/" directory because it
# reads BOTH raw and filtered matrices to estimate the soup profile.
# We therefore run SoupX BEFORE creating Seurat objects, and use the
# corrected count matrix as the input to CreateSeuratObject.

run_soupx <- function(outs_dir) {
  sc <- load10X(dataDir = outs_dir)
  sc <- autoEstCont(sc)
  adjustCounts(sc)
}

adj_matrices <- lapply(sample_dirs, run_soupx)
names(adj_matrices) <- samples

# Save the SoupX-corrected matrices for traceability
for (s in samples) {
  saveRDS(adj_matrices[[s]],
          file = file.path(out_dir, paste0(s, "_adj_matrix.rds")))
}


# ---------------------------------------------------------------------------
# 4. Build Seurat objects from SoupX-corrected counts
# ---------------------------------------------------------------------------
# NOTE: In the original draft, CreateSeuratObject() was called twice per
# sample -- once on the uncorrected Read10X matrix (with min.cells = 3,
# min.features = 200) and again on the SoupX-corrected matrix (with no
# filters). Because the second call OVERWRITES the first using the same
# variable name, only the second call defines the object that is actually
# used downstream, and the final object therefore has NO min.cells /
# min.features filter applied at object-construction time. Cell-level
# filtering is performed explicitly in section 7b
# (nFeature_RNA > 300, etc.).
#
# We replicate that exact behaviour here so that the published cell counts
# (sham = 9505, ADN = 9369, ADN_Stattic = 6286) are reproduced bit-for-bit.

create_seurat <- function(counts, project) {
  CreateSeuratObject(
    counts  = counts,
    project = project,
    assay   = "RNA"
  )
}

seurat_list <- mapply(
  create_seurat,
  counts  = adj_matrices,
  project = samples,
  SIMPLIFY = FALSE
)


# ---------------------------------------------------------------------------
# 5. QC metrics
# ---------------------------------------------------------------------------
# Mouse gene-symbol conventions are used (lowercase "mt-", capitalised
# "Rps/Rpl", "Hb"). For HUMAN data, switch to "^MT-", "^RP[SL]", "^HB[^P]".
#
# Hemoglobin pattern note:
#   The original draft used "^Hb[^(P)]". Inside a character class, the
#   parentheses are interpreted literally, not as grouping, so this matches
#   "any character that is not '(', ')', or 'P'". This works in practice
#   but is misleading. More importantly, mouse hemoglobin-pseudo gene names
#   such as Hbp1 use a LOWER-case 'p', so the upper-case [^P] would not
#   exclude them. The pattern has been corrected to "^Hb[^p]" for mouse.

add_qc <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^Hb[^p]")
  obj
}

seurat_list <- lapply(seurat_list, add_qc)


# ---------------------------------------------------------------------------
# 6. Doublet detection with scDblFinder
# ---------------------------------------------------------------------------
add_doublets <- function(obj) {
  sce <- Seurat::as.SingleCellExperiment(obj)
  set.seed(1234)
  sce <- scDblFinder(sce)
  obj$scDblFinder.class <- sce$scDblFinder.class
  obj$scDblFinder.score <- sce$scDblFinder.score
  obj
}

seurat_list <- lapply(seurat_list, add_doublets)

# Quick summary
lapply(seurat_list, function(x) table(x$scDblFinder.class))

# Save per-sample objects after QC annotation (before filtering)
for (s in samples) {
  saveRDS(seurat_list[[s]],
          file = file.path(out_dir, paste0(s, "_scDblFinder.rds")))
}


# ---------------------------------------------------------------------------
# 7. Cell filtering
# ---------------------------------------------------------------------------
# Thresholds were chosen by inspecting per-sample violin plots of
# nFeature_RNA, nCount_RNA and percent.mt (see Extended Data Fig. SX).
#
# IMPORTANT: percent.mt cut-off differs between samples. The ADN sample
# showed a tighter mitochondrial distribution and a stricter 5% cut-off
# was sufficient to remove low-quality cells, whereas the sham and
# ADN+Stattic samples required 20% to retain comparable cell numbers
# without admitting clearly stressed cells. This intentional difference
# is documented here for transparency.

qc_thresholds <- list(
  sham        = list(min_feat = 300, max_feat = 6000, max_count = 20000, max_mt = 20),
  ADN         = list(min_feat = 300, max_feat = 6000, max_count = 20000, max_mt =  5),
  ADN_Stattic = list(min_feat = 300, max_feat = 6000, max_count = 20000, max_mt = 20)
)


# ---------------------------------------------------------------------------
# 7a. Pre-filter merge (for the BEFORE-filtering QC violin plots)
# ---------------------------------------------------------------------------
# We build a merged-but-not-yet-filtered object so that the BEFORE-filtering
# violin plots in the manuscript QC figure (left half) are reproduced from
# the same grouping (`condition`) as the AFTER-filtering plots.

for (s in samples) seurat_list[[s]]$condition <- s

Kidney_prefilter <- merge(
  x            = seurat_list[[1]],
  y            = seurat_list[-1],
  add.cell.ids = samples,
  project      = "ADN_Stattic"
)
Kidney_prefilter$condition <- factor(Kidney_prefilter$condition, levels = samples)
Idents(Kidney_prefilter)   <- "condition"

saveRDS(Kidney_prefilter,
        file = file.path(out_dir, "Kidney_prefilter_merged.rds"))


# ---------------------------------------------------------------------------
# 7b. Apply per-sample filters
# ---------------------------------------------------------------------------
filter_cells <- function(obj, th) {
  subset(
    obj,
    subset = nFeature_RNA > th$min_feat &
             nFeature_RNA < th$max_feat &
             nCount_RNA   < th$max_count &
             percent.mt   < th$max_mt   &
             scDblFinder.class == "singlet"
  )
}

seurat_filt <- mapply(
  filter_cells,
  obj = seurat_list,
  th  = qc_thresholds,
  SIMPLIFY = FALSE
)

# Save per-sample filtered objects
for (s in samples) {
  saveRDS(seurat_filt[[s]],
          file = file.path(out_dir, paste0(s, "_filtered.rds")))
}


# ---------------------------------------------------------------------------
# 8. Post-filter merge
# ---------------------------------------------------------------------------
# In the original draft, merge() was called twice in succession, the second
# call silently overwriting the first with different add.cell.ids. Only the
# meaningful (named) merge is kept here.

for (s in samples) seurat_filt[[s]]$condition <- s

Kidney <- merge(
  x            = seurat_filt[[1]],
  y            = seurat_filt[-1],
  add.cell.ids = samples,
  project      = "ADN_Stattic"
)

# Lock the factor order so that figures are always plotted in
# sham -> ADN -> ADN_Stattic order downstream
Kidney$condition <- factor(Kidney$condition, levels = samples)
Idents(Kidney)   <- "condition"

table(Kidney$condition)
# Expected (this dataset):
#   sham         9505
#   ADN          9369
#   ADN_Stattic  6286

saveRDS(Kidney, file = file.path(out_dir, "Kidney_merged.rds"))


# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  Figure S3A  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================
# 8a. QC violin plots -- BEFORE vs AFTER filtering
# ---------------------------------------------------------------------------
# Reproduces the six violin panels of Figure S3A:
#   nFeature_RNA / nCount_RNA / percent.mt, grouped by `condition`
#   (sham, ADN, ADN_Stattic), shown before (left) and after (right) filtering.
# The threshold lines (< 6000, < 20000, < 20 / < 5 / < 20) and the connecting
# arrow in the manuscript figure were added in PowerPoint / Illustrator and
# are intentionally not drawn here.

fig_dir <- file.path(dirname(out_dir), "Figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

# --- Figure S3A LEFT: before filtering -------------------------------------
p_before <- VlnPlot(Kidney_prefilter, features = qc_features,
                    group.by = "condition", ncol = 3, pt.size = 0)

# --- Figure S3A RIGHT: after filtering -------------------------------------
p_after  <- VlnPlot(Kidney,            features = qc_features,
                    group.by = "condition", ncol = 3, pt.size = 0)

ggsave(file.path(fig_dir, "FigS3A_left_before_filter.pdf"),
       plot = p_before, width = 9, height = 3.2)
ggsave(file.path(fig_dir, "FigS3A_right_after_filter.pdf"),
       plot = p_after,  width = 9, height = 3.2)

# Combined panel matching the Figure S3A layout (6 plots in one row)
ggsave(file.path(fig_dir, "FigS3A_combined.pdf"),
       plot  = p_before | p_after,
       width = 18, height = 3.2)
# ===========================================================================
# <<<<<<<<<<<<<<<<<<<<<<<<<<  end Figure S3A  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ===========================================================================


# ---------------------------------------------------------------------------
# 9. Standard pre-integration workflow (Seurat v5 layered object)
# ---------------------------------------------------------------------------
# In Seurat v5 a merged object keeps per-sample count layers. Normalisation,
# variable-feature selection and scaling are layer-aware and work directly
# on the merged object; the layers are joined AFTER integration.

Kidney <- NormalizeData(Kidney)
Kidney <- FindVariableFeatures(Kidney)
Kidney <- ScaleData(Kidney)
Kidney <- RunPCA(Kidney)

ElbowPlot(Kidney, ndims = 30)

# Unintegrated embedding (used only for the before/after integration figure)
set.seed(1234)
Kidney <- FindNeighbors(Kidney, dims = 1:30, reduction = "pca")
Kidney <- FindClusters(Kidney, resolution = 1,
                       cluster.name = "unintegrated_clusters")
Kidney <- RunUMAP(Kidney, dims = 1:30, reduction = "pca",
                  reduction.name = "umap.unintegrated")

DimPlot(Kidney, reduction = "umap.unintegrated",
        group.by = c("condition", "unintegrated_clusters"))


# ---------------------------------------------------------------------------
# 10. Harmony integration
# ---------------------------------------------------------------------------
Kidney <- IntegrateLayers(
  object        = Kidney,
  method        = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose       = FALSE
)

# Re-join layers so that downstream DE / FindMarkers operate on a single
# count matrix.
Kidney[["RNA"]] <- JoinLayers(Kidney[["RNA"]])

set.seed(1234)
Kidney <- FindNeighbors(Kidney, reduction = "harmony", dims = 1:30)
Kidney <- FindClusters(Kidney, resolution = 1)
Kidney <- RunUMAP(Kidney, dims = 1:30, reduction = "harmony",
                  reduction.name = "umap")

# Visual sanity checks
# NOTE: the original draft tried to colour by `seurat_annotations`, but this
# column is only created in the next script (cell-type annotation). It has
# been replaced with `seurat_clusters` here.
DimPlot(Kidney, reduction = "umap",
        group.by = c("condition", "seurat_clusters"))
DimPlot(Kidney, reduction = "umap", split.by = "condition")


# ---------------------------------------------------------------------------
# 11. Save the integrated object for downstream analyses
# ---------------------------------------------------------------------------
saveRDS(Kidney,
        file = file.path(out_dir,
                         "Kidney_HarmonyIntegrated_dim30_res1.rds"))

# sessionInfo() should be saved alongside the RDS for full reproducibility
writeLines(capture.output(sessionInfo()),
           con = file.path(out_dir, "sessionInfo_01_QC.txt"))

###############################################################################
# End of script 01.
# Next: 02_CellType_Annotation.R
###############################################################################
