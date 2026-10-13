###############################################################################
# Script:      08_MacroFibro_CellChat_Figure5F_S11DEF.R
# Project:     Mouse kidney scRNA-seq (Sham / ADN / ADN + Stattic)
# Revision:    Current integrated manuscript, 2026-10-12
# Description: Macrophage-fibroblast communication and cross-compartment expression.
#
# Manuscript panels / source data:
#   Figures 5F; S11D-F.
#
# Input:       05 Macro; 07 Fibro; 06 PT CellChat; 04 GEC.
# Design:      One pooled library from two mice per condition (three libraries).
###############################################################################

# ---------------------------------------------------------------------------
# 1. Libraries
# ---------------------------------------------------------------------------
library(Seurat)
library(CellChat)
library(future)
library(ggplot2)
library(patchwork)
library(digest)

# ---------------------------------------------------------------------------
# 2. Configuration
# ---------------------------------------------------------------------------

# Run this R script with Rscript --vanilla or source() in a fresh R session.
locate_project <- function() {
  for(i in rev(seq_along(sys.frames()))) {
    f<-sys.frame(i)$ofile
    if(!is.null(f))return(dirname(normalizePath(f)))
  }
  arg<-grep("^--file=",commandArgs(),value=TRUE)
  if(length(arg)) {
    f<-sub("^--file=","",arg[1])
    # Some Rscript launchers encode spaces in --file as ~+~.
    if(!file.exists(f))f<-gsub("~+~"," ",f,fixed=TRUE)
    return(dirname(normalizePath(f,mustWork=TRUE)))
  }
  normalizePath(getwd()) # for selected-line execution, setwd() to this code folder
}
project_dir <- locate_project()
data_root <- Sys.getenv("KIDNEY_CELLRANGER_DIR",file.path(project_dir,"data","cellranger"))
processed_root <- Sys.getenv("KIDNEY_MOUSE_PROCESSED",file.path(project_dir,"data","processed"))
# Optional downloaded, annotated checkpoints. NULL checks the canonical name in processed_root.
# Set paths here when GEO filenames differ; annotations and embeddings are retained.
processed_rds <- lapply(c(atlas="KIDNEY_MOUSE_ATLAS_RDS",gec="KIDNEY_MOUSE_GEC_RDS",
  myeloid="KIDNEY_MOUSE_MYELOID_RDS",macrophage="KIDNEY_MOUSE_MACROPHAGE_RDS",
  fibroblast="KIDNEY_MOUSE_FIBROBLAST_RDS",pt_cellchat="KIDNEY_MOUSE_PT_RDS"),
  function(key) { p <- Sys.getenv(key, ""); if(nzchar(p)) p else NULL })
results_root <- Sys.getenv("KIDNEY_RESULTS_DIR", file.path(project_dir,"results"))
gene_list_dir <- file.path(project_dir,"gene_sets")
conditions <- c("sham","ADN","ADN_Stattic")
condition_cols <- c(sham="#62C15B",ADN="#E66342",ADN_Stattic="#3F6EFF")
# TRUE archives an incomplete stage before starting it again. Completed stages
# are never overwritten. Use this only after correcting/reviewing the failure.
restart_failed_stage <- FALSE
# Copy a token only AFTER reviewing the saved marker/QC table and proposed map.
# Cluster numbers alone are insufficient, even if the cluster count is identical.
reviewed_cluster_maps <- list()
# Reference mappings from the original study scripts. No external config files.
annotation_table <- utils::read.delim(text = "stage	cluster	cell_type
atlas	0	PST
atlas	1	Myeloid
atlas	2	Endo
atlas	3	LOH
atlas	4	Myeloid
atlas	5	PCT
atlas	6	CNT
atlas	7	Myeloid
atlas	8	T & NK
atlas	9	CNT
atlas	10	T & NK
atlas	11	Myeloid
atlas	12	PST
atlas	13	GEC
atlas	14	Myeloid
atlas	15	T & NK
atlas	16	IC
atlas	17	Myeloid
atlas	18	Endo
atlas	19	GEC
atlas	20	Fibroblast
atlas	21	Proliferating
atlas	22	B
atlas	23	Myeloid
atlas	24	LOH
atlas	25	LOH
atlas	26	Podo
atlas	27	IC
atlas	28	Mesangial
atlas	29	Myeloid
atlas	30	Fibroblast
atlas	31	Fibroblast
atlas	32	Myeloid
gec	0	GEC0
gec	1	GEC1
myeloid	0	Macrophage
myeloid	1	Macrophage
myeloid	2	Macrophage
myeloid	3	DC
myeloid	4	Granulocyte
myeloid	5	Monocyte
myeloid	6	DC
myeloid	7	Macrophage
myeloid	8	DC
myeloid	9	Macrophage
myeloid	10	DC
macrophage	0	Ccl4+ Macrophage
macrophage	1	Dlc1+ Macrophage
macrophage	2	C1qb+ Macrophage
macrophage	3	Monocyte
macrophage	4	Arg1+ Macrophage
macrophage	5	Cxcl10+ Macrophage
fibroblast	0	Myofibroblast
fibroblast	1	Proto-Myofibroblast
fibroblast	2	Pericyte
fibroblast	3	Quiescent Fibroblast
pt_cellchat	0	PST
pt_cellchat	1	PST
pt_cellchat	2	PCT
pt_cellchat	3	PCT
pt_cellchat	4	PST
",
    check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
marker_panels <- utils::read.delim(text = "panel	gene
atlas	Slc9a3
atlas	Slc4a4
atlas	Ppara
atlas	Acox1
atlas	Egfr
atlas	Slc34a1
atlas	Slc22a2
atlas	Slc22a8
atlas	Ggt1
atlas	Slc22a6
atlas	Slc7a13
atlas	Slc12a1
atlas	Umod
atlas	Aqp2
atlas	Slc12a3
atlas	Calb1
atlas	Atp6v1g3
atlas	Slc26a4
atlas	Atp6v0d2
atlas	Clic3
atlas	Nphs2
atlas	Nphs1
atlas	Sfrp2
atlas	Itga8
atlas	Gata3
atlas	Cd93
atlas	Igfbp3
atlas	Robo4
atlas	Ehd3
atlas	Tbx3
atlas	Col6a3
atlas	Dcn
atlas	Col1a1
atlas	Acta2
atlas	Ppp1r14a
atlas	Ccl12
atlas	Clec4a1
atlas	Cd68
atlas	S100a9
atlas	S100a8
atlas	Cxcr2
atlas	Mmp9
atlas	Nkg7
atlas	Icos
atlas	Cd3e
atlas	Cd3g
atlas	Themis
atlas	Il7r
atlas	Cd79a
atlas	Iglc2
atlas	Iglc3
atlas	Pou2af1
atlas	Mki67
atlas	Hist1h3c
macrophage	Chil3
macrophage	Vcan
macrophage	Ly6c2
macrophage	Arg1
macrophage	Vegfa
macrophage	Fn1
macrophage	Ccl4
macrophage	Cd83
macrophage	Atf3
macrophage	Dlc1
macrophage	Nrp1
macrophage	Psd3
macrophage	C1qb
macrophage	C1qa
macrophage	C1qc
macrophage	Cxcl10
macrophage	Ccl5
macrophage	Ifit3b
",
    check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
config <- list(.root=project_dir,output_dir=results_root,gene_set_dir=gene_list_dir,
  seed=1234L,conditions=conditions,condition_colors=condition_cols,
  mouse=list(annotation_table=annotation_table, marker_panels=marker_panels,
    processed_rds=processed_rds, processed_root=processed_root,
    reviewed_cluster_maps=reviewed_cluster_maps,
    dims=1:30,pca_seed=42L,umap_seed=42L,clustering_seed=0L,
    module_nbin=24L,module_ctrl=100L,module_seed=1L,dunn_adjust="bonferroni",
    proportion_permutations=1000L,proportion_fdr=0.05,proportion_log2fc=log2(1.5)))
set.seed(config$seed)
macro_rds <- NULL
fibro_rds <- NULL
pt_rds <- NULL
gec_rds <- NULL
selected_lr <- NULL # Optional character vector of CellChat interaction_name IDs in plot order.

stage_id <- "08"

# ---------------------------------------------------------------------------
# 3. Reproducibility and analysis helpers
# ---------------------------------------------------------------------------
need <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "),
                            ". See README.md.", call. = FALSE)
}

path_at <- function(root, x) {
  if (is.null(x)) return(NULL)
  if (grepl("^(~|/|[A-Za-z]:[/\\\\])", x)) path.expand(x) else file.path(root, x)
}

read_tsv <- function(file) {
  if (!file.exists(file)) stop("Missing input: ", file, call. = FALSE)
  read.delim(file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

check_columns <- function(x, required, label) {
  absent <- setdiff(required, names(x))
  if (length(absent)) stop(label, " lacks columns: ", paste(absent, collapse = ", "))
}

check_ids <- function(x, label) {
  if (anyNA(x) || any(!nzchar(trimws(x))) || anyDuplicated(x))
    stop(label, " must be nonempty, nonmissing, unique IDs.")
}

file_hash <- function(x) {
  need("digest")
  if (!file.exists(x)) stop("Cannot hash missing input: ", x)
  digest::digest(file = x, algo = "sha256")
}

hash_manifest <- function(paths, out, relative_to = NULL) {
  paths <- unique(paths)
  if (anyNA(paths) || any(!nzchar(paths))) stop("Invalid path in required input manifest.")
  missing <- paths[!file.exists(paths) | dir.exists(paths)]
  if (length(missing)) stop("Required files missing or not regular files: ", paste(missing, collapse=", "))
  labels <- normalizePath(paths, mustWork = TRUE)
  if (!is.null(relative_to)) labels <- substring(labels, nchar(normalizePath(relative_to)) + 2L)
  write_tsv(data.frame(file = labels, bytes = file.info(paths)$size,
                       sha256 = vapply(paths, file_hash, character(1))), out)
}

stage_dir <- function(config, stage) file.path(path_at(config$.root, config$output_dir), stage)

init_stage <- function(config, stage, inputs = character()) {
  out <- stage_dir(config, stage)
  if (dir.exists(out) && !file.exists(file.path(out,"SUCCESS")) && isTRUE(restart_failed_stage)) {
    archived <- paste0(out,"_failed_",format(Sys.time(),"%Y%m%d_%H%M%S"))
    if(file.exists(archived)||!file.rename(out,archived)) stop("Could not archive failed stage.")
  }
  if (dir.exists(out) || file.exists(out))
    stop("Stage directory already exists: ", out,
         ". Choose a new output_dir; failed-run files must not mix with a new run.")
  missing <- inputs[!file.exists(inputs) | dir.exists(inputs)]
  if (anyNA(inputs) || length(missing)) stop("Required stage inputs missing: ", paste(missing, collapse=", "))
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  if (!dir.create(out, showWarnings = FALSE)) stop("Could not claim new stage directory: ", out)
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_start.txt"))
  dput(config, file.path(out, "run_config.R"))
  hash_manifest(inputs, file.path(out, "input_sha256.tsv"))
  code <- c(list.files(config$.root,full.names=TRUE,pattern="\\.(R|txt)$"),
    list.files(file.path(config$.root,"gene_sets"),full.names=TRUE))
  hash_manifest(code, file.path(out, "code_sha256.tsv"), relative_to = config$.root)
  set.seed(config$seed)
  out
}

finish_stage <- function(out) {
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  loaded <- sort(loadedNamespaces())
  write_tsv(data.frame(package = loaded,
    version = vapply(loaded, function(p) as.character(getNamespaceVersion(p)), character(1)),
    path = vapply(loaded, function(p) if(p=="base") file.path(R.home(),"library/base") else getNamespaceInfo(asNamespace(p), "path"), character(1))),
            file.path(out, "loaded_packages.tsv"))
  files <- list.files(out, full.names=TRUE, recursive=TRUE)
  files <- files[!basename(files) %in% c("SUCCESS", "output_sha256.tsv")]
  hash_manifest(files, file.path(out, "output_sha256.tsv"), relative_to=out)
  writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), file.path(out, "SUCCESS"))
}

gene_sets <- function(config) {
  base <- path_at(config$.root, config$gene_set_dir)
  checks <- read_tsv(file.path(base, "SHA256SUMS.tsv"))
  check_columns(checks, c("file", "bytes", "sha256"), "Gene-set checksums")
  check_ids(checks$file, "Gene-set checksum filenames")
  if (any(basename(checks$file) != checks$file)) stop("Gene-set checksum paths must be filenames.")
  if (!all(c("gene_sets.tsv", "gene_set_manifest.tsv") %in% checks$file)) stop("Missing canonical gene-set checksums.")
  paths <- file.path(base, checks$file)
  if (any(!file.exists(paths))) stop("A checksummed gene-set file is missing.")
  if (any(file.info(paths)$size != checks$bytes) ||
      any(vapply(paths, file_hash, character(1)) != checks$sha256))
    stop("Gene-set export changed. Re-export and review Table S1; do not mix old aliases and new lists.")
  d <- read_tsv(file.path(base, "gene_sets.tsv"))
  manifest <- read_tsv(file.path(base, "gene_set_manifest.tsv"))
  check_columns(d, c("set_id", "gene", "order"), "gene_sets.tsv")
  check_columns(manifest, c("set_id", "n_genes"), "Gene-set manifest")
  check_ids(manifest$set_id, "Gene-set IDs")
  if (anyNA(d) || any(!nzchar(d$set_id)) || any(!nzchar(d$gene)) ||
      any(grepl("[[:space:]]", d$gene))) stop("Missing/invalid gene-set member.")
  if (anyDuplicated(paste(d$set_id, d$gene, sep = "::"))) stop("Duplicate gene within set.")
  if (!setequal(unique(d$set_id), manifest$set_id)) stop("Gene-set IDs differ from manifest.")
  result <- setNames(lapply(manifest$set_id, function(id) {
    z <- d[d$set_id == id, , drop=FALSE]
    if (!is.numeric(z$order) || anyDuplicated(z$order) ||
        !setequal(z$order, seq_len(nrow(z)))) stop("Invalid member order: ", id)
    z$gene[order(z$order)]
  }), manifest$set_id)
  if (any(lengths(result) != manifest$n_genes)) stop("Gene-set counts differ from manifest.")
  result
}

match_genes <- function(sets, available, out, context) {
  rows <- lapply(names(sets), function(nm) data.frame(set_id = nm, gene = sets[[nm]],
    present = sets[[nm]] %in% available, context = context))
  d <- do.call(rbind, rows)
  write_tsv(d, file.path(out, paste0(context, "_gene_matching.tsv")))
  result <- lapply(sets, function(g) intersect(g, available))
  if (any(lengths(result) == 0L)) stop("No genes matched in: ", paste(names(result)[lengths(result)==0L], collapse=", "))
  result
}

safe_name <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

save_plot <- function(p, out, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(out, paste0(name, ".pdf")), plot = p, width = width,
                 height = height, device = grDevices::pdf, useDingbats = FALSE)
}

require_completed <- function(config, stage) {
  if (!file.exists(file.path(stage_dir(config, stage), "SUCCESS")))
    stop("Upstream stage is incomplete: ", stage, ". A partial checkpoint is not a completed input.")
}
# A cluster number is not a cell identity. Bind manual marker review to the
# actual cell-to-cluster assignments AND the proposed annotation/exclusion map.

cluster_review_token <- function(cells, ids, map) {
  need("digest")
  check_ids(cells, "Cluster-review cell IDs")
  if (length(cells) != length(ids) || anyNA(ids)) stop("Invalid cluster assignments.")
  check_columns(map, c("cluster", "cell_type"), "Cluster review map")
  if (anyNA(map) || anyDuplicated(map$cluster) ||
      !setequal(unique(as.character(ids)), as.character(map$cluster)))
    stop("Cluster IDs differ from reference map; inspect markers and update the map.")
  ci <- order(cells, method="radix"); mi <- order(as.character(map$cluster), method="radix")
  digest::digest(list(cells=as.character(cells[ci]), clusters=as.character(ids[ci]),
    cluster=as.character(map$cluster[mi]), label=as.character(map$cell_type[mi])),
    algo="sha256", serializeVersion=2)
}

require_cluster_review <- function(cells, ids, map, config, stage, out) {
  token <- cluster_review_token(cells, ids, map)
  write_tsv(data.frame(cell=cells, cluster=as.character(ids)),
    file.path(out,paste0(stage,"_review_assignments.tsv")))
  write_tsv(map, file.path(out,paste0(stage,"_review_map.tsv")))
  write_tsv(data.frame(stage=stage, review_token=token),
    file.path(out,paste0(stage,"_review_token.tsv")))
  if (!identical(unname(config$mouse$reviewed_cluster_maps[[stage]]), token))
    stop("Marker review required for ", stage, ". Inspect saved markers/map/assignments, then record ",
      "the reviewed token in config$mouse$reviewed_cluster_maps[['",stage,"']]. ",
      "Set restart_failed_stage=TRUE for this incomplete stage. See README.md.")
  invisible(token)
}

# Only completed upstream results are accepted. Explicit external checkpoints
# may be supplied in the Configuration block; their hashes are recorded.
reference_input <- function(key, filename) {
  explicit <- config$mouse$processed_rds[[key]]
  if (!is.null(explicit)) {
    p <- path_at(project_dir,explicit)
    if (!file.exists(p)) stop("Configured processed RDS not found: ",p)
    return(p)
  }
  # Prefer exact original annotated names over generic checkpoint names.
  # Original working folders may also contain an unannotated Macro_subset.rds.
  # Explicit per-object paths above always take precedence.
  # Never select a checkpoint by a fuzzy filename match.
  original_names <- c(atlas="Stattic_after_levels2.rds",gec="GEC_subset_after_levels.rds",
    myeloid="Myeloid.rename.rds",macrophage="Macro.rename.rds",
    fibroblast="Fibro_subset_after_levels.rds",pt_cellchat="PT_subset_after_levels.rds")
  original <- unname(original_names[key])
  if(length(original)==1L && !is.na(original)) {
    p <- file.path(processed_root,original)
    if(file.exists(p)) return(p)
  }
  p <- file.path(processed_root,filename)
  if(file.exists(p)) p else NULL
}
validate_reference_subset <- function(o, parent, column, labels) {
  expected <- colnames(parent)[as.character(parent[[column]][,1]) %in% labels]
  if (!all(colnames(o) %in% expected))
    stop("Processed subset includes cells outside its expected parent lineage.")
  ix <- match(colnames(o),colnames(parent))
  if (!identical(as.character(o$condition),as.character(parent$condition[ix])))
    stop("Condition labels differ between processed subset and parent.")
  if (!"umap" %in% names(o@reductions)) stop("Processed subset lacks the published UMAP.")
  invisible(TRUE)
}

input_rds <- function(stage, file, override=NULL) {
  if (!is.null(override)) return(path_at(project_dir,override))
  # A deposited checkpoint can replace a computed upstream stage.
  path <- file.path(stage_dir(config,stage),file)
  if (!file.exists(file.path(stage_dir(config,stage),"SUCCESS"))) {
    key <- switch(file, "Stattic_annotated.rds"="atlas", "GEC_subset.rds"="gec",
      "Myeloid_subset.rds"="myeloid", "Macro_subset.rds"="macrophage",
      "Fibro_subset.rds"="fibroblast", "PT_CellChat_subset.rds"="pt_cellchat", "integrated")
    candidate <- reference_input(key,file)
    if (!is.null(candidate)) return(candidate)
  }
  require_completed(config,stage)
  if(!file.exists(path))stop("Missing upstream checkpoint: ",path)
  checks<-read_tsv(file.path(stage_dir(config,stage),"output_sha256.tsv"))
  expected<-checks$sha256[checks$file==file]
  if(length(expected)!=1L||file_hash(path)!=expected)stop("Upstream checkpoint changed after completion: ",path)
  path
}
read_mouse <- function(path, column=NULL, levels=NULL, require_normalized=TRUE) {
  o<-readRDS(path)
  if(!inherits(o,"Seurat")||!"RNA"%in%names(o@assays))stop("Expected a Seurat RNA object: ",path)
  check_ids(colnames(o),"Cell IDs")
  if(!"condition"%in%names(o[[]])||anyNA(o$condition)||!setequal(as.character(o$condition),conditions))
    stop("Object must contain sham, ADN, ADN_Stattic with no missing condition labels.")
  o$condition<-factor(as.character(o$condition),levels=conditions)
  SeuratObject::DefaultAssay(o)<-"RNA"
  if(require_normalized) {
    layers<-SeuratObject::Layers(o[["RNA"]])
    if(!all(c("counts","data")%in%layers)||any(grepl("^(counts|data)\\.",layers)))
      stop("Expected joined RNA counts and normalized data layers. Join/normalize upstream first.")
  }
  if(!is.null(column))o<-set_identity(o,column,levels)
  o
}
check_versions <- function(pins) {
  actual<-vapply(names(pins),function(p)as.character(utils::packageVersion(p)),character(1))
  if(any(actual!=pins))stop("Reference package mismatch: ",paste(names(pins)[actual!=pins],collapse=", "),
    ". Use the documented environment; do not label an unvalidated version change an exact reproduction.")
  if(as.character(getRversion())!="4.4.2")warning("Reference R version is 4.4.2; current R differs. Runtime versions are recorded in sessionInfo.txt.")
}
write_checkpoint <- function(o,filename)saveRDS(o,file.path(out,filename))
# Mouse scRNA-seq: one pooled library per condition. Cell-level tests are exploratory.
CELL_TYPES <- c("PCT", "PST", "LOH", "CNT", "IC", "Podo", "Mesangial", "Endo",
                "GEC", "Fibroblast", "Myeloid", "T & NK", "B", "Proliferating")
MACRO_TYPES <- c("Monocyte", "Arg1+ Macrophage", "Ccl4+ Macrophage", "Dlc1+ Macrophage",
                 "C1qb+ Macrophage", "Cxcl10+ Macrophage")
FIBRO_TYPES <- c("Quiescent Fibroblast", "Proto-Myofibroblast", "Myofibroblast", "Pericyte")

apply_annotation <- function(obj, config, stage, column, out) {
  map <- config$mouse$annotation_table
  map <- map[map$stage == stage, ]
  if (!nrow(map) || anyDuplicated(map$cluster)) stop("Invalid annotation map: ", stage)
  ids <- as.character(SeuratObject::Idents(obj))
  if (!setequal(unique(ids), as.character(map$cluster))) stop("Cluster IDs differ from the published ",
    stage, " map. Inspect marker outputs and revise annotations; never assign labels by position.")
  require_cluster_review(colnames(obj), ids, map, config, stage, out)
  obj[[column]] <- map$cell_type[match(ids, as.character(map$cluster))]
  SeuratObject::Idents(obj) <- column
  obj
}

set_identity <- function(o, col, order = NULL) {
  if (!col %in% names(o[[]])) stop("Missing annotation column: ", col)
  x <- as.character(o[[col]][,1])
  x[x == "Dlc+ Macrophage"] <- "Dlc1+ Macrophage"
  if (anyNA(x) || any(!nzchar(x))) stop("Missing cell annotations in ", col)
  if (!is.null(order) && !all(x %in% order)) stop("Unexpected labels in ", col)
  o[[col]] <- if (is.null(order)) x else factor(x, levels = order)
  SeuratObject::Idents(o) <- col
  o
}

cluster_subset <- function(o,config,resolution,out,name) {
  # Use inherited Harmony dimensions, as in the original scripts.
  # A new PCA on a subset is not used by this neighbor graph.
  if (!"harmony" %in% names(o@reductions)) stop("Subset lacks inherited Harmony dimensions.")
  o <- Seurat::FindNeighbors(o,reduction="harmony",dims=config$mouse$dims)
  o <- Seurat::FindClusters(o,resolution=resolution,random.seed=config$mouse$clustering_seed)
  o <- Seurat::RunUMAP(o,reduction="harmony",dims=config$mouse$dims,seed.use=config$mouse$umap_seed)
  mk <- Seurat::FindAllMarkers(o,only.pos=TRUE,min.pct=.25,logfc.threshold=1,assay="RNA",test.use="wilcox")
  write_tsv(mk,file.path(out,paste0(name,"_cluster_markers.tsv")))
  saveRDS(o,file.path(out,paste0(name,"_unannotated.rds")))
  o
}
incoming_strength <- function(edges,sources=c("PST","GEC")) {
  check_columns(edges,c("source","target","prob"),"CellChat edge table")
  d <- edges[edges$source %in% sources,,drop=FALSE]
  if(!nrow(d)) return(setNames(numeric(),character()))
  tapply(d$prob,as.character(d$target),sum)
}

incoming_comparison <- function(a,b,pseudocount=.01) {
  # Match the original inner join: only targets in both retained edge tables.
  targets <- intersect(names(a),names(b))
  x <- data.frame(target=targets,denominator=as.numeric(a[targets]),numerator=as.numeric(b[targets]))
  x$log2FC <- log2((x$numerator+pseudocount)/(x$denominator+pseudocount))
  x$rank_abs <- rank(-abs(x$log2FC),ties.method="average")
  x[order(x$rank_abs),,drop=FALSE]
}

relative_centrality <- function(x) {
  # Matches role-wise max scaling; an all-zero role has undefined importance.
  z <- sweep(x,1,apply(x,1,max),"/")
  z[!is.finite(z)] <- NA_real_
  z
}

combine_cc_objects <- function(objects,columns) {
  if(length(objects)!=length(columns)) stop("CellChat annotation input mismatch.")
  all_cells<-unlist(lapply(objects,colnames),use.names=FALSE)
  if(anyDuplicated(all_cells))stop("Overlapping cell IDs across communication compartments.")
  for(i in seq_along(objects)) {
    o <- set_identity(objects[[i]],columns[i])
    o$cc_cell_type <- as.character(SeuratObject::Idents(o))
    objects[[i]] <- o
  }
  o <- if(length(objects)==1L) objects[[1]] else merge(objects[[1]],y=objects[-1],add.cell.ids=names(objects))
  # Recompute only the expression layer consumed by CellChat. Integration and
  # reclustering of this union do not enter ligand-receptor probabilities.
  o <- Seurat::NormalizeData(o,normalization.method="LogNormalize",scale.factor=10000)
  if(inherits(o[["RNA"]],"Assay5")) o[["RNA"]] <- SeuratObject::JoinLayers(o[["RNA"]])
  o
}

run_cellchat_network <- function(o,cell_levels,name,min_cells,out) {
  oldplan<-future::plan();on.exit(future::plan(oldplan),add=TRUE)
  future::plan("sequential")
  database<-new.env(parent=emptyenv())
  utils::data("CellChatDB.mouse",package="CellChat",envir=database)
  utils::data("PPI.mouse",package="CellChat",envir=database)
  if(is.null(database$CellChatDB.mouse)||is.null(database$PPI.mouse))stop("Mouse database unavailable.")
  write_tsv(database$CellChatDB.mouse$interaction,file.path(out,paste0(name,"_database_interactions.tsv")))
  result<-setNames(vector("list",length(conditions)),conditions)
  for(cond in conditions) {
    o1<-subset(o,cells=colnames(o)[as.character(o$condition)==cond])
    lab<-as.character(o1$cc_cell_type)
    if(anyNA(lab)||!all(lab%in%cell_levels))stop("Unknown CellChat cell labels.")
    meta<-data.frame(labels=factor(lab,levels=cell_levels[cell_levels%in%lab]),row.names=colnames(o1))
    expr<-SeuratObject::GetAssayData(o1,assay="RNA",layer="data")
    if(!identical(colnames(expr),rownames(meta)))stop("CellChat metadata/expression order mismatch.")
    stem<-paste0("cellchat_",name,"_",cond)
    counts<-as.data.frame(table(meta$labels));names(counts)<-c("cell_type","n_cells")
    write_tsv(counts,file.path(out,paste0(stem,"_counts.tsv")))
    cc<-CellChat::createCellChat(object=expr,meta=meta,group.by="labels")
    cc@DB<-database$CellChatDB.mouse
    cc<-CellChat::subsetData(cc)
    cc<-CellChat::identifyOverExpressedGenes(cc)
    cc<-CellChat::identifyOverExpressedInteractions(cc)
    cc<-CellChat::smoothData(cc,adj=database$PPI.mouse)
    # raw.use=TRUE consumes normalized RNA signaling expression, not smoothed
    # PPI expression and not Harmony coordinates. Full CellChatDB.mouse is used.
    cc<-CellChat::computeCommunProb(cc,type="triMean",raw.use=TRUE,
      population.size=FALSE,nboot=100,seed.use=1L)
    cc<-CellChat::filterCommunication(cc,min.cells=min_cells)
    cc<-CellChat::computeCommunProbPathway(cc,thresh=.05)
    cc<-CellChat::aggregateNet(cc,thresh=.05)
    cc<-CellChat::netAnalysis_computeCentrality(cc,slot.name="netP",thresh=.05)
    saveRDS(cc,file.path(out,paste0(stem,"_final.rds")))
    # Original per-condition candidate probabilities/P values remain in @net.
    edges<-CellChat::subsetCommunication(cc,thresh=.05)
    write_tsv(edges,file.path(out,paste0(stem,"_retained_edges.tsv")))
    result[[cond]]<-cc
  }
  result
}
restricted_edges <- function(objects,sources,targets,pathways=NULL,threshold=.05) {
  rows<-list()
  for(cond in names(objects)) {
    cc<-objects[[cond]]
    # Extract once, then restrict by names. This handles absent source/target
    # populations without accidentally reassigning a numerical index.
    d<-CellChat::subsetCommunication(cc,thresh=threshold)
    d<-d[d$source%in%sources & d$target%in%targets,,drop=FALSE]
    if(!is.null(pathways))d<-d[d$pathway_name%in%pathways,,drop=FALSE]
    d$condition<-rep(cond,nrow(d));rows[[cond]]<-d
  }
  do.call(rbind,rows)
}
bubble_export <- function(objects,sources,targets,pathways,out,name,threshold=.05,pairs=NULL) {
  d<-restricted_edges(objects,sources,targets,pathways,threshold)
  if(!is.null(pairs))d<-d[d$interaction_name%in%pairs,,drop=FALSE]
  write_tsv(d,file.path(out,paste0(name,"_raw_edges.tsv")))
  if(!nrow(d)) {
    writeLines("No retained interactions for the requested groups/pathways.",file.path(out,paste0(name,"_unavailable.txt")))
    return(invisible(NULL))
  }
  merged<-CellChat::mergeCellChat(objects,add.names=names(objects))
  pair_df<-if(is.null(pairs))NULL else data.frame(interaction_name=pairs)
  p<-CellChat::netVisual_bubble(merged,sources.use=sources,targets.use=targets,
    comparison=1:3,signaling=if(is.null(pairs))pathways else NULL,pairLR.use=pair_df,
    thresh=threshold,angle.x=45,remove.isolate=FALSE,
    # CellChat 2.1.2 separator generation fails for a single source-target pair.
    line.on=nrow(unique(d[c("source","target")]))>1L)
  # Reserve space for long source-target labels at the left and bottom edges.
  p<-p+ggplot2::theme(plot.margin=ggplot2::margin(12,16,20,80,unit="pt"))
  save_plot(p,out,name,14,6)
  if(!is.null(p$data))write_tsv(p$data,file.path(out,paste0(name,"_scaled_plot_data.tsv")))
}
rank_export <- function(objects,sources,targets,out,name) {
  previous_warn<-getOption("warn");on.exit(options(warn=previous_warn),add=TRUE)
  d<-restricted_edges(objects,sources,targets)
  write_tsv(d,file.path(out,paste0(name,"_selected_raw_edges.tsv")))
  if(!nrow(d)) {
    writeLines("No retained source-target interactions; no information-flow plot.",file.path(out,paste0(name,"_unavailable.txt")))
    return(invisible(NULL))
  }
  merged<-CellChat::mergeCellChat(objects,add.names=names(objects))
  # All inferred pathways are compared. Scaling is the rankNet default; it is
  # not an absolute unit shared between separately inferred networks.
  z<-CellChat::rankNet(merged,mode="comparison",comparison=1:3,
    sources.use=sources,targets.use=targets,stacked=FALSE,do.stat=TRUE,
    show.raw=FALSE,return.data=TRUE,color.use=unname(condition_cols))
  write_tsv(z$signaling.contribution,file.path(out,paste0(name,"_raw_and_scaled_flow.tsv")))
  save_plot(z$gg.obj,out,name,7,8)
}

# ---------------------------------------------------------------------------
# 4. Analysis
# ---------------------------------------------------------------------------

check_versions(c(Seurat="5.3.0",CellChat="2.1.2"))
inputs<-c(macro=input_rds("05","Macro_subset.rds",macro_rds),
  fibro=input_rds("07","Fibro_subset.rds",fibro_rds),
  pt=input_rds("06","PT_CellChat_subset.rds",pt_rds),
  gec=input_rds("04","GEC_subset.rds",gec_rds))
config$selected_lr <- selected_lr
out<-init_stage(config,stage_id,inputs)
Macro<-read_mouse(inputs[["macro"]],"Macro_CellType",MACRO_TYPES)
Fibro<-read_mouse(inputs[["fibro"]],"Fibro_CellType",FIBRO_TYPES)
PT<-read_mouse(inputs[["pt"]],"PT_CellType",c("PCT","PST"))
GEC<-read_mouse(inputs[["gec"]],"GEC_CellType",c("GEC0","GEC1"))
mf_obj<-combine_cc_objects(list(macro=Macro,fibro=Fibro),c("Macro_CellType","Fibro_CellType"))
write_checkpoint(mf_obj,"Macro_Fibro_combined.rds")
mf<-run_cellchat_network(mf_obj,c(MACRO_TYPES,FIBRO_TYPES),"macro_fibro",1L,out)

# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>  Figure 5F: relative roles in the THBS pathway  <<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================

for(cond in conditions) {
  cc<-mf[[cond]]
  if(!"THBS"%in%names(cc@netP$centr)) {
    writeLines("No retained THBS pathway; centrality unavailable.",file.path(out,paste0("Fig5F_",cond,"_unavailable.txt")));next
  }
  raw<-do.call(rbind,cc@netP$centr$THBS[c("outdeg","indeg","flowbet","info")])
  rownames(raw)<-c("Sender","Receiver","Mediator","Influencer")
  rel<-relative_centrality(raw)
  write_tsv(data.frame(role=rownames(raw),raw,check.names=FALSE),file.path(out,paste0("Fig5F_",cond,"_raw.tsv")))
  write_tsv(data.frame(role=rownames(rel),rel,check.names=FALSE),file.path(out,paste0("Fig5F_",cond,"_relative.tsv")))
  if(!any(is.finite(rel)))next
  pdf(file.path(out,paste0("Fig5F_THBS_",cond,".pdf")),width=9,height=4)
  CellChat::netAnalysis_signalingRole_network(cc,signaling="THBS")
  dev.off()
}
# THBS is a CellChat pathway, whereas THBS1 is a ligand/protein. Role scores
# are max-normalized separately within each role and condition; they are not
# an absolute between-condition scale. An all-zero role is unavailable (NA).

# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>  Figures S11D-E: pathway flow and selected pathway LR interactions  <<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================

rank_export(mf,"Arg1+ Macrophage",FIBRO_TYPES,out,"FigS11D")
pairs<-NULL
if(!is.null(selected_lr)) {
  pairs<-selected_lr;check_ids(pairs,"Selected LR IDs")
  db<-mf[[1]]@DB$interaction
  if(!all(pairs%in%rownames(db)))stop("Selected LR IDs are absent from CellChatDB.mouse 2.1.2.")
  if(!all(db[pairs,"pathway_name"]%in%c("THBS","CypA","FN1")))stop("S11E selected LR vector contains another pathway.")
}
bubble_export(mf,"Arg1+ Macrophage",FIBRO_TYPES,c("THBS","CypA","FN1"),out,"FigS11E",pairs=pairs)
# Default: all retained LR rows within these three pathways. A historical
# manually selected subset/order can be set in the embedded selected_lr vector;
# this script does not infer row selection from an image or silently drop rows.

# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>  Figure S11F: cross-compartment expression  <<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================

GEC$expression_CellType<-"GEC" # unsplit expression display, distinct from 5E/S10D
all_obj<-combine_cc_objects(list(pt=PT,gec=GEC,macro=Macro,fibro=Fibro),
  c("PT_CellType","expression_CellType","Macro_CellType","Fibro_CellType"))
all_obj<-set_identity(all_obj,"cc_cell_type",c("PCT","PST","GEC",MACRO_TYPES,FIBRO_TYPES))
write_checkpoint(all_obj,"Ferroptotic_Macro_Fibro_combined.rds")
save_plot(Seurat::VlnPlot(all_obj,features=c("Spp1","Cd44","Thbs1"),split.by="condition",
  cols=condition_cols,pt.size=0,ncol=1),out,"FigS11F_Spp1_Cd44_Thbs1",20,10)
x<-SeuratObject::FetchData(all_obj,vars=c("Spp1","Cd44","Thbs1"),layer="data",clean=FALSE)
write_tsv(data.frame(cell=rownames(x),condition=all_obj$condition,cell_type=all_obj$cc_cell_type,x),file.path(out,"FigS11F_expression.tsv"))
writeLines(c("Relative centrality measures network roles and does not establish absolute pathway dominance.",
  "THBS1 expression is assessed across multiple cell populations.",
  "LR P values are within-condition permutation results, not condition-difference tests.",
  "Pathway flow and relative centrality tables are saved separately."),file.path(out,"CellChat_interpretation.txt"))

finish_stage(out)

###############################################################################
# End of script 08.
###############################################################################
