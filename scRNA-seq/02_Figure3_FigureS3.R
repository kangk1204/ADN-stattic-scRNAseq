###############################################################################
# Script:      02_Figure3_FigureS3.R
# Project:     Mouse kidney scRNA-seq (Sham / ADN / ADN + Stattic)
# Revision:    Current integrated manuscript, 2026-10-12
# Description: Atlas annotation, module scores and recovered-cell proportions.
#
# Manuscript panels / source data:
#   Figures 3A-F; S3B-D -- mouse single-cell analyses only.
#
# Input:       01/Kidney_HarmonyIntegrated_dim30_res1.rds; current Table S1 exports.
# Design:      One pooled library from two mice per condition (three libraries).
###############################################################################

# ---------------------------------------------------------------------------
# 1. Libraries
# ---------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(rstatix)
library(scProportionTest)
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
integrated_rds <- NULL # optional absolute path to an already integrated reference object

stage_id <- "02"

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

add_scores <- function(o,sets,config,out,label) {
  need("Seurat")
  match_genes(sets,rownames(o[["RNA"]]),out,label)
  for(nm in names(sets)) {
    g <- intersect(sets[[nm]],rownames(o[["RNA"]]))
    o <- Seurat::AddModuleScore(o,features=list(g),assay="RNA",slot="data",
      name=paste0(nm,"_score"),nbin=config$mouse$module_nbin,ctrl=config$mouse$module_ctrl,
      seed=config$mouse$module_seed,search=FALSE)
  }
  o
}

cell_tests <- function(d,adjust="bonferroni") {
  need("rstatix")
  if(length(unique(d$value)) < 2L) return(list(
    omnibus=data.frame(p=NA_real_,status="constant values; test undefined"),pairs=data.frame()))
  if(length(unique(d$condition)) < 2L) return(list(
    omnibus=data.frame(p=NA_real_,status="fewer than two groups"),pairs=data.frame()))
  list(omnibus=as.data.frame(rstatix::kruskal_test(d,value~condition)),
       pairs=as.data.frame(rstatix::dunn_test(d,value~condition,p.adjust.method=adjust)))
}

average_heatmap <- function(o,genes,group,order,out,name,cluster_labels=NULL) {
  need("pheatmap")
  genes <- unique(genes)
  available <- intersect(genes,rownames(o))
  write_tsv(data.frame(gene=genes,present=genes %in% available),file.path(out,paste0(name,"_genes.tsv")))
  # AverageExpression exponentiates log-normalized data before averaging.
  # Safe temporary group codes avoid Seurat converting underscores in labels.
  ids <- as.character(o[[group]][,1]); codes <- paste0("grp",seq_along(order))
  if(!setequal(unique(ids),order)) stop("Heatmap grouping mismatch: ",name)
  o$aggregation_group <- factor(codes[match(ids,order)],levels=codes)
  a <- as.matrix(Seurat::AverageExpression(o,features=available,group.by="aggregation_group",assays="RNA",layer="data",verbose=FALSE)$RNA)
  a <- a[available,codes,drop=FALSE]; colnames(a) <- order
  write_tsv(data.frame(gene=rownames(a),a,check.names=FALSE),file.path(out,paste0(name,"_linear_mean.tsv")))
  z <- t(scale(t(a))); keep <- apply(z,1,function(x) all(is.finite(x)))
  write_tsv(data.frame(gene=rownames(z),plotted=keep),file.path(out,paste0(name,"_row_scaling.tsv")))
  z <- z[keep,,drop=FALSE]
  if(!nrow(z)) stop("No variable heatmap rows.")
  ann <- if(is.null(cluster_labels)) NULL else data.frame(gene_cluster=cluster_labels[rownames(z)],row.names=rownames(z))
  pheatmap::pheatmap(z,cluster_rows=FALSE,cluster_cols=FALSE,annotation_row=ann,
    show_rownames=nrow(z)<100,filename=file.path(out,paste0(name,".pdf")),width=6,height=8)
}

proportion_exports <- function(o,column,config,out,label) {
  need("scProportionTest")
  counts <- as.data.frame(table(cell_type=o[[column]][,1],condition=o$condition))
  names(counts)[3] <- "n_cells"
  write_tsv(counts,file.path(out,paste0(label,"_counts.tsv")))
  for(pair in list(c("sham","ADN"),c("ADN","ADN_Stattic"))) {
    set.seed(config$seed)
    pt <- scProportionTest::sc_utils(o)
    pt <- scProportionTest::permutation_test(pt,cluster_identity=column,sample_identity="condition",
      sample_1=pair[1],sample_2=pair[2],n_permutations=config$mouse$proportion_permutations)
    stem <- paste(label,pair[1],pair[2],sep="_")
    saveRDS(pt,file.path(out,paste0(stem,".rds")))
    tab <- as.data.frame(pt@results$permutation)
    tab$significant <- tab$FDR < config$mouse$proportion_fdr &
      abs(tab$obs_log2FD) > config$mouse$proportion_log2fc
    write_tsv(tab,file.path(out,paste0(stem,".tsv")))
    save_plot(scProportionTest::permutation_plot(pt,
      FDR_threshold=config$mouse$proportion_fdr,
      log2FD_threshold=config$mouse$proportion_log2fc),out,stem,6,5)
  }
}


violin_exports <- function(o,features,config,out,label,split_col=NULL) {
  for(f in features) {
    if(!f%in%c(rownames(o),names(o[[]])))stop("Feature absent: ",f)
    dd<-SeuratObject::FetchData(o,vars=f,layer="data",clean=FALSE)
    ix<-match(rownames(dd),colnames(o))
    d<-data.frame(cell=rownames(dd),value=dd[[1]],condition=o$condition[ix])
    d$subset<-if(is.null(split_col))factor("all")else factor(as.character(o[[split_col]][ix,1]),levels=levels(o[[split_col]][,1]))
    if(anyNA(d)||any(!is.finite(d$value)))stop("Invalid violin input.")
    write_tsv(d,file.path(out,paste0(label,"_",f,"_values.tsv")))
    annotations<-list()
    for(s in unique(as.character(d$subset))) {
      z<-d[d$subset==s,,drop=FALSE];z$condition<-droplevels(z$condition)
      stat<-cell_tests(z,config$mouse$dunn_adjust)
      stem<-safe_name(paste(label,f,s,sep="_"))
      write_tsv(stat$omnibus,file.path(out,paste0(stem,"_KW.tsv")))
      write_tsv(stat$pairs,file.path(out,paste0(stem,"_Dunn.tsv")))
      if(nrow(stat$pairs)) {
        a<-stat$pairs;a$subset<-s
        a$x1<-match(a$group1,conditions);a$x2<-match(a$group2,conditions)
        span<-max(diff(range(z$value)),0.1)
        a$y<-max(z$value)+seq_len(nrow(a))*.14*span
        a$label<-ifelse(is.na(a$p.adj),"NA",ifelse(a$p.adj<.0001,"****",ifelse(a$p.adj<.001,"***",ifelse(a$p.adj<.01,"**",ifelse(a$p.adj<.05,"*","ns")))))
        annotations[[s]]<-a
      }
    }
    p<-ggplot2::ggplot(d,ggplot2::aes(as.numeric(condition),value,fill=condition,group=condition))+
      ggplot2::geom_violin(scale="width",trim=TRUE)+
      ggplot2::geom_boxplot(width=.12,outlier.shape=NA,fill="white")+
      ggplot2::facet_wrap(~subset,scales="free_y",nrow=1)+
      ggplot2::scale_x_continuous(breaks=1:3,labels=conditions)+
      ggplot2::scale_fill_manual(values=condition_cols)+ggplot2::theme_classic()+
      ggplot2::theme(legend.position="none",axis.text.x=ggplot2::element_text(angle=45,hjust=1))+
      ggplot2::labs(x=NULL,y=f,caption="Cell-level exploratory comparison of pooled preparations; not independent mice.")
    if(length(annotations)) {
      a<-do.call(rbind,annotations)
      a$subset<-factor(a$subset,levels=levels(d$subset))
      p<-p+ggplot2::geom_segment(data=a,ggplot2::aes(x=x1,xend=x2,y=y,yend=y),inherit.aes=FALSE)+
        ggplot2::geom_text(data=a,ggplot2::aes(x=(x1+x2)/2,y=y,label=label),inherit.aes=FALSE,vjust=-.25,size=3)
    }
    save_plot(p,out,safe_name(paste(label,f,sep="_")),if(is.null(split_col))4 else 13,4.5)
  }
}

# ---------------------------------------------------------------------------
# 4. Analysis
# ---------------------------------------------------------------------------

check_versions(c(Seurat="5.3.0",scProportionTest="0.0.0.9000"))
reference_atlas <- reference_input("atlas","Stattic_annotated.rds")
input<-if(!is.null(reference_atlas)) reference_atlas else input_rds("01","Kidney_HarmonyIntegrated_dim30_res1.rds",integrated_rds)
out<-init_stage(config,stage_id,c(input,
  file.path(gene_list_dir,"SHA256SUMS.tsv")))
Stattic<-read_mouse(input)
if ("CellType" %in% names(Stattic[[]])) {
  # Preserve deposited cell labels and the existing UMAP; no re-clustering.
  if (!"umap" %in% names(Stattic@reductions)) stop("Annotated atlas requires its stored UMAP.")
} else {
# Preserve the stored resolution-1 UMAP; change cluster resolution only.
Stattic<-Seurat::FindClusters(Stattic,resolution=1.8,random.seed=config$mouse$clustering_seed)
mk<-Seurat::FindAllMarkers(Stattic,only.pos=TRUE,min.pct=.25,logfc.threshold=1,test.use="wilcox",assay="RNA")
write_tsv(mk,file.path(out,"atlas_cluster_markers.tsv"))
write_checkpoint(Stattic,"atlas_unannotated.rds")
# Apply only after reviewing markers and the cell-to-cluster/map-specific token.
Stattic<-apply_annotation(Stattic,config,"atlas","CellType",out)
}
Stattic<-set_identity(Stattic,"CellType",CELL_TYPES)
write_checkpoint(Stattic,"Stattic_annotated.rds")
write_tsv(data.frame(cell=colnames(Stattic),Stattic[[]]),file.path(out,"cell_annotations.tsv"))
sets<-gene_sets(config)

# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>  Figures 3A-D  <<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================

save_plot(Seurat::DimPlot(Stattic,label=TRUE,pt.size=.01),out,"Fig3A_UMAP_celltype",7,6)
save_plot(Seurat::DimPlot(Stattic,split.by="condition",pt.size=.01),out,"Fig3B_UMAP_split_condition",15,5)
violin_exports(Stattic,"Stat3",config,out,"Fig3C")
violin_exports(Stattic,"Stat3",config,out,"Fig3D","CellType")

# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>  Figures 3E-F and S3D  <<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================

ids<-c("LipidOxidation","EpigeneticRegulation","InflammatoryResponse","ECMOrganization", "GOBP_FERROPTOSIS",paste0("group",1:4))
Stattic<-add_scores(Stattic,sets[ids],config,out,"atlas")
for(nm in ids) {
  panel<-if(nm=="GOBP_FERROPTOSIS")"Fig3F"else if(grepl("^group",nm))"FigS3D"else "Fig3E"
  save_plot(Seurat::FeaturePlot(Stattic,features=paste0(nm,"_score1"),pt.size=.01),out,paste0(panel,"_",nm,"_UMAP"),5,4)
}
# Current Table S1: GOBP_FERROPTOSIS, 32 genes; the same membership as bulk GSVA.
# AddModuleScore and GSVA are distinct scoring methods. This score is not a
# direct measurement of cell death. All atlas cell types are displayed.
save_plot(Seurat::VlnPlot(Stattic,features="GOBP_FERROPTOSIS_score1",pt.size=0)+
  ggplot2::geom_boxplot(width=.1,outlier.shape=NA),out,"Fig3F_GOBP_FERROPTOSIS_violin",12,4)
write_tsv(data.frame(cell=colnames(Stattic),condition=Stattic$condition,cell_type=Stattic$CellType,
  Stattic[[]][,paste0(ids,"_score1"),drop=FALSE]),file.path(out,"atlas_module_scores.tsv"))
write_checkpoint(Stattic,"Stattic_scored.rds")

# ===========================================================================
# >>>>>>>>>>>>>>>>>>>>>>>  Figures S3B-C  <<<<<<<<<<<<<<<<<<<<<<<<<
# ===========================================================================

markers<-config$mouse$marker_panels
average_heatmap(Stattic,markers$gene[markers$panel=="atlas"],"CellType",CELL_TYPES,out,"FigS3B_marker_heatmap")
proportion_exports(Stattic,"CellType",config,out,"FigS3C")

finish_stage(out)

###############################################################################
# End of script 02.
###############################################################################
