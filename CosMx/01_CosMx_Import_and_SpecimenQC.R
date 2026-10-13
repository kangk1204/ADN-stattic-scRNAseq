################################################################################
# 01_CosMx_Import_and_SpecimenQC
# Import the processed GEO object and audit all nine biopsy specimens
# Current manuscript: 2026-10-12. See the repository-root README.md.
# Run this R script as R code in a fresh R session.
################################################################################

####################  1. Libraries  ############################################
library(Matrix)
library(ggplot2)
library(patchwork)
library(digest)

####################  2. Configuration  ########################################
locate_project <- function ()
{
    for (i in rev(seq_along(sys.frames()))) {
        f <- sys.frame(i)$ofile
        if (!is.null(f))
            return(dirname(normalizePath(f)))
    }
    arg <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(arg)) {
        f <- sub("^--file=", "", arg[1])
        if (!file.exists(f))
            f <- gsub("~+~", " ", f, fixed = TRUE)
        return(dirname(normalizePath(f, mustWork = TRUE)))
    }
    normalizePath(getwd())
}


project_dir <- locate_project()
rds_file <- Sys.getenv("KIDNEY_COSMX_RDS",file.path(project_dir,"data","GSE253439_seurat_object_cosmx.rds"))
results_root <- Sys.getenv("KIDNEY_RESULTS_DIR", file.path(project_dir, "results"))
restart_failed_stage <- FALSE # TRUE archives only an incomplete run; completed runs stay protected.
cfg <- list(
  rds_file=rds_file, image_names=c("biopsy_1","biopsy_2","biopsy_3"),
  specimen_ids=paste0("Biopsy",1:9), expression_assay="SCT", expression_layer="data",
  expression_threshold=0, um_per_pixel=0.12,
  min_source_distance=1L, min_target_distance=1L, sparse_target_flag_below=5L,
  min_source_radial=5L, min_target_radial=5L,
  radius_um=seq(0,50,by=0.5), radial_correction="isotropic",
  observation_windows_rds=NULL, fov_rectangle_padding_um=0.12,
  save_cell_violin=TRUE, seed=1234L)
config <- list(.root=project_dir,output_dir=results_root,seed=cfg$seed,cosmx=cfg)
SOURCE_LEVELS <- c("PT","SPP1+ PT Injured")
TARGET_LEVELS <- c("STAT3+THBS1+ Macrophage","Myofibroblast")
ALL_LEVELS <- c(SOURCE_LEVELS,TARGET_LEVELS)
SOURCE_COLORS <- c("PT"="#FFA648","SPP1+ PT Injured"="#FF5E00")

####################  3. Reproducibility and analysis helpers  ###################
need <- function (packages)
{
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing))
        stop("Missing packages: ", paste(missing, collapse = ", "), ". See README.md.", call. = FALSE)
}

path_at <- function (root, x)
{
    if (is.null(x))
        return(NULL)
    if (grepl("^(~|/|[A-Za-z]:[/\\\\])", x))
        path.expand(x)
    else file.path(root, x)
}

read_tsv <- function (file)
{
    if (!file.exists(file))
        stop("Missing input: ", file, call. = FALSE)
    read.delim(file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
}

write_tsv <- function (x, path)
{
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

check_columns <- function (x, required, label)
{
    absent <- setdiff(required, names(x))
    if (length(absent))
        stop(label, " lacks columns: ", paste(absent, collapse = ", "))
}

check_ids <- function (x, label)
{
    if (anyNA(x) || any(!nzchar(trimws(x))) || anyDuplicated(x))
        stop(label, " must be nonempty, nonmissing, unique IDs.")
}

file_hash <- function (x)
{
    need("digest")
    if (!file.exists(x))
        stop("Cannot hash missing input: ", x)
    digest::digest(file = x, algo = "sha256")
}

hash_manifest <- function (paths, out, relative_to = NULL)
{
    paths <- unique(paths)
    if (anyNA(paths) || any(!nzchar(paths)))
        stop("Invalid path in required input manifest.")
    missing <- paths[!file.exists(paths) | dir.exists(paths)]
    if (length(missing))
        stop("Required files missing or not regular files: ", paste(missing, collapse = ", "))
    labels <- normalizePath(paths, mustWork = TRUE)
    if (!is.null(relative_to))
        labels <- substring(labels, nchar(normalizePath(relative_to)) + 2L)
    write_tsv(data.frame(file = labels, bytes = file.info(paths)$size, sha256 = vapply(paths, file_hash, character(1))), out)
}

stage_dir <- function (config, stage)
file.path(path_at(config$.root, config$output_dir), stage)

init_stage <- function (config, stage, inputs = character())
{
    out <- stage_dir(config, stage)
    if (dir.exists(out) && !file.exists(file.path(out, "SUCCESS")) && isTRUE(restart_failed_stage)) {
        archived <- paste0(out, "_failed_", format(Sys.time(), "%Y%m%d_%H%M%S"))
        if (file.exists(archived) || !file.rename(out, archived))
            stop("Could not archive failed stage.")
    }
    if (dir.exists(out) || file.exists(out))
        stop("Stage directory already exists: ", out, ". Choose a new output_dir; failed-run files must not mix with a new run.")
    missing <- inputs[!file.exists(inputs) | dir.exists(inputs)]
    if (anyNA(inputs) || length(missing))
        stop("Required stage inputs missing: ", paste(missing, collapse = ", "))
    dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
    if (!dir.create(out, showWarnings = FALSE))
        stop("Could not claim new stage directory: ", out)
    writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_start.txt"))
    dput(config, file.path(out, "run_config.R"))
    hash_manifest(inputs, file.path(out, "input_sha256.tsv"))
    code <- c(list.files(config$.root, full.names = TRUE, pattern = "\\.(R|txt)$"),
                  list.files(file.path(config$.root, "gene_sets"), full.names = TRUE))
    code <- code[!dir.exists(code)]
    hash_manifest(code, file.path(out, "code_sha256.tsv"), relative_to = config$.root)
    set.seed(config$seed)
    out
}

finish_stage <- function (out)
{
    writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
    loaded <- sort(loadedNamespaces())
    write_tsv(data.frame(package = loaded, version = vapply(loaded, function(p) as.character(getNamespaceVersion(p)), character(1)), path = vapply(loaded,
        function(p) if (p == "base")
            file.path(R.home(), "library/base")
        else getNamespaceInfo(asNamespace(p), "path"), character(1))), file.path(out, "loaded_packages.tsv"))
    files <- list.files(out, full.names = TRUE, recursive = TRUE)
    files <- files[!basename(files) %in% c("SUCCESS", "output_sha256.tsv")]
    hash_manifest(files, file.path(out, "output_sha256.tsv"), relative_to = out)
    writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), file.path(out, "SUCCESS"))
}


extraction_settings <- function(cfg) cfg[c("image_names","specimen_ids","expression_assay",
  "expression_layer","expression_threshold","um_per_pixel")]
input_checkpoint <- function(config,stage,file) {
  dir<-stage_dir(config,stage);p<-file.path(dir,file)
  if(!file.exists(file.path(dir,"SUCCESS")))stop("Upstream stage incomplete: ",stage)
  tab<-read_tsv(file.path(dir,"output_sha256.tsv"))
  hit<-which(tab$file==file)
  if(length(hit)!=1L||!file.exists(p)||file_hash(p)!=tab$sha256[hit])stop("Upstream checkpoint missing or changed: ",p)
  p
}
validate_cells <- function(d,cfg) {
  check_columns(d,c("cell","image","Specimen_ID","FOV","Annotation.Lvl2","x_px","y_px",
    "x_um","y_um","SPP1","STAT3","THBS1","fov_key","label"),"Extracted cells")
  check_ids(d$cell,"Spatial cell IDs")
  if(!setequal(unique(d$Specimen_ID),cfg$specimen_ids))stop("Specimen selection differs from cfg.")
  if(anyNA(d[,c("image","Specimen_ID","FOV","x_px","y_px","x_um","y_um")]))stop("Missing spatial identifiers/coordinates.")
  for(nm in c("x_px","y_px","x_um","y_um","SPP1","STAT3","THBS1"))
    if(!is.numeric(d[[nm]])||any(!is.finite(d[[nm]])))stop("Invalid numeric column: ",nm)
  if(any(abs(d$x_um-d$x_px*cfg$um_per_pixel)>1e-8)||any(abs(d$y_um-d$y_px*cfg$um_per_pixel)>1e-8))stop("Coordinate conversion differs from cfg.")
  if(!all(d$image%in%cfg$image_names))stop("Unknown image.")
  if(any(vapply(split(d$image,d$Specimen_ID),function(x)length(unique(x)),integer(1))!=1L))stop("One specimen spans multiple image coordinate frames; inspect before analysis.")
  if(any(!is.na(d$label)&!d$label%in%ALL_LEVELS))stop("Unknown analysis label.")
  expected<-rep(NA_character_,nrow(d));ct<-as.character(d$Annotation.Lvl2);thr<-cfg$expression_threshold
  expected[which(ct=="PT")]<-SOURCE_LEVELS[1]
  expected[which(ct=="PT Injured"&d$SPP1>thr)]<-SOURCE_LEVELS[2]
  expected[which(ct=="Macrophage"&d$STAT3>thr&d$THBS1>thr)]<-TARGET_LEVELS[1]
  expected[which(ct=="Myofibroblast")]<-TARGET_LEVELS[2]
  if(!identical(as.character(d$label),expected))stop("Analysis labels do not match expression/annotation rules.")
  if(!identical(as.character(d$fov_key),paste(d$image,d$Specimen_ID,d$FOV,sep="::")))stop("FOV keys mismatch.")
  invisible(TRUE)
}
save_panel <- function(p,out,name,width=8,height=6) {
  ggplot2::ggsave(file.path(out,paste0(name,".pdf")),plot=p,width=width,height=height,
    device=grDevices::pdf,useDingbats=FALSE)
  if(requireNamespace("ragg",quietly=TRUE))ggplot2::ggsave(file.path(out,paste0(name,".png")),
    plot=p,width=width,height=height,dpi=300,device=ragg::agg_png)
}
validate_cosmx_config <- function (config)
{
    check_ids(config$specimen_ids, "CosMx specimen IDs")
    check_ids(config$image_names, "CosMx image names")
    if (!length(config$specimen_ids) || !length(config$image_names))
        stop("CosMx selections cannot be empty.")
    if (!identical(config$radial_correction, "isotropic"))
        stop("This reference engine supports isotropic correction only.")
    r <- config$radius_um
    if (!is.numeric(r) || length(r) < 2L || anyNA(r) || any(!is.finite(r)) || r[1] != 0 || any(diff(r) <= 0) || max(abs(diff(r) - diff(r)[1])) >
        1e-08)
        stop("CosMx radii must be finite, equally spaced, increasing, and start at zero.")
    for (nm in c("um_per_pixel", "min_source_distance", "min_target_distance", "min_source_radial", "min_target_radial", "fov_rectangle_padding_um")) {
        v <- config[[nm]]
        if (length(v) != 1L || !is.numeric(v) || !is.finite(v) || v <= 0)
            stop("Invalid CosMx parameter: ", nm)
    }
    if (length(config$expression_threshold)!=1L || !is.numeric(config$expression_threshold) || !is.finite(config$expression_threshold)) stop("Invalid expression threshold.")
    for(nm in c("min_source_distance","min_target_distance","min_source_radial","min_target_radial","sparse_target_flag_below"))
      if(length(config[[nm]])!=1L||!is.finite(config[[nm]])||config[[nm]]<1||config[[nm]]!=as.integer(config[[nm]]))stop("Cell count threshold must be a positive integer: ",nm)
    invisible(TRUE)
}

extract_cells <- function (config)
{
    message("Reading RDS: ", config$rds_file)
    if (!file.exists(config$rds_file))
        stop("RDS not found. Check cfg$rds_file.")
    obj <- readRDS(config$rds_file)
    meta <- attr(obj, "meta.data")
    images <- attr(obj, "images")
    assays <- attr(obj, "assays")
    if (is.null(meta) || !is.data.frame(meta)) stop("Missing metadata table.")
    check_ids(rownames(meta), "Metadata cell IDs")
    required <- c("Specimen_ID", "FOV", "Annotation.Lvl2")
    if (!all(required %in% names(meta)))
        stop("Required metadata absent: ", paste(setdiff(required, names(meta)), collapse = ", "))
    if (!all(config$image_names %in% names(images)))
        stop("Requested spatial image absent.")
    aa <- assays[[config$expression_assay]]
    if (is.null(aa))
        stop("Assay absent: ", config$expression_assay)
    mat <- attr(aa, config$expression_layer)
    if (is.null(mat)) {
        if (!requireNamespace("SeuratObject", quietly = TRUE))
            stop("Assay5 input needs SeuratObject; this is not needed for the supplied Assay/SCTAssay RDS.")
        mat <- SeuratObject::LayerData(aa, layer = config$expression_layer)
    }
    check_ids(rownames(mat), "Assay gene IDs")
    check_ids(colnames(mat), "Assay cell IDs")
    genes <- c("SPP1", "STAT3", "THBS1")
    if (!all(genes %in% rownames(mat)))
        stop("A required gene is missing in the selected assay/layer.")
    coords <- do.call(rbind, lapply(config$image_names, function(im) {
        cc <- attr(images[[im]], "boundaries")[["centroids"]]
        xy <- attr(cc, "coords")
        ids <- attr(cc, "cells")
        if (is.null(xy) || ncol(xy) != 2L || nrow(xy) != length(ids))
            stop("Unexpected centroid structure: ", im)
        data.frame(cell = ids, image = im, x_px = xy[, 1], y_px = xy[, 2], row.names = NULL)
    }))
    if (anyDuplicated(coords$cell))
        stop("Duplicate cell IDs across selected images.")
    ix <- match(coords$cell, rownames(meta))
    if (anyNA(ix))
        stop("Coordinates contain cells absent from metadata.")
    keep <- intersect(c(required, "Slide_ID", "Sensor_ID", "Niche"), names(meta))
    d <- cbind(coords, meta[ix, keep, drop = FALSE])
    rownames(d) <- NULL
    for (v in intersect(c("Specimen_ID", "Annotation.Lvl2", "Slide_ID", "Sensor_ID", "Niche"), names(d))) d[[v]] <- as.character(d[[v]])
    d <- d[d$Specimen_ID %in% config$specimen_ids, , drop = FALSE]
    if (!all(config$specimen_ids %in% d$Specimen_ID))
        stop("Requested specimen absent: ", paste(setdiff(config$specimen_ids, d$Specimen_ID), collapse = ", "))
    if (anyNA(d[, c("Specimen_ID", "FOV", "x_px", "y_px")]))
        stop("Missing coordinates or specimen/FOV IDs.")
    if (any(!is.finite(d$x_px)) || any(!is.finite(d$y_px)))
        stop("Non-finite coordinates.")
    ix <- match(d$cell, colnames(mat))
    if (anyNA(ix))
        stop("Selected assay lacks some spatial cells.")
    expr <- t(as.matrix(mat[genes, ix, drop = FALSE]))
    if (anyNA(expr) || any(!is.finite(expr)))
        stop("Missing or non-finite marker expression.")
    d <- cbind(d, expr)
    d$x_um <- d$x_px * config$um_per_pixel
    d$y_um <- d$y_px * config$um_per_pixel
    d$fov_key <- paste(d$image, d$Specimen_ID, d$FOV, sep = "::")
    d$label <- NA_character_
    ct <- d$Annotation.Lvl2
    thr <- config$expression_threshold
    d$label[which(ct == "PT")] <- SOURCE_LEVELS[1]
    d$label[which(ct == "PT Injured" & d$SPP1 > thr)] <- SOURCE_LEVELS[2]
    d$label[which(ct == "Macrophage" & d$STAT3 > thr & d$THBS1 > thr)] <- TARGET_LEVELS[1]
    d$label[which(ct == "Myofibroblast")] <- TARGET_LEVELS[2]
    message("Extracted ", nrow(d), " cells from ", length(unique(d$Specimen_ID)), " specimens; marker criterion: ", config$expression_assay,
        "/", config$expression_layer, " > ", thr)
    rm(obj, assays, images, meta, mat, aa, expr)
    invisible(gc())
    d
}

cell_manifest <- function (d, config)
{
    do.call(rbind, lapply(config$specimen_ids, function(sp) {
        z <- d[d$Specimen_ID == sp, , drop = FALSE]
        data.frame(specimen = sp, image = paste(unique(z$image), collapse = ";"), n_all_cells = nrow(z), n_FOV = length(unique(z$fov_key)), n_PT = sum(z$label %in%
            SOURCE_LEVELS[1]), n_SPP1_PTinj = sum(z$label %in% SOURCE_LEVELS[2]), n_STAT3_THBS1_macrophage = sum(z$label %in% TARGET_LEVELS[1]),
            n_myofibroblast = sum(z$label %in% TARGET_LEVELS[2]), expression_assay = config$expression_assay, expression_layer = config$expression_layer)
    }))
}


####################  4. Analysis  #############################################

####################  Figure 6C/D: input and specimen audit  ####################
validate_cosmx_config(cfg)
out <- init_stage(config,"01",cfg$rds_file)
d <- extract_cells(cfg)
validate_cells(d,cfg)
manifest <- cell_manifest(d,cfg)
write.csv(manifest,file.path(out,"specimen_cell_counts.csv"),row.names=FALSE)
write.csv(as.data.frame(table(d$image,d$Specimen_ID)),file.path(out,"image_specimen_counts.csv"),row.names=FALSE)
write.csv(as.data.frame(table(d$Specimen_ID,d$Annotation.Lvl2)),file.path(out,"annotation_counts.csv"),row.names=FALSE)
if("Niche"%in%names(d))write.csv(as.data.frame(table(d$Specimen_ID,d$Niche)),file.path(out,"niche_counts.csv"),row.names=FALSE)
checkpoint <- list(schema_version=1L,settings=extraction_settings(cfg),
  source_sha256=read_tsv(file.path(out,"input_sha256.tsv"))$sha256[1],cells=d)
saveRDS(checkpoint,file.path(out,"CosMx_9biopsies_cells.rds"),compress="gzip")
# These are reference audit targets from the manuscript, never imposed counts.
audit <- data.frame(metric=c("specimens","cells","FOVs"),reference=c(9,180323,116),
  observed=c(length(unique(d$Specimen_ID)),nrow(d),length(unique(d$fov_key))))
audit$matches_reference<-audit$reference==audit$observed
write.csv(audit,file.path(out,"manuscript_count_comparison.csv"),row.names=FALSE)
if(!all(audit$matches_reference))warning("Counts differ from the study reference; input and selection details are recorded in the audit outputs.")
####################  Figure 6C: full-specimen centroid context  ################
cols<-c(SOURCE_COLORS,"STAT3+THBS1+ Macrophage"="#176F38","Myofibroblast"="#6A3D9A")
# Separate coordinate frames with a fixed aspect ratio in each panel.
plots<-lapply(cfg$specimen_ids,function(sp) {
  q<-d[d$Specimen_ID==sp,,drop=FALSE]
  ggplot(q,aes(x_um,y_um,color=label))+geom_point(size=.1)+coord_equal()+theme_void()+
    scale_color_manual(values=cols,na.value="grey90",limits=names(cols),drop=FALSE)+labs(title=sp,color=NULL)
})
save_panel(patchwork::wrap_plots(plots,ncol=3,guides="collect"),out,"Figure6C_specimen_centroid_context",12,12)
writeLines(c("Centroid context maps are not transcript-molecule or cell-boundary maps.",
  "Use script 02 for source-defined cropped segmentation/molecule panels.",
  "One image (e.g. biopsy_1) may contain more than one Specimen_ID.",
  "Marker positivity: SCT/data > 0 by default; not raw transcript counts."),file.path(out,"analysis_notes.txt"))
finish_stage(out)
