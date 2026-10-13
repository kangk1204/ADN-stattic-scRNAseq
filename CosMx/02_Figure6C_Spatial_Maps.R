################################################################################
# 02_Figure6C_Spatial_Maps
# Source-defined biopsy ROI segmentation and transcript molecule maps
# Current manuscript: 2026-10-12. See the repository-root README.md.
# Run this R script as R code in a fresh R session.
################################################################################

####################  1. Libraries  ############################################
library(Seurat)
library(SeuratObject)
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
# Crop coordinates from the original Figure 6C analysis (Seurat plot coordinates).
figure6c_regions <- utils::read.delim(text = "roi_id	image	x_min	x_max	y_min	y_max	molecule_size
biopsy1_overview	biopsy_1	-70000	-63000	-448000	-443000	0.5
biopsy1_detail	biopsy_1	-64500	-63000	-446500	-445500	1
",
    check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
molecule_limit<-1000L # explicit source ImageDimPlot default, not every transcript
allow_other_seurat_version<-FALSE

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

NICHE_COLORS<-c("PT"="#F45F5A","Tubular injury"="#8B2424","Fibrotic"="#3F0099",
  "Glomerular"="#17B682","Endothelia"="#17A6E6","CD"="#D1B2FF","LOH"="#F740CE")
CELL_COLORS<-c("PT"="#F39C12","PT Injured"="#670000","PT Inflammatory"="#99004C",
 "LOH-DCT"="#EAEAEA","LOH-DCT Injured"="#EAEAEA","LOH-DCT Inflammatory"="#EAEAEA",
 "CD Injured"="#EAEAEA","PEC"="#EAEAEA","Podocyte"="#EAEAEA","PC"="#EAEAEA","IC"="#EAEAEA",
 "SMC"="#EAEAEA","Mesangial Cell"="#EAEAEA","Endothelia"="#17AEBF","Endothelia Glomerular"="#2ECC71",
 "Fibroblast"="#8E6A3E","Myofibroblast"="#8E44AD","Monocyte"="#E84393","Macrophage"="#FF00DD",
 "cDC"="#EAEAEA","Mast Cell"="#EAEAEA","NK"="#EAEAEA","T Cell"="#EAEAEA","B Cell"="#EAEAEA","Plasma Cell"="#EAEAEA")
validate_rois<-function(r,cfg) {
 check_columns(r,c("roi_id","image","x_min","x_max","y_min","y_max","molecule_size"),"ROI table")
 check_ids(r$roi_id,"ROI IDs")
 if(!nrow(r)||any(!grepl("^[A-Za-z][A-Za-z0-9_]*$",r$roi_id))||any(r$roi_id%in%cfg$image_names))stop("Invalid ROI IDs.")
 if(any(!r$image%in%cfg$image_names))stop("ROI image outside selected biopsies.")
 for(nm in c("x_min","x_max","y_min","y_max","molecule_size"))
  if(!is.numeric(r[[nm]])||any(!is.finite(r[[nm]])))stop("Invalid ROI coordinate: ",nm)
 if(any(r$x_max<=r$x_min)|any(r$y_max<=r$y_min)|any(r$molecule_size<=0))stop("Invalid ROI bounds/size.")
 invisible(TRUE)
}

####################  4. Analysis  #############################################

####################  Figure 6C: fixed source-defined representative regions  ###
if(as.character(utils::packageVersion("Seurat"))!="5.3.0" && !allow_other_seurat_version)
  stop("Reference Figure 6C plotting uses Seurat 5.3.0; set allow_other_seurat_version only after reviewing differences.")
roi<-figure6c_regions; validate_rois(roi,cfg)
config$figure6c_regions <- roi
config$molecule_limit <- molecule_limit
checkpoint_file<-input_checkpoint(config,"01","CosMx_9biopsies_cells.rds")
out<-init_stage(config,"02",c(cfg$rds_file,checkpoint_file))
checkpoint<-readRDS(checkpoint_file)
if(!identical(checkpoint$settings,extraction_settings(cfg)))stop("01 extraction settings differ; use matching settings.")
if(!identical(checkpoint$source_sha256,read_tsv(file.path(out,"input_sha256.tsv"))$sha256[1]))stop("The original RDS changed after stage 01.")
d<-checkpoint$cells;validate_cells(d,cfg)
obj<-readRDS(cfg$rds_file)
if(!all(c("Niche","Annotation.Lvl2")%in%names(obj@meta.data)))stop("Niche/cell-type annotations required for Figure 6C.")
# Preserve author annotations; no reclustering or renormalization.
for(field in c("Niche","Annotation.Lvl2")) {
  observed<-as.character(obj@meta.data[d$cell,field]);palette<-if(field=="Niche")NICHE_COLORS else CELL_COLORS
  unknown<-setdiff(unique(observed[!is.na(observed)]),names(palette))
  if(length(unknown))stop("Update explicit color map for unexpected ",field,": ",paste(unknown,collapse=", "))
}
roi_manifest<-list()
for(i in seq_len(nrow(roi))) {
  r<-roi[i,,drop=FALSE];name<-r$roi_id;image<-r$image
  # Preserve the original Crop default coordinate convention explicitly.
  crop<-SeuratObject::Crop(obj[[image]],x=c(r$x_min,r$x_max),y=c(r$y_min,r$y_max),coords="plot")
  if(!"segmentation"%in%names(attr(crop,"boundaries")))stop("Segmentation boundary absent: ",name)
  crop<-SeuratObject::`DefaultBoundary<-`(crop,value="segmentation")
  if(!length(SeuratObject::Molecules(crop)))stop("Molecule coordinates absent: ",name)
  if(!all(c("SPP1","THBS1","STAT3")%in%SeuratObject::Features(crop,set=NA)))stop("Required transcript molecules absent: ",name)
  # Retain exactly the selected biopsy's cells; never draw another donor's cells.
  selected_ids<-intersect(SeuratObject::Cells(crop),d$cell)
  if(!length(selected_ids))stop("No retained cells in cropped FOV: ",name)
  actual<-d[match(selected_ids,d$cell),,drop=FALSE]
  if(length(unique(actual$Specimen_ID))!=1L)stop("Cropped segmentation crosses specimen boundaries: ",name)
  # Seurat 5.3.0 ImageDimPlot splits internal FOV names at underscores when
  # attaching molecule layers. Keep descriptive output IDs, but use a safe
  # internal FOV name so transcript overlays are not silently omitted.
  plot_fov<-paste0("roi",i)
  obj[[plot_fov]]<-crop
  write.csv(actual[,c("cell","image","Specimen_ID","FOV","Niche","Annotation.Lvl2")],
    file.path(out,paste0(name,"_displayed_cells.csv")),row.names=FALSE)
  # Store cropped spatial data for exact inspection of available molecules/boundaries.
  saveRDS(crop,file.path(out,paste0(name,"_cropped_FOV.rds")))
  draw <- function(field,genes=NULL,colors=NULL,mol_size=.5) {
    set.seed(cfg$seed) # fixed sampled molecule subset for each plot
    p<-Seurat::ImageDimPlot(obj,fov=plot_fov,boundaries="segmentation",group.by=field,
      cells=selected_ids,cols=if(field=="Niche")NICHE_COLORS else CELL_COLORS,
      molecules=genes,mols.cols=colors,mols.size=mol_size,nmols=molecule_limit,
      size=1,axes=TRUE,dark.background=FALSE,flip_xy=TRUE,combine=TRUE)+
      ggplot2::theme_bw()+ggplot2::theme(panel.grid=ggplot2::element_blank(),
        legend.text=ggplot2::element_text(color="grey10",size=9))
    if(length(genes)) {
      # Validate the actual plotting layers, not merely the gene names stored
      # in the FOV. A successful PDF save alone cannot verify an overlay.
      first_plot<-if(inherits(p,"patchwork"))p[[1L]]else p
      molecule_layers<-lapply(first_plot$layers,function(layer) {
        z<-layer$data
        if(is.data.frame(z)&&"molecule"%in%names(z))z else NULL
      })
      displayed<-do.call(rbind,molecule_layers)
      if(is.null(displayed)||!nrow(displayed)||
         !all(genes%in%as.character(displayed$molecule)))
        stop("Transcript overlay missing from actual plot: ",name," / ",field)
      write.csv(displayed,file.path(out,paste0(name,"_",field,"_",
        paste(genes,collapse="_"),"_displayed_molecules.csv")),row.names=FALSE)
    }
    p
  }
  save_panel(draw("Niche"),out,paste0("Figure6C_",name,"_Niche"),8,8)
  save_panel(draw("Niche","SPP1",c(SPP1="#FFFF00")),out,paste0("Figure6C_",name,"_Niche_SPP1"),8,8)
  save_panel(draw("Annotation.Lvl2"),out,paste0("Figure6C_",name,"_CellTypes"),8,8)
  save_panel(draw("Annotation.Lvl2",rev(c("SPP1","THBS1","STAT3")),
    c(SPP1="yellow",THBS1="black",STAT3="#00A500"),r$molecule_size),out,
    paste0("Figure6C_",name,"_CellTypes_SPP1_THBS1_STAT3"),8,8)
  roi_manifest[[name]]<-data.frame(r,specimen=unique(actual$Specimen_ID),
    n_displayed_cells=nrow(actual),nmols=molecule_limit,seed=cfg$seed,coords="plot",flip_xy=TRUE)
  obj[[plot_fov]]<-NULL
}
write.csv(do.call(rbind,roi_manifest),file.path(out,"Figure6C_ROI_manifest.csv"),row.names=FALSE)
writeLines(c("The two biopsy_1 crops reproduce source-script ROI coordinates, not a new significance-based ROI search.",
  "Verify final manuscript crop/composition against the supplied source regions.",
  "ImageDimPlot nmols=1000 is the source default; sampling is now explicitly seeded per plot.",
  "Transcript molecule overlays are not the SCT positivity gate used in Figure 6D.",
  "Regional co-occurrence does not itself prove single-cell co-expression or signaling."),file.path(out,"analysis_notes.txt"))
finish_stage(out)
