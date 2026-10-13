################################################################################
# 03_Figure6D_Specimenwise_Neighborhood
# Nine-specimen nearest-target distances and radial enrichment
# Current manuscript: 2026-10-12. See the repository-root README.md.
# Run this R script as R code in a fresh R session.
################################################################################

####################  1. Libraries  ############################################
library(Matrix)
library(RANN)
library(spatstat.geom)
library(spatstat.explore)
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

make_windows <- function (d, config, run_dir)
{
    supplied <- !is.null(config$observation_windows_rds)
    if (supplied) {
        windows <- readRDS(config$observation_windows_rds)
        if (!is.list(windows) || !all(config$specimen_ids %in% names(windows)))
            stop("Window RDS must be a named list covering all selected Specimen_IDs.")
        mode <- "user_supplied_window"
    }
    else {
        windows <- list()
        mode <- "estimated_union_of_all_cell_FOV_rectangles"
        message("Radial windows are estimated from all retained cells per FOV; window_QC.pdf displays the estimated observation windows.")
    }
    fov_rows <- list()
    window_rows <- list()
    for (sp in config$specimen_ids) {
        z <- d[d$Specimen_ID == sp, , drop = FALSE]
        fovs <- split(z, z$fov_key)
        rectangles <- lapply(fovs, function(f) {
            pad <- config$fov_rectangle_padding_um
            xr <- range(f$x_um) + c(-pad, pad)
            yr <- range(f$y_um) + c(-pad, pad)
            fov_rows[[f$fov_key[1]]] <<- data.frame(specimen = sp, fov_key = f$fov_key[1], n_all_cells = nrow(f), xmin_um = xr[1], xmax_um = xr[2],
                ymin_um = yr[1], ymax_um = yr[2])
            spatstat.geom::owin(xrange = xr, yrange = yr)
        })
        if (!supplied)
            windows[[sp]] <- Reduce(spatstat.geom::union.owin, rectangles)
        W <- windows[[sp]]
        if (!inherits(W, "owin") || !is.finite(spatstat.geom::area.owin(W)) || spatstat.geom::area.owin(W) <= 0)
            stop("Invalid observation window: ", sp)
        inside <- spatstat.geom::inside.owin(z$x_um, z$y_um, W)
        if (!all(inside))
            stop(sp, ": ", sum(!inside), " cells outside window. Check units and masks; cells are not silently dropped.")
        window_rows[[sp]] <- data.frame(specimen = sp, window_mode = mode, area_um2 = spatstat.geom::area.owin(W), n_FOV = length(fovs), n_all_cells = nrow(z))
    }
    saveRDS(windows, file.path(run_dir, "observation_windows_um.rds"))
    write.csv(do.call(rbind, fov_rows), file.path(run_dir, "FOV_coordinate_extents.csv"), row.names = FALSE)
    write.csv(do.call(rbind, window_rows), file.path(run_dir, "observation_window_manifest.csv"), row.names = FALSE)
    grDevices::pdf(file.path(run_dir, "window_QC.pdf"), width = 9, height = 7, onefile = TRUE)
    for (sp in config$specimen_ids) {
        z <- d[d$Specimen_ID == sp, , drop = FALSE]
        plot(windows[[sp]], main = paste(sp, "-", mode), xlab = "x (um)", ylab = "y (um)")
        graphics::points(z$x_um, z$y_um, pch = 16, cex = 0.15, col = grDevices::adjustcolor("grey40", 0.35))
        for (lab in ALL_LEVELS) {
            q <- z[z$label %in% lab, , drop = FALSE]
            colors <- c(SOURCE_COLORS, `STAT3+THBS1+ Macrophage` = "#176F38", Myofibroblast = "#6A3D9A")
            graphics::points(q$x_um, q$y_um, pch = 16, cex = if (lab == TARGET_LEVELS[1])
                0.8
            else 0.2, col = colors[[lab]])
        }
    }
    grDevices::dev.off()
    if (requireNamespace("ragg", quietly = TRUE)) {
        ragg::agg_png(file.path(run_dir, "window_QC_overview.png"), width = 1800, height = 1800, res = 150)
        graphics::par(mfrow = c(3, 3), mar = c(3, 3, 3, 1))
        for (sp in config$specimen_ids) {
            z <- d[d$Specimen_ID == sp, , drop = FALSE]
            plot(windows[[sp]], main = paste(sp, "|", length(unique(z$fov_key)), "FOVs"), xlab = "x (um)", ylab = "y (um)")
            graphics::points(z$x_um, z$y_um, pch = 16, cex = 0.12, col = grDevices::adjustcolor("#3A647A", 0.3))
        }
        grDevices::dev.off()
    }
    list(windows = windows, mode = mode)
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

distance_analysis <- function (d, config)
{
    rows <- list()
    eligibility <- list()
    summaries <- list()
    for (sp in config$specimen_ids) for (tar in TARGET_LEVELS) {
        z <- d[d$Specimen_ID == sp, , drop = FALSE]
        target <- z[z$label %in% tar, , drop = FALSE]
        n_src <- sapply(SOURCE_LEVELS, function(lab) sum(z$label %in% lab))
        eligible <- all(n_src >= config$min_source_distance) && nrow(target) >= config$min_target_distance
        reason <- if (eligible)
            "included"
        else if (nrow(target) < config$min_target_distance)
            "no_or_too_few_targets"
        else "missing_source_group"
        eligibility[[paste(sp, tar)]] <- data.frame(specimen = sp, target = tar, n_PT = n_src[1], n_SPP1_PTinj = n_src[2], n_target = nrow(target),
            eligible = eligible, status = reason, sparse_target = nrow(target) < config$sparse_target_flag_below)
        for (src in SOURCE_LEVELS) {
            q <- z[z$label %in% src, , drop = FALSE]
            if (!nrow(q))
                next
            rr <- data.frame(cell = q$cell, specimen = sp, source = src, target = tar, n_target = nrow(target), target_cell = NA_character_,
                target_specimen = NA_character_, distance_um = NA_real_, eligible_pair = eligible)
            if (nrow(target) >= config$min_target_distance) {
                near <- RANN::nn2(data = as.matrix(target[, c("x_um", "y_um")]), query = as.matrix(q[, c("x_um", "y_um")]), k = 1L)
                rr$distance_um <- near$nn.dists[, 1]
                rr$target_cell <- target$cell[near$nn.idx[, 1]]
                rr$target_specimen <- target$Specimen_ID[near$nn.idx[, 1]]
                stopifnot(all(rr$specimen == rr$target_specimen))
            }
            rows[[paste(sp, tar, src)]] <- rr
            summaries[[paste(sp, tar, src)]] <- data.frame(specimen = sp, source = src, target = tar, n_source = nrow(q), n_target = nrow(target),
                eligible_pair = eligible, median_distance_um = if (eligible)
                  median(rr$distance_um)
                else NA_real_, mean_distance_um = if (eligible)
                  mean(rr$distance_um)
                else NA_real_)
        }
    }
    cells <- if (length(rows))
        do.call(rbind, rows)
    else data.frame(cell = character(), specimen = character(), source = character(), target = character(), n_target = integer(), target_cell = character(),
        target_specimen = character(), distance_um = numeric(), eligible_pair = logical())
    summary <- if (length(summaries))
        do.call(rbind, summaries)
    else data.frame(specimen = character(), source = character(), target = character(), n_source = integer(), n_target = integer(), eligible_pair = logical(),
        median_distance_um = numeric(), mean_distance_um = numeric())
    tests <- list()
    per_specimen <- list()
    for (tar in TARGET_LEVELS) {
        ss <- summary[summary$target == tar & summary$eligible_pair, , drop = FALSE]
        a <- ss[ss$source == SOURCE_LEVELS[1], c("specimen", "median_distance_um")]
        b <- ss[ss$source == SOURCE_LEVELS[2], c("specimen", "median_distance_um")]
        names(a)[2] <- "median_PT_um"
        names(b)[2] <- "median_PTinj_um"
        paired <- merge(a, b, by = "specimen", all = FALSE)
        paired$target <- rep(tar, nrow(paired))
        paired$log2_PT_over_PTinj <- log2(paired$median_PT_um/paired$median_PTinj_um)
        per_specimen[[tar]] <- paired
        p <- NA_real_
        if (nrow(paired) >= 2L) {
            if (all(paired$median_PT_um == paired$median_PTinj_um))
                p <- 1
            else p <- stats::wilcox.test(paired$median_PTinj_um, paired$median_PT_um, paired = TRUE, alternative = "two.sided", exact = FALSE,
                correct = TRUE)$p.value
        }
        tests[[tar]] <- data.frame(target = tar, n_specimens = nrow(paired), n_specimens_PTinj_closer = sum(paired$median_PTinj_um < paired$median_PT_um),
            p_two_sided = p, test = "paired Wilcoxon signed-rank on specimen medians; exact=FALSE; correct=TRUE")
    }
    tests <- do.call(rbind, tests)
    tests$p_Holm <- stats::p.adjust(tests$p_two_sided, method = "holm")
    list(cells = cells, summaries = summary, eligibility = do.call(rbind, eligibility), paired = do.call(rbind, per_specimen), tests = tests)
}

radial_analysis <- function (d, windows, config)
{
    curves <- list()
    eligibility <- list()
    messages <- list()
    for (sp in config$specimen_ids) {
        z <- d[d$Specimen_ID == sp & !is.na(d$label), , drop = FALSE]
        pp <- spatstat.geom::ppp(z$x_um, z$y_um, window = windows[[sp]], marks = factor(z$label, levels = ALL_LEVELS), checkdup = TRUE)
        if (pp$n != nrow(z))
            stop("A point was dropped from the observation window: ", sp)
        for (tar in TARGET_LEVELS) {
            ni <- sapply(SOURCE_LEVELS, function(s) sum(z$label == s))
            nt <- sum(z$label == tar)
            eligible <- all(ni >= config$min_source_radial) && nt >= config$min_target_radial
            er <- data.frame(specimen = sp, target = tar, n_PT = ni[1], n_SPP1_PTinj = ni[2], n_target = nt, included = FALSE, reason = "too_few_source_or_target_cells")
            if (!eligible) {
                eligibility[[paste(sp, tar)]] <- er
                next
            }
            pair_curves <- list()
            ok <- TRUE
            for (src in SOURCE_LEVELS) {
                message("Radial: ", sp, "; ", src, " -> ", tar)
                warnings <- character()
                g <- tryCatch(withCallingHandlers(spatstat.explore::pcfcross(pp, i = src, j = tar, r = config$radius_um, correction = config$radial_correction),
                  warning = function(w) {
                    warnings <<- c(warnings, conditionMessage(w))
                    invokeRestart("muffleWarning")
                  }), error = function(e) e)
                if (length(warnings))
                  messages[[length(messages) + 1L]] <- data.frame(specimen = sp, source = src, target = tar, message = paste(warnings, collapse = " | "))
                if (inherits(g, "error")) {
                  ok <- FALSE
                  er$reason <- paste("pcfcross_error", conditionMessage(g))
                  break
                }
                v <- g$iso
                if (is.null(v))
                  stop("Isotropic g(r) column missing.")
                v[!is.finite(v) | v < 0] <- NA_real_
                lg <- rep(NA_real_, length(v))
                good <- !is.na(v) & v > 0
                lg[good] <- log2(v[good])
                pair_curves[[src]] <- data.frame(specimen = sp, source = src, target = tar, r_um = g$r, g_r = v, log2_g = lg)
            }
            if (ok) {
                valid <- is.finite(pair_curves[[1]]$g_r) & is.finite(pair_curves[[2]]$g_r) & config$radius_um > 0
                if (!any(valid)) {
                  er$reason <- "no_radius_with_both_curves_available"
                }
                else {
                  er$included <- TRUE
                  er$reason <- "included"
                  for (src in SOURCE_LEVELS) curves[[paste(sp, tar, src)]] <- pair_curves[[src]]
                }
            }
            eligibility[[paste(sp, tar)]] <- er
        }
    }
    if (!length(curves))
        return(list(curves = data.frame(specimen = character(), source = character(), target = character(), r_um = numeric(), g_r = numeric(),
            log2_g = numeric()), mean = data.frame(target = character(), r_um = numeric(), source = character(), n_specimens = integer(), mean_g_r = numeric(),
            log2_mean_g = numeric(), n_zero_estimates = integer()), contrast = data.frame(target = character(), r_um = numeric(), n_specimens = integer(),
            delta_log2_mean_g = numeric()), eligibility = do.call(rbind, eligibility), warnings = if (length(messages)) do.call(rbind, messages) else data.frame(specimen = character(),
            source = character(), target = character(), message = character())))
    curves <- do.call(rbind, curves)
    radial_mean <- list()
    contrast <- list()
    for (tar in TARGET_LEVELS) for (rr in config$radius_um[config$radius_um > 0]) {
        q <- curves[curves$target == tar & curves$r_um == rr, , drop = FALSE]
        a <- q[q$source == SOURCE_LEVELS[1], c("specimen", "g_r")]
        b <- q[q$source == SOURCE_LEVELS[2], c("specimen", "g_r")]
        names(a)[2] <- "g_PT"
        names(b)[2] <- "g_PTinj"
        pair <- merge(a, b, by = "specimen")
        pair <- pair[is.finite(pair$g_PT) & is.finite(pair$g_PTinj), , drop = FALSE]
        if (!nrow(pair))
            next
        av <- c(mean(pair$g_PT), mean(pair$g_PTinj))
        for (k in 1:2) radial_mean[[paste(tar, rr, k)]] <- data.frame(target = tar, r_um = rr, source = SOURCE_LEVELS[k], n_specimens = nrow(pair),
            mean_g_r = av[k], log2_mean_g = if (av[k] > 0)
                log2(av[k])
            else NA_real_, n_zero_estimates = sum(pair[[k + 1L]] == 0))
        contrast[[paste(tar, rr)]] <- data.frame(target = tar, r_um = rr, n_specimens = nrow(pair), delta_log2_mean_g = if (all(av > 0))
            log2(av[2]/av[1])
        else NA_real_)
    }
    list(curves = curves, mean = do.call(rbind, radial_mean), contrast = do.call(rbind, contrast), eligibility = do.call(rbind, eligibility),
        warnings = if (length(messages)) do.call(rbind, messages) else data.frame(specimen = character(), source = character(), target = character(),
            message = character()))
}

plot_outputs <- function (dist, radial, window_mode, config, run_dir)
{
    base_theme <- theme_classic(base_size = 11) + theme(axis.text = element_text(color = "black"), plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9), legend.position = "top", legend.title = element_blank())
    save_plot <- function(p, name, width, height) {
        ggsave(file.path(run_dir, paste0(name, ".pdf")), p, width = width, height = height, units = "in", device = grDevices::pdf, useDingbats = FALSE)
        if (requireNamespace("ragg", quietly = TRUE))
            ggsave(file.path(run_dir, paste0(name, ".png")), p, width = width, height = height, units = "in", dpi = 300, device = ragg::agg_png)
        else ggsave(file.path(run_dir, paste0(name, ".png")), p, width = width, height = height, units = "in", dpi = 300)
    }
    title_for <- function(t) if (t == TARGET_LEVELS[1])
        "STAT3+ THBS1+ macrophage"
    else "Myofibroblast"
    distance_panels <- radial_panels <- list()
    for (tar in TARGET_LEVELS) {
        s <- dist$summaries[dist$summaries$target == tar & dist$summaries$eligible_pair, , drop = FALSE]
        s$source <- factor(s$source, levels = SOURCE_LEVELS)
        tt <- dist$tests[dist$tests$target == tar, , drop = FALSE]
        ptext <- if (is.finite(tt$p_Holm))
            paste0("paired specimen medians; n = ", tt$n_specimens, "\nWilcoxon, Holm P = ", format.pval(tt$p_Holm, digits = 3, eps = 0.001))
        else paste0("paired specimen medians; n = ", tt$n_specimens, "; test unavailable")
        distance_panels[[tar]] <- ggplot(s, aes(source, median_distance_um, group = specimen)) + geom_line(color = "grey65", linewidth = 0.55) +
            geom_point(aes(fill = source), shape = 21, color = "grey20", size = 2.7) + scale_fill_manual(values = SOURCE_COLORS, drop = FALSE) +
            guides(fill = "none") + scale_x_discrete(labels = c("PT", "SPP1+\ninjured PT"), drop = FALSE) + scale_y_continuous(limits = c(0,
            NA), expand = expansion(mult = c(0.02, 0.12))) + labs(title = title_for(tar), subtitle = ptext, x = NULL, y = "Median nearest-observed distance (um)") +
            base_theme + theme(legend.position = "none")
        if (!nrow(s))
            distance_panels[[tar]] <- ggplot() + annotate("text", x = 0, y = 0, label = "No eligible specimens; test unavailable") + labs(title = title_for(tar)) +
                theme_void()
        q <- radial$mean[radial$mean$target == tar, , drop = FALSE]
        q$source <- factor(q$source, levels = SOURCE_LEVELS)
        ns <- sort(unique(q$n_specimens))
        ns_text <- if (!length(ns))
            "0; unavailable"
        else if (length(ns) == 1L)
            as.character(ns)
        else paste(range(ns), collapse = "--")
        radial_panels[[tar]] <- ggplot(q, aes(r_um, log2_mean_g, color = source)) + geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
            geom_line(linewidth = 0.85, na.rm = TRUE) + scale_color_manual(values = SOURCE_COLORS, drop = FALSE) + scale_x_continuous(limits = range(config$radius_um),
            breaks = pretty(range(config$radius_um), n = 5)) + labs(title = title_for(tar), subtitle = paste0("equal-specimen mean g(r); n = ",
            ns_text), x = "Distance from source centroid (um)", y = "log2 mean g(r)") + base_theme
    }
    for (tar in TARGET_LEVELS) if (!any(radial$mean$target == tar))
        radial_panels[[tar]] <- ggplot() + annotate("text", x = 0, y = 0, label = "No eligible radial estimates") + labs(title = title_for(tar)) +
            theme_void()
    main <- patchwork::wrap_plots(distance_panels[[1]], radial_panels[[1]], distance_panels[[2]], radial_panels[[2]], ncol = 2, guides = "collect") +
        patchwork::plot_annotation(title = paste0("Figure 6D | Spatial proximity across ", length(config$specimen_ids), " biopsy specimens"),
            caption = if (window_mode == "user_supplied_window")
                "Distance: specimen-level paired test. Radial curves: descriptive; supplied observation windows."
            else "Distance: specimen-level paired test. Radial curves: descriptive; FOV windows estimated from retained cell coordinates.")
    main <- main & theme(legend.position = "top")
    save_plot(main, "Figure6D_specimenwise", 9.5, 8.2)
    if (config$save_cell_violin && any(dist$cells$eligible_pair & is.finite(dist$cells$distance_um))) {
        v <- dist$cells[dist$cells$eligible_pair & is.finite(dist$cells$distance_um), , drop = FALSE]
        v$source <- factor(v$source, levels = SOURCE_LEVELS)
        p <- ggplot(v, aes(source, distance_um, fill = source)) + geom_violin(scale = "width", trim = TRUE, linewidth = 0.35) + geom_boxplot(width = 0.13,
            outlier.shape = NA, fill = "white", linewidth = 0.35) + facet_wrap(~target, scales = "free_y") + scale_fill_manual(values = SOURCE_COLORS) +
            scale_x_discrete(labels = c("PT", "SPP1+\ninjured PT")) + labs(title = "Exploratory cell-level distance distributions", subtitle = "Distances computed within specimens; no pooled-cell significance test",
            x = NULL, y = "Nearest observed target distance (um)") + base_theme + theme(legend.position = "none")
        save_plot(p, "Figure6D_cell_violin_exploratory", 8, 4.6)
        p2 <- p + facet_wrap(~specimen + target, scales = "free_y", ncol = 3) + labs(title = "Per-specimen source-cell distance distributions",
            subtitle = "No observed target = no distance distribution")
        save_plot(p2, "distance_violin_by_specimen", 11, 15)
    }
    if (!nrow(radial$curves)) {
        writeLines("Radial estimates unavailable: inspect radial_eligibility.csv. Missing estimates are not zero enrichment.", file.path(run_dir,
            "radial_unavailable.txt"))
        return(invisible(main))
    }
    c <- radial$curves
    c$source <- factor(c$source, levels = SOURCE_LEVELS)
    p <- ggplot(c[c$r_um > 0, ], aes(r_um, log2_g, color = source)) + geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") + geom_line(linewidth = 0.6,
        na.rm = TRUE) + facet_wrap(~specimen + target, scales = "free_y", ncol = 3) + scale_color_manual(values = SOURCE_COLORS) + labs(title = "Radial enrichment by specimen",
        subtitle = "Zero/nonpositive estimates are not assigned an arbitrary log pseudocount", x = "Radius (um)", y = "log2 g(r)") + base_theme
    save_plot(p, "radial_curves_by_specimen", 11, 15)
    invisible(main)
}


plot_manuscript_layout<-function(dist,radial,cfg,out) {
  panels<-list();labels<-list()
  for(tar in TARGET_LEVELS) {
    v<-dist$cells[dist$cells$target==tar&dist$cells$eligible_pair&is.finite(dist$cells$distance_um),,drop=FALSE]
    tt<-dist$tests[dist$tests$target==tar,,drop=FALSE]
    txt<-if(is.finite(tt$p_Holm))paste0("Specimen-level paired Wilcoxon; n = ",tt$n_specimens,
      "; Holm P = ",format.pval(tt$p_Holm,digits=3,eps=.0001))else "Specimen-level test unavailable"
    v$source<-factor(v$source,levels=SOURCE_LEVELS)
    p<-if(nrow(v))ggplot(v,aes(source,distance_um,fill=source))+
      geom_violin(scale="width",trim=TRUE,linewidth=.35)+geom_boxplot(width=.13,outlier.shape=NA,fill="white",linewidth=.35)+
      scale_fill_manual(values=SOURCE_COLORS)+scale_x_discrete(drop=FALSE,labels=c("PT","SPP1+ injured PT"))+
      labs(title=tar,subtitle=txt,x=NULL,y="Nearest observed target distance (um)")+theme_classic()+theme(legend.position="none") else
      ggplot()+annotate("text",x=0,y=0,label="No eligible distance distributions")+theme_void()
    q<-radial$mean[radial$mean$target==tar,,drop=FALSE]
    ns<-sort(unique(q$n_specimens));nt<-if(!length(ns))"0" else if(length(ns)==1L)as.character(ns) else paste(range(ns),collapse="-")
    rp<-if(nrow(q))ggplot(q,aes(r_um,log2_mean_g,color=source))+geom_hline(yintercept=0,linetype="dotted")+
      geom_line(linewidth=.8,na.rm=TRUE)+scale_color_manual(values=SOURCE_COLORS)+
      scale_x_continuous(limits=range(cfg$radius_um))+
      labs(title=tar,subtitle=paste0("Equal-specimen mean g(r); n = ",nt),x="Radius (um)",y="log2 mean g(r)",color=NULL)+theme_classic() else
      ggplot()+annotate("text",x=0,y=0,label="No eligible radial estimates")+theme_void()
    panels<-c(panels,list(p,rp))
    labels[[tar]]<-data.frame(target=tar,n_specimens_distance=tt$n_specimens,p_Holm=tt$p_Holm,
      n_source_cells_displayed=nrow(v),radial_n_min=if(length(ns))min(ns)else 0,
      radial_n_max=if(length(ns))max(ns)else 0)
  }
  p<-patchwork::wrap_plots(panels,ncol=2,guides="collect")+
    patchwork::plot_annotation(caption="Violins display cell distributions; P values compare specimen medians. Radial curves are descriptive.")
  save_panel(p,out,"Figure6D_cell_violin_specimen_statistics",11,8)
  write.csv(do.call(rbind,labels),file.path(out,"Figure6D_panel_statistics.csv"),row.names=FALSE)
}

####################  4. Analysis  #############################################

####################  Figure 6D: within-specimen nearest-target distances  ######
validate_cosmx_config(cfg)
checkpoint_file<-input_checkpoint(config,"01","CosMx_9biopsies_cells.rds")
checkpoint<-readRDS(checkpoint_file)
if(!identical(checkpoint$schema_version,1L)||!identical(checkpoint$settings,extraction_settings(cfg)))
  stop("Checkpoint settings differ from this script; do not mix expression layers/specimens/units.")
d<-checkpoint$cells;validate_cells(d,cfg)
window_file<-path_at(project_dir,cfg$observation_windows_rds)
cfg$observation_windows_rds<-window_file;config$cosmx<-cfg
out<-init_stage(config,"03",c(checkpoint_file,window_file))
write.csv(cell_manifest(d,cfg),file.path(out,"specimen_cell_counts.csv"),row.names=FALSE)
win<-make_windows(d,cfg,out)
dist<-distance_analysis(d,cfg)
for(nm in c("summaries","eligibility","paired","tests"))write.csv(dist[[nm]],
  file.path(out,paste0("distance_",nm,".csv")),row.names=FALSE)
con<-gzfile(file.path(out,"distance_per_source_cell.csv.gz"),"wt")
tryCatch(write.csv(dist$cells,con,row.names=FALSE),finally=close(con))
####################  Figure 6D: specimen-wise radial enrichment  ##############
radial<-radial_analysis(d,win$windows,cfg)
for(nm in c("curves","mean","contrast","eligibility","warnings"))write.csv(radial[[nm]],
  file.path(out,paste0("radial_",nm,".csv")),row.names=FALSE)
plot_outputs(dist,radial,win$mode,cfg,out)
plot_manuscript_layout(dist,radial,cfg,out)
checks<-data.frame(check=c("all_valid_distances_within_specimen","selected_specimens"),
  value=c(all(dist$cells$specimen[is.finite(dist$cells$distance_um)]==dist$cells$target_specimen[is.finite(dist$cells$distance_um)]),
  setequal(unique(d$Specimen_ID),cfg$specimen_ids)))
write.csv(checks,file.path(out,"completion_checks.csv"),row.names=FALSE)
if(!all(checks$value))stop("Spatial consistency check failed.")
writeLines(c("SCT/data >0; supplied author annotations retained.",
  "Distance queries confined to Specimen_ID; never join cells from different specimens.",
  "Two-sided paired Wilcoxon on specimen medians; exact=FALSE; correct=TRUE; Holm across two targets.",
  "Cell violins are descriptive; no pooled-cell P value is used.",
  "Radial curves: per-specimen isotropic pcfcross, default smoothing; equal-specimen mean g(r) then log2.",
  paste("Observation windows:",win$mode),
  "Estimated all-cell FOV rectangles are not validated tissue masks; retain this Methods description.",
  "No radial hypothesis test; unavailable values are not zero enrichment."),file.path(out,"analysis_notes.txt"))
finish_stage(out)
