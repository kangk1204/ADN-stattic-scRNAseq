# ============================================================================
# 01_Figure2E_GSVA — three gene sets / five conditions / 25 kidney samples
# Executable R code; run this R script with Rscript or source().
# Start a fresh R session. Original source: GSVAscore_2.R.
# ============================================================================

# 1. Configuration -------------------------------------------------------------
config <- list(
    tpm_csv = Sys.getenv("KIDNEY_TPM_CSV","data/kidney_GSVA.csv"),
    # NULL uses exact group headers (including repeated headers) or group + replicate suffix.
    # If headers are unrelated sample IDs, define a data.frame here with the columns
    # source_column, source_header, sample_id, group. No separate metadata file is read.
    sample_map = NULL,
    output_dir = Sys.getenv("KIDNEY_RESULTS_DIR", "results"),
    groups = c("Sham", "5D ADN", "5D ADN_Stattic", "14D ADN", "14D ADN_Stattic"),
    expected_n_per_group = 5L,
    set_files = c(GOBP_FERROPTOSIS = "gene_sets/GOBP_FERROPTOSIS.txt",
                  GOBP_LIPID_OXIDATION = "gene_sets/GOBP_LIPID_OXIDATION.txt",
                  WP_GLUTATHIONE_METABOLISM = "gene_sets/WP_GLUTATHIONE_METABOLISM.txt"),
    expected_gsva_version = "1.48.3",
    show_ns = FALSE,  # All pairwise results are exported regardless of display
    seed = 1234L
)
restart_failed_stage <- FALSE  # TRUE archives only an incomplete run before retry

# 2. Libraries and helpers -----------------------------------------------------
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


normalize_bulk_group <- function(x, groups) {
    # Exact aliases only; there is no column-position or fuzzy group inference.
    aliases <- c("Sham"="Sham", "sham"="Sham",
      "5D ADN"="5D ADN", "5D_ADN"="5D ADN", "5D.ADN"="5D ADN",
      "5D ADN_Stattic"="5D ADN_Stattic", "5D ADN + Stattic"="5D ADN_Stattic",
      "5D_ADN_Stattic"="5D ADN_Stattic", "5D.ADN_Stattic"="5D ADN_Stattic",
      "5D.ADN.Stattic"="5D ADN_Stattic",
      "14D ADN"="14D ADN", "14D_ADN"="14D ADN", "14D.ADN"="14D ADN",
      "14D ADN_Stattic"="14D ADN_Stattic", "14D ADN + Stattic"="14D ADN_Stattic",
      "14D_ADN_Stattic"="14D ADN_Stattic", "14D.ADN_Stattic"="14D ADN_Stattic",
      "14D.ADN.Stattic"="14D ADN_Stattic")
    y <- unname(aliases[x])
    # Supports explicit replicate suffixes such as Sham_1 or repaired duplicate names Sham.1.
    missing <- is.na(y)
    base <- sub("[._ ](?:rep)?[0-9]+$", "", x[missing], perl=TRUE)
    y[missing] <- unname(aliases[base])
    if (anyNA(y) || any(!y %in% groups))
        stop("TPM headers do not uniquely identify study groups: ",
             paste(unique(x[is.na(y)]), collapse=", "),
             ". Define config$sample_map inside this script using the GEO sample annotations; groups are never guessed from column positions.")
    y
}

resolve_bulk_samples <- function(headers, config) {
    m <- config$sample_map
    if (is.null(m)) {
        groups <- normalize_bulk_group(headers[-1], config$groups)
        # Unique source IDs are retained. Repeated group headers receive column-keyed
        # analysis IDs; these are not claimed to be original animal identifiers.
        ids <- headers[-1]
        if (anyDuplicated(ids)) ids <- paste0("TPM_column_",seq.int(2L,length(headers)))
        m <- data.frame(source_column=seq.int(2L,length(headers)), source_header=headers[-1],
                        sample_id=ids, group=groups, stringsAsFactors=FALSE)
    }
    check_columns(m, c("source_column", "source_header", "sample_id", "group"), "Embedded sample mapping")
    check_ids(as.character(m$sample_id), "Analysis sample IDs")
    if (!is.numeric(m$source_column) || anyNA(m$source_column) ||
        any(m$source_column != as.integer(m$source_column)) || anyDuplicated(m$source_column))
        stop("source_column must contain unique integer column positions; the gene column is column 1.")
    if (!setequal(m$source_column, seq.int(2L, length(headers))))
        stop("Embedded sample mapping must cover every TPM sample column exactly once.")
    actual_header <- headers[m$source_column]
    if (anyNA(m$source_header) || any(as.character(m$source_header) != actual_header))
        stop("TPM headers differ from config$sample_map. Check the source file and embedded mapping.")
    if (anyNA(m$group) || any(!m$group %in% config$groups)) stop("Unexpected group in sample mapping.")
    m$sample_id <- as.character(m$sample_id)
    m$group <- factor(m$group, levels=config$groups)
    counts <- table(m$group)
    if (any(counts != config$expected_n_per_group))
        stop("Expected ", config$expected_n_per_group, " samples in each study group; found: ", paste(counts,collapse=", "))
    m$source_header_actual <- actual_header
    m
}

read_bulk_csv <- function(config) {
    input <- path_at(config$.root, config$tpm_csv)
    if (!file.exists(input)) stop("Missing TPM input: ", input)
    x <- utils::read.csv(input, header=TRUE, row.names=NULL, check.names=FALSE,
                         colClasses="character", fileEncoding="UTF-8-BOM")
    if (!nrow(x) || ncol(x)<3L) stop("Expected mouse gene symbols in column 1 and TPM samples in the remaining columns.")
    m <- resolve_bulk_samples(names(x),config)
    genes <- x[[1]]
    check_ids(genes,"TPM gene symbols")
    if (any(genes!=trimws(genes))) stop("Whitespace in gene symbols; check the input.")
    vals <- lapply(m$source_column,function(j) suppressWarnings(as.numeric(x[[j]])))
    mat <- do.call(cbind,vals)
    if (anyNA(mat) || any(!is.finite(mat)) || any(mat<0))
        stop("TPM values must be numeric, finite, nonnegative and nonmissing. Invalid values are not replaced by zero.")
    dimnames(mat) <- list(genes,m$sample_id)
    list(tpm=mat,samples=m)
}

read_sets <- function(config) {
    paths <- vapply(config$set_files, function(f) path_at(config$.root, f), character(1))
    sets <- lapply(paths, function(file) {
        x <- read_tsv(file)
        if (!identical(names(x), "gene")) stop("Gene-set files must have exactly one column headed gene: ", basename(file))
        genes <- x$gene
        check_ids(genes, paste0(basename(file), " genes"))
        if (any(genes != trimws(genes))) stop("Whitespace in gene set: ", basename(file))
        genes
    })
    sets
}

match_sets <- function(sets, expr, out) {
    # Mirrors GSVA 1.48.3 .filterFeatures for the gsva method: SD < 1e-10 is treated as constant.
    sd_gene <- DelayedMatrixStats::rowSds(expr)
    retained <- is.finite(sd_gene) & sd_gene >= 1e-10
    if (sum(retained) < 2L) stop("Fewer than two variable genes in the input expression matrix.")
    write_tsv(data.frame(gene = rownames(expr), sd_log2_TPM_plus1 = sd_gene, retained_for_GSVA = retained),
              file.path(out, "expression_gene_filter.tsv"))
    detail <- do.call(rbind, lapply(names(sets), function(nm) {
        genes <- sets[[nm]]
        i <- match(genes, rownames(expr))
        status <- rep("not_in_TPM", length(genes))
        status[!is.na(i)] <- ifelse(retained[i[!is.na(i)]], "matched_variable", "matched_constant_removed")
        data.frame(gene_set = nm, gene = genes, status = status)
    }))
    write_tsv(detail, file.path(out, "gene_set_matching_details.tsv"))
    summary <- do.call(rbind, lapply(names(sets), function(nm) {
        d <- detail[detail$gene_set == nm, ]
        data.frame(gene_set = nm, n_input = nrow(d), n_matched_TPM = sum(d$status != "not_in_TPM"),
                   n_used_after_constant_filter = sum(d$status == "matched_variable"))
    }))
    write_tsv(summary, file.path(out, "gene_set_matching_summary.tsv"))
    if (any(summary$n_used_after_constant_filter < 1L))
        stop("A required gene set has no variable matched genes. Check gene IDs/species and gene_set_matching_summary.tsv.")
    filtered <- expr[retained, , drop = FALSE]
    # ALL variable background genes are kept, not just genes belonging to these three sets.
    list(expr = filtered, summary = summary)
}

calculate_scores <- function(expr, sets, config) {
    actual <- as.character(utils::packageVersion("GSVA"))
    if (actual != config$expected_gsva_version || actual != "1.48.3")
        stop("This reproduction targets GSVA 1.48.3; installed: ", actual,
             ". Do not silently switch algorithms/APIs. Run in the matching environment.")
    GSVA::gsva(expr, sets, method = "gsva", kcdf = "Gaussian", tau = 1,
               mx.diff = TRUE, abs.ranking = FALSE, min.sz = 1L, max.sz = Inf,
               parallel.sz = 1L, BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
}

significance_label <- function(p) {
    if (anyNA(p) || any(!is.finite(p)) || any(p < 0 | p > 1)) stop("Invalid P values for significance labels.")
    ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
}

score_statistics <- function(long, sets, groups) {
    a <- t <- list()
    for (nm in sets) {
        d <- long[long$category == nm, , drop = FALSE]
        d$group <- factor(d$group, levels = groups)
        model <- stats::aov(score ~ group, data = d)
        if (stats::df.residual(model) <= 0 || sum(stats::residuals(model)^2) <= 0)
            stop("ANOVA has no residual variance or degrees of freedom for ", nm, ".")
        # Retain the requested rstatix workflow. Also expose unrounded base-R values:
        # rstatix formats some P values; stats::TukeyHSD is the underlying test.
        aa <- as.data.frame(rstatix::anova_test(d, score ~ group))
        exact_a <- summary(model)[[1]]
        aa$p_unrounded <- exact_a["group", "Pr(>F)"]
        aa$F_unrounded <- exact_a["group", "F value"]
        a[[nm]] <- cbind(category = nm, aa)
        tt <- as.data.frame(rstatix::tukey_hsd(d, score ~ group))
        exact_t <- stats::TukeyHSD(model, which = "group", ordered = FALSE, conf.level = 0.95)$group
        idx <- match(paste(tt$group2, tt$group1, sep = "-"), rownames(exact_t))
        if (anyNA(idx) || nrow(tt) != choose(length(groups), 2)) stop("Unexpected Tukey contrast mapping.")
        tt$p.adj.rstatix <- tt$p.adj
        tt$p.adj <- exact_t[idx, "p adj"]
        tt$estimate <- exact_t[idx, "diff"]
        tt$conf.low <- exact_t[idx, "lwr"]
        tt$conf.high <- exact_t[idx, "upr"]
        tt$p.adj.signif <- significance_label(tt$p.adj)
        t[[nm]] <- cbind(category = nm, tt)
    }
    list(anova = do.call(rbind, a), tukey = do.call(rbind, t))
}

bracket_data <- function(long, tukey, show_ns = FALSE) {
    # Use unrounded adjusted P values, and keep all ten pairs per set in the output table.
    visible <- tukey[if (show_ns) rep(TRUE, nrow(tukey)) else tukey$p.adj < 0.05, , drop = FALSE]
    if (!nrow(visible)) return(visible)
    parts <- lapply(unique(as.character(visible$category)), function(nm) {
        tt <- visible[visible$category == nm, , drop = FALSE]
        x <- long$score[as.character(long$category) == nm]
        span <- diff(range(x))
        if (span <= 0) span <- max(1, abs(x[1]))
        tt$y.position <- max(x) + span * (0.15 + 0.13 * seq_len(nrow(tt)))
        tt
    })
    do.call(rbind, parts)
}

plot_scores <- function(long, tukey, config, out) {
    fill_cols <- setNames(c("#B2D1A8", "#FCEFCE", "#D2EEFD", "#F9D6AF", "#B3D4FC"), config$groups)
    dot_cols <- setNames(c("#3E8D27", "#F8D786", "#8DD4FB", "#F09837", "#4294F7"), config$groups)
    # shape 21 supports black borders independently of the dot fill.
    p <- ggplot2::ggplot(long, ggplot2::aes(group, score)) +
        ggplot2::geom_boxplot(ggplot2::aes(fill = group), outlier.shape = NA, width = 0.6,
                              color = "grey30", linewidth = 0.4) +
        ggplot2::scale_fill_manual(values = fill_cols) + ggnewscale::new_scale_fill() +
        ggplot2::geom_point(ggplot2::aes(fill = group), shape = 21, color = "black", stroke = 0.6,
                            size = 1.9, position = ggplot2::position_jitter(width = 0.10, seed = 1)) +
        ggplot2::scale_fill_manual(values = dot_cols) +
        ggplot2::facet_wrap(~category, scales = "free_y") +
        ggplot2::labs(x = NULL, y = "GSVA score") + ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
    brackets <- bracket_data(long, tukey, config$show_ns)
    write_tsv(brackets, file.path(out, "plotted_Tukey_brackets.tsv"))
    if (nrow(brackets)) {
        brackets$category <- factor(brackets$category, levels = levels(long$category))
        p <- p + ggpubr::stat_pvalue_manual(brackets, label = "p.adj.signif", xmin = "group1", xmax = "group2",
                                           y.position = "y.position", tip.length = 0.01,
                                           hide.ns = FALSE, inherit.aes = FALSE)
    }
    ggplot2::ggsave(file.path(out, "Figure2E_GSVA_boxplot.pdf"), p, width = 12, height = 5.5,
                    device = grDevices::pdf, useDingbats = FALSE)
    ggplot2::ggsave(file.path(out, "Figure2E_GSVA_boxplot.png"), p, width = 12, height = 5.5,
                    dpi = 300, bg = "white")
    invisible(p)
}

run_gsva <- function(config) {
    need(c("GSVA", "BiocParallel", "DelayedMatrixStats", "rstatix", "ggplot2", "ggpubr", "ggnewscale", "digest"))
    actual <- as.character(utils::packageVersion("GSVA"))
    if (actual != "1.48.3") stop("Expected GSVA 1.48.3; installed: ", actual, ". See README.md.")
    if (getRversion() != "4.4.2") warning("Reported manuscript R: 4.4.2; current: ", getRversion(), ". Actual versions are recorded.")
    dat <- read_bulk_csv(config)
    sets <- read_sets(config)
    inputs <- c(path_at(config$.root, config$tpm_csv),
                vapply(config$set_files, function(f) path_at(config$.root, f), character(1)))
    out <- init_stage(config, "01_Figure2E_GSVA", inputs)
    write_tsv(dat$samples, file.path(out, "sample_mapping_used.tsv"))
    write_tsv(data.frame(group = config$groups, n_samples = as.integer(table(dat$samples$group))),
              file.path(out, "sample_counts.tsv"))
    write_tsv(data.frame(sample_id = colnames(dat$tpm), total_TPM = colSums(dat$tpm)), file.path(out, "sample_TPM_totals.tsv"))
    expr <- log2(dat$tpm + 1)
    used <- match_sets(sets, expr, out)
    gsva_warnings <- character()
    scores <- withCallingHandlers(calculate_scores(used$expr, sets, config),
                warning = function(w) { gsva_warnings <<- c(gsva_warnings, conditionMessage(w)) })
    writeLines(if (length(gsva_warnings)) unique(gsva_warnings) else "No warnings during GSVA score calculation.",
               file.path(out, "GSVA_warnings.txt"))
    if (!is.matrix(scores) || !setequal(rownames(scores), names(sets)) ||
        !setequal(colnames(scores), dat$samples$sample_id) || any(!is.finite(scores)))
        stop("GSVA result has missing sets/samples or invalid values; inspect the input and warnings.")
    scores <- scores[names(sets), dat$samples$sample_id, drop = FALSE]
    long <- do.call(rbind, lapply(names(sets), function(nm)
        data.frame(sample_id = dat$samples$sample_id, group = dat$samples$group,
                   category = nm, score = as.numeric(scores[nm, ]))))
    long$group <- factor(long$group, levels = config$groups)
    long$category <- factor(long$category, levels = names(sets))
    stats <- score_statistics(long, names(sets), config$groups)
    write_tsv(long, file.path(out, "three_module_GSVA_scores_long.tsv"))
    wide <- data.frame(sample_id = dat$samples$sample_id, group = dat$samples$group, t(scores), check.names = FALSE)
    utils::write.csv(wide, file.path(out, "three_module_GSVA_scores.csv"), row.names = FALSE)
    utils::write.csv(stats$anova, file.path(out, "three_module_GSVA_ANOVA.csv"), row.names = FALSE)
    utils::write.csv(stats$tukey, file.path(out, "three_module_GSVA_TukeyHSD.csv"), row.names = FALSE)
    saveRDS(list(scores = scores, samples = dat$samples, gene_sets = sets,
                 matched_gene_counts = used$summary, anova = stats$anova, tukey = stats$tukey),
            file.path(out, "GSVA_analysis.rds"))
    plot_scores(long, stats$tukey, config, out)
    writeLines(c("Analysis unit: an input bulk RNA-seq sample; study design specifies one mouse per sample.",
                 "One-way ANOVA and all 10 Tukey comparisons were performed separately for each gene set.",
                 "Tukey correction is within each gene set, not a global adjustment across the three gene sets.",
                 "P.adj contains the unrounded stats::TukeyHSD value; p.adj.rstatix retains the rstatix-reported value.",
                 "A hidden nonsignificant bracket is not a missing test: all 30 comparisons are exported.",
                 "Scores are gene-set expression summaries, not direct measurements of ferroptotic cell death."),
               file.path(out, "ANALYSIS_NOTES.txt"))
    finish_stage(out)
    message("Completed GSVA analysis: ", out)
    invisible(out)
}

# 3. Analysis -----------------------------------------------------------------
config$.root <- locate_project()
run_gsva(config)
