# ============================================================================
# 01_KPMP_Atlas_FigureS4
# KPMP atlas / Figure S4
# Public, sequential workflow — 2026-10-12
# This R script contains executable R code; use Rscript or source().
# Start a fresh R session. Edit Configuration; do not edit Helpers routinely.
# Input data and software dependencies must be obtained separately; see README.md.
# ============================================================================

# 1. Configuration -------------------------------------------------------------
config <- list(
    input_rds = Sys.getenv("KIDNEY_KPMP_RDS","data/HumanKidney.rds"),  # Original consortium-processed atlas
    output_dir = Sys.getenv("KIDNEY_RESULTS_DIR", "results"),
    donor_id_column = "donor_id",             # Consortium-supplied donor identifier; never infer donor IDs from cell barcodes
    lineage_map = utils::read.delim(text = "cell_type	lineage
podocyte	epithelial cell
epithelial cell of proximal tubule	epithelial cell
parietal epithelial cell	epithelial cell
papillary tips cell	epithelial cell
kidney connecting tubule epithelial cell	epithelial cell
kidney distal convoluted tubule epithelial cell	epithelial cell
kidney loop of Henle thick ascending limb epithelial cell	epithelial cell
kidney loop of Henle thin ascending limb epithelial cell	epithelial cell
kidney loop of Henle thin descending limb epithelial cell	epithelial cell
kidney collecting duct principal cell	epithelial cell
kidney collecting duct intercalated cell	epithelial cell
endothelial cell	stromal cell
kidney interstitial fibroblast	stromal cell
renal interstitial pericyte	stromal cell
mast cell	immune cell
T cell	immune cell
B cell	immune cell
plasma cell	immune cell
plasmacytoid dendritic cell, human	immune cell
conventional dendritic cell	immune cell
mononuclear phagocyte	immune cell
kidney interstitial alternatively activated macrophage	immune cell
non-classical monocyte	immune cell
lymphocyte	immune cell
neutrophil	immune cell
neural cell	neural cell
",
    check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = ""),
    source_S4B = utils::read.delim(text = "order	feature_id
1	ENSG00000231079
2	ENSG00000004799
3	ENSG00000096060
4	ENSG00000109906
5	ENSG00000116285
6	ENSG00000229807
7	ENSG00000119138
8	ENSG00000100234
9	ENSG00000138798
10	ENSG00000118985
11	ENSG00000101443
12	ENSG00000008517
13	ENSG00000142669
14	ENSG00000051523
15	ENSG00000196126
16	ENSG00000034510
17	ENSG00000105640
18	ENSG00000118785
19	ENSG00000117984
20	ENSG00000134419
21	ENSG00000074803
22	ENSG00000184831
23	ENSG00000179399
24	ENSG00000152953
25	ENSG00000036828
26	ENSG00000151704
27	ENSG00000182621
28	ENSG00000164946
29	ENSG00000214128
30	ENSG00000113889
",
    check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = ""),
    # Study-curated S4D list: shared/high-in-other-group and duplicate markers were excluded.
    # Preserve this recorded selection; no unreported numeric specificity cutoff is inferred.
    source_S4D = utils::read.delim(text = "order	feature_id
1	ENSG00000116285
2	ENSG00000197375
3	ENSG00000138798
4	ENSG00000204531
5	ENSG00000166016
6	ENSG00000022567
7	ENSG00000066813
8	ENSG00000171234
9	ENSG00000233154
10	ENSG00000155974
11	ENSG00000101443
12	ENSG00000118785
13	ENSG00000185499
14	ENSG00000125144
15	ENSG00000137673
16	ENSG00000162896
17	ENSG00000105518
18	ENSG00000213625
19	ENSG00000196754
20	ENSG00000113361
21	ENSG00000074803
22	ENSG00000187783
23	ENSG00000036828
24	ENSG00000113889
25	ENSG00000166912
26	ENSG00000006453
27	ENSG00000152527
28	ENSG00000091136
29	ENSG00000138036
30	ENSG00000231535
31	ENSG00000104368
32	ENSG00000188517
33	ENSG00000124440
34	ENSG00000131634
35	ENSG00000077943
36	ENSG00000067715
37	ENSG00000204301
38	ENSG00000185666
39	ENSG00000116962
40	ENSG00000130300
41	ENSG00000165949
42	ENSG00000169908
43	ENSG00000140937
44	ENSG00000164692
45	ENSG00000175745
46	ENSG00000177469
47	ENSG00000163430
48	ENSG00000277443
49	ENSG00000139174
50	ENSG00000183853
51	ENSG00000178031
52	ENSG00000249669
53	ENSG00000138190
54	ENSG00000101333
55	ENSG00000101384
56	ENSG00000011295
57	ENSG00000083535
58	ENSG00000028203
59	ENSG00000110066
60	ENSG00000138443
61	ENSG00000186684
62	ENSG00000114948
63	ENSG00000185634
64	ENSG00000277363
65	ENSG00000259916
66	ENSG00000164318
67	ENSG00000078018
68	ENSG00000225914
69	ENSG00000143341
70	ENSG00000250159
71	ENSG00000009709
72	ENSG00000145794
73	ENSG00000253976
74	ENSG00000158887
75	ENSG00000226476
76	ENSG00000258637
77	ENSG00000168530
78	ENSG00000054938
79	ENSG00000018625
80	ENSG00000164122
81	ENSG00000163873
82	ENSG00000125851
83	ENSG00000182050
84	ENSG00000145242
85	ENSG00000187260
86	ENSG00000167601
87	ENSG00000101096
88	ENSG00000124225
89	ENSG00000065361
90	ENSG00000076706
91	ENSG00000177575
92	ENSG00000079215
93	ENSG00000260314
94	ENSG00000141506
95	ENSG00000185811
96	ENSG00000121966
97	ENSG00000136250
98	ENSG00000155926
99	ENSG00000134516
100	ENSG00000084070
101	ENSG00000090382
102	ENSG00000167286
103	ENSG00000169442
104	ENSG00000211772
105	ENSG00000126264
106	ENSG00000102879
107	ENSG00000128340
108	ENSG00000204482
109	ENSG00000136490
110	ENSG00000158869
111	ENSG00000173744
112	ENSG00000082996
113	ENSG00000107581
114	ENSG00000133657
115	ENSG00000198265
116	ENSG00000278311
117	ENSG00000129515
118	ENSG00000133961
119	ENSG00000189079
120	ENSG00000198162
",
    check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = ""),
    enrichment_tsv = NULL,              # Optional original S4E export; format in README.md
    expected_seurat_version = "5.3.0",
    allow_version_mismatch = FALSE,
    seed = 42L
)
restart_failed_stage <- FALSE           # TRUE archives only an incomplete stage before retry

# 2. Libraries and helpers -----------------------------------------------------
# Packages are checked below; no automatic installation, download, or setwd().
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


cell_tests <- function (d, adjust = "bonferroni")
{
    need("rstatix")
    if (length(unique(d$value)) < 2L)
        return(list(omnibus = data.frame(p = NA_real_, status = "constant values; test undefined"), pairs = data.frame()))
    if (length(unique(d$condition)) < 2L)
        return(list(omnibus = data.frame(p = NA_real_, status = "fewer than two groups"), pairs = data.frame()))
    list(omnibus = as.data.frame(rstatix::kruskal_test(d, value ~ condition)), pairs = as.data.frame(rstatix::dunn_test(d, value ~ condition,
        p.adjust.method = adjust)))
}

# Explicit metadata and feature handling; no network annotation lookup.
disease_levels <- c("Healthy_Reference", "AKI", "CKD")
lineage_levels <- c("epithelial cell", "stromal cell", "neural cell", "immune cell")
disease_colors <- c(Healthy_Reference = "#C9C9C9", AKI = "#FBA425", CKD = "#E66342")
marker_ids <- c(CD163 = "ENSG00000177575", MRC1 = "ENSG00000260314",
                SPP1 = "ENSG00000118785", FN1 = "ENSG00000115414",
                STAT3 = "ENSG00000168610", ITGAX = "ENSG00000140678",
                CD44 = "ENSG00000026508", THBS1 = "ENSG00000137801")

map_condition <- function(x) {
    key <- c("Normal Reference" = "Healthy_Reference",
             "Acute Kidney Injury (AKI)" = "AKI",
             "Chronic Kidney Disease" = "CKD",
             "Hypertensive CKD (H-CKD)" = "CKD",
             "COVID-AKI" = NA_character_,
             "Diabetic Kidney Disease (DKD)" = NA_character_)
    if (anyNA(x) || any(!x %in% names(key)))
        stop("Unrecognized/missing condition.long values: ",
             paste(unique(x[is.na(x) | !x %in% names(key)]), collapse = ", "),
             ". Check the input release; do not silently drop new conditions.")
    factor(unname(key[as.character(x)]), levels = disease_levels)
}

map_lineage <- function(x, mapping) {
    check_columns(mapping, c("cell_type", "lineage"), "lineage_map.tsv")
    check_ids(mapping$cell_type, "Mapped cell types")
    if (anyNA(mapping$lineage) || any(!mapping$lineage %in% lineage_levels))
        stop("Invalid broad lineage in lineage_map.tsv.")
    idx <- match(x, mapping$cell_type)
    if (anyNA(idx)) stop("Unmapped cell_type labels: ", paste(unique(x[is.na(idx)]), collapse = ", "))
    factor(mapping$lineage[idx], levels = lineage_levels)
}

resolve_features <- function(features, requested, aliases = NULL) {
    check_ids(features, "RNA feature IDs")
    # Ensembl version suffixes may vary; symbols are accepted only when unambiguous.
    clean <- sub("\\.[0-9]+$", "", features)
    vapply(seq_along(requested), function(i) {
        candidates <- requested[i]
        if (!is.null(aliases)) candidates <- c(candidates, aliases[i])
        hit <- which(features %in% candidates | clean %in% candidates)
        if (length(hit) != 1L)
            stop("Feature ", requested[i], " matched ", length(hit), " RNA rows. Check annotation; no duplicate collapse was performed.")
        features[hit]
    }, character(1), USE.NAMES = FALSE)
}

check_seurat_version <- function(config) {
    actual <- as.character(utils::packageVersion("Seurat"))
    if (actual != config$expected_seurat_version && !isTRUE(config$allow_version_mismatch))
        stop("Seurat ", config$expected_seurat_version, " is the manuscript version; installed: ", actual,
             ". Use the matching environment or explicitly set allow_version_mismatch = TRUE and review differences.")
    if (getRversion() != "4.4.2") warning("Manuscript R version: 4.4.2; running ", getRversion(), ". Runtime versions will be recorded.")
    future::plan(future::sequential)
    options(mc.cores = 1L)
}

prepare_rna <- function(object) {
    if (!inherits(object, "Seurat") || !"RNA" %in% names(object@assays))
        stop("Expected a processed Seurat object with an RNA assay.")
    SeuratObject::DefaultAssay(object) <- "RNA"
    layers <- SeuratObject::Layers(object[["RNA"]], search = NA)
    dat <- grep("^data($|\\.)", layers, value = TRUE)
    if (!length(dat)) stop("No normalized RNA data layer. Supply the processed atlas; this workflow does not infer a normalization recipe.")
    if (length(dat) > 1L || dat != "data") {
        if ("data" %in% dat && length(dat) > 1L)
            stop("Both joined and split RNA data layers are present. Resolve ambiguous input layers first.")
        object <- SeuratObject::JoinLayers(object, assay = "RNA", layers = "data", new = "data")
    }
    mat <- SeuratObject::LayerData(object, assay = "RNA", layer = "data")
    if (!identical(colnames(mat), colnames(object)))
        stop("RNA data columns do not match the full object cell order.")
    check_ids(colnames(object), "Cell IDs")
    check_ids(rownames(mat), "RNA feature IDs")
    rr <- range(mat)
    if (any(!is.finite(rr)) || rr[1] < 0)
        stop("RNA data contain negative/nonfinite values; expected nonnegative log-normalized expression.")
    object
}

validate_umap <- function(object, reduction = "umap") {
    if (!reduction %in% names(object@reductions)) stop("Missing reduction: ", reduction)
    emb <- SeuratObject::Embeddings(object, reduction = reduction)
    if (ncol(emb) < 2 || !identical(rownames(emb), colnames(object)) || any(!is.finite(emb[, 1:2])))
        stop("Invalid UMAP coordinates or cell order.")
    invisible(TRUE)
}

count_cells <- function(metadata, fields) {
    check_columns(metadata, fields, "metadata")
    z <- as.data.frame(table(metadata[, fields, drop = FALSE]), stringsAsFactors = FALSE)
    names(z)[ncol(z)] <- "n_cells"
    z
}

donor_manifest <- function(metadata, donor_column = NULL) {
    if (is.null(donor_column)) return(NULL)
    check_columns(metadata, c(donor_column, "disease_group"), "metadata")
    ids <- as.character(metadata[[donor_column]])
    if (anyNA(ids) || any(!nzchar(trimws(ids)))) stop("Missing donor IDs; do not replace them with cell IDs.")
    disease <- as.character(metadata$disease_group)
    if (any(vapply(split(disease, ids), function(z) length(unique(z)), integer(1)) != 1L))
        stop("Donor IDs span disease groups. Verify IDs and cohort metadata before donor counting.")
    ans <- aggregate(rep(1L, length(ids)), list(donor_id = ids, disease_group = disease), sum)
    names(ans)[3] <- "n_cells"
    ans
}

export_metadata <- function(object, out, donor_column = NULL) {
    md <- object[[]]
    write_tsv(data.frame(column = names(md), n_unique = vapply(md, function(z) length(unique(z)), integer(1))),
              file.path(out, "metadata_columns.tsv"))
    write_tsv(count_cells(md, "disease_group"), file.path(out, "cells_by_disease.tsv"))
    write_tsv(count_cells(md, c("subclass.full", "disease_group")), file.path(out, "cells_by_subclass_and_disease.tsv"))
    dm <- donor_manifest(md, donor_column)
    if (is.null(dm)) {
        writeLines(c("Donor ID column was not specified. No donor counts or donor-level inference were generated.",
                     "Review metadata_columns.tsv and set donor_id_column to a verified participant ID column."),
                   file.path(out, "DONOR_ID_NOT_CONFIGURED.txt"))
    } else {
        write_tsv(dm, file.path(out, "donor_cell_counts.tsv"))
        write_tsv(data.frame(disease_group = disease_levels,
                            n_donors = vapply(disease_levels, function(g) sum(dm$disease_group == g), integer(1))),
                  file.path(out, "donors_by_disease.tsv"))
    }
}

input_checkpoint <- function(config) {
    p <- path_at(config$.root, config$atlas_rds)
    directory <- dirname(p)
    if (!file.exists(file.path(directory, "SUCCESS"))) stop("Run script 01 successfully before scripts 02-04.")
    m <- read_tsv(file.path(directory, "output_sha256.tsv"))
    i <- match(basename(p), m$file)
    if (is.na(i) || file_hash(p) != m$sha256[i]) stop("Atlas checkpoint checksum mismatch.")
    p
}

run_markers <- function(object, group_column, logfc, seed, out, stem) {
    ids <- droplevels(factor(object[[]][[group_column]]))
    if (anyNA(ids) || nlevels(ids) < 2L || any(table(ids) < 3L))
        stop("Marker groups must contain at least 3 cells and at least two nonmissing groups: ", group_column)
    SeuratObject::Idents(object) <- ids
    groups <- levels(ids)
    counts <- as.integer(table(ids))
    write_tsv(data.frame(group = groups, n_in_group = counts, n_in_rest = length(ids) - counts,
                        comparator = "All remaining cells in this analysis object"),
              file.path(out, paste0(stem, "_comparisons.tsv")))
    warnings <- character()
    markers <- withCallingHandlers(
        Seurat::FindAllMarkers(object, assay = "RNA", slot = "data", test.use = "wilcox",
                              min.pct = 0.25, logfc.threshold = logfc, only.pos = TRUE,
                              return.thresh = 0.01, max.cells.per.ident = Inf,
                              random.seed = seed, base = 2, densify = FALSE, verbose = TRUE),
        warning = function(w) { warnings <<- c(warnings, conditionMessage(w)) })
    writeLines(if (length(warnings)) unique(warnings) else "No warnings.", file.path(out, paste0(stem, "_warnings.txt")))
    # Seurat may warn rather than stop when a DE comparison fails. Do not mark such a stage successful.
    if (any(grepl("not performed|failed|not joined|error", warnings, ignore.case = TRUE)))
        stop("Marker analysis reported a failed/skipped test; see the warning log.")
    if (!nrow(markers)) stop("No markers were returned; review warnings/input instead of plotting an empty result.")
    check_columns(markers, c("gene", "cluster", "avg_log2FC", "p_val", "p_val_adj"), "FindAllMarkers output")
    markers$cluster <- as.character(markers$cluster)
    write_tsv(markers, file.path(out, paste0(stem, "_all_returned.tsv")))
    markers
}

top_markers <- function(markers, groups, n = 10L) {
    check_columns(markers, c("gene", "cluster", "avg_log2FC"), "markers")
    z <- markers[is.finite(markers$avg_log2FC) & markers$avg_log2FC > 0, , drop = FALSE]
    parts <- lapply(groups, function(g) {
        d <- z[as.character(z$cluster) == g, , drop = FALSE]
        d <- d[order(-d$avg_log2FC, d$gene, method = "radix"), , drop = FALSE]
        d <- d[!duplicated(d$gene), , drop = FALSE]
        if (nrow(d) < n) stop("Fewer than ", n, " positive markers for ", g, "; review the release/results.")
        d <- head(d, n); d$rank_within_group <- seq_len(n); d
    })
    do.call(rbind, parts)
}

curated_lineage_markers <- function(markers, reference, groups, n = 10L) {
    check_columns(reference, c("order", "feature_id"), "Curated S4D list")
    reference <- reference[order(reference$order), , drop = FALSE]
    if (nrow(reference) != length(groups) * n || anyDuplicated(reference$order) ||
        anyDuplicated(reference$feature_id)) stop("Invalid curated S4D membership/order.")
    cluster <- rep(groups, each = n)
    index <- match(paste(cluster, reference$feature_id), paste(markers$cluster, markers$gene))
    if (anyNA(index)) stop("Curated S4D genes are absent from their recomputed positive-marker groups.")
    selected <- markers[index, , drop = FALSE]
    selected$rank_by_log2FC <- vapply(seq_len(nrow(selected)), function(i) {
        d <- markers[markers$cluster == selected$cluster[i], , drop = FALSE]
        d <- d[order(-d$avg_log2FC, d$gene, method = "radix"), , drop = FALSE]
        match(selected$gene[i], d$gene)
    }, integer(1))
    selected$rank_within_group <- rep(seq_len(n), length(groups))
    selected$selection <- "Study-curated list after excluding shared/high-in-other-group and duplicate markers"
    selected
}

mean_linear_expression <- function(logdata, features, labels, groups) {
    if (ncol(logdata) != length(labels) || anyNA(labels)) stop("Misaligned expression/group labels.")
    z <- logdata[features, , drop = FALSE]
    result <- vapply(groups, function(g) {
        keep <- which(labels == g)
        if (!length(keep)) stop("No cells for group ", g)
        # Equivalent to AverageExpression(layer = 'data'): mean in non-log space.
        as.numeric(Matrix::rowMeans(expm1(z[, keep, drop = FALSE])))
    }, numeric(length(features)))
    dimnames(result) <- list(features, groups)
    result
}

save_plot <- function(plot, path, width, height) {
    ggplot2::ggsave(path, plot = plot, width = width, height = height, units = "in", limitsize = FALSE,
                   device = grDevices::pdf, useDingbats = FALSE)
}

marker_heatmap <- function(object, selected, group_column, groups, out, stem) {
    # Repeated genes selected by different groups occupy repeated rows, preserving selection provenance.
    features <- unique(selected$gene)
    mat <- SeuratObject::LayerData(object, assay = "RNA", layer = "data")
    avg <- mean_linear_expression(mat, features, as.character(object[[]][[group_column]]), groups)
    write_tsv(data.frame(feature_id = rownames(avg), avg, check.names = FALSE), file.path(out, paste0(stem, "_mean_linear_expression.tsv")))
    z <- t(scale(t(avg)))
    constant <- rownames(z)[apply(z, 1, function(r) any(!is.finite(r)))]
    write_tsv(data.frame(feature_id = constant), file.path(out, paste0(stem, "_constant_features.tsv")))
    # A constant row has no relative enrichment; display as zero, report it explicitly.
    if (length(constant)) z[constant, ] <- 0
    selected$row_id <- seq_len(nrow(selected))
    z <- z[selected$gene, , drop = FALSE]
    df <- data.frame(row_id = rep(selected$row_id, times = length(groups)),
                     group = factor(rep(groups, each = nrow(selected)), levels = groups),
                     row_z = as.vector(z))
    df$row_id <- factor(df$row_id, levels = rev(selected$row_id))
    row_labels <- setNames(selected$gene, selected$row_id)
    p <- ggplot2::ggplot(df, ggplot2::aes(group, row_id, fill = row_z)) + ggplot2::geom_tile() +
        ggplot2::scale_y_discrete(labels = row_labels) +
        ggplot2::scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0) +
        ggplot2::labs(x = NULL, y = NULL, fill = "Row z-score") + ggplot2::theme_minimal(base_size = 9) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank())
    save_plot(p, file.path(out, paste0(stem, ".pdf")), max(7, length(groups) * 0.7), max(6, nrow(selected) * 0.13))
    write_tsv(selected, file.path(out, paste0(stem, "_selected_markers.tsv")))
    write_tsv(df, file.path(out, paste0(stem, "_row_z_scores.tsv")))
}

compare_source_features <- function(selected, reference, out, stem, curated = FALSE) {
    ref <- reference
    check_columns(ref, c("order", "feature_id"), "Original fixed feature list")
    ref <- ref[order(ref$order), ]
    a <- sub("\\.[0-9]+$", "", selected$gene)
    b <- ref$feature_id
    n <- max(length(a), length(b))
    length(a) <- length(b) <- n
    d <- data.frame(position = seq_len(n), recomputed_feature = a, source_fixed_feature = b,
                    same_at_position = !is.na(a) & !is.na(b) & a == b)
    write_tsv(d, file.path(out, paste0(stem, "_source_list_comparison.tsv")))
    if (!all(d$same_at_position)) {
        if (curated) message(stem, ": unfiltered top-10 and study-curated selections differ; both are exported separately.")
        else warning(stem, ": recomputed top markers differ from the original fixed list. Differences are listed in the exported reference-comparison table.")
    }
}

subset_recluster <- function(atlas, config, out) {
    md <- atlas[[]]
    keep <- rownames(md)[as.character(md$subclass.full) %in% config$subclasses]
    if (length(keep) <= config$npcs + 1L) stop("Too few cells for the requested subset/PCA.")
    object <- subset(atlas, cells = keep)
    object <- prepare_rna(object)
    features <- SeuratObject::VariableFeatures(object[["RNA"]])
    features <- features[features %in% rownames(object[["RNA"]])]
    if (length(features) <= config$npcs) stop("Insufficient inherited RNA variable features. Supply the original processed atlas with its variable-feature annotation.")
    write_tsv(data.frame(feature_id = features), file.path(out, "inherited_RNA_variable_features.tsv"))
    object <- Seurat::ScaleData(object, assay = "RNA", features = features, verbose = FALSE)
    object <- Seurat::RunPCA(object, assay = "RNA", features = features, npcs = config$npcs,
                             seed.use = config$seed, verbose = FALSE)
    if (ncol(SeuratObject::Embeddings(object, "pca")) < config$npcs) stop("PCA returned fewer dimensions than requested.")
    object <- Seurat::RunUMAP(object, reduction = "pca", dims = seq_len(config$npcs),
                              umap.method = "uwot", metric = "cosine", n.neighbors = 30L,
                              n.components = 2L, min.dist = 0.3, spread = 1,
                              seed.use = config$seed, verbose = FALSE)
    object <- Seurat::FindNeighbors(object, reduction = "pca", dims = seq_len(config$npcs),
                                    k.param = 20L, nn.method = "annoy", n.trees = 50L,
                                    prune.SNN = 1/15, verbose = FALSE)
    object <- Seurat::FindClusters(object, graph.name = "RNA_snn", resolution = config$resolution,
                                   algorithm = 1L, random.seed = config$cluster_seed, verbose = FALSE)
    object$analysis_cluster <- as.character(SeuratObject::Idents(object))
    validate_umap(object)
    write_tsv(count_cells(object[[]], c("analysis_cluster", "disease_group")), file.path(out, "cluster_cell_counts.tsv"))
    object
}

spp1_cluster_state <- function(markers, clusters, spp1_feature, alpha = 0.05) {
    check_columns(markers, c("gene", "cluster", "p_val_adj", "avg_log2FC"), "subset markers")
    high <- unique(as.character(markers$cluster[markers$gene == spp1_feature &
                          !is.na(markers$p_val_adj) & markers$p_val_adj < alpha & markers$avg_log2FC > 0]))
    data.frame(cluster = clusters, SPP1_state = ifelse(clusters %in% high, "SPP1-high", "Other clusters"))
}

umap_plot <- function(object, column, split = NULL, title = NULL) {
    p <- Seurat::DimPlot(object, reduction = "umap", group.by = column, split.by = split,
                        pt.size = 0.1, label = FALSE, raster = FALSE)
    if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
    p
}

expression_panels <- function(object, out, stem, donor_column, dunn_adjust) {
    fid <- resolve_features(rownames(object[["RNA"]]), marker_ids["SPP1"], "SPP1")
    write_tsv(data.frame(symbol = "SPP1", requested_id = unname(marker_ids["SPP1"]), actual_id = fid),
              file.path(out, "SPP1_feature_mapping.tsv"))
    mat <- SeuratObject::LayerData(object, assay = "RNA", layer = "data")
    md <- object[[]]
    d <- data.frame(cell_id = colnames(object), condition = md$disease_group, value = as.numeric(mat[fid, ]))
    d$condition <- factor(d$condition, levels = disease_levels)
    tests <- cell_tests(d, adjust = dunn_adjust)
    write_tsv(tests$omnibus, file.path(out, paste0(stem, "_cell_level_Kruskal_Wallis.tsv")))
    write_tsv(tests$pairs, file.path(out, paste0(stem, "_cell_level_Dunn.tsv")))
    if (!is.null(donor_column)) {
        donor_manifest(md, donor_column) # validate mapping before creating summaries
        d$donor_id <- as.character(md[[donor_column]])
        med <- aggregate(d$value, list(donor_id = d$donor_id, condition = d$condition), median)
        names(med)[3] <- "median_log_expression"
        write_tsv(med, file.path(out, paste0(stem, "_donor_medians_descriptive.tsv")))
    }
    write_tsv(d, file.path(out, paste0(stem, "_cell_expression.tsv")))
    writeLines(c("Kruskal-Wallis and Dunn tests use cells, not independent donors.",
                 paste0("Dunn adjustment within the three disease comparisons: ", dunn_adjust, "."),
                 "No donor-level hypothesis test is performed by this script."), file.path(out, "STATISTICAL_UNIT.txt"))
    p <- ggplot2::ggplot(d, ggplot2::aes(condition, value, fill = condition)) +
        ggplot2::geom_violin(scale = "width", trim = TRUE) +
        ggplot2::scale_fill_manual(values = disease_colors, drop = FALSE) +
        ggplot2::labs(x = NULL, y = "SPP1 log-normalized expression") + ggplot2::theme_classic() +
        ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
    save_plot(p, file.path(out, paste0(stem, "_violin.pdf")), 5, 4)
    # Global limits shared across disease panels, not three separately rescaled gradients.
    upper <- max(d$value)
    if (upper <= 0) upper <- 1
    p <- Seurat::FeaturePlot(object, features = fid, reduction = "umap", split.by = "disease_group",
                             order = TRUE, keep.scale = "feature", min.cutoff = 0, max.cutoff = upper,
                             cols = c("grey90", "#542788"), pt.size = 0.15, raster = FALSE)
    save_plot(p, file.path(out, paste0(stem, "_SPP1_by_disease.pdf")), 12, 4)
    invisible(fid)
}

plot_supplied_enrichment <- function(file, groups, out) {
    if (is.null(file)) {
        writeLines(c("Figure S4E was not generated because enrichment_tsv is NULL.",
                     "Enrichment is performed externally; this script plots a supplied enrichment-results table.",
                     "Marker tables are exported separately; a supplied enrichment table determines the displayed terms and adjusted P values.",
                     "To plot verified results, set enrichment_tsv to a TSV with term, group, p_adjusted columns; see README.md."),
                   file.path(out, "FigureS4E_NOT_GENERATED.txt"))
        return(invisible(NULL))
    }
    d <- read_tsv(file)
    check_columns(d, c("term", "group", "p_adjusted"), "Enrichment export")
    if (!nrow(d) || anyNA(d$term) || any(!nzchar(d$term)) || any(!d$group %in% groups) || anyNA(d$group))
        stop("Enrichment term/group labels are missing or unrecognized.")
    if (!is.numeric(d$p_adjusted) || any(!is.finite(d$p_adjusted)) || any(d$p_adjusted <= 0 | d$p_adjusted > 1))
        stop("Use original nonzero adjusted P values in (0,1]; rounded zeros cannot be transformed without a documented display floor.")
    if (anyDuplicated(d[c("term", "group")])) stop("Duplicate term/group rows in enrichment export.")
    terms <- unique(d$term)
    grid <- expand.grid(term = terms, group = groups, stringsAsFactors = FALSE)
    grid <- merge(grid, d[c("term", "group", "p_adjusted")], by = c("term", "group"), all.x = TRUE, sort = FALSE)
    grid$minus_log10_adjusted_p <- -log10(grid$p_adjusted)
    grid$term <- factor(grid$term, levels = rev(terms))
    grid$group <- factor(grid$group, levels = groups)
    p <- ggplot2::ggplot(grid, ggplot2::aes(group, term, fill = minus_log10_adjusted_p)) +
        ggplot2::geom_tile() + ggplot2::scale_fill_gradient(low = "white", high = "#B2182B", na.value = "grey85") +
        ggplot2::labs(x = NULL, y = NULL, fill = "-log10(adj. P)") + ggplot2::theme_minimal(base_size = 9) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank())
    save_plot(p, file.path(out, "FigureS4E_supplied_enrichment.pdf"), 11, max(5, length(terms) * 0.25))
    write_tsv(grid, file.path(out, "FigureS4E_plotted_values.tsv"))
    writeLines("Grey cells are absent from the supplied export; they were not assigned P = 1 or interpreted as untested/non-significant without evidence.",
               file.path(out, "FigureS4E_missing_cells.txt"))
}

# 3. Analysis -----------------------------------------------------------------
config$.root <- locate_project()
need(c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork", "digest", "future"))
check_seurat_version(config)
input <- path_at(config$.root, config$input_rds)
lineage_map <- config$lineage_map
source_B <- config$source_S4B
source_D <- config$source_S4D
enrichment_file <- path_at(config$.root, config$enrichment_tsv)
out <- init_stage(config, "01_Atlas", c(input, enrichment_file))

atlas <- readRDS(input)
atlas <- prepare_rna(atlas)
check_columns(atlas[[]], c("condition.long", "cell_type", "subclass.full"), "KPMP metadata")
validate_umap(atlas)
mapped <- map_condition(as.character(atlas$condition.long))
write_tsv(data.frame(condition = as.character(atlas$condition.long), disease_group = as.character(mapped)) |>
              unique(), file.path(out, "condition_mapping.tsv"))
write_tsv(count_cells(atlas[[]], "condition.long"), file.path(out, "cells_before_condition_filter.tsv"))
keep <- colnames(atlas)[!is.na(mapped)]
atlas$disease_group <- mapped
atlas <- subset(atlas, cells = keep)
atlas <- prepare_rna(atlas)
atlas$lineage <- map_lineage(as.character(atlas$cell_type), lineage_map)
group_order <- as.vector(t(outer(lineage_levels, disease_levels, paste, sep = "__")))
atlas$lineage_disease <- factor(paste(atlas$lineage, atlas$disease_group, sep = "__"), levels = group_order)
if (any(table(atlas$lineage_disease) < 3L)) stop("One or more manuscript lineage-disease groups have fewer than 3 cells.")
export_metadata(atlas, out, config$donor_id_column)
write_tsv(count_cells(atlas[[]], c("lineage", "disease_group")), file.path(out, "cells_by_lineage_and_disease.tsv"))
expected <- c(Healthy_Reference = 107701L, AKI = 65171L, CKD = 47776L)
observed <- as.integer(table(factor(atlas$disease_group, levels = disease_levels)))
count_audit <- data.frame(disease_group = disease_levels, manuscript_n_cells = unname(expected),
                         observed_n_cells = observed, matches = observed == unname(expected))
write_tsv(count_audit, file.path(out, "manuscript_cell_count_comparison.tsv"))
if (!all(count_audit$matches)) warning("Cell counts differ from the study reference; the cohort audit lists observed and reference counts. Cell selection follows the configured inclusion criteria.")

# Preserve the supplied atlas UMAP. No new global normalization/integration/QC is performed.
save_plot(umap_plot(atlas, "disease_group", split = "disease_group"),
          file.path(out, "FigureS4A_atlas_by_disease.pdf"), 12, 4)
save_plot(umap_plot(atlas, "lineage", split = "disease_group"),
          file.path(out, "FigureS4C_lineages_by_disease.pdf"), 13, 4)

# Figure S4B: disease vs all other retained diseases.
mk_disease <- run_markers(atlas, "disease_group", logfc = 0.25, seed = config$seed, out, "S4B_disease")
top_disease <- top_markers(mk_disease, disease_levels, n = 10L)
marker_heatmap(atlas, top_disease, "disease_group", disease_levels, out, "FigureS4B")
compare_source_features(top_disease, source_B, out, "FigureS4B")

# Figure S4D: each COMBINED lineage-disease group vs ALL other combined groups.
# This is NOT a within-lineage pairwise disease comparison and NOT a projection of mouse bulk gene groups.
mk_lineage <- run_markers(atlas, "lineage_disease", logfc = 0.25, seed = config$seed, out, "S4D_lineage_disease")
top_lineage <- top_markers(mk_lineage, group_order, n = 10L)
# Preserve both the unfiltered ranking and the manuscript's recorded curation.
marker_heatmap(atlas, top_lineage, "lineage_disease", group_order, out, "FigureS4D_unfiltered_top10")
compare_source_features(top_lineage, source_D, out, "FigureS4D", curated = TRUE)
curated_lineage <- curated_lineage_markers(mk_lineage, source_D, group_order)
marker_heatmap(atlas, curated_lineage, "lineage_disease", group_order, out, "FigureS4D")
writeLines(c("FigureS4D uses the 120-gene study-curated list embedded in source_S4D.",
             "Positive markers were ranked by average log2FC; shared/high-in-other-group and duplicate genes were excluded during study curation.",
             "The recorded list preserves the final membership and order; no numeric specificity threshold is inferred.",
             "FigureS4D_unfiltered_top10 shows the automatic ranking before this curation.",
             "Selected-marker exports retain recomputed fold changes and pre-curation ranks."),
           file.path(out, "FigureS4D_selection_notes.txt"))

# Marker-list exports support external enrichment analysis; Figure S4E requires a supplied enrichment-results table.
write_tsv(mk_lineage[!is.na(mk_lineage$p_val_adj) & mk_lineage$p_val_adj < 0.05, , drop = FALSE],
          file.path(out, "lineage_disease_positive_markers_adjusted_p_lt_0.05.tsv"))
write_tsv(data.frame(feature_id = rownames(atlas[["RNA"]])), file.path(out, "RNA_feature_universe.tsv"))
plot_supplied_enrichment(enrichment_file, group_order, out)

saveRDS(atlas, file.path(out, "KPMP_filtered_atlas.rds"))
finish_stage(out)
message("Completed KPMP atlas analysis: ", out)
