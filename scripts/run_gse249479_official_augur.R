#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: run_gse249479_official_augur.R INPUT_DIR OUTPUT_DIR EXPECTED_COMMIT")
input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
expected_commit <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required <- c("Augur", "Matrix", "data.table", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("missing R packages: ", paste(missing, collapse = ", "))
if (as.character(packageVersion("Augur")) != "1.0.3") stop("unexpected Augur version")

targets <- file.path(output_dir, c(
  "official_augur_PBS_vs_TNF.rds", "official_augur_PBS_vs_LPS.rds",
  "official_augur_state_auc.csv", "official_augur_subsample_auc.csv",
  "official_augur_versions.json", "official_augur_session_info.txt",
  "official_augur_run_metadata.json"
))
if (any(file.exists(targets))) stop("official Augur result already exists; refusing to overwrite")

metadata <- data.table::fread(file.path(input_dir, "cells.csv"), data.table = FALSE)
genes <- data.table::fread(file.path(input_dir, "genes.csv"), data.table = FALSE)
connection <- gzfile(file.path(input_dir, "expression_genes_by_cells.mtx.gz"), "rb")
on.exit(close(connection), add = TRUE)
expr <- Matrix::readMM(connection)
close(connection)
on.exit(NULL, add = FALSE)
expr <- methods::as(expr, "dgCMatrix")
rownames(expr) <- genes$gene_id
colnames(expr) <- metadata$cell_id
if (!identical(colnames(expr), metadata$cell_id)) stop("cell order mismatch")
if (!identical(dim(expr), c(3000L, nrow(metadata)))) stop("unexpected expression dimensions")
if (!all(metadata$inference_status == "descriptive_only")) stop("inference status mismatch")

state_rows <- list()
subsample_rows <- list()
durations <- list()
for (treated in c("TNF", "LPS")) {
  keep <- metadata$condition %in% c("PBS", treated)
  contrast <- paste0("PBS_vs_", treated)
  meta_contrast <- droplevels(metadata[keep, c("cell_id", "condition", "state_detailed")])
  rownames(meta_contrast) <- meta_contrast$cell_id
  start <- proc.time()
  # No scientific parameter overrides: invoke the official documented defaults.
  result <- Augur::calculate_auc(
    expr[, keep, drop = FALSE], meta_contrast,
    label_col = "condition", cell_type_col = "state_detailed"
  )
  elapsed <- unname((proc.time() - start)[["elapsed"]])
  durations[[contrast]] <- elapsed
  auc <- as.data.frame(result$AUC)
  names(auc) <- c("state", "official_augur_auc")
  auc$contrast <- contrast
  auc$inference_status <- "descriptive_only"
  auc$uncertainty_status <- "computational_cell_resampling_stability_only"
  state_rows[[contrast]] <- auc
  subs <- as.data.frame(result$results)
  subs <- subs[subs$metric == "roc_auc", , drop = FALSE]
  names(subs)[names(subs) == "cell_type"] <- "state"
  subs$contrast <- contrast
  subs$inference_status <- "descriptive_only"
  subs$uncertainty_status <- "computational_cell_resampling_stability_only"
  subsample_rows[[contrast]] <- subs
  saveRDS(result, file.path(output_dir, paste0("official_augur_", contrast, ".rds")), compress = "xz")
  rm(result); gc()
}

state_frame <- do.call(rbind, state_rows)
subsample_frame <- do.call(rbind, subsample_rows)
utils::write.csv(state_frame, file.path(output_dir, "official_augur_state_auc.csv"), row.names = FALSE)
utils::write.csv(subsample_frame, file.path(output_dir, "official_augur_subsample_auc.csv"), row.names = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x
package_names <- unique(c("Augur", tools::package_dependencies("Augur", recursive = TRUE)[[1]]))
versions <- lapply(package_names, function(package) {
  if (!requireNamespace(package, quietly = TRUE)) return(NULL)
  description <- utils::packageDescription(package)
  list(package = package, version = as.character(packageVersion(package)), remote_sha = description$RemoteSha %||% NA_character_, remote_repo = description$RemoteRepo %||% NA_character_)
})
versions <- Filter(Negate(is.null), versions)
version_record <- list(
  inference_status = "descriptive_only",
  R = R.version.string,
  platform = R.version$platform,
  official_augur_source = "https://github.com/neurorestore/Augur",
  official_augur_expected_commit = expected_commit,
  packages = versions
)
jsonlite::write_json(version_record, file.path(output_dir, "official_augur_versions.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
capture.output(sessionInfo(), file = file.path(output_dir, "official_augur_session_info.txt"))
run_metadata <- list(
  inference_status = "descriptive_only",
  biological_replicates = FALSE,
  uncertainty_status = "computational_cell_resampling_stability_only",
  official_function = "Augur::calculate_auc",
  parameter_overrides = list(label_col = "condition", cell_type_col = "state_detailed"),
  scientific_defaults = list(
    n_subsamples = 50, subsample_size = 20, folds = 3, min_cells = 20,
    var_quantile = 0.5, feature_perc = 0.5, n_threads = 4,
    select_var = TRUE, augur_mode = "default", classifier = "rf",
    rf_params = list(trees = 100, mtry = 2, min_n = NULL, importance = "accuracy")
  ),
  dimensions_genes_by_cells = dim(expr),
  elapsed_seconds = durations
)
jsonlite::write_json(run_metadata, file.path(output_dir, "official_augur_run_metadata.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(jsonlite::toJSON(run_metadata, auto_unbox = TRUE, pretty = TRUE, null = "null"), "\n")
