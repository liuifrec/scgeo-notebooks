#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1, 2)) stop("usage: diagnose_gse249479_official_augur.R OUTPUT_DIR [--versions-only]")
versions_only <- length(args) == 2 && identical(args[[2]], "--versions-only")
if (length(args) == 2 && !versions_only) stop("unknown option: ", args[[2]])
output_dir <- normalizePath(args[[1]], mustWork = TRUE)
feature_path <- file.path(output_dir, "official_augur_variance_selected_features.csv")
version_csv <- file.path(output_dir, "official_augur_dependency_versions.csv")
version_json <- file.path(output_dir, "official_augur_dependency_versions.json")
if (!versions_only && any(file.exists(c(feature_path, version_csv, version_json)))) stop("diagnostic output exists; refusing to overwrite")

required <- c("Augur", "data.table", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("missing packages: ", paste(missing, collapse = ", "))

feature_rows <- list()
if (!versions_only) {
  for (contrast in c("PBS_vs_TNF", "PBS_vs_LPS")) {
    result <- readRDS(file.path(output_dir, paste0("official_augur_", contrast, ".rds")))
    for (state in unique(as.character(result$cell_types))) {
      state_matrix <- result$X[, result$cell_types == state, drop = FALSE]
      selected <- Augur:::select_variance(state_matrix, var_quantile = 0.5, filter_negative_residuals = FALSE)
      feature_rows[[paste(contrast, state)]] <- data.frame(
        contrast = contrast,
        state = state,
        gene_id = rownames(selected),
        official_selection_method = "LOESS_residual_of_mean_over_sd_after_1pct_tail_trim",
        inference_status = "descriptive_only"
      )
      rm(state_matrix, selected); gc()
    }
    rm(result); gc()
  }
  utils::write.csv(do.call(rbind, feature_rows), feature_path, row.names = FALSE)
}

installed <- installed.packages()
dependencies <- tools::package_dependencies("Augur", db = installed, recursive = TRUE)[["Augur"]]
packages <- unique(c("Augur", dependencies, "jsonlite", "renv"))
packages <- packages[packages %in% rownames(installed)]
version_frame <- data.frame(
  package = packages,
  version = installed[packages, "Version"],
  library = installed[packages, "LibPath"],
  stringsAsFactors = FALSE
)
version_frame$source <- vapply(packages, function(package) {
  desc <- utils::packageDescription(package)
  if (identical(desc$RemoteType, "github") && !is.null(desc$RemoteUsername) && !is.null(desc$RemoteRepo)) {
    paste0(desc$RemoteUsername, "/", desc$RemoteRepo, "@", desc$RemoteSha)
  } else {
    "installed_CRAN_Bioconductor_or_system_package"
  }
}, character(1))
version_frame$inference_status <- "descriptive_only"
utils::write.csv(version_frame, version_csv, row.names = FALSE)
jsonlite::write_json(
  list(
    inference_status = "descriptive_only",
    R = R.version.string,
    platform = R.version$platform,
    packages = split(version_frame, seq_len(nrow(version_frame)))
  ),
  version_json, auto_unbox = TRUE, pretty = TRUE
)
cat("Recorded", nrow(version_frame), "R package versions and", length(feature_rows), "state/contrast feature sets.\n")
