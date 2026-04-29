script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_idx <- grep(file_arg, args, fixed = TRUE)

  if (length(file_idx) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", args[file_idx[1]], fixed = TRUE))))
  }

  normalizePath(getwd())
}

root <- normalizePath(file.path(script_dir(), ".."), mustWork = FALSE)

if (!file.exists(file.path(root, "R", "load.R"))) {
  root <- normalizePath(getwd())
}

source(file.path(root, "R", "load.R"))
source_project(root)

sim <- simulate_shapley_experiment(n = 40, p = 4, seed = 10)

stopifnot(is.matrix(sim$x))
stopifnot(all(dim(sim$x) == dim(sim$truth$cell_outlier)))
stopifnot(length(sim$truth$row_outlier) == nrow(sim$x))

if (requireNamespace("ShapleyOutlier", quietly = TRUE)) {
  result <- run_shapley_original(sim$x, estimator = "classical")
  metrics <- evaluate_result(result, sim$truth)
  stopifnot(is.data.frame(metrics))
  stopifnot(nrow(metrics) > 0)
  stopifnot(max(abs(rowSums(result$phi) - result$md2)) < 1e-8)
} else {
  message("Skipping ShapleyOutlier smoke test because the package is not installed.")
}

message("Smoke test passed.")
