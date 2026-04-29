script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_idx <- grep(file_arg, args, fixed = TRUE)

  if (length(file_idx) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", args[file_idx[1]], fixed = TRUE))))
  }

  normalizePath(getwd())
}

root <- script_dir()
source(file.path(root, "R", "load.R"))
source_project(root)

missing <- check_packages(c("ShapleyOutlier"))
if (length(missing) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing, collapse = ", "),
    "\nRun: Rscript scripts/install_packages.R",
    call. = FALSE
  )
}

set.seed(2026)

sim <- simulate_shapley_experiment(
  n = 150,
  p = 6,
  rho = 0.55,
  row_fraction = 0.06,
  cell_fraction = 0.03,
  magnitude = 6,
  seed = 2026
)

results <- list(
  original_classical = run_shapley_original(
    sim$x,
    estimator = "classical",
    q = 0.99,
    method_label = "Original: classical covariance"
  )
)

if (requireNamespace("robustbase", quietly = TRUE)) {
  results$extension_mcd <- run_robust_shapley(
    sim$x,
    estimator = "mcd",
    q = 0.99,
    method_label = "Extension: robust MCD covariance"
  )
}

metrics <- compare_methods(results, sim$truth)

dir.create(file.path(root, "outputs", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "outputs", "figures"), recursive = TRUE, showWarnings = FALSE)

write.csv(
  metrics,
  file.path(root, "outputs", "tables", "simulation_metrics.csv"),
  row.names = FALSE
)

saveRDS(
  list(simulation = sim, results = results, metrics = metrics),
  file.path(root, "outputs", "simulation_run.rds")
)

print(metrics)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::ggsave(
    filename = file.path(root, "outputs", "figures", "md_scores_original.png"),
    plot = plot_md_scores(results$original_classical),
    width = 7,
    height = 4,
    dpi = 150
  )

  ggplot2::ggsave(
    filename = file.path(root, "outputs", "figures", "shapley_heatmap_original.png"),
    plot = plot_shapley_heatmap(results$original_classical),
    width = 7,
    height = 5,
    dpi = 150
  )
}
