plot_shapley_heatmap <- function(result_or_phi, title = NULL) {
  require_package("ggplot2")

  phi <- if (is.list(result_or_phi) && !is.null(result_or_phi$phi)) {
    result_or_phi$phi
  } else {
    result_or_phi
  }

  phi <- as.matrix(phi)
  long <- as_long_matrix(phi, value_name = "phi")

  ggplot2::ggplot(long, ggplot2::aes(x = variable, y = factor(observation), fill = phi)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
    ggplot2::labs(
      title = title %||% "Shapley-Mahalanobis contributions",
      x = "Variable",
      y = "Observation",
      fill = "Contribution"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

plot_observation_contributions <- function(result, observation = NULL, title = NULL) {
  require_package("ggplot2")

  if (is.null(observation)) {
    observation <- which.max(result$md2)
  }

  dat <- data.frame(
    variable = colnames(result$phi),
    phi = as.numeric(result$phi[observation, ]),
    row.names = NULL
  )
  dat <- dat[order(abs(dat$phi), decreasing = TRUE), ]
  dat$variable <- factor(dat$variable, levels = rev(dat$variable))

  ggplot2::ggplot(dat, ggplot2::aes(x = variable, y = phi, fill = phi > 0)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac")) +
    ggplot2::labs(
      title = title %||% paste("Observation", observation, "contributions"),
      x = "Variable",
      y = "Shapley contribution"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

plot_md_scores <- function(result, title = NULL) {
  require_package("ggplot2")

  dat <- data.frame(
    observation = seq_along(result$md2),
    md2 = result$md2,
    row_outlier = result$row_outlier
  )

  ggplot2::ggplot(dat, ggplot2::aes(x = observation, y = md2, color = row_outlier)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = result$cutoff, linetype = "dashed", color = "#444444") +
    ggplot2::scale_color_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#444444")) +
    ggplot2::labs(
      title = title %||% "Squared Mahalanobis distances",
      x = "Observation",
      y = "Squared Mahalanobis distance",
      color = "Flagged"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

plot_method_comparison <- function(metrics,
                                   metric = "f1",
                                   scope = "cell",
                                   title = NULL) {
  require_package("ggplot2")

  dat <- metrics[metrics$metric == metric & metrics$metric_scope == scope, , drop = FALSE]

  ggplot2::ggplot(dat, ggplot2::aes(x = method, y = value, fill = estimator)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title %||% paste(scope, metric, "by method"),
      x = "Method",
      y = metric,
      fill = "Estimator"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
