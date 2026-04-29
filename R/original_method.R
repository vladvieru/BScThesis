extract_phi <- function(shapley_object, template_x = NULL) {
  phi <- if (is.list(shapley_object) && !is.null(shapley_object$phi)) {
    shapley_object$phi
  } else {
    shapley_object
  }

  phi <- as.matrix(phi)

  if (!is.null(template_x)) {
    template_x <- as_numeric_matrix(template_x)
    if (nrow(phi) == ncol(template_x) && ncol(phi) == 1 && nrow(template_x) == 1) {
      phi <- matrix(as.numeric(phi), nrow = 1)
    }

    dimnames(phi) <- dimnames(template_x)
  }

  phi
}

run_shapley_original <- function(x,
                                 mu = NULL,
                                 Sigma = NULL,
                                 estimator = c("classical", "mcd", "ogk"),
                                 q = 0.99,
                                 cells = NULL,
                                 shapley_method = "cellMCD",
                                 method_label = NULL) {
  require_package("ShapleyOutlier")

  estimator <- match.arg(estimator)
  x <- as_numeric_matrix(x)

  if (is.null(mu) || is.null(Sigma)) {
    pars <- estimate_location_scatter(x, estimator = estimator)
  } else {
    Sigma <- ensure_positive_definite(Sigma)
    pars <- list(
      mu = as.numeric(mu),
      Sigma = Sigma,
      Sigma_inv = solve(Sigma),
      estimator = "provided"
    )
  }

  shapley_object <- ShapleyOutlier::shapley(
    x = x,
    mu = pars$mu,
    Sigma = pars$Sigma,
    inverted = FALSE,
    method = shapley_method,
    check = TRUE,
    cells = cells
  )

  phi <- extract_phi(shapley_object, template_x = x)
  md2 <- mahalanobis_squared(x, pars$mu, Sigma_inv = pars$Sigma_inv)
  cutoff <- row_outlier_cutoff(ncol(x), q = q)

  out <- list(
    method = if (is.null(method_label)) "Original Shapley-Mahalanobis" else method_label,
    estimator = pars$estimator,
    x = x,
    mu = pars$mu,
    Sigma = pars$Sigma,
    Sigma_inv = pars$Sigma_inv,
    phi = phi,
    md2 = md2,
    q = q,
    cutoff = cutoff,
    row_outlier = md2 > cutoff,
    cells = cells,
    shapley_object = shapley_object
  )

  class(out) <- c("shapley_outlier_result", "list")
  out
}

print.shapley_outlier_result <- function(x, ...) {
  cat(x$method, "\n", sep = "")
  cat("Estimator: ", x$estimator, "\n", sep = "")
  cat("Observations: ", nrow(x$x), " | Variables: ", ncol(x$x), "\n", sep = "")
  cat("Row outliers at q = ", x$q, ": ", sum(x$row_outlier), "\n", sep = "")
  invisible(x)
}

run_scd_original <- function(x,
                             mu = NULL,
                             Sigma = NULL,
                             estimator = c("classical", "mcd", "ogk"),
                             q = 0.99,
                             step_size = 0.1,
                             min_deviation = 0,
                             max_iter = 1000,
                             cells = NULL) {
  require_package("ShapleyOutlier")

  estimator <- match.arg(estimator)
  x <- as_numeric_matrix(x)

  if (is.null(mu) || is.null(Sigma)) {
    pars <- estimate_location_scatter(x, estimator = estimator)
  } else {
    Sigma <- ensure_positive_definite(Sigma)
    pars <- list(mu = as.numeric(mu), Sigma = Sigma, Sigma_inv = solve(Sigma), estimator = "provided")
  }

  fit <- ShapleyOutlier::SCD(
    x = x,
    mu = pars$mu,
    Sigma = pars$Sigma,
    Sigma_inv = pars$Sigma_inv,
    step_size = step_size,
    min_deviation = min_deviation,
    max_iter = max_iter,
    q = q,
    cells = cells
  )

  list(
    method = "SCD cellwise detection",
    estimator = pars$estimator,
    fit = fit,
    phi = extract_phi(fit, template_x = x),
    x_original = x,
    x_imputed = fit$x
  )
}

run_moe_original <- function(x,
                             mu = NULL,
                             Sigma = NULL,
                             estimator = c("classical", "mcd", "ogk"),
                             q = 0.99,
                             step_size = 0.1,
                             min_deviation = 0,
                             local = TRUE,
                             max_iter = 1000,
                             cells = NULL) {
  require_package("ShapleyOutlier")

  estimator <- match.arg(estimator)
  x <- as_numeric_matrix(x)

  if (is.null(mu) || is.null(Sigma)) {
    pars <- estimate_location_scatter(x, estimator = estimator)
  } else {
    Sigma <- ensure_positive_definite(Sigma)
    pars <- list(mu = as.numeric(mu), Sigma = Sigma, Sigma_inv = solve(Sigma), estimator = "provided")
  }

  fit <- ShapleyOutlier::MOE(
    x = x,
    mu = pars$mu,
    Sigma = pars$Sigma,
    Sigma_inv = pars$Sigma_inv,
    step_size = step_size,
    min_deviation = min_deviation,
    local = local,
    max_iter = max_iter,
    q = q,
    cells = cells
  )

  list(
    method = "MOE cellwise detection",
    estimator = pars$estimator,
    fit = fit,
    phi = extract_phi(fit, template_x = x),
    x_original = x,
    x_imputed = fit$x
  )
}
