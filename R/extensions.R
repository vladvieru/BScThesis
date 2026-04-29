run_robust_shapley <- function(x,
                               estimator = c("mcd", "ogk"),
                               q = 0.99,
                               cells = NULL,
                               method_label = NULL) {
  estimator <- match.arg(estimator)

  run_shapley_original(
    x = x,
    estimator = estimator,
    q = q,
    cells = cells,
    method_label = if (is.null(method_label)) {
      paste("Extension: robust", toupper(estimator), "covariance")
    } else {
      method_label
    }
  )
}

run_interaction_extension <- function(x,
                                      observation = NULL,
                                      mu = NULL,
                                      Sigma = NULL,
                                      estimator = c("classical", "mcd", "ogk")) {
  require_package("ShapleyOutlier")

  estimator <- match.arg(estimator)
  x <- as_numeric_matrix(x)

  if (is.null(mu) || is.null(Sigma)) {
    pars <- estimate_location_scatter(x, estimator = estimator)
  } else {
    Sigma <- ensure_positive_definite(Sigma)
    pars <- list(mu = as.numeric(mu), Sigma = Sigma, Sigma_inv = solve(Sigma), estimator = "provided")
  }

  if (is.null(observation)) {
    md2 <- mahalanobis_squared(x, pars$mu, Sigma_inv = pars$Sigma_inv)
    observation <- which.max(md2)
  }

  interaction <- ShapleyOutlier::shapley_interaction(
    x = as.numeric(x[observation, ]),
    mu = pars$mu,
    Sigma = pars$Sigma
  )

  dimnames(interaction) <- list(colnames(x), colnames(x))

  list(
    method = "Extension: pairwise Shapley interactions",
    observation = observation,
    estimator = pars$estimator,
    interaction = interaction
  )
}

summarise_group_contributions <- function(phi, groups) {
  phi <- as.matrix(phi)

  if (is.list(groups)) {
    out <- vapply(
      groups,
      function(vars) rowSums(phi[, vars, drop = FALSE]),
      FUN.VALUE = numeric(nrow(phi))
    )
  } else {
    if (is.null(names(groups))) {
      if (length(groups) != ncol(phi)) {
        stop("Unnamed groups must have one entry per variable.", call. = FALSE)
      }
      names(groups) <- colnames(phi)
    }

    groups <- groups[colnames(phi)]
    out <- vapply(
      unique(groups),
      function(group) rowSums(phi[, groups == group, drop = FALSE]),
      FUN.VALUE = numeric(nrow(phi))
    )
  }

  rownames(out) <- rownames(phi)
  out
}

run_bootstrap_stability <- function(x,
                                    B = 50,
                                    sample_fraction = 0.8,
                                    estimator = c("classical", "mcd", "ogk"),
                                    q = 0.99,
                                    seed = NULL) {
  require_package("ShapleyOutlier")

  estimator <- match.arg(estimator)
  x <- as_numeric_matrix(x)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- nrow(x)
  p <- ncol(x)
  sample_size <- max(p + 2, floor(n * sample_fraction))
  sample_size <- min(sample_size, n)

  phi_array <- array(NA_real_, dim = c(n, p, B), dimnames = list(rownames(x), colnames(x), NULL))
  row_flag_array <- matrix(FALSE, nrow = n, ncol = B)

  for (b in seq_len(B)) {
    idx <- sample.int(n, sample_size, replace = TRUE)
    pars <- estimate_location_scatter(x[idx, , drop = FALSE], estimator = estimator)

    shapley_object <- ShapleyOutlier::shapley(
      x = x,
      mu = pars$mu,
      Sigma = pars$Sigma,
      inverted = FALSE,
      check = TRUE
    )

    phi_array[, , b] <- extract_phi(shapley_object, template_x = x)

    md2 <- mahalanobis_squared(x, pars$mu, Sigma_inv = pars$Sigma_inv)
    row_flag_array[, b] <- flag_row_outliers(md2, p = p, q = q)
  }

  list(
    method = "Extension: bootstrap contribution stability",
    estimator = estimator,
    B = B,
    phi_mean = apply(phi_array, c(1, 2), mean),
    phi_sd = apply(phi_array, c(1, 2), stats::sd),
    row_outlier_rate = rowMeans(row_flag_array),
    phi_draws = phi_array
  )
}

default_method_specs <- function(include_robust = TRUE) {
  specs <- list(
    list(name = "original_classical", type = "original", estimator = "classical")
  )

  if (isTRUE(include_robust)) {
    specs <- c(
      specs,
      list(
        list(name = "extension_mcd", type = "robust", estimator = "mcd"),
        list(name = "extension_ogk", type = "robust", estimator = "ogk")
      )
    )
  }

  specs
}

run_method_grid <- function(x, specs = default_method_specs(), q = 0.99) {
  results <- list()

  for (spec in specs) {
    if (identical(spec$type, "original")) {
      results[[spec$name]] <- run_shapley_original(
        x,
        estimator = spec$estimator,
        q = q,
        method_label = spec$name
      )
    } else if (identical(spec$type, "robust")) {
      results[[spec$name]] <- run_robust_shapley(
        x,
        estimator = spec$estimator,
        q = q,
        method_label = spec$name
      )
    } else {
      stop("Unknown method type: ", spec$type, call. = FALSE)
    }
  }

  results
}
