as_numeric_matrix <- function(x, name = "x") {
  if (is.data.frame(x)) {
    numeric_cols <- vapply(x, is.numeric, FUN.VALUE = logical(1))
    if (!all(numeric_cols)) {
      stop(name, " must contain only numeric columns.", call. = FALSE)
    }
    x <- as.matrix(x)
  }

  if (is.atomic(x) && is.null(dim(x))) {
    x <- matrix(as.numeric(x), nrow = 1)
  }

  if (!is.matrix(x)) {
    stop(name, " must be a numeric vector, matrix, or data frame.", call. = FALSE)
  }

  if (!is.numeric(x)) {
    stop(name, " must be numeric.", call. = FALSE)
  }

  storage.mode(x) <- "double"

  if (any(!is.finite(x))) {
    stop(name, " must not contain NA, NaN, or infinite values.", call. = FALSE)
  }

  if (is.null(colnames(x))) {
    colnames(x) <- paste0("V", seq_len(ncol(x)))
  }

  if (is.null(rownames(x))) {
    rownames(x) <- as.character(seq_len(nrow(x)))
  }

  x
}

ensure_positive_definite <- function(Sigma, eps = 1e-8, max_iter = 8) {
  Sigma <- as.matrix(Sigma)
  storage.mode(Sigma) <- "double"

  if (nrow(Sigma) != ncol(Sigma)) {
    stop("Sigma must be square.", call. = FALSE)
  }

  Sigma <- (Sigma + t(Sigma)) / 2

  for (i in seq_len(max_iter + 1)) {
    attempt <- try(chol(Sigma), silent = TRUE)
    if (!inherits(attempt, "try-error")) {
      return(Sigma)
    }

    diag(Sigma) <- diag(Sigma) + eps * 10^(i - 1)
  }

  stop("Sigma is not positive definite, even after diagonal regularisation.", call. = FALSE)
}

invert_spd <- function(Sigma) {
  solve(ensure_positive_definite(Sigma))
}

estimate_location_scatter <- function(x, estimator = c("classical", "mcd", "ogk")) {
  estimator <- match.arg(estimator)
  x <- as_numeric_matrix(x)

  if (nrow(x) < 2) {
    stop("At least two observations are required to estimate covariance.", call. = FALSE)
  }

  if (estimator == "classical") {
    mu <- colMeans(x)
    Sigma <- stats::cov(x)
  } else if (estimator == "mcd") {
    require_package("robustbase")
    fit <- robustbase::covMcd(x)
    mu <- fit$center
    Sigma <- fit$cov
  } else if (estimator == "ogk") {
    require_package("robustbase")
    fit <- robustbase::covOGK(x, sigmamu = robustbase::s_Qn)
    mu <- fit$center
    Sigma <- fit$cov
  }

  Sigma <- ensure_positive_definite(Sigma)

  list(
    mu = as.numeric(mu),
    Sigma = Sigma,
    Sigma_inv = solve(Sigma),
    estimator = estimator
  )
}

mahalanobis_squared <- function(x, mu, Sigma = NULL, Sigma_inv = NULL) {
  x <- as_numeric_matrix(x)

  if (is.null(Sigma_inv)) {
    if (is.null(Sigma)) {
      stop("Provide Sigma or Sigma_inv.", call. = FALSE)
    }
    Sigma_inv <- solve(ensure_positive_definite(Sigma))
  }

  centered <- sweep(x, 2, as.numeric(mu), FUN = "-")
  as.numeric(rowSums((centered %*% Sigma_inv) * centered))
}

row_outlier_cutoff <- function(p, q = 0.99) {
  stats::qchisq(q, df = p)
}

flag_row_outliers <- function(md2, p, q = 0.99) {
  md2 > row_outlier_cutoff(p, q)
}

as_long_matrix <- function(mat, value_name = "value") {
  mat <- as.matrix(mat)

  out <- data.frame(
    observation = rep(seq_len(nrow(mat)), times = ncol(mat)),
    variable = rep(colnames(mat), each = nrow(mat)),
    value = as.vector(mat),
    row.names = NULL
  )

  names(out)[names(out) == "value"] <- value_name
  out
}

top_abs_mask <- function(values, top_k) {
  values <- as.matrix(values)
  mask <- matrix(FALSE, nrow(values), ncol(values), dimnames = dimnames(values))

  if (top_k <= 0) {
    return(mask)
  }

  top_k <- min(top_k, length(values))
  keep <- order(abs(values), decreasing = TRUE)[seq_len(top_k)]
  mask[keep] <- TRUE
  mask
}
