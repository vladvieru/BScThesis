make_correlation_matrix <- function(p,
                                    rho = 0.5,
                                    structure = c("exchangeable", "ar1", "block")) {
  structure <- match.arg(structure)

  if (structure == "exchangeable") {
    Sigma <- matrix(rho, p, p)
    diag(Sigma) <- 1
  } else if (structure == "ar1") {
    idx <- seq_len(p)
    Sigma <- outer(idx, idx, function(i, j) rho^abs(i - j))
  } else if (structure == "block") {
    Sigma <- matrix(rho / 3, p, p)
    block <- rep(seq_len(ceiling(p / 3)), each = 3, length.out = p)
    Sigma[outer(block, block, "==")] <- rho
    diag(Sigma) <- 1
  }

  ensure_positive_definite(Sigma)
}

simulate_clean_data <- function(n,
                                p,
                                mu = rep(0, p),
                                Sigma = make_correlation_matrix(p),
                                seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  Sigma <- ensure_positive_definite(Sigma)
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  x <- sweep(z %*% chol(Sigma), 2, mu, FUN = "+")

  colnames(x) <- paste0("V", seq_len(p))
  rownames(x) <- as.character(seq_len(n))
  x
}

inject_cellwise_outliers <- function(x,
                                     fraction = 0.03,
                                     magnitude = 6,
                                     direction = c("random", "positive", "negative"),
                                     seed = NULL) {
  direction <- match.arg(direction)
  x <- as_numeric_matrix(x)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- nrow(x)
  p <- ncol(x)
  k <- if (fraction <= 0) 0 else max(1, round(fraction * n * p))
  truth <- matrix(FALSE, n, p, dimnames = dimnames(x))

  if (k == 0) {
    return(list(x = x, truth = truth, index = integer()))
  }

  index <- sample.int(n * p, k, replace = FALSE)
  rows <- ((index - 1) %% n) + 1
  cols <- ((index - 1) %/% n) + 1

  scales <- apply(x, 2, stats::sd)
  scales[!is.finite(scales) | scales == 0] <- 1

  signs <- switch(
    direction,
    random = sample(c(-1, 1), k, replace = TRUE),
    positive = rep(1, k),
    negative = rep(-1, k)
  )

  for (i in seq_len(k)) {
    x[rows[i], cols[i]] <- x[rows[i], cols[i]] + signs[i] * magnitude * scales[cols[i]]
    truth[rows[i], cols[i]] <- TRUE
  }

  list(x = x, truth = truth, index = index)
}

inject_rowwise_outliers <- function(x,
                                    fraction = 0.05,
                                    variables = NULL,
                                    magnitude = 6,
                                    direction = c("positive", "negative", "random"),
                                    seed = NULL) {
  direction <- match.arg(direction)
  x <- as_numeric_matrix(x)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- nrow(x)
  p <- ncol(x)
  k <- if (fraction <= 0) 0 else max(1, round(fraction * n))

  row_truth <- rep(FALSE, n)
  cell_truth <- matrix(FALSE, n, p, dimnames = dimnames(x))

  if (k == 0) {
    return(list(x = x, row_truth = row_truth, cell_truth = cell_truth, rows = integer()))
  }

  rows <- sample.int(n, k, replace = FALSE)

  if (is.null(variables)) {
    variables <- seq_len(max(1, ceiling(p / 3)))
  }

  scales <- apply(x, 2, stats::sd)
  scales[!is.finite(scales) | scales == 0] <- 1

  signs <- switch(
    direction,
    random = sample(c(-1, 1), k, replace = TRUE),
    positive = rep(1, k),
    negative = rep(-1, k)
  )

  for (i in seq_along(rows)) {
    x[rows[i], variables] <- x[rows[i], variables] + signs[i] * magnitude * scales[variables]
    row_truth[rows[i]] <- TRUE
    cell_truth[rows[i], variables] <- TRUE
  }

  list(x = x, row_truth = row_truth, cell_truth = cell_truth, rows = rows)
}

simulate_shapley_experiment <- function(n = 200,
                                        p = 8,
                                        rho = 0.5,
                                        correlation = c("exchangeable", "ar1", "block"),
                                        row_fraction = 0.05,
                                        cell_fraction = 0.03,
                                        magnitude = 6,
                                        seed = 1) {
  correlation <- match.arg(correlation)

  Sigma <- make_correlation_matrix(p, rho = rho, structure = correlation)
  x_clean <- simulate_clean_data(n, p, Sigma = Sigma, seed = seed)

  row_injected <- inject_rowwise_outliers(
    x_clean,
    fraction = row_fraction,
    magnitude = magnitude,
    seed = seed + 1
  )

  cell_injected <- inject_cellwise_outliers(
    row_injected$x,
    fraction = cell_fraction,
    magnitude = magnitude,
    seed = seed + 2
  )

  cell_truth <- row_injected$cell_truth | cell_injected$truth
  row_truth <- row_injected$row_truth | rowSums(cell_truth) > 0

  list(
    x_clean = x_clean,
    x = cell_injected$x,
    Sigma = Sigma,
    truth = list(
      row_outlier = row_truth,
      cell_outlier = cell_truth,
      rowwise_rows = row_injected$rows,
      cellwise_index = cell_injected$index
    )
  )
}
