safe_divide <- function(numerator, denominator) {
  ifelse(denominator == 0, NA_real_, numerator / denominator)
}

confusion_counts <- function(predicted, truth) {
  predicted <- as.logical(predicted)
  truth <- as.logical(truth)

  if (length(predicted) != length(truth)) {
    stop("predicted and truth must have the same length.", call. = FALSE)
  }

  c(
    tp = sum(predicted & truth),
    fp = sum(predicted & !truth),
    tn = sum(!predicted & !truth),
    fn = sum(!predicted & truth)
  )
}

classification_metrics <- function(predicted, truth) {
  counts <- confusion_counts(predicted, truth)

  precision <- safe_divide(counts[["tp"]], counts[["tp"]] + counts[["fp"]])
  recall <- safe_divide(counts[["tp"]], counts[["tp"]] + counts[["fn"]])
  specificity <- safe_divide(counts[["tn"]], counts[["tn"]] + counts[["fp"]])
  f1 <- safe_divide(2 * precision * recall, precision + recall)
  accuracy <- safe_divide(counts[["tp"]] + counts[["tn"]], sum(counts))

  c(
    counts,
    precision = precision,
    recall = recall,
    specificity = specificity,
    f1 = f1,
    accuracy = accuracy
  )
}

evaluate_row_detection <- function(result, truth_rows) {
  classification_metrics(result$row_outlier, truth_rows)
}

evaluate_cell_detection <- function(phi,
                                    truth_cells,
                                    predicted_cells = NULL,
                                    top_k = NULL,
                                    threshold = NULL) {
  phi <- as.matrix(phi)
  truth_cells <- as.matrix(truth_cells)

  if (!all(dim(phi) == dim(truth_cells))) {
    stop("phi and truth_cells must have the same dimensions.", call. = FALSE)
  }

  if (is.null(predicted_cells)) {
    if (!is.null(threshold)) {
      predicted_cells <- abs(phi) > threshold
    } else {
      if (is.null(top_k)) {
        top_k <- sum(truth_cells)
      }
      predicted_cells <- top_abs_mask(phi, top_k = top_k)
    }
  }

  metrics <- classification_metrics(as.vector(predicted_cells), as.vector(truth_cells))
  ranks <- rank(-abs(as.vector(phi)), ties.method = "average")
  true_ranks <- ranks[as.vector(truth_cells)]

  c(
    metrics,
    mean_true_rank = if (length(true_ranks) == 0) NA_real_ else mean(true_ranks),
    median_true_rank = if (length(true_ranks) == 0) NA_real_ else stats::median(true_ranks)
  )
}

evaluate_result <- function(result, truth) {
  row_metrics <- evaluate_row_detection(result, truth$row_outlier)
  cell_metrics <- evaluate_cell_detection(result$phi, truth$cell_outlier)

  data.frame(
    method = result$method,
    estimator = result$estimator,
    metric_scope = c(rep("row", length(row_metrics)), rep("cell", length(cell_metrics))),
    metric = c(names(row_metrics), names(cell_metrics)),
    value = as.numeric(c(row_metrics, cell_metrics)),
    row.names = NULL
  )
}

compare_methods <- function(results, truth) {
  out <- do.call(
    rbind,
    lapply(results, evaluate_result, truth = truth)
  )

  rownames(out) <- NULL
  out
}
