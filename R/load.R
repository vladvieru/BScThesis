find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)

  repeat {
    if (file.exists(file.path(current, "DESCRIPTION")) &&
        dir.exists(file.path(current, "R"))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find project root. Run from inside the BScThesis project.")
    }

    current <- parent
  }
}

source_project <- function(root = find_project_root()) {
  files <- c(
    "packages.R",
    "utils.R",
    "original_method.R",
    "extensions.R",
    "simulation.R",
    "evaluation.R",
    "plotting.R"
  )

  for (file in files) {
    source(file.path(root, "R", file), local = .GlobalEnv)
  }

  invisible(root)
}
