project_packages <- function(include_optional = TRUE) {
  core <- c("ShapleyOutlier")
  optional <- c("ggplot2", "robustbase", "rmarkdown", "knitr")

  if (isTRUE(include_optional)) {
    unique(c(core, optional))
  } else {
    core
  }
}

check_packages <- function(packages = project_packages()) {
  packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
}

require_package <- function(package) {
  missing <- check_packages(package)

  if (length(missing) > 0) {
    stop(
      "Package '", package, "' is required. ",
      "Run: Rscript scripts/install_packages.R",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

install_project_packages <- function(packages = project_packages(), repos = "https://cloud.r-project.org") {
  missing <- check_packages(packages)

  if (length(missing) == 0) {
    message("All requested packages are already installed.")
    return(invisible(character()))
  }

  install.packages(missing, repos = repos)
  invisible(missing)
}
