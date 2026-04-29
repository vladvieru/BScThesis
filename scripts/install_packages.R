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

if (!file.exists(file.path(root, "R", "packages.R"))) {
  root <- normalizePath(getwd())
}

source(file.path(root, "R", "packages.R"))
install_project_packages(project_packages(include_optional = TRUE))
