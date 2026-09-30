configure_project_library <- function(root = getwd()) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  activate <- file.path(root, "renv", "activate.R")
  if (file.exists(activate)) {
    source(activate)
  }
  invisible(.libPaths())
}

