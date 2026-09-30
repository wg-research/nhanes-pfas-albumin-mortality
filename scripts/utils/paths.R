find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, "Database"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the project root from: ", start)
    }
    current <- parent
  }
}

project_paths <- function(root = find_project_root()) {
  list(
    root = root,
    raw = file.path(root, "Database", "nhanes"),
    mortality = file.path(root, "Database", "nhanes", "linked_mortality"),
    processed = file.path(root, "data_rebuild"),
    results = file.path(root, "results", "latest_release"),
    audit = file.path(root, "results", "latest_release", "data_audit"),
    logs = file.path(root, "results", "latest_release", "logs")
  )
}

ensure_output_directories <- function(paths) {
  output_dirs <- unname(unlist(paths[c("processed", "results", "audit", "logs")]))
  invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

