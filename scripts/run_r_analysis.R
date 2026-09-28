#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
project_root <- dirname(dirname(script_path))

data_dir <- if (length(args)) args[[1]] else file.path(project_root, "data", "local")
if (!dir.exists(data_dir)) {
  stop("Data directory does not exist: ", data_dir, call. = FALSE)
}
data_dir <- normalizePath(data_dir)

stages <- list.files(
  file.path(project_root, "scripts", "R"),
  pattern = "^[0-9]{2}_.+[.]R$",
  full.names = TRUE
)
if (length(args) > 1) {
  if (length(args) > 2) {
    stop("Pass at most one final R stage after the data directory.", call. = FALSE)
  }
  final_stage <- match(args[[2]], basename(stages))
  if (is.na(final_stage)) {
    stop("Requested final R stage was not found: ", args[[2]], call. = FALSE)
  }
  stages <- stages[seq_len(final_stage)]
}

source(file.path(project_root, "R", "setup.R"))
old_working_directory <- setwd(data_dir)
on.exit(setwd(old_working_directory), add = TRUE)

for (stage in sort(stages)) {
  message("Running ", basename(stage))
  sys.source(stage, envir = globalenv())
}
