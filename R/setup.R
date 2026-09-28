# Shared setup for the behavior-preserving R analysis scripts.

hnscc_packages <- c(
  "arrow", "data.table", "dplyr", "ENmix", "ggplot2", "ggrepel",
  "grid", "gridExtra", "htmltools", "kableExtra", "knitr", "magick",
  "parallel", "pheatmap", "readr", "stringr", "survival", "survminer",
  "tidyr", "wateRmelon", "webshot", "writexl"
)

hnscc_load_packages <- function(packages = hnscc_packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      ". See README.md for installation instructions.",
      call. = FALSE
    )
  }

  invisible(lapply(packages, library, character.only = TRUE))
}

hnscc_load_packages()
knitr::opts_chunk$set(echo = FALSE, message = FALSE)
