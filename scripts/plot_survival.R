#!/usr/bin/env Rscript

# Generate the corrected Kaplan-Meier plot using only base R and survival.

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
project_root <- dirname(dirname(script_path))

input_file <- if (length(args)) args[[1]] else file.path(project_root, "data", "local", "clinical_All_CDEs.txt")
output_file <- if (length(args) > 1) args[[2]] else file.path(project_root, "results", "figures", "Kaplan_Meier_Survival_Curve_Corrected.png")
max_days <- if (length(args) > 2) as.numeric(args[[3]]) else 2500

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The R package 'survival' is required.", call. = FALSE)
}

clinical <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
clinical <- as.data.frame(t(clinical), stringsAsFactors = FALSE)
colnames(clinical) <- clinical[1, ]
clinical <- clinical[-1, , drop = FALSE]

survival_data <- clinical[
  clinical$hpv_status %in% c("negative", "positive") &
    clinical$tumor_tissue_site == "head and neck",
  c("hpv_status", "vital_status", "days_to_death", "days_to_last_followup"),
  drop = FALSE
]
survival_data$days_to_death <- suppressWarnings(as.numeric(survival_data$days_to_death))
survival_data$days_to_last_followup <- suppressWarnings(as.numeric(survival_data$days_to_last_followup))
survival_data$time <- ifelse(
  survival_data$vital_status == "dead",
  survival_data$days_to_death,
  survival_data$days_to_last_followup
)
survival_data$event <- survival_data$vital_status == "dead"
survival_data <- survival_data[
  is.finite(survival_data$time) & survival_data$time >= 0,
  ,
  drop = FALSE
]

fit <- survival::survfit(
  survival::Surv(time, event) ~ hpv_status,
  data = survival_data,
  conf.type = "log"
)
comparison <- survival::survdiff(
  survival::Surv(time, event) ~ hpv_status,
  data = survival_data
)
p_value <- pchisq(comparison$chisq, df = length(comparison$n) - 1, lower.tail = FALSE)

step_coordinates <- function(times, values, end_time) {
  keep <- times <= end_time
  times <- times[keep]
  values <- values[keep]
  if (!length(times)) return(list(x = c(0, end_time), y = c(1, 1)))
  before <- c(1, head(values, -1))
  list(
    x = c(0, rep(times, each = 2), end_time),
    y = c(1, as.vector(rbind(before, values)), tail(values, 1))
  )
}

strata_sizes <- fit$strata
strata_ends <- cumsum(strata_sizes)
strata_starts <- c(1, head(strata_ends, -1) + 1)
colors <- c("hpv_status=negative" = "#F8766D", "hpv_status=positive" = "#00BFC4")

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
png(output_file, width = 3000, height = 1800, res = 300, bg = "white")
par(mar = c(5.2, 5.5, 2.8, 1.5), mgp = c(3.2, 0.9, 0), las = 1)
plot(
  NA,
  xlim = c(0, max_days), ylim = c(0, 1.05),
  xlab = "Days", ylab = "Survival Probability",
  axes = FALSE, cex.lab = 1.25
)
abline(h = seq(0, 1, 0.125), v = seq(0, max_days, 500), col = "#EBEBEB", lwd = 1)
axis(1, at = seq(0, max_days, 500))
axis(2, at = seq(0, 1, 0.25), labels = sprintf("%.2f", seq(0, 1, 0.25)))
box(bty = "l", col = "#777777")

for (index in seq_along(strata_sizes)) {
  rows <- strata_starts[index]:strata_ends[index]
  label <- names(strata_sizes)[index]
  color <- colors[[label]]
  lower <- step_coordinates(fit$time[rows], fit$lower[rows], max_days)
  upper <- step_coordinates(fit$time[rows], fit$upper[rows], max_days)
  curve <- step_coordinates(fit$time[rows], fit$surv[rows], max_days)

  polygon(
    c(lower$x, rev(upper$x)), c(lower$y, rev(upper$y)),
    border = NA, col = adjustcolor(color, alpha.f = 0.25)
  )
  lines(curve$x, curve$y, col = color, lwd = 3)

  censor_rows <- rows[fit$n.censor[rows] > 0 & fit$time[rows] <= max_days]
  if (length(censor_rows)) {
    points(fit$time[censor_rows], fit$surv[censor_rows], pch = 3, col = color, cex = 0.8, lwd = 1.5)
  }
}

legend(
  "top", inset = c(0, -0.01), horiz = TRUE, bty = "n",
  legend = c("hpv_status=negative", "hpv_status=positive"),
  col = unname(colors), lwd = 3, pch = 3, title = "Strata", cex = 1.05
)
text(max_days * 0.05, 0.20, labels = paste0("p = ", format.pval(p_value, digits = 2)), adj = 0, cex = 1.35)
dev.off()

five_year <- summary(fit, times = 1825, extend = TRUE)
print(data.frame(
  hpv_status = sub("hpv_status=", "", five_year$strata),
  survival_5_years = five_year$surv,
  lower_95 = five_year$lower,
  upper_95 = five_year$upper
))
cat("Log-rank p-value:", format.pval(p_value, digits = 4), "\n")
cat("Plot written to:", normalizePath(output_file), "\n")
