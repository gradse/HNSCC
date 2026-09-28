# Analysis stage from archive/HNSCC_Analysis.Rmd (source lines 3195-1000000000).
# Run through scripts/run_r_analysis.R.

# Source notebook line 3195
clinical_all_CDEs <- fread("clinical_All_CDEs.txt", sep = "\t", header = TRUE)

# Transpose the clinical data
clinical_transposed <- t(clinical_all_CDEs)

clinical_transposed <- as.data.frame(clinical_transposed)

colnames(clinical_transposed) <- clinical_transposed[1,]

clinical_transposed <- clinical_transposed[-1,]

# Select necessary variables
clinical_focused_for_age <- clinical_transposed[, c(9, 12, 21:23, 46, 48, 54, 114)]

# Filter data
clinical_focused_filtered_for_age <- clinical_focused_for_age %>%
  filter(hpv_status != "indeterminate", tumor_tissue_site == "head and neck")

clinical_focused_filtered_for_age$age_at_initial_pathologic_diagnosis <- as.numeric(clinical_focused_filtered_for_age$age_at_initial_pathologic_diagnosis)

# Source notebook line 3219
# Ensure the hpv_status values are consistent
clinical_focused_filtered_for_age$hpv_status <- ifelse(clinical_focused_filtered_for_age$hpv_status %in% c("HPV(+)", "positive"), "HPV(+)", 
                                                       ifelse(clinical_focused_filtered_for_age$hpv_status %in% c("HPV(-)", "negative"), "HPV(-)", NA))

# Calculate the summary statistics
clinical_summary <- clinical_focused_filtered_for_age %>%
  group_by(hpv_status) %>%
  summarise(
    Mean_Age = round(mean(age_at_initial_pathologic_diagnosis, na.rm = TRUE), 2),
    Median_Age = median(age_at_initial_pathologic_diagnosis, na.rm = TRUE),
    Range_Age = paste0(min(age_at_initial_pathologic_diagnosis, na.rm = TRUE), " - ", max(age_at_initial_pathologic_diagnosis, na.rm = TRUE)),
    Female_Count = sum(gender == "female", na.rm = TRUE),
    Male_Count = sum(gender == "male", na.rm = TRUE),
    Total_Count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Female_Percentage = round(Female_Count / Total_Count * 100, 2),
    Male_Percentage = round(Male_Count / Total_Count * 100, 2)
  )

# Extract subdivision counts as a data frame and calculate percentages
subdivision_counts <- clinical_focused_filtered_for_age %>%
  group_by(hpv_status, anatomic_neoplasm_subdivision) %>%
  summarise(Count = n(), .groups = "drop") %>%
  spread(hpv_status, Count, fill = 0) %>%
  mutate(
    Total = `HPV(+)` + `HPV(-)`,
    HPV_positive_percentage = round(`HPV(+)` / sum(`HPV(+)`, na.rm = TRUE) * 100, 2),
    HPV_negative_percentage = round(`HPV(-)` / sum(`HPV(-)`, na.rm = TRUE) * 100, 2)
  )

# Prepare the table content
table_content <- data.frame(
  Category = c("Age", "", "", "Gender (%)", "", "Subdivision (%)", rep("", nrow(subdivision_counts) - 1)),
  Subcategory = c("Mean", "Median", "Range", "Female", "Male", subdivision_counts$anatomic_neoplasm_subdivision),
  `HPV(+)` = c(clinical_summary$Mean_Age[clinical_summary$hpv_status == "HPV(+)"],
               clinical_summary$Median_Age[clinical_summary$hpv_status == "HPV(+)"],
               clinical_summary$Range_Age[clinical_summary$hpv_status == "HPV(+)"],
               paste0(clinical_summary$Female_Count[clinical_summary$hpv_status == "HPV(+)"], " (", clinical_summary$Female_Percentage[clinical_summary$hpv_status == "HPV(+)"], ")"),
               paste0(clinical_summary$Male_Count[clinical_summary$hpv_status == "HPV(+)"], " (", clinical_summary$Male_Percentage[clinical_summary$hpv_status == "HPV(+)"], ")"),
               paste0(as.integer(subdivision_counts$`HPV(+)`), " (", subdivision_counts$HPV_positive_percentage, ")")),
  `HPV(-)` = c(clinical_summary$Mean_Age[clinical_summary$hpv_status == "HPV(-)"],
               clinical_summary$Median_Age[clinical_summary$hpv_status == "HPV(-)"],
               clinical_summary$Range_Age[clinical_summary$hpv_status == "HPV(-)"],
               paste0(clinical_summary$Female_Count[clinical_summary$hpv_status == "HPV(-)"], " (", clinical_summary$Female_Percentage[clinical_summary$hpv_status == "HPV(-)"], ")"),
               paste0(clinical_summary$Male_Count[clinical_summary$hpv_status == "HPV(-)"], " (", clinical_summary$Male_Percentage[clinical_summary$hpv_status == "HPV(-)"], ")"),
               paste0(as.integer(subdivision_counts$`HPV(-)`), " (", subdivision_counts$HPV_negative_percentage, ")"))
)

# Display the table
table_content %>%
  kable("html", col.names = c("", "", "HPV(+)", "HPV(-)")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))

# Source notebook line 3283
# Ensure the hpv_status values are consistent
clinical_focused_filtered_for_age$hpv_status <- ifelse(clinical_focused_filtered_for_age$hpv_status %in% c("HPV(+)", "positive"), "HPV(+)", 
                                                       ifelse(clinical_focused_filtered_for_age$hpv_status %in% c("HPV(-)", "negative"), "HPV(-)", NA))

# Calculate the summary statistics
clinical_summary <- clinical_focused_filtered_for_age %>%
  group_by(hpv_status) %>%
  summarise(
    Mean_Age = round(mean(age_at_initial_pathologic_diagnosis, na.rm = TRUE), 2),
    Median_Age = round(median(age_at_initial_pathologic_diagnosis, na.rm = TRUE), 2),
    Range_Age = paste0(min(age_at_initial_pathologic_diagnosis, na.rm = TRUE), " - ", max(age_at_initial_pathologic_diagnosis, na.rm = TRUE)),
    Female_Count = sum(gender == "female", na.rm = TRUE),
    Male_Count = sum(gender == "male", na.rm = TRUE),
    Total_Count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Female_Percentage = round(Female_Count / Total_Count * 100, 2),
    Male_Percentage = round(Male_Count / Total_Count * 100, 2)
  )

# Extract subdivision counts as a data frame and calculate percentages
subdivision_counts <- clinical_focused_filtered_for_age %>%
  group_by(hpv_status, anatomic_neoplasm_subdivision) %>%
  summarise(Count = n(), .groups = "drop") %>%
  spread(hpv_status, Count, fill = 0) %>%
  mutate(
    Total = `HPV(+)` + `HPV(-)`,
    HPV_positive_percentage = round(`HPV(+)` / sum(`HPV(+)`, na.rm = TRUE) * 100, 2),
    HPV_negative_percentage = round(`HPV(-)` / sum(`HPV(-)`, na.rm = TRUE) * 100, 2)
  )

# Prepare the table content
table_content <- data.frame(
  Category = c("Age", "", "", "Sex", "", "Subdivision", rep("", nrow(subdivision_counts) - 1)),
  Subcategory = c("Mean", "Median", "Range", "Female", "Male", subdivision_counts$anatomic_neoplasm_subdivision),
  `HPV(+)` = c(clinical_summary$Mean_Age[clinical_summary$hpv_status == "HPV(+)"],
               formatC(clinical_summary$Median_Age[clinical_summary$hpv_status == "HPV(+)"], format = "f", digits = 2),
               clinical_summary$Range_Age[clinical_summary$hpv_status == "HPV(+)"],
               paste0(clinical_summary$Female_Count[clinical_summary$hpv_status == "HPV(+)"], " (", formatC(clinical_summary$Female_Percentage[clinical_summary$hpv_status == "HPV(+)"], format = "f", digits = 2), "%)"),
               paste0(clinical_summary$Male_Count[clinical_summary$hpv_status == "HPV(+)"], " (", formatC(clinical_summary$Male_Percentage[clinical_summary$hpv_status == "HPV(+)"], format = "f", digits = 2), "%)"),
               paste0(as.integer(subdivision_counts$`HPV(+)`), " (", formatC(subdivision_counts$HPV_positive_percentage, format = "f", digits = 2), "%)")),
  `HPV(-)` = c(clinical_summary$Mean_Age[clinical_summary$hpv_status == "HPV(-)"],
               formatC(clinical_summary$Median_Age[clinical_summary$hpv_status == "HPV(-)"], format = "f", digits = 2),
               clinical_summary$Range_Age[clinical_summary$hpv_status == "HPV(-)"],
               paste0(clinical_summary$Female_Count[clinical_summary$hpv_status == "HPV(-)"], " (", formatC(clinical_summary$Female_Percentage[clinical_summary$hpv_status == "HPV(-)"], format = "f", digits = 2), "%)"),
               paste0(clinical_summary$Male_Count[clinical_summary$hpv_status == "HPV(-)"], " (", formatC(clinical_summary$Male_Percentage[clinical_summary$hpv_status == "HPV(-)"], format = "f", digits = 2), "%)"),
               paste0(as.integer(subdivision_counts$`HPV(-)`), " (", formatC(subdivision_counts$HPV_negative_percentage, format = "f", digits = 2), "%)"))
)

# Display the table with improved formatting
html_table <- table_content %>%
  kable("html", col.names = c("", "", "HPV(+)", "HPV(-)"), align = c("l", "l", "r", "r")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE, position = "left") %>%
  column_spec(1, bold = TRUE, border_left = TRUE) %>%
  column_spec(2, bold = TRUE) %>%
  column_spec(3, border_right = TRUE) %>%
  column_spec(4, border_right = TRUE) %>%
  row_spec(0, bold = TRUE) %>%
  kable_classic(full_width = FALSE, html_font = "Cambria")

# Add CSS for vertical lines
html_table <- gsub('<table class="kable_table', '<table class="kable_table" style="border-collapse: collapse; border-left: 2px solid black; border-right: 2px solid black;', html_table)

html_table <- gsub('<tbody>', '<tbody style="border-left: 2px solid black; border-right: 2px solid black;">', html_table)
html_table <- gsub('<th style="', '<th style="border-left: 2px solid black; border-right: 2px solid black; ', html_table)
html_table <- gsub('<td style="', '<td style="border-left: 2px solid black; border-right: 2px solid black; ', html_table)

# Ensure no vertical border between Category and Subcategory
html_table <- gsub('<td style="border-right: 2px solid black; ', '<td style="', html_table, fixed = TRUE)

html_table

save_kable(html_table, "table1_new.html")

# Use webshot to save the HTML as a PNG with specified dpi
webshot("table1_new.html", "table1_new.png", zoom = 3, vwidth = 1000, vheight = 800)

img <- image_read("table1_new.png")

# Rescale the image to ensure it is 100 dpi
img <- image_scale(img, "1000x800")

# Save the image with 100 dpi
image_write(img, path = "table1_new_100dpi.png", density = 800)

# Read the image
img <- image_read("table1_new.png")

# Optionally, rescale the image (adjust the dimensions if necessary)
# For example, scaling to a desired size while maintaining proportions
img <- image_scale(img, "1500x1200")

# Save the image with the specified DPI
image_write(img, path = "table1_new_700dpi.png", density = 700)


# Source notebook line 3389
clinical_all_CDEs <- fread("clinical_All_CDEs.txt", sep = "\t", header = TRUE)

# Transpose the clinical data
clinical_transposed <- t(clinical_all_CDEs)

clinical_transposed <- as.data.frame(clinical_transposed)

colnames(clinical_transposed) <- clinical_transposed[1,]

clinical_transposed <- clinical_transposed[-1,]

# Select necessary variables
clinical_focused_for_death <- clinical_transposed[, c(3, 4, 9, 46, 54, 114)]

# Filter data
clinical_focused_filtered_for_death <- clinical_focused_for_death %>%
  filter(hpv_status != "indeterminate", tumor_tissue_site == "head and neck")

# Convert days_to_death to numeric
clinical_focused_filtered_for_death$days_to_death <- as.numeric(clinical_focused_filtered_for_death$days_to_death)

# Define 5-year survival (1825 days)
five_years <- 1825

# Calculate 5-year survival status
clinical_focused_filtered_for_death$survived_5_years <- ifelse(clinical_focused_filtered_for_death$vital_status == "alive" | clinical_focused_filtered_for_death$days_to_death >= five_years, TRUE, FALSE)

# Calculate 5-year survival rate for HPV positive and negative groups
hpv_positive_5yr_survival_rate <- mean(clinical_focused_filtered_for_death$survived_5_years[clinical_focused_filtered_for_death$hpv_status == "positive"], na.rm = TRUE)
hpv_negative_5yr_survival_rate <- mean(clinical_focused_filtered_for_death$survived_5_years[clinical_focused_filtered_for_death$hpv_status == "negative"], na.rm = TRUE)

# Print the survival rates
cat("5-year survival rate for HPV positive group:", hpv_positive_5yr_survival_rate, "\n")
cat("5-year survival rate for HPV negative group:", hpv_negative_5yr_survival_rate, "\n")





# Calculate survival time and event status
clinical_focused_filtered_for_death$survival_time <- ifelse(
  clinical_focused_filtered_for_death$vital_status == "alive",
  five_years,
  clinical_focused_filtered_for_death$days_to_death
)

clinical_focused_filtered_for_death$event_observed <- clinical_focused_filtered_for_death$vital_status == "dead"

# Create the survival object
surv_object <- Surv(time = clinical_focused_filtered_for_death$survival_time, event = clinical_focused_filtered_for_death$event_observed)

# Fit the Kaplan-Meier curve
fit <- survfit(surv_object ~ hpv_status, data = clinical_focused_filtered_for_death)

km_plot <- ggsurvplot(
  fit,
  data = clinical_focused_filtered_for_death,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  ggtheme = theme_minimal(),
  xlab = "Days",
  ylab = "Survival Probability",
  xlim = c(0, 2500),
)

km_plot

# Save the plot
ggsave("Kaplan_Meier_Survival_Curve_Truncated.png", plot = km_plot$plot, dpi = 300, width = 10, height = 6)

# Source notebook line 3466
# Load the image
img <- image_read("network1.jpg")

img_scaled <- image_scale(img, "1024x1024")

# Save the image
image_write(img_scaled, path = "network1_scaled_300dpi.jpg", density = 300)

