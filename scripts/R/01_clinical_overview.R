# Analysis stage from archive/HNSCC_Analysis.Rmd (source lines 89-226).
# Run through scripts/run_r_analysis.R.

# Source notebook line 89
# Load data for clinical, methylation, and mRNA 
clinical_all_CDEs <- fread("clinical_All_CDEs.txt", sep = "\t", header = TRUE)

# Open the Parquet file and read it into a dataframe (methylation data)
methyl_level3 <- read_parquet("methylation_level3_within_bioassay.parquet")

# Source notebook line 98
# Transpose the clinical data
clinical_transposed <- t(clinical_all_CDEs)

clinical_transposed <- as.data.frame(clinical_transposed)

colnames(clinical_transposed) <- clinical_transposed[1,]

clinical_transposed <- clinical_transposed[-1,]

# Select necessary variables
clinical_focused <- clinical_transposed[, c(9, 12, 21:23, 46, 48, 54, 114)]

# Save the cleaned dataframe to a new Parquet file
write_parquet(clinical_focused, "clinical_focused.parquet")

# Filter data
clinical_focused_filtered <- clinical_focused %>%
  filter(!is.na(clinical_m), !is.na(clinical_n), !is.na(clinical_stage), hpv_status != "indeterminate", tumor_tissue_site == "head and neck")

# Source notebook line 120
clinical_export <- clinical_focused[, c(1, 2, 6, 8, 9)]

write.csv(clinical_export, file = "clinical_data.csv", row.names = FALSE)

# Source notebook line 131
subdivision <- as.factor(clinical_focused_filtered$anatomic_neoplasm_subdivision)

is_hpv <- as.factor(clinical_focused_filtered$hpv_status)

table1 <- table(subdivision, is_hpv)

table1_margin <- addmargins(table1)

kable(table1_margin, format = "html",
      caption = "Cross-Tabulation of Anatomic Neoplasm Subdivision and HPV Status") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                    full_width = FALSE)

# Source notebook line 147
# Distant metastasis (M)
# MX: Metastasis cannot be measured.
# M0: Cancer has not spread to other parts of the body.
# M1: Cancer has spread to other parts of the body.

clinical_m <- as.factor(clinical_focused_filtered$clinical_m)

table2 <- table(clinical_m, is_hpv)

table2_margin <- addmargins(table2)

kable(table2_margin, format = "html", 
      caption = "Cross-Tabulation of clinical_m and HPV Status") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                    full_width = FALSE)

# Source notebook line 166
# Regional lymph nodes (N)
# NX: Cancer in nearby lymph nodes cannot be measured.
# N0: There is no cancer in nearby lymph nodes.
# N1, N2, N3: Refers to the number and location of lymph nodes that contain cancer. The higher the number after the N, the more lymph nodes that contain cancer.

clinical_n <- as.factor(clinical_focused_filtered$clinical_n)

table3 <- table(clinical_n, is_hpv)

table3_margin <- addmargins(table3)

kable(table3_margin, format = "html", 
      caption = "Cross-Tabulation of clinical_n and HPV Status") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                    full_width = FALSE)

# Source notebook line 185
clinical_stg <- as.factor(clinical_focused_filtered$clinical_stage)

table4 <- table(clinical_stg, is_hpv)

table4_margin <- addmargins(table4)

kable(table4_margin, format = "html", 
      caption = "Cross-Tabulation of clinical_stage and HPV Status") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                    full_width = FALSE)

# Source notebook line 199
gender <- as.factor(clinical_focused_filtered$gender)

table5 <- table(gender, is_hpv)

table5_margin <- addmargins(table5)

kable(table5_margin, format = "html", 
      caption = "Cross-Tabulation of gender and HPV Status") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                    full_width = FALSE)

# Source notebook line 213
histological_type <- as.factor(clinical_focused_filtered$histological_type)

table6 <- table(histological_type, is_hpv)

table6_margin <- addmargins(table6)

kable(table6_margin, format = "html", 
      caption = "Cross-Tabulation of histological_type and HPV Status") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                    full_width = FALSE)

