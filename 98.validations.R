# After 0.mrio_prep.R + 1.mrio_convert.R:
# fish_exio_idx <- grep("fish|aqua", EXIO_sect, ignore.case = TRUE)

# For each RoW region, what share of mass comes from the outlier countries?
outliers <- c("ISL", "NZL", "NAM", "KIR")
oi <- FABIO_reg$ISO %in% outliers

for (r in unique(FABIO_reg$EXIOBASE_code[oi])) {
  countries_in_r <- FABIO_reg$ISO[FABIO_reg$EXIOBASE_code == r]

  cat(r, "-> top countries by fish mass:\n")
  df_fish %>% filter(iso3c %in% countries_in_r) %>% 
    select(iso3c, area, tonne, pop) %>%
    arrange(desc(tonne)) %>%
    print()


#   mass_by_country <- colSums(FABIO_x_in_EXIO[FABIO_reg$EXIOBASE_code == r, fish_idx_l])
# #   mass_by_country <- FABIO_x[fish_idx]
#   cat(r, "-> top countries by fish mass:\n")
#   print(sort(mass_by_country, decreasing=TRUE)[1:5])
}


# After fp_food and fp_nonfood are derived, we can do some validations for fishing:

# Make an index for FABIO region - EXIO nonfood sectors
io_nonfood = expand.grid(sect = exio_nonfood_sectors, iso3c = regions$iso3c) %>% select(iso3c, sect)

# Validate: find the biggest cell from fp_food[[2]] and check if it corresponds to expected source and destination
largest_cell_index = which(fp_food[[2]] == max(fp_food[[2]]), arr.ind = TRUE)
print(paste("Largest cell in fp_food[[2]] is at row", largest_cell_index[1], "and column", largest_cell_index[2]))
print(paste("This corresponds to country", io[largest_cell_index[1],c("iso3c", "item")], "and sector", regions$area[largest_cell_index[2]]))

largest_cell_index = which(fp_nonfood[[2]] == max(fp_nonfood[[2]]), arr.ind = TRUE)
print(paste("Largest cell in fp_food[[2]] is at row", largest_cell_index[1], "and column", largest_cell_index[2]))
print(paste0("This corresponds to country ", io_nonfood[largest_cell_index[1],"iso3c"], "'s '", 
             io_nonfood[largest_cell_index[1],"sect"], 
             "' going to ", regions$area[largest_cell_index[2]]))

# Debug: Find a country block for KIR and get the ten biggest cell values in fp_food[[2]] for that block and show their product and destination countries
kir_block_indices = which(io$iso3c == "KIR")
kir_block_values = fp_food[[2]][kir_block_indices, ]
kir_block_largest_indices = order(kir_block_values, decreasing = TRUE)[1:10]
print("Top 10 largest cells in KIR block of fp_food[[2]]:")
for (index in kir_block_largest_indices) {
  row_col = arrayInd(index, dim(kir_block_values))
  print(paste("Cell at row", row_col[1], "and column", row_col[2], "with value", kir_block_values[index]))
  print(paste("This corresponds to product", io[kir_block_indices[row_col[1]], c("iso3c", "item")], "and destination country", regions$area[row_col[2]]))
}

# Debug: Find a country block for ISL and get the ten biggest cell values in fp_food[[2]] for that block and show their product and destination countries
isl_block_indices = which(io$iso3c == "ISL")
isl_block_values = fp_food[[2]][isl_block_indices, ]
isl_block_largest_indices = order(isl_block_values, decreasing = TRUE)[1:10]
print("Top 10 largest cells in ISL block of fp_food[[2]]:")
for (index in isl_block_largest_indices) {
  row_col = arrayInd(index, dim(isl_block_values))
  print(paste("Cell at row", row_col[1], "and column", row_col[2], "with value", isl_block_values[index]))
  print(paste("This corresponds to product", io[isl_block_indices[row_col[1]], c("iso3c", "item")], "and destination country", regions$area[row_col[2]]))
}

# Debug: Find a country block for KIR and get the ten biggest cell values in fp_nonfood[[2]] for that block and show their product and destination countries
kir_block_indices = which(io_nonfood$iso3c == "KIR")
kir_block_values = fp_nonfood[[2]][kir_block_indices, ]
kir_block_largest_indices = order(kir_block_values, decreasing = TRUE)[1:10]
print("Top 10 largest cells in KIR block of fp_nonfood[[2]]:")
for (index in kir_block_largest_indices) {
  row_col = arrayInd(index, dim(kir_block_values))
  print(paste("Cell at row", row_col[1], "and column", row_col[2], "with value", kir_block_values[index]))
  print(paste("This corresponds to sector", io_nonfood[kir_block_indices[row_col[1]], c("sect")], "and destination country", regions$area[row_col[2]]))
}

# Nutrient supply validation: FAOSTAT vs FABIO ####
library(ggrepel)

# Load FAOSTAT nutrient file (contains both kcal [4003] and protein [4004])
fao_raw = read.csv("data/FAOSTAT_data_en_5-28-2026 protein.csv", check.names = FALSE)
fao_m49_to_iso3c = fao_raw %>%
  distinct(`Area Code (M49)`, Area) %>%
  mutate(iso3c = countrycode::countrycode(`Area Code (M49)`,
                                          origin = "un", destination = "iso3c",
                                          warn = FALSE))

fao_by_indicator = fao_raw %>%
  left_join(fao_m49_to_iso3c, by = c("Area Code (M49)", "Area")) %>%
  filter(!is.na(iso3c)) %>%
  group_by(iso3c, Area, `Indicator Code`, Unit) %>%
  summarise(value_cap_d = sum(Value, na.rm = TRUE), .groups = "drop")

validate_nutrient = function(fao_indicator_code, fao_col, fabio_mat, fabio_col, unit_label) {
  fao = fao_by_indicator %>%
    filter(`Indicator Code` == fao_indicator_code) %>%
    select(iso3c, Area, fao_val = value_cap_d)

  fabio = country_summary(agg_country_footprint(fabio_mat)) %>%
    mutate(fabio_val = (domestic_per_capita + import_per_capita) / 1e6 / 365)

  cmp = fao %>%
    inner_join(fabio %>% select(country, fabio_val), by = c("iso3c" = "country")) %>%
    left_join(regions %>% select(iso3c, continent), by = "iso3c") %>%
    mutate(diff_pct = (fabio_val - fao_val) / fao_val * 100)

  cat(sprintf("\n--- %s ---\n", unit_label))
  cat(sprintf("Countries matched: %d\n", nrow(cmp)))
  cat(sprintf("Correlation:       %.3f\n", cor(cmp$fao_val, cmp$fabio_val, use = "complete")))
  cat(sprintf("Mean %% diff:      %.1f%%\n", mean(cmp$diff_pct, na.rm = TRUE)))
  cat(sprintf("RMSE:              %.2f %s\n",
              sqrt(mean((cmp$fabio_val - cmp$fao_val)^2, na.rm = TRUE)), unit_label))
  cat("Top 15 outliers by |% diff|:\n")
  print(cmp %>% arrange(desc(abs(diff_pct))) %>%
          select(iso3c, Area, fao_val, fabio_val, diff_pct) %>% head(15))

  ggplot(cmp, aes(x = fao_val, y = fabio_val, color = continent, label = iso3c)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.8, size = 2) +
    ggrepel::geom_text_repel(size = 2.5, max.overlaps = 20, show.legend = FALSE) +
    scale_color_brewer(palette = "Set2") +
    labs(x = paste0("FAOSTAT ", unit_label),
         y = paste0("FABIO ", unit_label, " (domestic + import)"),
         color = "Continent",
         title = paste0(unit_label, ": FAOSTAT vs FABIO (", year, ")"),
         caption = "Dashed line = 1:1") +
    theme_minimal() +
    theme(legend.position = "right")
}

p_kcal_validation = validate_nutrient("4003", "kcal_cap_d", FABIO_y_hh_cal,
                                      "fabio_kcal_cap_d", "kcal/cap/day")
print(p_kcal_validation)
ggsave("results/validation_kcal_faostat_vs_fabio.pdf",
       p_kcal_validation, width = 10, height = 8)

p_pro_validation = validate_nutrient("4004", "pro_cap_d", FABIO_y_hh_pro,
                                     "fabio_pro_cap_d", "g protein/cap/day")
print(p_pro_validation)
ggsave("results/validation_protein_faostat_vs_fabio.pdf",
       p_pro_validation, width = 10, height = 8)


# Check mass consistency
sum(FABIO_x) == sum(exio_mass_x) == sum(FABIO_x_in_EXIO)
# about 15% loss potentially due to unmapped mass for Fodder crops & Grazing

# Total final energy in non-food sectors (TJ)  
sum(indir_sat_exio[[1]]) #== sum(p_fabio_exio %*% (t(total_intensity_fabio[[1]]) * colSums(FABIO_x_in_EXIO)))
== sum(FABIO_x_in_EXIO %*% t(total_intensity_fabio[[1]]) ) #== sum(FABIO_x * total_intensity_fabio[[1]] %*% t(p_fabio_exio))
