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

