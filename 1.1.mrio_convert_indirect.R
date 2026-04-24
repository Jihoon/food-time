#### 1. Define EXIO sectoral divides and extract non-food satellites ####

# Note: This is a prep step to derive indirect energy and labor use satellites for non-food sectors only.
#       

i_exio_bio_sectors = which(colSums(prod_map) > 0) 
# Remove some bio-sectors that are not edible.
bio_nonfood = c("Wool, silk-worm cocoons", "Chemicals nec", "Plant-based fibers")
i_exio_food_sectors = setdiff(i_exio_bio_sectors, 
                              which(EXIO_sect %in% bio_nonfood))
# Make a boolean vector for food sectors
idx_food = colSums(prod_map) > 0
idx_food[which(EXIO_sect %in% bio_nonfood)] = FALSE
idx_nonfood = !idx_food 

# Add test sum(idx_food) + sum(idx_nonfood) == 200

# Sector names
exio_food_sectors = EXIO_reg$sector[i_exio_food_sectors]
exio_nonfood_sectors = EXIO_reg$sector[setdiff(1:200, i_exio_food_sectors)] # 175 rows
exio_bio_sectors = EXIO_reg$sector[i_exio_bio_sectors]
  
# Derive indirect energy and labor satellites (F mtx) for non-food sectors (per ton)
ene_net_nonfood = matrix(ene_net*idx_nonfood, nrow=1)
lab_male_nonfood = matrix(lab_male*idx_nonfood, nrow=1) 
lab_female_nonfood = matrix(lab_female*idx_nonfood, nrow=1)

x_inv = 1/EXIO_x
x_inv[!is.finite(x_inv)] = 0

# S matrix (energy/labor intensity by $)
S_ene_net_nonfood = as.vector(ene_net_nonfood) * as.vector(x_inv)
S_lab_male_nonfood = as.vector(lab_male_nonfood) * as.vector(x_inv)
S_lab_female_nonfood = as.vector(lab_female_nonfood) * as.vector(x_inv)

# Prepare global total bio-sector FD (final demand) matrix
EXIO_Y_food = rowSums(EXIO_Y) * idx_food
EXIO_x_food = EXIO_x * idx_food

# Calculate total food sector footprint matrices (D) for energy and labor in non-food sectors 
# (S_nonfood * L * x_food) 8575x9800
# Multiply x_food to get producers' footprints
indir_sat_exio = list(
  sat_en = Diagonal(x=S_ene_net_nonfood)[EXIO_reg$sector %in% exio_nonfood_sectors,] %*% 
    EXIO_L %*% Diagonal(x=EXIO_x_food),
  sat_hr_m = Diagonal(x=S_lab_male_nonfood)[EXIO_reg$sector %in% exio_nonfood_sectors,] %*% 
    EXIO_L %*% Diagonal(x=EXIO_x_food),
  sat_hr_f = Diagonal(x=S_lab_female_nonfood)[EXIO_reg$sector %in% exio_nonfood_sectors,] %*% 
    EXIO_L %*% Diagonal(x=EXIO_x_food)
)
  
# # Convert FABIO_y to EXIO food sector mass vector
# l = convert_mass_vecs(exio_vec_monetary=matrix(EXIO_Y_food), fabio_vec_mass=rowSums(FABIO_y))
# exio_mass_y = l[[1]]  # Mass vector for EXIO food sectors (len 9800)
# FABIO_y_in_EXIO = l[[2]]  # Mass mapping (23001x37400) - 187 countries
# rm(l)

# Derive intensity-based total non-food energy/labor footprint matrices (J/kg or hr/kg) 
# (intensity by x (mass)) 8575x9800
# non-food columns are zero
# total_intensity_exio_by_mass = list(
#   d_en = sweep(indir_sat_exio[["sat_en"]], 2, exio_mass_x, "/"),
#   d_hr_m = sweep(indir_sat_exio[["sat_hr_m"]], 2, exio_mass_x, "/"),
#   d_hr_f = sweep(indir_sat_exio[["sat_hr_f"]], 2, exio_mass_x, "/")
# )
total_intensity_exio_by_mass = lapply(indir_sat_exio, function(d) {
  mat = sweep(d, 2, exio_mass_x, "/")
  mat[!is.finite(mat)] = 0 # Replace NaN and Inf with 0
  # Make mat to dgCMatrix
  mat = as(mat, "CsparseMatrix") # Convert to sparse matrix format to save memory
  return (mat)
})
names(total_intensity_exio_by_mass) <- c("en", "hr_m", "hr_f")

# Validate: Check column sums of total_intensity_exio_by_mass are zero for food sectors and non-zero for non-food sectors
total_intensity_exio_by_mass_check = lapply(total_intensity_exio_by_mass, function(d) {
  col_sums = colSums(d)
  print(summary(col_sums))
  food_sector_sums = col_sums[idx_food]
  nonfood_sector_sums = col_sums[idx_nonfood]
  return(list(food_sector_sums=food_sector_sums, nonfood_sector_sums=nonfood_sector_sums))
})

# Validate: Print column/row names of the biggest non-zero values in total_intensity_exio_by_mass to check if
# they correspond to expected non-food sectors and countries.
EXIO_reg_nonfood = EXIO_reg[EXIO_reg$sector %in% exio_nonfood_sectors,]
for (key in names(total_intensity_exio_by_mass)) {
  d = total_intensity_exio_by_mass[[key]]
  # Find the largest cell values in the matrix and print their corresponding sector and country names
  col_sums = colSums(d)
  top_indices = order(col_sums, decreasing = TRUE)[1]
  print(paste("Top food sectors for", key, ":"))
  print(EXIO_reg[top_indices,])
  
  top_rows = order(d[, top_indices], decreasing = TRUE)[1:10] # top 10 countries for the top non-food sector
  print(EXIO_reg_nonfood[top_rows,] %>% mutate(int = d[top_rows, top_indices]))
}


# The intensity (per masS) columns above can be reordered to match FABIO countries using:
# reorder_countries_to_FABIO(...)
# 8575x37400 (massive...)
total_intensity_fabio = list(
  d_en = reorder_countries_to_FABIO(total_intensity_exio_by_mass[["en"]]),
  d_hr_m = reorder_countries_to_FABIO(total_intensity_exio_by_mass[["hr_m"]]),
  d_hr_f = reorder_countries_to_FABIO(total_intensity_exio_by_mass[["hr_f"]])
)

# Each row in total_intensity_fabio can be multiplied by 
#       colSums(FABIO_y_in_EXIO) to get total energy/labor footprint for that EXIO non-food sector/origin country.
# And then, this vector can mapped to FABIO food sector/country by using p_fabio_exio.

# This gives an intensity matrix (in the same dim as satellite (32725x23001)).
l_int_i <- lapply(total_intensity_fabio, function(d) {
  FP_trans = p_fabio_exio %*% (t(d) * colSums(FABIO_x_in_EXIO))
  intensity = t(FP_trans / FABIO_x )
  intensity[!is.finite(intensity)] = 0
  intensity = reorder_countries_to_FABIO(intensity, direction=1)
  return (intensity) # per ton
})
names(l_int_i) <- c("en", "hr_m", "hr_f")
saveRDS(l_int_i, file = paste0("data/FABIO_exio_satellites_nonfood_", year, ".rds"))




