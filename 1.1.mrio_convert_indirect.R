#### 1. Define EXIO sectoral divides and extract non-food satellites ####

# Note: This is a prep step to derive indirect energy and labor use satellites for non-food sectors only.
#       

i_exio_bio_sectors = which(colSums(prod_map) > 0) 
# Remove wool and chemicals from food sectors
i_exio_food_sectors = setdiff(i_exio_bio_sectors, c(15, 90))

# Sector names
exio_food_sectors = EXIO_reg$sector[i_exio_food_sectors]
exio_nonfood_sectors = EXIO_reg$sector[setdiff(1:200, i_exio_food_sectors)] # 174 rows
exio_bio_sectors = EXIO_reg$sector[i_exio_bio_sectors]
  
# Derive indirect energy and labor satellites (F mtx) for non-food sectors (per ton)
ene_net_nonfood = matrix(ene_net*(colSums(prod_map)==0), nrow=1)
lab_male_nonfood = matrix(lab_male*(colSums(prod_map)==0), nrow=1) 
lab_female_nonfood = matrix(lab_female*(colSums(prod_map)==0), nrow=1)

x_inv = 1/EXIO_x
x_inv[!is.finite(x_inv)] = 0

# S matrix (energy/labor intensity by $)
S_ene_net_nonfood = as.vector(ene_net_nonfood) * as.vector(x_inv)
S_lab_male_nonfood = as.vector(lab_male_nonfood) * as.vector(x_inv)
S_lab_female_nonfood = as.vector(lab_female_nonfood) * as.vector(x_inv)

# Prepare global total bio-sector FD (final demand) matrix
EXIO_Y_food = rowSums(EXIO_Y) * (colSums(prod_map)>0)
EXIO_x_food = EXIO_x * (colSums(prod_map)>0)

# Calculate total food sector footprint matrices (D) for energy and labor in non-food sectors 
# (S_nonfood * L * x_food) 8526x9800
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
# (intensity by x (mass)) 8526x9800
# total_intensity_exio_by_mass = list(
#   d_en = sweep(indir_sat_exio[["sat_en"]], 2, exio_mass_x, "/"),
#   d_hr_m = sweep(indir_sat_exio[["sat_hr_m"]], 2, exio_mass_x, "/"),
#   d_hr_f = sweep(indir_sat_exio[["sat_hr_f"]], 2, exio_mass_x, "/")
# )
total_intensity_exio_by_mass = lapply(indir_sat_exio, function(d) {
  mat = sweep(d, 2, exio_mass_x, "/")
  mat[!is.finite(mat)] = 0 # Replace NaN and Inf with 0
  return (mat)
})
names(total_intensity_exio_by_mass) <- c("en", "hr_m", "hr_f")


# The intensity (per masS) columns above can be reordered to match FABIO countries using:
# reorder_countries_to_FABIO(...)
# 8526x37400 (massive...)
total_intensity_fabio = list(
  d_en = reorder_countries_to_FABIO(total_intensity_exio_by_mass[["en"]]),
  d_hr_m = reorder_countries_to_FABIO(total_intensity_exio_by_mass[["hr_m"]]),
  d_hr_f = reorder_countries_to_FABIO(total_intensity_exio_by_mass[["hr_f"]])
)

# Each row in total_intensity_fabio can be multiplied by 
#       colSums(FABIO_y_in_EXIO) to get total energy/labor footprint for that EXIO non-food sector/origin country.
# And then, this vector can mapped to FABIO food sector/country by using p_fabio_exio.

# m1 = t(t(total_intensity_fabio[["d_en"]])*colSums(FABIO_x_in_EXIO)) #8526x37400
# sat_nonfood_FABIO = p_fabio_exio %*% t(m1) #8526x23001

l_int_i <- lapply(total_intensity_fabio, function(d) {
  FP_trans = p_fabio_exio %*% (t(d) * colSums(FABIO_x_in_EXIO))
  intensity = t(FP_trans / FABIO_x )
  intensity[!is.finite(intensity)] = 0
  return (intensity) # per ton
})
names(l_int_i) <- c("en", "hr_m", "hr_f")
saveRDS(l_int_i, file = paste0("data/FABIO_exio_satellites_nonfood_", year, ".rds"))

# This gives an intensity matrix (in the same dim as satellite (8526x23001)).



