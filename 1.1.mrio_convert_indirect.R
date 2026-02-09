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

# Calculate total food sector footprint matrices (D) for energy and labor in non-food sectors 
# (S_nonfood * L * y_food) 
total_food_footprint_exio = list(
  D_en = Diagonal(x=S_ene_net_nonfood)[EXIO_reg$sector %in% exio_nonfood_sectors,] %*% 
    EXIO_L %*% Diagonal(x=EXIO_Y_food),
  D_hr_m = Diagonal(x=S_lab_male_nonfood)[EXIO_reg$sector %in% exio_nonfood_sectors,] %*% 
    EXIO_L %*% Diagonal(x=EXIO_Y_food),
  D_hr_f = Diagonal(x=S_lab_female_nonfood)[EXIO_reg$sector %in% exio_nonfood_sectors,] %*% 
    EXIO_L %*% Diagonal(x=EXIO_Y_food)
)
  
# Convert FABIO_y to EXIO food sector mass vector
l = convert_mass_vecs(exio_vec_monetary=matrix(EXIO_Y_food), fabio_vec_mass=rowSums(FABIO_y))
exio_mass_y = l[[1]]  # Mass vector for EXIO food sectors (len 9800)
FABIO_y_in_EXIO = l[[2]]  # Mass mapping (23001x37400) - 187 countries
rm(l)

# Derive intensity-based total non-food energy/labor footprint matrices (J/kg or hr/kg) 
# (intensity by y (mass))
total_intensity_exio_by_mass = list(
  d_en = sweep(total_food_footprint_exio[["D_en"]], 2, exio_mass_y, "/"),
  d_hr_m = sweep(total_food_footprint_exio[["D_hr_m"]], 2, exio_mass_y, "/"),
  d_hr_f = sweep(total_food_footprint_exio[["D_hr_f"]], 2, exio_mass_y, "/")
)

# The intensity (per masS) columns above can be reordered to match FABIO countries using:
# reorder_countries_to_FABIO(...)
# n_col = 37400
total_intensity_fabio = list(
  d_en = reorder_countries_to_FABIO(total_intensity_exio_by_mass[["d_en"]]),
  d_hr_m = reorder_countries_to_FABIO(total_intensity_exio_by_mass[["d_hr_m"]]),
  d_hr_f = reorder_countries_to_FABIO(total_intensity_exio_by_mass[["d_hr_f"]])
)

# Each row in total_intensity_fabio can be multiplied by 
#       colSums(FABIO_y_in_EXIO) to get total energy/labor footprint for that EXIO non-food sector/origin country.
# And then, this vector can mapped to FABIO food sector/country by using p_fabio_exio.

# m1 = t(t(total_intensity_fabio[["d_en"]])*colSums(FABIO_y_in_EXIO)) #8526x37400
# sat_nonfood_FABIO = p_fabio_exio %*% t(m1) #8526x23001

nonfood_int_FABIO <- lapply(total_intensity_fabio, function(d) {
  FP_trans = p_fabio_exio %*% (t(d) * colSums(FABIO_y_in_EXIO))
  intensity = t(FP_trans / FABIO_x )
  return (intensity) # per ton
})
names(nonfood_int_FABIO) <- c("en", "hr_m", "hr_f")

# This gives an intensity matrix (in the same dim as satellite (8526x23001)).
# Then multiply L*y to get FABIO footprint for non-food sectors.


#### Footprinting analysis ####

# Do nonfood_int_FABIO[[i]] %*% FABIO_L %*% FABIO_y_hh to get total non-food footprint for household consumption by FABIO country/sector.







# Obsolete

# # M_ene_net_nonfood = Diagonal(x=S_ene_net_nonfood) %*% EXIO_L
# # M_lab_male_nonfood = Diagonal(x=S_lab_male_nonfood) %*% EXIO_L
# # M_lab_female_nonfood = Diagonal(x=S_lab_female_nonfood) %*% EXIO_L
# 
# indir_sat_exio = list(
#   sat_en = Diagonal(x=S_ene_net_nonfood) %*% EXIO_L,
#   sat_hr_m = Diagonal(x=S_lab_male_nonfood) %*% EXIO_L,
#   sat_hr_f = Diagonal(x=S_lab_female_nonfood) %*% EXIO_L) 
# 
# 
# # Derive D matrices (footprint based on $ expenditure) for household consumption by region and sector
# D_list = list()
# for (s in names(indir_sat_exio)) {
#   D_list_satellite = list()
#   for (r in 1:n_reg_EXIO) {
#     # Final hh food expenditure
#     EXIO_Y_hh_region = EXIO_Y[,((r-1)*7+1)]
#     EXIO_Y_hh_region[EXIO_reg$sector %in% exio_food_sectors == FALSE] = 0
#     
#     # Calculate non-food footprint matrix for the region (ignore target country here. i.e. care only about consumer country)
#     D_hh_region = indir_sat_exio[[s]] %*% Diagonal(x=EXIO_Y_hh_region)
#     D_hh_region = rowsum(t(as.matrix(D_hh_region)), group = EXIO_reg$sector, reorder=FALSE)
#     D_hh_region = t(D_hh_region) # 9800x200
#     
#     D_list_satellite[[r]] = D_hh_region
#   }
#   D_list[[s]] = do.call(cbind, D_list_satellite)
#   print(paste("Completed satellite:", s))
#   print(head(colSums(D_list[[s]]))) # This has non-zero values only for food sectors.
# }
# 
# # Divide each D_list matrix by EXIO_mass (derived from FABIO Y_food) to get intensity-based D matrices
# FABIO_exio_y = Diagonal(x=FABIO_x) %*% p_fabio_exio
#   