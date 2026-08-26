library(tidyverse)
library(Matrix)


#### 1. Conversions to FABIO mass vector (123) to EXIO classifications (200, still mass) ####

# Prepare a conversion mapping (qualitative) from FABIO to EXIO sector mapping (for 187 countries)
p_fabio_exio = block_diag_repeat(as.matrix(prod_map), nrow(FABIO_reg)) 
  
# Debugging: Find the indices of FABIO/EXIO sectors that are considered fish sectors 
# Find fish-related FABIO rows among the 23,001 FABIO sectors/countries
# (regex is what Claude suggested, but "fish" is enough actually.)
fish_idx <- grep("fish|aqua|crustac|mollusc|seafood", io$item, 
                ignore.case = TRUE)
# Make a named vector for the 37,400 entries (187x200)
long_idx = (expand.grid(prod = EXIO_sect, ISO=regions$iso3c) %>%
  mutate(ISO_prod = paste0(ISO, "_", prod)))$ISO_prod
fish_idx_l <- grep("fish|aqua|crustac|mollusc|seafood", long_idx, 
                ignore.case = TRUE)

df_fish = io %>% mutate(tonne = FABIO_x) %>% slice(fish_idx) %>%
  left_join(pop_y %>% select(iso3c, pop), by = c("iso3c")) 

convert_mass_vecs <- function (exio_vec_monetary=EXIO_x, fabio_vec_mass=FABIO_x) {
  # Convert mass vectors from FABIO classification to EXIO classification, 
  #   based on the mapping p_fabio_exio
  
  # exio_vec_monetary: EXIO monetary vector (len 9800)
  # fabio_vec_mass: FABIO mass vector (len 23001)
  # exio_vec_monetary is used to derive row-wise shares for mapping mtx.
  # Then, fabio_vec_mass is multiplied row-wise to get a vector in mass unit 
  # under EXIO classification using the mapping mtx.
  
  # CHECK: Deriving the row (monetary) shares based on total x_i is just one way.
  #        Not necessarily good for rows with multiple cells with 1
  #        And any better ways?
  
  # Normalize mat_map by (monetary) row sums to get a proper mapping matrix
  mat_map = p_fabio_exio %*% Diagonal(x=reorder_countries_to_FABIO(exio_vec_monetary)) 
  p_inv = 1/Matrix::rowSums(mat_map)
  p_inv[!is.finite(p_inv)] = 0
  
  # mat_map has row-wise shares of monetary flows by sector.
  mat_map = Diagonal(x=p_inv) %*% mat_map
  
  # Multiply diag(FABIO_x) and mat_map to get a mass vector in EXIO sector classification
  FABIO_mass_in_EXIO = Diagonal(x=fabio_vec_mass) %*% mat_map
  v_mass = colSums(FABIO_mass_in_EXIO) # mass vec in EXIO sectors (187 countries = len 37400)
  
  # Aggregate FABIO regions to EXIO regions and reorder
  exio_mass = reorder_countries_to_EXIO(v_mass) # mass vec in EXIO sec/reg (9800)

  ### Diagnostic checks

  # Check the mass split for fish rows
  FABIO_mass_in_EXIO_fish <- FABIO_mass_in_EXIO[fish_idx, ]
  # Where does fish mass land in EXIO?
  fish_exio_sums <- colSums(FABIO_mass_in_EXIO_fish)
  top_fish_exio <- sort(fish_exio_sums, decreasing = TRUE)[1:20]

  # idx contains the original positions
  idx <- order(fish_exio_sums, decreasing = TRUE)[1:20]

  cat("Top EXIO sectors receiving fish mass:\n")
  print(data.frame(sector = long_idx[idx], mass = fish_exio_sums[idx]))

  return(list(exio_mass, FABIO_mass_in_EXIO))
}

l = convert_mass_vecs()
exio_mass_x = l[[1]]  # EXIO mass vector (len 9800)
FABIO_x_in_EXIO = l[[2]]  # EXIO mass vector (23001x37400)
rm(l)

# TEST: Validate sum(FABIO_x)-sum(FABIO_x[rowSums(p_fabio_exio)==0]) == sum(exio_mass_x) == sum(v_mass)




#### 2. Convert EXIO energy/labor satellites (and intensities) to FABIO classification ####

# Make EXIO satellite accounts in the same formatting 
dir_sat_exio = list(
  sat_en = ene_net,
  sat_hr_m = lab_male,
  sat_hr_f = lab_female) 

# Get direct energy/labor (hr) satellites in FABIO sectors
# dim = (23001, 37400), which can be row/column-summed to get vectors.
dir_sat_FAB = convert_intensities(dir_sat_exio)

# TEST: sum(dir_sat_FAB[[2]]) ~ sum(matrix(lab_male*idx_food, nrow=1)) 
# TEST: sum(dir_sat_FAB[[3]]) ~ sum(matrix(lab_female*idx_food, nrow=1)) 


# Save FABIO-based intensities (per tonne of product) by summing the satellite accounts and dividing by FABIO_x (mass vector in FABIO classification))
l_int_d <- lapply(dir_sat_FAB, function(d) {
  v = rowSums(d) / FABIO_x 
  v[!is.finite(v)] = 0
  return (v) 
})
names(l_int_d) <- c("en", "hr_m", "hr_f")

# # # Alternative way: directly calculate intensities by summing the satellite accounts and dividing by FABIO_x (mass vector in FABIO classification))
# # FABIO satellite account totals (intensities per tonne of product)
# FABIO_en_int_d = rowSums(dir_sat_FAB$sat_en_FAB) / FABIO_x 
# FABIO_hr_m_int_d = rowSums(dir_sat_FAB$sat_hr_m_FAB) / FABIO_x 
# FABIO_hr_f_int_d = rowSums(dir_sat_FAB$sat_hr_f_FAB) / FABIO_x 
# 
# # Replace NaN and Inf with 0
# FABIO_en_int_d[!is.finite(FABIO_en_int_d)] = 0
# FABIO_hr_m_int_d[!is.finite(FABIO_hr_m_int_d)] = 0
# FABIO_hr_f_int_d[!is.finite(FABIO_hr_f_int_d)] = 0
# 
# # Save FABIO-based EXIO satellite accounts and intensities
# l_int = list(
#   en_int_d = FABIO_en_int_d,
#   hr_m_int_d = FABIO_hr_m_int_d,
#   hr_f_int_d = FABIO_hr_f_int_d
# )

saveRDS(l_int_d, file = paste0("data/FABIO_exio_satellites_food_", year, ".rds"))



#### 3. Calculate consumption-based calorie/protein footprints ####

# Calculate consumption-based matrices
# Obs: sum(FABIO_y_hh)/sum(FABIO_x_hh) gives only 5.6%.
# Similarly, sum(FABIO_y)/sum(FABIO_x) gives only 5.9%.
# => Makes sense because x includes all production, double-counting intermediate flows.

# Import nutrient coefficients by FABIO products
coeff_cal <- fread(file.path(FABIO_path,"nutrient_coefficients.csv")) %>%
  select(kcal_per_kg) %>%
  # concatenate the vector 187 times 
  slice(rep(1:n(), times = nrow(FABIO_reg))) %>%
  # set kcal_per_kg=0 if the product is "Live animals" (comm_group == "Live animals")
  # 'Live animals' need to be removed from nutrient flow calculation because they are not consumed as such.
  mutate(kcal_per_kg = ifelse(items$comm_group[match(io$item, items$item)] == "Live animals", 0, kcal_per_kg))
  
coeff_pro <- fread(file.path(FABIO_path,"nutrient_coefficients_protein.csv")) %>%
  mutate(protein_per_kg = 10*protein_per_100g) %>% select(protein_per_kg) %>%
  # concatenate the vector 187 times 
  slice(rep(1:n(), times = nrow(FABIO_reg))) %>%
  # set protein_per_kg=0 if the product is "Live animals" (comm_group == "Live animals")
  mutate(protein_per_kg = ifelse(items$comm_group[match(io$item, items$item)] == "Live animals", 0, protein_per_kg))


# Consumption total
FABIO_y_hh_cal <- Matrix::Diagonal(x=coeff_cal$kcal_per_kg) %*% FABIO_y_hh * 1000 # kcal
FABIO_y_hh_pro <- Matrix::Diagonal(x=coeff_pro$protein_per_kg) %*% FABIO_y_hh * 1000 # g

# Production total
FABIO_x_hh_cal <- Matrix::Diagonal(x=coeff_cal$kcal_per_kg) %*% FABIO_x_hh * 1000 # kcal
FABIO_x_hh_pro <- Matrix::Diagonal(x=coeff_pro$protein_per_kg) %*% FABIO_x_hh * 1000 # g

# tend to have larger values for poorer ctys; Smaller values for richer (and smaller) ctys.
# might be because LICs are relying more on local supply chain, while HICs are more integrated to global chain.
df_stat = data.frame(
  country = regions$iso3c,
  total_cal_eff = colSums(FABIO_y_hh_cal)/colSums(FABIO_x_hh_cal),
  total_protein_eff = colSums(FABIO_y_hh_pro)/colSums(FABIO_x_hh_pro)
)




  