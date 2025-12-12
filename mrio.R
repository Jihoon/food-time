library(mrio)
library(tidyverse)
library(Matrix)

source("exiobase.R")

year = 2020  # Both labor and energy satellites available
type = "pxp"
n_col <- ifelse(type == "ixi", 7989, 9802)

EXIO_path = paste0("H:/MyDocuments/Data/EXIOBASE3/IOT_", year, "_", type)
EXIO_path_e = paste0(EXIO_path, "/energy/F.txt")
EXIO_path_l = paste0(EXIO_path, "/employment/F.txt")

FABIO_path = "H:/MyDocuments/Data/FABIO/input/"

# Read EXIO energy and labor satellite accounts
ene = as.matrix(data.table::fread(EXIO_path_e, select = 2:(n_col-1), skip = 3,
                                  header = F))
lab = as.matrix(data.table::fread(EXIO_path_l, select = 2:(n_col-1), skip = 3,
                                  header = F))
EXIO_x = as.matrix(data.table::fread(paste0(EXIO_path, "/x.txt"), select = 3, 
                                     header = T))

# Take the 'net energy'
# Note: "A net energy use account includes the energy content of all energy products 
# before being combusted or used for non-energy purposes, as well as the pseudo renewable energy 
# (before it is transformed into electricity or heat). So, it is not the output of 
# the electricity and heat generating industries, but only the input to it. 
# In other words, it is the final energy use, but with the inputs into heat and 
# electricity generation instead of the outputs." - Rasul et al. (2024) 
ene_net = matrix(ene[4,])

# Take the total male/female labor hours
lab_male = matrix(colSums(lab[c(7,9,11),]))
lab_female = matrix(colSums(lab[c(8,10,12),]))

# # Derive EXIO energy and labor intensity (by $) -> Not compatible with FABIO
# exio_en_int = ene_net / EXIO_x
# exio_hr_m_int = lab_male / EXIO_x
# exio_hr_f_int = lab_female / EXIO_x
# # Replace NaN and Inf with 0
# exio_en_int[!is.finite(exio_en_int)] = 0
# exio_hr_m_int[!is.finite(exio_hr_m_int)] = 0
# exio_hr_f_int[!is.finite(exio_hr_f_int)] = 0



# read "product_concordance" sheet from "fabio-exiobase-v2.xlsx" under 'data' folder
prod_map = readxl::read_xlsx("data/fabio-exiobase.xlsx", 
                             sheet = "product_concordance", 
                             range = "E4:GV126",
                             col_names = FALSE) %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) # Replace NAs with 0s

FABIO_x = readRDS(file.path(FABIO_path,"X.rds"))[,as.character(year)]

reg_map = readxl::read_xlsx("data/fabio-exiobase.xlsx", 
                            sheet = "regions_concordance", 
                            col_names = TRUE) 

FABIO_reg = readxl::read_xlsx(paste0(FABIO_path, "fabio_classifications_v2.xlsx"), 
                              sheet = "Countries") %>% select(-area) %>%
  rename(ISO = `iso3c`, FAO_code = `area_code`) %>%
  left_join(reg_map) %>%
  # replace NA cells with values from Country=="RoW Europe" (Lichtenstein, Monaco, Andorra, San Marino etc.)
  mutate(EXIOBASE_code = as.numeric(ifelse(EXIOBASE_code=="NA",
                                           47, EXIOBASE_code)),
         EXIOBASE = ifelse(EXIOBASE=="NA",
                           "RoW Europe", EXIOBASE))
# Note: RoW countries' mean GDP/cap is close to that of Italy (<- Perplexity)

#   # replace NA cells with values from Country=="Italy"
#   mutate(EXIOBASE_code = as.numeric(ifelse(EXIOBASE_code=="NA",
#                                 EXIOBASE_code[Country=="Italy"], EXIOBASE_code)),
#          EXIOBASE = ifelse(EXIOBASE=="NA",
#                            EXIOBASE[Country=="Italy"], EXIOBASE))
# # Note: RoW countries' mean GDP/cap is close to that of Italy (<- Perplexity)



# Re-order the EXIO intensity vectors (9800) to match FABIO region order based on FABIO_reg mapping
reorder_countries <- function(int_vector) {
  n_sect = 200
  n_countries = 49
  
  # # Make an index vector by repeating 200 times
  # exio_reg_idx = rep(1:n_countries, each=n_sect)
  
  M = matrix(int_vector, nrow = n_countries, ncol = n_sect, byrow = TRUE)
  M_reordered <- M[FABIO_reg$EXIOBASE_code, , drop = FALSE]
  
  # Flatten back to original vector shape
  data_reordered <- as.vector(t(M_reordered))
  
  return(data_reordered)
}

# exio_en_int_ord = reorder_countries(exio_en_int)
# exio_hr_m_int_ord = reorder_countries(exio_hr_m_int)
# exio_hr_f_int_ord = reorder_countries(exio_hr_f_int)

# Re-order the FABIO mass X vectors to match EXIO region order based on FABIO_reg mapping
reorder_countries_mass <- function(m_vector) {
  n_sect = 200
  n_countries = 187
  
  # # Make an index vector by repeating 200 times
  # exio_reg_idx = rep(1:n_countries, each=n_sect)
  r = FABIO_reg %>% arrange(EXIOBASE_code)
  r_idx = r %>% row_number()
  
  M = matrix(m_vector, nrow = n_countries, ncol = n_sect, byrow = TRUE)
  M_reordered = M[r_idx, , drop = FALSE]
  
  M_exio = rowsum(M_reordered, group = r$EXIOBASE_code)
  
  # Flatten back to original vector shape
  data_reordered <- as.vector(t(M_exio)) # 9800 mass vector to be used for intensity calcultion (J/kg)
  
  return(data_reordered)
}



# Make a block diagonal matrix having prod_map as diagonals, repeated for each FABIO region
block_diag_repeat <- function(M, n) {
  Msp <- as(M, "dgCMatrix") 
  kronecker(Diagonal(n), Msp)
}

p_fabio_exio = block_diag_repeat(as.matrix(prod_map), nrow(FABIO_reg)) 
# Normalize p_fabio_exio by (monetary) row sums to get a proper mapping matrix
a = p_fabio_exio %*% Diagonal(x=reorder_countries(EXIO_x)) 
a_inv = 1/Matrix::rowSums(a)
a_inv[!is.finite(a_inv)] = 0
a = Diagonal(x=a_inv) %*% a

# Multiply diag(FABIO_x) with p_fabio_exio to get a mass vector in EXIO classification
FABIO_exio_x = Diagonal(x=FABIO_x) %*% a
v_mass = colSums(FABIO_exio_x)
exio_mass = reorder_countries_mass(v_mass) # mass vec in EXIO sec/reg

# Intensity by mass (J/kg or hr/kg)
exio_en_int_mass = ene_net / exio_mass
exio_hr_m_int_mass = lab_male / exio_mass
exio_hr_f_int_mass = lab_female / exio_mass
# Replace NaN and Inf with 0
exio_en_int_mass[!is.finite(exio_en_int_mass)] = 0
exio_hr_m_int_mass[!is.finite(exio_hr_m_int_mass)] = 0
exio_hr_f_int_mass[!is.finite(exio_hr_f_int_mass)] = 0

# Re-order intensity vectors to match FABIO region order
exio_en_intm_ord = reorder_countries(exio_en_int_mass)
exio_hr_m_intm_ord = reorder_countries(exio_hr_m_int_mass)
exio_hr_f_intm_ord = reorder_countries(exio_hr_f_int_mass)

# Finally, calculate FABIO-based EXIO energy and labor use matrices
FABIO_en = FABIO_exio_x %*% Diagonal(x=exio_en_intm_ord)
FABIO_hr_m = FABIO_exio_x %*% Diagonal(x=exio_hr_m_intm_ord)
FABIO_hr_f = FABIO_exio_x %*% Diagonal(x=exio_hr_f_intm_ord)

