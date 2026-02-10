library(mrio)
library(tidyverse)
library(Matrix)

source("99.utils.R")

#### 1. Conversions to FABIO mass vector (123) to EXIO classifications (200, still mass) ####

# Prepare a conversion mapping (qualitative) from FABIO to EXIO sector mapping (for 187 countries)
p_fabio_exio = block_diag_repeat(as.matrix(prod_map), nrow(FABIO_reg)) 
  
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

  return(list(exio_mass, FABIO_mass_in_EXIO))
}

l = convert_mass_vecs()
exio_mass_x = l[[1]]  # EXIO mass vector (len 9800)
FABIO_x_in_EXIO = l[[2]]  # EXIO mass vector (len 9800)
rm(l)

# TEST: Validate sum(FABIO_x)-sum(FABIO_x[rowSums(p_fabio_exio)==0]) == sum(exio_mass) == sum(v_mass)




#### 2. Convert EXIO energy/labor satellites (and intensities) to FABIO classification ####

# Make EXIO satellite accounts in the same formatting 
dir_sat_exio = list(
  sat_en = ene_net,
  sat_hr_m = lab_male,
  sat_hr_f = lab_female) 

# Get direct energy/labor (hr) satellites in FABIO sectors
# dim = (23001, 37400), which can be row/column-summed to get vectors.
dir_sat_FAB = convert_intensities(dir_sat_exio)

# FABIO satellite account totals (intensities per tonne of product)
FABIO_en_int_d = rowSums(dir_sat_FAB$sat_en_FAB) / FABIO_x 
FABIO_hr_m_int_d = rowSums(dir_sat_FAB$sat_hr_m_FAB) / FABIO_x 
FABIO_hr_f_int_d = rowSums(dir_sat_FAB$sat_hr_f_FAB) / FABIO_x 

# Replace NaN and Inf with 0
FABIO_en_int_d[!is.finite(FABIO_en_int_d)] = 0
FABIO_hr_m_int_d[!is.finite(FABIO_hr_m_int_d)] = 0
FABIO_hr_f_int_d[!is.finite(FABIO_hr_f_int_d)] = 0

# Save FABIO-based EXIO satellite accounts and intensities
l_int = list(
  en_int_d = FABIO_en_int_d,
  hr_m_int_d = FABIO_hr_m_int_d,
  hr_f_int_d = FABIO_hr_f_int_d
)
saveRDS(l_int, file = paste0("data/FABIO_exio_satellites_", year, ".rds"))



#### 3. Calculate consumption-based calorie/protein footprints ####

# Calculate consumption-based matrices
FABIO_x_hh = FABIO_L %*% FABIO_y_hh

# Import nutrient coefficients by FABIO products
coeff_cal <- fread(file.path(FABIO_path,"nutrient_coefficients.csv")) %>%
  select(kcal_per_kg) %>%
  # concatenate the vector 187 times 
  slice(rep(1:n(), times = nrow(FABIO_reg)))

coeff_pro <- fread(file.path(FABIO_path,"nutrient_coefficients_protein.csv")) %>%
  mutate(protein_per_kg = 10*protein_per_100g) %>% select(protein_per_kg) %>%
  # concatenate the vector 187 times 
  slice(rep(1:n(), times = nrow(FABIO_reg)))

# Consumption total
FABIO_y_hh_cal <- Matrix::Diagonal(x=coeff_cal$kcal_per_kg) %*% FABIO_y_hh * 1000 # kcal
FABIO_y_hh_pro <- Matrix::Diagonal(x=coeff_pro$protein_per_kg) %*% FABIO_y_hh * 1000 # g

# Production total
FABIO_x_hh_cal <- Matrix::Diagonal(x=coeff_cal$kcal_per_kg) %*% FABIO_x_hh * 1000 # kcal
FABIO_x_hh_pro <- Matrix::Diagonal(x=coeff_pro$protein_per_kg) %*% FABIO_x_hh * 1000 # g

# Smaller values for poorer ctys; larger values for richer (and smaller) ctys
df_stat = data.frame(
  country = regions$iso3c,
  total_cal_eff = colSums(FABIO_y_hh_cal)/colSums(FABIO_x_hh_cal),
  total_protein_eff = colSums(FABIO_y_hh_pro)/colSums(FABIO_x_hh_pro)
)



#### 4. Energy and Labor footprint ####

# Total production-based footprint per country
energy_fp = Matrix::Diagonal(x=FABIO_en_int_d) %*% FABIO_x_hh #TJ = EJ/10^6 (sum = 22.5 EJ)
# Note: "In the United States, food production uses 10.11 quadrillion Btu annually" = 10.7 EJ
hr_m_fp = Matrix::Diagonal(x=FABIO_hr_m_int_d) %*% FABIO_x_hh #M.hr
hr_f_fp = Matrix::Diagonal(x=FABIO_hr_f_int_d) %*% FABIO_x_hh #M.hr 

colnames(energy_fp) = colnames(hr_m_fp) = colnames(hr_f_fp) = regions$iso3c

# Country-wise consumption-based footprint 
library(gt)
consumption = "food"
country = "KOR"
extension = "hr_m_int_d" #"en_int_d", "hr_f_int_d"

for (country in regions$iso3c) {
  
  Y_country <- FABIO_y[, which(fd$iso3c == country)]
  # Y_country <- Yi[, fd$iso3c == country]
  colnames(Y_country) <- fd$fd[fd$iso3c == country]
  
  pop = subset(countrypops, country_code_3 == country & year == yr)$population
  if (country %in% setdiff(regions$iso3c, unique(countrypops$country_code_3))) {pop=NA}
  
  print(paste("Population =", pop))
  
  # if(spread_stocks){
  #   stock_ratio <- Y_country[, "stock_addition"] / (rowSums(Y_country) - Y_country[, "stock_addition"])
  #   stock_ratio[!is.finite(stock_ratio)] <- 0
  #   Y_country <- as.data.table(as.matrix(Y_country))
  #   Y_country[, `:=`(food = food * (1 + stock_ratio),
  #                    other = other * (1 + stock_ratio),
  #                    tourist = tourist * (1 + stock_ratio),
  #                    unspecified = unspecified * (1 + stock_ratio),
  #                    stock_addition = 0)]
  # }
  
  int <- l_int[[extension]]
  MP <- int * FABIO_L
  
  # Initialize empty matrix to store results (row: exporter, col: importer)
  mat <- matrix(0, nrow=nrreg, ncol = nrreg, byrow=TRUE)
  rownames(mat) <- colnames(mat) <- sort(regions$iso3c)
  
  # Calculate footprints (energy & labor hr)
  country_consump = as.vector(as.matrix(Y_country[, consumption]))
  FP <- t(t(MP) * country_consump) # <= dim (23001x23001)
  colnames(FP) <- rownames(FP) <- paste0(io$iso3c, "_", io$item)
  FP <- as(FP, "TsparseMatrix")
  
  # Calculate calorie & protein production/consumption 
  x_country = t(t(FABIO_L) * country_consump) # mass flows
  # FP_cal = sweep(x_country, 2, 1000*coeff_cal, "*") 
  FP_cal = t(t(x_country) * as.vector(coeff_cal$kcal_per_kg)) * 1000 # kcal flows
  FP_pro = t(t(x_country) * as.vector(coeff_pro$protein_per_kg)) * 1000 # g-protein flows
  colnames(FP_cal) <- rownames(FP_cal) <- colnames(FP_pro) <- rownames(FP_pro) <- paste0(io$iso3c, "_", io$item)
  FP_cal <- as(FP_cal, "TsparseMatrix")
  FP_pro <- as(FP_pro, "TsparseMatrix")
  
  cal_consum = sum(country_consump*coeff_cal$kcal_per_kg)/365/pop*1000 #kcal/cap/day
  pro_consum = sum(country_consump*coeff_pro$protein_per_kg)/365/pop*1000 #kcal/cap/day
  
  results <- data.table(origin=rownames(FP)[FP@i + 1], target=colnames(FP)[FP@j + 1], 
                        value = FP@x, # M.hr
                        value_pcap = FP@x / pop * 1e6, # hr/cap (for labor)
                        cal = FP_cal@x,
                        cal_pcap = FP_cal@x / pop, # kcal per capita
                        pro = FP_pro@x,
                        pro_pcap = FP_pro@x / pop  # g-protein per capita
                        ) 
  results[,`:=`(country_consumer = country,
                year = year,
                indicator = extension,
                country_origin = substr(origin,1,3),
                item_origin = substr(origin,5,100),
                country_target = substr(target,1,3),
                item_target = substr(target,5,100))]
  
  results[,`:=`(group_origin = items$comm_group[match(results$item_origin,items$item)],
                group_target = items$comm_group[match(results$item_target,items$item)],
                continent_origin = regions$continent[match(results$country_origin, regions$iso3c)])]
  
  results$continent_origin[results$country_origin==country] <- country
  
  # Remove "Live animals" from the nutrient flows
  # CHECK: Need to be in energy/labor footprints?
  results = results %>%
    filter(group_origin != "Live animals")
  
  # print(paste(country, "ratio of FD vs. production:", 
  #             round(sum(as.matrix(Y_country)[, consumption]) / sum(results$value),4)))
  
  data_tot <- results %>%
    group_by(item_target, country_origin) %>%
    filter(value != 0) %>%
    summarise(value = (sum(value))) %>%
    spread(country_origin, value, fill = 0) %>% # Al
    # Add a row with column sums
    ungroup() %>%
    bind_rows(summarise(., item_target = "Total", across(-item_target, sum)))
  
  data.table::fwrite(data_tot, file=paste0("output/FABIO_", country,"_", year, "_", extension, "_", consumption,".csv"), sep=",")
  
  data_imp_tot = tail(data_tot, 1) %>% rename(importer = item_target) %>% mutate(importer = country)
  
  # Fill mat with data_imp_tot where rownames(mat) match names(data_imp_tot), and colnames(mat) match data_imp_tot$importer
  mat[rownames(mat) %in% names(data_imp_tot), data_imp_tot$importer[1]] <- 
    as.numeric(data_imp_tot[1, names(data_imp_tot) %in% rownames(mat)])
}

# Total Mcal trade between countries
saveRDS(mat, file.path("data/calorie_trade_mat.rds"))
  
  

# TODO
# - Consider per-capita results
# - Consider which metric to present
# - Consider how to handle energy & hour vs. calorie & protein
  
  
  