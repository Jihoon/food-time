library(tidyverse)
library(Matrix)

# Population
library(gt)
data(countrypops)

# countrypops (from the countrypops package) is missing Taiwan (TWN) --
# common for population datasets that follow UN-membership conventions.
# Supplement it from a hand-added CSV (region given as 2-letter "TW") so
# pop-normalized figures don't silently drop Taiwan -- previously this sent
# it down the pop=NA / dropped-from-plots path (see the
# setdiff(regions$iso3c, unique(countrypops$country_code_3)) check in
# 2.analyze_result.R).
taiwan_pop <- read.csv("data/Taiwan-population.csv") %>%
  transmute(country_name = "Taiwan", country_code_2 = "TW", country_code_3 = "TWN",
            year = Year, population = Population)
countrypops <- bind_rows(countrypops, taiwan_pop)

source("99.utils.R")

# Key parameters for config
year = yr = 2020  # Year when both labor and energy satellites are available
type = "pxp"

# Set up paths for data loading
EXIO_path = paste0("H:/MyDocuments/Data/EXIOBASE3/IOT_", year, "_", type)
EXIO_path_e = paste0(EXIO_path, "/energy/F.txt")
EXIO_path_l = paste0(EXIO_path, "/employment/F.txt")

FABIO_path = "H:/MyDocuments/Data/FABIO/input/"


#### 1. FABIO import ####

# read "product_concordance" sheet from "fabio-exiobase-v2.xlsx" under 'data' folder
# dim = 123x200
prod_map = readxl::read_xlsx("data/fabio-exiobase-v2.xlsx", 
                             sheet = "product_concordance", 
                             range = "E4:GV126",
                             col_names = FALSE) %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) # Replace NAs with 0s

FABIO_x = readRDS(file.path(FABIO_path,"X.rds"))[,as.character(year)]
FABIO_y = readRDS(file.path(FABIO_path,"Y.rds"))[[as.character(year)]]
FABIO_L = readRDS(file.path(FABIO_path, paste0(year, "_L_value.rds")))
FABIO_y_hh = FABIO_y[,grep("food", colnames(FABIO_y))] # 23001x187
FABIO_x_hh = FABIO_L %*% FABIO_y_hh

# Make FABIO-EXIOBASE region mapping 
reg_map = readxl::read_xlsx("data/fabio-exiobase.xlsx", 
                            sheet = "regions_concordance", 
                            col_names = TRUE) 
FABIO_reg = readxl::read_xlsx(paste0(FABIO_path, "fabio_classifications_v2.xlsx"), 
                              sheet = "Countries") %>% select(-area) %>%
  rename(ISO = `iso3c`, FAO_code = `area_code`) %>%
  left_join(reg_map) %>%
  # replace NA cells with values from Country=="RoW Europe" (Lichtenstein, Monaco, Andorra, San Marino etc.)
  mutate(EXIOBASE_code =
           as.numeric(ifelse(EXIOBASE_code=="NA", 47, EXIOBASE_code)),
         EXIOBASE =
           ifelse(EXIOBASE=="NA", "RoW Europe", EXIOBASE))
# Note: RoW countries' mean GDP/cap is close to that of Italy (<- Perplexity)

# Data-quality fix: code 41 (TWN's own individually-modeled EXIO region) is
# mislabeled "RoW Asia and Pacific" in the source concordance -- collides
# with the genuine RoW Asia & Pacific aggregate (code 45, ~37 countries),
# breaking anything that treats region names as unique (e.g. factor levels)
# and wrongly flagging Taiwan as a RoW aggregate member (is_row/row_countries
# match on grepl("RoW", EXIOBASE)). Taiwan is one of EXIOBASE's individually-
# modeled economies; correct its label.
stopifnot(identical(sort(unique(FABIO_reg$EXIOBASE[FABIO_reg$EXIOBASE_code == 41])), "RoW Asia and Pacific"))
FABIO_reg$EXIOBASE[FABIO_reg$EXIOBASE_code == 41] <- "Taiwan"

# Load FABIO metadata
library(data.table)
library(tidyselect)

regions <- fread(file=file.path(FABIO_path,"regions.csv"))
items <- fread(file=file.path(FABIO_path,"items.csv"))
nrreg <- nrow(regions)
nrcom <- nrow(items)
io <- fread(file.path(FABIO_path,"io_labels.csv"))
fd <- fread(file=file.path(FABIO_path,"fd_labels.csv"))



#### 2. EXIO import ####

n_col <- ifelse(type == "ixi", 7989, 9802)

# Read EXIO energy and labor direct satellite accounts
# Energy in TJ
# Labor in M.hour
ene = as.matrix(data.table::fread(EXIO_path_e, select = 2:(n_col-1), skip = 3,
                                  header = F))
lab = as.matrix(data.table::fread(EXIO_path_l, select = 2:(n_col-1), skip = 3,
                                  header = F))

# Read EXIO main matrices
EXIO_x = as.matrix(data.table::fread(paste0(EXIO_path, "/x.txt"), select = 3, 
                                     header = T))
EXIO_Y = as.matrix(data.table::fread(paste0(EXIO_path, "/Y.txt"), skip = 3, 
                                     select = 3:(7*49+2), header = F))
EXIO_reg = data.table::fread(paste0(EXIO_path, "/unit.txt"), header = T)
n_reg_EXIO = length(unique(EXIO_reg$region)) # = 49
EXIO_sect = unique(EXIO_reg$sector)

# Take the 'net energy' or 'final energy'?
# Note: "A net energy use account includes the energy content of all energy products 
# before being combusted or used for non-energy purposes, as well as the pseudo renewable energy 
# (before it is transformed into electricity or heat). So, it is not the output of 
# the electricity and heat generating industries, but only the input to it. 
# In other words, it is the final energy use, but with the inputs into heat and 
# electricity generation instead of the outputs." - Rasul et al. (2024) 
# Currently taking 'final energy' (ene[2,]) 
ene_net = matrix(ene[2,])

# Take the total male/female labor hours (direct)
lab_male = matrix(colSums(lab[c(7,9,11),]))
lab_female = matrix(colSums(lab[c(8,10,12),]))


#### 2.1. Population import ####
# Population data for per capita calculations
pop_y = subset(countrypops, year == yr) %>%
  select(iso3c = country_code_3, pop = population)

pop_y = pop_y %>% right_join(regions, by = "iso3c") %>%
  select(iso3c, pop)


#### 3. Prepare/Read indirect satellite through Leontief inverse 'L' ####
# DO ONLY ONCE AND SAVE THE RESULTING L MATRIX
# 
# # read matrices
# EXIO_A = as.matrix(data.table::fread(paste0(EXIO_path, "/A.txt"), select = 3:n_col, skip = 3,
#                                      header = F)) %>%
#   as("dgCMatrix")
# I <- diag(ncol(EXIO_A))
# EXIO_L <- Matrix::solve(I - EXIO_A) # inv_A <- chol2inv(chol(A)) faster?
# 
# # Save Leontief inverse L
# saveRDS(EXIO_L, file = paste0("data/EXIO_L_", year, "_", type, ".rds"))

# Load Leontief inverse L
EXIO_L = readRDS(file = paste0("data/EXIO_L_", year, "_", type, ".rds"))

# Note: This is a prep step to derive indirect energy and labor use satellites for non-food sectors only.       

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



#### 4. Read Fajzel data ####

path_GHD = "H:/MyDocuments/Data/GlobalHumanDay/"

# 0 Identify available country/year combinations in the data ####
# Read file names from the inputData folder under path_GHD
file_names = list.files(path = paste0(path_GHD, "inputData/"), pattern = "*.csv", full.names = TRUE)

# Use regex to extract the country names from three letters coming after "inputData/" and before "_TUS"
countries = str_extract(file_names, "(?<=inputData/)[A-Z]{3}(?=_TUS)")

# Use regex to extract years from the integers in file_names containing "TUS"
years = str_extract(file_names, "(?<=TUS_)[0-9]{4}")

# 4.1. Read the raw country data ####
# df_ghd = read_csv(paste0(path_GHD, "outputData/country_specific_results.csv"))
df_ghd_gender = read_csv(paste0(path_GHD, "outputData/gender_split/TUS_only_M24_251017.csv")) %>%
  filter(subcategory %in% c("processing", "preparation", "growth_collection")) %>%
  mutate(subcategory = case_when(
    subcategory == "processing" ~ "processing_non.econ",
    subcategory == "preparation" ~ "preparation_non.econ",
    subcategory == "growth_collection" ~ "growth_collection_non.econ"
  ))

# Non-economic food time by gender
df_ghd_gender = df_ghd_gender %>% rename(footprint_type = subcategory) %>%
  # Pivot longer to have a column "type" with values "hr_m" and "hr_f", and a column "per_capita_value" with the corresponding hours
  pivot_longer(cols = c("maleTotalHours", "femaleTotalHours"), names_to = "type", values_to = "per_capita_value") %>%
  mutate(type = ifelse(type == "maleTotalHours", "hr_m", "hr_f")) %>% 
  select(country, type, footprint_type, per_capita_value) %>%
  drop_na() # Drop rows with NA values in per_capita_value

# 4.2. Read the restaurant data ####
# Source has empl_restaurants counts but a blank mean_hrs_restaurants for a few
# countries (JPN, KOR, MNG, TUN as of this data vintage) -- drop those rather than
# carrying NA per_capita_value through to the plots (df_water_energy below does the
# same for its own source).
df_restaurant = read_csv(paste0(path_GHD, "outputData/restaurants_accommodations_time.csv")) %>%
  select(country, hr_m = mean_hrs_restaurants, hr_f = mean_hrs_restaurants) %>%
  # stack df_restaurant twice with type "hr_m" and "hr_f" to match the structure of df_ghd_gender
  pivot_longer(cols = c(hr_m, hr_f), names_to = "type", values_to = "per_capita_value") %>%
  mutate(footprint_type = "preparation_econ") %>%
  select(country, type, footprint_type, per_capita_value) %>%
  drop_na()

# 4.2b. Read water and firewood collection time (no gender split; apply the same value to both genders) ####
df_water_energy = read_csv("data/water_firewood.csv") %>%
  rename(water_non.econ = water_supply, energy_non.econ = forestry) %>%
  pivot_longer(cols = c(water_non.econ, energy_non.econ), names_to = "footprint_type", values_to = "per_capita_value") %>%
  mutate(hr_m = per_capita_value, hr_f = per_capita_value) %>%
  select(country, footprint_type, hr_m, hr_f) %>%
  pivot_longer(cols = c(hr_m, hr_f), names_to = "type", values_to = "per_capita_value") %>%
  drop_na()

# 4.2c. Read IEA residential energy end-use data — household energy for food-related appliances
# (cooking, refrigeration/freezing, dish washing), GJ per capita, 2020.
# CORRECTION has duplicate TC/X_TC rows per end-use (temperature-corrected vs not); these
# appliance end-uses aren't weather-sensitive so the two are identical — keep only "TC" to
# avoid double-counting. ####
df_iea_residential = read_csv("data/IEA - EEI RESIDENTIAL.csv") %>%
  mutate(country = countrycode::countrycode(`Country/Region`, origin = "country.name", destination = "iso3c")) %>%
  filter(!is.na(country), UNIT == "GJ_CAP", TIME_PERIOD == 2020, CORRECTION == "TC")

df_household_energy_food = df_iea_residential %>%
  filter(ENDUSE %in% c("FREEZER", "REFFREEZ", "REFRIG", "DISH_WASH", "COOKING")) %>%
  group_by(country) %>%
  summarise(energy_gj_cap = sum(OBS_VALUE, na.rm = TRUE), .groups = "drop")

# 4.2d. Availability check: COOKING vs TOTAL residential energy, for countries reporting both ####
df_cooking_total = df_iea_residential %>%
  filter(ENDUSE %in% c("COOKING", "TOTAL")) %>%
  select(country, ENDUSE, OBS_VALUE) %>%
  pivot_wider(names_from = ENDUSE, values_from = OBS_VALUE) %>%
  filter(!is.na(COOKING), !is.na(TOTAL)) %>%
  mutate(cooking_total_ratio = COOKING / TOTAL)

# 4.2e. By-country breakdown: COOKING, DISH_WASH, refrigeration/freezing (FREEZER+REFFREEZ+REFRIG), TOTAL ####
df_iea_breakdown = df_iea_residential %>%
  filter(ENDUSE %in% c("COOKING", "DISH_WASH", "FREEZER", "REFFREEZ", "REFRIG", "TOTAL")) %>%
  mutate(ENDUSE = ifelse(ENDUSE %in% c("FREEZER", "REFFREEZ", "REFRIG"), "REFRIG_ALL", ENDUSE)) %>%
  group_by(country, ENDUSE) %>%
  summarise(OBS_VALUE = sum(OBS_VALUE, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = ENDUSE, values_from = OBS_VALUE) %>%
  mutate(FOOD_ENERGY_TOTAL = COOKING + DISH_WASH + REFRIG_ALL) %>%
  select(country, COOKING, DISH_WASH, REFRIG_ALL, FOOD_ENERGY_TOTAL, TOTAL) 

# 4.3. Combine the dataframes ####
df_ghd_combined = bind_rows(df_ghd_gender, df_restaurant, df_water_energy)

cty_ghd = unique((df_ghd_gender %>% filter(footprint_type=="preparation_non.econ"))$country)
