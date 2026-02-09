library(mrio)
library(tidyverse)
library(Matrix)

# load some utility functions
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

# Take the 'net energy'
# Note: "A net energy use account includes the energy content of all energy products 
# before being combusted or used for non-energy purposes, as well as the pseudo renewable energy 
# (before it is transformed into electricity or heat). So, it is not the output of 
# the electricity and heat generating industries, but only the input to it. 
# In other words, it is the final energy use, but with the inputs into heat and 
# electricity generation instead of the outputs." - Rasul et al. (2024) 
ene_net = matrix(ene[2,])

# Take the total male/female labor hours (direct)
lab_male = matrix(colSums(lab[c(7,9,11),]))
lab_female = matrix(colSums(lab[c(8,10,12),]))



#### Prepare indirect satellite through Leontief inverse 'L' ####
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
