#### 4. Export direct food-sector labor hours by country ####
#
# Extracts the idx_food entries of lab_male/lab_female (direct EXIO labor
# satellites, M.hr) for each EXIOBASE region and writes a long-format CSV.
#
# Run this script in the same R session after script 0, which defines
# lab_male, lab_female, idx_food, and EXIO_reg.

library(dplyr)
library(data.table)

idx_food_full = rep(idx_food, times = nrow(EXIO_reg) / length(idx_food))

food_labor = EXIO_reg %>%
  mutate(lab_male = as.vector(lab_male),
         lab_female = as.vector(lab_female)) %>%
  filter(idx_food_full) %>%
  select(region, sector, lab_male, lab_female)

fwrite(food_labor, file = "output/food_sector_labor_by_country.csv")
