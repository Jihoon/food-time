#### 2.3.export_time_split.R -- per-capita provisioning time split by country ####
#
# Writes results/provisioning_time_split_by_country.csv: hours per person per
# day (women + men, averaged over the whole population), consumption basis,
# for the directly-modelled countries with GHD household time (Fig. 1 sample).
# Source AFTER 2.analyze_result.R; computes no new footprints.
#
#   unpaid         *_non.econ (preparation, processing, growth/collection,
#                  energy, water)
#   paid_domestic  domestic food + non-food sector hours, plus preparation_econ
#   paid_imported  all import_cont_* components, food + non-food
#   total          unpaid + paid_domestic + paid_imported
#   export         hours residents supply to other countries' consumption;
#                  reported separately and left out of total, so totals sit
#                  below Fig. 1's positive bars

library(tidyverse)

if (!exists("b_direct")) stop("Source 2.analyze_result.R first -- missing object: b_direct")

time_split = b_direct %>%
  filter(type %in% c("hr_f", "hr_m")) %>%
  mutate(ft = as.character(footprint_type),
         component = case_when(
           grepl("_non\\.econ$", ft)                                   ~ "unpaid",
           ft %in% c("preparation_econ", "domestic_per_capita",
                     "domestic_per_capita_nf")                         ~ "paid_domestic",
           grepl("^import_cont_", ft)                                  ~ "paid_imported",
           ft %in% c("export_per_capita", "export_per_capita_nf")      ~ "export",
           TRUE                                                        ~ NA_character_))

unmapped = unique(time_split$ft[is.na(time_split$component)])
if (length(unmapped) > 0) stop("Unmapped footprint types: ", paste(unmapped, collapse = ", "))

# abs(): import components are stored positive here but guard against a sign flip upstream
time_split_by_country = time_split %>%
  group_by(country = as.character(country), component) %>%
  summarise(hr = sum(abs(per_capita_value), na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = component, values_from = hr, values_fill = 0) %>%
  mutate(total = unpaid + paid_domestic + paid_imported,
         share_unpaid        = unpaid / total,
         share_paid_domestic = paid_domestic / total,
         share_paid_imported = paid_imported / total) %>%
  select(country, unpaid, paid_domestic, paid_imported, total,
         share_unpaid, share_paid_domestic, share_paid_imported, export) %>%
  arrange(-total)

write.csv(time_split_by_country, "results/provisioning_time_split_by_country.csv", row.names = FALSE)
