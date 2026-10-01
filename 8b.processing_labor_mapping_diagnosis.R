#### Diagnosis: EXIOBASE processing-sector labor landing on FABIO raw production ####
#
# convert_mass_vecs() / convert_intensities() (1.mrio_convert.R, 99.utils.R) split each
# FABIO product's PRODUCTION tonnage across the EXIOBASE sectors it is mapped to, then give
# each sector's full hours back to the FABIO products mapped to it. Where FABIO has no
# processing step for an EXIOBASE processing sector, that sector is mapped onto raw or
# "and products" items (FAO Food Balance Sheet, raw-commodity equivalents), so ALL of its
# labor -- including labor processing imported raw materials -- lands on the country's own
# farm or fishing output:
#   Fish products      <- Fish, Seafood (FABIO process: Fishing only)
#   Food products nec  <- ~40 crop, coffee/cocoa/tea, spice, egg, honey items
#   Dairy products     <- Milk - Excluding Butter, Butter, Ghee (dairy husbandry)
#   Processed rice     <- Rice
# Processing FABIO does model (oil extraction, sugar, beverages, slaughtering) maps onto
# products made where the processing happens, imported inputs included; those sectors are
# shown for comparison.
#
# For each EXIOBASE-modeled country this reports:
#   1. the share of its food-production hours that comes from each processing sector;
#   2. for the four suspect sectors: the sector's hours per tonne of the country's own
#      production of the items mapped to it (what the pipeline uses), against hours per tonne
#      of household use of those items from all origins (a processed-supply benchmark), and
#      the import share of that use. A large ratio with a high import share = inflation.
#
# Run in the session where 0.mrio_prep.R and 1.mrio_convert.R ran (needs dir_sat_FAB,
# EXIO_sect, prod_map, FABIO_x, FABIO_y_hh, io, regions, nrreg, nrcom, row_countries).

library(tidyverse)
library(Matrix)

PROJ <- "h:/MyDocuments/Projects/REMASS/food time"
setwd(PROJ)

need <- c("dir_sat_FAB", "EXIO_sect", "prod_map", "FABIO_x", "FABIO_y_hh", "regions", "nrreg", "nrcom")
miss <- need[!vapply(need, exists, logical(1))]
if (length(miss)) stop("Run 0.mrio_prep.R and 1.mrio_convert.R in this session first; missing: ",
                       paste(miss, collapse = ", "))
n_sec <- length(EXIO_sect)
hr_sat <- dir_sat_FAB[[2]] + dir_sat_FAB[[3]]        # FABIO rows x (EXIO sector x FABIO country), M.hr (hr_m + hr_f)
stopifnot(nrow(hr_sat) == nrreg * nrcom, ncol(hr_sat) == n_sec * nrreg)

suspect <- c("Fish products", "Food products nec", "Dairy products", "Processed rice")
modeled <- c("products of Vegetable oils and fats", "Sugar", "Beverages", "Products of meat cattle",
             "Products of meat pigs", "Products of meat poultry", "Meat products nec")
sec_of <- function(nm) {
  i <- which(tolower(trimws(EXIO_sect)) == tolower(nm))
  if (length(i) != 1) stop("EXIO sector not found uniquely: ", nm,
                           " | candidates: ", paste(grep(substr(nm, 1, 8), EXIO_sect, value = TRUE, ignore.case = TRUE), collapse = "; "))
  i
}
sec_idx <- vapply(c(suspect, modeled), sec_of, integer(1))
pm <- as.matrix(prod_map)                              # 123 FABIO items x 200 EXIO sectors (0/1)
items_of <- function(s) which(pm[, s] != 0)            # FABIO item positions mapped to EXIO sector s

cty <- regions$iso3c
rc  <- if (exists("row_countries")) row_countries else character(0)
targets <- which(!cty %in% rc)

# household use by item and consuming country, all origins (use_mat) and own-origin (dom_mat), tonnes
item_of_row <- rep(seq_len(nrcom), times = nrreg)      # FABIO rows: country-major, item within country
agg <- sparseMatrix(i = item_of_row, j = seq_along(item_of_row), x = 1, dims = c(nrcom, nrreg * nrcom))
use_mat <- as.matrix(agg %*% FABIO_y_hh)                # nrcom x nrreg
dom_mat <- vapply(seq_len(nrreg), function(j) as.numeric(FABIO_y_hh[((j - 1) * nrcom + 1):(j * nrcom), j]),
                  numeric(nrcom))

res <- map_dfr(targets, function(j) {
  rows <- ((j - 1) * nrcom + 1):(j * nrcom)            # FABIO production in country j
  cols <- ((j - 1) * n_sec + 1):(j * n_sec)            # country j's EXIO sector columns
  blk  <- hr_sat[rows, cols, drop = FALSE]
  total_hr <- sum(blk)
  map_dfr(names(sec_idx), function(nm) {
    s  <- sec_idx[[nm]]
    it <- items_of(s)
    prod_t  <- sum(FABIO_x[rows][it])                  # own production of the mapped items, tonnes
    use_all <- sum(use_mat[it, j])                     # household use of those items, all origins
    use_imp <- use_all - sum(dom_mat[it, j])           # ... of which imported
    sec_hr <- sum(blk[, s])
    tibble(country = cty[j], sector = nm, suspect = nm %in% suspect,
           sector_Mhr = sec_hr, share_of_food_hours = sec_hr / total_hr,
           prod_t = prod_t, use_t = use_all, import_share_of_use = ifelse(use_all > 0, use_imp / use_all, NA_real_),
           hr_per_t_on_production = ifelse(prod_t > 0, sec_hr * 1e6 / prod_t, NA_real_),
           hr_per_t_on_use        = ifelse(use_all > 0, sec_hr * 1e6 / use_all, NA_real_))
  })
}) %>% mutate(inflation = hr_per_t_on_production / hr_per_t_on_use)

cat("\n---- 1. Share of each country's food-production hours from processing sectors ----\n")
shares <- res %>% group_by(country) %>%
  summarise(suspect_4 = sum(share_of_food_hours[suspect]), modeled = sum(share_of_food_hours[!suspect]),
            fish_products = share_of_food_hours[sector == "Fish products"],
            food_nec = share_of_food_hours[sector == "Food products nec"],
            dairy = share_of_food_hours[sector == "Dairy products"], .groups = "drop") %>%
  arrange(desc(suspect_4))
print(shares %>% mutate(across(where(is.numeric), ~ round(.x, 2))), n = Inf)

cat("\n---- 2. Suspect sectors: hours per tonne on own production vs on household use (all origins) ----\n")
cat("(inflation = per-tonne on production / per-tonne on use; with a high import share it marks labor\n",
    " applied to imported raw materials being loaded onto domestic output)\n")
for (nm in suspect) {
  cat("\n", nm, ":\n", sep = "")
  print(res %>% filter(sector == nm, sector_Mhr > 0) %>%
          arrange(desc(inflation)) %>%
          transmute(country, `share of food hrs` = round(share_of_food_hours, 2),
                    `import share of use` = round(import_share_of_use, 2),
                    `hr/t on production` = round(hr_per_t_on_production),
                    `hr/t on use` = round(hr_per_t_on_use),
                    inflation = round(inflation, 1)) %>%
          slice_head(n = 12),
        n = Inf, width = Inf)
}
cat("\nCorrelation across countries, suspect sectors: log(inflation) vs import share of use:\n")
print(res %>% filter(suspect, is.finite(inflation), inflation > 0, !is.na(import_share_of_use)) %>%
        group_by(sector) %>%
        summarise(n = n(), r = cor(log(inflation), import_share_of_use),
                  median_inflation = median(inflation), .groups = "drop") %>%
        mutate(across(where(is.double), ~ round(.x, 2))))

write_csv(res, "results/processing_labor_mapping_diagnosis.csv")
cat("\nWrote results/processing_labor_mapping_diagnosis.csv\n")
