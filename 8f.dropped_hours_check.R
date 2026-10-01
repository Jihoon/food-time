#### Are any EXIOBASE food-sector hours / energy lost in the mapping to FABIO? ####
#
# convert_intensities() (99.utils.R) computes each EXIOBASE sector's per-tonne intensity as
# sector satellite / exio_mass_x and sets it to 0 when the sector received no FABIO tonnage
# (exio_mass_x == 0). That sector's hours and energy then vanish from every footprint. This
# compares, per EXIOBASE region x food sector, the satellite totals with what survives in
# dir_sat_FAB. It runs the TEST lines left commented in 1.mrio_convert.R (lines 92-93) and
# extends them to energy and to a region x sector breakdown. 98.validations.R checks mass only.
#
# Run in the session where 0.mrio_prep.R and 1.mrio_convert.R ran (needs lab_male, lab_female,
# ene_net, idx_food, EXIO_sect, exio_mass_x, dir_sat_FAB, reorder_countries_to_EXIO, EXIO_reg).

library(tidyverse)
library(Matrix)

need <- c("lab_male", "lab_female", "ene_net", "idx_food", "EXIO_sect", "exio_mass_x", "dir_sat_FAB",
          "reorder_countries_to_EXIO")
miss <- need[!vapply(need, exists, logical(1))]
if (length(miss)) stop("Run 0.mrio_prep.R and 1.mrio_convert.R in this session first; missing: ", paste(miss, collapse = ", "))

n_sec <- length(EXIO_sect)
n_reg <- length(as.numeric(lab_male)) / n_sec; stopifnot(n_reg == round(n_reg))
reg_names <- if (exists("EXIO_reg")) unique(EXIO_reg$region) else paste0("R", seq_len(n_reg))
food <- rep(as.logical(idx_food), times = n_reg)                       # EXIO order: sector within region

exio <- tibble(region = rep(reg_names, each = n_sec), sector = rep(EXIO_sect, times = n_reg),
               food = food, mass = as.numeric(exio_mass_x),
               hr_exio = as.numeric(lab_male) + as.numeric(lab_female), en_exio = as.numeric(ene_net))
# what survives: dir_sat_FAB columns (EXIO sector x FABIO country) summed and aggregated to EXIO regions
mapped_hr <- reorder_countries_to_EXIO(colSums(dir_sat_FAB[[2]]) + colSums(dir_sat_FAB[[3]]))
mapped_en <- reorder_countries_to_EXIO(colSums(dir_sat_FAB[[1]]))
stopifnot(length(mapped_hr) == nrow(exio))
exio <- exio %>% mutate(hr_mapped = as.numeric(mapped_hr), en_mapped = as.numeric(mapped_en))

f <- exio %>% filter(food)
tot <- f %>% summarise(across(c(hr_exio, hr_mapped, en_exio, en_mapped), sum))
cat(sprintf("\nFood sectors, all regions: hours mapped %.1f%% of EXIOBASE; energy mapped %.1f%%\n",
            100 * tot$hr_mapped / tot$hr_exio, 100 * tot$en_mapped / tot$en_exio))
lost <- f %>% mutate(hr_lost = hr_exio - hr_mapped, en_lost = en_exio - en_mapped)
cat(sprintf("Region-sector cells with satellite > 0 but no FABIO mass: %d (hours there: %.2f%% of food total)\n",
            sum(lost$mass == 0 & lost$hr_exio > 0),
            100 * sum(lost$hr_exio[lost$mass == 0]) / tot$hr_exio))

cat("\n---- Hours lost by sector (share of that sector's EXIOBASE hours) ----\n")
print(lost %>% group_by(sector) %>%
        summarise(hr_exio = sum(hr_exio), lost_pct = 100 * sum(hr_lost) / sum(hr_exio),
                  en_lost_pct = 100 * sum(en_lost) / sum(en_exio), .groups = "drop") %>%
        filter(abs(lost_pct) > 0.5 | abs(en_lost_pct) > 0.5) %>% arrange(desc(lost_pct)) %>%
        mutate(across(c(lost_pct, en_lost_pct), ~ round(.x, 1))),
      n = Inf, width = Inf)

cat("\n---- Regions losing the largest share of their food-sector hours ----\n")
print(lost %>% group_by(region) %>%
        summarise(lost_pct = 100 * sum(hr_lost) / sum(hr_exio), en_lost_pct = 100 * sum(en_lost) / sum(en_exio),
                  lost_cells = paste(sector[mass == 0 & hr_exio > 0], collapse = "; "), .groups = "drop") %>%
        arrange(desc(lost_pct)) %>% slice_head(n = 15) %>%
        mutate(across(c(lost_pct, en_lost_pct), ~ round(.x, 1))),
      n = Inf, width = Inf)
cat("\nRead: near 100% mapped and few empty cells -> nothing material is dropped. Losses concentrated",
    "\n      in a sector or region show where EXIOBASE labor has no FABIO production to attach to.\n")
dir.create("results", showWarnings = FALSE)
write_csv(lost, "results/dropped_hours_check.csv")
cat("Wrote results/dropped_hours_check.csv\n")
