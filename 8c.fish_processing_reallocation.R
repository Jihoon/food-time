#### Option A: attach fish-processing labor to where fish is eaten, not where it is caught ####
#
# 8b showed that EXIOBASE "Fish products" (processing) labor lands on each country's own
# fish catch, because FABIO models fishing but no fish processing. Countries that process
# mostly imported fish get inflated hours per tonne of catch (CHE 31x, AUT 19x, BEL 11x,
# SVK 10x), and in Malta fish processing is half of all food-production hours.
#
# Option A (a bound, not a pipeline change): treat processing as happening in the country
# where the fish is used.
#   current:  processing hours of country c ride on c's catch, through FABIO's supply chain,
#             to whoever eats it:   cur_j = sum_c  p_c * x_hh[(c, fish), j]
#             with p_c = processing hours / catch (hr per tonne caught)
#   option A: processing rate per tonne PROCESSED in c, r_c = processing hours /
#             (catch + fish imported for any final use), applied to the fish households in c
#             eat from all origins:   new_j = r_j * sum_o y_hh[(o, fish), j]
# Exported processed fish is still charged to domestic consumers (A's known bias; option C,
# EXIOBASE's own trade in processed fish, would fix that). Imports used as feed (FABIO Z) are
# left out of the denominator.
#
# Reports, per EXIO-modeled country: paid food hours per capita (household food demand) before
# and after, and the change in CF_paid. Run in the session where 0, 1 and 2.analyze_result.R
# ran (needs dir_sat_FAB, EXIO_sect, FABIO_x, FABIO_y, FABIO_y_hh, FABIO_x_hh, io, regions,
# nrreg, nrcom, pop_data_yr, tradeoff_protein_econlabor_consump).

library(tidyverse)
library(Matrix)

PROJ <- "h:/MyDocuments/Projects/REMASS/food time"
setwd(PROJ)

need <- c("dir_sat_FAB", "EXIO_sect", "FABIO_x", "FABIO_y", "FABIO_y_hh", "FABIO_x_hh", "io", "regions",
          "nrreg", "nrcom", "pop_data_yr", "tradeoff_protein_econlabor_consump")
miss <- need[!vapply(need, exists, logical(1))]
if (length(miss)) stop("Run 0, 1 and 2.analyze_result.R in this session first; missing: ", paste(miss, collapse = ", "))

cty <- regions$iso3c
n_sec <- length(EXIO_sect)
f <- which(io$item[1:nrcom] == "Fish, Seafood"); stopifnot(length(f) == 1)
s_proc <- which(tolower(trimws(EXIO_sect)) == "fish products"); stopifnot(length(s_proc) == 1)
fish_rows <- (seq_len(nrreg) - 1) * nrcom + f                         # fish row of each producing country

# processing hours per producing country (M.hr): fish row x that country's "Fish products" column
hr_sat <- dir_sat_FAB[[2]] + dir_sat_FAB[[3]]
proc_Mhr <- vapply(seq_len(nrreg), function(c) hr_sat[fish_rows[c], (c - 1) * n_sec + s_proc], numeric(1))
catch_t  <- FABIO_x[fish_rows]
p_c <- ifelse(catch_t > 0, proc_Mhr * 1e6 / catch_t, 0)               # hr per tonne caught (current)

# FABIO final demand: 6 use categories per country, country-major (fd_labels.csv order)
fd <- read_csv("H:/MyDocuments/Data/FABIO/input/fd_labels.csv", show_col_types = FALSE)
stopifnot(ncol(FABIO_y) == nrow(fd), all(unique(fd$iso3c) == cty))
fd_country <- match(fd$iso3c, cty)
y_fish <- as.matrix(FABIO_y[fish_rows, ])                               # producing country x fd columns, tonnes
final_use_by <- t(rowsum(t(y_fish), fd_country))                        # producing country x consuming country
imports_final_t <- colSums(final_use_by) - diag(final_use_by)          # fish imported for any final use
r_c <- ifelse(catch_t + imports_final_t > 0, proc_Mhr * 1e6 / (catch_t + imports_final_t), 0)  # hr per tonne processed

hh_fish_t <- colSums(as.matrix(FABIO_y_hh[fish_rows, ]))               # household fish use per consuming country
cur_hr <- as.numeric(p_c %*% as.matrix(FABIO_x_hh[fish_rows, ]))       # current: via catch, through the supply chain
new_hr <- r_c * hh_fish_t                                               # option A

pop <- pop_data_yr %>% select(iso3c, population)
paid_now <- tradeoff_protein_econlabor_consump %>%
  group_by(country) %>%
  summarise(paid_hr_cap_day = sum(hr_per_cap_day), cf_paid = sum(hr_per_50g_protein), .groups = "drop")

out <- tibble(country = cty, catch_t, imports_final_t, hh_fish_t, proc_Mhr,
              hr_per_t_caught = p_c, hr_per_t_processed = r_c, cur_hr, new_hr) %>%
  left_join(pop, by = c("country" = "iso3c")) %>%
  inner_join(paid_now, by = "country") %>%
  filter(!country %in% (if (exists("row_countries")) row_countries else character(0)), population > 0) %>%
  mutate(delta_hr_cap_day = (new_hr - cur_hr) / population / 365,
         paid_after = paid_hr_cap_day + delta_hr_cap_day,
         pct_change = 100 * delta_hr_cap_day / paid_hr_cap_day,
         cf_paid_after = cf_paid * paid_after / paid_hr_cap_day)

cat(sprintf("\nGlobal processing hours embodied in household fish: current %.0f M.hr, option A %.0f M.hr\n",
            sum(cur_hr) / 1e6, sum(new_hr) / 1e6))
cat("(A counts the household share of processing on all fish handled; the current method counts\n",
    " the household share of processing loaded onto catch. They need not match.)\n")
cat("\n---- Paid food hours per capita (household food) before and after, largest changes ----\n")
print(out %>% arrange(desc(abs(pct_change))) %>%
        transmute(country, `catch t` = round(catch_t), `imported t` = round(imports_final_t),
                  `hr/t caught` = round(hr_per_t_caught), `hr/t processed` = round(hr_per_t_processed),
                  `paid hr/cap/day` = round(paid_hr_cap_day, 3), `after` = round(paid_after, 3),
                  `% change` = round(pct_change, 1),
                  `CF_paid` = round(cf_paid, 3), `CF_paid after` = round(cf_paid_after, 3)),
      n = 20, width = Inf)
cat(sprintf("\nAcross %d countries: median |%% change| %.1f%%, %d countries change by more than 5%%.\n",
            nrow(out), median(abs(out$pct_change)), sum(abs(out$pct_change) > 5)))
cat(sprintf("Rank correlation of CF_paid before vs after: %.3f\n",
            cor(out$cf_paid, out$cf_paid_after, method = "spearman")))

write_csv(out, "results/fish_processing_option_A.csv")
cat("Wrote results/fish_processing_option_A.csv\n")
