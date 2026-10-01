#### Does the allocation of food-manufacturing labor change consumption-based paid hours? ####
#
# EXIOBASE processing sectors that FABIO does not model as separate steps are mapped onto
# raw or "and products" items (8b): Food products nec (~40 crop, coffee/cocoa/tea, spice,
# egg and honey items), Dairy products (milk, butter), Processed rice (rice). Their labor
# therefore rides on the country's own raw production and travels with those raw crops
# through FABIO's supply chain to whoever eats them. In rich countries Food products nec
# alone is 25-46% of food-production hours (8b).
#
# Two opposite assumptions bound the consumption-based result:
#   current   manufacturing labor follows the raw crops (pipeline as is):
#               cur_j = sum_c sum_i p_ci * x_hh[(c, i), j],  p_ci = hours / tonne produced
#   option A  manufacturing labor stays with the manufacturing country's own consumers:
#               new_j = r_j * sum_i sum_o y_hh[(o, i), j],
#               r_c = hours / (production + imports for final use) of the mapped items
# Reality lies in between (processed food is traded too). A drops processing done abroad on
# imported processed food and processing for export (see 8c: fish, Norway), so it is a bound,
# not a correction. If CF_paid is similar under both, results do not hinge on the allocation.
#
# Run in the session where 0, 1 and 2.analyze_result.R ran (same objects as 8c, plus prod_map).

library(tidyverse)
library(Matrix)

PROJ <- "h:/MyDocuments/Projects/REMASS/food time"
setwd(PROJ)

need <- c("dir_sat_FAB", "EXIO_sect", "prod_map", "FABIO_x", "FABIO_y", "FABIO_y_hh", "FABIO_x_hh", "io",
          "regions", "nrreg", "nrcom", "pop_data_yr", "tradeoff_protein_econlabor_consump")
miss <- need[!vapply(need, exists, logical(1))]
if (length(miss)) stop("Run 0, 1 and 2.analyze_result.R in this session first; missing: ", paste(miss, collapse = ", "))

cty <- regions$iso3c
n_sec <- length(EXIO_sect)
pm <- as.matrix(prod_map)
hr_sat <- dir_sat_FAB[[2]] + dir_sat_FAB[[3]]
fd <- read_csv("H:/MyDocuments/Data/FABIO/input/fd_labels.csv", show_col_types = FALSE)
stopifnot(ncol(FABIO_y) == nrow(fd), all(unique(fd$iso3c) == cty))
fd_country <- match(fd$iso3c, cty)

allocate <- function(sector) {
  s  <- which(tolower(trimws(EXIO_sect)) == tolower(sector)); stopifnot(length(s) == 1)
  it <- which(pm[, s] != 0)                                              # FABIO items mapped to the sector
  rows <- as.vector(outer(it, (seq_len(nrreg) - 1) * nrcom, "+"))        # (country, item) rows, all countries
  prod_c <- rep(seq_len(nrreg), each = length(it))                       # producing country of each row
  h <- vapply(seq_along(rows), function(k) hr_sat[rows[k], (prod_c[k] - 1) * n_sec + s], numeric(1)) * 1e6  # hours
  x <- FABIO_x[rows]
  p <- ifelse(x > 0, h / x, 0)                                           # hr per tonne produced
  cur <- as.numeric(p %*% as.matrix(FABIO_x_hh[rows, ]))                 # current, per consuming country
  H_c    <- tapply(h, prod_c, sum)                                       # sector hours by country
  prod_t <- tapply(x, prod_c, sum)
  fu <- t(rowsum(t(as.matrix(FABIO_y[rows, ])), fd_country))             # rows x consuming country, final use
  fu_by <- rowsum(fu, prod_c)                                            # producing x consuming country
  imp_t <- colSums(fu_by) - diag(fu_by)
  r <- ifelse(prod_t + imp_t > 0, H_c / (prod_t + imp_t), 0)             # hr per tonne processed
  hh_t <- colSums(as.matrix(FABIO_y_hh[rows, ]))                         # household use of the items, all origins
  tibble(country = cty, sector = sector, cur = cur, new = r * hh_t, own_hours = H_c)
}

sectors <- c("Food products nec", "Dairy products", "Processed rice")
alloc <- map_dfr(sectors, allocate)
cat("\nGlobal hours embodied in household food, current vs option A (M.hr):\n")
print(alloc %>% group_by(sector) %>% summarise(current = sum(cur) / 1e6, option_A = sum(new) / 1e6,
                                               .groups = "drop") %>% mutate(across(where(is.double), round)))

pop <- pop_data_yr %>% select(iso3c, population)
paid_now <- tradeoff_protein_econlabor_consump %>% group_by(country) %>%
  summarise(paid = sum(hr_per_cap_day), cf_paid = sum(hr_per_50g_protein), .groups = "drop")
rc <- if (exists("row_countries")) row_countries else character(0)

delta <- alloc %>% group_by(country) %>%
  summarise(d_food_nec = sum((new - cur)[sector == "Food products nec"]),
            d_all3 = sum(new - cur), .groups = "drop") %>%
  left_join(pop, by = c("country" = "iso3c")) %>% inner_join(paid_now, by = "country") %>%
  filter(!country %in% rc, population > 0) %>%
  mutate(paid_A_food_nec = paid + d_food_nec / population / 365,
         paid_A_all3     = paid + d_all3 / population / 365,
         pct_food_nec = 100 * (paid_A_food_nec / paid - 1),
         pct_all3     = 100 * (paid_A_all3 / paid - 1),
         cf_paid_A_all3 = cf_paid * paid_A_all3 / paid)

cat("\n---- Paid food hours per capita: current vs option A, largest changes ----\n")
print(delta %>% arrange(desc(abs(pct_all3))) %>%
        transmute(country, `paid hr/cap/day` = round(paid, 3),
                  `A: food nec` = round(paid_A_food_nec, 3), `% food nec` = round(pct_food_nec, 1),
                  `A: all three` = round(paid_A_all3, 3), `% all three` = round(pct_all3, 1),
                  `CF_paid` = round(cf_paid, 3), `CF_paid A` = round(cf_paid_A_all3, 3)),
      n = 20, width = Inf)
cat(sprintf("\nAcross %d countries (all three sectors): median |%% change| %.1f%%, %d countries beyond 10%%, %d beyond 25%%.\n",
            nrow(delta), median(abs(delta$pct_all3)), sum(abs(delta$pct_all3) > 10), sum(abs(delta$pct_all3) > 25)))
cat(sprintf("Rank correlation of CF_paid, current vs option A: %.3f; of paid hr/cap/day: %.3f\n",
            cor(delta$cf_paid, delta$cf_paid_A_all3, method = "spearman"),
            cor(delta$paid, delta$paid_A_all3, method = "spearman")))
cat("\nRead: small changes and rank correlations near 1 -> results do not hinge on where manufacturing",
    "\n      labor is attached. Large changes for some countries -> those countries' paid hours depend on",
    "\n      an allocation the data cannot settle; report the range or build EXIOBASE-based processing (option C).\n")
write_csv(delta, "results/food_manufacturing_allocation_bounds.csv")
cat("Wrote results/food_manufacturing_allocation_bounds.csv\n")
