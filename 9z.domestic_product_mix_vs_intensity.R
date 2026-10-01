#### Domestic food systems: product mix or per-tonne intensity? ####
#
# 9y (last block) showed that the import/domestic intensity ratio is ~85-95% driven
# by the DOMESTIC side: per tonne of food, domestic paid hours vary ~13-fold across
# countries and energy ~16-fold, while imported food looks alike everywhere. This
# script asks where that domestic spread comes from, for the domestically produced
# and consumed stream (Y_dom, 2.analyze_result.R section 1.3), domestic effort only,
# food sectors' direct intensities:
#
#   hours per tonne of final domestic food
#     = throughput (tonnes produced per tonne of final food; e.g. feed behind meat)
#     x mix        (which products are produced, at reference per-tonne intensities)
#     x intensity  (how many hours each product takes per tonne in this country,
#                   at a reference mix)
#   plus an interaction term, since mix and intensity are indexed separately.
#
# Reference intensity per product: geometric mean across countries producing it.
# Reference mix: mean production share across countries. Countries: EXIO-modeled
# (non-RoW) with population. Same for energy.
#
# Run in the pipeline session (after 0, 1, 2.analyze_result.R), which holds
# FABIO_L, FABIO_y_hh / Y_dom (or X_dom), l_int_d, io, items, regions, coeff_cal,
# nrreg, nrcom, row_countries.

library(tidyverse)
library(Matrix)

PROJ <- "h:/MyDocuments/Projects/REMASS/food time"
setwd(PROJ)

need <- c("l_int_d", "io", "items", "regions", "coeff_cal", "nrreg", "nrcom", "row_countries")
miss <- need[!vapply(need, exists, logical(1))]
if (length(miss)) stop("Run 0, 1 and 2.analyze_result.R in this session first; missing: ", paste(miss, collapse = ", "))
if (!exists("X_dom")) {
  if (!exists("Y_dom")) stop("Y_dom / X_dom missing: run 2.analyze_result.R (section 1.3) in this session.")
  message("Computing X_dom = FABIO_L %*% Y_dom ...")
  X_dom <- FABIO_L %*% Y_dom
}
stopifnot(nrow(X_dom) == nrreg * nrcom, length(l_int_d$hr_m) == nrreg * nrcom)

cty   <- regions$iso3c
edible <- as.numeric(coeff_cal$kcal_per_kg > 0)              # same mask as mass_origin_df
prod_name <- items$item[match(io$item, items$item)][1:nrcom]  # product names, one country block
# intensities per tonne of each (country, product) row: hr/t (M.hr -> hr) and MJ/t (TJ -> MJ)
int_hr <- (l_int_d$hr_m + l_int_d$hr_f) * 1e6
int_en <- l_int_d$en * 1e6

targets <- which(!cty %in% row_countries)
long <- map_dfr(targets, function(j) {
  rows <- ((j - 1) * nrcom + 1):(j * nrcom)                   # production in the consuming country itself
  final_t <- sum(Y_dom[rows, j] * edible[rows])               # final domestic food, tonnes (edible)
  tibble(country = cty[j], product = prod_name, final_t = final_t,
         t = as.numeric(X_dom[rows, j]),                      # tonnes of each product produced for this stream
         i_hr = int_hr[rows], i_en = int_en[rows])
}) %>% filter(t > 0, final_t > 0)
# FABIO records live animals in heads, not tonnes. Their hours and energy (animal
# rearing) stay in the totals, but their head counts must not count as tonnes: they are
# left out of throughput, the production-share weights and the product-effects model.
live_animals <- items$item[items$comm_group == "Live animals"]
long <- long %>% mutate(is_live = product %in% live_animals)
cat(sprintf("%d live-animal products kept in hours/energy totals, excluded from tonnage and weights.\n",
            length(live_animals)))

# reference intensity per product: geometric mean over countries producing it with a positive intensity
ref <- long %>% filter(!is_live) %>% group_by(product) %>%
  summarise(ref_hr = exp(mean(log(i_hr[i_hr > 0]))), ref_en = exp(mean(log(i_en[i_en > 0]))),
            .groups = "drop")
long <- long %>% group_by(country) %>%
  mutate(w = ifelse(is_live, NA_real_, t / sum(t[!is_live]))) %>% ungroup() %>%
  left_join(ref, by = "product")
ref_mix <- long %>% filter(!is_live) %>% group_by(product) %>%
  summarise(wbar = sum(w) / n_distinct(long$country), .groups = "drop")

decomp <- long %>%
  left_join(ref_mix, by = "product") %>%
  group_by(country) %>%
  summarise(final_t = first(final_t), thru_t = sum(t[!is_live]),       # tonnes only, no head counts
            hr = sum(i_hr * t), en = sum(i_en * t),                      # all hours / energy, incl. rearing
            mix_hr = sum(w * ref_hr, na.rm = TRUE), mix_en = sum(w * ref_en, na.rm = TRUE),
            .groups = "drop")
# A first version built the intensity index as each country's own per-tonne intensities
# at a reference mix. Minor products with extreme per-tonne values then dominated it
# (India ~18,800 hr/t against 1-23 elsewhere) and the factors did not multiply back to
# the total. Instead, estimate a two-way model on per-tonne intensities, weighted by
# tonnes:  log i_cp = country_c + product_p. The country effect is how labor- (energy-)
# heavy a country's production is across the products it makes ("how"); the product
# effects, weighted by the country's own production shares, give the mix factor
# ("what"). Then  h_c (per throughput tonne) ~ exp(country_c) x sum_p w_cp exp(product_p),
# and the remainder (country-specific product deviations) is reported.
# Grass, fodder and other feed crops dominate tonnage, so the mix is also reported on
# final-food-relevant products only (items with kcal_per_kg > 0 in the FABIO coefficients).
edible_prod <- prod_name[edible[1:nrcom] > 0]

fe_factors <- function(k) {
  i_col <- paste0("i_", k)
  d <- long %>% filter(!is_live, .data[[i_col]] > 0)   # per-head intensities are not per tonne
  fit <- lm(reformulate(c("0 + country", "product"), paste0("log(", i_col, ")")), d, weights = t)
  cf <- coef(fit)
  ce <- tibble(country = sort(unique(d$country)),
               c_eff = unname(cf[paste0("country", sort(unique(d$country)))]))
  pe <- tibble(product = sort(unique(d$product))) %>%
    mutate(p_eff = coalesce(unname(cf[paste0("product", product)]), 0))   # baseline product = 0
  d %>% left_join(pe, by = "product") %>%
    group_by(country) %>%
    summarise(h = sum(.data[[i_col]] * t) / sum(t),                         # per throughput tonne
              mix     = sum(t * exp(p_eff)) / sum(t),
              mix_edible = sum((t * exp(p_eff))[product %in% edible_prod]) / sum(t[product %in% edible_prod]),
              .groups = "drop") %>%
    left_join(ce, by = "country") %>%
    mutate(what = k, level = exp(c_eff))
}
ff <- bind_rows(fe_factors("hr"), fe_factors("en"))

decomp <- decomp %>%
  mutate(hr_per_final_t = hr / final_t, en_per_final_t = en / final_t, throughput = thru_t / final_t)
fac <- ff %>% left_join(decomp %>% select(country, throughput, hr_per_final_t, en_per_final_t), by = "country") %>%
  mutate(per_final_t = ifelse(what == "hr", hr_per_final_t, en_per_final_t),
         # defined against the full per-final-tonne figure, so the four parts multiply back exactly
         remainder = per_final_t / (throughput * level * mix))

share_tab <- fac %>% group_by(what) %>% group_modify(function(g, key) {
  y <- log(g$per_final_t)
  parts <- list(throughput = log(g$throughput), level = log(g$level), mix = log(g$mix),
                remainder = log(g$remainder))
  # exact additive split of log(per final tonne) = throughput + level + mix + remainder;
  # each part's share of the variance = cov(part, y) / var(y), summing to 1
  sh <- vapply(parts, function(p) cov(p, y) / var(y), numeric(1))
  tibble(n = nrow(g), `max/min per final tonne` = round(exp(max(y) - min(y)), 1),
         `share: throughput` = sh[["throughput"]], `share: country level (how)` = sh[["level"]],
         `share: product mix (what)` = sh[["mix"]], `share: remainder` = sh[["remainder"]],
         `R2 mix (edible products only)` = summary(lm(y ~ log(g$mix_edible)))$r.squared)
}) %>% ungroup() %>% mutate(across(where(is.double), ~ round(.x, 2)),
                            what = ifelse(what == "hr", "hours", "energy"))
cat("\n---- Domestic stream: what carries the cross-country spread in hours / energy per tonne of final food ----\n")
cat("(shares of the variance of log per-final-tonne intensity; they sum to 1 across the four parts)\n")
print(share_tab, n = Inf, width = Inf)

# sanity: this direct, food-sector, domestic-effort measure against 9y's per-tonne domestic intensity
if (file.exists("results/import_premium_per_gram_vs_per_tonne.csv")) {
  p9y <- read_csv("results/import_premium_per_gram_vs_per_tonne.csv", show_col_types = FALSE)
  chk <- decomp %>% inner_join(p9y %>% filter(what == "paid") %>% select(country, per_t_dom), by = "country")
  cat(sprintf("\nCorrelation (log) with 9y's domestic paid hours per tonne (all sectors and effort origins): r = %.2f, n = %d\n",
              cor(log(chk$hr_per_final_t), log(chk$per_t_dom)), nrow(chk)))
}

# which products carry each country's hours, for the extremes
extremes <- decomp %>% arrange(hr_per_final_t) %>% slice(c(1:4, (n() - 3):n())) %>% pull(country)
cat("\n---- Top products by share of domestic hours, extreme countries (hours per tonne vs reference) ----\n")
print(long %>% filter(country %in% extremes) %>%
        group_by(country) %>% mutate(hr_share = i_hr * t / sum(i_hr * t)) %>%
        slice_max(hr_share, n = 4) %>% ungroup() %>%
        transmute(country, product, `share of hours` = round(hr_share, 2), `hr/t` = round(i_hr),
                  `ref hr/t` = round(ref_hr), `share of tonnes` = round(w, 2)),
      n = Inf, width = Inf)

cat("\nRead: 'country level (how)' large -> the same products take very different hours / energy per",
    "\n      tonne across countries (technology, farm structure); 'product mix (what)' large -> countries",
    "\n      differ in which products they produce for themselves; 'throughput' -> upstream tonnes (feed,",
    "\n      grass) per tonne of final food; 'remainder' -> country-specific deviations for particular",
    "\n      products. A negative share means that part offsets the others.\n")
write_csv(fac, "results/domestic_mix_vs_intensity.csv")
cat("Wrote results/domestic_mix_vs_intensity.csv\n")
