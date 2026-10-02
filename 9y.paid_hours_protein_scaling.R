#### Does path b measure efficiency, or hours that do not scale with protein? ####
#
# CF = hours / protein, so log CF = log H - log P. Path b regresses log P on
# log GDP and log CF: protein sits on both sides. A negative b follows
# automatically if embodied hours per capita (H) grow less than one-for-one
# with protein supply (P) at a given income -- the "lower bars are wider"
# pattern in the protein mosaics, where rectangle area (hours) varies less
# than width (protein).
#
# Protein varies little across countries (~75-130 g/cap/day) while CF varies
# ~ninefold, so most CF variation is hours. The shared denominator still biases
# path b: if hours were unrelated to protein (net of GDP), the partial r of CF
# with protein would be -sd(P) / sqrt(sd(H)^2 + sd(P)^2), on GDP residuals.
# That null r is printed next to the observed partial r.
#
# Direct check: log H ~ log P + log GDP. The coefficient on log P is the
# protein elasticity of embodied hours, net of income.
# On GDP residuals, cov(CF, P) = (e - 1) var(P), so path b < 0 whenever e < 1.
#   e < 0     fewer hours per capita go with more protein: a link beyond the ratio
#   e ~ 0     hours do not move with protein; path b is the denominator arithmetic
#   0 < e < 1 hours rise with protein, but less than proportionally: the extra
#             protein costs fewer hours per gram than the average (diet mix,
#             sourcing), and the ratio turns that into a negative path b
#   e >= 1    no negative path b
# Same samples as Table 2 (9x), so each row pairs with a Table 2 row.
# H is rebuilt from the CF and protein columns: CF is per 50 g protein,
# protein is g/cap/day, so H = CF * P / 50 (hours or MJ per capita per day).
# CF_paid,dom is left out: its denominator is domestically produced protein,
# not total protein, so H cannot be rebuilt from log_protein.

library(tidyverse)

PROJ <- "h:/MyDocuments/Projects/REMASS/food time"
setwd(PROJ)

if (!exists("cf_paid_df38") || !exists("df2_collapsed")) source("9x.table2_mediation_summary.R")

X <- "log_gdp_pcap"; Y <- "log_protein"

scaling_row <- function(d, m, variant, sample) {
  d <- d[complete.cases(d[c(X, m, Y)]), ]
  d$log_H <- d[[m]] + d[[Y]] - log(50)
  fh <- lm(reformulate(c(Y, X), "log_H"), d)
  s <- summary(fh)$coefficients
  e <- s[Y, "Estimate"]; se <- s[Y, "Std. Error"]; dfr <- fh$df.residual
  fb <- lm(reformulate(c(X, m), Y), d)
  # spread of hours vs protein once income is taken out: the mosaic "area vs width" comparison
  rH <- resid(lm(reformulate(X, "log_H"), d)); rP <- resid(lm(reformulate(X, Y), d))
  rC <- resid(lm(reformulate(X, m), d))
  tibble(variant, sample, n = nrow(d),
         elasticity = e, ci_lo = e - qt(0.975, dfr) * se, ci_hi = e + qt(0.975, dfr) * se,
         p_vs_0 = s[Y, "Pr(>|t|)"],
         path_b = coef(fb)[[m]], partial_r_observed = cor(rC, rP),
         partial_r_null = -sd(rP) / sqrt(var(rH) + var(rP)),
         sd_ratio_H_to_P = sd(rH) / sd(rP))
}

scal <- bind_rows(
  scaling_row(df,             "log_cf",        "CF_time",   "EXIO-only"),
  scaling_row(df_collapsed,   "log_cf",        "CF_time",   "RoW-collapsed"),
  scaling_row(cf_paid_df33,   "log_cf_paid",   "CF_paid",   "EXIO-only"),
  scaling_row(cf_paid_df38,   "log_cf_paid",   "CF_paid",   "RoW-collapsed"),
  scaling_row(df2,            "log_cf_energy", "CF_energy", "EXIO-only"),
  scaling_row(df2_collapsed,  "log_cf_energy", "CF_energy", "RoW-collapsed")
)

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
cat("\n---- Protein elasticity of embodied hours (energy for CF_energy), net of GDP ----\n")
print(scal %>%
        transmute(variant, sample, n,
                  elasticity = sprintf("%.2f [%.2f, %.2f]", elasticity, ci_lo, ci_hi),
                  `p (=0)` = fmt_p(p_vs_0),
                  path_b = sprintf("%.2f", path_b),
                  `partial r` = sprintf("%.2f", partial_r_observed),
                  `null r` = sprintf("%.2f", partial_r_null),
                  `sd(H)/sd(P)` = sprintf("%.2f", sd_ratio_H_to_P)),
      n = Inf, width = Inf)
cat("\nRead: elasticity < 0 (CI below 0) -> fewer hours per capita go with more protein, net of income;",
    "\n      elasticity ~ 0 -> path b is mostly the shared protein denominator;",
    "\n      0 < elasticity < 1 -> extra protein costs fewer hours per gram than the average;",
    "\n      partial r vs null r -> how much of the observed correlation the ratio alone would produce.",
    "\n      sd(H)/sd(P) -> spread of per-capita hours relative to protein, both net of GDP.\n")

dir.create("results", showWarnings = FALSE)
write_csv(scal, "results/paid_hours_protein_scaling.csv")
cat("Wrote results/paid_hours_protein_scaling.csv\n")


#### Hours vs protein with and without the income control, by kind of hours ####
# Across the sample, higher-protein countries spend fewer hours (richer countries
# do both); net of GDP the relationship above turned positive for all-work hours.
# This splits all-work hours into paid and unpaid and shows both correlations
# side by side. CF_time and CF_paid share the denominator (total protein
# consumed), so unpaid = all-work - paid per capita.

hours_split <- function(d_all, d_paid, sample) {
  d <- inner_join(d_all %>% select(country, log_cf, all_of(c(X, Y))),
                  d_paid %>% select(country, log_cf_paid), by = "country") %>%
    filter(if_all(where(is.numeric), is.finite))
  stopifnot(nrow(d) > 10)
  P <- exp(d[[Y]])
  h <- tibble(all_work = exp(d$log_cf) * P / 50, paid = exp(d$log_cf_paid) * P / 50) %>%
    mutate(unpaid = all_work - paid)
  stopifnot(all(h$unpaid > 0))
  map_dfr(c("all_work", "unpaid", "paid"), function(k) {
    lh <- log(h[[k]])
    raw <- cor.test(lh, d[[Y]])
    pc  <- partial_cor(lh, d[[Y]], d[[X]])                # from 9.capability_set_income_control.R
    el  <- summary(lm(lh ~ d[[Y]] + d[[X]]))$coefficients[2, ]
    tibble(sample, n = nrow(d), hours = k,
           mean_hr_per_day = mean(h[[k]]),
           r_raw = unname(raw$estimate), p_raw = raw$p.value,
           r_net_gdp = pc$r, p_net_gdp = pc$p,
           elasticity_net_gdp = el[["Estimate"]], p_elasticity = el[["Pr(>|t|)"]])
  })
}

split <- bind_rows(hours_split(df, cf_paid_df33, "EXIO-only"),
                   hours_split(df_collapsed, cf_paid_df38, "RoW-collapsed"))
cat("\n---- log hours per capita vs log protein: raw, and net of GDP ----\n")
print(split %>%
        transmute(sample, n, hours, `hr/cap/day` = sprintf("%.2f", mean_hr_per_day),
                  `r raw` = sprintf("%.2f (%s)", r_raw, fmt_p(p_raw)),
                  `r net of GDP` = sprintf("%.2f (%s)", r_net_gdp, fmt_p(p_net_gdp)),
                  `elasticity net of GDP` = sprintf("%.2f (%s)", elasticity_net_gdp, fmt_p(p_elasticity))),
      n = Inf, width = Inf)
write_csv(split, "results/hours_protein_raw_vs_net.csv")
cat("Wrote results/hours_protein_raw_vs_net.csv\n")


#### Do hours follow total food (calories) rather than protein? ####
# Paid hours and energy are attached per tonne of each food product (Methods),
# so they could follow how much food is eaten rather than protein as such.
# Household (unpaid) hours come from time-use surveys, not tonnes. This block
# asks, net of GDP, whether each kind of hours follows food volume (calories),
# protein, or neither, and where protein differences at a given income come
# from: eating more food, or a more protein-rich diet.
#   kcal elasticity | GDP              ~ 1: hours rise in proportion to food volume;
#                                      ~ 0: hours do not follow food volume
#   protein elasticity | kcal, GDP     ~ 0: protein adds nothing once volume is known
#   decomposition of protein variation (net of GDP) into calories and protein
#   share of calories, with their covariance: a negative covariance means the
#   extra calories in high-calorie diets are protein-poor
# Calories and protein are highly collinear here, so read the "| kcal" and
# "| protein" columns together, not one alone. Calories are consumption-based
# (domestic + import) from the same FABIO flows as protein (2.analyze_result.R,
# kcal_consumption). EXIO-only sample (n = 33).

if (!exists("tradeoff_kcal_allwork_consump")) {
  stop("Calorie tables missing: run 2.analyze_result.R in this session, or re-run it and ",
       "2.2.save_objects_for_9.R so data/objects_for_9_2020.RData carries them.")
}
kcal_df <- tradeoff_kcal_allwork_consump %>%
  filter(!is_row) %>%
  group_by(country) %>%
  summarise(kcal_per_cap_day = first(kcal_per_cap_day), .groups = "drop")

d_k <- df %>% select(country, log_cf, all_of(c(X, Y))) %>%
  inner_join(cf_paid_df33 %>% select(country, log_cf_paid), by = "country") %>%
  inner_join(df2 %>% select(country, log_cf_energy), by = "country") %>%
  inner_join(kcal_df, by = "country") %>%
  filter(if_all(where(is.numeric), is.finite), kcal_per_cap_day > 0)
stopifnot(nrow(d_k) > 10)
P <- exp(d_k[[Y]])
d_k <- d_k %>%
  mutate(log_kcal   = log(kcal_per_cap_day),
         log_share  = .data[[Y]] - log_kcal,                            # protein per kcal
         log_all    = log(exp(log_cf) * P / 50),
         log_paid   = log(exp(log_cf_paid) * P / 50),
         log_unpaid = log(exp(log_cf) * P / 50 - exp(log_cf_paid) * P / 50),
         log_energy = log(exp(log_cf_energy) * P / 50))

coef_row <- function(fit, term) {
  s <- summary(fit)$coefficients
  sprintf("%.2f (%s)", s[term, "Estimate"], fmt_p(s[term, "Pr(>|t|)"]))
}
kcal_tab <- map_dfr(c(all_work = "log_all", unpaid = "log_unpaid", paid = "log_paid", energy = "log_energy"),
                    function(h) {
  f1 <- lm(reformulate(c("log_kcal", X), h), d_k)
  f2 <- lm(reformulate(c("log_kcal", Y, X), h), d_k)
  tibble(`kcal elasticity | GDP`                = coef_row(f1, "log_kcal"),
         `protein elasticity | kcal, GDP`       = coef_row(f2, Y),
         `kcal elasticity | protein, GDP`       = coef_row(f2, "log_kcal"),
         `R2 (kcal+GDP)` = sprintf("%.2f", summary(f1)$r.squared))
}, .id = "hours")
cat(sprintf("\n---- Hours per capita vs calories and protein, net of GDP (EXIO-only, n = %d) ----\n", nrow(d_k)))
print(kcal_tab, n = Inf, width = Inf)

# Where does protein variation at a given income come from: more food, or a protein-richer mix?
# log P = log kcal + log(P/kcal); split its GDP-residual variance accordingly.
rk <- resid(lm(reformulate(X, "log_kcal"), d_k))
rs <- resid(lm(reformulate(X, "log_share"), d_k))
vP <- var(rk + rs)
cat(sprintf(paste0("\nProtein variation net of GDP: calories %.0f%%, protein share of calories %.0f%%, ",
                   "covariance %.0f%% (sums to 100)\n",
                   "Correlation of calories with protein share, net of GDP: r = %.2f\n"),
            100 * var(rk) / vP, 100 * var(rs) / vP, 100 * 2 * cov(rk, rs) / vP, cor(rk, rs)))
cat("\nRead: kcal elasticity ~ 1 -> those hours rise in proportion to food volume at a given income;",
    "\n      kcal elasticity ~ 0 -> they do not follow food volume (product mix and sourcing matter more);",
    "\n      protein elasticity | kcal ~ 0 -> protein adds nothing beyond food volume.",
    "\n      Decomposition: a calorie share above 100% with a negative covariance means protein",
    "\n      differences come from eating more food whose extra calories are protein-poor;",
    "\n      a large protein-share term means they come from a more protein-rich diet.",
    "\n      Calories and protein are collinear here: read the '| kcal' and '| protein' columns together.\n")

write_csv(kcal_tab, "results/hours_kcal_protein.csv")
cat("Wrote results/hours_kcal_protein.csv\n")


#### Paid hours and energy by origin: which side fails to follow food supply? ####
# Split paid food-sector hours and embodied energy per capita by the origin of the
# final food product (domestically produced vs imported), the same split as protein
# (summary_pro_df_long, cat domestic/import) and food mass (mass_origin_df). This is
# the product-origin dimension, not where each hour was worked. EXIO-only sample (n = 33).
#   1. Does each origin's paid hours / energy follow that origin's protein, net of GDP?
#      (log H_o ~ log P_o + log GDP; elasticity ~ 1 = follows, ~ 0 = does not)
#   2. Which origin carries the cross-country spread in total paid hours / energy at a
#      given income? Levels are additive (H = H_dom + H_imp), so the GDP-residual
#      variance of H splits exactly into the two origins' contributions.
#   3. Does the import share of protein explain paid hours / energy per capita,
#      net of GDP and total protein?

# Split by the origin of the FINAL product (effort_consumption_df's protein_source),
# summing over where the effort itself happened, food sectors only (the scope of
# CF_paid). An earlier version took summary_food_df_long, whose domestic/import
# split is by where the work happened (2.analyze_result.R, section 1.3): totals
# matched, but hours were set against protein and mass split the other way.
origin_long <- effort_consumption_df %>%
  mutate(country = region_to_iso_1to1[as.character(exio_region)]) %>%
  filter(!is.na(country), sector == "food", type %in% c("hr_m", "hr_f", "en")) %>%
  mutate(origin = ifelse(protein_source == "domestic", "dom", "imp"),
         what   = ifelse(type == "en", "energy", "paid"),
         value  = ifelse(type == "en", per_capita_value / 365, per_capita_value)) %>%  # MJ or hr / cap / day
  group_by(country, what, origin) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

pro_origin <- summary_pro_df_long %>%
  filter(cat %in% c("domestic", "import")) %>%
  mutate(country = as.character(country), origin = ifelse(cat == "domestic", "dom", "imp")) %>%
  group_by(country, origin) %>%
  summarise(P = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

d_o <- origin_long %>%
  pivot_wider(names_from = c(what, origin), values_from = value) %>%
  inner_join(pro_origin %>% pivot_wider(names_from = origin, values_from = P, names_prefix = "P_"), by = "country") %>%
  inner_join(df %>% select(country, all_of(X)), by = "country") %>%
  mutate(paid = paid_dom + paid_imp, energy = energy_dom + energy_imp, P = P_dom + P_imp,
         import_share_P = P_imp / P)
stopifnot(nrow(d_o) == nrow(df))
# consistency with the CF_paid sample: rebuilt paid hours must match CF_paid * P / 50
chk <- d_o %>% inner_join(cf_paid_df33 %>% transmute(country, h_cf = exp(log_cf_paid + log_protein) / 50), by = "country")
cat(sprintf("\nConsistency: max |paid_dom + paid_imp - CF_paid*P/50| = %.2e hr/cap/day (n = %d)\n",
            max(abs(chk$paid - chk$h_cf)), nrow(chk)))

el <- function(y, x, extra = NULL) {
  f <- lm(reformulate(c(x, extra, X), y), d_o %>% mutate(across(c(paid_dom, paid_imp, energy_dom, energy_imp,
                                                                  paid, energy, P, P_dom, P_imp), log, .names = "log_{.col}")))
  s <- summary(f)$coefficients
  sprintf("%.2f (%s)", s[x, "Estimate"], fmt_p(s[x, "Pr(>|t|)"]))
}
contrib <- function(tot, dom, imp) {
  r  <- function(v) resid(lm(reformulate(X, v), d_o))
  rt <- r(tot)
  c(dom = cov(r(dom), rt) / var(rt), imp = cov(r(imp), rt) / var(rt))
}
origin_tab <- map_dfr(c("paid", "energy"), function(w) {
  cb <- contrib(w, paste0(w, "_dom"), paste0(w, "_imp"))
  tibble(what = w,
         `mean dom / imp per cap-day` = sprintf("%.3f / %.3f", mean(d_o[[paste0(w, "_dom")]]), mean(d_o[[paste0(w, "_imp")]])),
         `dom vs own protein`   = el(paste0("log_", w, "_dom"), "log_P_dom"),
         `imp vs own protein`   = el(paste0("log_", w, "_imp"), "log_P_imp"),
         `total vs total protein` = el(paste0("log_", w), "log_P"),
         `share of spread: dom / imp` = sprintf("%.0f%% / %.0f%%", 100 * cb[["dom"]], 100 * cb[["imp"]]),
         `import share of protein | P, GDP` = el(paste0("log_", w), "import_share_P", "log_P"))
})
cat(sprintf("\n---- Paid hours and energy per capita by product origin, net of GDP (EXIO-only, n = %d) ----\n", nrow(d_o)))
print(origin_tab, n = Inf, width = Inf)
cat("\nRead: 'vs own protein' ~ 1 -> that origin's hours/energy follow its protein; ~ 0 -> they do not.",
    "\n      'share of spread' -> which origin carries the cross-country variation in the total at a given",
    "\n      income (levels, sums to 100%). 'import share of protein' -> change in log hours/energy per",
    "\n      unit import share (0-1), holding total protein and GDP fixed.\n")
write_csv(origin_tab, "results/paid_energy_by_origin.csv")
cat("Wrote results/paid_energy_by_origin.csv\n")


#### Are per-gram intensities driven by how protein-dense the food basket is? ####
# In the import-effort trade-off plots (tradeoff_convfac_protein_*_import_effort_*),
# countries line up: import time and import energy per 50 g of imported protein rise
# together, ~8-fold across countries. Paid hours and energy are attached per tonne of
# ALL food, so a basket with more non-protein food per gram of protein (fruit, coffee,
# sugar, wine) sits high on both axes at once. Test: per origin, regress log hours (or
# energy) per 50 g protein on log kcal per g protein of the same basket.
#   elasticity ~ 1, high R2 -> per-gram intensity is mostly the basket's protein density
#   elasticity ~ 0          -> intensity differences come from elsewhere (sourcing, methods)
# All EXIO-modeled countries (not RoW members), no income control: this is descriptive,
# like the plots. Product origin as in the block above.

if (!exists("summary_kcal_df_long")) {
  stop("summary_kcal_df_long missing: run 2.analyze_result.R in this session, or re-run it and ",
       "2.2.save_objects_for_9.R so data/objects_for_9_2020.RData carries it.")
}
kcal_origin <- summary_kcal_df_long %>%
  filter(cat %in% c("domestic", "import")) %>%
  mutate(country = as.character(country), origin = ifelse(cat == "domestic", "dom", "imp")) %>%
  group_by(country, origin) %>%
  summarise(kcal = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

d_den <- origin_long %>%
  inner_join(pro_origin, by = c("country", "origin")) %>%
  inner_join(kcal_origin, by = c("country", "origin")) %>%
  filter(!country %in% row_countries, value > 0, P > 0, kcal > 0) %>%
  mutate(per_50g = value / P * 50, kcal_per_g = kcal / P)

den_tab <- d_den %>%
  group_by(what, origin) %>%
  group_modify(function(g, key) {
    f <- lm(log(per_50g) ~ log(kcal_per_g), g)
    s <- summary(f)$coefficients
    tibble(n = nrow(g),
           `kcal per g protein: range` = sprintf("%.0f-%.0f", min(g$kcal_per_g), max(g$kcal_per_g)),
           `per 50 g: max/min` = sprintf("%.1f", max(g$per_50g) / min(g$per_50g)),
           `elasticity on kcal/g` = sprintf("%.2f (%s)", s[2, "Estimate"], fmt_p(s[2, "Pr(>|t|)"])),
           R2 = sprintf("%.2f", summary(f)$r.squared))
  }) %>% ungroup()
cat("\n---- Per-gram intensity vs protein density of the basket, by origin (EXIO-modeled countries) ----\n")
print(den_tab, n = Inf, width = Inf)

# The plots' line: time and energy per gram moving together. If both are mostly
# basket density, hours per MJ should vary far less across countries than either axis.
ratio <- d_den %>% select(country, origin, what, per_50g) %>%
  pivot_wider(names_from = what, values_from = per_50g) %>%
  mutate(hr_per_mj = paid / energy)
cat("\nAcross countries, max/min of: time per 50 g | energy per 50 g | hours per MJ\n")
ratio %>% group_by(origin) %>%
  summarise(time = max(paid) / min(paid), energy = max(energy) / min(energy),
            hr_per_mj = max(hr_per_mj) / min(hr_per_mj), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ round(.x, 1))) %>% print()
cat("\nRead: elasticity ~ 1 with high R2 -> the per-gram spread is mostly how much non-protein food",
    "\n      comes with each gram of protein; hours per MJ varying much less than either axis means",
    "\n      time and energy are carried by the same tonnes.\n")
write_csv(den_tab, "results/intensity_vs_protein_density.csv")
cat("Wrote results/intensity_vs_protein_density.csv\n")


#### Same test with tonnes: food mass per gram of protein ####
# Calories per gram did not explain the per-gram spread. Hours and energy are
# attached per tonne, and water-heavy foods (vegetables, fruit, some dairy) are many
# tonnes but few calories, so mass per gram of protein is the measure that matches
# the numerator. mass_origin_df comes from 2.analyze_result.R; if this session ran
# an older version of it, rebuild it here from the same session objects.

if (!exists("mass_origin_df")) {
  if (!all(vapply(c("FABIO_y_hh", "coeff_cal", "agg_country_footprint", "country_summary"), exists, logical(1)))) {
    stop("mass_origin_df missing: run 2.analyze_result.R in this session, or re-run it and ",
         "2.2.save_objects_for_9.R so data/objects_for_9_2020.RData carries it.")
  }
  mass_origin_df <- country_summary(agg_country_footprint(
      Matrix::Diagonal(x = as.numeric(coeff_cal$kcal_per_kg > 0)) %*% FABIO_y_hh)) %>%
    transmute(country = as.character(country),
              dom = domestic_per_capita / 1e6 / 365 * 1000,   # kg/cap/day
              imp = import_per_capita   / 1e6 / 365 * 1000) %>%
    pivot_longer(c(dom, imp), names_to = "origin", values_to = "kg")
}

d_mass <- d_den %>%
  inner_join(mass_origin_df, by = c("country", "origin")) %>%
  filter(kg > 0) %>%
  mutate(kg_per_g = kg / P)

mass_tab <- d_mass %>%
  group_by(what, origin) %>%
  group_modify(function(g, key) {
    f1 <- lm(log(per_50g) ~ log(kg_per_g), g)
    f2 <- lm(log(per_50g) ~ log(kg_per_g) + log(kcal_per_g), g)
    s1 <- summary(f1)$coefficients; s2 <- summary(f2)$coefficients
    tibble(n = nrow(g),
           `kg per g protein: range` = sprintf("%.3f-%.3f", min(g$kg_per_g), max(g$kg_per_g)),
           `elasticity on kg/g` = sprintf("%.2f (%s)", s1[2, "Estimate"], fmt_p(s1[2, "Pr(>|t|)"])),
           R2 = sprintf("%.2f", summary(f1)$r.squared),
           `kg/g | kcal/g` = sprintf("%.2f (%s)", s2[2, "Estimate"], fmt_p(s2[2, "Pr(>|t|)"])),
           `R2 (kg + kcal)` = sprintf("%.2f", summary(f2)$r.squared))
  }) %>% ungroup()
cat("\n---- Per-gram intensity vs food mass per gram of protein, by origin (EXIO-modeled countries) ----\n")
print(mass_tab, n = Inf, width = Inf)
cat("\nRead: elasticity on kg/g ~ 1 with high R2 -> the per-gram spread is mostly how many tonnes of food",
    "\n      come with each gram of protein (water-heavy baskets). Low R2 here as well -> the spread",
    "\n      comes from sourcing (which producers) or per-tonne intensities, not basket mass.\n")
write_csv(mass_tab, "results/intensity_vs_mass_density.csv")
cat("Wrote results/intensity_vs_mass_density.csv\n")


#### Import premium per gram of protein vs per tonne of food ####
# Findings: 26 of 42 countries import protein at a higher time CF than they produce
# it at, and the import/domestic ratio rises with income (9q). The blocks above show
# import intensity per gram is mostly food mass per gram, so rich countries' premium
# could be WHAT they import (protein-poor, water-heavy food by the tonne) rather than
# HOW it is produced. Same construction as 9q (streams = domestically produced vs
# imported food; paid hours hr_m + hr_f and energy summed over food + non-food
# sectors and effort origin; countries with > 2% of protein imported), once per
# gram of protein (reproduces 9q) and once per tonne of food.
#   per-tonne ratio ~ 1, no income slope -> the premium is diet composition
#   per-tonne ratio > 1, rising with income -> imported food itself carries more hours

eff_stream <- effort_consumption_df %>%
  mutate(country = region_to_iso_1to1[as.character(exio_region)]) %>%
  filter(!is.na(country), type %in% c("hr_m", "hr_f", "en")) %>%
  mutate(origin = ifelse(protein_source == "domestic", "dom", "imp"),
         what   = ifelse(type == "en", "energy", "paid"),
         value  = ifelse(type == "en", per_capita_value / 365, per_capita_value)) %>%   # MJ or hr / cap / day
  group_by(country, what, origin) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  inner_join(pro_origin, by = c("country", "origin")) %>%
  inner_join(mass_origin_df, by = c("country", "origin")) %>%
  filter(P > 0, kg > 0) %>%
  mutate(per_g = value / P * 50, per_t = value / kg * 1000)             # per 50 g protein; per tonne of food

imp_share_9y <- pro_origin %>% pivot_wider(names_from = origin, values_from = P) %>%
  transmute(country, import_share = imp / (dom + imp))
prem <- eff_stream %>%
  select(country, what, origin, per_g, per_t) %>%
  pivot_wider(names_from = origin, values_from = c(per_g, per_t)) %>%
  mutate(ratio_g = per_g_imp / per_g_dom, ratio_t = per_t_imp / per_t_dom) %>%
  inner_join(imp_share_9y, by = "country") %>%
  filter(import_share > 0) %>%                      # no 2% floor: it only removed India (see 9q)
  inner_join(wdi_latest %>% select(country = iso3c, gdp_pcap_ppp), by = "country") %>%
  filter(is.finite(ratio_g), is.finite(ratio_t))

prem_tab <- prem %>%
  group_by(what) %>%
  group_modify(function(g, key) {
    sl <- function(r) {
      s <- summary(lm(log(g[[r]]) ~ log(g$gdp_pcap_ppp)))$coefficients
      sprintf("%.2f (%s)", s[2, "Estimate"], fmt_p(s[2, "Pr(>|t|)"]))
    }
    tibble(n = nrow(g),
           `per g: # ratio > 1` = sum(g$ratio_g > 1), `per g: median` = sprintf("%.2f", median(g$ratio_g)),
           `per g: income slope` = sl("ratio_g"),
           `per t: # ratio > 1` = sum(g$ratio_t > 1), `per t: median` = sprintf("%.2f", median(g$ratio_t)),
           `per t: income slope` = sl("ratio_t"))
  }) %>% ungroup()
cat("\n---- Import/domestic intensity ratio: per 50 g protein (as 9q) vs per tonne of food ----\n")
print(prem_tab, n = Inf, width = Inf)
cat("(per-gram columns should match 9q's stream ratios, now including India)\n")
cat("\nRead: per-tonne ratio near 1 with no income slope -> rich countries' import premium per gram is",
    "\n      diet composition (protein-poor food by the tonne); per-tonne ratio above 1 and rising with",
    "\n      income -> imported food itself carries more hours or energy per tonne.\n")
# Which side drives the per-tonne ratio across countries? log ratio = log imp - log dom, so
# its variance splits into a domestic part (var(dom) - cov) and an import part (var(imp) - cov).
cat("\n---- Per-tonne ratio: share of its cross-country variance from the domestic vs the import side ----\n")
print(prem %>% group_by(what) %>%
        summarise(n = n(),
                  dom_max_min = max(per_t_dom) / min(per_t_dom), imp_max_min = max(per_t_imp) / min(per_t_imp),
                  domestic_share = (var(log(per_t_dom)) - cov(log(per_t_imp), log(per_t_dom))) / var(log(ratio_t)),
                  import_share   = (var(log(per_t_imp)) - cov(log(per_t_imp), log(per_t_dom))) / var(log(ratio_t)),
                  .groups = "drop") %>%
        mutate(across(where(is.double), ~ round(.x, 2))))
write_csv(prem, "results/import_premium_per_gram_vs_per_tonne.csv")
cat("Wrote results/import_premium_per_gram_vs_per_tonne.csv (per-country ratios)\n")


#### Does the elasticity of hours / energy per person on protein differ by income? ####
# Table 3's specification applied to Table 4: log H ~ log P * log GDP. Reports e at India's
# and the Netherlands' income (delta method), the interaction's t-test p and its case-
# resampling bootstrap (the Methods' check for interactions at small n), and the countries
# with the largest influence (Cook's distance) on each interaction. Draft, low-income
# section footnote: no difference is significant; the lowest t-test p (paid hours, ~0.07)
# does not hold up under the bootstrap.

GDP_AT_9Y <- c(7187, 64862)
inc_row <- function(d, lh, what, sample, R = R_BOOT) {
  d <- d %>% mutate(lh = lh) %>% filter(is.finite(lh), is.finite(.data[[Y]]), is.finite(.data[[X]]))
  fml <- as.formula(paste("lh ~", Y, "*", X)); int <- paste0(Y, ":", X)
  fit <- lm(fml, d); s <- summary(fit)$coefficients; V <- vcov(fit)
  e_at <- function(g) {
    x <- log(g); est <- s[Y, 1] + s[int, 1] * x
    se <- sqrt(V[Y, Y] + x^2 * V[int, int] + 2 * x * V[Y, int])
    sprintf("%.2f (%s)", est, fmt_p(2 * pt(-abs(est / se), fit$df.residual)))
  }
  b <- replicate(R, coef(lm(fml, d[sample.int(nrow(d), replace = TRUE), ]))[[int]])
  cook <- sort(setNames(cooks.distance(fit), d$country), decreasing = TRUE)[1:3]
  tibble(what, sample, n = nrow(d), e_india = e_at(GDP_AT_9Y[1]), e_netherlands = e_at(GDP_AT_9Y[2]),
         interaction = sprintf("%.2f", s[int, 1]), t_p = fmt_p(s[int, 4]),
         boot_ci = sprintf("[%.2f, %.2f]", quantile(b, .025, na.rm = TRUE), quantile(b, .975, na.rm = TRUE)),
         boot_p = fmt_p(2 * min(mean(b < 0, na.rm = TRUE), mean(b > 0, na.rm = TRUE))),
         most_influential = paste(sprintf("%s %.2f", names(cook), cook), collapse = ", "))
}
inc_rows <- function(d_all, d_paid, d_en, sample) {
  a <- d_all %>% select(country, log_cf, all_of(c(Y, X))) %>%
    inner_join(d_paid %>% select(country, log_cf_paid), by = "country") %>%
    filter(if_all(where(is.numeric), is.finite))
  P <- exp(a[[Y]]); h_all <- exp(a$log_cf) * P / 50; h_paid <- exp(a$log_cf_paid) * P / 50
  e <- d_en %>% select(country, log_cf_energy, all_of(c(Y, X)))
  bind_rows(inc_row(a, log(h_all), "all work", sample),
            inc_row(a, log(h_all - h_paid), "unpaid", sample),
            inc_row(a, log(h_paid), "paid", sample),
            inc_row(e, e$log_cf_energy + e[[Y]] - log(50), "energy", sample))
}
set.seed(20261002)
inc_tab <- bind_rows(inc_rows(df, cf_paid_df33, df2, "EXIO-only"),
                     inc_rows(df_collapsed, cf_paid_df38, df2_collapsed, "RoW-collapsed"))
cat("\n---- Elasticity of hours / energy per person on protein, by income (net of GDP) ----\n")
print(inc_tab, n = Inf, width = Inf)
cat("\nRead: e at each income with delta-method p; 'interaction' tests whether e changes with income.",
    "\n      Trust the bootstrap over the t-test at this n; 'most_influential' shows which countries",
    "\n      the interaction rests on (Cook's distance).\n")
write_csv(inc_tab, "results/hours_protein_elasticity_by_income.csv")
cat("Wrote results/hours_protein_elasticity_by_income.csv\n")
