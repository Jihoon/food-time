#### Does CF limit protein supply independent of income, or is that riding on GDP? ####
#
# Motivating idea: daily time is a hard 24h budget, so (as with the stable
# travel-time-budget literature) societies settle on a roughly stable time
# budget for provisioning food too. Under a fixed budget, a country's time
# conversion factor (CF = hr / 50 g protein) caps the protein it can
# provision: low-CF (fast) societies get more protein per hour spent;
# high-CF (slow) societies get less, UNLESS income lets them buy more hours'
# worth of effort (e.g. via purchased/mechanized food supply). So CF and
# protein supply should correlate even after partialling out income --
# otherwise the apparent CF -> capability-set link would just be both
# variables riding on GDP.
#
# Run AFTER 2.analyze_result.R (reuses its `tradeoff_protein_allwork_consump`
# object, already in the R environment).
#
# Restricted to EXIO-individually-modeled countries (is_row == FALSE). Any
# FABIO country that falls in a RoW aggregate has its economic labor/energy
# intensity pasted from that aggregate's single EXIO value (see CLAUDE.md,
# "Architecture: The Bridging Problem" / reorder_countries_to_FABIO()) --
# its CF is not an independent estimate, so including it would just inflate
# n with copies of a handful of regional values.

library(tidyverse)
library(data.table)
library(WDI)

#### 1. Country-level CF and protein supply (EXIO-modeled countries only) ####
# Total labor (household + economic, both genders) per 50 g of protein
# actually consumed (domestic + import), and the matching protein supply.

cf_df <- tradeoff_protein_allwork_consump %>%
  filter(!is_row) %>%
  group_by(country) %>%
  summarise(hr_per_50g_protein = sum(hr_per_50g_protein, na.rm = TRUE),
            g_protein_per_cap_day = first(g_protein_per_cap_day),
            .groups = "drop")

#### 2. Income / labor-productivity controls (World Bank WDI) ####

wdi_raw <- WDI(country = "all",
               indicator = c(gdp_pcap_ppp   = "NY.GDP.PCAP.PP.KD",   # GDP per capita, PPP (constant 2021 int'l $)
                             gdp_per_worker = "SL.GDP.PCAP.EM.KD"),  # GDP per person employed (constant 2021 PPP $) -- labor-productivity proxy
               start = year - 5, end = year, extra = TRUE)

# WDI reports with a lag for some countries -- take each country's latest
# non-NA value within the window rather than requiring an exact `year` match.
latest_nonNA <- function(x, yr) {
  ok <- !is.na(x)
  if (!any(ok)) return(NA_real_)
  x[ok][which.max(yr[ok])]
}

wdi_latest <- wdi_raw %>%
  filter(region != "Aggregates") %>%   # drop WB's own region/income-group rollups
  group_by(iso3c) %>%
  summarise(gdp_pcap_ppp   = latest_nonNA(gdp_pcap_ppp, year),
            gdp_per_worker = latest_nonNA(gdp_per_worker, year),
            .groups = "drop")

#### 3. Merge ####

df <- cf_df %>%
  inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
  filter(hr_per_50g_protein > 0, g_protein_per_cap_day > 0) %>%
  mutate(log_cf         = log(hr_per_50g_protein),
         log_protein    = log(g_protein_per_cap_day),
         log_gdp_pcap   = log(gdp_pcap_ppp),
         log_gdp_worker = log(gdp_per_worker))

cat(sprintf("Countries with CF, protein supply, and >=1 income control: %d\n", nrow(df)))

#### 4. Partial correlation: CF vs. protein supply, controlling for income ####
# Partial correlation of x and y controlling for z = correlation of what's
# left of x and y after each is regressed on z. Equivalent to the standard
# Pearson partial-correlation formula for one control variable. Logged
# throughout since GDP, protein supply and CF are all right-skewed.

partial_cor <- function(x, y, z) {
  keep <- is.finite(x) & is.finite(y) & is.finite(z)
  rx <- resid(lm(x[keep] ~ z[keep]))
  ry <- resid(lm(y[keep] ~ z[keep]))
  ct <- cor.test(rx, ry)
  list(r = unname(ct$estimate), p = ct$p.value, n = sum(keep))
}

zero_order    <- cor.test(df$log_cf, df$log_protein)
pc_gdp_pcap   <- partial_cor(df$log_cf, df$log_protein, df$log_gdp_pcap)
pc_gdp_worker <- partial_cor(df$log_cf, df$log_protein, df$log_gdp_worker)

cat("\n---- CF (hr/50g protein) vs. protein supply (g/cap/day), log-log ----\n")
cat(sprintf("Zero-order r                              = %.3f (p = %.3g, n = %d)\n",
            zero_order$estimate, zero_order$p.value, sum(is.finite(df$log_cf) & is.finite(df$log_protein))))
cat(sprintf("Partial r | GDP per capita (PPP)           = %.3f (p = %.3g, n = %d)\n",
            pc_gdp_pcap$r, pc_gdp_pcap$p, pc_gdp_pcap$n))
cat(sprintf("Partial r | GDP per worker (labor prod.)   = %.3f (p = %.3g, n = %d)\n",
            pc_gdp_worker$r, pc_gdp_worker$p, pc_gdp_worker$n))
cat(sprintf("\n(context: r[CF, GDPpc] = %.3f | r[protein, GDPpc] = %.3f)\n",
            cor(df$log_cf, df$log_gdp_pcap, use = "complete.obs"),
            cor(df$log_protein, df$log_gdp_pcap, use = "complete.obs")))

# Optional cross-check against ppcor's closed-form partial correlation, if installed.
if (requireNamespace("ppcor", quietly = TRUE)) {
  chk <- ppcor::pcor.test(df$log_cf, df$log_protein, df$log_gdp_pcap)
  cat(sprintf("[ppcor cross-check] r = %.3f, p = %.3g\n", chk$estimate, chk$p.value))
}

#### 5. Visualize the partial relationship (GDP per capita as control) ####

df <- df %>%
  mutate(resid_cf      = resid(lm(log_cf ~ log_gdp_pcap)),
         resid_protein = resid(lm(log_protein ~ log_gdp_pcap)))

p_partial <- ggplot(df, aes(resid_cf, resid_protein)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "grey40", linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = country), size = 3, max.overlaps = 15) +
  labs(x = "CF residual (log hr/50g protein, GDP/cap partialled out)",
       y = "Protein supply residual (log g/cap/day, GDP/cap partialled out)",
       title = "Partial relationship: CF vs. protein supply, controlling for GDP per capita (PPP)") +
  theme_minimal()
print(p_partial)
ggsave("results/partial_cf_vs_protein_control_gdp.pdf", p_partial, width = 8, height = 6)

fwrite(df %>% select(country, hr_per_50g_protein, g_protein_per_cap_day, gdp_pcap_ppp, gdp_per_worker),
       "output/cf_protein_income_control_data.csv")

#### 6. Mediation decomposition: does GDP work THROUGH CF, or around it? ####
# The capability-set story puts CF on the causal path, not beside it as a
# confound: income -> mechanization -> lower CF -> more protein affordable
# within a fixed time budget. So instead of asking whether CF "survives"
# controlling for income, decompose income's total effect on protein supply
# (GDP -> protein) into:
#   direct effect   (ADE) : GDP -> protein supply, NOT through CF
#   indirect effect (ACME): GDP -> CF -> protein supply
# via the classic product-of-coefficients method (Baron & Kenny 1986), with
# a Sobel (1982) test for whether the indirect effect differs from zero.

mediation_decomp <- function(df, x_col, m_col, y_col) {
  keep <- complete.cases(df[c(x_col, m_col, y_col)])
  x <- df[[x_col]][keep]; m <- df[[m_col]][keep]; y <- df[[y_col]][keep]

  fit_total <- lm(y ~ x)      # total effect of X (income) on Y (protein supply)
  fit_med   <- lm(m ~ x)      # effect of X on the mediator M (CF)
  fit_out   <- lm(y ~ x + m)  # effect of M on Y, controlling for X -- also gives the direct effect of X

  a <- unname(coef(fit_med)["x"]);  se_a <- summary(fit_med)$coefficients["x", "Std. Error"]
  b <- unname(coef(fit_out)["m"]);  se_b <- summary(fit_out)$coefficients["m", "Std. Error"]

  total    <- unname(coef(fit_total)["x"])
  direct   <- unname(coef(fit_out)["x"])     # ADE
  indirect <- a * b                          # ACME: total - direct, algebraically identical for OLS
  se_indirect <- sqrt(b^2 * se_a^2 + a^2 * se_b^2)  # Sobel (1982) first-order delta-method SE
  z <- indirect / se_indirect

  list(n = length(x), total = total, direct = direct, indirect = indirect,
       prop_mediated = indirect / total, sobel_z = z, sobel_p = 2 * (1 - pnorm(abs(z))))
}

med_gdp_pcap   <- mediation_decomp(df, "log_gdp_pcap",   "log_cf", "log_protein")
med_gdp_worker <- mediation_decomp(df, "log_gdp_worker", "log_cf", "log_protein")

cat("\n---- Mediation decomposition: GDP -> CF -> protein supply ----\n")
cat(sprintf("Via GDP per capita (PPP):         total = %.3f | direct (ADE) = %.3f | indirect via CF (ACME) = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g), n = %d\n",
            med_gdp_pcap$total, med_gdp_pcap$direct, med_gdp_pcap$indirect,
            med_gdp_pcap$prop_mediated * 100, med_gdp_pcap$sobel_z, med_gdp_pcap$sobel_p, med_gdp_pcap$n))
cat(sprintf("Via GDP per worker (labor prod.): total = %.3f | direct (ADE) = %.3f | indirect via CF (ACME) = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g), n = %d\n",
            med_gdp_worker$total, med_gdp_worker$direct, med_gdp_worker$indirect,
            med_gdp_worker$prop_mediated * 100, med_gdp_worker$sobel_z, med_gdp_worker$sobel_p, med_gdp_worker$n))

# Optional cross-check with bootstrapped ACME/ADE (Imai et al.), if installed --
# handles the small-sample non-normality of the indirect effect's distribution
# better than the Sobel test's normal approximation above.
if (requireNamespace("mediation", quietly = TRUE)) {
  set.seed(1)
  fit_m <- lm(log_cf ~ log_gdp_pcap, data = df)
  fit_y <- lm(log_protein ~ log_gdp_pcap + log_cf, data = df)
  med_boot <- mediation::mediate(fit_m, fit_y, treat = "log_gdp_pcap", mediator = "log_cf",
                                  boot = TRUE, sims = 1000)
  cat("\n[mediation package cross-check, GDP per capita, bootstrapped]\n")
  print(summary(med_boot))
}

fwrite(data.frame(control = c("gdp_pcap_ppp", "gdp_per_worker"),
                   rbind(as.data.frame(med_gdp_pcap[c("n", "total", "direct", "indirect", "prop_mediated", "sobel_z", "sobel_p")]),
                         as.data.frame(med_gdp_worker[c("n", "total", "direct", "indirect", "prop_mediated", "sobel_z", "sobel_p")]))),
       "output/mediation_decomposition.csv")

#### 7. Robustness: collapse each RoW aggregate to ONE population-weighted point ####
# The all-countries run (n=54) pastes each RoW aggregate's single CF value
# across ~20 countries -- a handful of leverage points at the low-income tail,
# not ~20 independent ones. Dropping RoW entirely (n=33, EXIO-modeled
# countries only) instead discards the low/mid-income tail of the world.
# Middle ground: one row per RoW aggregate, built from a population-weighted
# average of protein supply and GDP across whichever of ITS FABIO members
# have usable GHD + WDI data (paired with the aggregate's one CF value).
#
# Caveat this does NOT fix: FABIO_reg assigns every non-EXIO-modeled country
# to one of 5 RoW aggregates, but most members lack GHD time-use data, WDI
# coverage, or both -- so even this collapsed point is built from a subset,
# not the aggregate's full membership, and that subset likely skews toward
# whichever countries have better statistical infrastructure (correlated
# with income). `pop_coverage_pct` reports how much of each aggregate's
# total population that subset represents, so this is a visible number
# rather than an assumption.

row_lookup <- FABIO_reg %>% transmute(country = ISO, exio_region = EXIOBASE)

row_full_pop <- row_lookup %>%
  filter(grepl("RoW", exio_region)) %>%
  left_join(pop_data_yr, by = c("country" = "iso3c")) %>%
  group_by(exio_region) %>%
  summarise(pop_total_region = sum(population, na.rm = TRUE), .groups = "drop")

row_collapsed <- tradeoff_protein_allwork_consump %>%
  filter(is_row) %>%
  group_by(country) %>%
  summarise(hr_per_50g_protein = sum(hr_per_50g_protein, na.rm = TRUE),
            g_protein_per_cap_day = first(g_protein_per_cap_day),
            .groups = "drop") %>%
  inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
  filter(hr_per_50g_protein > 0, g_protein_per_cap_day > 0) %>%
  left_join(row_lookup, by = "country") %>%
  left_join(pop_data_yr, by = c("country" = "iso3c")) %>%
  filter(!is.na(population), population > 0) %>%
  group_by(exio_region) %>%
  summarise(n_members_used = n(),
            pop_used = sum(population),
            hr_per_50g_protein    = unique(hr_per_50g_protein),  # identical within an aggregate by construction
            g_protein_per_cap_day = weighted.mean(g_protein_per_cap_day, population),
            gdp_pcap_ppp          = weighted.mean(gdp_pcap_ppp, population, na.rm = TRUE),
            gdp_per_worker        = weighted.mean(gdp_per_worker, population, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(row_full_pop, by = "exio_region") %>%
  mutate(pop_coverage_pct = pop_used / pop_total_region * 100)

cat("\n---- RoW-aggregate collapse: population coverage of the subset used ----\n")
print(row_collapsed %>% select(exio_region, n_members_used, pop_coverage_pct) %>%
        mutate(pop_coverage_pct = round(pop_coverage_pct, 1)))

df_collapsed <- bind_rows(
  df %>% select(country, hr_per_50g_protein, g_protein_per_cap_day, gdp_pcap_ppp, gdp_per_worker),
  row_collapsed %>% select(country = exio_region, hr_per_50g_protein, g_protein_per_cap_day, gdp_pcap_ppp, gdp_per_worker)
) %>%
  mutate(log_cf         = log(hr_per_50g_protein),
         log_protein    = log(g_protein_per_cap_day),
         log_gdp_pcap   = log(gdp_pcap_ppp),
         log_gdp_worker = log(gdp_per_worker))

cat(sprintf("\nCollapsed sample: %d EXIO-modeled countries + %d RoW-aggregate points = %d rows\n",
            nrow(df), nrow(row_collapsed), nrow(df_collapsed)))

zero_order_c    <- cor.test(df_collapsed$log_cf, df_collapsed$log_protein)
pc_gdp_pcap_c   <- partial_cor(df_collapsed$log_cf, df_collapsed$log_protein, df_collapsed$log_gdp_pcap)
pc_gdp_worker_c <- partial_cor(df_collapsed$log_cf, df_collapsed$log_protein, df_collapsed$log_gdp_worker)

cat("\n---- [RoW-collapsed] CF vs. protein supply, log-log ----\n")
cat(sprintf("Zero-order r                              = %.3f (p = %.3g, n = %d)\n",
            zero_order_c$estimate, zero_order_c$p.value, nrow(df_collapsed)))
cat(sprintf("Partial r | GDP per capita (PPP)           = %.3f (p = %.3g, n = %d)\n",
            pc_gdp_pcap_c$r, pc_gdp_pcap_c$p, pc_gdp_pcap_c$n))
cat(sprintf("Partial r | GDP per worker (labor prod.)   = %.3f (p = %.3g, n = %d)\n",
            pc_gdp_worker_c$r, pc_gdp_worker_c$p, pc_gdp_worker_c$n))

med_gdp_pcap_c   <- mediation_decomp(df_collapsed, "log_gdp_pcap",   "log_cf", "log_protein")
med_gdp_worker_c <- mediation_decomp(df_collapsed, "log_gdp_worker", "log_cf", "log_protein")

cat("\n---- [RoW-collapsed] Mediation decomposition: GDP -> CF -> protein supply ----\n")
cat(sprintf("Via GDP per capita (PPP):         total = %.3f | direct (ADE) = %.3f | indirect via CF (ACME) = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g), n = %d\n",
            med_gdp_pcap_c$total, med_gdp_pcap_c$direct, med_gdp_pcap_c$indirect,
            med_gdp_pcap_c$prop_mediated * 100, med_gdp_pcap_c$sobel_z, med_gdp_pcap_c$sobel_p, med_gdp_pcap_c$n))
cat(sprintf("Via GDP per worker (labor prod.): total = %.3f | direct (ADE) = %.3f | indirect via CF (ACME) = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g), n = %d\n",
            med_gdp_worker_c$total, med_gdp_worker_c$direct, med_gdp_worker_c$indirect,
            med_gdp_worker_c$prop_mediated * 100, med_gdp_worker_c$sobel_z, med_gdp_worker_c$sobel_p, med_gdp_worker_c$n))

fwrite(row_collapsed %>% select(exio_region, n_members_used, pop_used, pop_total_region, pop_coverage_pct,
                                hr_per_50g_protein, g_protein_per_cap_day, gdp_pcap_ppp, gdp_per_worker),
       "output/row_aggregate_collapse_coverage.csv")
