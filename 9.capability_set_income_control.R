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
