#### Two CFs, same denominator: total labor vs. domestic-only labor ####
#### going into this country's total protein consumption ####
#
# CF_1 (already established, script 9): all labor -- domestic-origin AND
# import-origin -- embodied in this country's total protein consumption,
# over that total consumption. "PS performance as a whole": how
# efficiently does everything that feeds this country, wherever the
# labor came from, convert time into its nutrition. This is
# tradeoff_protein_allwork_consump; nothing new to compute here.
#
# CF_2 (new): ONLY domestic-origin labor -- wherever it ends up embodied,
# including domestic labor exported, processed abroad, and re-imported
# (the footnote's second mechanism) -- over the SAME total protein
# consumption denominator. "PS performance with respect to domestic
# input": what the country's own labor buys toward its own nutrition,
# net of any assist from importing others' labor. Needs the
# effort_origin split (domestic vs import, of the labor itself) from
# effort_consumption_df, the same source 9q used for df_mosaic -- NOT
# tradeoff_protein_allwork, which splits by consumption origin
# (domestic/imported protein), not labor origin.
#
# Full run order, in one session:
#   1. 2.analyze_result.R       (effort_consumption_df, tradeoff_protein_
#                                 allwork_consump, region_to_iso_1to1)
#   2. 9.capability_set_income_control.R (wdi_latest, mediation_decomp(),
#                                 partial_cor(), row_lookup, row_full_pop,
#                                 pop_data_yr -- from its RoW-collapse block)
#   3. 9m.energy_time_dual_mediator.R    (builds df2 -- 9o needs this even
#                                 though 9s itself does not use df2 directly)
#   4. 9o.bootstrap_mediation_and_interaction.R (boot_stat(), report_boot(), R)
#   5. 9s.mediation_two_cf_comparison.R  (this file)
# Independent of 9q/9r otherwise -- rebuilds what it needs directly so it
# doesn't depend on their run order.

library(tidyverse)

# ---- Total protein consumption denominator (same for both CFs) ----
protein_total_9s <- tradeoff_protein_allwork_consump %>%
  filter(!is_row) %>%
  group_by(country) %>%
  summarise(g_protein_per_cap_day = first(g_protein_per_cap_day), .groups = "drop")

# ---- CF_2 numerator: domestic-origin labor, summed across BOTH
# consumption streams (domestic_per_capita and import_per_capita) and
# both sectors -- "all domestic effort," not restricted to where it ends
# up embodied. ----
domestic_effort_9s <- effort_consumption_df %>%
  mutate(country = region_to_iso_1to1[as.character(exio_region)]) %>%
  filter(type %in% c("hr_m", "hr_f"), effort_origin == "domestic") %>%
  group_by(country) %>%
  summarise(hr_per_cap_day_domestic_effort = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

cf2_df <- domestic_effort_9s %>%
  inner_join(protein_total_9s, by = "country") %>%
  mutate(hr_per_50g_protein_cf2 = hr_per_cap_day_domestic_effort / g_protein_per_cap_day * 50) %>%
  inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
  filter(hr_per_50g_protein_cf2 > 0, g_protein_per_cap_day > 0) %>%
  mutate(log_cf2       = log(hr_per_50g_protein_cf2),
         log_protein   = log(g_protein_per_cap_day),
         log_gdp_pcap  = log(gdp_pcap_ppp))

cat(sprintf("CF_2 (domestic-effort-only): n = %d\n", nrow(cf2_df)))

zero_order_2  <- cor.test(cf2_df$log_cf2, cf2_df$log_protein)
pc_gdp_pcap_2 <- partial_cor(cf2_df$log_cf2, cf2_df$log_protein, cf2_df$log_gdp_pcap)

cat("\n---- CF_2 vs. total protein supply, log-log ----\n")
cat(sprintf("Zero-order r                     = %.3f (p = %.3g, n = %d)\n",
            zero_order_2$estimate, zero_order_2$p.value, nrow(cf2_df)))
cat(sprintf("Partial r | GDP per capita (PPP) = %.3f (p = %.3g, n = %d)\n",
            pc_gdp_pcap_2$r, pc_gdp_pcap_2$p, pc_gdp_pcap_2$n))

med_cf2 <- mediation_decomp(cf2_df, "log_gdp_pcap", "log_cf2", "log_protein")

cat("\n---- Mediation decomposition, CF_2: GDP -> domestic-effort CF -> total protein supply ----\n")
cat(sprintf("total = %.3f | direct (ADE) = %.3f | indirect via CF (ACME) = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g), n = %d\n",
            med_cf2$total, med_cf2$direct, med_cf2$indirect,
            med_cf2$prop_mediated * 100, med_cf2$sobel_z, med_cf2$sobel_p, med_cf2$n))

cat("\n---- Path a for CF_2: does GDP predict domestic-effort CF? ----\n")
print(summary(lm(log_cf2 ~ log_gdp_pcap, data = cf2_df)))

cat("\n#### For comparison: CF_1, already established (script 9, unchanged) ####\n")
cat("path a:   a1 = -0.52 (n=33) / -0.46 (n=38, RoW-collapsed)\n")
cat("path b:   partial r = -0.13 (p=0.45, n=38); mediated share 16-24%, Sobel p>0.3 throughout\n")
cat("\nWhat to look at: does CF_2's path a differ meaningfully from CF_1's -0.52?\n")
cat("A much weaker CF_2-vs-income relationship would say income buys efficiency\n")
cat("mostly by mobilizing OTHER countries' labor (importing), not by making a\n")
cat("country's own labor more productive -- a real, reportable distinction, not\n")
cat("just a robustness check.\n")

#### RoW-collapsed (n=38) version, matching the "both country samples" ####
#### convention every other CF variant in this paper reports against.   ####
#
# FIX: effort_consumption_df is region-keyed throughout, including for
# the individually-modeled regions -- region_to_iso_1to1 "works" there
# only because those regions happen to map 1-to-1 to a country (a
# rename, not a real disaggregation). RoW aggregates have no such
# mapping because the underlying data genuinely doesn't resolve below
# the region level; the original (broken) version silently produced 0
# rows for the whole RoW block. Correct pattern, matching script 9's
# row_collapsed / 9p's row_collapsed_energy (which work because THEY
# start from tradeoff_protein_allwork_consump, already pasted to FABIO
# members upstream): pull one domestic-effort value per RoW region
# directly, then paste it onto each FABIO member via row_lookup
# (country -> exio_region) before population-weighting.

row_domestic_effort_region_9s <- effort_consumption_df %>%
  filter(type %in% c("hr_m", "hr_f"), effort_origin == "domestic") %>%
  mutate(exio_region = as.character(exio_region)) %>%
  group_by(exio_region) %>%
  summarise(hr_per_cap_day_domestic_effort = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

row_domestic_effort_9s <- row_lookup %>%
  filter(grepl("RoW", exio_region)) %>%
  inner_join(row_domestic_effort_region_9s, by = "exio_region")

row_protein_9s <- tradeoff_protein_allwork_consump %>%
  filter(is_row) %>%
  group_by(country) %>%
  summarise(g_protein_per_cap_day = first(g_protein_per_cap_day), .groups = "drop")

row_collapsed_cf2 <- row_domestic_effort_9s %>%
  inner_join(row_protein_9s, by = "country") %>%
  inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
  filter(hr_per_cap_day_domestic_effort > 0, g_protein_per_cap_day > 0) %>%
  left_join(pop_data_yr, by = c("country" = "iso3c")) %>%
  filter(!is.na(population), population > 0) %>%
  group_by(exio_region) %>%
  summarise(n_members_used = n(),
            pop_used = sum(population),
            hr_per_cap_day_domestic_effort = weighted.mean(hr_per_cap_day_domestic_effort, population),
            g_protein_per_cap_day = weighted.mean(g_protein_per_cap_day, population),
            gdp_pcap_ppp = weighted.mean(gdp_pcap_ppp, population, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(hr_per_50g_protein_cf2 = hr_per_cap_day_domestic_effort / g_protein_per_cap_day * 50) %>%
  left_join(row_full_pop, by = "exio_region") %>%
  mutate(pop_coverage_pct = pop_used / pop_total_region * 100)

cat("\n---- RoW-aggregate collapse (CF_2): population coverage of the subset used ----\n")
print(row_collapsed_cf2 %>% select(exio_region, n_members_used, pop_coverage_pct) %>%
        mutate(pop_coverage_pct = round(pop_coverage_pct, 1)))

cf2_df_collapsed <- bind_rows(
  cf2_df %>% select(country, hr_per_50g_protein_cf2, g_protein_per_cap_day, gdp_pcap_ppp),
  row_collapsed_cf2 %>% select(country = exio_region, hr_per_50g_protein_cf2, g_protein_per_cap_day, gdp_pcap_ppp)
) %>%
  mutate(log_cf2      = log(hr_per_50g_protein_cf2),
         log_protein  = log(g_protein_per_cap_day),
         log_gdp_pcap = log(gdp_pcap_ppp))

cat(sprintf("\nCollapsed CF_2 sample: %d EXIO-modeled + %d RoW-aggregate points = %d rows\n",
            nrow(cf2_df), nrow(row_collapsed_cf2), nrow(cf2_df_collapsed)))

pc_gdp_pcap_2c <- partial_cor(cf2_df_collapsed$log_cf2, cf2_df_collapsed$log_protein, cf2_df_collapsed$log_gdp_pcap)
cat(sprintf("\n[RoW-collapsed] Partial r | GDP per capita = %.3f (p = %.3g, n = %d)\n",
            pc_gdp_pcap_2c$r, pc_gdp_pcap_2c$p, pc_gdp_pcap_2c$n))

med_cf2_c <- mediation_decomp(cf2_df_collapsed, "log_gdp_pcap", "log_cf2", "log_protein")
cat(sprintf("[RoW-collapsed] Mediation: total = %.3f | direct = %.3f | indirect = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g), n = %d\n",
            med_cf2_c$total, med_cf2_c$direct, med_cf2_c$indirect,
            med_cf2_c$prop_mediated * 100, med_cf2_c$sobel_z, med_cf2_c$sobel_p, med_cf2_c$n))

cat("\n[RoW-collapsed] Path a for CF_2:\n")
print(round(summary(lm(log_cf2 ~ log_gdp_pcap, data = cf2_df_collapsed))$coefficients, 4))

#### Bootstrap check (reuses 9o's boot_stat()/report_boot() -- source 9o ####
#### first). ACME is a product of two coefficients; Sobel's normal      ####
#### approximation is known to fit that badly at small n (see 9o and    ####
#### "Is Energy CF the Missing Channel?" in the draft for why this      ####
#### matters and isn't just a formality).                               ####

acme_cf2_stat <- function(d) {
  fit_a <- lm(log_cf2 ~ log_gdp_pcap, data = d)
  fit_b <- lm(log_protein ~ log_gdp_pcap + log_cf2, data = d)
  unname(coef(fit_a)["log_gdp_pcap"]) * unname(coef(fit_b)["log_cf2"])
}
interaction_free_a_stat <- function(d) unname(coef(lm(log_cf2 ~ log_gdp_pcap, data = d))["log_gdp_pcap"])

point_acme_cf2   <- acme_cf2_stat(cf2_df)
point_acme_cf2_c <- acme_cf2_stat(cf2_df_collapsed)
point_a_cf2      <- interaction_free_a_stat(cf2_df)
point_a_cf2_c    <- interaction_free_a_stat(cf2_df_collapsed)

boot_acme_cf2   <- boot_stat(cf2_df, acme_cf2_stat, R)
boot_acme_cf2_c <- boot_stat(cf2_df_collapsed, acme_cf2_stat, R)
boot_a_cf2      <- boot_stat(cf2_df, interaction_free_a_stat, R)
boot_a_cf2_c    <- boot_stat(cf2_df_collapsed, interaction_free_a_stat, R)

report_boot(boot_acme_cf2, point_acme_cf2, "CF_2 ACME, EXIO-only (n=33)", sobel_p_for_comparison = med_cf2$sobel_p)
report_boot(boot_acme_cf2_c, point_acme_cf2_c, "CF_2 ACME, RoW-collapsed (n=38)", sobel_p_for_comparison = med_cf2_c$sobel_p)
report_boot(boot_a_cf2, point_a_cf2, "CF_2 path a, EXIO-only (n=33)")
report_boot(boot_a_cf2_c, point_a_cf2_c, "CF_2 path a, RoW-collapsed (n=38)")
