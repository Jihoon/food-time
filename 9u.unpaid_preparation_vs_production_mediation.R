#### Does unpaid provisioning time split into a "production" piece that ####
#### behaves like CF_paid, and a "meal transformation" piece that doesn't? ####
#
# Motivation: CF_paid's indirect effect on protein is significant (n=38,
# bootstrap-confirmed; 9t), while the all-work CF's is not (9/9o), at any
# sample size. An Alternative CF Definition (draft) offers a mechanism for
# this: unpaid time is disproportionately spent converting already-produced
# ingredients into food, not producing the protein those ingredients
# contain -- so it plausibly carries no nutrition-supply signal, unlike
# paid labor which spans production/processing/distribution.
#
# That mechanism is directly testable, not just plausible: the GHD source
# data (0.mrio_prep.R, summary_food_df_long_with_ghd) keeps THREE separate
# unpaid footprint_type categories -- "preparation_non.econ" (cooking/meal
# prep), "processing_non.econ" (e.g. milling, plucking, preserving -- turning
# a raw output into something edible, the same transformation logic as
# cooking, not production of it), and "growth_collection_non.econ" (growing,
# foraging, collecting food directly -- the only genuinely production-side
# unpaid activity) -- all the way through the pipeline.
# tradeoff_protein_allwork_consump just sums across all of them, which is
# why the split disappeared from every mediation script so far. This script
# rebuilds two narrower CFs straight from summary_food_df_long_with_ghd:
#
#   CF_unpaid_transform = hours, "preparation_non.econ" +
#                                 "processing_non.econ"         / total consumption
#   CF_unpaid_produce   = hours, "growth_collection_non.econ"   / total consumption
#
# Each numerator is that unpaid time component ALONE -- not summed with paid
# labor -- over the same total-consumption denominator (pro_consumption) as
# every other CF in Table 1. That mirrors how CF_paid isolates paid time
# alone and CF_energy isolates energy alone: one component per variant, same
# denominator throughout, so results are directly comparable across variants
# rather than each being its own bespoke ratio.
#
# If the mechanism in the draft is right, CF_unpaid_transform should look
# like CF_time (null path b) and CF_unpaid_produce should look more like
# CF_paid (a real indirect effect) -- which would also mean the operative
# distinction is "production- vs. transformation-adjacent activity," not
# literally "paid vs. unpaid," and the draft's framing would need revising
# to match.
#
# Run AFTER 2.analyze_result.R (summary_food_df_long_with_ghd, cty_ghd,
# pro_consumption, row_countries), 9.capability_set_income_control.R (df,
# wdi_latest, partial_cor(), mediation_decomp(), row_lookup, row_full_pop,
# pop_data_yr), and 9o.bootstrap_mediation_and_interaction.R (boot_stat(),
# report_boot(), R) in the same session.

library(tidyverse)

#### 0. Magnitude check first: is growth_collection big enough in this ####
#### sample to mean anything, or is it near-zero noise for most        ####
#### countries (subsistence food-growing/foraging is not universal)?   ####

unpaid_composition <- summary_food_df_long_with_ghd %>%
  filter(type %in% c("hr_m", "hr_f"),
         country %in% cty_ghd,
         footprint_type %in% c("preparation_non.econ", "processing_non.econ", "growth_collection_non.econ")) %>%
  mutate(country = as.character(country)) %>%
  group_by(country, footprint_type) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = footprint_type, values_from = hr_per_cap_day, values_fill = 0)

cat("---- Unpaid time composition: preparation vs. processing vs. growth/collection ----\n")
cat(sprintf("Countries with any GHD unpaid data: n = %d\n", nrow(unpaid_composition)))
cat(sprintf("Mean hr/cap/day -- preparation: %.3f | processing: %.3f | growth_collection: %.3f\n",
            mean(unpaid_composition$preparation_non.econ, na.rm = TRUE),
            mean(unpaid_composition$processing_non.econ, na.rm = TRUE),
            mean(unpaid_composition$growth_collection_non.econ, na.rm = TRUE)))
cat(sprintf("Countries with growth_collection > 0: %d (%.0f%%) | processing > 0: %d (%.0f%%) | preparation > 0: %d (%.0f%%)\n",
            sum(unpaid_composition$growth_collection_non.econ > 0), 100 * mean(unpaid_composition$growth_collection_non.econ > 0),
            sum(unpaid_composition$processing_non.econ > 0), 100 * mean(unpaid_composition$processing_non.econ > 0),
            sum(unpaid_composition$preparation_non.econ > 0), 100 * mean(unpaid_composition$preparation_non.econ > 0)))
cat("\nIf growth_collection is zero (or near-zero) for most countries, CF_unpaid_produce\n")
cat("below is being driven by a small subsistence-farming subset -- worth knowing before\n")
cat("reading its result as general rather than as evidence about that subset specifically.\n\n")

#### 1. CF_unpaid_transform: preparation + processing, n=33 (EXIO-modeled, non-RoW) ####

cf_unpaid_transform_df <- summary_food_df_long_with_ghd %>%
  filter(type %in% c("hr_m", "hr_f"),
         country %in% cty_ghd,
         footprint_type %in% c("preparation_non.econ", "processing_non.econ"),
         !(as.character(country) %in% row_countries)) %>%
  mutate(country = as.character(country)) %>%
  group_by(country) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  left_join(pro_consumption, by = "country") %>%
  inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
  filter(hr_per_cap_day > 0, g_protein_per_cap_day > 0) %>%
  mutate(hr_per_50g_protein = hr_per_cap_day / g_protein_per_cap_day * 50,
         log_cf_transform = log(hr_per_50g_protein),
         log_protein      = log(g_protein_per_cap_day),
         log_gdp_pcap     = log(gdp_pcap_ppp))

cat(sprintf("CF_unpaid_transform (preparation + processing, unpaid, total consumption): n = %d\n", nrow(cf_unpaid_transform_df)))

pc_transform <- partial_cor(cf_unpaid_transform_df$log_cf_transform, cf_unpaid_transform_df$log_protein, cf_unpaid_transform_df$log_gdp_pcap)
cat(sprintf("Partial r | GDP per capita = %.3f (p = %.3g, n = %d)\n", pc_transform$r, pc_transform$p, pc_transform$n))

med_transform <- mediation_decomp(cf_unpaid_transform_df, "log_gdp_pcap", "log_cf_transform", "log_protein")
cat(sprintf("Mediation: total = %.3f | direct = %.3f | indirect = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g), n = %d\n",
            med_transform$total, med_transform$direct, med_transform$indirect,
            med_transform$prop_mediated * 100, med_transform$sobel_z, med_transform$sobel_p, med_transform$n))
cat("Path a for CF_unpaid_transform:\n")
print(round(summary(lm(log_cf_transform ~ log_gdp_pcap, data = cf_unpaid_transform_df))$coefficients, 4))

acme_transform_stat <- function(d) {
  fit_a <- lm(log_cf_transform ~ log_gdp_pcap, data = d)
  fit_b <- lm(log_protein ~ log_gdp_pcap + log_cf_transform, data = d)
  unname(coef(fit_a)["log_gdp_pcap"]) * unname(coef(fit_b)["log_cf_transform"])
}
point_acme_transform <- acme_transform_stat(cf_unpaid_transform_df)
boot_acme_transform  <- boot_stat(cf_unpaid_transform_df, acme_transform_stat, R)
report_boot(boot_acme_transform, point_acme_transform, "CF_unpaid_transform single-mediator ACME, n=33", sobel_p_for_comparison = med_transform$sobel_p)

#### 2. CF_unpaid_produce: growth/collection only, n=33 ####

cf_unpaid_produce_df <- summary_food_df_long_with_ghd %>%
  filter(type %in% c("hr_m", "hr_f"),
         country %in% cty_ghd,
         footprint_type == "growth_collection_non.econ",
         !(as.character(country) %in% row_countries)) %>%
  mutate(country = as.character(country)) %>%
  group_by(country) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  left_join(pro_consumption, by = "country") %>%
  inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
  filter(hr_per_cap_day > 0, g_protein_per_cap_day > 0) %>%
  mutate(hr_per_50g_protein = hr_per_cap_day / g_protein_per_cap_day * 50,
         log_cf_produce = log(hr_per_50g_protein),
         log_protein    = log(g_protein_per_cap_day),
         log_gdp_pcap   = log(gdp_pcap_ppp))

cat(sprintf("\nCF_unpaid_produce (growth/collection only, unpaid, total consumption): n = %d\n", nrow(cf_unpaid_produce_df)))
cat("(n will likely be much smaller than CF_unpaid_transform's -- countries with zero\n")
cat("growth_collection are dropped by the hr_per_cap_day > 0 filter, and subsistence\n")
cat("growing/foraging is far from universal in this sample.)\n")

pc_produce <- partial_cor(cf_unpaid_produce_df$log_cf_produce, cf_unpaid_produce_df$log_protein, cf_unpaid_produce_df$log_gdp_pcap)
cat(sprintf("Partial r | GDP per capita = %.3f (p = %.3g, n = %d)\n", pc_produce$r, pc_produce$p, pc_produce$n))

med_produce <- mediation_decomp(cf_unpaid_produce_df, "log_gdp_pcap", "log_cf_produce", "log_protein")
cat(sprintf("Mediation: total = %.3f | direct = %.3f | indirect = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g), n = %d\n",
            med_produce$total, med_produce$direct, med_produce$indirect,
            med_produce$prop_mediated * 100, med_produce$sobel_z, med_produce$sobel_p, med_produce$n))
cat("Path a for CF_unpaid_produce:\n")
print(round(summary(lm(log_cf_produce ~ log_gdp_pcap, data = cf_unpaid_produce_df))$coefficients, 4))

acme_produce_stat <- function(d) {
  fit_a <- lm(log_cf_produce ~ log_gdp_pcap, data = d)
  fit_b <- lm(log_protein ~ log_gdp_pcap + log_cf_produce, data = d)
  unname(coef(fit_a)["log_gdp_pcap"]) * unname(coef(fit_b)["log_cf_produce"])
}
point_acme_produce <- acme_produce_stat(cf_unpaid_produce_df)
boot_acme_produce  <- boot_stat(cf_unpaid_produce_df, acme_produce_stat, R)
report_boot(boot_acme_produce, point_acme_produce, "CF_unpaid_produce single-mediator ACME, n=33", sobel_p_for_comparison = med_produce$sobel_p)

#### 3. RoW-collapsed (n=38-ish) versions of both, same pattern as 9p/9s/9t ####
#### -- NOTE: whether GHD unpaid time can even be collapsed this way      ####
#### depends on whether RoW-flagged FABIO members have their OWN GHD      ####
#### survey data (unlike paid labor, unpaid time has no EXIO-region       ####
#### value to paste -- if a RoW member has no national time-use survey,   ####
#### it simply has no row here, full stop). Check the printed n before    ####
#### trusting this section; it may come back small or even empty,        ####
#### especially for CF_unpaid_produce given section 0's coverage check.  ####

row_protein_transform_9u <- summary_food_df_long_with_ghd %>%
  filter(type %in% c("hr_m", "hr_f"),
         country %in% cty_ghd,
         footprint_type %in% c("preparation_non.econ", "processing_non.econ"),
         as.character(country) %in% row_countries) %>%
  mutate(country = as.character(country)) %>%
  group_by(country) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  left_join(pro_consumption, by = "country") %>%
  inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
  filter(hr_per_cap_day > 0, g_protein_per_cap_day > 0) %>%
  left_join(row_lookup, by = "country") %>%
  left_join(pop_data_yr, by = c("country" = "iso3c")) %>%
  filter(!is.na(population), population > 0) %>%
  group_by(exio_region) %>%
  summarise(n_members_used = n(),
            pop_used = sum(population),
            hr_per_cap_day_region = weighted.mean(hr_per_cap_day, population),
            g_protein_per_cap_day = weighted.mean(g_protein_per_cap_day, population),
            gdp_pcap_ppp = weighted.mean(gdp_pcap_ppp, population, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(hr_per_50g_protein = hr_per_cap_day_region / g_protein_per_cap_day * 50) %>%
  left_join(row_full_pop, by = "exio_region") %>%
  mutate(pop_coverage_pct = pop_used / pop_total_region * 100)

cat(sprintf("\n---- RoW-aggregate collapse (CF_unpaid_transform): %d of 5 regions have any GHD-covered members ----\n",
            nrow(row_protein_transform_9u)))
if (nrow(row_protein_transform_9u) > 0) {
  print(row_protein_transform_9u %>% select(exio_region, n_members_used, pop_coverage_pct) %>%
          mutate(pop_coverage_pct = round(pop_coverage_pct, 1)))
}

row_protein_produce_9u <- summary_food_df_long_with_ghd %>%
  filter(type %in% c("hr_m", "hr_f"),
         country %in% cty_ghd,
         footprint_type == "growth_collection_non.econ",
         as.character(country) %in% row_countries) %>%
  mutate(country = as.character(country)) %>%
  group_by(country) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  left_join(pro_consumption, by = "country") %>%
  inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
  filter(hr_per_cap_day > 0, g_protein_per_cap_day > 0) %>%
  left_join(row_lookup, by = "country") %>%
  left_join(pop_data_yr, by = c("country" = "iso3c")) %>%
  filter(!is.na(population), population > 0) %>%
  group_by(exio_region) %>%
  summarise(n_members_used = n(),
            pop_used = sum(population),
            hr_per_cap_day_region = weighted.mean(hr_per_cap_day, population),
            g_protein_per_cap_day = weighted.mean(g_protein_per_cap_day, population),
            gdp_pcap_ppp = weighted.mean(gdp_pcap_ppp, population, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(hr_per_50g_protein = hr_per_cap_day_region / g_protein_per_cap_day * 50) %>%
  left_join(row_full_pop, by = "exio_region") %>%
  mutate(pop_coverage_pct = pop_used / pop_total_region * 100)

cat(sprintf("\n---- RoW-aggregate collapse (CF_unpaid_produce): %d of 5 regions have any GHD-covered members ----\n",
            nrow(row_protein_produce_9u)))
if (nrow(row_protein_produce_9u) > 0) {
  print(row_protein_produce_9u %>% select(exio_region, n_members_used, pop_coverage_pct) %>%
          mutate(pop_coverage_pct = round(pop_coverage_pct, 1)))
}

if (nrow(row_protein_transform_9u) > 0) {
  cf_unpaid_transform_df38 <- bind_rows(
    cf_unpaid_transform_df %>% select(country, hr_per_50g_protein, g_protein_per_cap_day, gdp_pcap_ppp),
    row_protein_transform_9u %>% select(country = exio_region, hr_per_50g_protein, g_protein_per_cap_day, gdp_pcap_ppp)
  ) %>%
    mutate(log_cf_transform = log(hr_per_50g_protein), log_protein = log(g_protein_per_cap_day), log_gdp_pcap = log(gdp_pcap_ppp))

  cat(sprintf("\n[RoW-collapsed] CF_unpaid_transform: %d EXIO-modeled + %d RoW-aggregate = %d rows\n",
              nrow(cf_unpaid_transform_df), nrow(row_protein_transform_9u), nrow(cf_unpaid_transform_df38)))
  med_transform38 <- mediation_decomp(cf_unpaid_transform_df38, "log_gdp_pcap", "log_cf_transform", "log_protein")
  cat(sprintf("[RoW-collapsed] Mediation: total = %.3f | direct = %.3f | indirect = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g), n = %d\n",
              med_transform38$total, med_transform38$direct, med_transform38$indirect,
              med_transform38$prop_mediated * 100, med_transform38$sobel_z, med_transform38$sobel_p, med_transform38$n))
  point_acme_transform38 <- acme_transform_stat(cf_unpaid_transform_df38)
  boot_acme_transform38  <- boot_stat(cf_unpaid_transform_df38, acme_transform_stat, R)
  report_boot(boot_acme_transform38, point_acme_transform38, "CF_unpaid_transform single-mediator ACME, RoW-collapsed", sobel_p_for_comparison = med_transform38$sobel_p)
} else {
  cat("\n[RoW-collapsed] Skipped for CF_unpaid_transform -- no RoW-flagged FABIO members have their own GHD data.\n")
}

if (nrow(row_protein_produce_9u) > 0) {
  cf_unpaid_produce_df38 <- bind_rows(
    cf_unpaid_produce_df %>% select(country, hr_per_50g_protein, g_protein_per_cap_day, gdp_pcap_ppp),
    row_protein_produce_9u %>% select(country = exio_region, hr_per_50g_protein, g_protein_per_cap_day, gdp_pcap_ppp)
  ) %>%
    mutate(log_cf_produce = log(hr_per_50g_protein), log_protein = log(g_protein_per_cap_day), log_gdp_pcap = log(gdp_pcap_ppp))

  cat(sprintf("\n[RoW-collapsed] CF_unpaid_produce: %d EXIO-modeled + %d RoW-aggregate = %d rows\n",
              nrow(cf_unpaid_produce_df), nrow(row_protein_produce_9u), nrow(cf_unpaid_produce_df38)))
  med_produce38 <- mediation_decomp(cf_unpaid_produce_df38, "log_gdp_pcap", "log_cf_produce", "log_protein")
  cat(sprintf("[RoW-collapsed] Mediation: total = %.3f | direct = %.3f | indirect = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g), n = %d\n",
              med_produce38$total, med_produce38$direct, med_produce38$indirect,
              med_produce38$prop_mediated * 100, med_produce38$sobel_z, med_produce38$sobel_p, med_produce38$n))
  point_acme_produce38 <- acme_produce_stat(cf_unpaid_produce_df38)
  boot_acme_produce38  <- boot_stat(cf_unpaid_produce_df38, acme_produce_stat, R)
  report_boot(boot_acme_produce38, point_acme_produce38, "CF_unpaid_produce single-mediator ACME, RoW-collapsed", sobel_p_for_comparison = med_produce38$sobel_p)
} else {
  cat("\n[RoW-collapsed] Skipped for CF_unpaid_produce -- no RoW-flagged FABIO members have their own GHD data.\n")
}

cat("\n#### What to look at ####\n")
cat("1. The magnitude check (section 0): if growth_collection is near-zero for\n")
cat("   most countries, CF_unpaid_produce's result reflects a small subsistence-\n")
cat("   farming subset, not a general 'production-adjacent unpaid time' story --\n")
cat("   say so explicitly if so, and treat n and its bootstrap CI with extra\n")
cat("   caution.\n")
cat("2. If CF_unpaid_transform comes out null (like CF_time) and CF_unpaid_produce\n")
cat("   comes out significant (like CF_paid), that confirms the production-vs-\n")
cat("   transformation mechanism directly, and the draft's 'paid vs. unpaid'\n")
cat("   framing should be revised to 'production-adjacent vs. meal-\n")
cat("   transformation activity,' with paid/unpaid as a correlate, not the\n")
cat("   operative distinction.\n")
cat("3. If both come out null, or both come out significant, the mechanism\n")
cat("   proposed in the draft does not hold up and that paragraph needs to be\n")
cat("   walked back, not just caveated further.\n")
