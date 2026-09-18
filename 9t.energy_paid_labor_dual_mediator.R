#### CF_paid vs. CF_energy: the genuinely matched pair ####
#
# "Is Energy CF the Missing Channel?" (9m/9o) paired CF_energy (paid/
# market-economic energy only -- EXIOBASE cannot resolve household
# energy use, Limitations) against the all-work CF_time (paid + unpaid,
# Table 1). That comparison is not apples-to-apples on the labor side:
# CF_energy is structurally paid-only, CF_time is not. This script builds
# CF_paid -- paid labor only, domestic- AND import-origin, over the SAME
# total-consumption denominator as CF_energy -- so the pair differs only
# in labor vs. energy, nothing else. Table 1's proposed new row.
#
# tradeoff_protein_econlabor_consump mirrors tradeoff_protein_allwork_
# consump's construction exactly (same footprint_type summing, same
# total-consumption denominator via pro_consumption) but reads from
# summary_food_df_long (base EXIOBASE economic-labor data) instead of
# summary_food_df_long_with_ghd (the GHD-augmented, all-work version) --
# so it is paid-only by construction, not by a separate filter step.
#
# Run AFTER 2.analyze_result.R (builds tradeoff_protein_econlabor_consump,
# tradeoff_protein_allwork_consump, tradeoff_protein_allwork_energy or
# equivalent CF_energy source), 9.capability_set_income_control.R
# (wdi_latest, mediation_decomp(), partial_cor(), df -- reused below as
# the standard n=33 country list to restrict CF_paid's n=43 sample
# against -- and row_lookup, row_full_pop, pop_data_yr, from that
# script's own RoW-collapse block, needed for the n=38 version further
# down), 9o.bootstrap_mediation_and_interaction.R (boot_stat(),
# report_boot(), R), and 9p.energy_cf_row_collapsed.R (df2_collapsed --
# needed for section 4, the n=38 CF_paid+CF_energy pairing)
# in the same session.

library(tidyverse)

#### 1. CF_paid, single-mediator: does it reproduce the all-work story? ####

cf_paid_df <- tradeoff_protein_econlabor_consump %>%
  filter(!is_row) %>%
  group_by(country) %>%
  summarise(hr_per_50g_protein = sum(hr_per_50g_protein, na.rm = TRUE),
            g_protein_per_cap_day = first(g_protein_per_cap_day),
            .groups = "drop") %>%
  inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
  filter(hr_per_50g_protein > 0, g_protein_per_cap_day > 0) %>%
  mutate(log_cf_paid   = log(hr_per_50g_protein),
         log_protein   = log(g_protein_per_cap_day),
         log_gdp_pcap  = log(gdp_pcap_ppp))

cat(sprintf("CF_paid (paid labor, all-origin, total consumption): n = %d\n", nrow(cf_paid_df)))

pc_paid <- partial_cor(cf_paid_df$log_cf_paid, cf_paid_df$log_protein, cf_paid_df$log_gdp_pcap)
cat(sprintf("Partial r | GDP per capita = %.3f (p = %.3g, n = %d)\n", pc_paid$r, pc_paid$p, pc_paid$n))

med_paid <- mediation_decomp(cf_paid_df, "log_gdp_pcap", "log_cf_paid", "log_protein")
cat(sprintf("Mediation: total = %.3f | direct = %.3f | indirect = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g), n = %d\n",
            med_paid$total, med_paid$direct, med_paid$indirect,
            med_paid$prop_mediated * 100, med_paid$sobel_z, med_paid$sobel_p, med_paid$n))
cat("\nPath a for CF_paid:\n")
print(round(summary(lm(log_cf_paid ~ log_gdp_pcap, data = cf_paid_df))$coefficients, 4))

cat("\nFor comparison: all-work CF path a = -0.52; CF_paid,dom (SI-2) path a = -1.16.\n")
cat("CF_paid sits between the two on origin (all-origin, like all-work) and labor\n")
cat("type (paid-only, like CF_paid,dom) -- worth seeing where its own elasticity lands.\n")

#### 1b. Same test, restricted to the same n=33 EXIO-only sample as ####
#### everything else in this paper -- is the n=43 significance real, ####
#### or just the extra statistical power from 10 more countries?      ####

cf_paid_df33 <- cf_paid_df %>% filter(country %in% df$country)
cat(sprintf("\nCF_paid restricted to the standard n=33 sample: n = %d\n", nrow(cf_paid_df33)))

pc_paid33 <- partial_cor(cf_paid_df33$log_cf_paid, cf_paid_df33$log_protein, cf_paid_df33$log_gdp_pcap)
cat(sprintf("Partial r | GDP per capita = %.3f (p = %.3g, n = %d)\n", pc_paid33$r, pc_paid33$p, pc_paid33$n))

med_paid33 <- mediation_decomp(cf_paid_df33, "log_gdp_pcap", "log_cf_paid", "log_protein")
cat(sprintf("Mediation: total = %.3f | direct = %.3f | indirect = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g), n = %d\n",
            med_paid33$total, med_paid33$direct, med_paid33$indirect,
            med_paid33$prop_mediated * 100, med_paid33$sobel_z, med_paid33$sobel_p, med_paid33$n))
cat("\nPath a for CF_paid, n=33:\n")
print(round(summary(lm(log_cf_paid ~ log_gdp_pcap, data = cf_paid_df33))$coefficients, 4))

cat("\n#### What to look at ####\n")
cat("If partial r / mediation p stay roughly this significant at n=33, the n=43\n")
cat("result was real, not a sample-size artifact -- worth a bootstrap check next.\n")
cat("If they go null (like all-work, CF_paid,dom, and CF_dom-effort all did at\n")
cat("n=33), the n=43 significance was just extra power from the larger sample.\n")

#### 1c. Bootstrap the single-mediator (non-paired) CF_paid ACME at n=33 -- ####
#### the Sobel p (0.049) is right on the boundary, exactly the small-n     ####
#### case this session has repeatedly found Sobel unreliable for.         ####

acme_paid_single_stat <- function(d) {
  fit_a <- lm(log_cf_paid ~ log_gdp_pcap, data = d)
  fit_b <- lm(log_protein ~ log_gdp_pcap + log_cf_paid, data = d)
  unname(coef(fit_a)["log_gdp_pcap"]) * unname(coef(fit_b)["log_cf_paid"])
}
point_acme_paid_single <- acme_paid_single_stat(cf_paid_df33)
boot_acme_paid_single <- boot_stat(cf_paid_df33, acme_paid_single_stat, R)
report_boot(boot_acme_paid_single, point_acme_paid_single,
            "CF_paid single-mediator ACME, n=33 (not paired with CF_energy)",
            sobel_p_for_comparison = med_paid33$sobel_p)

#### 1d. RoW-collapsed (n=38) version, matching the standard "both       ####
#### country samples" convention every other CF variant here follows.   ####
#### Built from the n=33-restricted EXIO sample + 5 population-weighted ####
#### RoW regions (33+5=38), NOT the full n=43+RoW range, so it stays    ####
#### comparable to every other n=38 result in this paper. Same pattern  ####
#### as 9p/9s's RoW-collapse blocks, applied to CF_paid's source.       ####

row_protein_paid_9t <- tradeoff_protein_econlabor_consump %>%
  filter(is_row) %>%
  group_by(country) %>%
  summarise(hr_per_cap_day = sum(hr_per_cap_day, na.rm = TRUE),
            g_protein_per_cap_day = first(g_protein_per_cap_day),
            .groups = "drop") %>%
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

cat("\n---- RoW-aggregate collapse (CF_paid): population coverage of the subset used ----\n")
print(row_protein_paid_9t %>% select(exio_region, n_members_used, pop_coverage_pct) %>%
        mutate(pop_coverage_pct = round(pop_coverage_pct, 1)))

cf_paid_df38 <- bind_rows(
  cf_paid_df33 %>% select(country, hr_per_50g_protein, g_protein_per_cap_day, gdp_pcap_ppp),
  row_protein_paid_9t %>% select(country = exio_region, hr_per_50g_protein, g_protein_per_cap_day, gdp_pcap_ppp)
) %>%
  mutate(log_cf_paid  = log(hr_per_50g_protein),
         log_protein  = log(g_protein_per_cap_day),
         log_gdp_pcap = log(gdp_pcap_ppp))

cat(sprintf("\nCollapsed CF_paid sample: %d EXIO-modeled + %d RoW-aggregate = %d rows\n",
            nrow(cf_paid_df33), nrow(row_protein_paid_9t), nrow(cf_paid_df38)))

pc_paid38 <- partial_cor(cf_paid_df38$log_cf_paid, cf_paid_df38$log_protein, cf_paid_df38$log_gdp_pcap)
cat(sprintf("[RoW-collapsed] Partial r | GDP per capita = %.3f (p = %.3g, n = %d)\n",
            pc_paid38$r, pc_paid38$p, pc_paid38$n))

med_paid38 <- mediation_decomp(cf_paid_df38, "log_gdp_pcap", "log_cf_paid", "log_protein")
cat(sprintf("[RoW-collapsed] Mediation: total = %.3f | direct = %.3f | indirect = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g), n = %d\n",
            med_paid38$total, med_paid38$direct, med_paid38$indirect,
            med_paid38$prop_mediated * 100, med_paid38$sobel_z, med_paid38$sobel_p, med_paid38$n))

cat("\n[RoW-collapsed] Path a for CF_paid:\n")
print(round(summary(lm(log_cf_paid ~ log_gdp_pcap, data = cf_paid_df38))$coefficients, 4))

point_acme_paid38 <- acme_paid_single_stat(cf_paid_df38)
boot_acme_paid38 <- boot_stat(cf_paid_df38, acme_paid_single_stat, R)
report_boot(boot_acme_paid38, point_acme_paid38,
            "CF_paid single-mediator ACME, n=38 RoW-collapsed",
            sobel_p_for_comparison = med_paid38$sobel_p)

#### 2. CF_paid vs. CF_energy: correlation, then the parallel two-mediator ####
#### model, replacing CF_time with CF_paid throughout. ####

# NOTE: adjust the object/column names in this block to match whatever
# CF_energy actually reads from in 9m (cf_df_energy / df2) -- reusing that
# same energy source here, joined to cf_paid_df instead of the all-work df.
cf_paid_energy_df <- cf_paid_df %>%
  inner_join(df2 %>% select(country, log_cf_energy), by = "country")

cat(sprintf("\nCF_paid + CF_energy, joined: n = %d\n", nrow(cf_paid_energy_df)))

cor_paid_energy <- cor.test(cf_paid_energy_df$log_cf_paid, cf_paid_energy_df$log_cf_energy)
cat(sprintf("Correlation, CF_paid vs. CF_energy: r = %.3f (p = %.3g, n = %d)\n",
            cor_paid_energy$estimate, cor_paid_energy$p.value, nrow(cf_paid_energy_df)))
cat("(compare against CF_time vs. CF_energy: r = -0.53, p = 0.002, n = 33)\n")

parallel_mediation_decomp_paid <- function(df, x_col, m1_col, m2_col, y_col) {
  keep <- complete.cases(df[c(x_col, m1_col, m2_col, y_col)])
  x <- df[[x_col]][keep]; m1 <- df[[m1_col]][keep]; m2 <- df[[m2_col]][keep]; y <- df[[y_col]][keep]
  fit_total <- lm(y ~ x)
  fit_med1  <- lm(m1 ~ x)
  fit_med2  <- lm(m2 ~ x)
  fit_out   <- lm(y ~ x + m1 + m2)
  a1 <- unname(coef(fit_med1)["x"]); se_a1 <- summary(fit_med1)$coefficients["x", "Std. Error"]
  a2 <- unname(coef(fit_med2)["x"]); se_a2 <- summary(fit_med2)$coefficients["x", "Std. Error"]
  b1 <- unname(coef(fit_out)["m1"]); se_b1 <- summary(fit_out)$coefficients["m1", "Std. Error"]
  b2 <- unname(coef(fit_out)["m2"]); se_b2 <- summary(fit_out)$coefficients["m2", "Std. Error"]
  total <- unname(coef(fit_total)["x"]); direct <- unname(coef(fit_out)["x"])
  ind1 <- a1 * b1; se_ind1 <- sqrt(b1^2*se_a1^2 + a1^2*se_b1^2); z1 <- ind1/se_ind1
  ind2 <- a2 * b2; se_ind2 <- sqrt(b2^2*se_a2^2 + a2^2*se_b2^2); z2 <- ind2/se_ind2
  list(n = length(x), total = total, direct = direct,
       indirect_paid = ind1, prop_mediated_paid = ind1/total, sobel_z_paid = z1, sobel_p_paid = 2*(1-pnorm(abs(z1))),
       indirect_energy = ind2, prop_mediated_energy = ind2/total, sobel_z_energy = z2, sobel_p_energy = 2*(1-pnorm(abs(z2))))
}

med2_paid <- parallel_mediation_decomp_paid(cf_paid_energy_df, "log_gdp_pcap", "log_cf_paid", "log_cf_energy", "log_protein")
cat(sprintf("\nParallel dual-mediator (CF_paid + CF_energy): total = %.3f | direct = %.3f\n", med2_paid$total, med2_paid$direct))
cat(sprintf("  indirect via CF_paid   = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g)\n",
            med2_paid$indirect_paid, med2_paid$prop_mediated_paid*100, med2_paid$sobel_z_paid, med2_paid$sobel_p_paid))
cat(sprintf("  indirect via CF_energy = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g)\n",
            med2_paid$indirect_energy, med2_paid$prop_mediated_energy*100, med2_paid$sobel_z_energy, med2_paid$sobel_p_energy))

#### 3. Bootstrap both indirect effects (reuses 9o's boot_stat/report_boot) ####

acme_paid_stat <- function(d) {
  fit_a <- lm(log_cf_paid ~ log_gdp_pcap, data = d)
  fit_b <- lm(log_protein ~ log_gdp_pcap + log_cf_paid + log_cf_energy, data = d)
  unname(coef(fit_a)["log_gdp_pcap"]) * unname(coef(fit_b)["log_cf_paid"])
}
acme_energy_stat_paidpair <- function(d) {
  fit_a <- lm(log_cf_energy ~ log_gdp_pcap, data = d)
  fit_b <- lm(log_protein ~ log_gdp_pcap + log_cf_paid + log_cf_energy, data = d)
  unname(coef(fit_a)["log_gdp_pcap"]) * unname(coef(fit_b)["log_cf_energy"])
}
point_acme_paid   <- acme_paid_stat(cf_paid_energy_df)
point_acme_energy2 <- acme_energy_stat_paidpair(cf_paid_energy_df)
boot_acme_paid    <- boot_stat(cf_paid_energy_df, acme_paid_stat, R)
boot_acme_energy2 <- boot_stat(cf_paid_energy_df, acme_energy_stat_paidpair, R)

report_boot(boot_acme_paid, point_acme_paid, "ACME via CF_paid (paid-labor pair)", sobel_p_for_comparison = med2_paid$sobel_p_paid)
report_boot(boot_acme_energy2, point_acme_energy2, "ACME via CF_energy (paid-labor pair)", sobel_p_for_comparison = med2_paid$sobel_p_energy)

cat("\n#### What to look at ####\n")
cat("Compare against the all-work pair (9m/9o): CF_time ACME was null (Sobel\n")
cat("p=0.42, bootstrap CI included zero); CF_energy ACME was the strongest result\n")
cat("in the paper (Sobel p=0.124, bootstrap p=0.018, CI excluded zero). If CF_paid's\n")
cat("ACME here comes out differently from CF_time's -- e.g. reaching significance\n")
cat("where CF_time didn't -- that's evidence the energy relationship is specifically\n")
cat("about paid labor, not unpaid, resolving the question the Limitations-note\n")
cat("caveat currently leaves open. If CF_paid looks the same as CF_time, the paid/\n")
cat("unpaid split isn't where the action is.\n")

#### 4. RoW-collapsed (n=38) version of the CF_paid + CF_energy pair, ####
#### matching this paper's standard "both country samples" convention  ####
#### (as extended for the all-work pair in 9p section 4). Uses         ####
#### cf_paid_df38 (section 1d, above) joined against df2_collapsed     ####
#### (9p) instead of cf_paid_df/df2 -- same 5 RoW regions on both      ####
#### sides, built from the same row_lookup/pop_data_yr machinery, so   ####
#### the country keys line up directly.                                ####
#### Requires 9p run in this session (df2_collapsed).                  ####

cf_paid_energy_df38 <- cf_paid_df38 %>%
  inner_join(df2_collapsed %>% select(country, log_cf_energy), by = "country")

cat(sprintf("\n[RoW-collapsed] CF_paid + CF_energy, joined: n = %d\n", nrow(cf_paid_energy_df38)))

cor_paid_energy38 <- cor.test(cf_paid_energy_df38$log_cf_paid, cf_paid_energy_df38$log_cf_energy)
cat(sprintf("[RoW-collapsed] Correlation, CF_paid vs. CF_energy: r = %.3f (p = %.3g, n = %d)\n",
            cor_paid_energy38$estimate, cor_paid_energy38$p.value, nrow(cf_paid_energy_df38)))

med2_paid38 <- parallel_mediation_decomp_paid(cf_paid_energy_df38, "log_gdp_pcap", "log_cf_paid", "log_cf_energy", "log_protein")
cat(sprintf("\n[RoW-collapsed] Parallel dual-mediator (CF_paid + CF_energy): total = %.3f | direct = %.3f\n",
            med2_paid38$total, med2_paid38$direct))
cat(sprintf("  indirect via CF_paid   = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g)\n",
            med2_paid38$indirect_paid, med2_paid38$prop_mediated_paid*100, med2_paid38$sobel_z_paid, med2_paid38$sobel_p_paid))
cat(sprintf("  indirect via CF_energy = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g)\n",
            med2_paid38$indirect_energy, med2_paid38$prop_mediated_energy*100, med2_paid38$sobel_z_energy, med2_paid38$sobel_p_energy))

acme_paid_stat38 <- function(d) {
  fit_a <- lm(log_cf_paid ~ log_gdp_pcap, data = d)
  fit_b <- lm(log_protein ~ log_gdp_pcap + log_cf_paid + log_cf_energy, data = d)
  unname(coef(fit_a)["log_gdp_pcap"]) * unname(coef(fit_b)["log_cf_paid"])
}
acme_energy_stat_paidpair38 <- function(d) {
  fit_a <- lm(log_cf_energy ~ log_gdp_pcap, data = d)
  fit_b <- lm(log_protein ~ log_gdp_pcap + log_cf_paid + log_cf_energy, data = d)
  unname(coef(fit_a)["log_gdp_pcap"]) * unname(coef(fit_b)["log_cf_energy"])
}
point_acme_paid38    <- acme_paid_stat38(cf_paid_energy_df38)
point_acme_energy238 <- acme_energy_stat_paidpair38(cf_paid_energy_df38)
boot_acme_paid38b    <- boot_stat(cf_paid_energy_df38, acme_paid_stat38, R)
boot_acme_energy238  <- boot_stat(cf_paid_energy_df38, acme_energy_stat_paidpair38, R)

report_boot(boot_acme_paid38b, point_acme_paid38,
            sprintf("ACME via CF_paid (paid-labor pair), RoW-collapsed (n=%d)", nrow(cf_paid_energy_df38)),
            sobel_p_for_comparison = med2_paid38$sobel_p_paid)
report_boot(boot_acme_energy238, point_acme_energy238,
            sprintf("ACME via CF_energy (paid-labor pair), RoW-collapsed (n=%d)", nrow(cf_paid_energy_df38)),
            sobel_p_for_comparison = med2_paid38$sobel_p_energy)

cat("\n#### What to look at ####\n")
cat("Compare against section 3's n=33 result and against 9p section 4's n=38\n")
cat("all-work pair. Four cells total now exist for CF_energy's ACME: {n=33,\n")
cat("n=38} x {paired with CF_time, paired with CF_paid} -- if CF_energy stays\n")
cat("significant (bootstrap CI excludes zero) in all four, that is about as\n")
cat("robust as a single mediation finding gets in this dataset.\n")
