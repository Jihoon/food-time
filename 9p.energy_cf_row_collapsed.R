#### Does adding RoW (low-income) regions change the energy-CF interaction ####
#### and its bootstrap check? ####
#
# Motivation: the GDP x CF_energy interaction (9n) is identified mostly by
# the few low-income leverage points (India, Indonesia, RoW Africa -- see
# Findings/Discussion), and EXIO-only (n=33) under-represents exactly that
# part of the income distribution. A case-resampling bootstrap (9o) can
# read that sparsity as instability (wide CI) even when the underlying
# relationship is real -- more low-income data should stabilize both the
# original estimate AND the bootstrap CI, not just one or the other.
#
# This builds the RoW-collapsed (n=38) version of CF_energy, mirroring
# script 9's row_collapsed block exactly (same population-weighting, same
# row_lookup/pop_data_yr/row_full_pop), just adding mj_per_cap_day
# alongside hr_per_cap_day -- both come from the same source rows in
# tradeoff_protein_allwork_consump, so this is one extension, not two.
#
# Run AFTER 2.analyze_result.R, 9.capability_set_income_control.R (reuses
# row_lookup, row_full_pop, wdi_latest, pop_data_yr, df_collapsed),
# 9m.energy_time_dual_mediator.R (reuses df2), 9o's boot_stat()/report_boot()
# helper functions (source 9o first, or paste those two functions in).

#### 1. RoW-collapsed energy CF, mirroring script 9's row_collapsed exactly ####

row_collapsed_energy <- tradeoff_protein_allwork_consump %>%
  filter(is_row) %>%
  group_by(country) %>%
  summarise(hr_per_cap_day = sum(hr_per_cap_day, na.rm = TRUE),
            mj_per_cap_day = sum(mj_per_cap_day, na.rm = TRUE),
            g_protein_per_cap_day = first(g_protein_per_cap_day),
            .groups = "drop") %>%
  inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
  filter(hr_per_cap_day > 0, mj_per_cap_day > 0, g_protein_per_cap_day > 0) %>%
  left_join(row_lookup, by = "country") %>%
  left_join(pop_data_yr, by = c("country" = "iso3c")) %>%
  filter(!is.na(population), population > 0) %>%
  group_by(exio_region) %>%
  summarise(n_members_used = n(),
            pop_used = sum(population),
            hr_per_cap_day_region = weighted.mean(hr_per_cap_day, population),
            mj_per_cap_day_region = weighted.mean(mj_per_cap_day, population),
            g_protein_per_cap_day = weighted.mean(g_protein_per_cap_day, population),
            gdp_pcap_ppp          = weighted.mean(gdp_pcap_ppp, population, na.rm = TRUE),
            gdp_per_worker        = weighted.mean(gdp_per_worker, population, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(hr_per_50g_protein = hr_per_cap_day_region / g_protein_per_cap_day * 50,
         mj_per_50g_protein = mj_per_cap_day_region / g_protein_per_cap_day * 50) %>%
  left_join(row_full_pop, by = "exio_region") %>%
  mutate(pop_coverage_pct = pop_used / pop_total_region * 100)

cat("\n---- RoW-aggregate collapse: population coverage of the subset used ----\n")
print(row_collapsed_energy %>% select(exio_region, n_members_used, pop_coverage_pct) %>%
        mutate(pop_coverage_pct = round(pop_coverage_pct, 1)))

df2_collapsed <- bind_rows(
  df2 %>% select(country, hr_per_50g_protein, mj_per_50g_protein, g_protein_per_cap_day, gdp_pcap_ppp, gdp_per_worker),
  row_collapsed_energy %>% select(country = exio_region, hr_per_50g_protein, mj_per_50g_protein, g_protein_per_cap_day, gdp_pcap_ppp, gdp_per_worker)
) %>%
  mutate(log_cf         = log(hr_per_50g_protein),
         log_cf_energy  = log(mj_per_50g_protein),
         log_protein    = log(g_protein_per_cap_day),
         log_gdp_pcap   = log(gdp_pcap_ppp),
         log_gdp_worker = log(gdp_per_worker))

cat(sprintf("\nCollapsed energy-CF sample: %d EXIO-modeled + %d RoW-aggregate points = %d rows\n",
            nrow(df2), nrow(row_collapsed_energy), nrow(df2_collapsed)))

#### 1b. CF_time vs. CF_energy correlation, n=38 -- the draft currently  ####
#### only reports this at n=33 (9m: r=-0.53, p=0.002); this fills in the ####
#### RoW-collapsed counterpart so that number can be reported on its own ####
#### rather than as a combined range.                                   ####

cor_cf_time_energy_c <- cor.test(df2_collapsed$log_cf, df2_collapsed$log_cf_energy)
cat(sprintf("\n[RoW-collapsed] Correlation between CF_time and CF_energy (log-log): r = %.3f (p = %.3g, n = %d)\n",
            cor_cf_time_energy_c$estimate, cor_cf_time_energy_c$p.value, nrow(df2_collapsed)))
cat("(compare against EXIO-only, n=33: r = -0.53, p = 0.002)\n")

#### 2. GDP x CF_energy interaction, n=38 (compare against 9n's n=33 result) ####

fit_int_energy_c <- lm(log_protein ~ log_gdp_pcap * log_cf_energy, data = df2_collapsed)
print(round(summary(fit_int_energy_c)$coefficients, 4))
sobel_p_int_energy_c <- summary(fit_int_energy_c)$coefficients["log_gdp_pcap:log_cf_energy", "Pr(>|t|)"]
cat(sprintf("\nGDP x CF_energy interaction [RoW-collapsed, n=%d]: estimate = %.4f, p = %.4g\n",
            nrow(df2_collapsed),
            coef(fit_int_energy_c)["log_gdp_pcap:log_cf_energy"], sobel_p_int_energy_c))

#### 3. Bootstrap the same interaction on the n=38 sample (needs 9o's helpers) ####

interaction_energy_stat_c <- function(d) {
  fit <- lm(log_protein ~ log_gdp_pcap * log_cf_energy, data = d)
  unname(coef(fit)["log_gdp_pcap:log_cf_energy"])
}
point_int_energy_c <- interaction_energy_stat_c(df2_collapsed)
boot_int_energy_c  <- boot_stat(df2_collapsed, interaction_energy_stat_c, R)
report_boot(boot_int_energy_c, point_int_energy_c,
            sprintf("[5] GDP x CF_energy interaction, RoW-collapsed (n=%d) -- does more low-income data stabilize this?", nrow(df2_collapsed)),
            sobel_p_for_comparison = sobel_p_int_energy_c)

cat("\n#### What to look at ####\n")
cat("Compare [5]'s CI width against 9o's [4] (n=33, CI was [-0.31, 0.12]).\n")
cat("Narrower here would confirm the n=33 CI was inflated by leverage-point\n")
cat("sparsity, not just noise -- and if [5]'s CI excludes zero where [4]'s\n")
cat("didn't, that's a real, reportable update to the draft, not a reversal\n")
cat("of the bootstrap finding, an extension of it with better data.\n")

#### 3b. Marginal effect of CF_energy on protein at representative income ####
#### levels ($2,000 and $49,000/capita), same calculation as 9d's for      ####
#### CF_time -- the draft currently only states the CF_energy interaction  ####
#### is "positive at low income, negative at high income" qualitatively,   ####
#### without the matching numbers CF_time got. Both samples, for parity.   ####

marginal_effect_at_energy <- function(data, label, income_levels = c(2000, 49000)) {
  fit <- lm(log_protein ~ log_gdp_pcap * log_cf_energy, data = data)
  g2 <- unname(coef(fit)["log_cf_energy"])
  g3 <- unname(coef(fit)["log_gdp_pcap:log_cf_energy"])
  V  <- vcov(fit)
  cat(sprintf("\n==== %s (n = %d): marginal effect of CF_energy on protein, by income level ====\n", label, nrow(data)))
  for (x in income_levels) {
    lx <- log(x)
    me <- g2 + g3 * lx
    se <- sqrt(V["log_cf_energy", "log_cf_energy"] + lx^2 * V["log_gdp_pcap:log_cf_energy", "log_gdp_pcap:log_cf_energy"] +
                 2 * lx * V["log_cf_energy", "log_gdp_pcap:log_cf_energy"])
    cat(sprintf("  at $%s/capita: marginal effect = %.4f (SE %.4f)\n", format(x, big.mark = ","), me, se))
  }
  invisible(fit)
}

marginal_effect_at_energy(df2, "EXIO-only")
marginal_effect_at_energy(df2_collapsed, "RoW-collapsed")

#### 4. RoW-collapsed (n=38) version of 9m's PARALLEL dual-mediator ACME ####
#### decomposition (GDP -> {CF_time, CF_energy} -> protein), the source  ####
#### of "Is Energy CF the Missing Channel?"'s headline result. 9m only   ####
#### ever ran this at n=33; df2_collapsed above already has both        ####
#### hr_per_50g_protein (CF_time) and mj_per_50g_protein (CF_energy), so ####
#### no new data build is needed -- just re-run 9m's                    ####
#### parallel_mediation_decomp() on df2_collapsed instead of df2, then   ####
#### bootstrap both indirect effects the same way 9o did at n=33.        ####
#### Requires 9m's parallel_mediation_decomp() sourced in this session.  ####

med2_gdp_pcap_c <- parallel_mediation_decomp(df2_collapsed, "log_gdp_pcap", "log_cf", "log_cf_energy", "log_protein")

cat(sprintf("\n---- Parallel two-mediator decomposition, RoW-collapsed (n=%d): GDP -> {CF_time, CF_energy} -> protein ----\n",
            med2_gdp_pcap_c$n))
cat(sprintf("total = %.3f | direct (net of both mediators) = %.3f\n", med2_gdp_pcap_c$total, med2_gdp_pcap_c$direct))
cat(sprintf("indirect via CF_time   = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g) | b_time (net of CF_energy)   = %.3f (SE %.3f)\n",
            med2_gdp_pcap_c$indirect_time, med2_gdp_pcap_c$prop_mediated_time * 100,
            med2_gdp_pcap_c$sobel_z_time, med2_gdp_pcap_c$sobel_p_time, med2_gdp_pcap_c$b1, med2_gdp_pcap_c$se_b1))
cat(sprintf("indirect via CF_energy = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g) | b_energy (net of CF_time) = %.3f (SE %.3f)\n",
            med2_gdp_pcap_c$indirect_energy, med2_gdp_pcap_c$prop_mediated_energy * 100,
            med2_gdp_pcap_c$sobel_z_energy, med2_gdp_pcap_c$sobel_p_energy, med2_gdp_pcap_c$b2, med2_gdp_pcap_c$se_b2))

acme_time_stat_c <- function(d) {
  fit_a <- lm(log_cf ~ log_gdp_pcap, data = d)
  fit_b <- lm(log_protein ~ log_gdp_pcap + log_cf + log_cf_energy, data = d)
  unname(coef(fit_a)["log_gdp_pcap"]) * unname(coef(fit_b)["log_cf"])
}
acme_energy_stat_c <- function(d) {
  fit_a <- lm(log_cf_energy ~ log_gdp_pcap, data = d)
  fit_b <- lm(log_protein ~ log_gdp_pcap + log_cf + log_cf_energy, data = d)
  unname(coef(fit_a)["log_gdp_pcap"]) * unname(coef(fit_b)["log_cf_energy"])
}
point_acme_time_c   <- acme_time_stat_c(df2_collapsed)
point_acme_energy_c <- acme_energy_stat_c(df2_collapsed)
boot_acme_time_c    <- boot_stat(df2_collapsed, acme_time_stat_c, R)
boot_acme_energy_c  <- boot_stat(df2_collapsed, acme_energy_stat_c, R)

report_boot(boot_acme_time_c, point_acme_time_c,
            sprintf("[6] ACME via CF_time (all-work pair), RoW-collapsed (n=%d)", nrow(df2_collapsed)),
            sobel_p_for_comparison = med2_gdp_pcap_c$sobel_p_time)
report_boot(boot_acme_energy_c, point_acme_energy_c,
            sprintf("[7] ACME via CF_energy (all-work pair), RoW-collapsed (n=%d)", nrow(df2_collapsed)),
            sobel_p_for_comparison = med2_gdp_pcap_c$sobel_p_energy)

cat("\n#### What to look at ####\n")
cat("Compare [6]/[7] against 9o's n=33 dual-mediator bootstrap (CF_time: Sobel\n")
cat("p=0.42, bootstrap CI included zero -- null both ways; CF_energy: Sobel\n")
cat("p=0.124, bootstrap p=0.018, CI excluded zero -- the paper's most robust\n")
cat("single finding). If [7] still excludes zero at n=38, that finding holds up\n")
cat("with more low-income data, not just the EXIO-only leverage points. If [6]\n")
cat("stays null, CF_time's null result is confirmed, not just under-powered.\n")
