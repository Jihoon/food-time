#### Does energy CF mediate GDP -> protein alongside time CF, or is time CF ####
#### doing the work once energy CF is also in the model? ####
#
# Extends 9.capability_set_income_control.R's single-mediator Baron & Kenny
# decomposition (GDP -> CF_time -> protein) to a parallel two-mediator model
# (Hayes PROCESS model 4 logic): GDP -> {CF_time, CF_energy} -> protein,
# with both mediators entered simultaneously in the outcome model so each
# mediator's indirect effect is estimated net of the other, not treated as
# if it were the only path.
#
# Run AFTER 2.analyze_result.R (reuses tradeoff_protein_allwork_consump) AND
# AFTER 9.capability_set_income_control.R (reuses df, wdi_latest, `year`).

library(tidyverse)

#### 1. Extend cf_df with energy CF (mj_per_50g_protein), same aggregation
#### as the existing hr_per_50g_protein (summed across sector x origin
#### buckets per country) ####

cf_df_energy <- tradeoff_protein_allwork_consump %>%
  filter(!is_row) %>%
  group_by(country) %>%
  summarise(hr_per_50g_protein = sum(hr_per_50g_protein, na.rm = TRUE),
            mj_per_50g_protein = sum(mj_per_50g_protein, na.rm = TRUE),
            g_protein_per_cap_day = first(g_protein_per_cap_day),
            .groups = "drop")

df2 <- cf_df_energy %>%
  inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
  filter(hr_per_50g_protein > 0, mj_per_50g_protein > 0, g_protein_per_cap_day > 0) %>%
  mutate(log_cf         = log(hr_per_50g_protein),
         log_cf_energy  = log(mj_per_50g_protein),
         log_protein    = log(g_protein_per_cap_day),
         log_gdp_pcap   = log(gdp_pcap_ppp),
         log_gdp_worker = log(gdp_per_worker))

cat(sprintf("Countries with CF_time, CF_energy, protein supply, and >=1 income control: %d\n", nrow(df2)))

#### 2. How correlated are the two mediators? (matters for interpreting the ####
#### parallel model -- if they move together, the split between them is    ####
#### less stable than either indirect effect alone)                        ####

cor_cf_time_energy <- cor.test(df2$log_cf, df2$log_cf_energy)
cat(sprintf("\nCorrelation between CF_time and CF_energy (log-log): r = %.3f (p = %.3g, n = %d)\n",
            cor_cf_time_energy$estimate, cor_cf_time_energy$p.value, nrow(df2)))

#### 3. Parallel two-mediator decomposition ####
# total    : Y ~ X
# path a1  : M1 ~ X            (GDP -> CF_time)
# path a2  : M2 ~ X            (GDP -> CF_energy)
# path b1,b2,direct : Y ~ X + M1 + M2   (both mediators entered together --
#   b1/b2 are each mediator's effect on Y net of the other, and `direct` is
#   GDP's effect on protein net of BOTH mediators, not just one)
# indirect_1 = a1*b1, indirect_2 = a2*b2; total = direct + indirect_1 + indirect_2
#   (still an exact algebraic identity for OLS, same as the single-mediator case)

parallel_mediation_decomp <- function(df, x_col, m1_col, m2_col, y_col) {
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

  total  <- unname(coef(fit_total)["x"])
  direct <- unname(coef(fit_out)["x"])   # net of BOTH mediators now

  ind1 <- a1 * b1; se_ind1 <- sqrt(b1^2 * se_a1^2 + a1^2 * se_b1^2); z1 <- ind1 / se_ind1
  ind2 <- a2 * b2; se_ind2 <- sqrt(b2^2 * se_a2^2 + a2^2 * se_b2^2); z2 <- ind2 / se_ind2

  list(n = length(x), total = total, direct = direct,
       indirect_time = ind1, prop_mediated_time = ind1 / total,
       sobel_z_time = z1, sobel_p_time = 2 * (1 - pnorm(abs(z1))),
       indirect_energy = ind2, prop_mediated_energy = ind2 / total,
       sobel_z_energy = z2, sobel_p_energy = 2 * (1 - pnorm(abs(z2))),
       b1 = b1, se_b1 = se_b1, b2 = b2, se_b2 = se_b2)
}

med2_gdp_pcap   <- parallel_mediation_decomp(df2, "log_gdp_pcap",   "log_cf", "log_cf_energy", "log_protein")
med2_gdp_worker <- parallel_mediation_decomp(df2, "log_gdp_worker", "log_cf", "log_cf_energy", "log_protein")

cat("\n---- Parallel two-mediator decomposition: GDP -> {CF_time, CF_energy} -> protein supply ----\n")
for (m in list(list(label = "GDP per capita (PPP)",   r = med2_gdp_pcap),
                list(label = "GDP per worker (labor prod.)", r = med2_gdp_worker))) {
  r <- m$r
  cat(sprintf("\nVia %s (n = %d):\n", m$label, r$n))
  cat(sprintf("  total = %.3f | direct (net of both mediators) = %.3f\n", r$total, r$direct))
  cat(sprintf("  indirect via CF_time   = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g) | b_time (Y~..+CF_time, net of CF_energy)   = %.3f (SE %.3f)\n",
              r$indirect_time, r$prop_mediated_time * 100, r$sobel_z_time, r$sobel_p_time, r$b1, r$se_b1))
  cat(sprintf("  indirect via CF_energy = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g) | b_energy (Y~..+CF_energy, net of CF_time) = %.3f (SE %.3f)\n",
              r$indirect_energy, r$prop_mediated_energy * 100, r$sobel_z_energy, r$sobel_p_energy, r$b2, r$se_b2))
}

#### 4. What to look at in the output ####
# - If CF_time and CF_energy are highly correlated (step 2), b1 and b2 above
#   are each other's confound: a b_time that goes to ~0 (or flips sign) once
#   CF_energy is added, versus the single-mediator b1 from script 9, would
#   mean energy CF was doing (or masking) some of what looked like a time-CF
#   effect there -- worth comparing r$b1 here against script 9's single-
#   mediator b (implied by med_gdp_pcap$indirect / med_gdp_pcap$total * ... ,
#   or just re-run mediation_decomp(df, ..., "log_cf", "log_protein") side by
#   side).
# - indirect_time + indirect_energy + direct should equal total (algebraic
#   check, not a finding) -- if it visibly doesn't, something's off in the
#   merge (likely a country dropped by the mj_per_50g_protein > 0 filter that
#   wasn't dropped in script 9's df).
