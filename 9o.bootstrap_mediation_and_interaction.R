#### Bootstrap check: does the Sobel normal-approximation hold up at n=33-38? ####
#
# Every Sobel z/p reported so far (9, 9e, 9m) uses the classic delta-method
# SE: se_indirect = sqrt(b^2*se_a^2 + a^2*se_b^2), then a normal-approx p.
# That approximation is known to be weakest exactly where we're leaning on
# it hardest: small n, and indirect effects (products of two estimated
# coefficients) whose true sampling distribution is skewed, not normal.
# This does a straightforward case-resampling bootstrap (no `mediation`
# package dependency, same "implement by hand" choice 9c already made for
# cluster-robust SEs) for the results closest to conventional significance:
# the CF_energy x GDP interaction (9n, Sobel-adjacent p=0.049) and, for
# comparison, the CF_time x GDP interaction (9d, p=0.163) and both
# dual-mediator indirect effects (9m).
#
# Percentile CI is the primary output (does it exclude zero), with a
# symmetric two-sided bootstrap p (2*min(prop below 0, prop above 0)) for
# direct comparison against the already-reported Sobel p's. R=5000
# replications, case (row) resampling with replacement -- appropriate here
# since X is not fixed by design (observational country-level data), so
# resampling rows (not just residuals) is the right bootstrap for this
# design.
#
# Run AFTER 2.analyze_result.R, 9.capability_set_income_control.R,
# 9m.energy_time_dual_mediator.R in the same session (reuses df, df2).

set.seed(1)
R <- 5000

boot_stat <- function(data, statistic_fn, R = 5000) {
  n <- nrow(data)
  stats <- numeric(R)
  for (i in seq_len(R)) {
    idx <- sample.int(n, n, replace = TRUE)
    stats[i] <- statistic_fn(data[idx, , drop = FALSE])
  }
  stats
}

report_boot <- function(stats, point_estimate, label, sobel_p_for_comparison = NA) {
  ci <- quantile(stats, c(0.025, 0.975), na.rm = TRUE)
  boot_se <- sd(stats, na.rm = TRUE)
  prop_below0 <- mean(stats < 0, na.rm = TRUE)
  prop_above0 <- mean(stats > 0, na.rm = TRUE)
  boot_p <- 2 * min(prop_below0, prop_above0)
  cat(sprintf("\n%s\n", label))
  cat(sprintf("  point estimate = %.4f | bootstrap SE = %.4f | 95%% percentile CI = [%.4f, %.4f]\n",
              point_estimate, boot_se, ci[1], ci[2]))
  cat(sprintf("  bootstrap two-sided p = %.4f%s | CI excludes zero: %s\n",
              boot_p,
              if (!is.na(sobel_p_for_comparison)) sprintf(" (Sobel/asymptotic p was %.4f)", sobel_p_for_comparison) else "",
              ifelse(ci[1] > 0 || ci[2] < 0, "YES", "no")))
  invisible(list(ci = ci, se = boot_se, p = boot_p))
}

#### 1. Single-mediator indirect effect (CF_time only), GDP per capita, n=33 ####

acme_time_stat <- function(d) {
  fit_a <- lm(log_cf ~ log_gdp_pcap, data = d)
  fit_b <- lm(log_protein ~ log_gdp_pcap + log_cf, data = d)
  unname(coef(fit_a)["log_gdp_pcap"]) * unname(coef(fit_b)["log_cf"])
}
point_acme_time <- acme_time_stat(df)
boot_acme_time <- boot_stat(df, acme_time_stat, R)
report_boot(boot_acme_time, point_acme_time,
            "[1] Single-mediator ACME, CF_time (GDP per capita, n=33)",
            sobel_p_for_comparison = med_gdp_pcap$sobel_p)

#### 2. Dual-mediator indirect effects (CF_time & CF_energy), GDP per capita, n=33 ####

acme_dual_stat <- function(which_m) {
  function(d) {
    fit_m1 <- lm(log_cf ~ log_gdp_pcap, data = d)
    fit_m2 <- lm(log_cf_energy ~ log_gdp_pcap, data = d)
    fit_y  <- lm(log_protein ~ log_gdp_pcap + log_cf + log_cf_energy, data = d)
    if (which_m == "time") {
      unname(coef(fit_m1)["log_gdp_pcap"]) * unname(coef(fit_y)["log_cf"])
    } else {
      unname(coef(fit_m2)["log_gdp_pcap"]) * unname(coef(fit_y)["log_cf_energy"])
    }
  }
}
point_acme_dual_time   <- acme_dual_stat("time")(df2)
point_acme_dual_energy <- acme_dual_stat("energy")(df2)
boot_acme_dual_time   <- boot_stat(df2, acme_dual_stat("time"), R)
boot_acme_dual_energy <- boot_stat(df2, acme_dual_stat("energy"), R)
report_boot(boot_acme_dual_time, point_acme_dual_time,
            "[2a] Dual-mediator ACME, CF_time net of CF_energy (GDP per capita, n=33)",
            sobel_p_for_comparison = med2_gdp_pcap$sobel_p_time)
report_boot(boot_acme_dual_energy, point_acme_dual_energy,
            "[2b] Dual-mediator ACME, CF_energy net of CF_time (GDP per capita, n=33)",
            sobel_p_for_comparison = med2_gdp_pcap$sobel_p_energy)

#### 3. GDP x CF_time interaction (9d), GDP per capita, n=33 ####

interaction_time_stat <- function(d) {
  fit <- lm(log_protein ~ log_gdp_pcap * log_cf, data = d)
  unname(coef(fit)["log_gdp_pcap:log_cf"])
}
point_int_time <- interaction_time_stat(df)
boot_int_time  <- boot_stat(df, interaction_time_stat, R)
fit_int_time_orig <- lm(log_protein ~ log_gdp_pcap * log_cf, data = df)
sobel_p_int_time <- summary(fit_int_time_orig)$coefficients["log_gdp_pcap:log_cf", "Pr(>|t|)"]
report_boot(boot_int_time, point_int_time,
            "[3] GDP x CF_time interaction (GDP per capita, n=33)",
            sobel_p_for_comparison = sobel_p_int_time)

#### 4. GDP x CF_energy interaction (9n), GDP per capita, n=33 -- the borderline one ####

interaction_energy_stat <- function(d) {
  fit <- lm(log_protein ~ log_gdp_pcap * log_cf_energy, data = d)
  unname(coef(fit)["log_gdp_pcap:log_cf_energy"])
}
point_int_energy <- interaction_energy_stat(df2)
boot_int_energy  <- boot_stat(df2, interaction_energy_stat, R)
fit_int_energy_orig <- lm(log_protein ~ log_gdp_pcap * log_cf_energy, data = df2)
sobel_p_int_energy <- summary(fit_int_energy_orig)$coefficients["log_gdp_pcap:log_cf_energy", "Pr(>|t|)"]
report_boot(boot_int_energy, point_int_energy,
            "[4] GDP x CF_energy interaction (GDP per capita, n=33) -- the p=0.049 result",
            sobel_p_for_comparison = sobel_p_int_energy)

cat("\n#### What to look at ####\n")
cat("For [4] specifically: does the 95% percentile CI exclude zero? If yes,\n")
cat("the significant result survives a small-sample-appropriate bootstrap and\n")
cat("we can lean on it a bit more. If the CI straddles zero despite the\n")
cat("t-test p landing at 0.049, that's real evidence the normal approximation\n")
cat("was flattering this result and the hedge in the draft should get\n")
cat("stronger, not softer.\n")
