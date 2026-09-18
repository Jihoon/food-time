#### Does the CF_paid -> protein relationship (net of GDP) vary with income? ####
#
# Companion to 9d (GDP x CF_time) and 9n (GDP x CF_energy). Same 4-parameter
# specification (intercept, GDP, CF_paid, GDP:CF_paid), same marginal-effect
# evaluation points ($2,000 and $49,000 per capita), same case-resampling
# bootstrap on the interaction term as 9o. Both country samples, since 9t
# already builds cf_paid_df33 and cf_paid_df38.
#
# Motivation: the manuscript reports income heterogeneity for CF_time and
# CF_energy only. CF_paid is the one time-measure whose path b is not null,
# so whether its slope is steeper at low income is a distinct question and
# has not been asked.
#
# Run AFTER 9o.bootstrap_mediation_and_interaction.R (boot_stat, report_boot,
# R) and 9t.energy_paid_labor_dual_mediator.R (cf_paid_df33, cf_paid_df38)
# in the same session.

gdp_paid_interaction_test <- function(data, label, gdp_col = "log_gdp_pcap",
                                      eval_gdp = c(2000, 49000)) {
  fml <- as.formula(paste0("log_protein ~ ", gdp_col, " * log_cf_paid"))
  fit <- lm(fml, data = data)
  s   <- summary(fit)
  V   <- vcov(fit)
  b   <- coef(fit)

  int_name <- paste0(gdp_col, ":log_cf_paid")

  cat(sprintf("\n==== %s (control = %s, n = %d) ====\n", label, gdp_col, nrow(data)))
  print(round(s$coefficients, 4))

  cat(sprintf("\nInteraction (%s): does CF_paid's slope on protein change with income?\n", int_name))
  cat(sprintf("Estimate = %.4f, SE = %.4f, t = %.2f, p = %.3g\n",
              s$coefficients[int_name, "Estimate"], s$coefficients[int_name, "Std. Error"],
              s$coefficients[int_name, "t value"], s$coefficients[int_name, "Pr(>|t|)"]))

  cat("\nMarginal effect of log_cf_paid on log_protein at representative incomes:\n")
  for (g in eval_gdp) {
    x  <- log(g)
    me <- b["log_cf_paid"] + b[int_name] * x
    se <- sqrt(V["log_cf_paid", "log_cf_paid"] + x^2 * V[int_name, int_name] +
               2 * x * V["log_cf_paid", int_name])
    cat(sprintf("  GDP/cap = $%6d: dY/dM = %.3f (SE = %.3f, t = %.2f, p = %.3g)\n",
                g, me, se, me / se, 2 * pt(-abs(me / se), df = fit$df.residual)))
  }

  interaction_stat <- function(d) {
    f <- lm(fml, data = d)
    unname(coef(f)[int_name])
  }
  boot_int <- boot_stat(data, interaction_stat, R)
  report_boot(boot_int, unname(b[int_name]),
              sprintf("GDP x CF_paid interaction, %s (n = %d)", label, nrow(data)),
              sobel_p_for_comparison = s$coefficients[int_name, "Pr(>|t|)"])

  invisible(list(fit = fit, boot = boot_int))
}

cat("\n#### GDP x CF_paid interaction ####\n")
int_paid_33 <- gdp_paid_interaction_test(cf_paid_df33, "EXIO-only")
int_paid_38 <- gdp_paid_interaction_test(cf_paid_df38, "RoW-collapsed")

cat("\n#### What to look at ####\n")
cat("Compare with CF_time (9d/9o): slope -0.44 at $2,000 vs -0.04 at $49,000,\n")
cat("interaction p = 0.096 at n = 38, bootstrap CI includes zero. If CF_paid's\n")
cat("low-income slope is clearly more negative than its high-income one AND the\n")
cat("bootstrap CI on the interaction excludes zero, the low-income threshold\n")
cat("claim can be extended to CF_paid. If the CI straddles zero, it stays a\n")
cat("directed hypothesis, as for CF_time.\n")
