#### Does the CF_energy -> protein relationship (net of GDP) get stronger ####
#### at lower income, the same way CF_time's does? ####
#
# Companion to 9d (GDP x CF_time interaction). Same logic and same
# power-matched specification (4 parameters: intercept, GDP, CF_energy,
# GDP:CF_energy) -- NOT the joint dual-mediator interaction model
# (GDP*CF_time + GDP*CF_energy, 6 parameters), which at n=33 and with
# CF_time/CF_energy correlated at r=-0.53 (9m) would compound an already
# marginal power problem. Kept as a separate single-mediator-style test so
# it is directly comparable to 9d's result, not a weaker joint one.
#
# EXIO-only (n=33) only -- no RoW-collapsed (n=38) version here. That would
# need the same population-weighted regional aggregation script 9 does for
# CF_time's df_collapsed (row_lookup / pop_data_yr / row_full_pop, applied
# to a per-capita energy variable upstream of mj_per_50g_protein), which
# isn't built out for energy. Flagging as a gap, not done silently.
#
# Run AFTER 2.analyze_result.R, 9.capability_set_income_control.R, AND
# 9m.energy_time_dual_mediator.R in the same session (reuses df2).

gdp_energy_interaction_test <- function(data, label, gdp_col) {
  fml <- as.formula(paste0("log_protein ~ ", gdp_col, " * log_cf_energy"))
  fit <- lm(fml, data = data)
  s <- summary(fit)

  cat(sprintf("\n==== %s (control = %s, n = %d) ====\n", label, gdp_col, nrow(data)))
  print(round(s$coefficients, 4))

  int_name <- paste0(gdp_col, ":log_cf_energy")
  b    <- s$coefficients[int_name, "Estimate"]
  se   <- s$coefficients[int_name, "Std. Error"]
  tval <- s$coefficients[int_name, "t value"]
  pval <- s$coefficients[int_name, "Pr(>|t|)"]
  cat(sprintf("\nInteraction (%s): does CF_energy's slope on protein change with income?\n", int_name))
  cat(sprintf("Estimate = %.4f, SE = %.4f, t = %.2f, p = %.3g\n", b, se, tval, pval))

  invisible(fit)
}

cat("\n#### GDP x CF_energy interaction: EXIO-only (n=33) ####\n")
gdp_energy_interaction_test(df2, "EXIO-only", "log_gdp_pcap")
gdp_energy_interaction_test(df2, "EXIO-only", "log_gdp_worker")
