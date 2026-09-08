#### Does the CF -> protein relationship (net of GDP) get stronger at lower income? ####
#
# Companion to 9c (gender interaction). Same logic, but the moderator is GDP
# instead of gender: does the marginal effect of CF on protein supply change
# continuously with income level, rather than being one fixed global slope?
# A negative log_gdp_pcap:log_cf interaction would mean CF's slope on protein
# gets MORE negative (more binding) as GDP falls -- consistent with a fixed-
# time-budget constraint that bites at low income, and a ceiling effect
# (nutrition already adequate, further CF gains buy nothing) at high income.
#
# GDP is exogenous to this model (not the outcome variable), so a continuous
# interaction term is valid here -- unlike a protein-supply-level version of
# this same question, which would have to interact log_cf with a term built
# from the dependent variable itself (circular). This replaces that attempt.
#
# Run AFTER 9.capability_set_income_control.R in the same session (reuses its
# `df`, `df_collapsed`).

gdp_interaction_test <- function(data, label, gdp_col) {
  fml <- as.formula(paste0("log_protein ~ ", gdp_col, " * log_cf"))
  fit <- lm(fml, data = data)
  s <- summary(fit)

  cat(sprintf("\n==== %s (control = %s, n = %d) ====\n", label, gdp_col, nrow(data)))
  print(round(s$coefficients, 4))

  int_name <- paste0(gdp_col, ":log_cf")
  b    <- s$coefficients[int_name, "Estimate"]
  se   <- s$coefficients[int_name, "Std. Error"]
  tval <- s$coefficients[int_name, "t value"]
  pval <- s$coefficients[int_name, "Pr(>|t|)"]
  cat(sprintf("\nInteraction (%s): does CF's slope on protein change with income?\n", int_name))
  cat(sprintf("Estimate = %.4f, SE = %.4f, t = %.2f, p = %.3g\n", b, se, tval, pval))

  invisible(fit)
}

cat("\n#### GDP x CF interaction: EXIO-only (n=33) ####\n")
gdp_interaction_test(df, "EXIO-only", "log_gdp_pcap")
gdp_interaction_test(df, "EXIO-only", "log_gdp_worker")

cat("\n\n#### GDP x CF interaction: RoW-collapsed (n=38) ####\n")
gdp_interaction_test(df_collapsed, "RoW-collapsed", "log_gdp_pcap")
gdp_interaction_test(df_collapsed, "RoW-collapsed", "log_gdp_worker")
