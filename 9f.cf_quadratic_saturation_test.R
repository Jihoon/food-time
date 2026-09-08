#### Does CF's own effect on protein supply saturate, and does that curvature itself depend on GDP? ####
#
# Combines 9d (GDP x CF linear interaction: does the SLOPE change with
# income) and a plain quadratic-in-CF test (does the relationship flatten
# out at low CF, independent of income) into one nested comparison, plus the
# term that actually answers the question asked: GDP x CF^2 -- does the
# CURVATURE itself change with income, not just the slope?
#
# log_cf is CENTERED before squaring/interacting to reduce collinearity
# between the linear and quadratic terms -- log_cf and log_gdp_pcap are
# already strongly correlated (elasticity -0.52, R^2=0.67), so leaving it
# uncentered would compound that further.
#
# Fitted as three NESTED models (linear-only -> +quadratic -> +quadratic x
# GDP) compared via ANOVA, rather than read off individual t-stats in the
# fully saturated model -- at n=33-38 with 5-6 correlated regressors, the
# saturated model's own coefficients are likely to be numerically unstable
# even where the nested F-test is informative. Treat this as exploratory:
# there isn't much power left once curvature AND its interaction with GDP
# both have to be estimated off this sample.
#
# Run AFTER 9.capability_set_income_control.R in the same session (reuses
# `df`, `df_collapsed`).

cf_curvature_nested_test <- function(data, label, gdp_col) {
  data <- data %>% mutate(log_cf_c = log_cf - mean(log_cf))
  fml_linear <- as.formula(paste0("log_protein ~ ", gdp_col, " * log_cf_c"))
  fml_quad   <- as.formula(paste0("log_protein ~ ", gdp_col, " * log_cf_c + I(log_cf_c^2)"))
  fml_full   <- as.formula(paste0("log_protein ~ ", gdp_col, " * (log_cf_c + I(log_cf_c^2))"))

  fit_linear <- lm(fml_linear, data = data)
  fit_quad   <- lm(fml_quad,   data = data)
  fit_full   <- lm(fml_full,   data = data)

  cat(sprintf("\n==== %s (control = %s, n = %d) ====\n", label, gdp_col, nrow(data)))
  cat("-- Full model (GDP x CF + GDP x CF^2) coefficients --\n")
  print(round(summary(fit_full)$coefficients, 4))

  cat("\n-- Nested model comparison (linear -> +CF^2 main effect -> +CF^2 x GDP) --\n")
  print(anova(fit_linear, fit_quad, fit_full))

  invisible(list(fit_linear = fit_linear, fit_quad = fit_quad, fit_full = fit_full))
}

cat("\n#### CF curvature x GDP nested test ####\n")
cf_curvature_nested_test(df, "EXIO-only", "log_gdp_pcap")
cf_curvature_nested_test(df, "EXIO-only", "log_gdp_worker")
cf_curvature_nested_test(df_collapsed, "RoW-collapsed", "log_gdp_pcap")
cf_curvature_nested_test(df_collapsed, "RoW-collapsed", "log_gdp_worker")
