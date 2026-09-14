#### GDP x CF interaction, controlling for gender (not interacting with it) ####
#
# Leaner alternative to a full three-way GDP x CF x gender model. 9d tested
# the income-heterogeneity pattern on the combined-gender CF (n=33/38). This
# adds gender as an ADDITIVE control on the stacked male/female data (n=66/76,
# same construction as 9c) -- asks "does the income pattern hold once gender
# is controlled for", not "does the income pattern itself differ by gender"
# (that would be the log_gdp:log_cf:gender three-way term, deliberately left
# out here -- flagged as underpowered at this n when discussed).
#
# Country-clustered SEs are required for the same reason as in 9c: log_protein
# is identical for the male-row and female-row of the same country, so the
# two rows are not independent draws.
#
# Run AFTER 9c.gender_interaction_test.R in the same session (reuses `df_f`,
# `df_m`, `df_collapsed_f`, `df_collapsed_m`, `cluster_robust_se()`).

gdp_cf_gender_controlled <- function(df_f, df_m, label, gdp_col) {
  stacked <- bind_rows(
    df_f %>% mutate(gender = "Female"),
    df_m %>% mutate(gender = "Male")
  ) %>%
    mutate(gender = factor(gender, levels = c("Male", "Female")))

  fml <- as.formula(paste0("log_protein ~ ", gdp_col, " * log_cf + gender"))
  fit <- lm(fml, data = stacked)
  cr  <- cluster_robust_se(fit, stacked$country)

  b    <- coef(fit)
  se   <- cr$se
  tval <- b / se
  pval <- 2 * pt(-abs(tval), df = cr$G - 1)

  cat(sprintf("\n==== %s (control = %s, gender as additive control, n = %d, clusters = %d) ====\n",
              label, gdp_col, nrow(stacked), cr$G))
  print(data.frame(coef = round(b, 4), cluster_se = round(se, 4),
                    t = round(tval, 2), p = signif(pval, 3)))

  int_name <- paste0(gdp_col, ":log_cf")
  cat(sprintf("\nIncome-heterogeneity term (%s), controlling for gender:\n", int_name))
  cat(sprintf("Estimate = %.4f, cluster-robust SE = %.4f, t = %.2f, p = %.3g\n",
              b[int_name], se[int_name], tval[int_name], pval[int_name]))

  invisible(list(fit = fit, cr = cr))
}

cat("\n#### GDP x CF interaction, gender-controlled: EXIO-only (n=66 pooled) ####\n")
gdp_cf_gender_controlled(df_f, df_m, "EXIO-only", "log_gdp_pcap")
gdp_cf_gender_controlled(df_f, df_m, "EXIO-only", "log_gdp_worker")

cat("\n\n#### GDP x CF interaction, gender-controlled: RoW-collapsed (n=76 pooled) ####\n")
gdp_cf_gender_controlled(df_collapsed_f, df_collapsed_m, "RoW-collapsed", "log_gdp_pcap")
gdp_cf_gender_controlled(df_collapsed_f, df_collapsed_m, "RoW-collapsed", "log_gdp_worker")
