#### One-row gender test: does female CF predict protein differently from male CF? ####
#
# Alternative to 9c's stacked interaction model. 9c duplicates each country's
# log_protein across a female and a male row and then clusters SEs by country
# to undo the duplication. This asks the same question without the
# duplication: one row per country, both gender CFs on the right-hand side,
# and the gender difference is the contrast c - d.
#
#   log_protein = a + b*log_gdp + c*log_cf_f + d*log_cf_m + e
#
# c - d here plays the role of 9c's log_cf:genderFemale interaction term.
# Ordinary OLS SEs are valid (n independent rows), so no clustering.
#
# The limitation both forms share: log_cf_f and log_cf_m are both driven by
# income and by the country's total CF, so c and d are collinear and their
# difference is imprecise. The correlation between the two CFs is printed so
# the reader can see how much of the imprecision is that.
#
# Run AFTER 9b.gender_split_income_control.R in the same session (reuses its
# df_f, df_m, df_collapsed_f, df_collapsed_m).

horse_race <- function(df_f, df_m, label, gdp_col = "log_gdp_pcap") {
  one_row <- inner_join(
    df_f %>% select(country, log_protein, all_of(gdp_col), log_cf_f = log_cf),
    df_m %>% select(country, log_cf_m = log_cf),
    by = "country"
  )

  fml <- as.formula(paste0("log_protein ~ ", gdp_col, " + log_cf_f + log_cf_m"))
  fit <- lm(fml, data = one_row)
  V   <- vcov(fit)
  b   <- coef(fit)

  diff    <- b["log_cf_f"] - b["log_cf_m"]
  se_diff <- sqrt(V["log_cf_f", "log_cf_f"] + V["log_cf_m", "log_cf_m"] -
                  2 * V["log_cf_f", "log_cf_m"])
  t_diff  <- diff / se_diff
  p_diff  <- 2 * pt(-abs(t_diff), df = fit$df.residual)

  r_fm <- cor(one_row$log_cf_f, one_row$log_cf_m)

  cat(sprintf("\n==== %s (control = %s, n = %d) ====\n", label, gdp_col, nrow(one_row)))
  print(summary(fit)$coefficients)
  cat(sprintf("\ncor(log_cf_f, log_cf_m) = %.3f\n", r_fm))
  cat(sprintf("c - d (female minus male slope) = %.3f, SE = %.3f, t = %.2f, p = %.3g, df = %d\n",
              diff, se_diff, t_diff, p_diff, fit$df.residual))

  invisible(list(fit = fit, one_row = one_row, diff = diff, se_diff = se_diff,
                 t = t_diff, p = p_diff, r_fm = r_fm))
}

cat("\n#### One-row gender horse race: female vs. male CF, net of GDP ####\n")
hr_1 <- horse_race(df_f, df_m, "EXIO-only")
hr_2 <- horse_race(df_collapsed_f, df_collapsed_m, "RoW-collapsed")
