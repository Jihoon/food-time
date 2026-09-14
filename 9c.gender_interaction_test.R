#### Pooled interaction test: does the CF -> protein slope differ by gender? ####
#
# 9b ran the mediation/partial-correlation machinery separately on female-only
# (df_f / df_collapsed_f) and male-only (df_m / df_collapsed_m) subsamples --
# two n=33 (or n=38) fits, each underpowered on its own to detect a slope
# difference. This pools both genders into one regression per sample and
# tests the CF x gender interaction directly -- the correctly-powered way to
# ask "does this slope differ by gender" instead of eyeballing two separate
# p-values against each other.
#
# Wrinkle: log_protein (national protein supply) is IDENTICAL for the female
# and male row of the same country -- it's the same national total joined to
# both, not gender-specific consumption. So the two rows per country are NOT
# independent draws of y. Clustering standard errors by country is required
# for valid inference here, not a robustness nicety. Cluster-robust SE (CR1,
# matching Stata's default vce(cluster)) is implemented by hand below to
# avoid a new package dependency.
#
# Caveat: G=33 (EXIO-only) / G=38 (RoW-collapsed) clusters is on the low side
# for cluster-robust asymptotics (rule of thumb wants ~30-50+). Treat the
# p-values here as indicative, not exact -- a wild cluster bootstrap would be
# the next step up in rigor if this result is borderline and worth leaning on.
#
# Run AFTER 9b.gender_split_income_control.R in the same session (reuses its
# df_f, df_m, df_collapsed_f, df_collapsed_m).

cluster_robust_se <- function(fit, cluster) {
  X <- model.matrix(fit)
  u <- residuals(fit)
  n <- nrow(X); k <- ncol(X)
  G <- length(unique(cluster))
  meat <- matrix(0, k, k)
  for (g in unique(cluster)) {
    idx <- cluster == g
    Xg <- X[idx, , drop = FALSE]
    ug <- u[idx]
    score_g <- t(Xg) %*% ug
    meat <- meat + score_g %*% t(score_g)
  }
  bread <- solve(t(X) %*% X)
  V0 <- bread %*% meat %*% bread
  correction <- (G / (G - 1)) * ((n - 1) / (n - k))
  V <- correction * V0
  se <- sqrt(diag(V))
  names(se) <- colnames(X)
  list(V = V, se = se, G = G)
}

interaction_test <- function(df_f, df_m, label, gdp_col) {
  stacked <- bind_rows(
    df_f %>% mutate(gender = "Female"),
    df_m %>% mutate(gender = "Male")
  ) %>%
    mutate(gender = factor(gender, levels = c("Male", "Female")))  # Male = reference

  fml <- as.formula(paste0("log_protein ~ ", gdp_col, " + log_cf * gender"))
  fit <- lm(fml, data = stacked)
  cr  <- cluster_robust_se(fit, stacked$country)

  b    <- coef(fit)
  se   <- cr$se
  tval <- b / se
  pval <- 2 * pt(-abs(tval), df = cr$G - 1)

  cat(sprintf("\n==== %s (control = %s) ====\n", label, gdp_col))
  cat(sprintf("Clusters = %d, obs = %d\n", cr$G, nrow(stacked)))
  print(data.frame(coef = round(b, 3), cluster_se = round(se, 3),
                    t = round(tval, 2), p = signif(pval, 3)))

  int_name <- "log_cf:genderFemale"
  cat(sprintf("\nInteraction term -- does the CF->protein slope differ, female vs. male, net of %s?\n", gdp_col))
  cat(sprintf("Estimate = %.3f, cluster-robust SE = %.3f, t = %.2f, p = %.3g\n",
              b[int_name], se[int_name], tval[int_name], pval[int_name]))

  invisible(list(fit = fit, cr = cr, stacked = stacked))
}

cat("\n#### Pooled CF x gender interaction test ####\n")
res_1 <- interaction_test(df_f, df_m, "EXIO-only (n=66 pooled)", "log_gdp_pcap")
res_2 <- interaction_test(df_f, df_m, "EXIO-only (n=66 pooled)", "log_gdp_worker")
res_3 <- interaction_test(df_collapsed_f, df_collapsed_m, "RoW-collapsed (n=76 pooled)", "log_gdp_pcap")
res_4 <- interaction_test(df_collapsed_f, df_collapsed_m, "RoW-collapsed (n=76 pooled)", "log_gdp_worker")
