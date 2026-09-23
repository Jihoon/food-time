#### Table 2: every mediation estimate in the paper, computed one way, in one place ####
#
# The manuscript's mediation numbers were produced across 9, 9e, 9m, 9o, 9p, 9s
# and 9t, at different times, with bootstraps from separate runs. This script
# recomputes every Table 2 cell from the same session objects with one seed,
# so the table and the prose can be checked against a single output.
#
# Rows: single-mediator models for CF_time, CF_paid, CF_energy, CF_paid,dom,
# CF_domestic-effort; two-mediator models CF_time + CF_energy and
# CF_paid + CF_energy (one row per mediator). Each row at n = 33 (EXIO-only)
# and n = 38 (RoW-collapsed) where the variant exists at that sample.
# Columns: n, path a (GDP -> mediator elasticity, p), partial r (mediator ->
# protein, GDP fixed; single-mediator rows only), indirect effect a*b,
# share of GDP's total effect, Sobel p, bootstrap 95% CI and two-sided p.
#
# Sources its prerequisites itself (9, 9e, 9m, 9o, 9p, 9s, 9t, in dependency
# order). The upstream objects come from data/objects_for_9_2020.RData
# (2.2.save_objects_for_9.R) when present, else from running 2.analyze_result.R.
# Writes results/table2_mediation_summary.csv and, from the block at the end,
# results/table3_income_interaction.csv (GDP x CF interaction, Table 3). Re-running
# in the same session skips the prerequisite scripts.

library(tidyverse)

PROJ <- "h:/MyDocuments/Projects/REMASS/food time"
setwd(PROJ)

if (!exists("tradeoff_protein_allwork_consump")) {
  saved <- "data/objects_for_9_2020.RData"            # written by 2.2.save_objects_for_9.R
  if (file.exists(saved)) load(saved) else source("2.analyze_result.R")
}
needed <- c("df", "df_collapsed", "df_domestic_econ", "df2", "df2_collapsed", "cf2_df", "cf2_df_collapsed",
            "cf_paid_df33", "cf_paid_df38", "cf_paid_energy_df", "cf_paid_energy_df38", "partial_cor")
if (!all(vapply(needed, exists, logical(1)))) {   # skip the slow prerequisites when re-running in the same session
  for (f in c("9.capability_set_income_control.R",      # df, df_collapsed, wdi_latest, partial_cor(), mediation_decomp()
              "9e.domestic_econ_cf_mediation_test.R",   # df_domestic_econ (CF_paid,dom)
              "9m.energy_time_dual_mediator.R",         # df2 (CF_time + CF_energy, n = 33)
              "9o.bootstrap_mediation_and_interaction.R",  # boot_stat(), report_boot(), R; needs 9 and 9m
              "9p.energy_cf_row_collapsed.R",           # df2_collapsed; needs 9o
              "9s.mediation_two_cf_comparison.R",       # cf2_df, cf2_df_collapsed (domestic-effort CF)
              "9t.energy_paid_labor_dual_mediator.R")) {  # cf_paid_df33/38, cf_paid_energy_df/_df38; needs 9o, 9p
    cat("\n\n######## sourcing", f, "########\n")
    source(f)
  }
}

set.seed(20260921)  # after the sourced scripts, which set their own seeds
R_BOOT <- 5000
X <- "log_gdp_pcap"; Y <- "log_protein"

boot_ci <- function(d, stat_fn) {
  n <- nrow(d)
  s <- replicate(R_BOOT, stat_fn(d[sample.int(n, n, replace = TRUE), , drop = FALSE]))
  ci <- quantile(s, c(0.025, 0.975), na.rm = TRUE)
  c(lo = unname(ci[1]), hi = unname(ci[2]),
    p = 2 * min(mean(s < 0, na.rm = TRUE), mean(s > 0, na.rm = TRUE)))
}

single_row <- function(d, m, variant, sample) {
  d <- d[complete.cases(d[c(X, m, Y)]), ]
  fa <- lm(reformulate(X, m), d); fb <- lm(reformulate(c(X, m), Y), d); ft <- lm(reformulate(X, Y), d)
  a <- coef(fa)[[X]]; sa <- summary(fa)$coefficients[X, "Std. Error"]; pa <- summary(fa)$coefficients[X, "Pr(>|t|)"]
  b <- coef(fb)[[m]]; sb <- summary(fb)$coefficients[m, "Std. Error"]
  tot <- coef(ft)[[X]]
  pc <- partial_cor(d[[m]], d[[Y]], d[[X]])
  z <- a * b / sqrt(b^2 * sa^2 + a^2 * sb^2)
  bt <- boot_ci(d, function(s) coef(lm(reformulate(X, m), s))[[X]] * coef(lm(reformulate(c(X, m), Y), s))[[m]])
  tibble(model = "single", variant, mediator = variant, sample, n = nrow(d),
         path_a = a, path_a_p = pa, partial_r = pc$r, partial_r_p = pc$p,
         indirect = a * b, share = a * b / tot, sobel_p = 2 * (1 - pnorm(abs(z))),
         boot_lo = bt[["lo"]], boot_hi = bt[["hi"]], boot_p = bt[["p"]])
}

dual_rows <- function(d, m1, m2, name1, name2, sample) {
  d <- d[complete.cases(d[c(X, m1, m2, Y)]), ]
  ft <- lm(reformulate(X, Y), d); fo <- lm(reformulate(c(X, m1, m2), Y), d); tot <- coef(ft)[[X]]
  map_dfr(list(c(m1, name1), c(m2, name2)), function(mm) {
    m <- mm[1]
    fa <- lm(reformulate(X, m), d)
    a <- coef(fa)[[X]]; sa <- summary(fa)$coefficients[X, "Std. Error"]; pa <- summary(fa)$coefficients[X, "Pr(>|t|)"]
    b <- coef(fo)[[m]]; sb <- summary(fo)$coefficients[m, "Std. Error"]
    z <- a * b / sqrt(b^2 * sa^2 + a^2 * sb^2)
    bt <- boot_ci(d, function(s) coef(lm(reformulate(X, m), s))[[X]] * coef(lm(reformulate(c(X, m1, m2), Y), s))[[m]])
    tibble(model = paste(name1, "+", name2), variant = paste(name1, "+", name2), mediator = mm[2], sample, n = nrow(d),
           path_a = a, path_a_p = pa, partial_r = NA_real_, partial_r_p = NA_real_,
           indirect = a * b, share = a * b / tot, sobel_p = 2 * (1 - pnorm(abs(z))),
           boot_lo = bt[["lo"]], boot_hi = bt[["hi"]], boot_p = bt[["p"]])
  })
}

tab <- bind_rows(
  single_row(df,                "log_cf",        "CF_time",            "EXIO-only"),
  single_row(df_collapsed,      "log_cf",        "CF_time",            "RoW-collapsed"),
  single_row(cf_paid_df33,      "log_cf_paid",   "CF_paid",            "EXIO-only"),
  single_row(cf_paid_df38,      "log_cf_paid",   "CF_paid",            "RoW-collapsed"),
  single_row(df2,               "log_cf_energy", "CF_energy",          "EXIO-only"),
  single_row(df2_collapsed,     "log_cf_energy", "CF_energy",          "RoW-collapsed"),
  single_row(df_domestic_econ,  "log_cf",        "CF_paid,dom",        "EXIO-only"),
  single_row(cf2_df,            "log_cf2",       "CF_domestic-effort", "EXIO-only"),
  single_row(cf2_df_collapsed,  "log_cf2",       "CF_domestic-effort", "RoW-collapsed"),
  dual_rows(df2,                  "log_cf",      "log_cf_energy", "CF_time", "CF_energy", "EXIO-only"),
  dual_rows(df2_collapsed,        "log_cf",      "log_cf_energy", "CF_time", "CF_energy", "RoW-collapsed"),
  dual_rows(cf_paid_energy_df,    "log_cf_paid", "log_cf_energy", "CF_paid", "CF_energy", "EXIO-only"),
  dual_rows(cf_paid_energy_df38,  "log_cf_paid", "log_cf_energy", "CF_paid", "CF_energy", "RoW-collapsed")
)

dir.create("results", showWarnings = FALSE)
write_csv(tab, "results/table2_mediation_summary.csv")

fmt_p <- function(p) ifelse(is.na(p), "—", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
print(tab %>%
        transmute(variant, mediator, sample, n,
                  path_a = sprintf("%.2f (%s)", path_a, fmt_p(path_a_p)),
                  partial_r = ifelse(is.na(partial_r), "—", sprintf("%.2f (%s)", partial_r, fmt_p(partial_r_p))),
                  indirect = sprintf("%.3f", indirect), share = sprintf("%.1f%%", 100 * share),
                  sobel_p = fmt_p(sobel_p),
                  boot_ci = sprintf("[%.3f, %.3f]", boot_lo, boot_hi), boot_p = fmt_p(boot_p)),
      n = Inf, width = Inf)
cat("\nWrote results/table2_mediation_summary.csv (", nrow(tab), "rows )\n")


#### Table 3: income interaction — does the mediator's slope on protein change with income? ####
# log_protein ~ log_gdp_pcap * log_M, per CF variant and sample (the 9d/9n/9w specification).
# Reports the interaction coefficient (t-test p, bootstrap CI and p) and the marginal slope
# dY/dM = b_M + b_int * log(GDP) at two observed incomes (India, Netherlands), with
# delta-method SE and p.
# Own seed, so re-running this file leaves Table 2's bootstrap draws unchanged.

set.seed(20260922)
# Evaluate the marginal slope at observed incomes rather than round numbers:
# India ($7,187) is the poorest country in the EXIO-only sample, so $2,000 (used
# earlier) extrapolated well below the data; the Netherlands ($64,862) anchors the
# top. Both are inside the range of the RoW-collapsed sample too.
GDP_AT <- c(7187, 64862)

interaction_row <- function(d, m, variant, sample) {
  d <- d[complete.cases(d[c(X, m, Y)]), ]
  fml <- as.formula(paste(Y, "~", X, "*", m))
  fit <- lm(fml, d); s <- summary(fit)$coefficients; V <- vcov(fit)
  int <- paste0(X, ":", m)
  bt <- boot_ci(d, function(b) coef(lm(fml, b))[[int]])
  me <- map_dfr(GDP_AT, function(g) {
    x <- log(g); est <- s[m, "Estimate"] + s[int, "Estimate"] * x
    se <- sqrt(V[m, m] + x^2 * V[int, int] + 2 * x * V[m, int])
    tibble(gdp = g, me = est, me_se = se, me_p = 2 * pt(-abs(est / se), fit$df.residual))
  })
  tibble(variant, sample, n = nrow(d),
         interaction = s[int, "Estimate"], interaction_se = s[int, "Std. Error"], interaction_p = s[int, "Pr(>|t|)"],
         boot_lo = bt[["lo"]], boot_hi = bt[["hi"]], boot_p = bt[["p"]],
         gdp_low = GDP_AT[1], slope_low = me$me[1], slope_low_p = me$me_p[1],
         gdp_high = GDP_AT[2], slope_high = me$me[2], slope_high_p = me$me_p[2])
}

tab3 <- bind_rows(
  interaction_row(df,            "log_cf",        "CF_time",   "EXIO-only"),
  interaction_row(df_collapsed,  "log_cf",        "CF_time",   "RoW-collapsed"),
  interaction_row(cf_paid_df33,  "log_cf_paid",   "CF_paid",   "EXIO-only"),
  interaction_row(cf_paid_df38,  "log_cf_paid",   "CF_paid",   "RoW-collapsed"),
  interaction_row(df2,           "log_cf_energy", "CF_energy", "EXIO-only"),
  interaction_row(df2_collapsed, "log_cf_energy", "CF_energy", "RoW-collapsed")
)
write_csv(tab3, "results/table3_income_interaction.csv")
print(tab3 %>%
        transmute(variant, sample, n,
                  interaction = sprintf("%.3f (%s)", interaction, fmt_p(interaction_p)),
                  boot_ci = sprintf("[%.3f, %.3f]", boot_lo, boot_hi), boot_p = fmt_p(boot_p),
                  slope_low = sprintf("%.2f (%s)", slope_low, fmt_p(slope_low_p)),
                  slope_high = sprintf("%.2f (%s)", slope_high, fmt_p(slope_high_p))),
      n = Inf, width = Inf)
cat("\nWrote results/table3_income_interaction.csv (", nrow(tab3), "rows )\n")
