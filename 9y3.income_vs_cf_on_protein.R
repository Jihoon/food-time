#### Income against time efficiency as explanations of protein supply ####
#
# How much of the cross-country variation in log protein supply does income
# explain, and what does each CF add once income is in the model?
#   M0: log P ~ log GDP
#   M1: log P ~ log GDP + log CF
# Reported: R2 of M0 and M1, the gain in R2 from adding CF, and the partial R2 of
# each term in M1 (t^2 / (t^2 + df)). The correlation of log CF with log GDP and
# the implied variance inflation factor show how much independent variation CF
# has left once income is taken out, i.e. how much power the CF term has.
# Same samples as Table 2 (9x). The partial correlation of CF in M1 is path b;
# its sign carries the shared-denominator bias described in 9y, so this
# comparison is descriptive and is not a test of mediation.

library(tidyverse)

PROJ <- "h:/MyDocuments/Projects/REMASS/food time"
setwd(PROJ)

if (!exists("cf_paid_df38") || !exists("df2_collapsed")) source("9x.table2_mediation_summary.R")

X <- "log_gdp_pcap"; Y <- "log_protein"
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

one_row <- function(d, m, variant, sample) {
  d <- d[complete.cases(d[c(X, m, Y)]), ]
  f0 <- lm(reformulate(X, Y), d)
  f1 <- lm(reformulate(c(X, m), Y), d)
  s0 <- summary(f0)$coefficients; s1 <- summary(f1)$coefficients; dfr <- f1$df.residual
  pr2 <- function(term) s1[term, "t value"]^2 / (s1[term, "t value"]^2 + dfr)
  r_xm <- cor(d[[X]], d[[m]])
  tibble(variant, sample, n = nrow(d),
         income_elasticity = s0[X, 1], p_income_alone = s0[X, 4],
         R2_income_alone = summary(f0)$r.squared,
         R2_income_plus_CF = summary(f1)$r.squared,
         gain_in_R2 = summary(f1)$r.squared - summary(f0)$r.squared,
         partial_R2_income = pr2(X), partial_R2_CF = pr2(m),
         p_CF = s1[m, 4], cor_CF_income = r_xm, VIF = 1 / (1 - r_xm^2),
         share_of_CF_variance_left = 1 - r_xm^2)
}

tab <- bind_rows(
  one_row(df,            "log_cf",        "CF_time",   "EXIO-only"),
  one_row(df_collapsed,  "log_cf",        "CF_time",   "RoW-collapsed"),
  one_row(cf_paid_df33,  "log_cf_paid",   "CF_paid",   "EXIO-only"),
  one_row(cf_paid_df38,  "log_cf_paid",   "CF_paid",   "RoW-collapsed"),
  one_row(df2,           "log_cf_energy", "CF_energy", "EXIO-only"),
  one_row(df2_collapsed, "log_cf_energy", "CF_energy", "RoW-collapsed")
)

cat("\n---- Protein supply: income alone vs income + CF ----\n")
print(tab %>% transmute(variant, sample, n,
        `income elasticity` = sprintf("%.2f (%s)", income_elasticity, fmt_p(p_income_alone)),
        `R2 income` = sprintf("%.2f", R2_income_alone),
        `R2 income+CF` = sprintf("%.2f", R2_income_plus_CF),
        `gain in R2` = sprintf("%.3f", gain_in_R2),
        `partial R2 income` = sprintf("%.2f", partial_R2_income),
        `partial R2 CF` = sprintf("%.3f", partial_R2_CF),
        `p (CF)` = fmt_p(p_CF),
        `r(CF, income)` = sprintf("%.2f", cor_CF_income),
        VIF = sprintf("%.1f", VIF),
        `CF variance left` = sprintf("%.0f%%", 100 * share_of_CF_variance_left)),
      n = Inf, width = Inf)
cat("\nRead: 'R2 income' is what income alone explains of log protein; 'gain in R2' is what CF adds;",
    "\n      'CF variance left' is the share of log CF not explained by income, the variation the CF",
    "\n      coefficient is estimated from.\n")

dir.create("results", showWarnings = FALSE)
write_csv(tab, "results/income_vs_cf_on_protein.csv")
cat("Wrote results/income_vs_cf_on_protein.csv\n")
