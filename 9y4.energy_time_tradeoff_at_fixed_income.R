#### Is the energy-time tradeoff there once income is held fixed? ####
#
# Log CF_time and log CF_energy correlate negatively across countries, but both
# track income in opposite directions (CF_time falls, CF_energy rises), so the
# correlation could be income alone. This asks whether time per gram of protein
# still falls where energy per gram rises, at the same income.
#   A: log CF_time ~ log CF_energy                       (raw)
#   B: log CF_time ~ log CF_energy + log GDP             (income held fixed)
#   C: log CF_time ~ log CF_energy + log GDP + log P     (income and protein fixed)
#   D: log H ~ log E + log GDP                           (per person, no shared denominator)
# CF_time and CF_energy share the protein denominator, which by itself pushes their
# correlation UP; a negative association here is therefore not produced by the ratio.
# Model D uses hours and energy per person (H = CF * P / 50), the check that does not
# divide by protein at all. Time is run as all-work (CF_time) and as paid only
# (CF_paid, the like-for-like scope with EXIOBASE energy). Same samples as Table 2.
# Case-resampling bootstrap, R = 5,000.

library(tidyverse)

PROJ <- "h:/MyDocuments/Projects/REMASS/food time"
setwd(PROJ)

if (!exists("cf_paid_df38") || !exists("df2_collapsed")) source("9x.table2_mediation_summary.R")

X <- "log_gdp_pcap"; Y <- "log_protein"
R_B <- if (exists("R_BOOT")) R_BOOT else 5000
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

build <- function(d_time, d_paid, d_en) {
  d_time %>% select(country, log_cf, all_of(c(X, Y))) %>%
    inner_join(d_paid %>% select(country, log_cf_paid), by = "country") %>%
    inner_join(d_en %>% select(country, log_cf_energy), by = "country") %>%
    filter(if_all(where(is.numeric), is.finite))
}

fit_row <- function(d, y, xterm, controls, model, time_kind, sample) {
  d$.y <- y; d$.x <- xterm
  f <- reformulate(c(".x", controls), ".y")
  m <- lm(f, d); s <- summary(m)$coefficients; dfr <- m$df.residual
  b <- replicate(R_B, coef(lm(f, d[sample.int(nrow(d), replace = TRUE), ]))[[".x"]])
  tibble(sample, time_kind, model, n = nrow(d),
         coef = s[".x", 1], ci_lo = s[".x", 1] - qt(.975, dfr) * s[".x", 2],
         ci_hi = s[".x", 1] + qt(.975, dfr) * s[".x", 2], p = s[".x", 4],
         boot_lo = unname(quantile(b, .025, na.rm = TRUE)), boot_hi = unname(quantile(b, .975, na.rm = TRUE)),
         partial_r = s[".x", 3] / sqrt(s[".x", 3]^2 + dfr))
}

run_sample <- function(d, sample) {
  P <- exp(d[[Y]])
  kinds <- list(all_work = list(lcf = d$log_cf),
                paid     = list(lcf = d$log_cf_paid))
  bind_rows(lapply(names(kinds), function(k) {
    lcf <- kinds[[k]]$lcf
    lH  <- lcf + d[[Y]] - log(50)                       # log hours per person
    lE  <- d$log_cf_energy + d[[Y]] - log(50)           # log MJ per person
    bind_rows(
      fit_row(d, lcf, d$log_cf_energy, NULL,        "A raw",               k, sample),
      fit_row(d, lcf, d$log_cf_energy, X,           "B | income",          k, sample),
      fit_row(d, lcf, d$log_cf_energy, c(X, Y),     "C | income, protein", k, sample),
      fit_row(d, lH,  lE,              X,           "D per person | income", k, sample))
  }))
}

set.seed(20261123)
tab <- bind_rows(
  run_sample(build(df,           cf_paid_df33, df2),           "EXIO-only"),
  run_sample(build(df_collapsed, cf_paid_df38, df2_collapsed), "RoW-collapsed"))

cat("\n---- Energy-time tradeoff: coefficient on log energy (CF_energy, or MJ per person in D) ----\n")
print(tab %>% transmute(sample, time = time_kind, model, n,
        coef = sprintf("%.2f [%.2f, %.2f]", coef, ci_lo, ci_hi), p = fmt_p(p),
        `boot CI` = sprintf("[%.2f, %.2f]", boot_lo, boot_hi), `partial r` = sprintf("%.2f", partial_r)),
      n = Inf, width = Inf)
cat("\nRead: a negative coefficient = where energy per gram (or per person) is higher, time is lower.",
    "\n      B vs A shows how much of the raw tradeoff is income; C and D remove the shared protein",
    "\n      denominator. A coefficient near zero in B-D means the raw correlation was income alone.\n")

dir.create("results", showWarnings = FALSE)
write_csv(tab, "results/energy_time_tradeoff_at_fixed_income.csv")
cat("Wrote results/energy_time_tradeoff_at_fixed_income.csv\n")
