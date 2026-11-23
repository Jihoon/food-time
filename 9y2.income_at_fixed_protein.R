#### Income gradient of provisioning hours at the same protein supply ####
#
# The reverse of path b: among countries with the same protein supply, how much of
# the spread in hours per person does income explain?
#   log H_k ~ log GDP + log P,   k = all work, unpaid, paid (hours), energy (MJ)
# Because P is held fixed, the GDP coefficient is also the income elasticity of
# hours per gram of protein. Reported next to the unadjusted income slope
# (log H ~ log GDP), with the case-resampling bootstrap CI used elsewhere
# (9o, 9y) and the partial R2 of log GDP, t^2 / (t^2 + df).
# Same samples and H construction as 9y.paid_hours_protein_scaling.R:
# H = CF * P / 50, unpaid = all work - paid.

library(tidyverse)

PROJ <- "h:/MyDocuments/Projects/REMASS/food time"
setwd(PROJ)

if (!exists("cf_paid_df38") || !exists("df2_collapsed")) source("9x.table2_mediation_summary.R")

X <- "log_gdp_pcap"; Y <- "log_protein"
R_B <- if (exists("R_BOOT")) R_BOOT else 5000
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

gdp_row <- function(d, lh, what, sample) {
  d <- tibble(lh = lh, P = d[[Y]], G = d[[X]]) %>% filter(if_all(everything(), is.finite))
  fit <- lm(lh ~ G + P, d); s <- summary(fit)$coefficients; dfr <- fit$df.residual
  raw <- lm(lh ~ G, d); sr <- summary(raw)$coefficients
  b <- replicate(R_B, coef(lm(lh ~ G + P, d[sample.int(nrow(d), replace = TRUE), ]))[["G"]])
  tibble(what, sample, n = nrow(d),
         slope_fixed_P = s["G", 1], se = s["G", 2],
         ci_lo = s["G", 1] - qt(.975, dfr) * s["G", 2], ci_hi = s["G", 1] + qt(.975, dfr) * s["G", 2],
         p = s["G", 4],
         boot_lo = unname(quantile(b, .025, na.rm = TRUE)), boot_hi = unname(quantile(b, .975, na.rm = TRUE)),
         partial_R2 = s["G", 3]^2 / (s["G", 3]^2 + dfr),
         slope_unadjusted = sr["G", 1], p_unadjusted = sr["G", 4],
         protein_elasticity = s["P", 1])
}

rows <- function(d_all, d_paid, d_en, sample) {
  a <- d_all %>% select(country, log_cf, all_of(c(Y, X))) %>%
    inner_join(d_paid %>% select(country, log_cf_paid), by = "country") %>%
    filter(if_all(where(is.numeric), is.finite))
  P <- exp(a[[Y]]); h_all <- exp(a$log_cf) * P / 50; h_paid <- exp(a$log_cf_paid) * P / 50
  stopifnot(all(h_all - h_paid > 0))
  e <- d_en %>% select(country, log_cf_energy, all_of(c(Y, X))) %>% filter(if_all(where(is.numeric), is.finite))
  bind_rows(gdp_row(a, log(h_all), "all work", sample),
            gdp_row(a, log(h_all - h_paid), "unpaid", sample),
            gdp_row(a, log(h_paid), "paid", sample),
            gdp_row(e, e$log_cf_energy + e[[Y]] - log(50), "energy", sample))
}

set.seed(20261121)
tab <- bind_rows(rows(df, cf_paid_df33, df2, "EXIO-only"),
                 rows(df_collapsed, cf_paid_df38, df2_collapsed, "RoW-collapsed"))

cat("\n---- Income elasticity of hours (energy) per person at the same protein supply ----\n")
print(tab %>% transmute(sample, n, what,
        `GDP | protein` = sprintf("%.2f [%.2f, %.2f]", slope_fixed_P, ci_lo, ci_hi),
        p = fmt_p(p), `boot CI` = sprintf("[%.2f, %.2f]", boot_lo, boot_hi),
        `partial R2` = sprintf("%.2f", partial_R2),
        `GDP alone` = sprintf("%.2f (%s)", slope_unadjusted, fmt_p(p_unadjusted)),
        `protein el.` = sprintf("%.2f", protein_elasticity)),
      n = Inf, width = Inf)
cat("\nRead: 'GDP | protein' is the % change in hours per person (= per gram, protein fixed) for 1% more",
    "\n      income; 'partial R2' is the share of the remaining variation in log hours that income explains;",
    "\n      'GDP alone' is the slope without the protein control.\n")

dir.create("results", showWarnings = FALSE)
write_csv(tab, "results/income_at_fixed_protein.csv")
cat("Wrote results/income_at_fixed_protein.csv\n")
