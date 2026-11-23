# Domestic income slopes (domestic effort per 50 g of domestically produced protein).
# Run in the session where 2.analyze_result.R has been sourced. Writes results/domestic_income_slopes.csv.
suppressMessages(library(tidyverse))
need <- c("tradeoff_protein_totalecon_domestic_effort", "tradeoff_protein_foodecon_domestic_effort",
          "region_to_iso_1to1", "wdi_latest")
miss <- need[!sapply(need, exists)]
if (length(miss)) stop("Missing objects: ", paste(miss, collapse = ", "))

dom <- function(tab) {
  tab %>% filter(effort_label == "Domestic effort") %>%
    group_by(exio_region) %>%
    summarise(mj = first(mj_per_50g_protein), mj_sd = sd(mj_per_50g_protein),
              hr = sum(hr_per_50g_protein), .groups = "drop") %>%
    mutate(iso3c = unname(region_to_iso_1to1[as.character(exio_region)])) %>%
    filter(!is.na(iso3c), iso3c != "TWN", mj > 0, hr > 0) %>%
    inner_join(wdi_latest %>% select(iso3c, gdp_pcap_ppp), by = "iso3c")
}
fit <- function(d, y, label) {
  f <- lm(log(d[[y]]) ~ log(d$gdp_pcap_ppp)); s <- summary(f); ci <- confint(f)[2, ]
  tibble(series = label, n = nrow(d), slope = coef(f)[2], ci_lo = ci[1], ci_hi = ci[2],
         p = s$coefficients[2, 4], R2 = s$r.squared,
         IND = d[[y]][d$iso3c == "IND"][1], NLD = d[[y]][d$iso3c == "NLD"][1])
}
dT <- dom(tradeoff_protein_totalecon_domestic_effort)
dF <- dom(tradeoff_protein_foodecon_domestic_effort)
cat("max within-region sd of energy across gender rows (should be ~0):",
    max(c(dT$mj_sd, dF$mj_sd), na.rm = TRUE), "\n")
dom_slopes <- bind_rows(fit(dT, "mj", "domestic energy, food + non-food (MJ/50 g)"),
                        fit(dT, "hr", "domestic paid hours, food + non-food (hr/50 g)"),
                        fit(dF, "mj", "domestic energy, food only (MJ/50 g)"),
                        fit(dF, "hr", "domestic paid hours, food only (hr/50 g)"))
print(dom_slopes, width = Inf)
write.csv(dom_slopes, "results/domestic_income_slopes.csv", row.names = FALSE)
cat("Wrote results/domestic_income_slopes.csv\n")
