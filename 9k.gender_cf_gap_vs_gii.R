#### Does the gender CF gap track a country's Gender Inequality Index (GII), independent of GDP? ####
#
# 9c found no evidence that GDP moves the gendered distribution of
# provisioning time (GDP x gender interaction, p = 0.50-0.90). GII is a
# more direct, purpose-built measure of societal gender inequality (UNDP;
# reproductive health, empowerment, labor-market participation) than income
# is -- worth testing whether IT predicts the gap where income didn't,
# rather than concluding the gap is unrelated to gender inequality broadly.
#
# Gap defined as log(female CF) - log(male CF) per country (log-ratio,
# consistent with the log-log framing used throughout), using df_f/df_m
# from 9b.gender_split_income_control.R (EXIO-only, n=33).
#
# gii_lookup loading matches the precedent in 97.deprecate_reserve.r:
# data/gender_inequality_index_2020.csv, columns `code` (ISO3) and `gii`
# (UNDP Gender Inequality Index, higher = more unequal).
#
# Run AFTER 9b.gender_split_income_control.R in the same session (reuses
# `df_f`, `df_m`, `df`, `partial_cor()`).

gii_lookup <- read.csv("data/gender_inequality_index_2020.csv") %>%
  select(country = code, gii)

gap_df <- df_f %>%
  select(country, log_cf_f = log_cf) %>%
  inner_join(df_m %>% select(country, log_cf_m = log_cf), by = "country") %>%
  inner_join(df %>% select(country, gdp_pcap_ppp), by = "country") %>%
  inner_join(gii_lookup, by = "country") %>%
  mutate(gap = log_cf_f - log_cf_m,
         log_gdp_pcap = log(gdp_pcap_ppp))

cat(sprintf("\nGender CF gap with GII match: n = %d (vs. n = %d in the gender-split sample)\n",
            nrow(gap_df), nrow(df_f)))

r_gdp_gii <- cor(gap_df$log_gdp_pcap, gap_df$gii)
cat(sprintf("cor(log_gdp_pcap, gii) = %.3f\n", r_gdp_gii))

cat("\n---- Gap (log female CF - log male CF) vs. GII ----\n")
zero_order <- cor.test(gap_df$gap, gap_df$gii)
cat(sprintf("Zero-order r = %.3f (p = %.3g, n = %d)\n",
            zero_order$estimate, zero_order$p.value, nrow(gap_df)))

pc_gii <- partial_cor(gap_df$gap, gap_df$gii, gap_df$log_gdp_pcap)
cat(sprintf("Partial r | GDP per capita = %.3f (p = %.3g, n = %d)\n",
            pc_gii$r, pc_gii$p, pc_gii$n))

cat("\n---- Regression: gap ~ GDP + GII ----\n")
fit <- lm(gap ~ log_gdp_pcap + gii, data = gap_df)
print(summary(fit)$coefficients)

p_gap_gii <- ggplot(gap_df, aes(gii, gap)) +
  geom_point(alpha = 0.7) +
  ggrepel::geom_text_repel(aes(label = country), size = 3, max.overlaps = 20) +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
  labs(x = "UNDP Gender Inequality Index (higher = more unequal)",
       y = "log(female CF) - log(male CF)  [gender gap in provisioning time]",
       title = sprintf("Gender CF gap vs. Gender Inequality Index (n = %d)", nrow(gap_df))) +
  theme_minimal()
print(p_gap_gii)
ggsave("results/gender_cf_gap_vs_gii.pdf", p_gap_gii, width = 9, height = 6)

#### Zero-order relation with income, and how large the narrowing is ####
# The joint regression above gives GDP's coefficient NET of GII (p ~ 0.10); the
# manuscript needs the plain, zero-order relation too, reported like GII's.
# Size: gap is log(female CF / male CF), so exp(gap) is the female-to-male ratio
# of hours per gram. Fitted ratios at the 10th and 90th percentile of each
# predictor show how much the ratio moves across the sample.

cat("\n---- Gap vs. log GDP per capita, zero-order ----\n")
zero_order_gdp <- cor.test(gap_df$gap, gap_df$log_gdp_pcap)
cat(sprintf("Zero-order r = %.3f (p = %.3g, n = %d)\n",
            zero_order_gdp$estimate, zero_order_gdp$p.value, nrow(gap_df)))

size_line <- function(pred, label) {
  f <- lm(reformulate(pred, "gap"), data = gap_df)
  q <- quantile(gap_df[[pred]], c(0.1, 0.9))
  fit_q <- predict(f, newdata = setNames(data.frame(q), pred))
  cat(sprintf("%-14s slope = %.3f (p = %.3g) | fitted F/M ratio at P10 %s = %.2f, at P90 = %.2f\n",
              label, coef(f)[[pred]], summary(f)$coefficients[pred, "Pr(>|t|)"],
              label, exp(fit_q[1]), exp(fit_q[2])))
}
size_line("log_gdp_pcap", "log GDP")
size_line("gii", "GII")
cat(sprintf("Observed F/M ratio: median %.2f, range %.2f-%.2f (P10 %.2f, P90 %.2f)\n",
            median(exp(gap_df$gap)), min(exp(gap_df$gap)), max(exp(gap_df$gap)),
            quantile(exp(gap_df$gap), 0.1), quantile(exp(gap_df$gap), 0.9)))
print(gap_df %>% mutate(ratio = round(exp(gap), 2), gdp = round(gdp_pcap_ppp)) %>%
        select(country, gdp, gii, ratio) %>% arrange(gdp), n = Inf)

#### Which hours drive the narrowing? The F/M ratio by component ####
# The all-work gap above sums, per sex: unpaid household time (*_non.econ), paid
# domestic-origin hours (domestic_per_capita + preparation_econ = food service), and
# paid IMPORTED hours (import_per_capita), whose sex split is the exporting
# countries' workforce, not the consuming households'. Exports are excluded, as in
# tradeoff_protein_allwork_consump. Protein cancels in F/M, so each gap is simply
# log(women's hours / men's hours) for that component.

comp_types <- list(
  all_work        = c("preparation_non.econ", "processing_non.econ", "growth_collection_non.econ",
                      "energy_non.econ", "water_non.econ", "preparation_econ",
                      "domestic_per_capita", "import_per_capita"),
  unpaid          = c("preparation_non.econ", "processing_non.econ", "growth_collection_non.econ",
                      "energy_non.econ", "water_non.econ"),
  domestic_origin = c("preparation_non.econ", "processing_non.econ", "growth_collection_non.econ",
                      "energy_non.econ", "water_non.econ", "preparation_econ", "domestic_per_capita"),
  paid_domestic   = c("preparation_econ", "domestic_per_capita"),
  paid_imported   = c("import_per_capita")
)

hours_by_sex <- summary_food_df_long_with_ghd %>%
  filter(type %in% c("hr_f", "hr_m"), as.character(country) %in% gap_df$country) %>%
  mutate(country = as.character(country), footprint_type = as.character(footprint_type))

comp_gap <- map_dfr(names(comp_types), function(k) {
  hours_by_sex %>%
    filter(footprint_type %in% comp_types[[k]]) %>%
    group_by(country, type) %>% summarise(h = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = type, values_from = h) %>%
    transmute(country, component = k, hr_f, hr_m, gap_k = log(hr_f / hr_m))
}) %>%
  inner_join(gap_df %>% select(country, log_gdp_pcap, gii, gap), by = "country")

chk <- comp_gap %>% filter(component == "all_work")
cat(sprintf("\nSanity check: all_work gap reproduces the gap above -> max |diff| = %.2e\n",
            max(abs(chk$gap_k - chk$gap), na.rm = TRUE)))

cat("\n---- F/M ratio of hours, by component: relation with income and GII ----\n")
comp_summary <- comp_gap %>%
  filter(is.finite(gap_k)) %>%
  group_by(component) %>%
  group_modify(function(d, key) {
    ct_gdp <- cor.test(d$gap_k, d$log_gdp_pcap); ct_gii <- cor.test(d$gap_k, d$gii)
    f <- lm(gap_k ~ log_gdp_pcap, data = d); q <- quantile(d$log_gdp_pcap, c(0.1, 0.9))
    fr <- exp(predict(f, newdata = data.frame(log_gdp_pcap = q)))
    tibble(n = nrow(d), women_share_of_hours = median(d$hr_f / (d$hr_f + d$hr_m)),
           median_ratio = median(exp(d$gap_k)),
           r_gdp = unname(ct_gdp$estimate), p_gdp = ct_gdp$p.value,
           r_gii = unname(ct_gii$estimate), p_gii = ct_gii$p.value,
           fitted_ratio_P10_gdp = unname(fr[1]), fitted_ratio_P90_gdp = unname(fr[2]))
  }) %>% ungroup() %>%
  arrange(match(component, names(comp_types)))
print(comp_summary %>% mutate(across(where(is.numeric), ~ signif(., 3))), width = Inf)
write_csv(comp_summary, "results/gender_gap_by_component.csv")
cat("\nWrote results/gender_gap_by_component.csv\n")


#### Per-country F/M ratios by component, and the anchor countries used in the text ####
# The manuscript quotes observed countries rather than interpolated percentiles.
# India is the poorest country in this sample; the Netherlands is near the top of
# the income range (and the larger of the two countries bracketing the 90th percentile).

ratio_wide <- comp_gap %>%
  mutate(ratio = exp(gap_k)) %>%
  select(country, component, ratio, log_gdp_pcap) %>%
  pivot_wider(names_from = component, values_from = ratio) %>%
  mutate(gdp_pcap = round(exp(log_gdp_pcap))) %>%
  select(country, gdp_pcap, all_work, unpaid, domestic_origin, paid_domestic, paid_imported) %>%
  arrange(gdp_pcap)

cat("\n---- F/M ratio of hours by country and component (ascending income) ----\n")
print(ratio_wide %>% mutate(across(where(is.numeric), ~ round(., 2))), n = Inf)
write_csv(ratio_wide, "results/gender_gap_by_country_component.csv")

cat("\n---- Anchor countries ----\n")
print(ratio_wide %>% filter(country %in% c("IND", "NLD")) %>%
        mutate(across(where(is.numeric), ~ round(., 2))))
cat("\nWrote results/gender_gap_by_country_component.csv\n")


#### Population-weighted gender share: of ALL unpaid food hours in these countries, how much is women's? ####
# The medians above treat every country alike. Weighting each country's per-capita hours by
# its population gives the share of the hours actually worked across the sample -- the number
# the Abstract wants. pop_data_yr (iso3c, population) comes from 2.analyze_result.R.

pop_share <- comp_gap %>%
  inner_join(pop_data_yr, by = c("country" = "iso3c")) %>%
  mutate(H_f = hr_f * population, H_m = hr_m * population) %>%   # per-capita hours/day -> national hours/day
  group_by(component) %>%
  summarise(n = n(),
            pop_covered_bn = sum(population) / 1e9,
            women_share_pop_weighted = sum(H_f) / sum(H_f + H_m),
            ratio_pop_weighted = sum(H_f) / sum(H_m),
            women_share_median_country = median(hr_f / (hr_f + hr_m)),
            .groups = "drop") %>%
  arrange(match(component, names(comp_types)))

cat("\n---- Gender share of hours: population-weighted vs median country ----\n")
print(pop_share %>% mutate(across(where(is.numeric), ~ round(., 3))), width = Inf)
write_csv(pop_share, "results/gender_share_population_weighted.csv")

miss <- setdiff(gap_df$country, pop_data_yr$iso3c)
if (length(miss)) cat("Countries without population data (dropped):", paste(miss, collapse = ", "), "\n")
cat(sprintf("Population covered: %.2f bn across %d countries\n",
            pop_share$pop_covered_bn[1], pop_share$n[1]))
cat("\nWrote results/gender_share_population_weighted.csv\n")
