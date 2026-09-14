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
