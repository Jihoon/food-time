#### CF vs. GDP, both expressed per 50 g of domestically-consumed protein ####
#
# 9i put domestic-economic CF against GDP PER CAPITA -- an income level, a
# different kind of quantity than CF or the energy-time chart's own axes,
# both of which are embodied inputs PER 50g of domestically-consumed
# protein (see that chart's own subtitle). This rebuilds the x-axis to
# match: GDP per capita divided by domestic protein supply, same
# denominator convention as CF's own construction, so money, time, and
# energy are all on the same "embodied input per unit of protein" basis
# and directly comparable to tradeoff_convfac_protein_foodecon_domestic_
# effort_gender_range.
#
# Uses `domestic_protein_df` (built in 9e -- country-level
# g_protein_per_cap_day_domestic, the same protein denominator
# domestic-economic CF itself uses) rather than df_domestic_econ's
# g_protein_per_cap_day, which is TOTAL consumption protein (kept there
# deliberately for comparability with the earlier mediation tests, but the
# wrong denominator for this particular comparison).
#
# Population-weighted as in 9i, same reasoning (LUX-style microstate
# leverage).
#
# Run AFTER 9i.cf_domestic_econ_vs_gdp_popweighted.R in the same session
# (reuses `domestic_econ_cf_df`, `domestic_protein_df`, `wdi_latest`,
# `pop_data_yr`).

df_gdp_per_protein <- domestic_econ_cf_df %>%
  select(country, hr_per_50g_protein_domestic_econ) %>%
  inner_join(domestic_protein_df, by = "country") %>%
  inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
  inner_join(pop_data_yr %>% select(country = iso3c, population), by = "country") %>%
  filter(gdp_pcap_ppp > 0, g_protein_per_cap_day_domestic > 0) %>%
  mutate(gdp_per_50g_protein_domestic = gdp_pcap_ppp / g_protein_per_cap_day_domestic * 50,
         log_cf  = log(hr_per_50g_protein_domestic_econ),
         log_gdp_per_protein = log(gdp_per_50g_protein_domestic))

cat(sprintf("\nCF vs. GDP-per-50g-protein, domestic-only basis: n = %d\n", nrow(df_gdp_per_protein)))

p_gdp_per_protein <- ggplot(df_gdp_per_protein, aes(log_gdp_per_protein, log_cf, weight = population)) +
  geom_point(aes(size = population), alpha = 0.6) +
  ggrepel::geom_text_repel(aes(label = country), size = 3, max.overlaps = 20) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, color = "firebrick") +
  scale_size_continuous(name = "Population", labels = scales::comma) +
  labs(x = "log(GDP per capita / 50g domestic protein, PPP $)",
       y = "log(CF, hr / 50g protein) — domestic, paid labor only",
       title = sprintf("CF vs. GDP, both per 50g domestically-consumed protein, population-weighted (n = %d)",
                        nrow(df_gdp_per_protein)),
       subtitle = "Same basis as the energy-time tradeoff chart. Point size/fit weighting = population.") +
  theme_minimal()
print(p_gdp_per_protein)
ggsave("results/cf_vs_gdp_per_protein_popweighted.pdf", p_gdp_per_protein, width = 10, height = 7)
