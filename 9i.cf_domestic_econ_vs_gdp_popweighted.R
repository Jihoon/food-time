#### CF (domestic, paid-labor-only) vs. GDP, population-weighted: does the achievable-frontier shape show up here? ####
#
# Rebuild of 9h using the narrower domestic-economic-only CF from 9e
# (food-sector econ. labor, domestic effort) instead of the all-work CF --
# this is the measure that actually matches
# tradeoff_convfac_protein_foodecon_domestic_effort_gender_range's
# frontier-looking cluster, and it has a much steeper GDP elasticity (-1.16
# vs. -0.52) than the all-work version, so it's the better candidate for
# showing the same shape in GDP-space.
#
# Population-weighted (both the loess/linear fit AND point size) so a single
# small-population outlier like LUX -- which produced a spurious upward hook
# in 9h, identical in both country samples, the signature of a one-point
# artifact rather than a real pattern -- doesn't dominate the fit the way it
# does unweighted. For an "is this level broadly achievable" argument, a
# large, populous low-CF country should count for more than a wealthy
# microstate. Both weighted and unweighted versions are plotted so the
# effect of weighting itself is visible, not assumed.
#
# Run AFTER 9e.domestic_econ_cf_mediation_test.R in the same session
# (reuses `df_domestic_econ`; needs `pop_data_yr` from 2.analyze_result.R
# for population weights).

df_domestic_econ_pop <- df_domestic_econ %>%
  inner_join(pop_data_yr %>% select(country = iso3c, population), by = "country")

cat(sprintf("\nDomestic-economic CF with population data: n = %d (vs. n = %d without population match)\n",
            nrow(df_domestic_econ_pop), nrow(df_domestic_econ)))

p_weighted <- ggplot(df_domestic_econ_pop, aes(log_gdp_pcap, log_cf, weight = population)) +
  geom_point(aes(size = population), alpha = 0.6) +
  ggrepel::geom_text_repel(aes(label = country), size = 3, max.overlaps = 20) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, color = "firebrick") +
  scale_size_continuous(name = "Population", labels = scales::comma) +
  labs(x = "log(GDP per capita, PPP)",
       y = "log(CF, hr / 50g protein) — domestic, paid labor only",
       title = sprintf("Domestic economic-only CF vs. GDP per capita, population-weighted (n = %d)",
                        nrow(df_domestic_econ_pop)),
       subtitle = "Point size and fit weighting = population. Grey dashed = linear; red = loess") +
  theme_minimal()
print(p_weighted)
ggsave("results/cf_domestic_econ_vs_gdp_popweighted.pdf", p_weighted, width = 10, height = 7)

p_unweighted <- ggplot(df_domestic_econ_pop, aes(log_gdp_pcap, log_cf)) +
  geom_point(alpha = 0.7) +
  ggrepel::geom_text_repel(aes(label = country), size = 3, max.overlaps = 20) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, color = "steelblue") +
  labs(x = "log(GDP per capita, PPP)",
       y = "log(CF, hr / 50g protein) — domestic, paid labor only",
       title = sprintf("Same, UNWEIGHTED for comparison (n = %d)", nrow(df_domestic_econ_pop)),
       subtitle = "Grey dashed = linear; blue = loess") +
  theme_minimal()
print(p_unweighted)
ggsave("results/cf_domestic_econ_vs_gdp_unweighted.pdf", p_unweighted, width = 10, height = 7)
