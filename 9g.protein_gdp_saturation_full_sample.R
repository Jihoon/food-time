#### Big-n check: does protein supply itself saturate with GDP? ####
#
# Every CF-based test so far is capped at n=33 (EXIO-modeled) or n=38
# (RoW-collapsed) because CF needs Global Human Day time-use survey coverage.
# Protein supply alone doesn't need that -- summary_pro_df_long (built from
# FABIO, itself reconciled with FAO food balance sheets) covers the full
# ~187-country FABIO panel. This can't test whether CF mediates anything --
# there's no CF here -- but it CAN test, with far more power and many more
# low-income points, whether the underlying premise of the ceiling-effect
# story holds at all: does protein supply's response to income visibly
# flatten out at high GDP, using (almost) every country there is data for?
#
# Run AFTER 9.capability_set_income_control.R in the same session (reuses
# `wdi_latest`; needs `summary_pro_df_long` from 2.analyze_result.R).

protein_full <- summary_pro_df_long %>%
  filter(footprint_type %in% c("domestic_per_capita", "import_per_capita")) %>%
  mutate(country = as.character(country)) %>%
  group_by(country) %>%
  summarise(g_protein_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
  filter(g_protein_per_cap_day > 0, gdp_pcap_ppp > 0) %>%
  # BRN excluded: protein_per_cap_day ~0.0003 g/cap/day, a data artifact (not
  # a real value for any country) that sat as a high-leverage outlier at the
  # high-income end of the fit. Confirmed by inspecting the saved plot.
  filter(country != "BRN") %>%
  mutate(log_protein  = log(g_protein_per_cap_day),
         log_gdp_pcap = log(gdp_pcap_ppp))

cat(sprintf("\nFull FABIO/WDI overlap: n = %d countries (vs. n = 33 in every CF-based test so far)\n",
            nrow(protein_full)))

fit_linear <- lm(log_protein ~ log_gdp_pcap, data = protein_full)
fit_quad   <- lm(log_protein ~ poly(log_gdp_pcap, 2, raw = FALSE), data = protein_full)

cat("\n---- Linear: log(protein) ~ log(GDP/cap) ----\n")
print(summary(fit_linear)$coefficients)

cat("\n---- Quadratic: log(protein) ~ poly(log(GDP/cap), 2) ----\n")
print(summary(fit_quad)$coefficients)

cat("\n---- Does the quadratic term improve fit over linear? ----\n")
print(anova(fit_linear, fit_quad))

p_saturation <- ggplot(protein_full, aes(log_gdp_pcap, log_protein)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "firebrick", se = TRUE) +
  geom_smooth(method = "loess", color = "steelblue", se = FALSE, linetype = "dashed") +
  labs(x = "log(GDP per capita, PPP)", y = "log(Protein supply, g/cap/day)",
       title = sprintf("Does protein supply saturate with income? (n = %d, full FABIO/WDI overlap)",
                        nrow(protein_full)),
       subtitle = "Red = quadratic fit with 95% CI; blue dashed = loess (nonparametric) — agreement between them is the check") +
  theme_minimal()
print(p_saturation)
ggsave("results/protein_gdp_saturation_full_sample.pdf", p_saturation, width = 9, height = 6)
