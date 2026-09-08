#### CF vs. GDP per capita, log-log: is there an achievable-frontier / elbow shape? ####
#
# Motivating question: tradeoff_convfac_protein_foodecon_domestic_effort_
# gender_range shows a dense low-income-elasticity cluster of mostly rich
# countries at low time-per-50g-protein, and a small set of high-CF outliers
# (IND, IDN, CHN, MEX, RoW Africa) well above them -- an apparent achievable
# frontier rather than a smooth continuum. Does the same shape show up
# plotting CF directly against GDP per capita? If so, that supports a
# normative "this CF level is already demonstrated as achievable, poorer
# countries should be able to reach it" argument, distinct from (and less
# fragile than) the weak GDP x CF statistical interaction (9d).
#
# Linear fit (constant elasticity, what the a-path regression already
# assumes) vs. loess (nonparametric, no assumption) overlaid -- if they
# diverge with the loess flattening at high income, that IS the elbow shape;
# if they track each other closely, the relationship really is closer to a
# constant elasticity throughout, and the frontier reading from the energy-
# time chart would need a different visual (e.g. that chart itself) to make
# the normative point, not this one.
#
# Run AFTER 9.capability_set_income_control.R in the same session (reuses
# `df`, `df_collapsed`).

p_cf_gdp <- function(data, label) {
  p <- ggplot(data, aes(log_gdp_pcap, log_cf)) +
    geom_point(alpha = 0.7) +
    ggrepel::geom_text_repel(aes(label = country), size = 3, max.overlaps = 20) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed") +
    geom_smooth(method = "loess", se = FALSE, color = "firebrick") +
    labs(x = "log(GDP per capita, PPP)", y = "log(CF, hr / 50g protein)",
         title = sprintf("CF vs. GDP per capita (%s, n = %d)", label, nrow(data)),
         subtitle = "Grey dashed = linear (constant elasticity); red = loess -- divergence at the high-income end is the elbow check") +
    theme_minimal()
  print(p)
  invisible(p)
}

p1 <- p_cf_gdp(df, "EXIO-only")
ggsave("results/cf_vs_gdp_exio_only.pdf", p1, width = 9, height = 7)

p2 <- p_cf_gdp(df_collapsed, "RoW-collapsed")
ggsave("results/cf_vs_gdp_row_collapsed.pdf", p2, width = 9, height = 7)
