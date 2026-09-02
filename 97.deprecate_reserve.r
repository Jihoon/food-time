
#### Protein g per hour: population-weighted inequality (Gini & Atkinson) ####

# Helper: population-weighted Gini coefficient
weighted_gini <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0 & x >= 0
  x <- x[ok]; w <- w[ok]
  n <- length(x)
  if (n == 0) return(NA_real_)
  ord  <- order(x)
  x    <- x[ord]; w <- w[ord]
  N    <- sum(w)
  cumw <- cumsum(w)
  L    <- cumsum(w * x) / sum(w * x)
  F_lo <- c(0, cumw[-n]) / N
  F_hi <- cumw / N
  area <- sum((F_hi - F_lo) * (c(0, L[-n]) + L) / 2)
  1 - 2 * area
}

# Helper: population-weighted Atkinson index
# epsilon = 0.5 (moderate aversion), 1.0 (log-utilitarian)
weighted_atkinson <- function(x, w, epsilon = 0.5) {
  ok <- !is.na(x) & !is.na(w) & w > 0 & x > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) == 0) return(NA_real_)
  N  <- sum(w)
  mu <- sum(w * x) / N
  ede <- if (epsilon == 1) {
    exp(sum(w * log(x)) / N)
  } else {
    (sum(w * x^(1 - epsilon)) / N)^(1 / (1 - epsilon))
  }
  1 - ede / mu
}

pop_data_ineq <- subset(countrypops, year == yr) %>%
  select(country = country_code_3, population)

# Derive protein_g_per_hr for GHD countries: sum all domestic time components
# (economic + non-economic) per country × gender, then invert to get g protein per hour.
# Each gender covers ~half the national population, so weight each row by population / 2.
protein_g_per_hr_df <- summary_time_protein %>%
  filter(country %in% cty_ghd, cat == "domestic") %>%
  group_by(country, type) %>%
  summarise(hr_per_50g_protein = sum(hr_per_50g_protein, na.rm = TRUE), .groups = "drop") %>%
  mutate(protein_g_per_hr = 50 / hr_per_50g_protein) %>%
  left_join(pop_data_ineq, by = "country") %>%
  mutate(weight = population / 2) %>%
  filter(!is.na(weight), is.finite(protein_g_per_hr), protein_g_per_hr > 0)

ineq_protein_df <- protein_g_per_hr_df %>%
  summarise(
    n_obs         = n(),
    n_countries   = n_distinct(country),
    mean_g_per_hr = weighted.mean(protein_g_per_hr, weight),
    gini          = weighted_gini(protein_g_per_hr, weight),
    atkinson_0.5  = weighted_atkinson(protein_g_per_hr, weight, epsilon = 0.5),
    atkinson_1.0  = weighted_atkinson(protein_g_per_hr, weight, epsilon = 1.0),
    atkinson_2.0  = weighted_atkinson(protein_g_per_hr, weight, epsilon = 2.0)
  )

print(ineq_protein_df)

## Global average protein g per hour: population-weighted mean of protein_g_per_hr across GHD countries
global_avg_protein_g_per_hr <- weighted.mean(protein_g_per_hr_df$protein_g_per_hr, protein_g_per_hr_df$weight)
print(global_avg_protein_g_per_hr) # = 89.53 g/hr
perfect_equality <-global_avg_protein_g_per_hr * (1-ineq_protein_df$atkinson_1.0) # = 71.62 g/hr at epsilon=1.0, which is the "equally distributed equivalent" level of protein g/hr if we had perfect equality across countries.






#### Attempts to incorporate additional country-level indicators into the tradeoff plots ####

# e.g. gender equity, agrifood-system GHG emissions, and agricultural energy-use GHG emissions.


# UNDP Gender Inequality Index (2020, via Our World in Data), inverted to a "gender
# equity index" (1 - GII) so a bigger bubble reads as more gender-equitable, not less.
gii_lookup = read.csv("data/gender_inequality_index_2020.csv") %>%
  transmute(country = code, gender_equity_index = 1 - gii) %>%
  filter(country %in% regions$iso3c)

# FAOSTAT agrifood-system total GHG emissions (2020, kt CO2eq AR5), converted to
# kg CO2eq per capita using the same population source as the rest of the pipeline.
ghg_agrifood_raw = read.csv("data/FAOSTAT_agrifood_ghg_totals_2020.csv", check.names = FALSE)
ghg_agrifood_lookup = ghg_agrifood_raw %>%
  mutate(`Area Code (M49)` = as.numeric(gsub("'", "", `Area Code (M49)`)),
         iso3c = countrycode::countrycode(`Area Code (M49)`,
                                          origin = "un", destination = "iso3c",
                                          warn = FALSE)) %>%
  filter(!is.na(iso3c)) %>%
  left_join(subset(countrypops, year == yr) %>% select(iso3c = country_code_3, population),
            by = "iso3c") %>%
  drop_na(population) %>%
  mutate(ghg_kg_co2eq_per_cap = Value * 1e6 / population) %>%  # kt -> kg
  select(country = iso3c, ghg_kg_co2eq_per_cap)

# FAO's total is production-based (domestic production = what's consumed domestically
# + what's exported; imports are excluded since they're another country's production).
# Split it into domestic vs. export portions using each country's domestic/export kcal
# shares as a proxy, assuming uniform emissions intensity per kcal produced.
kcal_domestic_export_share = summary_kcal_df_long %>%
  filter(footprint_type %in% c("domestic_per_capita", "export_per_capita")) %>%
  mutate(country = as.character(country)) %>%
  select(country, footprint_type, kcal_per_cap_day = per_capita_value) %>%
  pivot_wider(names_from = footprint_type, values_from = kcal_per_cap_day) %>%
  mutate(domestic_share = domestic_per_capita / (domestic_per_capita + export_per_capita),
         export_share   = export_per_capita   / (domestic_per_capita + export_per_capita))

ghg_agrifood_lookup = ghg_agrifood_lookup %>%
  left_join(kcal_domestic_export_share %>% select(country, domestic_share, export_share),
            by = "country") %>%
  mutate(ghg_domestic_kg_per_cap = ghg_kg_co2eq_per_cap * domestic_share,
         ghg_export_kg_per_cap   = ghg_kg_co2eq_per_cap * export_share)

# Domestic vs. export split, per capita, non-RoW countries only — reusable so the
# same treatment can be applied to any total/domestic/export triple of columns.
make_ghg_split_plot = function(lookup_df, domestic_col, export_col, total_col, y_lab, title) {
  ord = lookup_df %>%
    filter(!country %in% row_countries) %>%
    arrange(-.data[[total_col]]) %>%
    pull(country)

  df_long = lookup_df %>%
    filter(!country %in% row_countries) %>%
    select(country, all_of(c(domestic_col, export_col))) %>%
    pivot_longer(cols = c(domestic_col, export_col), names_to = "flow", values_to = "value") %>%
    mutate(flow = ifelse(flow == domestic_col, "Domestic", "Export"),
           flow = factor(flow, levels = c("Export", "Domestic")),
           country = factor(country, levels = ord)) %>%
    drop_na()

  ggplot(df_long, aes(x = country, y = value, fill = flow)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Domestic" = "#1f77b4", "Export" = "#ff7f0e")) +
    labs(x = "Country (ISO3)", y = y_lab, fill = "", title = title) +
    theme_minimal() +
    theme(legend.position = "top", axis.text.x = element_text(angle = 90, hjust = 1))
}

p_ghg_split_nonrow = make_ghg_split_plot(
  ghg_agrifood_lookup, "ghg_domestic_kg_per_cap", "ghg_export_kg_per_cap", "ghg_kg_co2eq_per_cap",
  "kg CO2eq / cap / yr",
  paste0("Agrifood-system total GHG emissions per capita (", year, ") — non-RoW"))

# FAOSTAT "Emissions from Energy use in agriculture" (on-farm energy use only — a
# narrower scope than the "Energy" line item inside Emissions Totals, which covers
# the whole agrifood supply chain). No CO2eq element is provided directly, so it's
# built from the raw CH4/CO2/N2O elements using IPCC AR5 GWP100 (no climate-carbon
# feedback) factors — the same convention FAOSTAT itself uses for its CO2eq(AR5) figures.
energy_ag_raw = read.csv("data/FAOSTAT_energy_in_agriculture_totals_2020.csv", check.names = FALSE)
energy_ag_lookup = energy_ag_raw %>%
  mutate(`Area Code (M49)` = as.numeric(gsub("'", "", `Area Code (M49)`)),
         iso3c = countrycode::countrycode(`Area Code (M49)`,
                                          origin = "un", destination = "iso3c",
                                          warn = FALSE)) %>%
  filter(!is.na(iso3c)) %>%
  select(iso3c, Element, Value) %>%
  pivot_wider(names_from = Element, values_from = Value, values_fn = sum) %>%
  mutate(co2eq_kt = `Emissions (CO2)` * 1 + `Emissions (CH4)` * 28 + `Emissions (N2O)` * 265) %>%
  left_join(subset(countrypops, year == yr) %>% select(iso3c = country_code_3, population),
            by = "iso3c") %>%
  drop_na(population) %>%
  mutate(energy_ag_kg_co2eq_per_cap = co2eq_kt * 1e6 / population) %>%
  select(country = iso3c, energy_ag_kg_co2eq_per_cap) %>%
  left_join(kcal_domestic_export_share %>% select(country, domestic_share, export_share), by = "country") %>%
  mutate(energy_ag_domestic_kg_per_cap = energy_ag_kg_co2eq_per_cap * domestic_share,
         energy_ag_export_kg_per_cap   = energy_ag_kg_co2eq_per_cap * export_share)

p_energy_ag_split_nonrow = make_ghg_split_plot(
  energy_ag_lookup, "energy_ag_domestic_kg_per_cap", "energy_ag_export_kg_per_cap", "energy_ag_kg_co2eq_per_cap",
  "kg CO2eq / cap / yr",
  paste0("Agricultural energy-use GHG emissions per capita (", year, ") — non-RoW"))

p_ghg_energy_combined = (p_ghg_split_nonrow | p_energy_ag_split_nonrow) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
print(p_ghg_energy_combined)

ggsave("results/ghg_and_energy_ag_domestic_vs_export_nonrow.pdf", p_ghg_energy_combined, width = 30, height = 9)

combined_pcap = bind_rows(
  tradeoff_pcap_econlabor %>% mutate(scope = "Economic"),
  tradeoff_pcap_allwork %>% filter(!is_row) %>% mutate(scope = "Economic + non-economic")
) %>%
  inner_join(gii_lookup, by = "country")

combined_protein50g = bind_rows(
  tradeoff_protein_econlabor %>% mutate(scope = "Economic"),
  tradeoff_protein_allwork %>% filter(!is_row) %>% mutate(scope = "Economic + non-economic")
) %>%
  inner_join(gii_lookup, by = "country")

p_left_overlay = overlay_scatter(
  combined_pcap, "mj_per_cap_day", "hr_per_cap_day", "gender_equity_index",
  "Energy (MJ/cap/day)", "Time (hr/cap/day)", "Gender equity index (1 - GII)",
  paste0("Energy vs. time per capita (", year, ") — Economic vs. Economic + non-economic, non-RoW"))

p_right_overlay = overlay_scatter(
  combined_protein50g, "mj_per_50g_protein", "hr_per_50g_protein", "gender_equity_index",
  "Energy (MJ / 50 g protein)", "Time (hr / 50 g protein)", "Gender equity index (1 - GII)",
  paste0("Energy vs. time to provision 50 g protein (", year, ") — Economic vs. Economic + non-economic, non-RoW"))

p_tradeoff_protein_overlay = (p_left_overlay | p_right_overlay) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")
print(p_tradeoff_protein_overlay)

ggsave("results/tradeoff_protein_pcap_vs_50g_nonrow_overlay.pdf", p_tradeoff_protein_overlay, width = 22, height = 8)
