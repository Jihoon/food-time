#### Is the domestic/import CF gap much bigger for time than for energy? ####
#### Consistent with wage arbitrage adding to the labor-time gap beyond ####
#### whatever pure technology/resource differences show up in energy too? ####
#
# Energy is a globally-traded commodity (prices converge across borders,
# net of transport/tariffs); labor is not (migration is restricted, wages
# are not arbitraged the way commodity prices are). If the domestic/import
# CF ratio is much larger for time than for energy, that's consistent with
# something beyond pure technology/resource differences driving the time
# gap specifically -- the wage-differential mechanism already discussed
# (International Embodied Labor and Trade) but without wage data of our
# own. NOT proof of that mechanism: labor and energy sit on different
# production margins (mechanization trades one for the other), so part of
# any gap-size difference could reflect substitution elasticities rather
# than price arbitrage. Report both ratios plainly and let the reader see
# the caveat, don't oversell it.
#
# df_mosaic_9q / df_mosaic_energy_9q: one row per country x footprint_type
# (domestic_per_capita / import_per_capita) x sector (food/non-food) x
# effort_origin (domestic/import of the labor or energy itself).
# g_protein_per_cap_day is constant within each footprint_type block (the
# protein quantity in that consumption stream), so summing hr_per_cap_day
# (or gj_per_cap_day) across sector x effort_origin within a footprint_type
# and dividing by that stream's protein quantity gives the CF of
# domestically-consumed vs. imported protein as a whole -- the level this
# question is actually asked at, not the domestic/import-of-the-labor axis
# the footnote in Findings/Discussion is about (a different split).
#
# Run AFTER whatever builds df_mosaic_9q and df_mosaic_energy_9q in
# 2.analyze_result.R (both already in the R environment once that's run).

library(tidyverse)

stream_cf <- function(df, effort_col, unit_scale = 1) {
  df %>%
    group_by(country, footprint_type) %>%
    summarise(effort_stream = sum(.data[[effort_col]], na.rm = TRUE) * unit_scale,
              g_protein_stream = first(g_protein_per_cap_day),
              .groups = "drop") %>%
    filter(g_protein_stream > 0) %>%
    mutate(cf_stream = effort_stream / g_protein_stream * 50) %>%
    select(country, footprint_type, cf_stream) %>%
    pivot_wider(names_from = footprint_type, values_from = cf_stream,
                names_prefix = "cf_")
}


pro_mosaic_9q = summary_pro_df_long %>%
  filter(footprint_type %in% c("domestic_per_capita", "import_per_capita")) %>%
  mutate(country = as.character(country)) %>%
  select(country, footprint_type, g_protein_per_cap_day = per_capita_value)

pro_width_9q = pro_mosaic_9q %>%
  mutate(source = ifelse(footprint_type == "domestic_per_capita", "Domestic", "Import"),
         source = factor(source, levels = c("Domestic", "Import"))) %>%
  arrange(country, source) %>%
  group_by(country) %>%
  mutate(xmax = cumsum(g_protein_per_cap_day),
         xmin = xmax - g_protein_per_cap_day) %>%
  ungroup() %>%
  select(country, footprint_type, source, xmin, xmax, g_protein_per_cap_day)

hr_mosaic_sector_9q = effort_consumption_df %>%
  mutate(country = region_to_iso_1to1[as.character(exio_region)]) %>%
  filter(type %in% c("hr_m", "hr_f")) %>%
  mutate(footprint_type = ifelse(protein_source == "domestic", "domestic_per_capita", "import_per_capita")) %>%
  group_by(country, footprint_type, sector, effort_origin) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

df_mosaic_9q = hr_mosaic_sector_9q %>%
  left_join(pro_width_9q, by = c("country", "footprint_type")) %>%
  drop_na() %>%
  mutate(min_per_g = hr_per_cap_day * 60 / g_protein_per_cap_day,
         sector = factor(sector, levels = c("food", "non-food")),
         effort_label = ifelse(effort_origin == "domestic", "Domestic effort", "Import effort"),
         fill_grp = factor(paste(sector, effort_label, sep = " | "), levels = mosaic_fill_levels)) %>%
  arrange(country, source, sector, effort_label) %>%
  group_by(country, source) %>%
  mutate(ymax = cumsum(min_per_g),
         ymin = ymax - min_per_g) %>%
  ungroup() %>%
  mutate(min_per_cap_day = hr_per_cap_day * 60,
         label_min = ifelse(round(min_per_cap_day) == 0, NA_character_,
                            paste0(round(min_per_cap_day), " min")),
         fits_inside = (xmax - xmin) >= 0.05 * max(xmax, na.rm = TRUE) &
                       (ymax - ymin) >= 0.06 * max(ymax, na.rm = TRUE),
         label_y = ifelse(fits_inside, (ymin + ymax) / 2, ymax + 0.03 * max(ymax, na.rm = TRUE)))


en_mosaic_sector_9q = effort_consumption_df %>%
  mutate(country = region_to_iso_1to1[as.character(exio_region)]) %>%
  filter(type == "en") %>%
  mutate(footprint_type = ifelse(protein_source == "domestic", "domestic_per_capita", "import_per_capita"),
         gj_per_cap_day = per_capita_value / 365 / 1000) %>%
  select(country, footprint_type, sector, effort_origin, gj_per_cap_day)

df_mosaic_energy_9q = en_mosaic_sector_9q %>%
  left_join(pro_width_9q, by = c("country", "footprint_type")) %>%
  drop_na() %>%
  mutate(mj_per_g = gj_per_cap_day * 1000 / g_protein_per_cap_day,
         sector = factor(sector, levels = c("food", "non-food")),
         effort_label = ifelse(effort_origin == "domestic", "Domestic effort", "Import effort"),
         fill_grp = factor(paste(sector, effort_label, sep = " | "), levels = mosaic_fill_levels)) %>%
  arrange(country, source, sector, effort_label) %>%
  group_by(country, source) %>%
  mutate(ymax = cumsum(mj_per_g),
         ymin = ymax - mj_per_g) %>%
  ungroup() %>%
  mutate(gj_per_year = gj_per_cap_day * 365,
         label_gj = ifelse(round(gj_per_year, 1) == 0, NA_character_,
                            paste0(round(gj_per_year, 1), " GJ/yr")),
         # If a rectangle is too small (in either dimension) for the label to plausibly
         # fit inside it, place the label just above the rectangle instead.
         fits_inside = (xmax - xmin) >= 0.05 * max(xmax, na.rm = TRUE) &
                       (ymax - ymin) >= 0.06 * max(ymax, na.rm = TRUE),
         label_y = ifelse(fits_inside, (ymin + ymax) / 2, ymax + 0.03 * max(ymax, na.rm = TRUE)))


# time: hr_per_cap_day is already in hours -> CF in hr/50g, no rescale
cf_time <- stream_cf(df_mosaic_9q, "hr_per_cap_day") %>%
  rename(cf_time_domestic = cf_domestic_per_capita, cf_time_import = cf_import_per_capita) %>%
  mutate(ratio_time = cf_time_import / cf_time_domestic)

# energy: gj_per_cap_day -> MJ (x1000) -> CF in MJ/50g, matching mj_per_50g_protein elsewhere
cf_energy <- stream_cf(df_mosaic_energy_9q, "gj_per_cap_day", unit_scale = 1000) %>%
  rename(cf_energy_domestic = cf_domestic_per_capita, cf_energy_import = cf_import_per_capita) %>%
  mutate(ratio_energy = cf_energy_import / cf_energy_domestic)

comparison <- cf_time %>%
  select(country, cf_time_domestic, cf_time_import, ratio_time) %>%
  inner_join(cf_energy %>% select(country, cf_energy_domestic, cf_energy_import, ratio_energy),
             by = "country") %>%
  mutate(ratio_of_ratios = ratio_time / ratio_energy)

cat("\n---- Domestic/import CF, time vs. energy, per country ----\n")
print(comparison %>% mutate(across(where(is.numeric), ~round(., 3))), n = Inf)

# Countries with negligible import protein give an unstable/extreme ratio
# (dividing by ~0) -- exclude those explicitly rather than let them
# silently dominate a summary statistic; report how many were dropped.
imp_share <- df_mosaic_9q %>%
  distinct(country, footprint_type, g_protein_per_cap_day) %>%
  pivot_wider(names_from = footprint_type, values_from = g_protein_per_cap_day,
              names_prefix = "g_") %>%
  mutate(import_share = g_import_per_capita / (g_domestic_per_capita + g_import_per_capita))

meaningful <- imp_share %>% filter(import_share > 0.02) %>% pull(country)
comp_filtered <- comparison %>% filter(country %in% meaningful)

cat(sprintf("\n%d of %d countries have >2%% of protein supply from imports; ratio summary uses those.\n",
            length(meaningful), nrow(comparison)))
cat(sprintf("\nMedian domestic/import CF ratio -- time: %.2f | energy: %.2f\n",
            median(comp_filtered$ratio_time, na.rm = TRUE),
            median(comp_filtered$ratio_energy, na.rm = TRUE)))
cat(sprintf("Mean domestic/import CF ratio   -- time: %.2f | energy: %.2f\n",
            mean(comp_filtered$ratio_time, na.rm = TRUE),
            mean(comp_filtered$ratio_energy, na.rm = TRUE)))
cat(sprintf("Median ratio-of-ratios (time ratio / energy ratio): %.2f\n",
            median(comp_filtered$ratio_of_ratios, na.rm = TRUE)))

cat("\n#### What to look at ####\n")
cat("If ratio_time is consistently and substantially larger than ratio_energy\n")
cat("across most countries (ratio_of_ratios >> 1), that's the pattern worth\n")
cat("writing up -- paste the full per-country table and the summary lines back.\n")

#### Does the ratio's DIRECTION (above/below 1) track income? ####
# Rich countries importing at a *higher* CF than their own domestic
# production, and poor countries importing at a *lower* CF, is a distinct
# claim from "the time ratio is bigger than the energy ratio" above -- a
# systematic reversal by income, not just a magnitude difference. Test it
# directly against log GDP per capita (reuses wdi_latest from script 9)
# rather than eyeballing a handful of countries.

comp_gdp <- comp_filtered %>%
  inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
  mutate(log_gdp_pcap = log(gdp_pcap_ppp),
         log_ratio_time = log(ratio_time),
         log_ratio_energy = log(ratio_energy))

cat("\n---- Sorted by income: does ratio_time / ratio_energy cross 1 in order? ----\n")
print(comp_gdp %>%
        arrange(gdp_pcap_ppp) %>%
        select(country, gdp_pcap_ppp, ratio_time, ratio_energy) %>%
        mutate(across(c(gdp_pcap_ppp, ratio_time, ratio_energy), ~round(., 2))),
      n = Inf)

fit_time   <- lm(log_ratio_time ~ log_gdp_pcap, data = comp_gdp)
fit_energy <- lm(log_ratio_energy ~ log_gdp_pcap, data = comp_gdp)

cat("\n---- Does log(import/domestic CF ratio) rise with log GDP per capita? ----\n")
cat("(a positive, significant slope here is the formal version of the\n")
cat(" rich-import-high / poor-import-low pattern)\n")
print(round(summary(fit_time)$coefficients, 4))
cat(sprintf("time:   slope = %.4f, p = %.4g, n = %d\n",
            coef(fit_time)["log_gdp_pcap"],
            summary(fit_time)$coefficients["log_gdp_pcap", "Pr(>|t|)"],
            nrow(comp_gdp)))
print(round(summary(fit_energy)$coefficients, 4))
cat(sprintf("energy: slope = %.4f, p = %.4g, n = %d\n",
            coef(fit_energy)["log_gdp_pcap"],
            summary(fit_energy)$coefficients["log_gdp_pcap", "Pr(>|t|)"],
            nrow(comp_gdp)))

n_time_above1   <- sum(comp_gdp$ratio_time > 1, na.rm = TRUE)
n_energy_above1 <- sum(comp_gdp$ratio_energy > 1, na.rm = TRUE)
cat(sprintf("\n%d/%d countries have ratio_time > 1 (imports less time-efficient than domestic); %d/%d for energy.\n",
            n_time_above1, nrow(comp_gdp), n_energy_above1, nrow(comp_gdp)))

#### Is this an income effect, or a Europe-specific / bloc-specific one? ####
# A hand-coded "Europe vs. non-Europe" split isn't the right test: it lumps
# rich non-European countries (USA, Canada, Australia, Japan, South Korea)
# in with poorer ones (Brazil, China, Indonesia), which is not the
# comparison this needs, and doesn't match how the paper itself already
# frames this argument (Conclusions: "the North's short provisioning day"
# -- standard unequal-exchange vocabulary, not "Europe"). Rather than
# invent a Global-North/South boundary by hand (Russia, China, Turkey,
# and Eastern Europe are all genuinely contested cases in that
# literature), use the World Bank's own income-group classification --
# already pulled in script 9's `wdi_raw` (WDI(..., extra = TRUE) returns
# an `income` column), external and non-arbitrary.

income_lookup <- wdi_raw %>%
  filter(region != "Aggregates", !is.na(income)) %>%
  group_by(iso3c) %>%
  slice_max(year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(country = iso3c, income_group = income)

comp_income <- comp_gdp %>%
  inner_join(income_lookup, by = "country") %>%
  mutate(income_group = factor(income_group,
           levels = c("Low income", "Lower middle income",
                      "Upper middle income", "High income")))

cat("\n---- N, median, mean ratio by World Bank income group ----\n")
print(comp_income %>%
        group_by(income_group) %>%
        summarise(n = n(),
                   median_ratio_time = median(ratio_time, na.rm = TRUE),
                   median_ratio_energy = median(ratio_energy, na.rm = TRUE),
                   mean_ratio_time = mean(ratio_time, na.rm = TRUE),
                   mean_ratio_energy = mean(ratio_energy, na.rm = TRUE),
                   .groups = "drop"))

cat("\n---- Income slope WITHIN High-income and WITHIN non-High-income separately ----\n")
comp_income2 <- comp_income %>%
  mutate(income_bin = ifelse(income_group == "High income", "High income", "Not high income"))
fit_time_hi     <- lm(log_ratio_time ~ log_gdp_pcap, data = filter(comp_income2, income_bin == "High income"))
fit_time_lo     <- lm(log_ratio_time ~ log_gdp_pcap, data = filter(comp_income2, income_bin == "Not high income"))
fit_energy_hi   <- lm(log_ratio_energy ~ log_gdp_pcap, data = filter(comp_income2, income_bin == "High income"))
fit_energy_lo   <- lm(log_ratio_energy ~ log_gdp_pcap, data = filter(comp_income2, income_bin == "Not high income"))
cat(sprintf("time,   High-income only:     slope = %.4f, p = %.4g, n = %d\n",
            coef(fit_time_hi)["log_gdp_pcap"], summary(fit_time_hi)$coefficients["log_gdp_pcap","Pr(>|t|)"], nrow(filter(comp_income2, income_bin=="High income"))))
cat(sprintf("time,   Not-high-income only: slope = %.4f, p = %.4g, n = %d\n",
            coef(fit_time_lo)["log_gdp_pcap"], summary(fit_time_lo)$coefficients["log_gdp_pcap","Pr(>|t|)"], nrow(filter(comp_income2, income_bin=="Not high income"))))
cat(sprintf("energy, High-income only:     slope = %.4f, p = %.4g, n = %d\n",
            coef(fit_energy_hi)["log_gdp_pcap"], summary(fit_energy_hi)$coefficients["log_gdp_pcap","Pr(>|t|)"], nrow(filter(comp_income2, income_bin=="High income"))))
cat(sprintf("energy, Not-high-income only: slope = %.4f, p = %.4g, n = %d\n",
            coef(fit_energy_lo)["log_gdp_pcap"], summary(fit_energy_lo)$coefficients["log_gdp_pcap","Pr(>|t|)"], nrow(filter(comp_income2, income_bin=="Not high income"))))

cat("\n---- Sanity check: where do rich non-European countries land? ----\n")
print(comp_income %>% filter(country %in% c("USA","CAN","AUS","JPN","KOR")) %>%
        select(country, gdp_pcap_ppp, income_group, ratio_time, ratio_energy))
