#### Domestic-origin vs import-origin effort per gram: the "diagonal" ratio ####
#
# 9q compares protein STREAMS: all effort (any origin) behind domestically
# produced protein, per gram, vs all effort behind imported protein, per
# gram. The wage-differential argument in "International Embodied Labor and
# Trade" is about effort ORIGIN: hours worked abroad vs hours worked at home.
# Effort origin has no protein denominator of its own (footnote 2), so the
# closest defensible test uses only the diagonal cells of the mosaic:
#
#   cf_diag_domestic = domestic-origin effort in the domestic stream / domestic grams
#   cf_diag_import   = import-origin effort in the import stream   / imported grams
#
# dropping the off-diagonal cells (imported inputs or machinery in domestic
# food; domestic effort in food exported for processing and re-imported).
# If ratio_diag tracks 9q's ratio_time closely, the stream result stands in
# for the origin argument; if it diverges, the off-diagonal cells matter and
# the prose should say which comparison it is making.
#
# Run AFTER 9q.mosaic_domestic_import_cf_ratio.R in the same session
# (reuses df_mosaic_9q, df_mosaic_energy_9q, comparison, meaningful, comp_gdp).

diag_cf <- function(df, effort_col, unit_scale = 1) {
  df %>%
    mutate(stream_origin = ifelse(footprint_type == "domestic_per_capita", "domestic", "import")) %>%
    filter(effort_origin == stream_origin) %>%
    group_by(country, footprint_type) %>%
    summarise(effort = sum(.data[[effort_col]], na.rm = TRUE) * unit_scale,
              g = first(g_protein_per_cap_day), .groups = "drop") %>%
    filter(g > 0) %>%
    mutate(cf = effort / g * 50) %>%
    select(country, footprint_type, cf) %>%
    pivot_wider(names_from = footprint_type, values_from = cf, names_prefix = "cf_")
}

diag_time <- diag_cf(df_mosaic_9q, "hr_per_cap_day") %>%
  rename(cf_diag_time_dom = cf_domestic_per_capita, cf_diag_time_imp = cf_import_per_capita) %>%
  mutate(ratio_diag_time = cf_diag_time_imp / cf_diag_time_dom)

diag_energy <- diag_cf(df_mosaic_energy_9q, "gj_per_cap_day", unit_scale = 1000) %>%
  rename(cf_diag_en_dom = cf_domestic_per_capita, cf_diag_en_imp = cf_import_per_capita) %>%
  mutate(ratio_diag_energy = cf_diag_en_imp / cf_diag_en_dom)

# share of each stream's effort that sits on the diagonal (how much 9q and this agree by construction)
diag_share <- df_mosaic_9q %>%
  mutate(stream_origin = ifelse(footprint_type == "domestic_per_capita", "domestic", "import"),
         on_diag = effort_origin == stream_origin) %>%
  group_by(country, footprint_type) %>%
  summarise(diag_share = sum(hr_per_cap_day[on_diag], na.rm = TRUE) / sum(hr_per_cap_day, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = footprint_type, values_from = diag_share, names_prefix = "diag_share_")

side_by_side <- comparison %>%
  select(country, ratio_time, ratio_energy) %>%
  inner_join(diag_time %>% select(country, ratio_diag_time), by = "country") %>%
  inner_join(diag_energy %>% select(country, ratio_diag_energy), by = "country") %>%
  left_join(diag_share, by = "country") %>%
  filter(country %in% meaningful)

cat("\n---- Stream ratio (9q) vs diagonal ratio, per country (imports > 2%) ----\n")
print(side_by_side %>% mutate(across(where(is.numeric), ~round(., 2))), n = Inf)

cat(sprintf("\ncor(ratio_time, ratio_diag_time) = %.3f | cor(ratio_energy, ratio_diag_energy) = %.3f\n",
            cor(side_by_side$ratio_time, side_by_side$ratio_diag_time, use = "complete.obs"),
            cor(side_by_side$ratio_energy, side_by_side$ratio_diag_energy, use = "complete.obs")))
cat(sprintf("median diagonal share of stream effort: domestic stream %.2f, import stream %.2f\n",
            median(side_by_side$diag_share_domestic_per_capita, na.rm = TRUE),
            median(side_by_side$diag_share_import_per_capita, na.rm = TRUE)))
cat(sprintf("%d/%d countries have ratio_diag_time > 1 (9q: ratio_time > 1 in %d)\n",
            sum(side_by_side$ratio_diag_time > 1, na.rm = TRUE), nrow(side_by_side),
            sum(side_by_side$ratio_time > 1, na.rm = TRUE)))

diag_gdp <- side_by_side %>%
  inner_join(comp_gdp %>% select(country, log_gdp_pcap), by = "country") %>%
  mutate(log_ratio_diag_time = log(ratio_diag_time),
         log_ratio_diag_energy = log(ratio_diag_energy)) %>%
  filter(is.finite(log_ratio_diag_time), is.finite(log_ratio_diag_energy))

cat("\n---- Does log(diagonal import/domestic CF ratio) rise with log GDP per capita? ----\n")
for (v in c("log_ratio_diag_time", "log_ratio_diag_energy")) {
  f <- lm(as.formula(paste(v, "~ log_gdp_pcap")), data = diag_gdp)
  s <- summary(f)$coefficients["log_gdp_pcap", ]
  cat(sprintf("%-22s slope = %.4f, SE = %.4f, p = %.3g, n = %d\n", v, s[1], s[2], s[4], nrow(diag_gdp)))
}
cat("(9q stream version: time slope = 0.7745, p = 0.009; energy slope = -0.0972, p = 0.63)\n")

#### Column cut: foreign share of the paid hours behind ALL consumed protein ####
# Same denominator on both sides (total protein), so this is a share of
# labor, not a per-gram efficiency. It quantifies "richer countries draw on
# labor embodied in other countries' exports"; the row/diagonal ratios above
# quantify whether those hours are less efficient.

foreign_share <- df_mosaic_9q %>%
  group_by(country, effort_origin) %>%
  summarise(hr = sum(hr_per_cap_day, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = effort_origin, values_from = hr, names_prefix = "hr_") %>%
  mutate(foreign_share_paid = hr_import / (hr_domestic + hr_import)) %>%
  inner_join(comp_gdp %>% select(country, log_gdp_pcap), by = "country")

cat("\n---- Foreign share of paid hours behind all consumed protein ----\n")
print(foreign_share %>% arrange(log_gdp_pcap) %>%
        mutate(gdp_pcap = round(exp(log_gdp_pcap)), foreign_share_paid = round(foreign_share_paid, 3)) %>%
        select(country, gdp_pcap, foreign_share_paid), n = Inf)
cat(sprintf("median foreign share = %.3f | min = %.3f (%s) | max = %.3f (%s)\n",
            median(foreign_share$foreign_share_paid),
            min(foreign_share$foreign_share_paid), foreign_share$country[which.min(foreign_share$foreign_share_paid)],
            max(foreign_share$foreign_share_paid), foreign_share$country[which.max(foreign_share$foreign_share_paid)]))
f_share <- lm(qlogis(foreign_share_paid) ~ log_gdp_pcap, data = foreign_share)
s <- summary(f_share)$coefficients["log_gdp_pcap", ]
cat(sprintf("logit(foreign share) on log GDP per capita: slope = %.3f, SE = %.3f, p = %.3g, n = %d\n",
            s[1], s[2], s[4], nrow(foreign_share)))
r_share <- cor.test(foreign_share$foreign_share_paid, foreign_share$log_gdp_pcap)
cat(sprintf("Pearson r (share, log GDP) = %.3f, p = %.3g\n", r_share$estimate, r_share$p.value))
