#### Gender-specific version of script 9: does CF's protein-supply link (net of GDP) differ by gender? ####
#
# Companion to 9.capability_set_income_control.R -- same logic, but instead of
# summing hr_m + hr_f into one combined CF per country (as cf_df does there),
# keeps each gender's hours separate. tradeoff_protein_allwork_consump already
# carries one row per (country, type) with type in {hr_m, hr_f}, so this needs
# no new upstream computation -- just a different filter.
#
# Run AFTER 9.capability_set_income_control.R in the same session (reuses its
# `wdi_latest`, `partial_cor()`, and `mediation_decomp()`).

gender_cf_df <- function(type_filter) {
  tradeoff_protein_allwork_consump %>%
    filter(!is_row, type == type_filter) %>%
    distinct(country, .keep_all = TRUE) %>%
    select(country, hr_per_50g_protein, g_protein_per_cap_day) %>%
    inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
    filter(hr_per_50g_protein > 0, g_protein_per_cap_day > 0) %>%
    mutate(log_cf         = log(hr_per_50g_protein),
           log_protein    = log(g_protein_per_cap_day),
           log_gdp_pcap   = log(gdp_pcap_ppp),
           log_gdp_worker = log(gdp_per_worker))
}

df_f <- gender_cf_df("hr_f")
df_m <- gender_cf_df("hr_m")

run_gender_block <- function(df, label) {
  cat(sprintf("\n==== %s (n = %d) ====\n", label, nrow(df)))

  zero_order    <- cor.test(df$log_cf, df$log_protein)
  pc_gdp_pcap   <- partial_cor(df$log_cf, df$log_protein, df$log_gdp_pcap)
  pc_gdp_worker <- partial_cor(df$log_cf, df$log_protein, df$log_gdp_worker)

  cat(sprintf("Zero-order r                             = %.3f (p = %.3g)\n",
              zero_order$estimate, zero_order$p.value))
  cat(sprintf("Partial r | GDP per capita (PPP)         = %.3f (p = %.3g, n = %d)\n",
              pc_gdp_pcap$r, pc_gdp_pcap$p, pc_gdp_pcap$n))
  cat(sprintf("Partial r | GDP per worker (labor prod.) = %.3f (p = %.3g, n = %d)\n",
              pc_gdp_worker$r, pc_gdp_worker$p, pc_gdp_worker$n))

  med_pcap   <- mediation_decomp(df, "log_gdp_pcap",   "log_cf", "log_protein")
  med_worker <- mediation_decomp(df, "log_gdp_worker", "log_cf", "log_protein")

  cat(sprintf("Mediation via GDP/cap:    total = %.3f | ADE = %.3f | ACME = %.3f (%.1f%%, Sobel z = %.2f, p = %.3g)\n",
              med_pcap$total, med_pcap$direct, med_pcap$indirect,
              med_pcap$prop_mediated * 100, med_pcap$sobel_z, med_pcap$sobel_p))
  cat(sprintf("Mediation via GDP/worker: total = %.3f | ADE = %.3f | ACME = %.3f (%.1f%%, Sobel z = %.2f, p = %.3g)\n",
              med_worker$total, med_worker$direct, med_worker$indirect,
              med_worker$prop_mediated * 100, med_worker$sobel_z, med_worker$sobel_p))

  invisible(list(df = df, zero_order = zero_order, pc_gdp_pcap = pc_gdp_pcap,
                 pc_gdp_worker = pc_gdp_worker, med_pcap = med_pcap, med_worker = med_worker))
}

res_f <- run_gender_block(df_f, "FEMALE hours only")
res_m <- run_gender_block(df_m, "MALE hours only")

#### Visualize both genders' partial relationships on one plot ####

plot_df <- bind_rows(
  df_f %>% mutate(resid_cf      = resid(lm(log_cf ~ log_gdp_pcap)),
                   resid_protein = resid(lm(log_protein ~ log_gdp_pcap)),
                   gender = "Female"),
  df_m %>% mutate(resid_cf      = resid(lm(log_cf ~ log_gdp_pcap)),
                   resid_protein = resid(lm(log_protein ~ log_gdp_pcap)),
                   gender = "Male")
)

p_gender_partial <- ggplot(plot_df, aes(resid_cf, resid_protein, color = gender)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  labs(x = "CF residual (log hr/50g protein, GDP/cap partialled out)",
       y = "Protein supply residual (log g/cap/day, GDP/cap partialled out)",
       title = "Partial relationship: CF vs. protein supply by gender, controlling for GDP per capita (PPP)",
       color = NULL) +
  theme_minimal()
print(p_gender_partial)
ggsave("results/partial_cf_vs_protein_control_gdp_by_gender.pdf", p_gender_partial, width = 9, height = 6)

#### RoW-inclusive version (n=38): same gender split, collapsing each RoW ####
#### aggregate to one population-weighted point, as script 9 section 7 does ####
#### for the combined-gender case. Reuses `row_lookup`, `row_full_pop`, and  ####
#### `wdi_latest` from 9.capability_set_income_control.R (already in the     ####
#### session after running its section 7).                                  ####

gender_row_collapsed <- function(type_filter) {
  tradeoff_protein_allwork_consump %>%
    filter(is_row, type == type_filter) %>%
    distinct(country, .keep_all = TRUE) %>%
    select(country, hr_per_cap_day, g_protein_per_cap_day) %>%
    inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
    filter(hr_per_cap_day > 0, g_protein_per_cap_day > 0) %>%
    left_join(row_lookup, by = "country") %>%
    left_join(pop_data_yr, by = c("country" = "iso3c")) %>%
    filter(!is.na(population), population > 0) %>%
    group_by(exio_region) %>%
    summarise(n_members_used = n(),
              pop_used = sum(population),
              hr_per_cap_day_region = weighted.mean(hr_per_cap_day, population),
              g_protein_per_cap_day = weighted.mean(g_protein_per_cap_day, population),
              gdp_pcap_ppp          = weighted.mean(gdp_pcap_ppp, population, na.rm = TRUE),
              gdp_per_worker        = weighted.mean(gdp_per_worker, population, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(hr_per_50g_protein = hr_per_cap_day_region / g_protein_per_cap_day * 50) %>%
    left_join(row_full_pop, by = "exio_region") %>%
    mutate(pop_coverage_pct = pop_used / pop_total_region * 100)
}

row_collapsed_f <- gender_row_collapsed("hr_f")
row_collapsed_m <- gender_row_collapsed("hr_m")

df_collapsed_f <- bind_rows(
  df_f %>% select(country, hr_per_50g_protein, g_protein_per_cap_day, gdp_pcap_ppp, gdp_per_worker),
  row_collapsed_f %>% select(country = exio_region, hr_per_50g_protein, g_protein_per_cap_day, gdp_pcap_ppp, gdp_per_worker)
) %>%
  mutate(log_cf = log(hr_per_50g_protein), log_protein = log(g_protein_per_cap_day),
         log_gdp_pcap = log(gdp_pcap_ppp), log_gdp_worker = log(gdp_per_worker))

df_collapsed_m <- bind_rows(
  df_m %>% select(country, hr_per_50g_protein, g_protein_per_cap_day, gdp_pcap_ppp, gdp_per_worker),
  row_collapsed_m %>% select(country = exio_region, hr_per_50g_protein, g_protein_per_cap_day, gdp_pcap_ppp, gdp_per_worker)
) %>%
  mutate(log_cf = log(hr_per_50g_protein), log_protein = log(g_protein_per_cap_day),
         log_gdp_pcap = log(gdp_pcap_ppp), log_gdp_worker = log(gdp_per_worker))

cat(sprintf("\n[RoW-collapsed] Female sample: %d EXIO-modeled + %d RoW-aggregate points = %d rows\n",
            nrow(df_f), nrow(row_collapsed_f), nrow(df_collapsed_f)))
cat(sprintf("[RoW-collapsed] Male sample:   %d EXIO-modeled + %d RoW-aggregate points = %d rows\n",
            nrow(df_m), nrow(row_collapsed_m), nrow(df_collapsed_m)))

res_f_collapsed <- run_gender_block(df_collapsed_f, "FEMALE hours only [RoW-collapsed, n=38]")
res_m_collapsed <- run_gender_block(df_collapsed_m, "MALE hours only [RoW-collapsed, n=38]")

plot_df_collapsed <- bind_rows(
  df_collapsed_f %>% mutate(resid_cf      = resid(lm(log_cf ~ log_gdp_pcap)),
                             resid_protein = resid(lm(log_protein ~ log_gdp_pcap)),
                             gender = "Female"),
  df_collapsed_m %>% mutate(resid_cf      = resid(lm(log_cf ~ log_gdp_pcap)),
                             resid_protein = resid(lm(log_protein ~ log_gdp_pcap)),
                             gender = "Male")
)

p_gender_partial_collapsed <- ggplot(plot_df_collapsed, aes(resid_cf, resid_protein, color = gender)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  labs(x = "CF residual (log hr/50g protein, GDP/cap partialled out)",
       y = "Protein supply residual (log g/cap/day, GDP/cap partialled out)",
       title = "Partial relationship: CF vs. protein supply by gender, controlling for GDP per capita (PPP) [RoW-collapsed, n=38]",
       color = NULL) +
  theme_minimal()
print(p_gender_partial_collapsed)
ggsave("results/partial_cf_vs_protein_control_gdp_by_gender_row_collapsed.pdf", p_gender_partial_collapsed, width = 9, height = 6)
