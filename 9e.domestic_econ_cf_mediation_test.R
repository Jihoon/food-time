#### Alternative CF definition: economic (paid) labor only, domestic-only ####
#
# Motivating question: the well-documented link in the literature is
# agricultural productivity -> nutrition, using paid/economic labor. Whether
# TIME-EFFICIENCY specifically -- including household unpaid time, and net of
# trade -- carries independent explanatory weight is what this paper actually
# tests, and it hasn't been asked before. This script isolates the more
# conventional, narrower measure for comparison: CF built from domestic PAID
# labor only (food + non-food sectors, both genders, no household time) over
# domestically-produced-and-consumed protein (excludes imports on both the
# labor and the protein side) -- as opposed to the all-work, all-origin CF
# used in scripts 9/9b, which is this paper's actual object of interest.
#
# Y (protein supply) and the GDP controls are kept IDENTICAL to every earlier
# mediation test (df's g_protein_per_cap_day, gdp_pcap_ppp, gdp_per_worker) --
# only the mediator (M = CF) changes, so ADE/ACME magnitudes are directly
# comparable to the all-work results already reported.
#
# Restricted to EXIO-individually-modeled countries (excludes RoW
# aggregates), matching cf_df's n=33 sample in
# 9.capability_set_income_control.R, for the same reason given there (a
# RoW-collapsed version could be added the same way 9b/9d did, if wanted).
#
# Run AFTER 9.capability_set_income_control.R in the same session (reuses its
# `df`, `wdi_latest`, `partial_cor()`, `mediation_decomp()`; needs
# `effort_consumption_df`, `summary_pro_df_long`, `region_to_iso_1to1` from
# 2.analyze_result.R).

domestic_protein_df <- summary_pro_df_long %>%
  filter(footprint_type == "domestic_per_capita") %>%
  transmute(country = as.character(country), g_protein_per_cap_day_domestic = per_capita_value)

domestic_econ_cf_df <- effort_consumption_df %>%
  filter(type %in% c("hr_m", "hr_f"), protein_source == "domestic", effort_origin == "domestic",
         !grepl("RoW", exio_region)) %>%
  mutate(country = region_to_iso_1to1[as.character(exio_region)]) %>%
  group_by(country) %>%
  summarise(hr_per_cap_day_domestic_econ = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  inner_join(domestic_protein_df, by = "country") %>%
  filter(hr_per_cap_day_domestic_econ > 0, g_protein_per_cap_day_domestic > 0) %>%
  mutate(hr_per_50g_protein_domestic_econ = hr_per_cap_day_domestic_econ / g_protein_per_cap_day_domestic * 50)

df_domestic_econ <- domestic_econ_cf_df %>%
  select(country, hr_per_50g_protein_domestic_econ) %>%
  inner_join(df %>% select(country, g_protein_per_cap_day, gdp_pcap_ppp, gdp_per_worker), by = "country") %>%
  mutate(log_cf         = log(hr_per_50g_protein_domestic_econ),
         log_protein    = log(g_protein_per_cap_day),
         log_gdp_pcap   = log(gdp_pcap_ppp),
         log_gdp_worker = log(gdp_per_worker))

cat(sprintf("\nDomestic-economic-only CF: n = %d (vs. n = %d in the all-work version)\n",
            nrow(df_domestic_econ), nrow(df)))

zero_order_de    <- cor.test(df_domestic_econ$log_cf, df_domestic_econ$log_protein)
pc_gdp_pcap_de   <- partial_cor(df_domestic_econ$log_cf, df_domestic_econ$log_protein, df_domestic_econ$log_gdp_pcap)
pc_gdp_worker_de <- partial_cor(df_domestic_econ$log_cf, df_domestic_econ$log_protein, df_domestic_econ$log_gdp_worker)

cat("\n---- Domestic-economic-only CF vs. protein supply, log-log ----\n")
cat(sprintf("Zero-order r                             = %.3f (p = %.3g, n = %d)\n",
            zero_order_de$estimate, zero_order_de$p.value, nrow(df_domestic_econ)))
cat(sprintf("Partial r | GDP per capita (PPP)         = %.3f (p = %.3g, n = %d)\n",
            pc_gdp_pcap_de$r, pc_gdp_pcap_de$p, pc_gdp_pcap_de$n))
cat(sprintf("Partial r | GDP per worker (labor prod.) = %.3f (p = %.3g, n = %d)\n",
            pc_gdp_worker_de$r, pc_gdp_worker_de$p, pc_gdp_worker_de$n))

med_pcap_de   <- mediation_decomp(df_domestic_econ, "log_gdp_pcap",   "log_cf", "log_protein")
med_worker_de <- mediation_decomp(df_domestic_econ, "log_gdp_worker", "log_cf", "log_protein")

cat("\n---- Mediation decomposition (domestic-economic-only CF): GDP -> CF -> protein supply ----\n")
cat(sprintf("Via GDP per capita (PPP):         total = %.3f | direct (ADE) = %.3f | indirect via CF (ACME) = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g), n = %d\n",
            med_pcap_de$total, med_pcap_de$direct, med_pcap_de$indirect,
            med_pcap_de$prop_mediated * 100, med_pcap_de$sobel_z, med_pcap_de$sobel_p, med_pcap_de$n))
cat(sprintf("Via GDP per worker (labor prod.): total = %.3f | direct (ADE) = %.3f | indirect via CF (ACME) = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g), n = %d\n",
            med_worker_de$total, med_worker_de$direct, med_worker_de$indirect,
            med_worker_de$prop_mediated * 100, med_worker_de$sobel_z, med_worker_de$sobel_p, med_worker_de$n))

cat("\n---- Path a for this CF definition: does GDP predict domestic-economic-only CF? ----\n")
print(summary(lm(log_cf ~ log_gdp_pcap, data = df_domestic_econ)))
