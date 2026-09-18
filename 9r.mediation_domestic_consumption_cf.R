#### Primary mediation test, rebuilt on domestically-consumed CF ####
#### (tradeoff_protein_allwork, not tradeoff_protein_allwork_consump) ####
#
# Switches the mediator from "labor embodied in total consumption
# (domestic + imported food)" to "labor embodied in domestically-consumed
# food only" -- matches what the mediation test is conceptually about
# (Methods: Mediation Model Specification note added alongside this), and
# matches that unpaid labor, the larger and stickier component of
# all-work CF, cannot be imported in the first place.
#
# tradeoff_protein_allwork has one row per country x gender (type =
# hr_m/hr_f), already restricted to domestically-consumed protein
# (cat == "domestic", i.e. footprint_type == "domestic_per_capita") and
# GHD-covered countries. Sum across gender to get the all-work total per
# country, mirroring how the original cf_df summed across sector x origin
# buckets in the _consump version. Note this does NOT strip
# foreign-sourced labor embedded in domestic production (imported inputs/
# machinery) -- it only drops the imported-*food* consumption stream.
# That is the intended scope (see Methods note), not an oversight.
#
# Y (protein supply / nutrition outcome) stays TOTAL consumption
# (domestic + imported), pulled from tradeoff_protein_allwork_consump --
# NOT domestic-only. Nutrition doesn't care where the food came from;
# only M (the mediator) has a principled reason to be domestic-only. An
# earlier version of this script pulled Y from tradeoff_protein_allwork
# too, which silently made Y domestic-only protein supply and produced a
# spurious-looking reversal (total effect went negative, partial r
# flipped to a strongly significant -0.65) -- both artifacts of GDP
# correlating with import-reliance, not a real finding. Fixed here.
#
# Run AFTER 2.analyze_result.R (builds tradeoff_protein_allwork and
# tradeoff_protein_allwork_consump) and 9.capability_set_income_control.R
# (reuses wdi_latest, mediation_decomp(), partial_cor()) in the same
# session.

cf_df_domestic <- tradeoff_protein_allwork %>%
  filter(!is_row) %>%
  group_by(country) %>%
  summarise(hr_per_50g_protein = sum(hr_per_50g_protein, na.rm = TRUE),
            .groups = "drop")

# Y (nutrition outcome) stays TOTAL protein supply -- domestic + imported
# -- the same variable used everywhere else in this paper. Nutrition
# doesn't care where the food came from; only M (the mediator) has a
# principled reason to be domestic-only (Methods note: unpaid labor can't
# be imported). Pull total protein supply from the _consump object, NOT
# from tradeoff_protein_allwork's own g_protein_per_cap_day, which is
# domestic-only and would silently change what Y means if used here.
protein_total <- tradeoff_protein_allwork_consump %>%
  filter(!is_row) %>%
  group_by(country) %>%
  summarise(g_protein_per_cap_day = first(g_protein_per_cap_day), .groups = "drop")

df_domestic <- cf_df_domestic %>%
  inner_join(protein_total, by = "country") %>%
  inner_join(wdi_latest, by = c("country" = "iso3c")) %>%
  filter(hr_per_50g_protein > 0, g_protein_per_cap_day > 0) %>%
  mutate(log_cf       = log(hr_per_50g_protein),
         log_protein  = log(g_protein_per_cap_day),
         log_gdp_pcap = log(gdp_pcap_ppp))

cat(sprintf("Countries with domestic-consumption CF, protein supply, GDP: %d\n", nrow(df_domestic)))

zero_order_d  <- cor.test(df_domestic$log_cf, df_domestic$log_protein)
pc_gdp_pcap_d <- partial_cor(df_domestic$log_cf, df_domestic$log_protein, df_domestic$log_gdp_pcap)

cat("\n---- Domestic-consumption CF vs. protein supply, log-log ----\n")
cat(sprintf("Zero-order r                     = %.3f (p = %.3g, n = %d)\n",
            zero_order_d$estimate, zero_order_d$p.value, nrow(df_domestic)))
cat(sprintf("Partial r | GDP per capita (PPP) = %.3f (p = %.3g, n = %d)\n",
            pc_gdp_pcap_d$r, pc_gdp_pcap_d$p, pc_gdp_pcap_d$n))

med_domestic <- mediation_decomp(df_domestic, "log_gdp_pcap", "log_cf", "log_protein")

cat("\n---- Mediation decomposition (domestic-consumption CF): GDP -> CF -> protein supply ----\n")
cat(sprintf("total = %.3f | direct (ADE) = %.3f | indirect via CF (ACME) = %.3f (%.1f%% of total, Sobel z = %.2f, p = %.3g), n = %d\n",
            med_domestic$total, med_domestic$direct, med_domestic$indirect,
            med_domestic$prop_mediated * 100, med_domestic$sobel_z, med_domestic$sobel_p, med_domestic$n))

cat("\n---- Path a: does GDP predict domestic-consumption CF? ----\n")
print(summary(lm(log_cf ~ log_gdp_pcap, data = df_domestic)))

cat("\n#### For comparison, the _consump (total-consumption) version already in the draft ####\n")
cat("path a:   a1 = -0.52 (n=33) / -0.46 (n=38, RoW-collapsed)\n")
cat("path b:   partial r = -0.13 (p=0.45, n=38); mediated share 16-24%, Sobel p>0.3 throughout\n")
cat("\nIf df_domestic's n differs a lot from 33-38, or the path a/b numbers move\n")
cat("substantially, that's the signal this needs to become the new primary result\n")
cat("(with everything downstream -- robustness checks, energy dual-mediator,\n")
cat("interaction tests, gender split -- re-run to match), not just a footnote.\n")
