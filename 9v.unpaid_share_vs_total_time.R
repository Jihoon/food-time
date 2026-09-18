#### Does the unpaid share of provisioning time shrink where total time is ####
#### short, as the Nested Time Substitution hypothesis (Conclusions,       ####
#### Future work) claims -- "where total provisioning time is short, the   ####
#### unpaid share is smaller and the embodied energy higher"?             ####
#
# That specific share claim has never been tested anywhere in this paper's
# scripts -- it is carried-over hypothesis text. It may also be in tension
# with something the paper DOES establish: unpaid time is repeatedly
# described elsewhere as "sticky," barely moving with income, while total
# CF_time (all-work) falls substantially with income (elasticity ~ -0.46)
# and CF_paid's own elasticity is similar in magnitude to the all-work
# total's (~ -0.44 to -0.52), not obviously steeper. If unpaid stays flat
# while the total shrinks, unpaid's SHARE of a shrinking total should rise,
# not fall -- the opposite of the claim. This script checks directly
# rather than reasons about it secondhand.
#
# Construction: CF_paid and the primary (all-work) CF share the exact same
# protein denominator (g_protein_per_cap_day) per country, so the ratio
# CF_paid / CF_time_allwork equals paid_hours / allwork_hours directly --
# the g_protein_per_cap_day term cancels. No new hours data needs to be
# pulled; this is just a ratio of two already-built CF variants.
#
# Run AFTER 2.analyze_result.R, 9.capability_set_income_control.R (df,
# wdi_latest), and 9t.energy_paid_labor_dual_mediator.R (cf_paid_df) in the
# same session.

library(tidyverse)

share_df <- df %>%
  select(country, hr_per_50g_protein_allwork = hr_per_50g_protein,
         g_protein_per_cap_day, gdp_pcap_ppp, gdp_per_worker) %>%
  inner_join(cf_paid_df %>% select(country, hr_per_50g_protein_paid = hr_per_50g_protein),
             by = "country") %>%
  mutate(paid_share   = hr_per_50g_protein_paid / hr_per_50g_protein_allwork,
         unpaid_share = 1 - paid_share,
         log_cf_allwork = log(hr_per_50g_protein_allwork),
         log_gdp_pcap   = log(gdp_pcap_ppp))

cat(sprintf("Countries with both all-work and paid-only CF (n = %d)\n", nrow(share_df)))
cat(sprintf("Unpaid share: mean = %.3f, median = %.3f, range [%.3f, %.3f]\n",
            mean(share_df$unpaid_share), median(share_df$unpaid_share),
            min(share_df$unpaid_share), max(share_df$unpaid_share)))

cat("\n---- Sorted by total CF_time (ascending: shortest provisioning time first) ----\n")
print(share_df %>% arrange(hr_per_50g_protein_allwork) %>%
        select(country, hr_per_50g_protein_allwork, unpaid_share, gdp_pcap_ppp) %>%
        mutate(across(c(hr_per_50g_protein_allwork, unpaid_share), ~round(.x, 3)),
               gdp_pcap_ppp = round(gdp_pcap_ppp, 0)),
      n = Inf)

#### The actual test: does unpaid_share correlate with total CF_time, or ####
#### with income? The claim needs a NEGATIVE correlation with CF_time     ####
#### (short time -> smaller unpaid share) or a NEGATIVE correlation with  ####
#### log_gdp_pcap in the other direction (higher income -> smaller share, ####
#### since higher income means shorter CF_time) -- same prediction, two   ####
#### ways of stating it.                                                  ####

cor_share_cf <- cor.test(share_df$unpaid_share, share_df$log_cf_allwork)
cor_share_gdp <- cor.test(share_df$unpaid_share, share_df$log_gdp_pcap)

cat("\n---- Does unpaid share track total time or income? ----\n")
cat(sprintf("unpaid_share vs. log(CF_time, all-work): r = %.3f (p = %.3g, n = %d)\n",
            cor_share_cf$estimate, cor_share_cf$p.value, nrow(share_df)))
cat("  Claim predicts a POSITIVE r here (shorter total time = smaller unpaid share\n")
cat("  = these two move together in the same direction on a positive scale).\n")
cat(sprintf("unpaid_share vs. log(GDP per capita): r = %.3f (p = %.3g, n = %d)\n",
            cor_share_gdp$estimate, cor_share_gdp$p.value, nrow(share_df)))
cat("  Claim predicts a NEGATIVE r here (higher income = smaller unpaid share).\n")

cat("\n#### What to look at ####\n")
cat("If unpaid_share vs. CF_time comes back negative (not positive) and/or\n")
cat("unpaid_share vs. GDP comes back positive (not negative), the claim in\n")
cat("Conclusions/Future work is backwards on this dataset, consistent with\n")
cat("the 'unpaid time is sticky' finding stated elsewhere -- fix the text,\n")
cat("don't just soften it.\n")
