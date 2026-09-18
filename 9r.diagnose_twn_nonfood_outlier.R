#### Diagnose: TWN shows the world's biggest domestic non-food ENERGY   ####
#### footprint and biggest imported non-food LABOR footprint in the    ####
#### mosaic. Is this a genuine economic pattern (Taiwan's energy-      ####
#### intensive, trade-heavy non-food economy) or a data/index bug?     ####
#
# "Domestic effort" / "import effort" here follow the effort_origin
# convention used throughout 2.analyze_result.R and 9q: domestic = the
# non-food EXIO region doing the work equals the consuming FABIO country's
# own EXIO region (diagonal block); import = any other region (off-diagonal).
# This is crossed with the *separate* domestic/import-protein bucket
# (whether the food product itself was domestically-traded or imported,
# per FABIO_y_hh) via fp_domcons / fp_impcons.
#
# Reuses fp_domcons/fp_impcons (data/fp_domcons_2020.rds, data/fp_impcons_2020.rds
# -- already computed and saved by 2.analyze_result.R section 1.3) and objects
# left in the R environment by a full run of 0.mrio_prep.R through
# 2.analyze_result.R (regions, FABIO_reg, n_reg_EXIO, n_nf, exio_nonfood_sectors,
# region_name_of_index, region_population). Run AFTER 2.analyze_result.R.

library(tidyverse)
library(Matrix)

fp_domcons = readRDS("data/fp_domcons_2020.rds")
fp_impcons = readRDS("data/fp_impcons_2020.rds")

own_rows_fun    = function(reg) ((reg - 1) * n_nf + 1):(reg * n_nf)
import_rows_fun = function(reg) setdiff(seq_len(n_reg_EXIO * n_nf), own_rows_fun(reg))

# Per-country, per-capita-per-day value for one non-food metric matrix (8575 x 187),
# restricted to either the country's own-region rows (domestic effort) or every
# other region's rows (import effort), for every FABIO country -- to confirm/refute
# "biggest in the world" and show the runner-up.
compare_across_countries <- function(mat, rows_fun, metric_label) {
  vals = sapply(seq_len(nrreg), function(jj) {
    iso = regions$iso3c[jj]
    reg = FABIO_reg$EXIOBASE_code[FABIO_reg$ISO == iso]
    pop = region_population[reg]
    if (is.na(pop) || pop == 0) return(NA_real_)
    sum(mat[rows_fun(reg), jj]) / pop * 1e6 / 365
  })
  df = data.frame(country = regions$iso3c, value = vals) %>% arrange(desc(value))
  cat("\nTop 10 countries --", metric_label, "--\n")
  print(head(df, 10))
  df
}

en_domestic_effort_both_buckets = fp_domcons$nonfood$en + fp_impcons$nonfood$en
hr_import_effort_both_buckets   = fp_domcons$nonfood$hr_m + fp_domcons$nonfood$hr_f +
                                   fp_impcons$nonfood$hr_m + fp_impcons$nonfood$hr_f

rank_energy = compare_across_countries(en_domestic_effort_both_buckets, own_rows_fun,
                                        "domestic non-food ENERGY per cap/day (both protein buckets)")
rank_labor  = compare_across_countries(hr_import_effort_both_buckets, import_rows_fun,
                                        "imported non-food LABOR per cap/day (both protein buckets)")

# ---- Drill into TWN specifically: which origin EXIO region / which EXIO ----
# ---- non-food sector is driving each number.                            ----
diagnose_country <- function(country_iso) {
  j = which(regions$iso3c == country_iso)
  own_region = FABIO_reg$EXIOBASE_code[FABIO_reg$ISO == country_iso]
  own_rows = own_rows_fun(own_region)
  import_rows = import_rows_fun(own_region)
  pop = region_population[own_region]

  summarize_block <- function(mat, rows, label) {
    total = sum(mat[rows, j])
    cat(sprintf("\n%s: %.4g (per-cap/day: %.4g)\n", label, total, total / pop * 1e6 / 365))

    by_region = sapply(seq_len(n_reg_EXIO), function(r) {
      rr = intersect(own_rows_fun(r), rows)
      if (length(rr) == 0) return(0)
      sum(mat[rr, j])
    })
    df_r = data.frame(region = region_name_of_index, value = by_region) %>%
      filter(value != 0) %>% arrange(desc(value)) %>% mutate(share_pct = value / sum(value) * 100)
    cat("By origin EXIO region (top 8):\n"); print(head(df_r, 8))

    by_sector = sapply(seq_len(n_nf), function(s) {
      rr = intersect(seq(s, n_reg_EXIO * n_nf, by = n_nf), rows)
      if (length(rr) == 0) return(0)
      sum(mat[rr, j])
    })
    df_s = data.frame(sector = exio_nonfood_sectors, value = by_sector) %>%
      filter(value != 0) %>% arrange(desc(value)) %>% mutate(share_pct = value / sum(value) * 100)
    cat("By EXIO non-food sector (top 8):\n"); print(head(df_s, 8))

    invisible(list(total = total, by_region = df_r, by_sector = df_s))
  }

  cat("\n========", country_iso, ": domestic non-food ENERGY (own-region effort) ========\n")
  r1 = summarize_block(fp_domcons$nonfood$en, own_rows, "domestic-protein bucket x domestic-effort (energy)")
  r2 = summarize_block(fp_impcons$nonfood$en, own_rows, "import-protein bucket x domestic-effort (energy)")

  cat("\n========", country_iso, ": imported non-food LABOR, hr_m+hr_f (off-region effort) ========\n")
  hr_domcons = fp_domcons$nonfood$hr_m + fp_domcons$nonfood$hr_f
  hr_impcons = fp_impcons$nonfood$hr_m + fp_impcons$nonfood$hr_f
  r3 = summarize_block(hr_domcons, import_rows, "domestic-protein bucket x import-effort (labor)")
  r4 = summarize_block(hr_impcons, import_rows, "import-protein bucket x import-effort (labor)")

  list(dom_energy_domprotein = r1, dom_energy_impprotein = r2,
       imp_labor_domprotein = r3, imp_labor_impprotein = r4)
}

res_twn = diagnose_country("TWN")

# Sanity comparator: a large, food-import-heavy but non-outlier economy, to see
# whether TWN's breakdown looks qualitatively different (concentrated in one
# sector/region -- possible artifact) or just scaled up (genuine pattern).
res_kor = diagnose_country("KOR")


#### Rate vs. volume: is TWN's domestic-energy number a genuine rate outlier, ####
#### or just a lot of embodied food mass moving through an ordinary rate?     ####
#
# footprint = (TJ per tonne, mass-weighted across every item embodied in this
# country's domestic-protein consumption) x (tonnes of embodied mass itself).
# A RATE outlier would point at the convert_mass_vecs()/exio_mass_x bridging
# step (99.utils.R convert_intensities(), 1.mrio_convert.R convert_mass_vecs())
# -- i.e. TWN's EXIO "food" sectors reporting monetary output too small relative
# to the physical FABIO tonnage routed through them, inflating the per-tonne
# intensity. A VOLUME outlier instead says TWN's own non-food industries are
# genuinely working with (or processing/re-exporting under the "domestic
# consumption" label) an unusually large physical food mass per capita --
# a supply-chain-structure story, not a bridging artifact.
#
# Needs X_dom (embodied mass, any producer country, feeding this country's
# domestically-traded food consumption -- 23001 x 187, built in
# 2.analyze_result.R section 1.3 alongside fp_domcons/fp_impcons). If it isn't
# in the environment (e.g. fp_domcons/fp_impcons were loaded from the saved
# .rds without re-running section 1.3), rebuild it from Y_dom, which is cheap
# once FABIO_L and Y_dom already exist.
if (!exists("X_dom")) {
  stopifnot(exists("FABIO_L"), exists("Y_dom"))
  X_dom = FABIO_L %*% Y_dom
}

compare_rate_volume <- function(country_iso) {
  i = which(regions$iso3c == country_iso)
  own_region = FABIO_reg$EXIOBASE_code[FABIO_reg$ISO == country_iso]
  own_rows = own_rows_fun(own_region)
  pop = region_population[own_region]

  mass_i = as.vector(X_dom[, i])                       # tonnes, any producer, embodied in domestic-protein consumption
  total_mass = sum(mass_i)
  footprint = sum(fp_domcons$nonfood$en[own_rows, i])  # TJ, domestic-effort energy, domestic-protein bucket
  rate = footprint / total_mass                        # TJ per tonne, mass-weighted average

  cat(sprintf("%-4s: mass = %8.1f kg/cap/day | rate = %8.4g TJ/tonne | footprint = %6.2f MJ/cap/day\n",
              country_iso, total_mass / pop * 1e6 / 365 * 1000, rate, footprint / pop * 1e6 / 365))
  invisible(list(total_mass = total_mass, rate = rate, footprint = footprint, pop = pop))
}

cat("\n==== Rate x volume decomposition: domestic-protein bucket, domestic-effort energy ====\n")
rv_twn = compare_rate_volume("TWN")
rv_kor = compare_rate_volume("KOR")
rv_usa = compare_rate_volume("USA")
rv_chn = compare_rate_volume("CHN")

cat(sprintf("\nTWN / KOR -- mass ratio: %.2fx, rate ratio: %.2fx (product ~= footprint ratio %.2fx)\n",
            (rv_twn$total_mass / rv_twn$pop) / (rv_kor$total_mass / rv_kor$pop),
            rv_twn$rate / rv_kor$rate,
            (rv_twn$footprint / rv_twn$pop) / (rv_kor$footprint / rv_kor$pop)))
cat(sprintf("TWN / USA -- mass ratio: %.2fx, rate ratio: %.2fx (product ~= footprint ratio %.2fx)\n",
            (rv_twn$total_mass / rv_twn$pop) / (rv_usa$total_mass / rv_usa$pop),
            rv_twn$rate / rv_usa$rate,
            (rv_twn$footprint / rv_twn$pop) / (rv_usa$footprint / rv_usa$pop)))

cat("\nRead this as: if the MASS ratio alone is already close to the footprint\n")
cat("ratio, TWN's own consumption/production of embodied food mass is doing\n")
cat("most of the work (a volume story, e.g. re-processing/re-export activity\n")
cat("counted as domestic consumption). If the RATE ratio is doing most of the\n")
cat("work instead, the per-tonne intensity itself is inflated for TWN --\n")
cat("go back to exio_mass_x and EXIO_x_food (1.1.mrio_convert_indirect.R) for\n")
cat("TWN's own EXIO region (41) and check whether its food-sector monetary\n")
cat("output is implausibly small relative to the FABIO mass mapped onto it.\n")

# RESULT (first run): mass ratio TWN/KOR = 0.83x, TWN/USA = 0.62x -- TWN
# actually moves LESS embodied food mass per capita than either comparator.
# The rate ratio (11.71x vs KOR, 6.67x vs USA) accounts for essentially all of
# the footprint gap. So this is a RATE outlier, not a volume story -- go to
# stage 1 of the bridging (convert_mass_vecs()'s exio_mass_x, before the
# FABIO-item reallocation in l_int_i) to see whether the inflation is already
# present there.


#### Stage-1 check: is exio_mass_x (mass mapped onto TWN's own EXIO food-  ####
#### sector columns) too small relative to the embodied energy the        ####
#### Leontief inverse attributes to those columns' economic output?       ####
#
# total_intensity_exio_by_mass (1.1.mrio_convert_indirect.R) = indir_sat_exio
# / exio_mass_x, computed at EXIO-sector resolution (9800 columns) BEFORE
# l_int_i reallocates it to FABIO items via FABIO_x_in_EXIO. Recomputing just
# this ratio for each country's own food-sector columns isolates whether the
# rate anomaly already exists at this first bridging step, or only appears
# after the second (FABIO-item) reallocation.
#
# idx_food is a 200-sector pattern (same food sectors flagged in every EXIO
# region, since prod_map's food classification isn't region-specific) that
# R recycles across the 9800-length EXIO_x/exio_mass_x/indir_sat_exio columns
# (49 regions x 200 sectors, region-major order) -- see 0.mrio_prep.R.
stopifnot(exists("indir_sat_exio"), exists("exio_mass_x"), exists("idx_food"))

food_cols_of_region <- function(region_code) ((region_code - 1) * 200) + which(idx_food)

stage1_check <- function(country_iso) {
  region_code = FABIO_reg$EXIOBASE_code[FABIO_reg$ISO == country_iso]
  cols = food_cols_of_region(region_code)

  mass_tonnes = exio_mass_x[cols]
  embodied_TJ = Matrix::colSums(indir_sat_exio$sat_en[, cols, drop = FALSE])
  rate = embodied_TJ / mass_tonnes
  rate[!is.finite(rate)] = 0

  df = data.frame(sector = EXIO_reg$sector[which(idx_food)],
                   mass_tonnes = mass_tonnes, embodied_TJ = embodied_TJ, rate_TJ_per_tonne = rate) %>%
    arrange(desc(embodied_TJ))

  cat("\n====", country_iso, "-- EXIO food-sector columns, stage-1 bridging ====\n")
  cat(sprintf("Total mass: %.4g tonnes | Total embodied TJ: %.4g | Mass-weighted rate: %.4g TJ/tonne\n",
              sum(mass_tonnes), sum(embodied_TJ), sum(embodied_TJ) / sum(mass_tonnes)))
  print(head(df %>% mutate(across(where(is.numeric), ~signif(., 4))), 10))
  invisible(df)
}

s1_twn = stage1_check("TWN")
s1_kor = stage1_check("KOR")
s1_usa = stage1_check("USA")

cat("\nIf TWN's stage-1 mass-weighted rate is already ~10x KOR's/~7x USA's here,\n")
cat("the anomaly is in convert_mass_vecs()'s monetary-share bridging itself\n")
cat("(exio_mass_x too small for TWN's food columns, or EXIO's own reported\n")
cat("monetary output for those columns is miscalibrated for Taiwan). If the\n")
cat("stage-1 rates look comparable across countries and the gap only opens up\n")
cat("in l_int_i's later FABIO-item reallocation, look at FABIO_x_in_EXIO's\n")
cat("split specifically for TWN's items instead (1.1.mrio_convert_indirect.R,\n")
cat("the FP_trans = t(FABIO_x_in_EXIO %*% t(d)) step).\n")
