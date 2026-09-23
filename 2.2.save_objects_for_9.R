#### Save everything the 9* analysis scripts need from 0-2, so they can run without re-running 2 ####
#
# Run once, at the end of a session in which 0.mrio_prep.R, 1*.R and
# 2.analyze_result.R have been run. Afterwards, any 9* script can start from a
# fresh R session with:
#
#   setwd("h:/MyDocuments/Projects/REMASS/food time")
#   load(paste0("data/objects_for_9_", 2020, ".RData"))
#
# (9x does this itself when the objects are missing.)
#
# The list comes from a static scan of the 9* scripts for names they use but
# define upstream. Not covered: 9r.diagnose_twn_nonfood_outlier.R, a one-off
# diagnostic that needs the full MRIO matrices (FABIO_L, X_dom, Y_dom, ...);
# run it after 2.analyze_result.R as before.

objects_for_9 <- c(
  # 0.mrio_prep.R: settings and lookups
  "year", "yr", "type", "FABIO_reg", "countries", "cty_ghd", "countrypops",
  # 2.analyze_result.R: analysis tables
  "tradeoff_protein_allwork_consump",     # 9, 9b, 9k, 9m, 9p, 9r.mediation, 9s, ...
  "tradeoff_protein_allwork",             # 9r.mediation_domestic_consumption_cf
  "tradeoff_protein_econlabor_consump",   # 9t (CF_paid)
  "summary_pro_df_long",                  # 9e, 9q
  "summary_food_df_long",
  "summary_food_df_long_with_ghd",        # 9u
  "effort_consumption_df",                # 9e, 9q, 9s
  "region_to_iso_1to1",                   # 9e, 9q, 9s
  "row_countries",                        # 9u
  "pop_data_yr",                          # RoW-collapsed samples (9, 9b, 9p, 9s, 9t)
  "mosaic_fill_levels"                    # 9q
)

missing <- objects_for_9[!vapply(objects_for_9, exists, logical(1))]
if (length(missing)) stop("Not in the session (run 0-2 first): ", paste(missing, collapse = ", "))

out <- paste0("data/objects_for_9_", year, ".RData")
save(list = objects_for_9, file = out, compress = "xz")
cat(sprintf("Saved %d objects to %s (%.1f MB)\n", length(objects_for_9), out, file.size(out) / 2^20))
