#### Verify: has the RoW-duplication fix actually landed, and what's the RoW ####
#### share of "import effort" now (post-fix)?                                  ####
#
# This script originally checked whether two FABIO countries sharing a RoW
# EXIO region had byte-identical footprint rows -- the smoking gun for the
# duplication bug fixed in 1.1.mrio_convert_indirect.R / 2.analyze_result.R
# (see chat: l_int_i used to paste each EXIO region's row-block onto every
# FABIO country mapped to it, so summing origin rows -- as
# agg_country_footprint()/country_summary()'s colSums(mat)-diag(mat) did --
# counted a RoW region's true total once PER MEMBER COUNTRY instead of once).
#
# That check no longer applies: l_int_i's origin/row axis is now natively
# EXIO-region-resolved (49 rows, not 187 FABIO-country-pasted rows), so
# there's only ONE row per region to begin with -- nothing left to compare
# for duplication. This version instead (a) confirms the data actually has
# the new shape, and (b) re-runs the "RoW share of import effort" diagnostic
# that originally exposed the bug (83-93% RoW share), which should now come
# out far lower and more economically plausible.

library(tidyverse)
library(Matrix)
library(data.table)

year = 2020
FABIO_path = "H:/MyDocuments/Data/FABIO/input/"

regions <- fread(file.path(FABIO_path, "regions.csv"))
nrreg <- nrow(regions)

reg_map <- readxl::read_xlsx("data/fabio-exiobase.xlsx", sheet = "regions_concordance", col_names = TRUE)
FABIO_reg <- readxl::read_xlsx(file.path(FABIO_path, "fabio_classifications_v2.xlsx"), sheet = "Countries") %>%
  select(-area) %>% rename(ISO = iso3c, FAO_code = area_code) %>% left_join(reg_map) %>%
  mutate(EXIOBASE_code = as.numeric(ifelse(EXIOBASE_code == "NA", 47, EXIOBASE_code)),
         EXIOBASE = ifelse(EXIOBASE == "NA", "RoW Europe", EXIOBASE))

# Data-quality fix: code 41 (TWN's own individually-modeled EXIO region) is
# mislabeled "RoW Asia and Pacific" in the source concordance -- see
# 7.import_effort_hypothesis_test.R for the diagnosis. Correct its label.
stopifnot(identical(sort(unique(FABIO_reg$EXIOBASE[FABIO_reg$EXIOBASE_code == 41])), "RoW Asia and Pacific"))
FABIO_reg$EXIOBASE[FABIO_reg$EXIOBASE_code == 41] <- "Taiwan"

stopifnot(identical(regions$iso3c, FABIO_reg$ISO))
row_countries <- FABIO_reg$ISO[grepl("RoW", FABIO_reg$EXIOBASE)]

n_reg_EXIO <- length(unique(FABIO_reg$EXIOBASE_code))
region_row_of_country <- FABIO_reg$EXIOBASE_code
region_members <- split(seq_len(nrreg), region_row_of_country)  # region index (as character) -> FABIO column indices
region_name_of_index <- FABIO_reg$EXIOBASE[match(seq_len(n_reg_EXIO), FABIO_reg$EXIOBASE_code)]

cat("\nLoading fp_impcons (~95MB)...\n")
fp_impcons <- readRDS(paste0("data/fp_impcons_", year, ".rds"))
mat <- fp_impcons$nonfood$hr_m   # should now be (n_reg_EXIO * n_nf) x 187, NOT (nrreg * n_nf) x 187

n_nf <- nrow(mat) / n_reg_EXIO
cat(sprintf("\nnrow(fp_impcons$nonfood$hr_m) = %d ; n_reg_EXIO = %d ; implied sectors/region = %.4f\n",
           nrow(mat), n_reg_EXIO, n_nf))
stopifnot(n_nf == round(n_nf))
cat("OK -- data has the fixed, EXIO-region-resolved shape: one row-block per region, not per FABIO country.\n")
cat(sprintf("(For reference, the OLD/buggy shape would have had nrow = %d = nrreg(%d) * %d sectors.)\n",
           nrreg * n_nf, nrreg, n_nf))

agg_exio_region_footprint <- function(mat) {
  mat_region = matrix(0, nrow = n_reg_EXIO, ncol = nrreg)
  nsect = nrow(mat) / n_reg_EXIO
  for (i in 1:n_reg_EXIO) mat_region[i, ] = colSums(mat[((i - 1) * nsect + 1):(i * nsect), ])
  rownames(mat_region) = region_name_of_index
  colnames(mat_region) = regions$iso3c
  mat_region
}

hr_nonfood_region = agg_exio_region_footprint(fp_impcons$nonfood$hr_m) + agg_exio_region_footprint(fp_impcons$nonfood$hr_f)  # 49 x 187, M.hr

row_idx = which(grepl("RoW", rownames(hr_nonfood_region)))

# For each importer, what share of its non-food import-effort hours comes
# from RoW-aggregate regions vs. named (individually-modeled) regions? This
# is the SAME diagnostic that originally found 83-93% -- now computed on the
# fixed (non-duplicated) data.
domestic_vec = hr_nonfood_region[cbind(region_row_of_country, seq_len(nrreg))]
hr_import_total_vec = colSums(hr_nonfood_region) - domestic_vec
row_hr_share_vec = colSums(hr_nonfood_region[row_idx, , drop = FALSE]) / hr_import_total_vec

df_row_share = data.frame(country = regions$iso3c, row_hr_share = row_hr_share_vec) %>%
  filter(!country %in% row_countries, is.finite(row_hr_share))

cat("\n---- Post-fix: RoW-aggregate share of non-food import-effort hours ----\n")
cat(sprintf("  median: %.1f%%  (range %.1f - %.1f%%)\n",
           median(df_row_share$row_hr_share) * 100,
           min(df_row_share$row_hr_share) * 100, max(df_row_share$row_hr_share) * 100))
cat("  (Pre-fix figure, for comparison: median 88.4%, range 83.2 - 93.2%.)\n")
print(df_row_share %>% arrange(desc(row_hr_share)) %>% head(10) %>% mutate(row_hr_share = round(row_hr_share * 100, 1)))
