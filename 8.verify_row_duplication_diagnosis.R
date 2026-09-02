#### Verify: are l_int_i's per-FABIO-country non-food rows literally duplicated ####
#### for countries sharing a RoW EXIO-region mapping, and does naive row-summing ####
#### (agg_country_footprint / country_summary's colSums(mat)-diag(mat)) inflate  ####
#### the "import effort" total by the region's member count?                    ####
#
# Cheap check -- reuses fp_impcons_<year>.rds (95MB, already used in
# 7.import_effort_hypothesis_test.R) instead of the full 3.7GB l_int_i file:
# if row-duplication is real, it should already be visible at the footprint
# level (fp_impcons$nonfood = l_int_i %*% X_imp), since d[Somalia,:] ==
# d[Ethiopia,:] (identical vectors) implies d[Somalia,:] %*% X[,j] ==
# d[Ethiopia,:] %*% X[,j] for ANY X -- no need to touch l_int_i itself.

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
stopifnot(identical(regions$iso3c, FABIO_reg$ISO))

# Region membership counts alone predict the inflation factor, if the
# duplication mechanism is real.
region_counts <- FABIO_reg %>% count(EXIOBASE, sort = TRUE)
cat("---- EXIO region membership counts (RoW aggregates only) ----\n")
print(region_counts %>% filter(grepl("RoW", EXIOBASE)))

cat("\nLoading fp_impcons (~95MB)...\n")
fp_impcons <- readRDS(paste0("data/fp_impcons_", year, ".rds"))
mat <- fp_impcons$nonfood$hr_m   # 32725 x 187, rows = origin country x non-food sector (blocks of n_nf per country)

n_nf <- nrow(mat) / nrreg
stopifnot(n_nf == round(n_nf))

country_block <- function(iso) {
  i <- match(iso, regions$iso3c)
  ((i - 1) * n_nf + 1):(i * n_nf)
}

# Pick one RoW region and one real importer to test the mechanism against.
test_region   <- "RoW Africa"
test_importer <- "USA"

members <- FABIO_reg$ISO[FABIO_reg$EXIOBASE == test_region]
j <- match(test_importer, regions$iso3c)
cat(sprintf("\n%s has %d member FABIO countries: %s%s\n",
           test_region, length(members), paste(head(members, 10), collapse = ", "),
           ifelse(length(members) > 10, ", ...", "")))

member_contrib <- sapply(members, function(iso) sum(mat[country_block(iso), j]))
cat(sprintf("\nEach %s member's OWN row-block contribution to %s's non-food import-effort hours\n(should be IDENTICAL across members if rows are duplicated):\n",
           test_region, test_importer))
print(round(member_contrib, 6))

all_identical <- isTRUE(all.equal(member_contrib, rep(member_contrib[1], length(member_contrib)), tolerance = 1e-8))
cat(sprintf("\nAll members identical? %s\n", all_identical))

naive_total <- sum(member_contrib)     # what colSums(mat)-style aggregation actually computes
true_value  <- member_contrib[1]       # the single underlying (shared) regional value

cat(sprintf("\nNaive summed total across all %d members (what the current pipeline computes): %.4f M.hr\n",
           length(members), naive_total))
cat(sprintf("Single deduplicated regional value (any one member's row-block):                 %.4f M.hr\n", true_value))
cat(sprintf("Inflation factor: %.2fx  (member count = %d)\n", naive_total / true_value, length(members)))
