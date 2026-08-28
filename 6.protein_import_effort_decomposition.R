#### Ad-hoc analysis: decompose the IMPORT-EFFORT time embodied in a country's ####
#### domestic protein-consumption bucket, by origin country/EXIO region and   ####
#### by commodity/sector.                                                     ####
#
# Standalone script (run with the project root as working directory) — reuses
# only lightweight metadata (regions/items/EXIO sector lists/FABIO_reg) plus
# the pre-computed fp_domcons_<year>.rds saved by 2.analyze_result.R section
# 1.3, rather than re-running the full pipeline.
#
# "Domestic protein-consumption bucket" = fp_domcons = compute_footprints(X_dom),
# i.e. the footprint of protein that a country sources through domestically-
# traded final goods (per FABIO_y_hh's own domestic/import split), NOT protein
# imported as a finished product. Within that bucket, "import effort" is labor
# that happened in a DIFFERENT country than the consumer (off-diagonal in the
# country x country footprint matrix) -- e.g. imported packaging materials,
# machinery, or trade/transport/other services embodied in getting that
# domestically-sourced food to the table. This is exactly what backs the
# "Import effort" rectangles in the mosaic plot at 2.analyze_result.R ~L1420.

library(tidyverse)
library(Matrix)
library(data.table)
library(gt)
data(countrypops)

year = 2020
yr = 2020
FABIO_path = "H:/MyDocuments/Data/FABIO/input/"
EXIO_path  = paste0("H:/MyDocuments/Data/EXIOBASE3/IOT_", year, "_pxp")

regions <- fread(file.path(FABIO_path, "regions.csv"))
items   <- fread(file.path(FABIO_path, "items.csv"))
nrreg <- nrow(regions)
nrcom <- nrow(items)

EXIO_reg <- fread(file.path(EXIO_path, "unit.txt"), header = TRUE)
EXIO_sect <- unique(EXIO_reg$sector)

prod_map <- readxl::read_xlsx("data/fabio-exiobase-v2.xlsx",
                              sheet = "product_concordance",
                              range = "E4:GV126",
                              col_names = FALSE) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

i_exio_bio_sectors <- which(colSums(prod_map) > 0)
bio_nonfood <- c("Wool, silk-worm cocoons", "Chemicals nec", "Plant-based fibers")
i_exio_food_sectors <- setdiff(i_exio_bio_sectors, which(EXIO_sect %in% bio_nonfood))

exio_nonfood_sectors <- EXIO_reg$sector[setdiff(1:200, i_exio_food_sectors)]
n_nf <- length(exio_nonfood_sectors)
stopifnot(n_nf == 175)

pop_y <- subset(countrypops, year == yr) %>%
  select(iso3c = country_code_3, pop = population)

reg_map <- readxl::read_xlsx("data/fabio-exiobase.xlsx",
                             sheet = "regions_concordance",
                             col_names = TRUE)
FABIO_reg <- readxl::read_xlsx(paste0(FABIO_path, "fabio_classifications_v2.xlsx"),
                               sheet = "Countries") %>% select(-area) %>%
  rename(ISO = `iso3c`, FAO_code = `area_code`) %>%
  left_join(reg_map) %>%
  mutate(EXIOBASE_code = as.numeric(ifelse(EXIOBASE_code=="NA", 47, EXIOBASE_code)),
         EXIOBASE = ifelse(EXIOBASE=="NA", "RoW Europe", EXIOBASE))

# Sanity check: regions/FABIO_reg must share the same 187-country ordering for
# the block-index arithmetic below (producer_country_food, origin_country_nf) to
# line up with fp_domcons' rows.
stopifnot(identical(regions$iso3c, FABIO_reg$ISO))
iso_to_exioregion <- setNames(FABIO_reg$EXIOBASE, FABIO_reg$ISO)

cat("Loading fp_domcons (large file)...\n")
fp_domcons <- readRDS(paste0("data/fp_domcons_", year, ".rds"))

producer_country_food <- rep(regions$iso3c, each = nrcom)
item_food             <- rep(items$item, times = nrreg)

origin_country_nf <- rep(regions$iso3c, each = n_nf)
sector_nf         <- rep(exio_nonfood_sectors, times = nrreg)

# For country `cty`: split fp_domcons (food + non-food, both genders summed) into
# domestic-effort (origin == cty) vs import-effort (origin != cty) rows, then
# break the import-effort part down by origin country, by EXIOBASE region
# (the real underlying resolution -- RoW regions paste one value to all FABIO
# countries within them, so raw per-country rankings are mostly a resolution
# artifact), and by commodity/sector.
decompose_import_effort <- function(cty, bucket_list = fp_domcons) {
  col_idx <- match(cty, regions$iso3c)
  pop_cty <- pop_y$pop[match(cty, pop_y$iso3c)]
  if (is.na(col_idx)) stop("country not found in regions")
  if (is.na(pop_cty)) stop("country not found in countrypops for year ", yr)

  food_hr <- bucket_list$food[["hr_m"]][, col_idx] + bucket_list$food[["hr_f"]][, col_idx]
  nf_hr   <- bucket_list$nonfood[["hr_m"]][, col_idx] + bucket_list$nonfood[["hr_f"]][, col_idx]

  df_food <- data.frame(sector = "food", origin = producer_country_food,
                        item = item_food, mhr = as.numeric(food_hr))
  df_nf   <- data.frame(sector = "non-food", origin = origin_country_nf,
                        item = sector_nf, mhr = as.numeric(nf_hr))

  df <- bind_rows(df_food, df_nf) %>%
    mutate(effort = ifelse(origin == cty, "domestic", "import"),
           min_per_cap_day = mhr * 1e6 / pop_cty / 365 * 60)

  totals <- df %>% group_by(sector, effort) %>%
    summarise(min_per_cap_day = sum(min_per_cap_day), .groups = "drop")
  import_total = sum(totals$min_per_cap_day[totals$effort == "import"])

  cat("\n====", cty, "-- import effort within domestic protein-consumption bucket ====\n")
  print(totals)
  cat("TOTAL import-effort min/cap/day:", round(import_total, 1), "\n")
  cat("TOTAL (domestic+import) min/cap/day:", round(sum(totals$min_per_cap_day), 1), "\n")

  by_origin <- df %>% filter(effort == "import") %>%
    group_by(origin) %>%
    summarise(min_per_cap_day = sum(min_per_cap_day), .groups = "drop") %>%
    arrange(desc(min_per_cap_day)) %>%
    mutate(pct = min_per_cap_day / sum(min_per_cap_day) * 100)

  by_exio_region <- df %>% filter(effort == "import") %>%
    mutate(exio_region = iso_to_exioregion[origin]) %>%
    group_by(exio_region) %>%
    summarise(min_per_cap_day = sum(min_per_cap_day), .groups = "drop") %>%
    arrange(desc(min_per_cap_day)) %>%
    mutate(pct = min_per_cap_day / sum(min_per_cap_day) * 100)
  cat("\nTop EXIOBASE regions for import effort (real underlying resolution):\n")
  print(head(by_exio_region, 10))

  by_item <- df %>% filter(effort == "import") %>%
    group_by(sector, item) %>%
    summarise(min_per_cap_day = sum(min_per_cap_day), .groups = "drop") %>%
    arrange(desc(min_per_cap_day)) %>%
    mutate(pct = min_per_cap_day / sum(min_per_cap_day) * 100)
  cat("\nTop commodities/sectors for import effort:\n")
  print(head(by_item, 15))

  list(totals = totals, import_total_min_per_cap_day = import_total,
       by_origin = by_origin, by_exio_region = by_exio_region, by_item = by_item)
}

countries_of_interest <- c("FRA", "AUS", "USA", "IND", "IDN")
results_domcons <- lapply(countries_of_interest, decompose_import_effort) %>%
  setNames(countries_of_interest)

saveRDS(results_domcons,
       file.path("results", paste0("protein_import_effort_domcons_", year, ".rds")))
