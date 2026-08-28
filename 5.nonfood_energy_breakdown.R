#### Ad-hoc analysis: decompose the DOMESTIC non-food-sector energy footprint ####
#### by EXIO non-food sector and by FABIO commodity, for a given country.     ####
#
# Standalone script (run with the project root as working directory) — reuses
# only what's needed from 0.mrio_prep.R / 1.mrio_convert.R / 1.1.mrio_convert_indirect.R
# rather than sourcing the full pipeline, since l_int_i is already saved to disk.
#
# "Domestic" here means the diagonal convention used throughout 2.analyze_result.R
# (agg_country_footprint()/country_summary()): producer country == consumer country.

library(tidyverse)
library(Matrix)
library(data.table)

year = yr = 2020
type = "pxp"
FABIO_path = "H:/MyDocuments/Data/FABIO/input/"
EXIO_path = paste0("H:/MyDocuments/Data/EXIOBASE3/IOT_", year, "_", type)

cat("Loading FABIO core objects...\n")
FABIO_x = readRDS(file.path(FABIO_path,"X.rds"))[,as.character(year)]
FABIO_y = readRDS(file.path(FABIO_path,"Y.rds"))[[as.character(year)]]
FABIO_L = readRDS(file.path(FABIO_path, paste0(year, "_L_value.rds")))
FABIO_y_hh = FABIO_y[,grep("food", colnames(FABIO_y))]

regions <- fread(file=file.path(FABIO_path,"regions.csv"))
items <- fread(file=file.path(FABIO_path,"items.csv"))
nrreg <- nrow(regions)
nrcom <- nrow(items)
io <- fread(file.path(FABIO_path,"io_labels.csv"))

reg_map = readxl::read_xlsx("data/fabio-exiobase.xlsx", sheet = "regions_concordance", col_names = TRUE)
FABIO_reg = readxl::read_xlsx(paste0(FABIO_path, "fabio_classifications_v2.xlsx"), sheet = "Countries") %>%
  select(-area) %>% rename(ISO = `iso3c`, FAO_code = `area_code`) %>% left_join(reg_map) %>%
  mutate(EXIOBASE_code = as.numeric(ifelse(EXIOBASE_code=="NA", 47, EXIOBASE_code)),
         EXIOBASE = ifelse(EXIOBASE=="NA", "RoW Europe", EXIOBASE))

# Sanity check: regions/FABIO_reg/io must all share the same 187-country ordering
# for the block-index arithmetic below (item_rows_i, nf_row_block) to line up.
stopifnot(identical(regions$iso3c, FABIO_reg$ISO))

# Nutrient coefficients (protein), same construction as 1.mrio_convert.R
coeff_pro <- fread(file.path(FABIO_path,"nutrient_coefficients_protein.csv")) %>%
  mutate(protein_per_kg = 10*protein_per_100g) %>% select(protein_per_kg) %>%
  slice(rep(1:n(), times = nrow(FABIO_reg))) %>%
  mutate(protein_per_kg = ifelse(items$comm_group[match(io$item, items$item)] == "Live animals", 0, protein_per_kg))

# EXIO region/sector info, needed only for non-food sector name labels
EXIO_reg = data.table::fread(paste0(EXIO_path, "/unit.txt"), header = T)
n_reg_EXIO = length(unique(EXIO_reg$region))
EXIO_sect = unique(EXIO_reg$sector)

prod_map = readxl::read_xlsx("data/fabio-exiobase-v2.xlsx",
                             sheet = "product_concordance",
                             range = "E4:GV126",
                             col_names = FALSE) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

i_exio_bio_sectors = which(colSums(prod_map) > 0)
bio_nonfood = c("Wool, silk-worm cocoons", "Chemicals nec", "Plant-based fibers")
i_exio_food_sectors = setdiff(i_exio_bio_sectors, which(EXIO_sect %in% bio_nonfood))
exio_nonfood_sectors = EXIO_reg$sector[setdiff(1:200, i_exio_food_sectors)]
n_nf = length(exio_nonfood_sectors)

cat("Loading l_int_i (indirect non-food satellites, per-tonne intensity)...\n")
l_int_i = readRDS(file = paste0("data/FABIO_exio_satellites_nonfood_", year, ".rds"))
stopifnot(dim(l_int_i$en) == c(nrreg*n_nf, nrow(io)))

# ---- Decompose the domestic non-food ENERGY footprint for one country ----
# l_int_i is a per-tonne intensity (TJ/tonne). To get an absolute footprint we
# multiply by the actual embodied production mass (tonnes) feeding that
# country's own consumption, restricted to items it produces itself (the
# "domestic" diagonal slice), then sum across items (-> by sector) or across
# sectors (-> by commodity).
compute_domestic_breakdown <- function(country_iso) {
  i = which(regions$iso3c == country_iso)
  if (length(i) == 0) stop("country not found in regions")

  # Embodied production (tonnes) feeding this country's total food consumption
  X_i = as.vector(FABIO_L %*% FABIO_y_hh[, i])

  item_rows_i = ((i-1)*nrcom + 1):(i*nrcom)          # FABIO rows produced by country i
  X_i_dom_items = X_i[item_rows_i]                    # tonnes of i's own items feeding i's own consumption

  # Protein-equivalent of that domestic mass (FABIO_x_hh_pro convention), used only
  # to express the total per gram of protein for context -- this is a cumulative,
  # production-weighted quantity (not final per-capita consumption), so treat the
  # resulting MJ/g ratio as a rough comparison, not a literal per-consumed-gram figure.
  dom_protein_g = sum(X_i_dom_items * coeff_pro$protein_per_kg[item_rows_i]) * 1000

  nf_row_block = ((i-1)*n_nf + 1):(i*n_nf)            # country i's own 175 non-food-sector rows
  A_i_dom = l_int_i$en[nf_row_block, item_rows_i, drop=FALSE]  # TJ/tonne, sector x item

  FP_dom = sweep(as.matrix(A_i_dom), 2, X_i_dom_items, `*`)    # TJ, sector x item
  total_TJ = sum(FP_dom)

  by_sector_df = data.frame(sector = exio_nonfood_sectors, TJ = rowSums(FP_dom)) %>%
    arrange(desc(TJ)) %>% mutate(share_pct = TJ/sum(TJ)*100)

  by_item_df = data.frame(item = items$item, TJ = colSums(FP_dom)) %>%
    arrange(desc(TJ)) %>% mutate(share_pct = TJ/sum(TJ)*100)

  cat("\n====", country_iso, "====\n")
  cat("Total domestic non-food energy footprint (TJ):", total_TJ, "\n")
  cat("MJ per g protein (domestic, non-food-sector energy only, rough):", total_TJ*1e6 / dom_protein_g, "\n")
  cat("\nTop 10 EXIO non-food sectors:\n"); print(head(by_sector_df, 10))
  cat("\nTop 10 FABIO commodities:\n"); print(head(by_item_df, 10))

  list(total_TJ = total_TJ, dom_protein_g = dom_protein_g, by_sector = by_sector_df, by_item = by_item_df)
}

res_usa = compute_domestic_breakdown("USA")
res_chn = compute_domestic_breakdown("CHN")

saveRDS(list(usa = res_usa, chn = res_chn), file.path("results", paste0("nonfood_energy_breakdown_domestic_", year, ".rds")))
