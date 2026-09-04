#### Result analysis ####


#### 1. Energy and Labor footprint ####

# Load economic (direct food sector) and non-economic (indirect non-food sector) satellite data
l_int_d <- readRDS(file = paste0("data/FABIO_exio_satellites_food_", year, ".rds"))
l_int_i <- readRDS(file = paste0("data/FABIO_exio_satellites_nonfood_RoW_", year, ".rds"))

n_nf = length(exio_nonfood_sectors) # EXIO non-food sectors

# Compute food-sector (direct) and non-food-sector (indirect) footprints for
# a given embodied-production matrix X (23001 x 187, rows = producer
# country/product, columns = consuming country). Used both for the full
# FABIO_x_hh (below) and, further down, for the domestic-/import-consumption
# restricted variants (X_dom/X_imp) that split effort by whether it feeds
# domestically- or internationally-traded final food products.
compute_footprints <- function(X) {
  fp_food_X <- lapply(l_int_d, function(d) Matrix::Diagonal(x=d) %*% X)

  # Footprint for non-food sectors is calculated by multiplying the intensity matrix
  # (in the same dim as satellite (32725x23001)) with the X vector (23001x187) to get a matrix of dimension (32725x187), and then summing every 187 rows to get a matrix of dimension (187x187) where rows are origin countries and columns are target countries.
  #
  # BUG (found via code review, see branch discussion): the loop below does NOT do
  # the general multiplication the comment above describes. l_int_i[[k]] (= d) is a
  # per-tonne-of-commodity-r intensity, fully general across origin country of the
  # non-food effort (row, via EXIO_L's cross-country linkages) and producer country
  # of commodity r (baked into column r) -- see 1.1.mrio_convert_indirect.R. It is
  # structurally the same kind of object as l_int_d (used above via a plain
  # Diagonal(x=d) %*% X), and the correct general footprint is likewise just
  # `d %*% X` (rows = non-food effort's origin country, columns = consuming
  # country j, matching country_summary()'s domestic/export/import formula).
  #
  # Instead, the loop forces the row-block index i (origin region of the non-food
  # industry, via A_i <- d[block i,]) to also select the SAME i as the consuming
  # country (via B_i, built from X[,i] only -- never X[,k] for k != i), then
  # regroups the leftover axis by the FOOD COMMODITY'S PRODUCER country instead of
  # leaving it as the general consuming-country axis. Concretely, B_i just scatters
  # the vector X[,i] into a 23001x187 matrix, zero except in the column matching
  # each row's own producer country (a one-hot "group by producer country" trick),
  # so A_i %*% B_i computes, for consumer i only: [s, producer q] = sum over
  # commodity rows r produced in q of A_i[s,r] * X[r,i].
  #
  # Net effect on country_summary()'s domestic/export/import split of the
  # resulting matrix:
  #  - domestic (diag[i,i]): narrower than intended -- only captures the
  #    triple-match case (non-food origin = food producer = food consumer = i),
  #    not "non-food effort in i embodied in i's total consumption, any food
  #    producer" (likely understated).
  #  - export (rowSums - diag): for fixed consumer i, sums over producer q != i
  #    -- that's "foreign-grown food consumed at home," an import concept
  #    mislabeled as export.
  #  - import (colSums - diag): for fixed producer q, sums across DIFFERENT
  #    countries' own-consumption buckets that happen to use food from q -- not
  #    q's import of anything; mixes unrelated countries together.
  #
  # This does not affect the food-sector branch above (fp_food_X), which already
  # uses the full X and is fine. Fixed below by replacing the loop with the
  # direct, general `d %*% X` -- old code kept commented out for reference.
  #
  # fp_nonfood_X <- lapply(l_int_i, function(d) {
  #   fp_list <- vector("list", nrreg)
  #   for (i in 1:nrreg) {
  #     print(paste("Processing region", i, "of", nrreg))
  #     A_i <- d[((i-1)*n_nf + 1):(i*n_nf), ]
  #     # Build B_i as block diagonal: each column k gets X commodities for region k, source i
  #     B_i <- Matrix::Diagonal(x = as.vector(X[, i])) %*%
  #       kronecker(Matrix::Diagonal(nrreg), Matrix::Matrix(1, nrcom, 1))
  #     fp_list[[i]] <- A_i %*% B_i
  #   }
  #   do.call(rbind, fp_list) %>% as("CsparseMatrix")
  # })
  fp_nonfood_X <- lapply(l_int_i, function(d) (d %*% X) %>% as("CsparseMatrix"))

  list(food = fp_food_X, nonfood = fp_nonfood_X)
}

# Footprint summed at the FABIO country level
fp_res <- compute_footprints(FABIO_x_hh)
fp_food <- fp_res$food
fp_nonfood <- fp_res$nonfood

# Save the results
saveRDS(fp_food, file = paste0("data/footprint_food_", year, ".rds"))
saveRDS(fp_nonfood, file = paste0("data/footprint_nonfood_", year, ".rds"))



# Read the results back and do some validations
fp_food = readRDS(file = paste0("data/footprint_food_", year, ".rds"))
fp_nonfood = readRDS(file = paste0("data/footprint_nonfood_", year, ".rds"))


#### 1.1. Aggregate all at country level ####

# fp_food and fp_nonfood lists have three large matrices each, where rows are multiples of 187 (number of countries).
# For each matrix, make partial row sums by adding every 187 rows. This will give us a matrix of dimension (187 x 187), where rows are origin countries and columns are target countries.

agg_country_footprint <- function(mat) {
  mat_country = matrix(0, nrow=nrreg, ncol=nrreg)
  nsect = dim(mat)[1]/nrreg
  print(paste("Aggregating footprint matrix with", nsect, "sectors per country..."))

  for (i in 1:nrreg) {
    mat_country[i, ] = colSums(mat[((i-1)*nsect+1):(i*nsect), ])
  }
  rownames(mat_country) = colnames(mat_country) = regions$iso3c
  return(mat_country)
}

# fp_nonfood's rows are now EXIO-region-resolved (8575 = 49 regions * 175
# sectors), NOT FABIO-country-resolved (see 1.1.mrio_convert_indirect.R --
# l_int_i deliberately no longer pastes each region's row-block onto every
# member FABIO country, since that duplication inflated RoW regions'
# contribution by their member count whenever origin rows got summed --
# confirmed in 8.verify_row_duplication_diagnosis.R). agg_country_footprint()
# above would silently misbehave here (8575/187 isn't even an integer), so
# use this region-resolved counterpart instead: output is 49 (EXIO region)
# x 187 (FABIO consuming country), with EACH region counted exactly once.
agg_exio_region_footprint <- function(mat) {
  mat_region = matrix(0, nrow = n_reg_EXIO, ncol = nrreg)
  nsect = dim(mat)[1] / n_reg_EXIO
  stopifnot(nsect == round(nsect))
  print(paste("Aggregating non-food footprint matrix with", nsect, "sectors per EXIO region..."))

  for (i in 1:n_reg_EXIO) {
    mat_region[i, ] = colSums(mat[((i-1)*nsect+1):(i*nsect), ])
  }
  rownames(mat_region) = unique(EXIO_reg$region)  # 49 EXIO region codes, same order FABIO_reg$EXIOBASE_code indexes into
  colnames(mat_region) = regions$iso3c
  return(mat_region)
}

# Country -> region lookup: which of the 49 EXIO-region rows is FABIO country
# j's own origin (its 1:1 region if individually modeled, or its shared RoW
# aggregate's row if not) -- and the reverse: which FABIO countries (column
# indices) belong to each of the 49 regions.
region_row_of_country = FABIO_reg$EXIOBASE_code
stopifnot(identical(regions$iso3c, FABIO_reg$ISO))  # region_row_of_country must line up with mat_region's columns positionally
region_members = split(seq_len(nrreg), region_row_of_country)  # list keyed by region index (as character) -> column indices

region_pop = subset(countrypops, year == yr) %>% select(iso3c = country_code_3, population)
pop_by_country_idx = region_pop$population[match(regions$iso3c, region_pop$iso3c)]
region_population = sapply(seq_len(n_reg_EXIO), function(i) sum(pop_by_country_idx[region_members[[as.character(i)]]], na.rm = TRUE))
region_name_of_index = FABIO_reg$EXIOBASE[match(seq_len(n_reg_EXIO), FABIO_reg$EXIOBASE_code)]  # human-readable label per region

# ISO3 -> EXIO region name (for translating a list of FABIO countries into the
# deduplicated set of regions they belong to, e.g. when a hand-picked country
# list for some chart happens to include multiple RoW-mapped countries that
# share one region).
iso_to_region = setNames(FABIO_reg$EXIOBASE, FABIO_reg$ISO)
countries_to_regions = function(iso_codes) unique(iso_to_region[iso_codes])

# Reverse: EXIO region name -> ISO3, defined only for named (non-RoW) regions
# where the mapping is genuinely 1:1. Used when a curated country list for
# some chart is already known to be all-named (e.g. mosaic_isos below), so
# effort_consumption_df's exio_region can be translated straight back to the
# ISO3 codes the rest of that chart's code already expects.
region_to_iso_1to1 = FABIO_reg %>% filter(!grepl("RoW", EXIOBASE)) %>%
  distinct(EXIOBASE, ISO) %>% { setNames(.$ISO, .$EXIOBASE) }

# Continent label for each of the 49 EXIO regions -- for any continent-of-
# origin breakdown of non-food effort (which is region-resolved on its
# origin axis). A named (non-RoW) region's members are just itself, so this
# reduces to regions$continent. A RoW aggregate can span more than one FABIO
# continent label -- fall back to the region's own name, which is already a
# continent-scoped EXIOBASE grouping (RoW Africa/America/Asia and Pacific/
# Europe/Middle East). Used by both b_import (below) and make_continent_data().
region_continent_of_index = sapply(seq_len(n_reg_EXIO), function(i) {
  conts_i = unique(regions$continent[region_members[[as.character(i)]]])
  if (length(conts_i) == 1) conts_i else region_name_of_index[i]
})

# Collapse a matrix's row and/or column axis from FABIO-country resolution
# (187) to EXIO-region resolution (49), summing member countries within each
# region block. An axis already region-resolved (e.g. l_int_i's origin/row
# axis, collapsed at the source -- see 1.1.mrio_convert_indirect.R) is left
# as-is.
collapse_axis_to_region = function(mat, axis) {
  n = if (axis == 1) nrow(mat) else ncol(mat)
  if (n == n_reg_EXIO) return(mat)  # already region-resolved on this axis
  stopifnot(n == nrreg)
  if (axis == 1) {
    out = matrix(0, nrow = n_reg_EXIO, ncol = ncol(mat))
    for (i in seq_len(n_reg_EXIO)) out[i, ] = colSums(mat[region_members[[as.character(i)]], , drop = FALSE])
  } else {
    out = matrix(0, nrow = nrow(mat), ncol = n_reg_EXIO)
    for (i in seq_len(n_reg_EXIO)) out[, i] = rowSums(mat[, region_members[[as.character(i)]], drop = FALSE])
  }
  out
}

agg_to_region_matrix = function(mat) {
  m = collapse_axis_to_region(mat, axis = 1)
  m = collapse_axis_to_region(m, axis = 2)
  rownames(m) = colnames(m) = region_name_of_index
  m
}

# Region-level domestic/export/import summary, generalizing country_summary()
# to work identically for food (a 187x187 country matrix, both axes collapsed
# here) and non-food (l_int_i's already-49-row mat_region, only columns
# collapsed here) -- used ONLY for combined food+nonfood "total" presentation.
# summary_food (187-country, built separately below) stays the source of
# truth for food-only per-country displays, since food has no RoW ambiguity
# to begin with (see chat discussion: l_int_d's per-row diagonal multiply
# keeps each FABIO country's own production quantities distinct even though
# the underlying per-tonne RATE may be RoW-shared). For non-food, "RoW
# Africa's domestic" honestly means "all internal (within-region) activity
# across every member country" -- including e.g. Somalia's imports from
# Ethiopia -- and for the combined total, food is summed the same way (a
# legitimate aggregation of genuinely distinct country quantities, unlike the
# non-food duplication bug). No FABIO-country rows (SOM, ETH, ...) appear in
# either object; only the 49 EXIO regions.
region_summary = function(mat) {
  mat_r = agg_to_region_matrix(mat)
  df = data.frame(
    exio_region = rownames(mat_r),
    domestic = diag(mat_r),
    export = rowSums(mat_r) - diag(mat_r),
    import = colSums(mat_r) - diag(mat_r),
    population = region_population
  ) %>%
    mutate(domestic_per_capita = domestic / population * 1e6,
           export_per_capita = export / population * 1e6,
           import_per_capita = import / population * 1e6)
  return(df)
}

l_food_country = lapply(fp_food, agg_country_footprint)
l_nonfood_region = lapply(fp_nonfood, agg_exio_region_footprint)

# Validate: sum(l_food_country[[i]]) == sum(fp_food[[i]])

# # Validate: find the ten biggest cell from l_food_country[[2]] and check source and destination.
# find_top_cells(l_food_country[[2]], matrix_name = names(l_food_country)[2])
# find_top_cells(l_food_country[[3]], matrix_name = names(l_food_country)[3])
# find_top_cells(l_nonfood_country[[2]], matrix_name = names(l_nonfood_country)[2])
# find_top_cells(l_nonfood_country[[3]], matrix_name = names(l_nonfood_country)[3])



#### 1.2. Domestic/import/export summary by country ####

# For each of 187 countries, make summarized footprints of domestic, imported, and exported.
# Diagonal elements of the matrices give the domestic footprint, while off-diagonal row elements give the exported footprint for the row country, and off-diagonal column elements give the imported footprint for the column country.
country_summary <- function(mat) {
  df = data.frame(
    country = regions$iso3c,
    domestic = diag(mat),
    export = rowSums(mat) - diag(mat),
    import = colSums(mat) - diag(mat)
  )
  
  # Import population data and calculate per capita footprints
  pop_data = subset(countrypops, year == yr) %>% select(iso3c = country_code_3, population)
  df = df %>% left_join(pop_data, by = c("country" = "iso3c")) %>%
    # TJ to MJ, M.hr to hr, and then per capita (MJ/cap or hr/cap)
    mutate(domestic_per_capita = domestic / population *1e6,
           export_per_capita = export / population*1e6,
           import_per_capita = import / population*1e6) 
  
  return(df)
}

summary_food = lapply(l_food_country, country_summary)          # 187-country, full precision -- food-only displays
summary_food_region = lapply(l_food_country, region_summary)    # 49-region -- ONLY for combining with summary_nonfood below
summary_nonfood = lapply(l_nonfood_region, region_summary)       # 49-region (see region_summary()'s comment for why)

# Make per capita labor hour footprints daily by dividing by 365
for (i in names(summary_food[2:3])) {
  summary_food[[i]]        %>% mutate(across(ends_with("per_capita"), ~ .x / 365)) -> summary_food[[i]]
  summary_food_region[[i]] %>% mutate(across(ends_with("per_capita"), ~ .x / 365)) -> summary_food_region[[i]]
  summary_nonfood[[i]]     %>% mutate(across(ends_with("per_capita"), ~ .x / 365)) -> summary_nonfood[[i]]
}



#### 1.3. Domestic-effort vs. import-effort within domestic/import protein consumption ####

# summary_food/summary_nonfood above split effort by *where the work happened*
# (domestic_per_capita = producer country == consumer country), independent of
# whether the food product the work is embodied in was itself domestically-
# traded or imported (per FABIO_y_hh). That means "domestic" here doesn't line
# up with FABIO_y_hh's "domestic protein": e.g. labor spent domestically
# processing an imported raw ingredient counts as "domestic" above, blurring
# together two different things.
#
# Split FABIO_y_hh itself first, by whether the final good was domestically-
# traded (Y_dom) or imported (Y_imp) — its existing diagonal/off-diagonal
# structure, nothing new — then trace *each* piece separately through
# FABIO_L and the same intensities used above. This gives, for each country's
# domestically-traded food consumption and (separately) its imported food
# consumption, a further domestic-effort/import-effort split — matched
# exactly to FABIO_y_hh's own domestic/import split, since it's the same
# consumption slice on both sides. Because Y_dom + Y_imp = FABIO_y_hh exactly
# and FABIO_L is linear, the four resulting matrices sum back to
# l_food_country/l_nonfood_country exactly (see TEST below).
row_country_idx = rep(seq_len(nrreg), each = nrcom)  # producing-country index for each FABIO_y_hh row

y_hh_t = as(FABIO_y_hh, "TsparseMatrix")
is_domestic_flow = row_country_idx[y_hh_t@i + 1] == (y_hh_t@j + 1)

Y_dom = Matrix::sparseMatrix(i = y_hh_t@i[is_domestic_flow] + 1, j = y_hh_t@j[is_domestic_flow] + 1,
                             x = y_hh_t@x[is_domestic_flow], dims = dim(FABIO_y_hh))
Y_imp = Matrix::sparseMatrix(i = y_hh_t@i[!is_domestic_flow] + 1, j = y_hh_t@j[!is_domestic_flow] + 1,
                             x = y_hh_t@x[!is_domestic_flow], dims = dim(FABIO_y_hh))

X_dom = FABIO_L %*% Y_dom  # embodied production feeding each country's domestically-traded food consumption
X_imp = FABIO_L %*% Y_imp  # embodied production feeding each country's imported food consumption

fp_domcons = compute_footprints(X_dom)
fp_impcons = compute_footprints(X_imp)

# Save the results
saveRDS(fp_domcons, file = paste0("data/fp_domcons_", year, ".rds"))
saveRDS(fp_impcons, file = paste0("data/fp_impcons_", year, ".rds"))

# Save the results
fp_domcons = readRDS(file = paste0("data/fp_domcons_", year, ".rds"))
fp_impcons = readRDS(file = paste0("data/fp_impcons_", year, ".rds"))

l_food_country_domcons   = lapply(fp_domcons$food,    agg_country_footprint)
l_food_country_impcons   = lapply(fp_impcons$food,    agg_country_footprint)
l_nonfood_region_domcons = lapply(fp_domcons$nonfood, agg_exio_region_footprint)
l_nonfood_region_impcons = lapply(fp_impcons$nonfood, agg_exio_region_footprint)

# TEST: for each metric, l_food_country_domcons[[m]] + l_food_country_impcons[[m]] ~ l_food_country[[m]]
#       (and likewise for nonfood, at region resolution)

# effort_consumption_df (below) combines food+nonfood, so -- per chat
# discussion -- both go through region_summary() here (49-row, exio_region-
# keyed), the same as summary_food_region/summary_nonfood above. Nothing else
# in the file consumes summary_food_domcons/summary_food_impcons at 187-
# country resolution, so there's no separate country-level version to keep.
summary_food_domcons    = lapply(l_food_country_domcons,   region_summary)
summary_food_impcons    = lapply(l_food_country_impcons,   region_summary)
summary_nonfood_domcons = lapply(l_nonfood_region_domcons, region_summary)
summary_nonfood_impcons = lapply(l_nonfood_region_impcons, region_summary)

for (i in names(summary_food_domcons[2:3])) {
  summary_food_domcons[[i]]    %>% mutate(across(ends_with("per_capita"), ~ .x / 365)) -> summary_food_domcons[[i]]
  summary_food_impcons[[i]]    %>% mutate(across(ends_with("per_capita"), ~ .x / 365)) -> summary_food_impcons[[i]]
  summary_nonfood_domcons[[i]] %>% mutate(across(ends_with("per_capita"), ~ .x / 365)) -> summary_nonfood_domcons[[i]]
  summary_nonfood_impcons[[i]] %>% mutate(across(ends_with("per_capita"), ~ .x / 365)) -> summary_nonfood_impcons[[i]]
}

# Tidy 8-category table: {domestic, import} protein-consumption bucket x
# {food, non-food} sector x {domestic, import} effort-origin, per country x
# metric (hr_m/hr_f/en). Bar widths (protein supply) should use the existing
# y_hh-based domestic/import protein totals unchanged; this table only
# supplies the height-stacking within each protein bucket.
effort_consumption_df = bind_rows(
  bind_rows(lapply(names(summary_food_domcons), function(m)
    summary_food_domcons[[m]] %>% mutate(type = m, sector = "food", protein_source = "domestic"))),
  bind_rows(lapply(names(summary_food_impcons), function(m)
    summary_food_impcons[[m]] %>% mutate(type = m, sector = "food", protein_source = "import"))),
  bind_rows(lapply(names(summary_nonfood_domcons), function(m)
    summary_nonfood_domcons[[m]] %>% mutate(type = m, sector = "non-food", protein_source = "domestic"))),
  bind_rows(lapply(names(summary_nonfood_impcons), function(m)
    summary_nonfood_impcons[[m]] %>% mutate(type = m, sector = "non-food", protein_source = "import")))
) %>%
  select(exio_region, type, sector, protein_source, domestic_per_capita, import_per_capita) %>%
  pivot_longer(cols = c(domestic_per_capita, import_per_capita),
               names_to = "effort_origin", values_to = "per_capita_value") %>%
  mutate(effort_origin = ifelse(effort_origin == "domestic_per_capita", "domestic", "import"),
         is_row = grepl("RoW", exio_region))



# Stack vertically all elements in each list after adding a column "type" filled with the name the element
summary_food_df = bind_rows(lapply(names(summary_food),
                                   function(x) summary_food[[x]] %>% mutate(type = x))) %>%
    drop_na(domestic, export, import)
summary_nonfood_df = bind_rows(lapply(names(summary_nonfood),
                                      function(x) summary_nonfood[[x]] %>% mutate(type = x))) %>%
    drop_na(domestic, export, import)

# Flag countries mapped to aggregate "Rest of" EXIO regions (labour intensities not country-specific)
row_countries      <- FABIO_reg$ISO[grepl("RoW", FABIO_reg$EXIOBASE)]
summary_food_df    <- summary_food_df    %>% mutate(is_row = country %in% row_countries)
# summary_nonfood_df is region-resolved (region_summary()) -- every row IS either one named
# country's own EXIO region or one of the 5 RoW aggregates, so is_row is just whether the region
# name itself is a RoW aggregate (no per-country membership check needed/possible here).
summary_nonfood_df <- summary_nonfood_df %>% mutate(is_row = grepl("RoW", exio_region))

# Order of domestic hours (by female)
sum_ord = (summary_food_df %>%
    filter(type %in% c("hr_f")) %>%
    group_by(country) %>%
    summarise(d = sum(domestic_per_capita)) %>%
    arrange(-d))$country

# Same idea, but keyed by exio_region for the (49-row) non-food summary.
sum_ord_region = (summary_nonfood_df %>%
    filter(type %in% c("hr_f")) %>%
    group_by(exio_region) %>%
    summarise(d = sum(domestic_per_capita)) %>%
    arrange(-d))$exio_region

# Make summary_food_df and summary_nonfood_df in long format for plotting
summary_food_df_long = summary_food_df %>%
  select(-(domestic:population)) %>%
  pivot_longer(cols = c("domestic_per_capita", "export_per_capita", "import_per_capita"),
               names_to = "footprint_type", values_to = "per_capita_value") %>%
  # Order countries by sum of domestic_per_capita of type "hr_m" and "hr_f"
  mutate(country = factor(country, levels = sum_ord),
         footprint_type = factor(footprint_type,
                                 levels = c("preparation_non.econ", "processing_non.econ", "growth_collection_non.econ",
                                            "preparation_econ",
                                            "export_per_capita", "import_per_capita", "domestic_per_capita")))
summary_nonfood_df_long = summary_nonfood_df %>%
  select(-(domestic:population)) %>%
  pivot_longer(cols = c("domestic_per_capita", "export_per_capita", "import_per_capita"),
               names_to = "footprint_type", values_to = "per_capita_value") %>%
  mutate(exio_region = factor(exio_region, levels = sum_ord_region),
         footprint_type = factor(footprint_type, levels = c("export_per_capita", "import_per_capita", "domestic_per_capita")))

# Region-resolved (49-row) counterpart of summary_food_df/summary_food_df_long,
# built from summary_food_region -- exists ONLY so it can be bind_rows'd with
# summary_nonfood_df_long for combined food+nonfood totals (see chat: summing
# every RoW-Africa member's own, individually-correct food footprint into one
# "RoW Africa" row is a legitimate aggregation, unlike the non-food
# double-counting bug). summary_food/summary_food_df/summary_food_df_long
# (187-country) remain the source of truth for food-only per-country displays
# and are untouched by this.
summary_food_region_df = bind_rows(lapply(names(summary_food_region),
                                          function(x) summary_food_region[[x]] %>% mutate(type = x))) %>%
    drop_na(domestic, export, import) %>%
    mutate(is_row = grepl("RoW", exio_region))

summary_food_region_df_long = summary_food_region_df %>%
  select(-(domestic:population)) %>%
  pivot_longer(cols = c("domestic_per_capita", "export_per_capita", "import_per_capita"),
               names_to = "footprint_type", values_to = "per_capita_value") %>%
  # Same region ordering as summary_nonfood_df_long, so the two line up when combined.
  mutate(exio_region = factor(exio_region, levels = sum_ord_region),
         footprint_type = factor(footprint_type, levels = c("export_per_capita", "import_per_capita", "domestic_per_capita")))


# TEST: compare food-sector (direct) energy footprint from MRIO (summary_food[["en"]], MJ/cap/year)
# against IEA household food-related appliance energy (df_iea_breakdown, GJ/cap/year), both
# converted to MJ/cap/day ####
test_energy_mrio_iea = summary_food[["en"]] %>%
  select(country, mrio_domestic_mj_cap_day = domestic_per_capita,
         mrio_export_mj_cap_day = export_per_capita,
         mrio_import_mj_cap_day = import_per_capita) %>%
  mutate(across(starts_with("mrio_"), ~ .x / 365)) %>%  # MJ/cap/year -> MJ/cap/day
  inner_join(
    df_iea_breakdown %>%
      mutate(across(c(COOKING, DISH_WASH, REFRIG_ALL, TOTAL), ~ .x * 1000 / 365,  # GJ/cap/year -> MJ/cap/day
                    .names = "iea_{.col}_mj_cap_day")) %>%
      select(country, starts_with("iea_")),
    by = "country"
  )


# Vertically stack df_ghd_gender to summary_food_df_long, and then plot again with the same function plot_countries. This will add the non-economic food time to the existing labor footprint
summary_food_df_long_with_ghd = bind_rows(summary_food_df_long, df_ghd_combined) %>%
  arrange(country, type, footprint_type) %>%
  mutate(is_row = as.character(country) %in% row_countries,
         country = factor(country, levels = sum_ord),
         footprint_type = factor(footprint_type,
                                 levels = c("preparation_non.econ", "processing_non.econ", "growth_collection_non.econ",
                                            "energy_non.econ", "water_non.econ",
                                            "preparation_econ",
                                            "export_per_capita", "import_per_capita", "domestic_per_capita")))

# Plot the results (all countries)

# Factor inputs (labor + energy)
# Have three main colors for footprint_type; and two alpha values for type (hr_m and hr_f)
# "export_per_capita" is placed on top of "domestic_per_capita" for each type.
p_hr_f = plot_countries(summary_food_df_long_with_ghd %>% filter(type %in% c("hr_f")),
               "Daily time footprint per capita (hr/cap/day)",
               paste0("Female time footprint per capita by country (", year, ")"))

p_hr_m = plot_countries(summary_food_df_long_with_ghd %>% filter(type %in% c("hr_m")),
               "Daily time footprint per capita (hr/cap/day)",
               paste0("Male time footprint per capita by country (", year, ")"))
               
p_en = plot_countries(summary_food_df_long %>% filter(type %in% c("en")),
               "Food sector energy footprint per capita (GJ/cap/yr)",
               paste0("Energy footprint per capita by country (", year, ")"))

# plot_countries() hardcodes aes(x=country, ...); summary_nonfood_df_long is
# region-resolved (exio_region), so rename for these standalone nonfood-only
# displays -- each bar is really one of the 49 EXIO regions (a RoW aggregate
# for 5 of them), not a FABIO country, but plot_countries() itself is
# unchanged/shared with the food-sector plots above.
p_hr_f_nonfood = plot_countries(summary_nonfood_df_long %>% filter(type %in% c("hr_f")) %>% rename(country = exio_region),
               "Daily time footprint per capita (hr/cap/day)",
               paste0("Female nonfood time footprint per capita by EXIO region (", year, ")"))

p_hr_m_nonfood = plot_countries(summary_nonfood_df_long %>% filter(type %in% c("hr_m")) %>% rename(country = exio_region),
               "Daily time footprint per capita (hr/cap/day)",
               paste0("Male nonfood time footprint per capita by EXIO region (", year, ")"))

p_en_nonfood = plot_countries(summary_nonfood_df_long %>% filter(type %in% c("en")) %>% rename(country = exio_region),
               "Nonfood sector energy footprint per capita (GJ/cap/yr)",
               paste0("Energy footprint per capita by EXIO region (", year, ")"))

# Stack data from df_ghd_gender to the labor footprint plots by gender
# df_ghd_gender has two columns "maleTotalHours" and "femaleTotalHours" which are the total hours spent on food-related activities per capita per day by activity.
# Add these hours to the existing bars in p_hr_f and p_hr_m, with different alpha values for different activities (growth_collection, processing, preparation).
p_hr_f = p_hr_f + ylim(-1, 3.5) 
p_hr_m = p_hr_m + ylim(-1, 3.5) 


# Combine three panels vertically with patchwork and have only one shared x-axis 
library(patchwork)
p_combined  = p_hr_f / p_hr_m + plot_layout(guides = "collect") & theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# keep x labels only on bottom panel
p_combined[[2]] <- p_combined[[2]] +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.ticks.x = element_line())
print(p_combined)


# Plot countries with non-econ obs

# See only the countries with non.econ observations
partial_cty = summary_food_df_long_with_ghd %>% 
  filter(footprint_type == "preparation_non.econ") %>% pull(country) %>% unique()
partial_ord = (summary_food_df_long_with_ghd %>% filter(country %in% partial_cty) %>%
                 filter(type %in% c("hr_f"), footprint_type != "import_per_capita") %>%
                 group_by(country) %>% 
                 summarise(d = sum(per_capita_value, na.rm=TRUE)) %>% 
                 arrange(-d))$country
partial_df = summary_food_df_long_with_ghd %>% filter(country %in% partial_cty) %>%
  mutate(country = factor(country, levels = partial_ord)) %>%
  arrange(country)

fao_pro_raw = read.csv("data/FAOSTAT_data_en_5-28-2026 protein.csv", check.names = FALSE)
fao_pro_lookup = fao_pro_raw %>%
  filter(`Indicator Code` == "4004") %>%
  group_by(`Area Code (M49)`, Area) %>%
  summarise(pro_per_cap_day = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  mutate(iso3c = countrycode::countrycode(`Area Code (M49)`,
                                          origin = "un", destination = "iso3c",
                                          warn = FALSE)) %>%
  filter(!is.na(iso3c)) %>%
  select(country = iso3c, pro_per_cap_day)

pro_partial = country_summary(agg_country_footprint(FABIO_y_hh_pro)) %>%
  # mutate(pro_per_cap_day = (domestic_per_capita + import_per_capita) / 1e6 / 365) %>%
  left_join(fao_pro_lookup, by = "country") %>%
  filter(as.character(country) %in% as.character(partial_cty)) %>%
  mutate(country = factor(as.character(country), levels = as.character(partial_ord))) %>%
  select(country, pro_per_cap_day)
pro_max   <- max(pro_partial$pro_per_cap_day, na.rm = TRUE)
pro_scale <- 3.0 / pro_max

fao_und_raw = read.csv("data/FAOSTAT_undernourishment_3yravg_2016-2022.csv", check.names = FALSE)
fao_und_lookup = fao_und_raw %>%
  filter(Year == "2018-2020") %>%
  # FAO censors low prevalence as "<2.5" instead of a number; treat that as 2.5.
  mutate(Value = as.numeric(gsub("<", "", Value)),
         `Area Code (M49)` = as.numeric(gsub("'", "", `Area Code (M49)`)),
         iso3c = countrycode::countrycode(`Area Code (M49)`,
                                          origin = "un", destination = "iso3c",
                                          warn = FALSE)) %>%
  filter(!is.na(iso3c)) %>%
  select(country = iso3c, undernourishment_pct = Value)

und_partial = tibble(country = as.character(partial_cty)) %>%
  left_join(fao_und_lookup, by = "country") %>%
  mutate(country = factor(country, levels = as.character(partial_ord))) %>%
  select(country, undernourishment_pct)
und_max   <- max(und_partial$undernourishment_pct, na.rm = TRUE)
und_scale <- 3.0 / und_max

p1 = plot_countries(partial_df %>% filter(type %in% c("hr_f"), footprint_type != "import_per_capita"), "Resident time use per capita: Female (hr/day)", "") +
  geom_point(data = pro_partial, aes(x = country, y = pro_per_cap_day * pro_scale),
             color = "black", size = 2.5, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(-0.2, 3.5),
                     sec.axis = sec_axis(~ . / pro_scale, name = "g protein/cap/day"))
p2 = plot_countries(partial_df %>% filter(type %in% c("hr_m"), footprint_type != "import_per_capita"), "Resident time use per capita: Male (hr/day)", "") +
  geom_point(data = und_partial, aes(x = country, y = undernourishment_pct * und_scale),
             color = "black", size = 2.5, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(-0.2, 3.5),
                     sec.axis = sec_axis(~ . / und_scale, name = "Prevalence of undernourishment (%)"))
p3 = plot_countries(partial_df %>% filter(type %in% c("en")), "Energy footprint per capita (GJ/yr)", "")

# Protein supply vs. total food provisioning time, one point per country per gender
# (excludes export_per_capita, since that time provisions other countries' consumption)
time_partial = partial_df %>%
  filter(type %in% c("hr_m", "hr_f"),
         !footprint_type %in% c("export_per_capita", "import_per_capita")) %>%
  group_by(country, type, is_row) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

df_time_protein = time_partial %>%
  left_join(pro_partial, by = "country") %>%
  drop_na(pro_per_cap_day)

label_partial = df_time_protein %>%
  group_by(country) %>%
  slice_max(hr_per_cap_day, n = 1, with_ties = FALSE) %>%
  mutate(fontface = ifelse(is_row, "plain", "bold"))

p4 = ggplot(df_time_protein, aes(x = pro_per_cap_day, y = hr_per_cap_day)) +
  geom_line(aes(group = country, linetype = is_row, linewidth = is_row), color = "grey60") +
  geom_point(aes(color = type), size = 2.5) +
  ggrepel::geom_text_repel(
    data = label_partial,
    aes(label = country, fontface = I(fontface)), size = 3, max.overlaps = 20, show.legend = FALSE) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
  scale_color_manual(values = c(hr_f = "#ca2323", hr_m = "#1f77b4"),
                      labels = c(hr_f = "Female", hr_m = "Male")) +
  scale_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dotted"),
                         labels = c("FALSE" = "Directly modeled (EXIO)", "TRUE" = "Aggregated 'Rest of World' region"),
                         name = "Data coverage") +
  scale_linewidth_manual(values = c("FALSE" = 1.1, "TRUE" = 0.4), guide = "none") +
  # scale_x_continuous(limits = c(0, NA)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Protein supply (g/cap/day)", y = "Food provisioning time (hr/cap/day)",
       color = "Gender") +
  theme_minimal()
dev.new(width = 16, height = 6)
print(p4)

# Same protein-vs-time scatter, restricted to non-RoW countries only
df_time_protein_nonrow = df_time_protein %>% filter(!is_row)

label_partial_nonrow = df_time_protein_nonrow %>%
  group_by(country) %>%
  slice_max(hr_per_cap_day, n = 1, with_ties = FALSE)

p4_nonrow = ggplot(df_time_protein_nonrow, aes(x = pro_per_cap_day, y = hr_per_cap_day)) +
  geom_line(aes(group = country), color = "grey60", linewidth = 1.1) +
  geom_point(aes(color = type), size = 2.5) +
  ggrepel::geom_text_repel(
    data = label_partial_nonrow,
    aes(label = country), size = 3, fontface = "bold", max.overlaps = 20, show.legend = FALSE) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
  scale_color_manual(values = c(hr_f = "#ca2323", hr_m = "#1f77b4"),
                      labels = c(hr_f = "Female", hr_m = "Male")) +
  # scale_x_continuous(limits = c(0, NA)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Protein supply (g/cap/day)", y = "Food provisioning time (hr/cap/day)",
       color = "Gender") +
  theme_minimal()
dev.new(width = 16, height = 6)
print(p4_nonrow)

# Same as p4_nonrow, but with the y-axis as time per 50 g protein (hr/50g)
# instead of total time per capita per day -- i.e. each country's own time is
# divided by its own protein supply, so the metric is a rate rather than a total.
df_time_protein_nonrow_50g = df_time_protein_nonrow %>%
  mutate(hr_per_50g_protein = hr_per_cap_day / pro_per_cap_day * 50)

label_partial_nonrow_50g = df_time_protein_nonrow_50g %>%
  group_by(country) %>%
  slice_max(hr_per_50g_protein, n = 1, with_ties = FALSE)

p4_nonrow_50g = ggplot(df_time_protein_nonrow_50g, aes(x = pro_per_cap_day, y = hr_per_50g_protein)) +
  geom_line(aes(group = country), color = "grey60", linewidth = 1.1) +
  geom_point(aes(color = type), size = 2.5) +
  ggrepel::geom_text_repel(
    data = label_partial_nonrow_50g,
    aes(label = country), size = 3, fontface = "bold", max.overlaps = 20, show.legend = FALSE) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
  scale_color_manual(values = c(hr_f = "#ca2323", hr_m = "#1f77b4"),
                      labels = c(hr_f = "Female", hr_m = "Male")) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Protein supply (g/cap/day)", y = "Food provisioning time per 50 g protein (hr/50g)",
       color = "Gender") +
  theme_minimal()
dev.new(width = 16, height = 6)
print(p4_nonrow_50g)

# Same as p4_nonrow_50g, but adds indirect (non-food-sector) domestic labor --
# transport, packaging, etc. embodied in domestically-produced food -- on top
# of the direct food-sector + household time already in hr_per_cap_day. Non-
# food is EXIO-region-resolved (summary_nonfood_df), so for these non-RoW
# countries each has its own 1:1 region (iso_to_region); domestic_per_capita
# there is already hr/cap/day, same units/convention as the food side (see
# country_summary()/region_summary()). Both the food-only and food+indirect
# values are kept and plotted (open vs. filled point, joined by a thin colored
# segment) so the added indirect component stays visible per point, not just
# folded silently into a bigger total.
nonfood_domestic_nonrow = summary_nonfood_df %>%
  filter(type %in% c("hr_m", "hr_f"), !is_row) %>%
  select(exio_region, type, hr_per_cap_day_nonfood = domestic_per_capita)

df_time_protein_nonrow_indirect = df_time_protein_nonrow_50g %>%
  mutate(exio_region = iso_to_region[as.character(country)]) %>%
  left_join(nonfood_domestic_nonrow, by = c("exio_region", "type")) %>%
  mutate(hr_per_cap_day_nonfood    = replace_na(hr_per_cap_day_nonfood, 0),
         hr_per_50g_protein_direct = hr_per_50g_protein,
         hr_per_50g_protein_total  = (hr_per_cap_day + hr_per_cap_day_nonfood) / pro_per_cap_day * 50)

label_partial_nonrow_indirect = df_time_protein_nonrow_indirect %>%
  group_by(country) %>%
  slice_max(hr_per_50g_protein_total, n = 1, with_ties = FALSE)

df_components_nonrow_indirect = df_time_protein_nonrow_indirect %>%
  select(country, type, pro_per_cap_day,
         direct = hr_per_50g_protein_direct, total = hr_per_50g_protein_total) %>%
  pivot_longer(c(direct, total), names_to = "component", values_to = "hr_per_50g_protein") %>%
  mutate(component = factor(component, levels = c("direct", "total"),
                            labels = c("Food (direct)", "+ Non-food (indirect)")))

p4_nonrow_50g_indirect = ggplot(df_components_nonrow_indirect,
                                aes(x = pro_per_cap_day, y = hr_per_50g_protein)) +
  geom_line(data = df_time_protein_nonrow_indirect,
            aes(x = pro_per_cap_day, y = hr_per_50g_protein_total, group = country),
            color = "grey60", linewidth = 1.1, inherit.aes = FALSE) +
  geom_segment(data = df_time_protein_nonrow_indirect,
              aes(x = pro_per_cap_day, xend = pro_per_cap_day,
                  y = hr_per_50g_protein_direct, yend = hr_per_50g_protein_total, color = type),
              linewidth = 0.7, alpha = 0.5, inherit.aes = FALSE) +
  geom_point(aes(color = type, shape = component), size = 2.5) +
  ggrepel::geom_text_repel(
    data = label_partial_nonrow_indirect,
    aes(x = pro_per_cap_day, y = hr_per_50g_protein_total, label = country),
    size = 3, fontface = "bold", max.overlaps = 20, show.legend = FALSE, inherit.aes = FALSE) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
  scale_color_manual(values = c(hr_f = "#ca2323", hr_m = "#1f77b4"),
                      labels = c(hr_f = "Female", hr_m = "Male")) +
  scale_shape_manual(values = c("Food (direct)" = 1, "+ Non-food (indirect)" = 16)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Protein supply (g/cap/day)",
       y = "Domestic food provisioning time per 50 g protein (hr/50g)",
       color = "Gender", shape = "Effort component") +
  theme_minimal()
dev.new(width = 16, height = 6)
print(p4_nonrow_50g_indirect)

# Table version of the same food/non-food split, by non-RoW country and gender
# -- this time also separating food-sector paid (economic: domestic_per_capita
# + preparation_econ, the FABIO/EXIO labor satellite + GHD's paid food-service
# time) from unpaid (household non-econ: preparation/processing/growth_
# collection, GHD's time-use time), and adding total daily protein supply as
# its own column (the un-normalized denominator, alongside the hr/50g rates).
food_time_split_nonrow = partial_df %>%
  filter(!is_row, type %in% c("hr_m", "hr_f"),
         !footprint_type %in% c("export_per_capita", "import_per_capita")) %>%
  mutate(pay_type = ifelse(grepl("non\\.econ", footprint_type), "unpaid_food_hr", "paid_food_hr"),
         country = as.character(country)) %>%
  group_by(country, type, pay_type) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = pay_type, values_from = hr_per_cap_day, values_fill = 0)

table_food_nonfood_50g = food_time_split_nonrow %>%
  left_join(pro_partial %>% mutate(country = as.character(country)), by = "country") %>%
  mutate(exio_region = iso_to_region[country]) %>%
  left_join(nonfood_domestic_nonrow, by = c("exio_region", "type")) %>%
  mutate(hr_per_cap_day_nonfood = replace_na(hr_per_cap_day_nonfood, 0),
         gender = ifelse(type == "hr_f", "Female", "Male"),
         unpaid_food_hr_per_50g = unpaid_food_hr / pro_per_cap_day * 50,
         paid_food_hr_per_50g   = paid_food_hr   / pro_per_cap_day * 50,
         nonfood_hr_per_50g     = hr_per_cap_day_nonfood / pro_per_cap_day * 50,
         nonfood_share = nonfood_hr_per_50g / (paid_food_hr_per_50g + nonfood_hr_per_50g) * 100,
         total_hr_per_50g       = unpaid_food_hr_per_50g + paid_food_hr_per_50g + nonfood_hr_per_50g) %>%
  select(country, gender, protein_g_per_cap_day = pro_per_cap_day,
         unpaid_food_hr_per_50g, paid_food_hr_per_50g, nonfood_hr_per_50g, nonfood_share, total_hr_per_50g) %>%
  arrange(desc(total_hr_per_50g))

gt_food_nonfood_50g = table_food_nonfood_50g %>%
  gt() %>%
  fmt_number(columns = protein_g_per_cap_day, decimals = 1) %>%
  fmt_number(columns = c(unpaid_food_hr_per_50g, paid_food_hr_per_50g, nonfood_hr_per_50g, total_hr_per_50g), decimals = 3) %>%
  cols_label(country = "Country", gender = "Gender",
             protein_g_per_cap_day = "Protein supply (g/cap/day)",
             unpaid_food_hr_per_50g = "Food, unpaid (household), hr/50g",
             paid_food_hr_per_50g = "Food, paid (economic), hr/50g",
             nonfood_hr_per_50g = "Non-food, paid (indirect), hr/50g",
             nonfood_share = "Non-food share of total paid, %",
             total_hr_per_50g = "Total, hr/50g") %>%
  tab_header(title = "Domestic effort per 50 g protein: unpaid food, paid food, and non-food sector",
             subtitle = paste0("Non-RoW countries, ", year))
print(gt_food_nonfood_50g)

# Domestic-only variant: domestic + household time (same as above) vs. FABIO's
# domestic-only protein supply (instead of FAO's domestic+import national total)
pro_domestic_partial = country_summary(agg_country_footprint(FABIO_y_hh_pro)) %>%
  mutate(across(ends_with("per_capita"), ~ .x / 1e6 / 365)) %>%
  filter(as.character(country) %in% as.character(partial_cty)) %>%
  mutate(country = factor(as.character(country), levels = as.character(partial_ord))) %>%
  select(country, pro_per_cap_day_domestic = domestic_per_capita)

df_time_protein_domestic = time_partial %>%
  left_join(pro_domestic_partial, by = "country") %>%
  drop_na(pro_per_cap_day_domestic)

label_partial_domestic = df_time_protein_domestic %>%
  group_by(country) %>%
  slice_max(hr_per_cap_day, n = 1, with_ties = FALSE) %>%
  mutate(fontface = ifelse(is_row, "plain", "bold"))

p4_domestic = ggplot(df_time_protein_domestic, aes(x = pro_per_cap_day_domestic, y = hr_per_cap_day)) +
  geom_line(aes(group = country, linetype = is_row, linewidth = is_row), color = "grey60") +
  geom_point(aes(color = type), size = 2.5) +
  ggrepel::geom_text_repel(
    data = label_partial_domestic,
    aes(label = country, fontface = I(fontface)), size = 3, max.overlaps = 20, show.legend = FALSE) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
  scale_color_manual(values = c(hr_f = "#ca2323", hr_m = "#1f77b4"),
                      labels = c(hr_f = "Female", hr_m = "Male")) +
  scale_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dotted"),
                         labels = c("FALSE" = "Directly modeled (EXIO)", "TRUE" = "Aggregated 'Rest of World' region"),
                         name = "Data coverage") +
  scale_linewidth_manual(values = c("FALSE" = 1.1, "TRUE" = 0.4), guide = "none") +
  scale_x_continuous(limits = c(0, NA)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Domestic protein supply (g/cap/day, FABIO)",
       y = "Domestic + household food provisioning time (hr/cap/day)",
       color = "Gender") +
  theme_minimal()
dev.new(width = 16, height = 6)
print(p4_domestic)

p_combined_partial  = p1 / p2 + plot_layout(guides = "collect") & theme(legend.position = "top",
                                                                        axis.text.x = element_blank(),
                                                                        axis.ticks.x = element_blank(),
                                                                        text              = element_text(size = 15),
                                                                        axis.text.y       = element_text(size = 22),
                                                                        axis.title        = element_text(size = 15),
                                                                        legend.text       = element_text(size = 22),
                                                                        legend.title      = element_text(size = 15)
)
p_combined_partial[[2]] <- p_combined_partial[[2]] +
  theme(axis.text.x = element_text(angle=90, hjust=1, size = 22),
        axis.ticks.x = element_line())
print(p_combined_partial)

# Save it to a PDF
ggsave("results/footprint_all_countries.pdf", p_combined, width = 18, height = 12)
ggsave("results/footprint_partial_countries.pdf", p_combined_partial, width = 18, height = 12)
ggsave("results/protein_vs_time.pdf", p4, width = 16, height = 6)
ggsave("results/protein_vs_time_nonrow.pdf", p4_nonrow, width = 16, height = 6)
ggsave("results/protein_vs_time_per50g_nonrow.pdf", p4_nonrow_50g, width = 16, height = 6)
ggsave("results/protein_vs_time_per50g_nonrow_indirect.pdf", p4_nonrow_50g_indirect, width = 12, height = 6)
ggsave("results/protein_vs_time_domestic.pdf", p4_domestic, width = 16, height = 6)



# Same dual-axis figure, restricted to non-RoW countries only
nonrow_ord = partial_ord[!as.character(partial_ord) %in% row_countries]

partial_df_nonrow = partial_df %>% filter(!is_row) %>%
  mutate(country = factor(as.character(country), levels = as.character(nonrow_ord)))

pro_partial_nonrow = pro_partial %>%
  filter(as.character(country) %in% as.character(nonrow_ord)) %>%
  mutate(country = factor(as.character(country), levels = as.character(nonrow_ord)))
pro_max_nonrow   <- max(pro_partial_nonrow$pro_per_cap_day, na.rm = TRUE)
pro_scale_nonrow <- 3.0 / pro_max_nonrow

und_partial_nonrow = und_partial %>%
  filter(as.character(country) %in% as.character(nonrow_ord)) %>%
  mutate(country = factor(as.character(country), levels = as.character(nonrow_ord)))
und_max_nonrow   <- max(und_partial_nonrow$undernourishment_pct, na.rm = TRUE)
und_scale_nonrow <- 3.0 / und_max_nonrow

p1_nonrow = plot_countries(partial_df_nonrow %>% filter(type %in% c("hr_f"), footprint_type != "import_per_capita"), "Resident time use per capita: Female (hr/day)", "") +
  geom_point(data = pro_partial_nonrow, aes(x = country, y = pro_per_cap_day * pro_scale_nonrow),
             color = "black", size = 2.5, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(-0.2, 3.5),
                     sec.axis = sec_axis(~ . / pro_scale_nonrow, name = "g protein/cap/day"))
p2_nonrow = plot_countries(partial_df_nonrow %>% filter(type %in% c("hr_m"), footprint_type != "import_per_capita"), "Resident time use per capita: Male (hr/day)", "") +
  geom_point(data = und_partial_nonrow, aes(x = country, y = undernourishment_pct * und_scale_nonrow),
             color = "black", size = 2.5, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(-0.2, 3.5),
                     sec.axis = sec_axis(~ . / und_scale_nonrow, name = "Prevalence of undernourishment (%)"))

p_combined_partial_nonrow  = p1_nonrow / p2_nonrow + plot_layout(guides = "collect") & theme(legend.position = "top",
                                                                        axis.text.x = element_blank(),
                                                                        axis.ticks.x = element_blank(),
                                                                        text              = element_text(size = 15),
                                                                        axis.text.y       = element_text(size = 22),
                                                                        axis.title        = element_text(size = 15),
                                                                        legend.text       = element_text(size = 22),
                                                                        legend.title      = element_text(size = 15)
)
p_combined_partial_nonrow[[2]] <- p_combined_partial_nonrow[[2]] +
  theme(axis.text.x = element_text(angle=90, hjust=1, size = 22),
        axis.ticks.x = element_line())
print(p_combined_partial_nonrow)

ggsave("results/footprint_partial_countries_nonrow update.pdf", p_combined_partial_nonrow, width = 18, height = 12)

# Non-food sector: combined female + male + energy, ordered by non-food female time.
# Region-resolved (see plot_countries() note above) -- ordered/renamed at point of use.
nonfood_ord = (summary_nonfood_df_long %>%
  filter(type == "hr_m", footprint_type == "domestic_per_capita") %>%
  arrange(-per_capita_value))$exio_region

p_hr_f_nf = plot_countries(
  summary_nonfood_df_long %>% filter(type == "hr_f") %>%
    mutate(exio_region = factor(exio_region, levels = nonfood_ord)) %>% rename(country = exio_region),
  "Daily time footprint per capita (hr/cap/day)",
  paste0("Female non-food time footprint per capita by EXIO region (", year, ")"))

p_hr_m_nf = plot_countries(
  summary_nonfood_df_long %>% filter(type == "hr_m") %>%
    mutate(exio_region = factor(exio_region, levels = nonfood_ord)) %>% rename(country = exio_region),
  "Daily time footprint per capita (hr/cap/day)",
  paste0("Male non-food time footprint per capita by EXIO region (", year, ")"))

p_en_nf = plot_countries(
  summary_nonfood_df_long %>% filter(type == "en") %>%
    mutate(exio_region = factor(exio_region, levels = nonfood_ord)) %>% rename(country = exio_region),
  "Non-food sector energy footprint per capita (GJ/cap/yr)",
  paste0("Non-food energy footprint per capita by EXIO region (", year, ")"))

p_combined_nonfood = p_hr_f_nf / p_hr_m_nf / p_en_nf +
  plot_layout(guides = "collect") & theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
p_combined_nonfood[[3]] <- p_combined_nonfood[[3]] +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.x = element_line())

ggsave("results/footprint_nonfood_countries.pdf", p_combined_nonfood, width = 18, height = 18)

# Plot directly-modelled EXIO countries only (non-RoW)
direct_ord = (summary_food_df_long_with_ghd %>%
  filter(!is_row, type == "hr_f") %>%
  group_by(country) %>%
  summarise(d = sum(per_capita_value, na.rm = TRUE)) %>%
  arrange(-d))$country

b = summary_food_df_long_with_ghd %>%
  filter(!is_row, footprint_type != "import_per_capita") %>%
  mutate(country = factor(country, levels = direct_ord))

# Add non-food-sector domestic/export time (labor in packaging, transport, etc.
# that supports food provisioning) into the same bars. footprint_type gets a
# "_nf" suffix so plot_countries() can color and stack it separately from the
# matching food component. (Non-food import is added below, disaggregated by
# continent, same as food import.)
#
# summary_nonfood_df_long is region-resolved (exio_region); restricted to
# !is_row here, EXIO region <-> FABIO country is 1:1, so translate exio_region
# back to its matching ISO3 code (region_to_iso_1to1, defined near the top
# alongside iso_to_region) -- not just rename the column -- to align with
# direct_ord/b's country-keyed x-axis. A plain rename would leave `country`
# holding region NAMES ("China") that never match direct_ord's ISO3 codes
# ("CHN"), silently sending every bar's x position to NA instead of erroring.
b_nonfood = summary_nonfood_df_long %>%
  filter(!is_row, type %in% c("hr_f", "hr_m"), footprint_type != "import_per_capita") %>%
  mutate(country = factor(region_to_iso_1to1[as.character(exio_region)], levels = direct_ord),
         footprint_type = paste0(as.character(footprint_type), "_nf"))

# Continent-of-origin breakdown of the import component -----------------------------
# For each consuming country, split the import total (colSums(mat) - diag(mat)) by the
# continent of the producing country -- the same continent grouping (regions$continent)
# and imp_by_cont logic used by make_continent_data()'s single-country spotlight plots,
# just computed for every country at once instead of one focal iso.
import_by_continent = function(mat) {
  cty  = colnames(mat)
  cont = regions$continent[match(cty, regions$iso3c)]

  bind_rows(lapply(seq_along(cty), function(j) {
    not_focal = seq_along(cty) != j
    imp_by_cont = tapply(mat[not_focal, j], cont[not_focal], sum, na.rm = TRUE)
    data.frame(country = cty[j], continent = names(imp_by_cont), import = as.numeric(imp_by_cont))
  }))
}

# Non-food counterpart: mat is 49 (EXIO region) x 187 (FABIO country) --
# origin is region-resolved, so exclude the focal country's own REGION (not
# just its own column) and group by region_continent_of_index rather than
# regions$continent.
import_by_continent_nonfood = function(mat) {
  cty = colnames(mat)
  bind_rows(lapply(seq_along(cty), function(j) {
    own_region_idx = region_row_of_country[j]
    not_focal_region = seq_len(n_reg_EXIO) != own_region_idx
    imp_by_cont = tapply(mat[not_focal_region, j], region_continent_of_index[not_focal_region], sum, na.rm = TRUE)
    data.frame(country = cty[j], continent = names(imp_by_cont), import = as.numeric(imp_by_cont))
  }))
}

pop_data_yr = subset(countrypops, year == yr) %>% select(iso3c = country_code_3, population)

# Long df of import_per_capita by continent-of-origin for hr_m/hr_f, in the same units as
# summary_food_df_long/summary_nonfood_df_long's import_per_capita (M.hr -> hr/cap/day).
# footprint_type = "import_cont_<Continent>" (food) / "..._nf" (non-food); plot_countries()
# recognizes both the "import" and "_nf" naming to color/stack/negate them correctly.
make_import_continent_df = function(l_country, sector_suffix) {
  bind_rows(lapply(c("hr_m", "hr_f"), function(m) {
    mat = l_country[[m]]
    imp_df = if (nrow(mat) == nrreg) import_by_continent(mat) else import_by_continent_nonfood(mat)
    imp_df %>%
      filter(country %in% direct_ord) %>%
      left_join(pop_data_yr, by = c("country" = "iso3c")) %>%
      mutate(country = factor(country, levels = direct_ord),
             type = m,
             is_row = FALSE,
             per_capita_value = import / population * 1e6 / 365,
             footprint_type = paste0("import_cont_", continent, sector_suffix)) %>%
      select(country, type, is_row, footprint_type, per_capita_value)
  }))
}

b_import = bind_rows(make_import_continent_df(l_food_country, ""),
                      make_import_continent_df(l_nonfood_region, "_nf"))
import_levels = sort(unique(b_import$footprint_type))

b_direct = bind_rows(b %>% mutate(footprint_type = as.character(footprint_type)), b_nonfood, b_import) %>%
  mutate(footprint_type = factor(footprint_type, levels = c(
    "preparation_non.econ", "processing_non.econ", "growth_collection_non.econ",
    "energy_non.econ", "water_non.econ",
    "preparation_econ",
    "domestic_per_capita_nf", "export_per_capita_nf",
    "domestic_per_capita", "export_per_capita",
    import_levels
  )))

# Share y-axis limits between the female/male panels, sized to the combined
# food + non-food stack (positive) and food + non-food import (negative)
stack_totals_direct = b_direct %>%
  filter(type %in% c("hr_m", "hr_f")) %>%
  mutate(part = ifelse(grepl("^import_cont_", footprint_type), "neg", "pos")) %>%
  group_by(country, type, part) %>%
  summarise(total = sum(per_capita_value, na.rm = TRUE), .groups = "drop")
y_max_direct =  max(stack_totals_direct$total[stack_totals_direct$part == "pos"], na.rm = TRUE) * 1.1
y_min_direct = -max(stack_totals_direct$total[stack_totals_direct$part == "neg"], na.rm = TRUE) * 1.1

p_hr_f_direct = plot_countries(b_direct %>% filter(type == "hr_f"), "Female time footprint per capita (hr/day)", "") + ylim(y_min_direct, y_max_direct)
p_hr_m_direct_base = plot_countries(b_direct %>% filter(type == "hr_m"), "Male time footprint per capita (hr/day)", "") + ylim(y_min_direct, y_max_direct)

# LUX (male) has by far the largest import-effort bar (drives y_min_direct itself --
# see the y-axis scaling fix above) and its continent-of-origin makeup isn't obvious
# from color alone, so label each import_cont_* segment right next to its own block,
# in that block's own fill color (via plot_countries()'s attached "fill_scheme").
# Segment y-ranges are derived analytically rather than via ggplot_build(): with
# position_stack() on negative values, the FIRST level of footprint_type ends up
# farthest from zero and the LAST level ends up closest to zero (verified
# empirically) -- ordering by import_levels and cumulative-summing from the total
# reproduces that exactly.
lux_idx = which(direct_ord == "LUX")
fill_scheme = attr(p_hr_m_direct_base, "fill_scheme")

lux_segments = b_direct %>%
  filter(country == "LUX", type == "hr_m", grepl("^import_cont_", as.character(footprint_type))) %>%
  mutate(footprint_type = factor(as.character(footprint_type), levels = import_levels)) %>%
  arrange(footprint_type) %>%
  mutate(ymax = cumsum(per_capita_value) - sum(per_capita_value),
         ymin = ymax - per_capita_value,
         y_mid = (ymin + ymax) / 2,
         continent = gsub("_nf$", "", gsub("^import_cont_", "", as.character(footprint_type))),
         sector = ifelse(grepl("_nf$", as.character(footprint_type)), "non-food", "food"),
         label = paste0(continent, ": ", round(per_capita_value, 2), "h"),
         color = fill_scheme[as.character(footprint_type)]) %>%
  filter(per_capita_value >= 0.15)  # skip slivers too thin for non-overlapping labels

p_hr_m_direct = p_hr_m_direct_base +
  { if (length(lux_idx) == 1 && nrow(lux_segments) > 0)
      geom_text(data = lux_segments, aes(x = lux_idx + 0.6, y = y_mid, label = label),
                color = lux_segments$color, hjust = 0, size = 3, fontface = "bold", inherit.aes = FALSE) }

p_combined_direct = p_hr_f_direct / p_hr_m_direct + plot_layout(guides = "collect") & theme(
  legend.position = "top",
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()
)
p_combined_direct[[2]] <- p_combined_direct[[2]] +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.x = element_line())
print(p_combined_direct)

ggsave(paste0("results/footprint_direct_countries_incl_nonfood update.pdf"), p_combined_direct, width = 18, height = 12)





# Nutrient plotting ####
df_nutri = list("kcal" = country_summary(agg_country_footprint(FABIO_y_hh_cal)),
                "protein" = country_summary(agg_country_footprint(FABIO_y_hh_pro))
)

for (i in 1:length(df_nutri)) {
  df_nutri[[i]] %>% mutate(across(ends_with("per_capita"), ~ .x / 1e6 / 365)) -> df_nutri[[i]]
}

summary_kcal_df_long = df_nutri[["kcal"]] %>% 
  select(-c(population, domestic, export, import)) %>%
  pivot_longer(cols = c("domestic_per_capita", "export_per_capita", "import_per_capita"), 
               names_to = "footprint_type", values_to = "per_capita_value") %>%
  # Order countries by sum of domestic_per_capita of type "hr_m" and "hr_f"
  mutate(country = factor(country, levels = sum_ord),
         footprint_type = factor(footprint_type, 
                                 levels = c("export_per_capita", "import_per_capita", "domestic_per_capita"))) %>%
  drop_na() %>%
  mutate(cat = case_when(
    footprint_type == "export_per_capita" ~ "export",
    footprint_type == "import_per_capita" ~ "import",
    .default = "domestic"
  ))
summary_pro_df_long = df_nutri[["protein"]] %>%
  select(-c(population, domestic, export, import)) %>%
  pivot_longer(cols = c("domestic_per_capita", "export_per_capita", "import_per_capita"), 
               names_to = "footprint_type", values_to = "per_capita_value") %>%
  # Order countries by sum of domestic_per_capita of type "hr_m" and "hr_f"
  mutate(country = factor(country, levels = sum_ord),
         footprint_type = factor(footprint_type, 
                                 levels = c("export_per_capita", "import_per_capita", "domestic_per_capita"))) %>%
  drop_na() %>%
  mutate(cat = case_when(
    footprint_type == "export_per_capita" ~ "export",
    footprint_type == "import_per_capita" ~ "import",
    .default = "domestic"
  ))

# Region-resolved (49-row) counterparts, built the same way as
# summary_food_region_df_long -- protein/kcal data has no RoW ambiguity of its
# own (it's purely FABIO-native, agg_country_footprint(FABIO_y_hh_pro/cal),
# never touches EXIOBASE's region-pasted satellites), but it needs to be
# summed to the same 49-region resolution as summary_nonfood_df_long/
# summary_food_region_df_long wherever it's paired with combined food+nonfood
# effort (e.g. the protein-per-hour conversion-factor charts). summary_kcal_
# df_long/summary_pro_df_long (187-country) remain untouched for food-only
# nutrient displays.
df_nutri_region = list(
  kcal    = region_summary(agg_country_footprint(FABIO_y_hh_cal)),
  protein = region_summary(agg_country_footprint(FABIO_y_hh_pro))
)
for (i in 1:length(df_nutri_region)) {
  df_nutri_region[[i]] %>% mutate(across(ends_with("per_capita"), ~ .x / 1e6 / 365)) -> df_nutri_region[[i]]
}

summary_kcal_region_df_long = df_nutri_region[["kcal"]] %>%
  select(-c(population, domestic, export, import)) %>%
  pivot_longer(cols = c("domestic_per_capita", "export_per_capita", "import_per_capita"),
               names_to = "footprint_type", values_to = "per_capita_value") %>%
  mutate(exio_region = factor(exio_region, levels = sum_ord_region),
         footprint_type = factor(footprint_type,
                                 levels = c("export_per_capita", "import_per_capita", "domestic_per_capita"))) %>%
  drop_na() %>%
  mutate(cat = case_when(
    footprint_type == "export_per_capita" ~ "export",
    footprint_type == "import_per_capita" ~ "import",
    .default = "domestic"
  ))
summary_pro_region_df_long = df_nutri_region[["protein"]] %>%
  select(-c(population, domestic, export, import)) %>%
  pivot_longer(cols = c("domestic_per_capita", "export_per_capita", "import_per_capita"),
               names_to = "footprint_type", values_to = "per_capita_value") %>%
  mutate(exio_region = factor(exio_region, levels = sum_ord_region),
         footprint_type = factor(footprint_type,
                                 levels = c("export_per_capita", "import_per_capita", "domestic_per_capita"))) %>%
  drop_na() %>%
  mutate(cat = case_when(
    footprint_type == "export_per_capita" ~ "export",
    footprint_type == "import_per_capita" ~ "import",
    .default = "domestic"
  ))

p_kcal = plot_countries(summary_kcal_df_long, "Daily kcal supply per capita (kcal/cap/day)", "kcal")
p_protein = plot_countries(summary_pro_df_long, "Daily protein supply per capita (g/cap/day)", "protein")
ggsave("results/kcal_supply.pdf",    p_kcal,    width = 18, height = 6)
ggsave("results/protein_supply.pdf", p_protein, width = 18, height = 6)


# Domestic vs. import conversion factors: protein per hour and protein per unit energy,
# using total (food-sector + non-food-sector) economic hours/energy — no household time.
#
# Combines food+nonfood effort, so operates at EXIO-region resolution (49
# regions) throughout, per chat discussion -- non-food (and, to match it,
# protein/kcal here) can't honestly be attributed below region granularity for
# RoW aggregates. `countries` (a list of FABIO ISO3 codes, as callers already
# have them) is translated to its deduplicated set of EXIO regions; any
# RoW-mapped countries in the input that share a region collapse to one bar.
# Named (non-RoW) countries are unaffected in substance -- their region IS
# just them, under a name instead of an ISO3 code.
make_protein_conversion_plot = function(countries, country_order = countries, show_protein_dots = FALSE) {
  regions_in   = countries_to_regions(countries)
  region_order = countries_to_regions(country_order)

  hr_conv_total = bind_rows(summary_food_region_df_long, summary_nonfood_df_long) %>%
    filter(type %in% c("hr_m", "hr_f"), footprint_type %in% c("domestic_per_capita", "import_per_capita"),
           exio_region %in% regions_in) %>%
    mutate(exio_region = as.character(exio_region)) %>%
    group_by(exio_region, footprint_type) %>%
    summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

  en_conv_total = bind_rows(summary_food_region_df_long, summary_nonfood_df_long) %>%
    filter(type == "en", footprint_type %in% c("domestic_per_capita", "import_per_capita"),
           exio_region %in% regions_in) %>%
    mutate(exio_region = as.character(exio_region)) %>%
    group_by(exio_region, footprint_type) %>%
    summarise(gj_per_cap_day = sum(per_capita_value, na.rm = TRUE) / 365 / 1000, .groups = "drop")

  pro_conv_total = summary_pro_region_df_long %>%
    filter(footprint_type %in% c("domestic_per_capita", "import_per_capita"),
           exio_region %in% regions_in) %>%
    mutate(exio_region = as.character(exio_region)) %>%
    select(exio_region, footprint_type, g_protein_per_cap_day = per_capita_value)

  df_conv = hr_conv_total %>%
    left_join(en_conv_total, by = c("exio_region", "footprint_type")) %>%
    left_join(pro_conv_total, by = c("exio_region", "footprint_type")) %>%
    drop_na() %>%
    mutate(source = ifelse(footprint_type == "domestic_per_capita", "Domestic", "Import"),
           time_conv   = g_protein_per_cap_day / hr_per_cap_day,   # g protein / hr
           energy_conv = g_protein_per_cap_day / gj_per_cap_day,   # g protein / GJ
           is_row = grepl("RoW", exio_region))

  df_conv_long = df_conv %>%
    select(exio_region, source, is_row, time_conv, energy_conv) %>%
    pivot_longer(cols = c(time_conv, energy_conv), names_to = "factor_type", values_to = "value") %>%
    mutate(exio_region = factor(exio_region, levels = as.character(region_order)),
           source = factor(source, levels = c("Domestic", "Import")))

  dot_df = df_conv %>%
    mutate(exio_region = factor(exio_region, levels = as.character(region_order)),
           source = factor(source, levels = c("Domestic", "Import")))
  pro_max = max(dot_df$g_protein_per_cap_day, na.rm = TRUE)

  p_conv_time = ggplot(df_conv_long %>% filter(factor_type == "time_conv"),
                       aes(x = exio_region, y = value, fill = source, alpha = is_row)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Domestic" = "#1f77b4", "Import" = "#ff7f0e")) +
    scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 0.9), guide = "none") +
    labs(x = "EXIO region", y = "g protein / hr", fill = "") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p_conv_energy = ggplot(df_conv_long %>% filter(factor_type == "energy_conv"),
                         aes(x = exio_region, y = value, fill = source, alpha = is_row)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Domestic" = "#1f77b4", "Import" = "#ff7f0e")) +
    scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 0.9), guide = "none") +
    labs(x = "EXIO region", y = "g protein / GJ", fill = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.ticks.x = element_line())

  if (show_protein_dots) {
    time_max = max((df_conv_long %>% filter(factor_type == "time_conv"))$value, na.rm = TRUE)
    pro_scale_time = time_max / pro_max
    p_conv_time = p_conv_time +
      geom_point(data = dot_df, aes(x = exio_region, y = g_protein_per_cap_day * pro_scale_time, group = source),
                 position = position_dodge(width = 0.9), color = "black", size = 2, inherit.aes = FALSE) +
      scale_y_continuous(sec.axis = sec_axis(~ . / pro_scale_time, name = "g protein/cap/day"))

    energy_max = max((df_conv_long %>% filter(factor_type == "energy_conv"))$value, na.rm = TRUE)
    pro_scale_energy = energy_max / pro_max
    p_conv_energy = p_conv_energy +
      geom_point(data = dot_df, aes(x = exio_region, y = g_protein_per_cap_day * pro_scale_energy, group = source),
                 position = position_dodge(width = 0.9), color = "black", size = 2, inherit.aes = FALSE) +
      scale_y_continuous(sec.axis = sec_axis(~ . / pro_scale_energy, name = "g protein/cap/day"))
  }

  p_conv_time / p_conv_energy + plot_layout(guides = "collect") &
    theme(legend.position = "top")
}

# All countries shown in the dual-axis figure (footprint_partial_countries.pdf)
p_conv_combined = make_protein_conversion_plot(partial_cty, partial_ord)
print(p_conv_combined)
ggsave("results/protein_conversion_domestic_vs_import.pdf", p_conv_combined, width = 18, height = 10)


# Derive time conversion factors
# For each country, divide the total food-related time footprint (hr/cap/day) by the total kcal supply per capita (kcal/cap/day) to get hr/kcal. Do this separately for domestic, export, and import footprints.
summary_time_kcal = summary_food_df_long_with_ghd %>% 
  filter(type %in% c("hr_m", "hr_f")) %>%
  mutate(cat = case_when(
    footprint_type == "export_per_capita" ~ "export",
    footprint_type == "import_per_capita" ~ "import",
    .default = "domestic"
  )) %>%
  select(country, type, cat, footprint_type, per_capita_value) %>%
  # pivot_wider(names_from = footprint_type, values_from = per_capita_value) %>%
  left_join(summary_kcal_df_long %>% select(country, cat, footprint_type, per_capita_value),
              # pivot_wider(names_from = footprint_type, values_from = per_capita_value),
            by = c("country", "cat"), suffix = c("_time", "_kcal")) %>% ungroup() %>%
  mutate(hr_per_2000kcal = per_capita_value_time / per_capita_value_kcal * 2e3)


# Separate view for all countries vs. those with preparation_non.econ
# All countries
df_convfac_kcal_econlabor = summary_time_kcal %>% filter(grepl("per_capita", footprint_type_time)) 

# Countries with non-economic observations - this will be a subset of df_convfac_kcal_econlabor but with additional rows for non-economic footprint types (preparation_non.econ, processing_non.econ, growth_collection_non.econ)
df_convfac_kcal_nonecon = summary_time_kcal %>% filter(country %in% cty_ghd) 

# Plot the distribution of hr/kcal conversion factors by country and by type (domestic, export, import) with different colors for type and different facets for hr_m and hr_f.
# Order countries by domestic_hr_per_2000kcal of hr_f
v_ord_econlabor = (df_convfac_kcal_econlabor %>%   
           filter(type %in% c("hr_f"), cat=="domestic") %>% 
           # group_by(country) %>% 
           arrange(-hr_per_2000kcal))$country
v_ord_alltime = (df_convfac_kcal_nonecon %>% 
               filter(type %in% c("hr_f"), cat=="domestic") %>% 
               group_by(country) %>% summarise(d = sum(hr_per_2000kcal, na.rm=TRUE)) %>%
               arrange(-d))$country

# Plot only the economic conversion factors first, and then add the non-economic ones as a separate plot. This way we can see the difference between the two and also avoid having too many bars in one plot.
p_conversion_kcal_econlabor = ggplot(df_convfac_kcal_econlabor %>% 
                                  # For certain export, this can lead to very high hr_per_2000kcal values which can make the plot hard to read. So we can filter out those extreme values for better visualization.
                                  filter(footprint_type_time=="domestic_per_capita") %>% 
                                  select(country, type, footprint_type_time, hr_per_2000kcal) %>%
                                  mutate(country = factor(country, levels = v_ord_econlabor)),
                                aes(x=country, y=hr_per_2000kcal, fill=footprint_type_time)) +
  geom_bar(stat="identity", position="stack") +
  facet_wrap(~type, ncol=1, scales = "fixed") +
  labs(x="Country (ISO3)", y="Time per 2000 kcal (hr/2000kcal)", fill="Footprint type",
       title=paste0("Time per 2000 kcal conversion factors by country (", year, ") - Economic time only")) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_fill_manual(values=c("domestic_per_capita"="#1f77b4", "export_per_capita"="#2ca02c", "import_per_capita"="#ff7f0e")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# And have different facets for hr_m and hr_f
p_conversion_kcal_nonecon = ggplot(df_convfac_kcal_nonecon %>%
                                select(country, type, footprint_type_time, hr_per_2000kcal) %>%
                                mutate(country = factor(country, levels = v_ord_alltime)),
                              aes(x=country, y=hr_per_2000kcal, fill=footprint_type_time)) +
  geom_bar(stat="identity", position="stack") +
  facet_wrap(~type, ncol=1, scales = "fixed") +
  labs(x="Country (ISO3)", y="Time per 2000 kcal (hr/2000kcal)", fill="Footprint type",
       title=paste0("Time per 2000 kcal conversion factors by country (", year, ")")) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_fill_manual(values=c("domestic_per_capita"="#1f77b4", "export_per_capita"="#2ca02c", "import_hr_per_capita"="#ff7f0e",
                             "preparation_non.econ" = "#b41f87",
                             "processing_non.econ" = "#f542f5", 
                             "growth_collection_non.econ" = "#bc36dd", 
                             "preparation_econ" = "#1f2eb4")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/conversion_kcal_econlabor.pdf",    p_conversion_kcal_econlabor,    width = 18, height = 8)
ggsave("results/conversion_kcal_nonecon.pdf", p_conversion_kcal_nonecon, width = 18, height = 8)

# v_ord.tot = (summary_time_kcal %>% drop_na() %>%
#            filter(type %in% c("hr_f")) %>% 
#            group_by(country) %>% summarise(d = sum(hr_per_2000kcal, na.rm=TRUE)) %>%
#            arrange(-hr_per_2000kcal))$country
# p_conversion_kcal = ggplot(summary_time_kcal %>% drop_na() %>%
#                              select(country, type, per_capita_value_time) %>%
#                              pivot_longer(cols = starts_with(c("domestic", "export", "import")), 
#                                           names_to = "footprint_type", values_to = "hr_per_2000kcal") %>% 
#                              mutate(country = factor(country, levels = v_ord.tot)),
#                            aes(x=country, y=hr_per_2000kcal, fill=footprint_type)) +
#   geom_bar(stat="identity", position="stack") +
#   facet_wrap(~type, ncol=1, scales = "fixed") +
#   labs(x="Country (ISO3)", y="Time per 2000 kcal (hr/2000kcal)", fill="Footprint type",
#        title=paste0("Time per 2000 kcal conversion factors by country (", year, ")")) +
#   theme_minimal() +
#   theme(legend.position = "top") +
#   scale_fill_manual(values=c("domestic_hr_per_2000kcal"="#1f77b4", "export_hr_per_2000kcal"="#2ca02c", "import_hr_per_2000kcal"="#ff7f0e")) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))      


# Do the same for protein (long format, mirrors summary_time_kcal)
summary_time_protein = summary_food_df_long_with_ghd %>%
  filter(type %in% c("hr_m", "hr_f")) %>%
  mutate(cat = case_when(
    footprint_type == "export_per_capita" ~ "export",
    footprint_type == "import_per_capita" ~ "import",
    .default = "domestic"
  )) %>%
  select(country, type, cat, footprint_type, per_capita_value) %>%
  left_join(summary_pro_df_long %>% select(country, cat, footprint_type, per_capita_value),
            by = c("country", "cat"), suffix = c("_time", "_protein")) %>%
  ungroup() %>%
  mutate(hr_per_50g_protein = per_capita_value_time / per_capita_value_protein * 50)

v_ord_protein_econlabor = (summary_time_protein %>%
  filter(type == "hr_f", cat == "domestic", grepl("per_capita", footprint_type_time)) %>%
  arrange(-hr_per_50g_protein))$country

df_convfac_protein_econlabor    = summary_time_protein %>% filter(grepl("per_capita", footprint_type_time))
df_convfac_protein_nonecon = summary_time_protein %>% filter(country %in% cty_ghd)

v_ord_protein_alltime = (df_convfac_protein_nonecon %>% 
               filter(type %in% c("hr_f"), cat=="domestic") %>% 
               group_by(country) %>% summarise(d = sum(hr_per_50g_protein, na.rm=TRUE)) %>%
               arrange(-d))$country

p_conversion_protein_econlabor = ggplot(
  df_convfac_protein_econlabor %>%
    filter(footprint_type_time == "domestic_per_capita") %>%
    select(country, type, footprint_type_time, hr_per_50g_protein) %>%
    mutate(country = factor(country, levels = v_ord_protein_econlabor)),
  aes(x = country, y = hr_per_50g_protein, fill = footprint_type_time)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~type, ncol = 1, scales = "fixed") +
  labs(x = "Country (ISO3)", y = "Time per 50 g protein (hr/50g)", fill = "Footprint type",
       title = paste0("Time per 50 g protein by country (", year, ") - Economic time only")) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("domestic_per_capita" = "#1f77b4",
                               "export_per_capita"   = "#2ca02c",
                               "import_per_capita"   = "#ff7f0e")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_conversion_protein_nonecon = ggplot(
  df_convfac_protein_nonecon %>%
    select(country, type, footprint_type_time, hr_per_50g_protein) %>%
    mutate(country = factor(country, levels = v_ord_alltime)),
  aes(x = country, y = hr_per_50g_protein, fill = footprint_type_time)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~type, ncol = 1, scales = "fixed") +
  labs(x = "Country (ISO3)", y = "Time per 50 g protein (hr/50g)", fill = "Footprint type",
       title = paste0("Time per 50 g protein by country (", year, ")")) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("domestic_per_capita"        = "#1f77b4",
                               "export_per_capita"          = "#2ca02c",
                               "import_per_capita"          = "#ff7f0e",
                               "preparation_non.econ"       = "#ffa6a6",
                               "processing_non.econ"        = "#fc4a4a",
                               "growth_collection_non.econ" = "#ce0303",
                               "preparation_econ"           = "#8610ca")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/conversion_protein_econlabor.pdf",    p_conversion_protein_econlabor,    width = 18, height = 8)
ggsave("results/conversion_protein_nonecon.pdf", p_conversion_protein_nonecon, width = 18, height = 8)


#### Paid vs. unpaid food provisioning time by country and gender (scatter) ####

# "Paid" = economic (EXIOBASE-derived) labor embodied in food provisioning: food-sector labor
# (excludes non-food-sector/indirect labor — packaging, transport, etc.), plus restaurant/
# food-service time (preparation_econ), split into domestic (domestic_per_capita, already net
# of exported labor: = territorial production consumed at home, i.e. territorial production
# minus what's exported and excluded here) and import. "Unpaid" = household non-economic time
# from GHD: food preparation/processing/growth-collection plus household water and firewood/
# energy collection (see commit 96d7e3c).

# Paid time is pooled across gender (total economic labor) since labor-market hours aren't
# attributed to the unpaid worker's own gender; only unpaid (household) time is gender-split.
paid_domestic = summary_food_df_long_with_ghd %>%
  filter(type %in% c("hr_m", "hr_f"),
         footprint_type %in% c("domestic_per_capita", "preparation_econ")) %>%
  mutate(country = as.character(country)) %>%
  group_by(country) %>%
  summarise(paid_value = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  mutate(source = "Domestic")

paid_import = summary_food_df_long_with_ghd %>%
  filter(type %in% c("hr_m", "hr_f"), footprint_type == "import_per_capita") %>%
  mutate(country = as.character(country)) %>%
  group_by(country) %>%
  summarise(paid_value = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  mutate(source = "Import")

unpaid = summary_food_df_long_with_ghd %>%
  filter(type %in% c("hr_m", "hr_f"),
         footprint_type %in% c("preparation_non.econ", "processing_non.econ",
                                "growth_collection_non.econ", "energy_non.econ", "water_non.econ")) %>%
  mutate(country = as.character(country)) %>%
  group_by(country, type) %>%
  summarise(unpaid_value = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

# One row per country x gender x source (domestic/import); paid_value (gender-pooled) is
# repeated across the Female/Male facets of the same source. Restricted to non-RoW countries
# (directly modeled EXIO regions only, not aggregated "Rest of" regions).
df_paid_unpaid = bind_rows(paid_domestic, paid_import) %>%
  left_join(unpaid, by = "country") %>%
  filter(country %in% cty_ghd, !country %in% row_countries) %>%
  drop_na(unpaid_value) %>%
  mutate(source = factor(source, levels = c("Domestic", "Import")),
         type = factor(type, levels = c("hr_f", "hr_m"), labels = c("Female", "Male")))

# R^2 of paid ~ unpaid, computed separately per facet (source x gender)
r2_paid_unpaid = df_paid_unpaid %>%
  group_by(source, type) %>%
  summarise(r2 = summary(lm(paid_value ~ unpaid_value))$r.squared, .groups = "drop") %>%
  mutate(label = sprintf("R² = %.2f", r2))

p_paid_unpaid = ggplot(df_paid_unpaid, aes(x = unpaid_value, y = paid_value)) +
  geom_point(color = "#1f77b4", size = 2.5) +
  # Linear fit line: female facets only
  geom_smooth(data = df_paid_unpaid %>% filter(type == "Female"),
              method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  ggrepel::geom_text_repel(aes(label = country), size = 3, max.overlaps = 20, show.legend = FALSE) +
  # R^2 label: all facets
  geom_text(data = r2_paid_unpaid, aes(x = Inf, y = Inf, label = label),
             hjust = 1.1, vjust = 1.5, inherit.aes = FALSE, size = 3.5, color = "black") +
  facet_grid(source ~ type) +
  scale_x_continuous(limits = c(0, NA)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Unpaid (household) time (hr/cap/day)", y = "Paid (economic) time, total (hr/cap/day)",
       title = paste0("Paid (total) vs. unpaid food provisioning time by country and gender (", year, ") — non-RoW")) +
  theme_minimal()
dev.new(width = 12, height = 10)
print(p_paid_unpaid)

ggsave("results/paid_vs_unpaid_time update.pdf", p_paid_unpaid, width = 12, height = 10)


#### Country spotlight ####

# Combines food+nonfood, so the non-food (and, to match it, nutrient) side
# reads from the region-resolved tables -- see chat discussion. iso is
# translated to its own EXIO region; for a named (non-RoW) country this is a
# 1:1, lossless lookup, so current call sites (all named countries) are
# unaffected in substance.
make_spotlight_data = function(iso) {
  region = iso_to_region[[iso]]
  time_energy = bind_rows(
    summary_food_region_df_long %>% mutate(sector = "food sector"),
    summary_nonfood_df_long     %>% mutate(sector = "non-food sector")
  ) %>%
    filter(as.character(exio_region) == region) %>%
    mutate(footprint_type = as.character(footprint_type),
           per_capita_value = ifelse(type == "en", per_capita_value / 1e3, per_capita_value))

  nutrition = bind_rows(
    summary_kcal_region_df_long %>% mutate(type = "kcal",    sector = "—"),
    summary_pro_region_df_long  %>% mutate(type = "protein", sector = "—")
  ) %>%
    filter(as.character(exio_region) == region) %>%
    select(type, footprint_type, per_capita_value, sector) %>%
    mutate(footprint_type = as.character(footprint_type))

  bind_rows(time_energy, nutrition) %>%
    mutate(
      flow = case_when(
        grepl("export", footprint_type) ~ "export",
        grepl("import", footprint_type) ~ "import",
        TRUE ~ "domestic"
      ),
      metric = factor(type,
        levels = c("hr_f", "hr_m", "en", "kcal", "protein"),
        labels = c("Female labor", "Male labor", "Energy", "kcal", "Protein"))
    )
}

mat_kcal_country = agg_country_footprint(FABIO_y_hh_cal)  # 187×187, kcal
mat_pro_country  = agg_country_footprint(FABIO_y_hh_pro)  # 187×187, g protein

# region_continent_of_index is defined near the top-level region lookups
# (shared with b_import's make_import_continent_df, above).

make_continent_data = function(iso) {
  pop = subset(countrypops, year == yr & country_code_3 == iso)$population
  region_idx_iso = region_row_of_country[match(iso, regions$iso3c)]

  # Food sector: l_food_country is still 187x187 (country-resolved, no RoW
  # ambiguity of its own) -- unchanged logic.
  df_food = bind_rows(lapply(names(l_food_country), function(metric) {
    mat  = l_food_country[[metric]]
    cty  = colnames(mat)
    cont = regions$continent[match(cty, regions$iso3c)]
    not_focal = cty != iso

    exp_by_cont = tapply(mat[iso, not_focal],  cont[not_focal], sum, na.rm = TRUE)
    imp_by_cont = tapply(mat[not_focal, iso],  cont[not_focal], sum, na.rm = TRUE)

    denom = if (metric == "en") pop / 1e3 else pop / 1e6 * 365

    bind_rows(
      data.frame(continent = "Domestic",          value = mat[iso, iso] / denom,              flow = "domestic"),
      data.frame(continent = names(exp_by_cont), value = as.numeric(exp_by_cont) / denom, flow = "export"),
      data.frame(continent = names(imp_by_cont), value = as.numeric(imp_by_cont) / denom, flow = "import")
    ) %>% mutate(metric = metric)
  })) %>% mutate(sector = "food sector")

  # Non-food: l_nonfood_region is 49 (EXIO region) x 187 (FABIO country).
  # Export: iso's own region row, grouped by DESTINATION continent (still a
  # clean per-FABIO-country lookup). Import: every OTHER region's
  # contribution to iso's own column, grouped by the ORIGIN region's
  # continent (region_continent_of_index, above).
  df_nonfood = bind_rows(lapply(names(l_nonfood_region), function(metric) {
    mat = l_nonfood_region[[metric]]
    cty = colnames(mat)
    cont = regions$continent[match(cty, regions$iso3c)]
    not_focal = cty != iso
    other_region_idx = setdiff(seq_len(n_reg_EXIO), region_idx_iso)

    exp_by_cont = tapply(mat[region_idx_iso, not_focal], cont[not_focal], sum, na.rm = TRUE)
    imp_by_cont = tapply(mat[other_region_idx, iso], region_continent_of_index[other_region_idx], sum, na.rm = TRUE)

    denom = if (metric == "en") pop / 1e3 else pop / 1e6 * 365

    bind_rows(
      data.frame(continent = "Domestic",          value = mat[region_idx_iso, iso] / denom,   flow = "domestic"),
      data.frame(continent = names(exp_by_cont), value = as.numeric(exp_by_cont) / denom, flow = "export"),
      data.frame(continent = names(imp_by_cont), value = as.numeric(imp_by_cont) / denom, flow = "import")
    ) %>% mutate(metric = metric)
  })) %>% mutate(sector = "non-food sector")

  df_el = bind_rows(df_food, df_nonfood) %>%
    mutate(metric = factor(metric,
      levels = c("hr_f", "hr_m", "en"),
      labels = c("Female labor", "Male labor", "Energy")))

  df_nut = bind_rows(lapply(
    list(kcal = mat_kcal_country, protein = mat_pro_country),
    function(mat) {
      cty  = colnames(mat)
      cont = regions$continent[match(cty, regions$iso3c)]
      not_focal = cty != iso

      exp_by_cont = tapply(mat[iso, not_focal],  cont[not_focal], sum, na.rm = TRUE)
      imp_by_cont = tapply(mat[not_focal, iso],  cont[not_focal], sum, na.rm = TRUE)

      denom = pop * 365  # → kcal/cap/day or g/cap/day

      bind_rows(
        data.frame(continent = "Domestic",          value = mat[iso, iso] / denom,              flow = "domestic"),
        data.frame(continent = names(exp_by_cont), value = as.numeric(exp_by_cont) / denom, flow = "export"),
        data.frame(continent = names(imp_by_cont), value = as.numeric(imp_by_cont) / denom, flow = "import")
      )
    }
  ), .id = "metric") %>%
    mutate(sector = "—",
           metric = factor(metric,
             levels = c("kcal", "protein"),
             labels = c("kcal", "Protein")))

  bind_rows(df_el, df_nut)
}

plot_spotlight = function(iso) {
  df = make_spotlight_data(iso)
  ggplot(df, aes(x = flow, y = per_capita_value, fill = sector)) +
    geom_col(width = 0.6, position = "stack") +
    facet_wrap(~metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c("food sector" = "#4e9af1", "non-food sector" = "#f1884e",
                                 "—" = "#aaaaaa")) +
    labs(x = NULL, y = NULL, fill = "Sector",
         title = paste0(iso, " food system footprint per capita (", year, ")")) +
    theme_minimal() +
    theme(legend.position = "top", strip.text = element_text(size = 10))
}

plot_continent = function(iso) {
  df = make_continent_data(iso)

  cont_levels = sort(setdiff(unique(df$continent), "Domestic"))
  fill_vals = c(
    setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(cont_levels)), cont_levels),
    "Domestic" = "grey65"
  )

  base_theme = list(
    scale_fill_manual(values = fill_vals, name = "Continent"),
    labs(x = NULL, y = NULL),
    theme_minimal(),
    theme(strip.text = element_text(size = 10))
  )

  p_el = df %>%
    filter(sector != "—") %>%
    ggplot(aes(x = flow, y = value, fill = continent)) +
    geom_col(width = 0.6, position = "stack") +
    facet_grid(metric ~ sector, scales = "free_y") +
    base_theme

  p_nut = df %>%
    filter(sector == "—") %>%
    ggplot(aes(x = flow, y = value, fill = continent)) +
    geom_col(width = 0.6, position = "stack") +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    base_theme

  (p_el / p_nut) +
    plot_layout(heights = c(3, 2), guides = "collect") &
    theme(legend.position = "right") &
    plot_annotation(title = paste0(iso, " food trade flows by continent (", year, ")"))
}

for (iso in c("BRA", "ZAF", "DNK", "KOR")) {
# for (iso in c("USA", "CHN", "IND", "AUT")) {
  ggsave(paste0("results/", tolower(iso), "_spotlight update.pdf"),  plot_spotlight(iso), width = 14, height = 5)
  ggsave(paste0("results/", tolower(iso), "_continent update.pdf"),  plot_continent(iso), width = 12, height = 14)
}

# Combined three-country comparison: absolute footprints + domestic/import shares
make_select_country_data = function(isos = c("USA", "AUT", "CHN", "IND", "BRA"), include_export = FALSE) {
  flows = if (include_export) c("domestic", "import", "export") else c("domestic", "import")
  bind_rows(lapply(isos, function(iso) {
    make_spotlight_data(iso) %>% mutate(country = iso)
  })) %>%
    filter(flow %in% flows) %>%
    mutate(country = factor(country, levels = isos),
           flow    = factor(flow, levels = flows))
}

plot_select_country = function(isos = c("USA", "AUT", "CHN", "IND", "BRA")) {
  df = make_select_country_data(isos)

  fill_vals = c(
    "domestic | food sector"     = "#4e9af1",
    "domestic | non-food sector" = "#f1884e",
    "domestic | —"               = "#2ca02c",
    "import | food sector"       = "#91c4f8",
    "import | non-food sector"   = "#f8ba91",
    "import | —"                 = "#9467bd"
  )
  fill_labels = c(
    "domestic | food sector"     = "Domestic: food",
    "domestic | non-food sector" = "Domestic: non-food",
    "domestic | —"               = "Domestic",
    "import | food sector"       = "Import: food",
    "import | non-food sector"   = "Import: non-food",
    "import | —"                 = "Import"
  )
  flow_cols = c("domestic" = "#2ca02c", "import" = "#9467bd")

  df_abs = df %>%
    mutate(fill_grp = factor(paste(flow, sector, sep = " | "), levels = names(fill_vals))) %>%
    group_by(country, fill_grp, metric) %>%
    summarise(value = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

  p_abs = ggplot(df_abs, aes(x = country, y = value, fill = fill_grp)) +
    geom_col(width = 0.6, position = "stack") +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = fill_vals, labels = fill_labels, name = NULL) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(legend.position = "right", strip.text = element_text(size = 9))

  df_share = df %>%
    group_by(country, flow, metric) %>%
    summarise(value = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
    group_by(country, metric) %>%
    mutate(share = value / sum(value) * 100) %>%
    ungroup()

  p_share = ggplot(df_share, aes(x = country, y = share, fill = flow)) +
    geom_col(width = 0.6, position = "stack") +
    facet_wrap(~ metric, nrow = 1) +
    scale_fill_manual(values = flow_cols, name = "Flow") +
    scale_y_continuous(labels = scales::percent_format(scale = 1),
                       breaks = c(0, 25, 50, 75, 100)) +
    labs(x = NULL, y = "Share (%)") +
    theme_minimal() +
    theme(legend.position = "right", strip.text = element_blank())

  (p_abs / p_share) +
    plot_layout(heights = c(3, 1)) +
    plot_annotation(title = paste0("Selected country comparison: domestic vs import footprints (", year, ")"))
}

p_select_country = plot_select_country()
ggsave("results/selected_country_comparison.pdf", p_select_country, width = 18, height = 10)

# Country-arranged comparison with dual y-axis
plot_select_country_dual = function(isos = c("USA", "AUT", "CHN", "IND", "BRA", "ZAF"),
                                    mode       = c("labor_energy", "kcal_protein"),
                                    bar_width  = 0.6,
                                    show_strip = TRUE) {
  mode = match.arg(mode)
  df = make_select_country_data(isos, include_export = TRUE)

  fill_vals = c(
    "domestic | food sector"     = "#4e9af1",
    "domestic | non-food sector" = "#f1884e",
    "domestic | —"               = "#2ca02c",
    "import | food sector"       = "#91c4f8",
    "import | non-food sector"   = "#f8ba91",
    "import | —"                 = "#9467bd",
    "export | food sector"       = "#c0392b",
    "export | non-food sector"   = "#f1948a",
    "export | —"                 = "#922b21"
  )
  fill_labels = c(
    "domestic | food sector"     = "Domestic: food",
    "domestic | non-food sector" = "Domestic: non-food",
    "domestic | —"               = "Domestic",
    "import | food sector"       = "Import: food",
    "import | non-food sector"   = "Import: non-food",
    "import | —"                 = "Import",
    "export | food sector"       = "Export: food",
    "export | non-food sector"   = "Export: non-food",
    "export | —"                 = "Export"
  )

  df_abs = df %>%
    mutate(fill_grp = factor(paste(flow, sector, sep = " | "), levels = names(fill_vals))) %>%
    group_by(country, fill_grp, metric) %>%
    summarise(value = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
    mutate(value = ifelse(grepl("^export", as.character(fill_grp)), -value, value))

  if (mode == "labor_energy") {
    primary_pat   = "labor";       secondary_pat = "Energy"
    y_primary_lab = "Paid embodied labor only (hr/cap/day)"
    y_second_lab  = "Embodied Energy (GJ/cap/yr)"
  } else {
    primary_pat   = "kcal"; secondary_pat = "Protein"
    y_primary_lab = "Food supply (kcal/cap/day)"
    y_second_lab  = "Protein supply (g/cap/day)"
  }

  df_sel = df_abs %>% filter(grepl(paste(primary_pat, secondary_pat, sep = "|"), as.character(metric)))

  totals = df_sel %>% group_by(country, metric) %>% summarise(tot = sum(value), .groups = "drop")
  scale_fac = max(totals$tot[grepl(primary_pat,   as.character(totals$metric))], na.rm = TRUE) /
              max(totals$tot[grepl(secondary_pat, as.character(totals$metric))], na.rm = TRUE)

  df_scaled = df_sel %>%
    mutate(y_plot = ifelse(grepl(secondary_pat, as.character(metric)), value * scale_fac, value))

  ggplot(df_scaled, aes(x = metric, y = y_plot, fill = fill_grp)) +
    geom_col(width = bar_width, position = "stack") +
    facet_wrap(~country, nrow = 1) +
    scale_fill_manual(values = fill_vals, labels = fill_labels, name = NULL) +
    scale_y_continuous(
      name     = y_primary_lab,
      sec.axis = sec_axis(~ . / scale_fac, name = y_second_lab)
    ) +
    labs(x = NULL) +
    theme_minimal() +
    theme(legend.position   = "right",
          strip.text        = if (show_strip) element_text(size = 13) else element_blank(),
          axis.text         = element_text(size = 12),
          axis.text.x       = element_text(size = 12, angle = 45, hjust = 1),
          axis.title        = element_text(size = 15),
          axis.title.y.right = element_text(color = "grey50", size = 15),
          panel.spacing     = unit(1.2, "cm"))
}

# labor_energy: 3 metrics, bar_width = 0.85; kcal_protein: 2 metrics, bar_width = 0.57 → same physical width
# strip shown only on bottom plot (top of its panel = seam between plots)
p_select_labor_energy = plot_select_country_dual(mode = "labor_energy", bar_width = 0.85, show_strip = FALSE)
p_select_kcal_protein  = plot_select_country_dual(mode = "kcal_protein", bar_width = 0.57, show_strip = TRUE)
p_select_combined = (p_select_labor_energy / p_select_kcal_protein) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.text = element_text(size = 12))
# Kept under a separate filename since selected_country_dual_axis.pdf below now
# points at the mosaic/marimekko version instead.
ggsave("results/selected_country_bars.pdf", p_select_combined, width = 15, height = 10)

# Marimekko-style mosaic: one facet per country, two x-spans (Domestic/Import
# protein, width = protein supply g/cap/day for that y_hh-based bucket, see
# pro_width below), each split into four STACKED rectangles: food/non-food
# sector x domestic/import effort-origin (from effort_consumption_df, section
# 1.3 above — the domestic-effort/import-effort split computed on the *same*
# y_hh-defined consumption slice as the bucket it sits in, unlike the
# footprint_type-matched join this mosaic used before). Height = conversion
# factor (effort / protein for that bucket), so rectangle AREA = per-capita
# effort for that sector x effort-origin x bucket.
mosaic_isos = c("USA", "FRA", "AUS", "CHN", "IDN", "IND", "BRA", "ZAF")

mosaic_fill_levels = c("food | Domestic effort", "food | Import effort",
                       "non-food | Domestic effort", "non-food | Import effort")
mosaic_fill_colors = c("food | Domestic effort" = "#4e9af1", "food | Import effort" = "#91c4f8",
                       "non-food | Domestic effort" = "#f1884e", "non-food | Import effort" = "#f8ba91")

pro_mosaic = summary_pro_df_long %>%
  filter(footprint_type %in% c("domestic_per_capita", "import_per_capita"),
         country %in% mosaic_isos) %>%
  mutate(country = as.character(country)) %>%
  select(country, footprint_type, g_protein_per_cap_day = per_capita_value)

pro_width = pro_mosaic %>%
  mutate(source = ifelse(footprint_type == "domestic_per_capita", "Domestic", "Import"),
         country = factor(country, levels = mosaic_isos),
         source = factor(source, levels = c("Domestic", "Import"))) %>%
  arrange(country, source) %>%
  group_by(country) %>%
  mutate(xmax = cumsum(g_protein_per_cap_day),
         xmin = xmax - g_protein_per_cap_day) %>%
  ungroup() %>%
  select(country, footprint_type, source, xmin, xmax, g_protein_per_cap_day)

# effort_consumption_df is exio_region-keyed; mosaic_isos are all named
# (non-RoW) countries, so translate straight back to ISO3 (1:1, lossless) to
# keep the rest of this mosaic's country-keyed joins/faceting unchanged.
hr_mosaic_sector = effort_consumption_df %>%
  mutate(country = region_to_iso_1to1[as.character(exio_region)]) %>%
  filter(type %in% c("hr_m", "hr_f"), country %in% mosaic_isos) %>%
  mutate(footprint_type = ifelse(protein_source == "domestic", "domestic_per_capita", "import_per_capita")) %>%
  group_by(country, footprint_type, sector, effort_origin) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

df_mosaic = hr_mosaic_sector %>%
  left_join(pro_width, by = c("country", "footprint_type")) %>%
  drop_na() %>%
  mutate(min_per_g = hr_per_cap_day * 60 / g_protein_per_cap_day,
         sector = factor(sector, levels = c("food", "non-food")),
         effort_label = ifelse(effort_origin == "domestic", "Domestic effort", "Import effort"),
         fill_grp = factor(paste(sector, effort_label, sep = " | "), levels = mosaic_fill_levels)) %>%
  arrange(country, source, sector, effort_label) %>%
  group_by(country, source) %>%
  mutate(ymax = cumsum(min_per_g),
         ymin = ymax - min_per_g) %>%
  ungroup() %>%
  mutate(min_per_cap_day = hr_per_cap_day * 60,
         label_min = ifelse(round(min_per_cap_day) == 0, NA_character_,
                            paste0(round(min_per_cap_day), " min")),
         fits_inside = (xmax - xmin) >= 0.05 * max(xmax, na.rm = TRUE) &
                       (ymax - ymin) >= 0.06 * max(ymax, na.rm = TRUE),
         label_y = ifelse(fits_inside, (ymin + ymax) / 2, ymax + 0.03 * max(ymax, na.rm = TRUE)),
         # left_join can coerce the factor "country" column back to character;
         # re-assert factor order explicitly so facets follow mosaic_isos exactly.
         country = factor(as.character(country), levels = mosaic_isos))

# Highest total (domestic + import) protein supply among the selected countries;
# drawn on every facet except USA's own (since USA is the country that sets it),
# with a "USA" label so it's clear whose value the line marks.
protein_max = max(df_mosaic$xmax, na.rm = TRUE)
vline_df = tibble(country = factor(setdiff(mosaic_isos, "USA"), levels = mosaic_isos))

title_theme = theme(plot.title = element_text(size = 26, hjust = 0.5, face = "bold",
                                              margin = margin(b = 20)))

p_mosaic = ggplot(df_mosaic) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_grp),
            color = "white", linewidth = 0.4) +
  geom_text(aes(x = (xmin + xmax) / 2, y = label_y, label = label_min),
            color = "black", size = 4, fontface = "bold", na.rm = TRUE) +
  geom_vline(data = vline_df, aes(xintercept = protein_max),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text(data = vline_df, aes(x = protein_max, y = Inf, label = "USA"),
            color = "black", size = 4, fontface = "bold", hjust = 1.1, vjust = 1.5) +
  facet_wrap(~country, ncol = 3) +
  scale_fill_manual(values = mosaic_fill_colors, name = "") +
  labs(x = "Protein supply (g/cap/day)", y = "Time conversion factor (min / g protein)",
       title = "Daily time embodied in domestic/imported protein provisioning") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.8, "lines"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 18, color = "black", face = "bold")) +
  title_theme
print(p_mosaic)

# Energy counterpart of the mosaic above: same daily protein width as the time mosaic
# (for visual alignment), y-axis = MJ / g protein. The rectangle's literal geometric
# area is technically GJ/day, but the printed label deliberately shows the annualized
# GJ/yr figure instead, for readability.
en_mosaic_sector = effort_consumption_df %>%
  mutate(country = region_to_iso_1to1[as.character(exio_region)]) %>%
  filter(type == "en", country %in% mosaic_isos) %>%
  mutate(footprint_type = ifelse(protein_source == "domestic", "domestic_per_capita", "import_per_capita"),
         gj_per_cap_day = per_capita_value / 365 / 1000) %>%
  select(country, footprint_type, sector, effort_origin, gj_per_cap_day)

df_mosaic_energy = en_mosaic_sector %>%
  left_join(pro_width, by = c("country", "footprint_type")) %>%
  drop_na() %>%
  mutate(mj_per_g = gj_per_cap_day * 1000 / g_protein_per_cap_day,
         sector = factor(sector, levels = c("food", "non-food")),
         effort_label = ifelse(effort_origin == "domestic", "Domestic effort", "Import effort"),
         fill_grp = factor(paste(sector, effort_label, sep = " | "), levels = mosaic_fill_levels)) %>%
  arrange(country, source, sector, effort_label) %>%
  group_by(country, source) %>%
  mutate(ymax = cumsum(mj_per_g),
         ymin = ymax - mj_per_g) %>%
  ungroup() %>%
  mutate(gj_per_year = gj_per_cap_day * 365,
         label_gj = ifelse(round(gj_per_year, 1) == 0, NA_character_,
                            paste0(round(gj_per_year, 1), " GJ/yr")),
         # If a rectangle is too small (in either dimension) for the label to plausibly
         # fit inside it, place the label just above the rectangle instead.
         fits_inside = (xmax - xmin) >= 0.05 * max(xmax, na.rm = TRUE) &
                       (ymax - ymin) >= 0.06 * max(ymax, na.rm = TRUE),
         label_y = ifelse(fits_inside, (ymin + ymax) / 2, ymax + 0.03 * max(ymax, na.rm = TRUE)),
         # left_join can coerce the factor "country" column back to character;
         # re-assert factor order explicitly so facets follow mosaic_isos exactly.
         country = factor(as.character(country), levels = mosaic_isos))

# Same reference line as the time mosaic (same underlying protein widths), skipping
# USA's own facet and labeled "USA" on every other facet.
protein_max_energy = max(df_mosaic_energy$xmax, na.rm = TRUE)
vline_df_energy = tibble(country = factor(setdiff(mosaic_isos, "USA"), levels = mosaic_isos))

p_mosaic_energy = ggplot(df_mosaic_energy) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_grp),
            color = "white", linewidth = 0.4) +
  geom_text(aes(x = (xmin + xmax) / 2, y = label_y, label = label_gj),
            color = "black", size = 4, fontface = "bold", na.rm = TRUE) +
  geom_vline(data = vline_df_energy, aes(xintercept = protein_max_energy),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text(data = vline_df_energy, aes(x = protein_max_energy, y = Inf, label = "USA"),
            color = "black", size = 4, fontface = "bold", hjust = 1.1, vjust = 1.5) +
  facet_wrap(~country, ncol = 3) +
  scale_fill_manual(values = mosaic_fill_colors, name = "") +
  labs(x = "Protein supply (g/cap/day)", y = "Energy conversion factor (MJ / g protein)",
       title = "Yearly energy embodied in domestic/imported protein provisioning") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.8, "lines"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 18, color = "black", face = "bold")) +
  title_theme
print(p_mosaic_energy)

# Combine time + energy side by side with one shared, centered, bottom legend
# (both plots use the same fill scale, so guides="collect" deduplicates it).
# A blank plot_spacer() between them creates the gap; each mosaic spans 3 facet
# columns (ncol = 3) within its "1" width unit, so a gap of 1/12 = (1/4) / 3
# comes out to a quarter of one facet's width.
p_mosaic_combined = (p_mosaic | plot_spacer() | p_mosaic_energy) +
  plot_layout(guides = "collect", widths = c(1, 1/12, 1)) &
  theme(legend.position = "bottom")
print(p_mosaic_combined)

ggsave("results/protein_mosaic_selected_combined update.pdf", p_mosaic_combined, width = 27, height = 11)

# Same mosaic, but for every directly-modeled (non-RoW) country instead of the six selected ones
mosaic_isos_nonrow = setdiff(regions$iso3c, row_countries)

pro_mosaic_nonrow = summary_pro_df_long %>%
  filter(footprint_type %in% c("domestic_per_capita", "import_per_capita"),
         country %in% mosaic_isos_nonrow) %>%
  mutate(country = as.character(country)) %>%
  select(country, footprint_type, g_protein_per_cap_day = per_capita_value)

# Order countries by domestic-only time/cap/day, descending. mosaic_isos_nonrow
# is all-named by construction (setdiff of every RoW-mapped country), so the
# combined food+nonfood total translates back to ISO3 losslessly.
hr_domestic_nonrow = bind_rows(summary_food_region_df_long, summary_nonfood_df_long) %>%
  mutate(country = region_to_iso_1to1[as.character(exio_region)]) %>%
  filter(type %in% c("hr_m", "hr_f"), footprint_type == "domestic_per_capita",
         country %in% mosaic_isos_nonrow) %>%
  group_by(country) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop")
mosaic_nonrow_ord = hr_domestic_nonrow %>% arrange(-hr_per_cap_day) %>% pull(country)

pro_width_nonrow = pro_mosaic_nonrow %>%
  mutate(source = ifelse(footprint_type == "domestic_per_capita", "Domestic", "Import"),
         country = factor(country, levels = mosaic_nonrow_ord),
         source = factor(source, levels = c("Domestic", "Import"))) %>%
  arrange(country, source) %>%
  group_by(country) %>%
  mutate(xmax = cumsum(g_protein_per_cap_day),
         xmin = xmax - g_protein_per_cap_day) %>%
  ungroup() %>%
  select(country, footprint_type, source, xmin, xmax, g_protein_per_cap_day)

hr_mosaic_nonrow_sector = effort_consumption_df %>%
  mutate(country = region_to_iso_1to1[as.character(exio_region)]) %>%
  filter(type %in% c("hr_m", "hr_f"), country %in% mosaic_isos_nonrow) %>%
  mutate(footprint_type = ifelse(protein_source == "domestic", "domestic_per_capita", "import_per_capita")) %>%
  group_by(country, footprint_type, sector, effort_origin) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

df_mosaic_nonrow = hr_mosaic_nonrow_sector %>%
  left_join(pro_width_nonrow, by = c("country", "footprint_type")) %>%
  drop_na() %>%
  mutate(hr_per_g = hr_per_cap_day / g_protein_per_cap_day,
         sector = factor(sector, levels = c("food", "non-food")),
         effort_label = ifelse(effort_origin == "domestic", "Domestic effort", "Import effort"),
         fill_grp = factor(paste(sector, effort_label, sep = " | "), levels = mosaic_fill_levels)) %>%
  arrange(country, source, sector, effort_label) %>%
  group_by(country, source) %>%
  mutate(ymax = cumsum(hr_per_g),
         ymin = ymax - hr_per_g) %>%
  ungroup() %>%
  mutate(min_per_cap_day = hr_per_cap_day * 60,
         label_min = ifelse(round(min_per_cap_day) == 0, NA_character_,
                            paste0(round(min_per_cap_day), " min")),
         fits_inside = (xmax - xmin) >= 0.05 * max(xmax, na.rm = TRUE) &
                       (ymax - ymin) >= 0.06 * max(ymax, na.rm = TRUE),
         label_y = ifelse(fits_inside, (ymin + ymax) / 2, ymax + 0.03 * max(ymax, na.rm = TRUE)),
         label_col = ifelse(fits_inside, "white", "black"))

p_mosaic_nonrow = ggplot(df_mosaic_nonrow) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_grp),
            color = "white", linewidth = 0.3) +
  geom_text(aes(x = (xmin + xmax) / 2, y = label_y, label = label_min, color = I(label_col)),
            size = 3, fontface = "bold", na.rm = TRUE) +
  facet_wrap(~country, ncol = 7) +
  scale_fill_manual(values = mosaic_fill_colors, name = "") +
  labs(x = "Protein supply (g/cap/day)", y = "Time conversion factor (hr / g protein)",
       title = "Daily time embodied in domestic/imported protein provisioning") +
  theme_minimal() +
  theme(legend.position = "top",
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 6))
print(p_mosaic_nonrow)

ggsave("results/protein_time_mosaic_nonrow update.pdf", p_mosaic_nonrow,
       width = 24, height = ceiling(length(mosaic_isos_nonrow) / 7) * 3.2, limitsize = FALSE)

# Energy counterpart of the non-RoW mosaic above: same daily protein width as the time
# mosaic (for alignment), y-axis = MJ / g protein, label deliberately shows annual GJ/yr.
en_mosaic_nonrow_sector = effort_consumption_df %>%
  mutate(country = region_to_iso_1to1[as.character(exio_region)]) %>%
  filter(type == "en", country %in% mosaic_isos_nonrow) %>%
  mutate(footprint_type = ifelse(protein_source == "domestic", "domestic_per_capita", "import_per_capita"),
         gj_per_cap_day = per_capita_value / 365 / 1000) %>%
  select(country, footprint_type, sector, effort_origin, gj_per_cap_day)

df_mosaic_nonrow_energy = en_mosaic_nonrow_sector %>%
  left_join(pro_width_nonrow, by = c("country", "footprint_type")) %>%
  drop_na() %>%
  mutate(mj_per_g = gj_per_cap_day * 1000 / g_protein_per_cap_day,
         sector = factor(sector, levels = c("food", "non-food")),
         effort_label = ifelse(effort_origin == "domestic", "Domestic effort", "Import effort"),
         fill_grp = factor(paste(sector, effort_label, sep = " | "), levels = mosaic_fill_levels)) %>%
  arrange(country, source, sector, effort_label) %>%
  group_by(country, source) %>%
  mutate(ymax = cumsum(mj_per_g),
         ymin = ymax - mj_per_g) %>%
  ungroup() %>%
  mutate(gj_per_year = gj_per_cap_day * 365,
         label_gj = ifelse(round(gj_per_year, 1) == 0, NA_character_,
                            paste0(round(gj_per_year, 1), " GJ/yr")),
         fits_inside = (xmax - xmin) >= 0.05 * max(xmax, na.rm = TRUE) &
                       (ymax - ymin) >= 0.06 * max(ymax, na.rm = TRUE),
         label_y = ifelse(fits_inside, (ymin + ymax) / 2, ymax + 0.03 * max(ymax, na.rm = TRUE)),
         label_col = ifelse(fits_inside, "white", "black"))

p_mosaic_nonrow_energy = ggplot(df_mosaic_nonrow_energy) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_grp),
            color = "white", linewidth = 0.3) +
  geom_text(aes(x = (xmin + xmax) / 2, y = label_y, label = label_gj, color = I(label_col)),
            size = 3, fontface = "bold", na.rm = TRUE) +
  facet_wrap(~country, ncol = 7) +
  scale_fill_manual(values = mosaic_fill_colors, name = "") +
  labs(x = "Protein supply (g/cap/day)", y = "Energy conversion factor (MJ / g protein)",
       title = "Yearly energy embodied in domestic/imported protein provisioning") +
  theme_minimal() +
  theme(legend.position = "top",
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 6))
print(p_mosaic_nonrow_energy)

ggsave("results/protein_energy_mosaic_nonrow update.pdf", p_mosaic_nonrow_energy,
       width = 24, height = ceiling(length(mosaic_isos_nonrow) / 7) * 3.2, limitsize = FALSE)

# Same protein conversion factor comparison (domestic vs. import), but for the five
# countries used in the dual-axis figure above (USA, AUT, CHN, IND, BRA)
p_conv_combined_selected = make_protein_conversion_plot(c("USA", "AUT", "CHN", "IND", "BRA", "ZAF"),
                                                         show_protein_dots = TRUE)
print(p_conv_combined_selected)
ggsave("results/protein_conversion_domestic_vs_import_selected.pdf", p_conv_combined_selected, width = 10, height = 8)


#### Energy-time tradeoff ####
library(ggrepel)

# Energy conversion factors (food sector, per unit nutrition)
# Energy is MJ/cap/year; divide by 365 for daily basis consistent with time (hr/cap/day)
summary_energy_kcal = summary_food_df_long %>%
  filter(type == "en") %>%
  left_join(summary_kcal_df_long %>% select(country, footprint_type, per_capita_value),
            by = c("country", "footprint_type"), suffix = c("_energy", "_kcal")) %>%
  ungroup() %>%
  mutate(mj_per_2000kcal = (per_capita_value_energy / 365) / per_capita_value_kcal * 2000)

summary_energy_protein = summary_food_df_long %>%
  filter(type == "en") %>%
  left_join(summary_pro_df_long %>% select(country, footprint_type, per_capita_value),
            by = c("country", "footprint_type"), suffix = c("_energy", "_protein")) %>%
  ungroup() %>%
  mutate(mj_per_50g_protein = (per_capita_value_energy / 365) / per_capita_value_protein * 50)

# Domestic tradeoff: economic time only (all countries)
tradeoff_kcal_econlabor = summary_time_kcal %>%
  filter(cat == "domestic", footprint_type_time == "domestic_per_capita") %>%
  select(country, type, hr_per_2000kcal) %>%
  left_join(summary_energy_kcal %>%
              filter(footprint_type == "domestic_per_capita") %>%
              select(country, mj_per_2000kcal, kcal_per_cap_day = per_capita_value_kcal),
            by = "country") %>%
  left_join(regions %>% select(iso3c, continent), by = c("country" = "iso3c")) %>%
  drop_na() %>%
  filter(!as.character(country) %in% row_countries)

tradeoff_protein_econlabor = summary_time_protein %>%
  filter(cat == "domestic", grepl("per_capita", footprint_type_time)) %>%
  select(country, type, hr_per_50g_protein) %>%
  left_join(summary_energy_protein %>%
              filter(footprint_type == "domestic_per_capita") %>%
              select(country, mj_per_50g_protein, g_protein_per_cap_day = per_capita_value_protein),
            by = "country") %>%
  left_join(regions %>% select(iso3c, continent), by = c("country" = "iso3c")) %>%
  drop_na() %>%
  filter(!as.character(country) %in% row_countries)

# Domestic tradeoff: total time (economic + non-economic), GHD countries only
# Sum hr_per_2000kcal across all domestic footprint types per country+gender
tradeoff_kcal_allwork = summary_time_kcal %>%
  filter(cat == "domestic", country %in% cty_ghd) %>%
  group_by(country, type) %>%
  summarise(hr_per_2000kcal = sum(hr_per_2000kcal, na.rm = TRUE), .groups = "drop") %>%
  left_join(summary_energy_kcal %>%
              filter(footprint_type == "domestic_per_capita") %>%
              select(country, mj_per_2000kcal, kcal_per_cap_day = per_capita_value_kcal),
            by = "country") %>%
  left_join(regions %>% select(iso3c, continent), by = c("country" = "iso3c")) %>%
  drop_na() %>%
  mutate(is_row = as.character(country) %in% row_countries)

tradeoff_protein_allwork = summary_time_protein %>%
  filter(cat == "domestic", country %in% cty_ghd) %>%
  group_by(country, type) %>%
  summarise(hr_per_50g_protein = sum(hr_per_50g_protein, na.rm = TRUE), .groups = "drop") %>%
  left_join(summary_energy_protein %>%
              filter(footprint_type == "domestic_per_capita") %>%
              select(country, mj_per_50g_protein, g_protein_per_cap_day = per_capita_value_protein),
            by = "country") %>%
  left_join(regions %>% select(iso3c, continent), by = c("country" = "iso3c")) %>%
  drop_na() %>%
  mutate(is_row = as.character(country) %in% row_countries)

# Scatter plot helper: energy (x) vs time (y), sized by nutrition supply per capita
tradeoff_scatter <- function(df, x_col, y_col, size_col, x_lab, y_lab, size_lab, title, label_size = 3.5) {
  has_row <- "is_row" %in% colnames(df)
  base_aes <- aes(x = .data[[x_col]], y = .data[[y_col]], size = .data[[size_col]],
                  color = continent, label = country)
  g <- ggplot(df, base_aes)
  if (has_row) {
    g <- g +
      geom_point(aes(alpha = is_row)) +
      ggrepel::geom_text_repel(aes(alpha = is_row), size = label_size, max.overlaps = 20, show.legend = FALSE) +
      scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 0.9), guide = "none")
  } else {
    g <- g +
      geom_point(alpha = 0.9) +
      ggrepel::geom_text_repel(size = label_size, max.overlaps = 20, show.legend = FALSE)
  }
  g +
    scale_size_continuous(range = c(1, 8)) +
    facet_wrap(~type, nrow = 1,
               labeller = labeller(type = c("hr_m" = "Male", "hr_f" = "Female"))) +
    labs(x = x_lab, y = y_lab, size = size_lab,
         color = "Continent", title = title) +
    theme_minimal() +
    theme(legend.position = "right",
          strip.text   = element_text(size = rel(1.5)),
          axis.title   = element_text(size = rel(1.4)),
          legend.text  = element_text(size = rel(1.3)),
          legend.title = element_text(size = rel(1.3)))
}

p_tradeoff_kcal_econlabor = tradeoff_scatter(
  tradeoff_kcal_econlabor %>% mutate(hr_per_2000kcal = hr_per_2000kcal * 60),
  "mj_per_2000kcal", "hr_per_2000kcal", "kcal_per_cap_day",
  "Energy (MJ / 2000 kcal)", "Time (min / 2000 kcal)", "kcal/cap/day",
  paste0("Energy vs. time to provision 2000 kcal (", year, ") — Econ. labor"))

p_tradeoff_kcal_allwork = tradeoff_scatter(
  tradeoff_kcal_allwork, "mj_per_2000kcal", "hr_per_2000kcal", "kcal_per_cap_day",
  "Energy (MJ / 2000 kcal)", "Time (hr / 2000 kcal)", "kcal/cap/day",
  paste0("Energy vs. time to provision 2000 kcal (", year, ") — Econ. + non-econ labor")) +
  coord_cartesian(ylim = c(NA, 3))

p_tradeoff_protein_econlabor = tradeoff_scatter(
  tradeoff_protein_econlabor %>% mutate(hr_per_50g_protein = hr_per_50g_protein * 60),
  "mj_per_50g_protein", "hr_per_50g_protein", "g_protein_per_cap_day",
  "Energy (MJ / 50 g protein)", "Time (min / 50 g protein)", "g protein/cap/day",
  paste0("Energy vs. time to provision 50 g protein (", year, ") — Econ. labor"))

p_tradeoff_protein_allwork = tradeoff_scatter(
  tradeoff_protein_allwork, "mj_per_50g_protein", "hr_per_50g_protein", "g_protein_per_cap_day",
  "Energy (MJ / 50 g protein)", "Time (hr / 50 g protein)", "g protein/cap/day",
  paste0("Energy vs. time to provision 50 g protein (", year, ") — Econ. + non-econ labor")) +
  coord_cartesian(ylim = c(NA, 6))

# Same domestic protein tradeoff (economic labor only), but collapsing the
# male/female facets into a single panel: each country becomes a vertical
# segment at its (gender-invariant) energy value, running from male hours to
# female hours of paid food-sector labor.
tradeoff_protein_econlabor_min = tradeoff_protein_econlabor %>%
  mutate(hr_per_50g_protein = hr_per_50g_protein * 60)

label_tradeoff_protein_econlabor_gender = tradeoff_protein_econlabor_min %>%
  group_by(country) %>%
  slice_max(hr_per_50g_protein, n = 1, with_ties = FALSE)

p_tradeoff_protein_econlabor_gender_range = ggplot(
  tradeoff_protein_econlabor_min,
  aes(x = mj_per_50g_protein, y = hr_per_50g_protein)) +
  geom_line(aes(group = country), color = "grey60", linewidth = 0.6) +
  geom_point(aes(color = type, size = g_protein_per_cap_day), alpha = 0.85) +
  ggrepel::geom_text_repel(
    data = label_tradeoff_protein_econlabor_gender,
    aes(label = country), size = 3.5, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c(hr_f = "#ca2323", hr_m = "#1f77b4"),
                      labels = c(hr_f = "Female", hr_m = "Male")) +
  scale_size_continuous(range = c(1, 8)) +
  labs(x = "Energy (MJ / 50 g protein)", y = "Time (min / 50 g protein)",
       color = "Gender", size = "g protein/cap/day",
       title = paste0("Energy vs. time to provision 50 g protein (", year, ") — Econ. labor, by gender")) +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title   = element_text(size = rel(1.4)),
        legend.text  = element_text(size = rel(1.1)),
        legend.title = element_text(size = rel(1.1)))

ggsave("results/tradeoff_convfac_protein_econlabor_gender_range.pdf",
       p_tradeoff_protein_econlabor_gender_range, width = 14, height = 9)

p_tradeoff_econlabor = (p_tradeoff_kcal_econlabor | p_tradeoff_protein_econlabor) +
  plot_layout(guides = "collect") & theme(legend.position = "right")

p_tradeoff_allwork = (p_tradeoff_kcal_allwork | p_tradeoff_protein_allwork) +
  plot_layout(guides = "collect") & theme(legend.position = "right")

ggsave(paste0("results/tradeoff_convfac_econlabor.pdf"), p_tradeoff_econlabor, width = 20, height = 8)
ggsave(paste0("results/tradeoff_convfac_allwork.pdf"), p_tradeoff_allwork, width = 20, height = 8)

# Econlabor scatter sized by non-economic time (GHD countries only)
nonecon_size_kcal = summary_time_kcal %>%
  filter(cat == "domestic", country %in% cty_ghd,
         grepl("non.econ", footprint_type_time)) %>%
  group_by(country, type) %>%
  summarise(nonecon_hr_per_2000kcal = sum(hr_per_2000kcal, na.rm = TRUE), .groups = "drop")

tradeoff_kcal_econlabor_noneconsize = tradeoff_kcal_econlabor %>%
  inner_join(nonecon_size_kcal, by = c("country", "type"))

nonecon_size_protein = summary_time_protein %>%
  filter(cat == "domestic", country %in% cty_ghd,
         grepl("non.econ", footprint_type_time)) %>%
  group_by(country, type) %>%
  summarise(nonecon_hr_per_50g_protein = sum(hr_per_50g_protein, na.rm = TRUE), .groups = "drop")

tradeoff_protein_econlabor_noneconsize = tradeoff_protein_econlabor %>%
  inner_join(nonecon_size_protein, by = c("country", "type"))

p_tradeoff_kcal_econlabor_noneconsize = tradeoff_scatter(
  tradeoff_kcal_econlabor_noneconsize %>% mutate(hr_per_2000kcal = hr_per_2000kcal * 60),
  "mj_per_2000kcal", "hr_per_2000kcal", "nonecon_hr_per_2000kcal",
  "Energy (MJ / 2000 kcal)", "Economic time (min / 2000 kcal)", "Non-econ. time (hr / 2000 kcal)",
  paste0("Energy vs. economic time to provision 2000 kcal (", year, ") — sized by non-econ. time"))

p_tradeoff_protein_econlabor_noneconsize = tradeoff_scatter(
  tradeoff_protein_econlabor_noneconsize %>% mutate(hr_per_50g_protein = hr_per_50g_protein * 60),
  "mj_per_50g_protein", "hr_per_50g_protein", "nonecon_hr_per_50g_protein",
  "Energy (MJ / 50 g protein)", "Economic time (min / 50 g protein)", "Non-econ. time (hr / 50 g protein)",
  paste0("Energy vs. economic time to provision 50 g protein (", year, ") — sized by non-econ. time"))

p_tradeoff_econlabor_noneconsize = (p_tradeoff_kcal_econlabor_noneconsize | p_tradeoff_protein_econlabor_noneconsize) +
  plot_layout(guides = "collect") & theme(legend.position = "right")

ggsave(paste0("results/tradeoff_convfac_econlabor_noneconsize update.pdf"), p_tradeoff_econlabor_noneconsize, width = 20, height = 8)

# Protein scatter sized by protein import share: import / (domestic + import)
protein_import_share = summary_pro_df_long %>%
  filter(cat %in% c("domestic", "import")) %>%
  mutate(country = as.character(country)) %>%
  group_by(country) %>%
  summarise(
    protein_import_share = sum(per_capita_value[cat == "import"], na.rm = TRUE) /
      sum(per_capita_value, na.rm = TRUE),
    .groups = "drop"
  )

# Female-only version of p4_nonrow_50g_indirect: one point per country (total
# domestic food + non-food hr/50g protein, already a single row per country
# for type == "hr_f" in df_time_protein_nonrow_indirect), dot size instead
# encodes protein_import_share rather than splitting direct vs. indirect.
df_female_50g_importshare = df_time_protein_nonrow_indirect %>%
  filter(type == "hr_f") %>%
  mutate(country = as.character(country)) %>%
  left_join(protein_import_share, by = "country") %>%
  left_join(table_food_nonfood_50g %>%
              filter(gender == "Female") %>%
              select(country, unpaid_food_hr_per_50g, paid_food_hr_per_50g, nonfood_hr_per_50g),
            by = "country") %>%
  drop_na(protein_import_share) %>%
  mutate(unpaid_share = unpaid_food_hr_per_50g / hr_per_50g_protein_total,
         paid_share = (paid_food_hr_per_50g + nonfood_hr_per_50g) / hr_per_50g_protein_total)

p4_female_50g_importshare = ggplot(df_female_50g_importshare,
                                    aes(x = pro_per_cap_day, y = hr_per_50g_protein_total)) +
  geom_point(aes(size = protein_import_share), color = "#ca2323", alpha = 0.8) +
  ggrepel::geom_text_repel(
    aes(label = country), size = 3, fontface = "bold", max.overlaps = 20, show.legend = FALSE) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
  scale_size_continuous(labels = scales::percent, range = c(1.5, 8)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Total protein supply (g/cap/day)",
       y = "Domestic food provisioning time per 50 g protein (hr/50g)",
       size = "Protein import share",
       title = "Female: paid & unpaid time per 50 g protein, sized by protein import share") +
  theme_minimal()
dev.new(width = 16, height = 6)
print(p4_female_50g_importshare)
ggsave(paste0("results/protein_vs_time_per50g_nonrow_indirect_importshare.pdf"), p4_female_50g_importshare, width = 15, height = 8)

# dev.new()
# ppp= ggplot(df_female_50g_importshare, aes(x = protein_import_share, y = unpaid_share)) +
#   geom_point(color = "#ca2323", size = 2.5) +
#   ggrepel::geom_text_repel(aes(label = country), size = 3, max.overlaps = 20) +
#   scale_x_continuous(labels = scales::percent) +
#   scale_y_continuous(labels = scales::percent) +
#   labs(x = "Protein import share", y = "Unpaid share of time per 50 g protein",
#        title = "Female: unpaid share vs. protein import share") +
#   theme_minimal()
# print(ppp)
# ggsave(paste0("results/ppp.pdf"), ppp, width = 15, height = 8)

# Same female country set/axes as p4_female_50g_importshare (x = protein
# supply, y = hr/50g protein), but each country is a stacked bar (unpaid food
# / paid food / paid non-food) instead of a single dot -- bar width is a
# fraction of the x range since pro_per_cap_day is continuous, not discrete.
df_female_stack = df_female_50g_importshare %>%
  select(country, pro_per_cap_day, unpaid_food_hr_per_50g, paid_food_hr_per_50g, nonfood_hr_per_50g) %>%
  pivot_longer(cols = c(unpaid_food_hr_per_50g, paid_food_hr_per_50g, nonfood_hr_per_50g),
               names_to = "component", values_to = "hr_per_50g_protein") %>%
  mutate(component = factor(component,
                             levels = c("unpaid_food_hr_per_50g", "paid_food_hr_per_50g", "nonfood_hr_per_50g"),
                             labels = c("Unpaid food (household)", "Paid food (economic)", "Paid non-food (indirect)")))

bar_width_female = diff(range(df_female_stack$pro_per_cap_day, na.rm = TRUE)) / 60

p4_female_50g_stack = ggplot(df_female_stack, aes(x = pro_per_cap_day, y = hr_per_50g_protein, fill = component)) +
  geom_col(aes(group = country),
           position = position_dodge2(width = bar_width_female, padding = 0.1, preserve = "single"),
           width = bar_width_female, alpha = 0.85, color = "white", linewidth = 0.15) +
  ggrepel::geom_text_repel(
    data = df_female_50g_importshare,
    aes(x = pro_per_cap_day, y = hr_per_50g_protein_total, label = country),
    size = 3, fontface = "bold", max.overlaps = 20, inherit.aes = FALSE) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Unpaid food (household)" = "#ca2323",
                                "Paid food (economic)" = "#f4a582",
                                "Paid non-food (indirect)" = "#4393c3")) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Total protein supply (g/cap/day)",
       y = "Domestic food provisioning time per 50 g protein (hr/50g)",
       fill = "Component",
       title = "Female: unpaid/paid-food/paid-non-food time per 50 g protein") +
  theme_minimal()
dev.new(width = 16, height = 6)
print(p4_female_50g_stack)
ggsave(paste0("results/protein_vs_time_per50g_nonrow_indirect_stack.pdf"), p4_female_50g_stack, width = 15, height = 8)

tradeoff_protein_econlabor_importshare = tradeoff_protein_econlabor_noneconsize %>%
  left_join(protein_import_share, by = "country") %>%
  drop_na(protein_import_share) %>%
  mutate(hr_per_50g_protein = hr_per_50g_protein * 60)

x_lim_importshare = range(tradeoff_protein_econlabor_importshare$mj_per_50g_protein, na.rm = TRUE)

p_importshare_top = tradeoff_scatter(
  tradeoff_protein_econlabor_importshare, "mj_per_50g_protein", "hr_per_50g_protein", "protein_import_share",
  "Energy (MJ / 50 g protein)", "Economic time (min / 50 g protein)", "Protein import share",
  paste0("Energy vs. economic time to provision 50 g protein (", year, ") — sized by protein import share")) +
  coord_cartesian(xlim = x_lim_importshare, ylim = c(42, 45)) +
  scale_y_continuous(breaks = c(42, 45)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        strip.text = element_text(size = rel(1.5)))
p_importshare_top$layers <- Filter(
  function(l) !inherits(l$geom, "GeomTextRepel"), p_importshare_top$layers)

p_importshare_bottom = tradeoff_scatter(
  tradeoff_protein_econlabor_importshare, "mj_per_50g_protein", "hr_per_50g_protein", "protein_import_share",
  "Energy (MJ / 50 g protein)", "Economic time (min / 50 g protein)", "Protein import share",
  NULL) +
  coord_cartesian(xlim = x_lim_importshare, ylim = c(0, 24)) +
  theme(strip.text = element_blank(), legend.position = "none")
for (i in seq_along(p_importshare_bottom$layers)) {
  if (inherits(p_importshare_bottom$layers[[i]]$geom, "GeomTextRepel")) {
    p_importshare_bottom$layers[[i]]$aes_params$size <- 3.5
    p_importshare_bottom$layers[[i]]$geom_params$point.padding <- 0.8
    p_importshare_bottom$layers[[i]]$geom_params$min.segment.length <- 0
  }
}

# heights ratio = top range / bottom range = 3 / 24 = 1:8 → same physical scale per unit
p_tradeoff_protein_econlabor_importshare = (p_importshare_top / p_importshare_bottom) +
  plot_layout(heights = c(1, 8), guides = "collect") &
  theme(legend.position = "right",
        axis.title = element_text(size = rel(1.4)),
        legend.text = element_text(size = rel(1.3)),
        legend.title = element_text(size = rel(1.3)))

ggsave(paste0("results/tradeoff_convfac_protein_importshare.pdf"), p_tradeoff_protein_econlabor_importshare, width = 15, height = 7)

tradeoff_protein_allwork_importshare = tradeoff_protein_allwork %>%
  mutate(country = as.character(country)) %>%
  left_join(protein_import_share, by = "country") %>%
  drop_na(protein_import_share)

p_tradeoff_protein_allwork_importshare = tradeoff_scatter(
  tradeoff_protein_allwork_importshare, "mj_per_50g_protein", "hr_per_50g_protein", "protein_import_share",
  "Energy (MJ / 50 g protein)", "Time (hr / 50 g protein)", "Protein import share",
  paste0("Energy vs. time to provision 50 g protein (", year, ") — Econ. + non-econ labor, sized by protein import share")) +
  coord_cartesian(ylim = c(NA, 6))

ggsave(paste0("results/tradeoff_convfac_protein_allwork_importshare.pdf"), p_tradeoff_protein_allwork_importshare, width = 12, height = 8)

# Consumption-based tradeoff: total (domestic + import) nutrients and hours
# Shared consumption totals: sum domestic + import, excluding exports
kcal_consumption = summary_kcal_df_long %>%
  filter(cat %in% c("domestic", "import")) %>%
  mutate(country = as.character(country)) %>%
  group_by(country) %>%
  summarise(kcal_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

pro_consumption = summary_pro_df_long %>%
  filter(cat %in% c("domestic", "import")) %>%
  mutate(country = as.character(country)) %>%
  group_by(country) %>%
  summarise(g_protein_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

en_consumption = summary_food_df_long %>%
  filter(type == "en", footprint_type %in% c("domestic_per_capita", "import_per_capita")) %>%
  mutate(country = as.character(country)) %>%
  group_by(country) %>%
  summarise(mj_per_cap_day = sum(per_capita_value / 365, na.rm = TRUE), .groups = "drop")

# Economic time (domestic + import), consumption-based kcal
tradeoff_kcal_econlabor_consump = summary_food_df_long %>%
  filter(type %in% c("hr_m", "hr_f"),
         footprint_type %in% c("domestic_per_capita", "import_per_capita")) %>%
  mutate(country = as.character(country)) %>%
  group_by(country, type) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  left_join(kcal_consumption, by = "country") %>%
  left_join(en_consumption, by = "country") %>%
  mutate(hr_per_2000kcal = hr_per_cap_day / kcal_per_cap_day * 2000,
         mj_per_2000kcal = mj_per_cap_day / kcal_per_cap_day * 2000) %>%
  left_join(regions %>% select(iso3c, continent), by = c("country" = "iso3c")) %>%
  drop_na() %>%
  mutate(is_row = as.character(country) %in% row_countries)

# Economic time (domestic + import), consumption-based protein
tradeoff_protein_econlabor_consump = summary_food_df_long %>%
  filter(type %in% c("hr_m", "hr_f"),
         footprint_type %in% c("domestic_per_capita", "import_per_capita")) %>%
  mutate(country = as.character(country)) %>%
  group_by(country, type) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  left_join(pro_consumption, by = "country") %>%
  left_join(en_consumption, by = "country") %>%
  mutate(hr_per_50g_protein = hr_per_cap_day / g_protein_per_cap_day * 50,
         mj_per_50g_protein = mj_per_cap_day / g_protein_per_cap_day * 50) %>%
  left_join(regions %>% select(iso3c, continent), by = c("country" = "iso3c")) %>%
  drop_na() %>%
  mutate(is_row = as.character(country) %in% row_countries)

# Total time (economic + household), consumption-based kcal — GHD countries only
# "_consump": Nutrients based on sum of domestic + import (excluding exports), i.e. consumption-based; time based on sum of domestic (all types) + import (economic only)
tradeoff_kcal_allwork_consump = summary_food_df_long_with_ghd %>%
  filter(type %in% c("hr_m", "hr_f"),
         country %in% cty_ghd,
         !footprint_type %in% c("export_per_capita")) %>%
  mutate(country = as.character(country)) %>%
  group_by(country, type) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  left_join(kcal_consumption, by = "country") %>%
  left_join(en_consumption, by = "country") %>%
  mutate(hr_per_2000kcal = hr_per_cap_day / kcal_per_cap_day * 2000,
         mj_per_2000kcal = mj_per_cap_day / kcal_per_cap_day * 2000) %>%
  left_join(regions %>% select(iso3c, continent), by = c("country" = "iso3c")) %>%
  drop_na() %>%
  mutate(is_row = as.character(country) %in% row_countries)

tradeoff_protein_allwork_consump = summary_food_df_long_with_ghd %>%
  filter(type %in% c("hr_m", "hr_f"),
         country %in% cty_ghd,
         !footprint_type %in% c("export_per_capita")) %>%
  mutate(country = as.character(country)) %>%
  group_by(country, type) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  left_join(pro_consumption, by = "country") %>%
  left_join(en_consumption, by = "country") %>%
  mutate(hr_per_50g_protein = hr_per_cap_day / g_protein_per_cap_day * 50,
         mj_per_50g_protein = mj_per_cap_day / g_protein_per_cap_day * 50) %>%
  left_join(regions %>% select(iso3c, continent), by = c("country" = "iso3c")) %>%
  drop_na() %>%
  mutate(is_row = as.character(country) %in% row_countries)

p_tradeoff_kcal_econlabor_consump = tradeoff_scatter(
  tradeoff_kcal_econlabor_consump %>% mutate(hr_per_2000kcal = hr_per_2000kcal * 60),
  "mj_per_2000kcal", "hr_per_2000kcal", "kcal_per_cap_day",
  "Energy (MJ / 2000 kcal)", "Time (min / 2000 kcal)", "kcal/cap/day",
  paste0("Energy vs. time to provision 2000 kcal (", year, ") — Economic, consumption-based"))

p_tradeoff_kcal_allwork_consump = tradeoff_scatter(
  tradeoff_kcal_allwork_consump, "mj_per_2000kcal", "hr_per_2000kcal", "kcal_per_cap_day",
  "Energy (MJ / 2000 kcal)", "Time (hr / 2000 kcal)", "kcal/cap/day",
  paste0("Energy vs. time to provision 2000 kcal (", year, ") — Total incl. household, consumption-based"))

p_tradeoff_protein_econlabor_consump = tradeoff_scatter(
  tradeoff_protein_econlabor_consump %>% mutate(hr_per_50g_protein = hr_per_50g_protein * 60),
  "mj_per_50g_protein", "hr_per_50g_protein", "g_protein_per_cap_day",
  "Energy (MJ / 50 g protein)", "Time (min / 50 g protein)", "g protein/cap/day",
  paste0("Energy vs. time to provision 50 g protein (", year, ") — Economic, consumption-based"))

p_tradeoff_protein_allwork_consump = tradeoff_scatter(
  tradeoff_protein_allwork_consump, "mj_per_50g_protein", "hr_per_50g_protein", "g_protein_per_cap_day",
  "Energy (MJ / 50 g protein)", "Time (hr / 50 g protein)", "g protein/cap/day",
  paste0("Energy vs. time to provision 50 g protein (", year, ") — Total incl. household, consumption-based"))

p_tradeoff_econlabor_consump = (p_tradeoff_kcal_econlabor_consump | p_tradeoff_protein_econlabor_consump) +
  plot_layout(guides = "collect") & theme(legend.position = "right")
p_tradeoff_allwork_consump = (p_tradeoff_kcal_allwork_consump | p_tradeoff_protein_allwork_consump) +
  plot_layout(guides = "collect") & theme(legend.position = "right")

ggsave(paste0("results/tradeoff_econlabor_consump.pdf"), p_tradeoff_econlabor_consump, width = 20, height = 8)
ggsave(paste0("results/tradeoff_allwork_consump.pdf"), p_tradeoff_allwork_consump, width = 20, height = 8)

# Per-capita scatter: energy (MJ/cap/day) vs time (hr/cap/day) — no nutrition normalisation

en_domestic = summary_food_df_long %>%
  filter(type == "en", footprint_type == "domestic_per_capita") %>%
  select(country, mj_per_cap_day = per_capita_value) %>%
  mutate(mj_per_cap_day = mj_per_cap_day / 365)

pro_domestic = summary_pro_df_long %>%
  filter(footprint_type == "domestic_per_capita") %>%
  select(country, g_protein_per_cap_day = per_capita_value)

pro_import = summary_pro_df_long %>%
  filter(footprint_type == "import_per_capita") %>%
  select(country, g_protein_per_cap_day = per_capita_value)

# Region-resolved counterparts, for joining with combined food+nonfood effort
# (which is region-keyed -- see chat discussion). pro_domestic/pro_import
# (187-country) stay as-is for food-only joins.
pro_domestic_region = summary_pro_region_df_long %>%
  filter(footprint_type == "domestic_per_capita") %>%
  select(exio_region, g_protein_per_cap_day = per_capita_value)

pro_import_region = summary_pro_region_df_long %>%
  filter(footprint_type == "import_per_capita") %>%
  select(exio_region, g_protein_per_cap_day = per_capita_value)

# Table: FABIO countries' household final protein consumption by source
# (domestic vs. import products, per FABIO_y_hh's own domestic/import split --
# see section 1.3 comment above for how this differs from "effort origin").
# df_nutri$protein's raw (non-per-capita) domestic/import columns are absolute
# country totals in grams protein/year (FABIO_y_hh_pro's native unit -- see
# 1.mrio_convert.R); converted here to g protein/capita/day using pop_data_yr.
protein_consumption_domestic_import = df_nutri[["protein"]] %>%
  transmute(country = as.character(country),
            domestic_protein_g = domestic,
            import_protein_g   = import) %>%
  left_join(pop_data_yr, by = c("country" = "iso3c")) %>%
  mutate(domestic_g_per_cap_day = domestic_protein_g / population / 365,
         import_g_per_cap_day   = import_protein_g / population / 365,
         total_g_per_cap_day    = domestic_g_per_cap_day + import_g_per_cap_day,
         import_share_pct       = import_g_per_cap_day / total_g_per_cap_day * 100) %>%
  filter(total_g_per_cap_day > 0) %>%
  left_join(fao_und_lookup, by = "country") %>%
  arrange(desc(total_g_per_cap_day)) %>%
  select(country, domestic_g_per_cap_day, import_g_per_cap_day, total_g_per_cap_day, import_share_pct,
         undernourishment_pct)

fwrite(protein_consumption_domestic_import, file = "output/protein_consumption_domestic_import_2020.csv")

gt_protein_consumption_domestic_import = protein_consumption_domestic_import %>%
  gt() %>%
  tab_header(title = paste0("Household protein consumption by source (", year, ")"),
             subtitle = "Domestic vs. imported products, FABIO countries") %>%
  cols_label(country = "Country",
             domestic_g_per_cap_day = "Domestic (g/cap/day)",
             import_g_per_cap_day = "Import (g/cap/day)",
             total_g_per_cap_day = "Total (g/cap/day)",
             import_share_pct = "Import share (%)",
             undernourishment_pct = "Undernourishment (%, 2018-2020)") %>%
  fmt_number(columns = c(domestic_g_per_cap_day, import_g_per_cap_day, total_g_per_cap_day), decimals = 1) %>%
  fmt_number(columns = c(import_share_pct, undernourishment_pct), decimals = 1)
print(gt_protein_consumption_domestic_import)

kcal_domestic = summary_kcal_df_long %>%
  filter(footprint_type == "domestic_per_capita") %>%
  mutate(country = as.character(country)) %>%
  select(country, kcal_per_cap_day = per_capita_value)

# Food only
tradeoff_pcap_econlabor = summary_food_df_long %>%
  filter(type %in% c("hr_m", "hr_f"), footprint_type == "domestic_per_capita") %>%
  select(country, type, hr_per_cap_day = per_capita_value) %>%
  inner_join(en_domestic, by = "country") %>%
  inner_join(pro_domestic, by = "country") %>%
  left_join(regions %>% select(iso3c, continent), by = c("country" = "iso3c")) %>%
  drop_na() %>%
  filter(!as.character(country) %in% row_countries)

tradeoff_pcap_allwork = summary_food_df_long_with_ghd %>%
  filter(type %in% c("hr_m", "hr_f"), country %in% cty_ghd,
         !footprint_type %in% c("export_per_capita", "import_per_capita")) %>%
  group_by(country, type) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  inner_join(en_domestic, by = "country") %>%
  inner_join(pro_domestic, by = "country") %>%
  left_join(regions %>% select(iso3c, continent), by = c("country" = "iso3c")) %>%
  drop_na() %>%
  mutate(is_row = as.character(country) %in% row_countries)

# Non-food sector only
# No need to consider non-economic time here since non-food sector time data is only economic (from EXIO); just sum across all domestic footprint types per country
en_domestic_nonfood = summary_nonfood_df_long %>%
  filter(type == "en", footprint_type == "domestic_per_capita") %>%
  mutate(country = as.character(country)) %>%
  group_by(country) %>%
  summarise(mj_per_cap_day = sum(per_capita_value, na.rm = TRUE) / 365, .groups = "drop")

tradeoff_pcap_nonfood_econlabor = summary_nonfood_df_long %>%
  filter(type %in% c("hr_m", "hr_f"), footprint_type == "domestic_per_capita") %>%
  mutate(country = as.character(country)) %>%
  group_by(country, type) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  inner_join(en_domestic_nonfood, by = "country") %>%
  inner_join(pro_domestic, by = "country") %>%
  left_join(regions %>% select(iso3c, continent), by = c("country" = "iso3c")) %>%
  drop_na() %>%
  filter(!as.character(country) %in% row_countries)

p_tradeoff_pcap_econlabor = tradeoff_scatter(
  tradeoff_pcap_econlabor %>% mutate(hr_per_cap_day = hr_per_cap_day * 60),
  "mj_per_cap_day", "hr_per_cap_day", "g_protein_per_cap_day",
  "Energy (MJ/cap/day)", "Time (min/cap/day)", "g protein/cap/day",
  paste0("Energy vs. time per capita (", year, ") — Food, Economic"))

p_tradeoff_pcap_allwork = tradeoff_scatter(
  tradeoff_pcap_allwork, "mj_per_cap_day", "hr_per_cap_day", "g_protein_per_cap_day",
  "Energy (MJ/cap/day)", "Time (hr/cap/day)", "g protein/cap/day",
  paste0("Energy vs. time per capita (", year, ") — Food, Economic + non-economic"))

p_tradeoff_pcap_nonfood_econlabor = tradeoff_scatter(
  tradeoff_pcap_nonfood_econlabor %>% mutate(hr_per_cap_day = hr_per_cap_day * 60),
  "mj_per_cap_day", "hr_per_cap_day", "g_protein_per_cap_day",
  "Energy (MJ/cap/day)", "Time (min/cap/day)", "g protein/cap/day",
  paste0("Energy vs. time per capita (", year, ") — Non-food, Economic"))

ggsave(paste0("results/tradeoff_pcap_econlabor update.pdf"),        p_tradeoff_pcap_econlabor,        width = 14, height = 6)
ggsave(paste0("results/tradeoff_pcap_allwork update.pdf"),        p_tradeoff_pcap_allwork,        width = 14, height = 6)
ggsave(paste0("results/tradeoff_pcap_nonfood_econlabor update.pdf"),  p_tradeoff_pcap_nonfood_econlabor,  width = 14, height = 6)

# Same as tradeoff_pcap_allwork.pdf, but hard-excluding RoW countries instead of just fading them
p_tradeoff_pcap_allwork_nonrow = tradeoff_scatter(
  tradeoff_pcap_allwork %>% filter(!is_row), "mj_per_cap_day", "hr_per_cap_day", "g_protein_per_cap_day",
  "Energy (MJ/cap/day)", "Time (hr/cap/day)", "g protein/cap/day",
  paste0("Energy vs. time per capita (", year, ") — Food, Economic + non-economic, non-RoW"))
ggsave("results/tradeoff_pcap_allwork_nonrow.pdf", p_tradeoff_pcap_allwork_nonrow, width = 14, height = 6)

# Protein energy-time tradeoff, side by side: left = total per-capita energy/time,
# right = per-50g-protein conversion factors, for non-RoW countries. Scope (which
# pcap/protein dataframe pair to draw from, e.g. econlabor vs. allwork) is a parameter
# so the same reused tradeoff_scatter() machinery can be pointed at either.
make_tradeoff_protein_side_by_side = function(pcap_df, protein_df, scope_label, filename, y_unit = c("hr", "min")) {
  y_unit = match.arg(y_unit)
  y_scale = if (y_unit == "min") 60 else 1

  # Hard-exclude RoW countries (rather than just fading them) to match "non-RoW countries"
  if ("is_row" %in% colnames(pcap_df))    pcap_df    = pcap_df    %>% filter(!is_row)
  if ("is_row" %in% colnames(protein_df)) protein_df = protein_df %>% filter(!is_row)

  pcap_df    = pcap_df    %>% mutate(hr_per_cap_day = hr_per_cap_day * y_scale)
  protein_df = protein_df %>% mutate(hr_per_50g_protein = hr_per_50g_protein * y_scale)

  p_left = tradeoff_scatter(
    pcap_df, "mj_per_cap_day", "hr_per_cap_day", "g_protein_per_cap_day",
    "Energy (MJ/cap/day)", paste0("Time (", y_unit, "/cap/day)"), "g protein/cap/day",
    paste0("Energy vs. time per capita (", year, ") — ", scope_label, ", non-RoW"))

  p_right = tradeoff_scatter(
    protein_df, "mj_per_50g_protein", "hr_per_50g_protein", "g_protein_per_cap_day",
    "Energy (MJ / 50 g protein)", paste0("Time (", y_unit, " / 50 g protein)"), "g protein/cap/day",
    paste0("Energy vs. time to provision 50 g protein (", year, ") — ", scope_label, ", non-RoW"))

  p_combined = (p_left | p_right) + plot_layout(guides = "collect") & theme(legend.position = "right")
  print(p_combined)
  ggsave(filename, p_combined, width = 20, height = 8)
  p_combined
}

# Food-sector, economic labor only (all non-RoW FABIO countries)
p_tradeoff_protein_side_by_side_econlabor = make_tradeoff_protein_side_by_side(
  tradeoff_pcap_econlabor, tradeoff_protein_econlabor, "Food, Economic",
  "results/tradeoff_protein_pcap_vs_50g_nonrow_econlabor.pdf", y_unit = "min")

# Food-sector, economic + non-economic (household) time — restricted to GHD-covered
# countries only, since household time isn't measured elsewhere; a much smaller set
# than the econlabor version above.
p_tradeoff_protein_side_by_side_allwork = make_tradeoff_protein_side_by_side(
  tradeoff_pcap_allwork, tradeoff_protein_allwork, "Food, Economic + non-economic",
  "results/tradeoff_protein_pcap_vs_50g_nonrow_allwork.pdf")

# Overlay econlabor and allwork on the same axes instead of separate panels. Energy
# (x) is identical between the two scopes (both draw domestic food-sector energy from
# the same source); only time (y) differs, since allwork adds household/non-economic
# time on top — so each country's two points are joined by a dotted VERTICAL line.
overlay_scatter = function(df, x_col, y_col, size_col, x_lab, y_lab, size_lab, title) {
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_line(aes(group = country), linetype = "dotted", color = "grey40", linewidth = 0.5) +
    geom_point(aes(size = .data[[size_col]], color = continent, shape = scope), alpha = 0.5) +
    ggrepel::geom_text_repel(data = df %>% filter(scope == "Economic + non-economic"),
                              aes(label = country), size = 3, max.overlaps = 15, show.legend = FALSE) +
    scale_size_continuous(range = c(1, 12)) +
    scale_shape_manual(values = c("Economic" = 16, "Economic + non-economic" = 17)) +
    guides(color = guide_legend(override.aes = list(size = 6))) +
    facet_wrap(~type, nrow = 1, labeller = labeller(type = c("hr_m" = "Male", "hr_f" = "Female"))) +
    labs(x = x_lab, y = y_lab, size = size_lab, color = "Continent", shape = "Scope", title = title) +
    theme_minimal() +
    theme(legend.position = "right",
          strip.text   = element_text(size = rel(1.5)),
          axis.title   = element_text(size = rel(1.4)),
          legend.text  = element_text(size = rel(1.3)),
          legend.title = element_text(size = rel(1.3)))
}


# Non-food: normalised by kcal and protein (non-food time/energy per unit nutrition)
tradeoff_kcal_nonfood_econlabor = tradeoff_pcap_nonfood_econlabor %>%
  inner_join(kcal_domestic, by = "country") %>%
  mutate(mj_per_2000kcal = mj_per_cap_day / kcal_per_cap_day * 2000,
         hr_per_2000kcal = hr_per_cap_day  / kcal_per_cap_day * 2000)

tradeoff_protein_nonfood_econlabor = tradeoff_pcap_nonfood_econlabor %>%
  mutate(mj_per_50g_protein = mj_per_cap_day / g_protein_per_cap_day * 50,
         hr_per_50g_protein = hr_per_cap_day  / g_protein_per_cap_day * 50)

p_tradeoff_kcal_nonfood_econlabor = tradeoff_scatter(
  tradeoff_kcal_nonfood_econlabor %>% mutate(hr_per_2000kcal = hr_per_2000kcal * 60),
  "mj_per_2000kcal", "hr_per_2000kcal", "kcal_per_cap_day",
  "Energy (MJ / 2000 kcal)", "Time (min / 2000 kcal)", "kcal/cap/day",
  paste0("Energy vs. time per 2000 kcal (", year, ") — Non-food, Economic"))

p_tradeoff_protein_nonfood_econlabor = tradeoff_scatter(
  tradeoff_protein_nonfood_econlabor %>% mutate(hr_per_50g_protein = hr_per_50g_protein * 60),
  "mj_per_50g_protein", "hr_per_50g_protein", "g_protein_per_cap_day",
  "Energy (MJ / 50 g protein)", "Time (min / 50 g protein)", "g protein/cap/day",
  paste0("Energy vs. time per 50 g protein (", year, ") — Non-food, Economic"))

p_tradeoff_nonfood_econlabor = (p_tradeoff_kcal_nonfood_econlabor | p_tradeoff_protein_nonfood_econlabor) +
  plot_layout(guides = "collect") & theme(legend.position = "right")

ggsave(paste0("results/tradeoff_convfac_nonfood_econlabor update.pdf"), p_tradeoff_nonfood_econlabor, width = 20, height = 8)

# Domestic tradeoff: total economic labor (food + non-food sectors combined),
# gender-collapsed — each country is a vertical segment at its (gender-invariant)
# total economic energy value, running from male hours to female hours of paid
# labor (food-sector + non-food-sector combined).
# Combines food+nonfood -> region-resolved (see chat discussion). country ->
# exio_region throughout; label uses exio_region (region names for RoW
# aggregates, otherwise identical in substance to the country it represents).
en_domestic_totalecon = bind_rows(summary_food_region_df_long, summary_nonfood_df_long) %>%
  filter(type == "en", footprint_type == "domestic_per_capita") %>%
  mutate(exio_region = as.character(exio_region)) %>%
  group_by(exio_region) %>%
  summarise(mj_per_cap_day = sum(per_capita_value, na.rm = TRUE) / 365, .groups = "drop")

tradeoff_pcap_totalecon = bind_rows(summary_food_region_df_long, summary_nonfood_df_long) %>%
  filter(type %in% c("hr_m", "hr_f"), footprint_type == "domestic_per_capita") %>%
  mutate(exio_region = as.character(exio_region)) %>%
  group_by(exio_region, type) %>%
  summarise(hr_per_cap_day = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  inner_join(en_domestic_totalecon, by = "exio_region") %>%
  inner_join(pro_domestic_region, by = "exio_region") %>%
  drop_na()

tradeoff_protein_totalecon = tradeoff_pcap_totalecon %>%
  mutate(mj_per_50g_protein = mj_per_cap_day / g_protein_per_cap_day * 50,
         hr_per_50g_protein = hr_per_cap_day  / g_protein_per_cap_day * 50)

label_tradeoff_protein_totalecon_gender = tradeoff_protein_totalecon %>%
  group_by(exio_region) %>%
  slice_max(hr_per_50g_protein, n = 1, with_ties = FALSE)

p_tradeoff_protein_totalecon_gender_range = ggplot(
  tradeoff_protein_totalecon,
  aes(x = mj_per_50g_protein, y = hr_per_50g_protein)) +
  geom_line(aes(group = exio_region), color = "grey60", linewidth = 0.6) +
  geom_point(aes(color = type, size = g_protein_per_cap_day), alpha = 0.85) +
  ggrepel::geom_text_repel(
    data = label_tradeoff_protein_totalecon_gender,
    aes(label = exio_region), size = 3.5, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c(hr_f = "#ca2323", hr_m = "#1f77b4"),
                      labels = c(hr_f = "Female", hr_m = "Male")) +
  scale_size_continuous(range = c(1, 8)) +
  labs(x = "Energy (MJ / 50 g protein)", y = "Time (hr / 50 g protein)",
       color = "Gender", size = "g protein/cap/day",
       title = paste0("Energy vs. time to provision 50 g protein (", year, ") — Total econ. labor (food + non-food), by gender")) +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title   = element_text(size = rel(1.4)),
        legend.text  = element_text(size = rel(1.1)),
        legend.title = element_text(size = rel(1.1)))

ggsave("results/tradeoff_convfac_protein_totalecon_gender_range update.pdf",
       p_tradeoff_protein_totalecon_gender_range, width = 14, height = 9)

# Same total economic (food + non-food) conversion factors, but for 50 g of
# *domestically-purchased/consumed* protein specifically (protein_source ==
# "domestic" in effort_consumption_df — the y_hh-based domestic/import
# protein-consumption split, matching pro_domestic), split into two panels by
# where the underlying energy/labor effort itself occurred: "Domestic effort"
# (produced and consumed in the same country) vs. "Import effort" (e.g. the
# labor/energy embodied in tuna caught abroad, for a domestically-canned can
# of tuna). See section 1.3 above for how effort_consumption_df is built.
# Combines food+nonfood -> region-resolved throughout (see chat discussion).
effort_totalecon_domestic = effort_consumption_df %>%
  filter(protein_source == "domestic") %>%
  mutate(exio_region = as.character(exio_region)) %>%
  group_by(exio_region, type, effort_origin, is_row) %>%
  summarise(per_capita_value = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

en_effort_totalecon_domestic = effort_totalecon_domestic %>%
  filter(type == "en") %>%
  select(exio_region, effort_origin, mj_per_cap_day = per_capita_value) %>%
  mutate(mj_per_cap_day = mj_per_cap_day / 365)

tradeoff_pcap_totalecon_domestic_effort = effort_totalecon_domestic %>%
  filter(type %in% c("hr_m", "hr_f")) %>%
  select(exio_region, type, effort_origin, hr_per_cap_day = per_capita_value) %>%
  left_join(en_effort_totalecon_domestic, by = c("exio_region", "effort_origin")) %>%
  inner_join(pro_domestic_region, by = "exio_region") %>%
  drop_na()

tradeoff_protein_totalecon_domestic_effort = tradeoff_pcap_totalecon_domestic_effort %>%
  mutate(mj_per_50g_protein = mj_per_cap_day / g_protein_per_cap_day * 50,
         hr_per_50g_protein = hr_per_cap_day  / g_protein_per_cap_day * 50,
         effort_label = factor(ifelse(effort_origin == "domestic", "Domestic effort", "Import effort"),
                                levels = c("Domestic effort", "Import effort")))

label_tradeoff_protein_totalecon_domestic_effort = tradeoff_protein_totalecon_domestic_effort %>%
  group_by(exio_region, effort_label) %>%
  slice_max(hr_per_50g_protein, n = 1, with_ties = FALSE)

p_tradeoff_protein_totalecon_domestic_effort = ggplot(
  tradeoff_protein_totalecon_domestic_effort,
  aes(x = mj_per_50g_protein, y = hr_per_50g_protein)) +
  geom_line(aes(group = exio_region), color = "grey60", linewidth = 0.6) +
  geom_point(aes(color = type, size = g_protein_per_cap_day), alpha = 0.85) +
  ggrepel::geom_text_repel(
    data = label_tradeoff_protein_totalecon_domestic_effort,
    aes(label = exio_region), size = 3.5, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c(hr_f = "#ca2323", hr_m = "#1f77b4"),
                      labels = c(hr_f = "Female", hr_m = "Male")) +
  scale_size_continuous(range = c(1, 8)) +
  # Extra expansion (esp. on the upper end) so the largest points -- up to 8mm
  # radius via scale_size_continuous above -- don't get clipped by the panel
  # border when they sit at/near the axis max (default 5% expansion is sized
  # for point *position*, not point *radius*).
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  facet_wrap(~effort_label, nrow = 1, scales = "free") +
  labs(x = "Energy (MJ / 50 g protein)", y = "Time (hr / 50 g protein)",
       color = "Gender", size = "g protein/cap/day",
       title = paste0("Energy vs. time per 50 g of domestically-consumed protein (", year,
                       ") — Total econ. labor (food + non-food), by effort origin")) +
  theme_minimal() +
  theme(legend.position = "right",
        strip.text   = element_text(size = rel(1.3)),
        axis.title   = element_text(size = rel(1.4)),
        legend.text  = element_text(size = rel(1.1)),
        legend.title = element_text(size = rel(1.1)))

ggsave("results/tradeoff_convfac_protein_totalecon_domestic_effort_gender_range update.pdf",
       p_tradeoff_protein_totalecon_domestic_effort, width = 20, height = 9)

# Same as above, but food sector only (no non-food sector labor/energy). Also
# region-resolved -- it's sliced from the same effort_consumption_df, which
# is region-keyed throughout (both food and nonfood domcons/impcons went
# through region_summary()); a separate country-resolved version isn't kept
# since nothing else needs it. Flag if you'd rather this one panel show full
# 187-country precision -- would need a parallel country-level effort table.
effort_foodecon_domestic = effort_consumption_df %>%
  filter(protein_source == "domestic", sector == "food") %>%
  mutate(exio_region = as.character(exio_region)) %>%
  group_by(exio_region, type, effort_origin, is_row) %>%
  summarise(per_capita_value = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

en_effort_foodecon_domestic = effort_foodecon_domestic %>%
  filter(type == "en") %>%
  select(exio_region, effort_origin, mj_per_cap_day = per_capita_value) %>%
  mutate(mj_per_cap_day = mj_per_cap_day / 365)

tradeoff_pcap_foodecon_domestic_effort = effort_foodecon_domestic %>%
  filter(type %in% c("hr_m", "hr_f")) %>%
  select(exio_region, type, effort_origin, hr_per_cap_day = per_capita_value) %>%
  left_join(en_effort_foodecon_domestic, by = c("exio_region", "effort_origin")) %>%
  inner_join(pro_domestic_region, by = "exio_region") %>%
  drop_na()

tradeoff_protein_foodecon_domestic_effort = tradeoff_pcap_foodecon_domestic_effort %>%
  mutate(mj_per_50g_protein = mj_per_cap_day / g_protein_per_cap_day * 50,
         hr_per_50g_protein = hr_per_cap_day  / g_protein_per_cap_day * 50,
         effort_label = factor(ifelse(effort_origin == "domestic", "Domestic effort", "Import effort"),
                                levels = c("Domestic effort", "Import effort")))

label_tradeoff_protein_foodecon_domestic_effort = tradeoff_protein_foodecon_domestic_effort %>%
  group_by(exio_region, effort_label) %>%
  slice_max(hr_per_50g_protein, n = 1, with_ties = FALSE)

p_tradeoff_protein_foodecon_domestic_effort = ggplot(
  tradeoff_protein_foodecon_domestic_effort,
  aes(x = mj_per_50g_protein, y = hr_per_50g_protein)) +
  geom_line(aes(group = exio_region), color = "grey60", linewidth = 0.6) +
  geom_point(aes(color = type, size = g_protein_per_cap_day), alpha = 0.85) +
  ggrepel::geom_text_repel(
    data = label_tradeoff_protein_foodecon_domestic_effort,
    aes(label = exio_region), size = 3.5, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c(hr_f = "#ca2323", hr_m = "#1f77b4"),
                      labels = c(hr_f = "Female", hr_m = "Male")) +
  scale_size_continuous(range = c(1, 8)) +
  # Extra expansion (esp. on the upper end) so the largest points -- up to 8mm
  # radius via scale_size_continuous above -- don't get clipped by the panel
  # border when they sit at/near the axis max (default 5% expansion is sized
  # for point *position*, not point *radius*).
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  facet_wrap(~effort_label, nrow = 1, scales = "free") +
  labs(x = "Energy (MJ / 50 g protein)", y = "Time (hr / 50 g protein)",
       color = "Gender", size = "g protein/cap/day",
       title = paste0("Energy vs. time per 50 g of domestically-consumed protein (", year,
                       ") — Food-sector econ. labor, by effort origin")) +
  theme_minimal() +
  theme(legend.position = "right",
        strip.text   = element_text(size = rel(1.3)),
        axis.title   = element_text(size = rel(1.4)),
        legend.text  = element_text(size = rel(1.1)),
        legend.title = element_text(size = rel(1.1)))

ggsave("results/tradeoff_convfac_protein_foodecon_domestic_effort_gender_range update.pdf",
       p_tradeoff_protein_foodecon_domestic_effort, width = 20, height = 9)

# Same as tradeoff_..._totalecon_domestic_effort above, but for 50 g of
# *imported* protein (protein_source == "import" in effort_consumption_df,
# matching pro_import), still split into "Domestic effort" / "Import effort"
# panels by where the underlying energy/labor itself occurred -- e.g. a can
# of tuna imported whole (import effort) vs. raw fish imported and canned
# domestically (domestic effort), both counted as imported protein here.
# Combines food+nonfood -> region-resolved throughout (see chat discussion).
effort_totalecon_import = effort_consumption_df %>%
  filter(protein_source == "import") %>%
  mutate(exio_region = as.character(exio_region)) %>%
  group_by(exio_region, type, effort_origin, is_row) %>%
  summarise(per_capita_value = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

en_effort_totalecon_import = effort_totalecon_import %>%
  filter(type == "en") %>%
  select(exio_region, effort_origin, mj_per_cap_day = per_capita_value) %>%
  mutate(mj_per_cap_day = mj_per_cap_day / 365)

tradeoff_pcap_totalecon_import_effort = effort_totalecon_import %>%
  filter(type %in% c("hr_m", "hr_f")) %>%
  select(exio_region, type, effort_origin, hr_per_cap_day = per_capita_value) %>%
  left_join(en_effort_totalecon_import, by = c("exio_region", "effort_origin")) %>%
  inner_join(pro_import_region, by = "exio_region") %>%
  drop_na()

tradeoff_protein_totalecon_import_effort = tradeoff_pcap_totalecon_import_effort %>%
  mutate(mj_per_50g_protein = mj_per_cap_day / g_protein_per_cap_day * 50,
         hr_per_50g_protein = hr_per_cap_day  / g_protein_per_cap_day * 50,
         effort_label = factor(ifelse(effort_origin == "domestic", "Domestic effort", "Import effort"),
                                levels = c("Domestic effort", "Import effort")))

label_tradeoff_protein_totalecon_import_effort = tradeoff_protein_totalecon_import_effort %>%
  group_by(exio_region, effort_label) %>%
  slice_max(hr_per_50g_protein, n = 1, with_ties = FALSE)

p_tradeoff_protein_totalecon_import_effort = ggplot(
  tradeoff_protein_totalecon_import_effort,
  aes(x = mj_per_50g_protein, y = hr_per_50g_protein)) +
  geom_line(aes(group = exio_region), color = "grey60", linewidth = 0.6) +
  geom_point(aes(color = type, size = g_protein_per_cap_day), alpha = 0.85) +
  ggrepel::geom_text_repel(
    data = label_tradeoff_protein_totalecon_import_effort,
    aes(label = exio_region), size = 3.5, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c(hr_f = "#ca2323", hr_m = "#1f77b4"),
                      labels = c(hr_f = "Female", hr_m = "Male")) +
  scale_size_continuous(range = c(1, 8)) +
  # Extra expansion (esp. on the upper end) so the largest points -- up to 8mm
  # radius via scale_size_continuous above -- don't get clipped by the panel
  # border when they sit at/near the axis max (default 5% expansion is sized
  # for point *position*, not point *radius*).
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  facet_wrap(~effort_label, nrow = 1, scales = "free") +
  labs(x = "Energy (MJ / 50 g protein)", y = "Time (hr / 50 g protein)",
       color = "Gender", size = "g protein/cap/day",
       title = paste0("Energy vs. time per 50 g of imported-consumed protein (", year,
                       ") — Total econ. labor (food + non-food), by effort origin")) +
  theme_minimal() +
  theme(legend.position = "right",
        strip.text   = element_text(size = rel(1.3)),
        axis.title   = element_text(size = rel(1.4)),
        legend.text  = element_text(size = rel(1.1)),
        legend.title = element_text(size = rel(1.1)))

ggsave("results/tradeoff_convfac_protein_totalecon_import_effort_gender_range update.pdf",
       p_tradeoff_protein_totalecon_import_effort, width = 20, height = 9)

# Same as above, but food sector only (no non-food sector labor/energy).
# Region-resolved for the same reason as effort_foodecon_domestic above.
effort_foodecon_import = effort_consumption_df %>%
  filter(protein_source == "import", sector == "food") %>%
  mutate(exio_region = as.character(exio_region)) %>%
  group_by(exio_region, type, effort_origin, is_row) %>%
  summarise(per_capita_value = sum(per_capita_value, na.rm = TRUE), .groups = "drop")

en_effort_foodecon_import = effort_foodecon_import %>%
  filter(type == "en") %>%
  select(exio_region, effort_origin, mj_per_cap_day = per_capita_value) %>%
  mutate(mj_per_cap_day = mj_per_cap_day / 365)

tradeoff_pcap_foodecon_import_effort = effort_foodecon_import %>%
  filter(type %in% c("hr_m", "hr_f")) %>%
  select(exio_region, type, effort_origin, hr_per_cap_day = per_capita_value) %>%
  left_join(en_effort_foodecon_import, by = c("exio_region", "effort_origin")) %>%
  inner_join(pro_import_region, by = "exio_region") %>%
  drop_na()

tradeoff_protein_foodecon_import_effort = tradeoff_pcap_foodecon_import_effort %>%
  mutate(mj_per_50g_protein = mj_per_cap_day / g_protein_per_cap_day * 50,
         hr_per_50g_protein = hr_per_cap_day  / g_protein_per_cap_day * 50,
         effort_label = factor(ifelse(effort_origin == "domestic", "Domestic effort", "Import effort"),
                                levels = c("Domestic effort", "Import effort")))

label_tradeoff_protein_foodecon_import_effort = tradeoff_protein_foodecon_import_effort %>%
  group_by(exio_region, effort_label) %>%
  slice_max(hr_per_50g_protein, n = 1, with_ties = FALSE)

p_tradeoff_protein_foodecon_import_effort = ggplot(
  tradeoff_protein_foodecon_import_effort,
  aes(x = mj_per_50g_protein, y = hr_per_50g_protein)) +
  geom_line(aes(group = exio_region), color = "grey60", linewidth = 0.6) +
  geom_point(aes(color = type, size = g_protein_per_cap_day), alpha = 0.85) +
  ggrepel::geom_text_repel(
    data = label_tradeoff_protein_foodecon_import_effort,
    aes(label = exio_region), size = 3.5, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c(hr_f = "#ca2323", hr_m = "#1f77b4"),
                      labels = c(hr_f = "Female", hr_m = "Male")) +
  scale_size_continuous(range = c(1, 8)) +
  # Extra expansion (esp. on the upper end) so the largest points -- up to 8mm
  # radius via scale_size_continuous above -- don't get clipped by the panel
  # border when they sit at/near the axis max (default 5% expansion is sized
  # for point *position*, not point *radius*).
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  facet_wrap(~effort_label, nrow = 1, scales = "free") +
  labs(x = "Energy (MJ / 50 g protein)", y = "Time (hr / 50 g protein)",
       color = "Gender", size = "g protein/cap/day",
       title = paste0("Energy vs. time per 50 g of imported-consumed protein (", year,
                       ") — Food-sector econ. labor, by effort origin")) +
  theme_minimal() +
  theme(legend.position = "right",
        strip.text   = element_text(size = rel(1.3)),
        axis.title   = element_text(size = rel(1.4)),
        legend.text  = element_text(size = rel(1.1)),
        legend.title = element_text(size = rel(1.1)))

ggsave("results/tradeoff_convfac_protein_foodecon_import_effort_gender_range update.pdf",
       p_tradeoff_protein_foodecon_import_effort, width = 20, height = 9)


# library(ggplot2)
# # Labor hours
# ggplot(summary_food_df %>% filter(type %in% c("hr_m", "hr_f")) %>%
#          # Order countries by sum of domestic_per_capita of type "hr_m" and "hr_f"
#           mutate(country = factor(country, levels = sum_ord)),
#        aes(x=country, y=domestic_per_capita, fill=type)) +
#   geom_bar(stat="identity", position="stack") +
#   labs(x="Country (ISO3)", y="Daily time footprint per capita (hr/cap/day)", alpha="gender",
#        title=paste0("Food-related time footprint per capita by country (", year, ")")) +
#   theme_minimal() +
#   theme(legend.position = "top") +
#   # Add alpha values for genders
#   scale_fill_manual(values=c("hr_m"="#1f77b4", "hr_f"="#2ca02c")) +
#   # Tilt x-axis labels
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#
#   # Add import_per_capita values as negative y axis bars by gender
#   geom_bar(data = summary_food_df %>% filter(type %in% c("hr_m", "hr_f")),
#            aes(x=country, y=-import_per_capita, fill=type), stat="identity", position="stack") +
#   scale_fill_manual(values=c("hr_m"="#1f77b4", "hr_f"="#2ca02c")) +
#
#   # Stack export_per_capita values as positive y axis bars by gender
#   geom_bar(data = summary_food_df %>% filter(type %in% c("hr_m", "hr_f")),
#            aes(x=country, y=export_per_capita, fill=type), stat="identity", position="stack") +
#   scale_fill_manual(values=c("hr_m"="#1f77b4", "hr_f"="#2ca02c")) +
#
#   labs(fill="hours by gender", title=paste0("Food-related time footprint per capita by country (", year, ")\n(positive: domestic+export, negative: import)"))



# Aggregate indirect (non-food sector) intensities to FABIO-product level
l_int_i_agg <- lapply(l_int_i, colSums)

# Country-wise consumption-based footprint

consumption = "food"

for (country in regions$iso3c) {

  Y_country <- FABIO_y[, which(fd$iso3c == country)]
  colnames(Y_country) <- fd$fd[fd$iso3c == country]

  pop = subset(countrypops, country_code_3 == country & year == yr)$population
  if (country %in% setdiff(regions$iso3c, unique(countrypops$country_code_3))) {pop=NA}

  for (extension in names(l_int_d)) {
    for (sector in c("food", "non_food")) {

      int <- if (sector == "food") l_int_d[[extension]] else l_int_i_agg[[extension]]

      print(paste("Calculating", extension, sector, "footprint for", country))
    
      # if(spread_stocks){
      #   stock_ratio <- Y_country[, "stock_addition"] / (rowSums(Y_country) - Y_country[, "stock_addition"])
      #   stock_ratio[!is.finite(stock_ratio)] <- 0
      #   Y_country <- as.data.table(as.matrix(Y_country))
      #   Y_country[, `:=`(food = food * (1 + stock_ratio),
      #                    other = other * (1 + stock_ratio),
      #                    tourist = tourist * (1 + stock_ratio),
      #                    unspecified = unspecified * (1 + stock_ratio),
      #                    stock_addition = 0)]
      # }
      
      MP <- int * FABIO_L
      
      # Initialize empty matrix to store results (row: exporter, col: importer)
      mat <- matrix(0, nrow=nrreg, ncol = nrreg, byrow=TRUE)
      rownames(mat) <- colnames(mat) <- sort(regions$iso3c)
      
      # Calculate footprints (energy & labor hr)
      country_consump = as.vector(as.matrix(Y_country[, consumption]))
      FP <- t(t(MP) * country_consump) # <= dim (23001x23001)
      colnames(FP) <- rownames(FP) <- paste0(io$iso3c, "_", io$item)
      FP <- as(FP, "TsparseMatrix")
      
      # Calculate calorie & protein production/consumption 
      x_country = t(t(FABIO_L) * country_consump) # mass flows
      
      # FP_cal = sweep(x_country, 2, 1000*coeff_cal, "*") 
      FP_cal = t(t(x_country) * as.vector(coeff_cal$kcal_per_kg)) * 1000 # kcal flows
      FP_pro = t(t(x_country) * as.vector(coeff_pro$protein_per_kg)) * 1000 # g-protein flows
      colnames(FP_cal) <- rownames(FP_cal) <- colnames(FP_pro) <- rownames(FP_pro) <- paste0(io$iso3c, "_", io$item)
      FP_cal <- as(FP_cal, "TsparseMatrix")
      FP_pro <- as(FP_pro, "TsparseMatrix")
      
      cal_consum = (country_consump*coeff_cal$kcal_per_kg)/365/pop*1000 #kcal/cap/day
      pro_consum = (country_consump*coeff_pro$protein_per_kg)/365/pop*1000 #kcal/cap/day

      # Make by-country (for both origin and target) for calorie and protein flows by partial summing each nrreg entries of cal_consum and pro_consum 
    
      cal_prod = (rowSums(x_country)*coeff_cal$kcal_per_kg)/365/pop*1000 #kcal/cap/day
      pro_prod = (rowSums(x_country)*coeff_pro$protein_per_kg)/365/pop*1000 #kcal/cap/day
      
      results <- data.table(origin=rownames(FP)[FP@i + 1], target=colnames(FP)[FP@j + 1], 
                            value = FP@x, # M.hr
                            value_pcap = FP@x / pop * 1e6, # hr/cap (for labor)
                            cal = FP_cal@x,
                            cal_pcap = FP_cal@x / pop, # kcal per capita
                            pro = FP_pro@x,
                            pro_pcap = FP_pro@x / pop  # g-protein per capita
      ) 
      results[,`:=`(country_consumer = country,
                    year = year,
                    indicator = extension,
                    sector = sector,
                    country_origin = substr(origin,1,3),
                    item_origin = substr(origin,5,100),
                    country_target = substr(target,1,3),
                    item_target = substr(target,5,100))]
      
      results[,`:=`(group_origin = items$comm_group[match(results$item_origin,items$item)],
                    group_target = items$comm_group[match(results$item_target,items$item)],
                    continent_origin = regions$continent[match(results$country_origin, regions$iso3c)])]
      
      results$continent_origin[results$country_origin==country] <- country
      
      # Remove "Live animals" from the nutrient flows
      # CHECK: Need to be in energy/labor footprints?
      results = results %>%
        filter(group_origin != "Live animals")
      
      # print(paste(country, "ratio of FD vs. production:", 
      #             round(sum(as.matrix(Y_country)[, consumption]) / sum(results$value),4)))
      
      data_tot <- results %>%
        group_by(item_target, country_origin) %>%
        filter(value != 0) %>%
        summarise(value = (sum(value))) %>%
        spread(country_origin, value, fill = 0) %>% # Al
        # Add a row with column sums
        ungroup() %>%
        bind_rows(summarise(., item_target = "Total", across(-item_target, sum)))
      
      data.table::fwrite(data_tot, file=paste0("output/FABIO_", country,"_", year, "_", extension, "_", sector, "_", consumption,".csv"), sep=",")

      data_imp_tot = tail(data_tot, 1) %>% rename(importer = item_target) %>% mutate(importer = country)

      # Fill mat with data_imp_tot where rownames(mat) match names(data_imp_tot), and colnames(mat) match data_imp_tot$importer
      mat[rownames(mat) %in% names(data_imp_tot), data_imp_tot$importer[1]] <-
        as.numeric(data_imp_tot[1, names(data_imp_tot) %in% rownames(mat)])
    }
  }
}

# Total Mcal trade between countries
saveRDS(mat, file.path("data/calorie_trade_mat.rds"))

mat = readRDS(file.path("data/calorie_trade_mat.rds"))


# TODO
# - Consider per-capita results
# - Consider which metric to present
# - Consider how to handle energy & hour vs. calorie & protein





# ====

# Format the output matrices
Y_cons = data.table(as.matrix(FABIO_y_hh_cal))
# rownames(Y_cons) <- paste0(io$iso3c, "_", io$item)
colnames(Y_cons) <- regions$iso3c
Y_cons[,`:=`(iso3c = io$iso3c,
             item_origin = io$item)]

Y_sq <- Y_cons %>%
  group_by(iso3c) %>%
  # Summarise total consumption by country
  dplyr::summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
  column_to_rownames(var = "iso3c") %>%
  # order columns alphabetically
  select(sort(peek_vars()))

# Total kcal trade from Y
# mat_y is basically the same as agg_country_footprint(FABIO_y_hh_cal).
mat_y = as.matrix(Y_sq) 
# Save total consumption kcal trade between countries
saveRDS(mat_y, file.path("data/calorie_trade_mat_cons.rds")) # in kcal


mat_y = readRDS(file.path("data/calorie_trade_mat_cons.rds")) # in kcal

# Summary matrix by country (by row)
mat_cons = get_mat_summary(mat_y)

mat_cons_net = get_net_trade(mat_y)
mat_cons_net[mat_cons_net < 3e9] <- 0

# # Circular plot(too heavy)
# library(circlize)
# chordDiagram(mat_cons_net, directional = 1, direction.type = c("arrows", "diffHeight"), )
# title(main = paste0("Net (consumption) calorie trade ", year, " (Mcal)"))


# Look into individual countries

cty = "USA"
imp = round(mat_y[,cty] / pop_y$pop[pop_y$iso3c==cty] / 365, 4)
ex = round(mat_y[cty,] / pop_y$pop[pop_y$iso3c==cty] / 365, 4)

# Domestic consumption
print(paste("Total consumption in", cty, ":", sum(imp, na.rm=TRUE), "kcal/capita/day"))
# Tot imports
print(paste("Total imports in", cty, ":", sum(imp[names(imp)!=cty], na.rm=TRUE), "kcal/capita/day"))
# Tot exports
print(paste("Total exports from", cty, ":", sum(ex[names(ex)!=cty], na.rm=TRUE), "kcal/capita/day"))



#### Sankey diagrams: kcal and protein trade by continent ####

library(networkD3)

# Build protein trade matrix (same structure as mat_y for kcal)
Y_pro_cons = data.table(as.matrix(FABIO_y_hh_pro))
colnames(Y_pro_cons) <- regions$iso3c
Y_pro_cons[, iso3c := io$iso3c]

mat_y_pro = Y_pro_cons %>%
  group_by(iso3c) %>%
  dplyr::summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
  column_to_rownames(var = "iso3c") %>%
  select(sort(peek_vars())) %>%
  as.matrix()  # g protein; rows = producer, cols = consumer country

# Aggregate country×country matrix to continent×continent (not used)
agg_to_continent <- function(mat) {
  cty  = colnames(mat)
  cont = regions$continent[match(cty, regions$iso3c)]
  conts = sort(unique(na.omit(cont)))
  mat_c = matrix(0, nrow = length(conts), ncol = length(conts),
                 dimnames = list(conts, conts))
  for (i in conts)
    for (j in conts)
      mat_c[i, j] = sum(mat[which(cont == i), which(cont == j)], na.rm = TRUE)
  mat_c
}

# Aggregate country×country matrix to country level, collapsing small countries into
# "Other {continent}" buckets. Diagonal (domestic flows) is zeroed out.
agg_to_country_sankey <- function(mat, n_top = 20) {
  cty   <- rownames(mat)
  cont  <- regions$continent[match(cty, regions$iso3c)]

  # Rank by total off-diagonal trade volume
  off   <- mat; diag(off) <- 0
  total <- rowSums(off) + colSums(off)
  top   <- names(sort(total, decreasing = TRUE))[seq_len(min(n_top, length(cty)))]

  labels <- ifelse(cty %in% top, cty,
                   paste0("Other ", ifelse(is.na(cont), "World", cont)))

  mat_agg        <- t(rowsum(t(rowsum(mat, group = labels)), group = labels))
  diag(mat_agg)  <- 0  # remove domestic / intra-group flows

  pop_lookup <- data.frame(iso3c = cty, label = labels) %>%
    left_join(pop_y, by = "iso3c") %>%
    group_by(label) %>%
    summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop")

  list(mat = mat_agg, pop = pop_lookup)
}

# Region-resolved counterpart, for non-food matrices (l_nonfood_region):
# nodes are EXIO regions (a RoW aggregate is one node, not split across its
# member countries) rather than FABIO countries. mat must already be a
# square 49x49 region-region matrix (agg_to_region_matrix(l_nonfood_region[[metric]])).
agg_to_region_sankey <- function(mat, n_top = 20) {
  reg  <- rownames(mat)
  cont <- region_continent_of_index[match(reg, region_name_of_index)]

  off   <- mat; diag(off) <- 0
  total <- rowSums(off) + colSums(off)
  top   <- names(sort(total, decreasing = TRUE))[seq_len(min(n_top, length(reg)))]

  labels <- ifelse(reg %in% top, reg,
                   paste0("Other ", ifelse(is.na(cont), "World", cont)))

  mat_agg        <- t(rowsum(t(rowsum(mat, group = labels)), group = labels))
  diag(mat_agg)  <- 0

  pop_lookup <- data.frame(exio_region = reg, label = labels, pop = region_population[match(reg, region_name_of_index)]) %>%
    group_by(label) %>%
    summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop")

  list(mat = mat_agg, pop = pop_lookup)
}

# Shared aggregation for two matrices, using the same country labels for both.
# Ensures prod and cons Sankey nodes are directly comparable.
agg_to_country_sankey_shared <- function(mat_prod, mat_cons, n_top = 20) {
  cty  <- rownames(mat_prod)
  cont <- regions$continent[match(cty, regions$iso3c)]

  # Rank by combined off-diagonal volume across both matrices
  off_prod <- mat_prod; diag(off_prod) <- 0
  off_cons <- mat_cons; diag(off_cons) <- 0
  total <- rowSums(off_prod) + colSums(off_prod) + rowSums(off_cons) + colSums(off_cons)
  top   <- names(sort(total, decreasing = TRUE))[seq_len(min(n_top, length(cty)))]

  labels <- ifelse(cty %in% top, cty,
                   paste0("Other ", ifelse(is.na(cont), "World", cont)))

  agg_mat <- function(m) {
    m_agg       <- t(rowsum(t(rowsum(m, group = labels)), group = labels))
    diag(m_agg) <- 0
    m_agg
  }

  pop_lookup <- data.frame(iso3c = cty, label = labels) %>%
    left_join(pop_y, by = "iso3c") %>%
    group_by(label) %>%
    summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop")

  list(mat_prod = agg_mat(mat_prod), mat_cons = agg_mat(mat_cons), pop = pop_lookup)
}

# Only the last stage of trade before final consumption is considered.
res_kcal_cty_cons = agg_to_country_sankey(mat_y);      
mat_kcal_cty_cons = res_kcal_cty_cons$mat; 
pop_kcal_cty_cons = res_kcal_cty_cons$pop
res_pro_cty_cons  = agg_to_country_sankey(mat_y_pro);  
mat_pro_cty_cons  = res_pro_cty_cons$mat;  
pop_pro_cty_cons  = res_pro_cty_cons$pop
# This reflects the flow of calories/protein from producer to consumer countries through the trade chain. This avoids double-counting of intermediate trade flows.
res_kcal_cty = agg_to_country_sankey(agg_country_footprint(FABIO_x_hh_cal)); 
mat_kcal_cty = res_kcal_cty$mat; 
pop_kcal_cty = res_kcal_cty$pop
res_pro_cty  = agg_to_country_sankey(agg_country_footprint(FABIO_x_hh_pro));  
mat_pro_cty  = res_pro_cty$mat;  
pop_pro_cty  = res_pro_cty$pop

# Helper: map a Sankey node label to its continent.
# Labels are either iso3c codes (top countries) or "Other {continent}".
label_to_continent <- function(labels) {
  other <- grepl("^Other ", labels)
  # Try ISO3 first (food Sankeys' node labels); fall back to EXIO region name
  # (non-food Sankeys' node labels -- e.g. "China", "RoW Africa" -- which
  # never match regions$iso3c).
  cont_iso    <- regions$continent[match(labels, regions$iso3c)]
  cont_region <- region_continent_of_index[match(labels, region_name_of_index)]
  cont <- ifelse(other, sub("^Other ", "", labels),
                 ifelse(!is.na(cont_iso), cont_iso, cont_region))
  ifelse(is.na(cont), "World", cont)
}

# Colour scale: one distinct colour per continent, shared across all Sankeys.
.cont_levs      <- sort(unique(c(na.omit(regions$continent), "World")))
.cont_pal       <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
sankey_cont_cols <- setNames(.cont_pal[seq_along(.cont_levs)], .cont_levs)
sankey_colour_scale <- htmlwidgets::JS(sprintf(
  'd3.scaleOrdinal().domain([%s]).range([%s])',
  paste0('"', .cont_levs, '"', collapse = ", "),
  paste0('"', .cont_pal[seq_along(.cont_levs)], '"', collapse = ", ")
))
rm(.cont_levs, .cont_pal)

# Convert continent×continent matrix to Sankey links/nodes
# Producer nodes sit on the left, consumer nodes on the right (avoids self-loop issue)
mat_to_sankey <- function(mat, scale, pop = NULL, pcap_label = NULL, digits = 1, pcap_scale = 1) {
  conts = rownames(mat)
  prod_names = conts[order(rowSums(mat), decreasing = TRUE)]
  cons_names = colnames(mat)[order(colSums(mat), decreasing = TRUE)]
  nodes = data.frame(
    name  = c(paste0(prod_names, " (prod)"), paste0(cons_names, " (cons)")),
    group = label_to_continent(c(prod_names, cons_names))
  )

  if (!is.null(pop) && !is.null(pcap_label)) {
    prod_pop = pop$pop[match(prod_names, pop$label)]
    cons_pop = pop$pop[match(cons_names, pop$label)]
    pcap = c(rowSums(mat)[prod_names] / prod_pop / 365 * pcap_scale,
             colSums(mat)[cons_names] / cons_pop / 365 * pcap_scale)
    nodes$tooltip = paste0(formatC(pcap, format = "f", digits = digits, big.mark = ","), " ", pcap_label)
  }

  links = as.data.frame(mat) %>%
    tibble::rownames_to_column("source_name") %>%
    pivot_longer(-source_name, names_to = "target_name", values_to = "value") %>%
    filter(value > 0) %>%
    mutate(
      source = match(paste0(source_name, " (prod)"), nodes$name) - 1L,
      target = match(paste0(target_name, " (cons)"), nodes$name) - 1L,
      value  = value / scale
    )
  list(nodes = nodes, links = as.data.frame(links))
}

# Build a Sankey where left nodes are sized by production (x) and right nodes by
# consumption (y). Supply-chain loss per producer fills the gap, making flows
# conserved at every node.
#
# mat_prod[A,B]: production in A attributed to B's consumption (Leontief-based)
# mat_cons[A,B]: actual household calories in B that originated from A (final demand Y)
# Loss per producer A = rowSums(mat_prod)[A] - rowSums(mat_cons)[A]
mat_to_sankey_with_loss <- function(mat_prod, mat_cons, scale,
                                    pop = NULL, pcap_label = NULL) {
  producers <- rownames(mat_prod)
  consumers <- colnames(mat_cons)

  # Order nodes by descending total volume
  prod_ord <- producers[order(rowSums(mat_prod), decreasing = TRUE)]
  cons_ord <- consumers[order(colSums(mat_cons), decreasing = TRUE)]
  loss_ord <- prod_ord  # one loss node per producer, same order

  nodes <- data.frame(
    name  = c(paste0(prod_ord, " (prod)"),
              paste0(cons_ord, " (cons)"),
              paste0(loss_ord, " (loss)")),
    group = c(label_to_continent(prod_ord),
              label_to_continent(cons_ord),
              label_to_continent(loss_ord))
  )

  if (!is.null(pop) && !is.null(pcap_label)) {
    prod_pop  <- pop$pop[match(prod_ord, pop$label)]
    cons_pop  <- pop$pop[match(cons_ord, pop$label)]
    loss_vals <- rowSums(mat_prod - mat_cons)[loss_ord]
    loss_pct  <- loss_vals / rowSums(mat_prod)[loss_ord] * 100
    nodes$tooltip <- c(
      paste0(formatC(rowSums(mat_prod)[prod_ord] / prod_pop / 365,
                     format = "f", digits = 1), " ", pcap_label, " (prod)"),
      paste0(formatC(colSums(mat_cons)[cons_ord] / cons_pop / 365,
                     format = "f", digits = 1), " ", pcap_label, " (cons)"),
      paste0(formatC(loss_pct, format = "f", digits = 1), "% supply-chain loss")
    )
  }

  # Links: producer → consumer (actual consumption values)
  links_cons <- as.data.frame(mat_cons) %>%
    tibble::rownames_to_column("source_name") %>%
    pivot_longer(-source_name, names_to = "target_name", values_to = "value") %>%
    filter(value > 0) %>%
    mutate(
      source = match(paste0(source_name, " (prod)"), nodes$name) - 1L,
      target = match(paste0(target_name, " (cons)"), nodes$name) - 1L,
      value  = value / scale
    )

  # Links: producer → its own loss node (remainder not reaching final consumption)
  loss_vals <- rowSums(mat_prod - mat_cons)
  links_loss <- data.frame(source_name = names(loss_vals), value = loss_vals) %>%
    filter(value > 0) %>%
    mutate(
      source = match(paste0(source_name, " (prod)"), nodes$name) - 1L,
      target = match(paste0(source_name, " (loss)"), nodes$name) - 1L,
      value  = value / scale
    )

  links <- bind_rows(
    links_cons %>% select(source, target, value),
    links_loss %>% select(source, target, value)
  )

  list(nodes = nodes, links = as.data.frame(links))
}

sankey_kcal      = mat_to_sankey(mat_kcal_cty,      scale = 1e12, pop = pop_kcal_cty,      pcap_label = "kcal/cap/day")
sankey_pro       = mat_to_sankey(mat_pro_cty,       scale = 1e9,  pop = pop_pro_cty,       pcap_label = "g protein/cap/day")
sankey_kcal_cons = mat_to_sankey(mat_kcal_cty_cons, scale = 1e12, pop = pop_kcal_cty_cons, pcap_label = "kcal/cap/day")
sankey_pro_cons  = mat_to_sankey(mat_pro_cty_cons,  scale = 1e9,  pop = pop_pro_cty_cons,  pcap_label = "g protein/cap/day")

add_pcap_tooltip <- function(p, sankey) {
  # sankeyNetwork drops all columns except name/group, so re-inject tooltip directly
  p$x$nodes$tooltip <- sankey$nodes$tooltip
  htmlwidgets::onRender(p, "
    function(el, x) {
      d3.select(el).selectAll('.node').select('title').each(function(d) {
        if (d.tooltip) d3.select(this).text(d3.select(this).text() + '\\n' + d.tooltip);
      });
    }
  ")
}

add_continent_legend <- function(p) {
  domain_json <- paste0('["', paste(names(sankey_cont_cols), collapse = '","'), '"]')
  range_json  <- paste0('["', paste(unname(sankey_cont_cols), collapse = '","'), '"]')
  htmlwidgets::onRender(p, sprintf("
    function(el, x) {
      var domain = %s;
      var range  = %s;
      var legW   = 120;
      var svg    = d3.select(el).select('svg');
      var mainG  = svg.select('g');

      // Shift the diagram right to make room for the legend on the left
      var t = mainG.attr('transform') || 'translate(0,0)';
      var m = t.match(/translate\\(\\s*([^,\\s]+)[,\\s]+([^)\\s]+)\\s*\\)/);
      var tx = m ? +m[1] : 0, ty = m ? +m[2] : 0;
      mainG.attr('transform', 'translate(' + (tx + legW) + ',' + ty + ')');
      svg.attr('width', (+svg.attr('width') || el.offsetWidth) + legW);

      // Continent legend
      var leg = svg.append('g')
        .attr('class', 'cont-legend')
        .attr('transform', 'translate(8, 10)');
      domain.forEach(function(d, i) {
        var g = leg.append('g').attr('transform', 'translate(0,' + i * 20 + ')');
        g.append('rect')
          .attr('width', 13).attr('height', 13)
          .attr('fill', range[i]).attr('opacity', 0.85);
        g.append('text')
          .attr('x', 18).attr('y', 11)
          .style('font-size', '12px').style('font-family', 'sans-serif')
          .text(d);
      });

      // Node total labels — vertical text centred on each block
      var fmt   = d3.format('.3s');
      var units = x.options.units || '';
      mainG.selectAll('.node').each(function(d) {
        if (d.dy < 30) return;
        var cx = d.dx / 2, cy = d.dy / 2;
        d3.select(this).append('text')
          .attr('class', 'node-total')
          .attr('transform', 'rotate(-90,' + cx + ',' + cy + ')')
          .attr('x', cx).attr('y', cy)
          .attr('dy', '0.35em')
          .attr('text-anchor', 'middle')
          .style('font-size', '9px')
          .style('font-family', 'sans-serif')
          .style('fill', '#333')
          .style('stroke', 'white')
          .style('stroke-width', '2.5px')
          .style('paint-order', 'stroke fill')
          .style('pointer-events', 'none')
          .text(fmt(d.value) + ' ' + units);
      });
    }
  ", domain_json, range_json))
}

p_sankey_kcal = sankeyNetwork(
  Links = sankey_kcal$links, Nodes = sankey_kcal$nodes,
  Source = "source", Target = "target", Value = "value", NodeID = "name",
  NodeGroup = "group", colourScale = sankey_colour_scale,
  sinksRight = FALSE, fontSize = 13, nodeWidth = 20, nodePadding = 10,
  units = "Tcal", iterations = 0
) %>% add_pcap_tooltip(sankey_kcal) %>% add_continent_legend()

p_sankey_pro = sankeyNetwork(
  Links = sankey_pro$links, Nodes = sankey_pro$nodes,
  Source = "source", Target = "target", Value = "value", NodeID = "name",
  NodeGroup = "group", colourScale = sankey_colour_scale,
  sinksRight = FALSE, fontSize = 13, nodeWidth = 20, nodePadding = 10,
  units = "kt protein", iterations = 0
) %>% add_pcap_tooltip(sankey_pro) %>% add_continent_legend()

p_sankey_kcal_cons = sankeyNetwork(
  Links = sankey_kcal_cons$links, Nodes = sankey_kcal_cons$nodes,
  Source = "source", Target = "target", Value = "value", NodeID = "name",
  NodeGroup = "group", colourScale = sankey_colour_scale,
  sinksRight = FALSE, fontSize = 13, nodeWidth = 20, nodePadding = 10,
  units = "Tcal", iterations = 0
) %>% add_pcap_tooltip(sankey_kcal_cons) %>% add_continent_legend()

p_sankey_pro_cons = sankeyNetwork(
  Links = sankey_pro_cons$links, Nodes = sankey_pro_cons$nodes,
  Source = "source", Target = "target", Value = "value", NodeID = "name",
  NodeGroup = "group", colourScale = sankey_colour_scale,
  sinksRight = FALSE, fontSize = 13, nodeWidth = 20, nodePadding = 10,
  units = "kt protein", iterations = 0
) %>% add_pcap_tooltip(sankey_pro_cons) %>% add_continent_legend()

# These are without  losses.
htmlwidgets::saveWidget(p_sankey_kcal, "results/sankey_kcal.html",     selfcontained = FALSE)
htmlwidgets::saveWidget(p_sankey_pro,  "results/sankey_protein.html",  selfcontained = FALSE)
htmlwidgets::saveWidget(p_sankey_kcal_cons, "results/sankey_kcal_cons.html",     selfcontained = FALSE)
htmlwidgets::saveWidget(p_sankey_pro_cons,  "results/sankey_protein_cons.html",  selfcontained = FALSE)

# Combined prod-x / cons-y Sankey with supply-chain loss nodes:
#   Left  node width = production calories attributed to food consumption (FABIO_x_hh_cal)
#   Right node width = actual household consumption calories (mat_y)
#   Loss  node       = the difference per producer (supply-chain + processing losses)
res_combined_kcal <- agg_to_country_sankey_shared(
  agg_country_footprint(FABIO_x_hh_cal), agg_country_footprint(FABIO_y_hh_cal))
res_combined_pro  <- agg_to_country_sankey_shared(
  agg_country_footprint(FABIO_x_hh_pro), agg_country_footprint(FABIO_y_hh_pro))

sankey_combined_kcal <- mat_to_sankey_with_loss(
  res_combined_kcal$mat_prod, res_combined_kcal$mat_cons,
  scale = 1e12, pop = res_combined_kcal$pop, pcap_label = "kcal/cap/day")
sankey_combined_pro  <- mat_to_sankey_with_loss(
  res_combined_pro$mat_prod,  res_combined_pro$mat_cons,
  scale = 1e9,  pop = res_combined_pro$pop,  pcap_label = "g protein/cap/day")

p_sankey_combined_kcal <- sankeyNetwork(
  Links = sankey_combined_kcal$links, Nodes = sankey_combined_kcal$nodes,
  Source = "source", Target = "target", Value = "value", NodeID = "name",
  NodeGroup = "group", colourScale = sankey_colour_scale,
  sinksRight = FALSE, fontSize = 13, nodeWidth = 20, nodePadding = 10,
  units = "Tcal", iterations = 0
) %>% add_pcap_tooltip(sankey_combined_kcal) %>% add_continent_legend()

p_sankey_combined_pro <- sankeyNetwork(
  Links = sankey_combined_pro$links, Nodes = sankey_combined_pro$nodes,
  Source = "source", Target = "target", Value = "value", NodeID = "name",
  NodeGroup = "group", colourScale = sankey_colour_scale,
  sinksRight = FALSE, fontSize = 13, nodeWidth = 20, nodePadding = 10,
  units = "kt protein", iterations = 0
) %>% add_pcap_tooltip(sankey_combined_pro) %>% add_continent_legend()

htmlwidgets::saveWidget(p_sankey_combined_kcal, "results/sankey_combined_kcal.html",
                        selfcontained = FALSE)
htmlwidgets::saveWidget(p_sankey_combined_pro,  "results/sankey_combined_protein.html",
                        selfcontained = FALSE)

# Food- and non-food-sector labor time Sankeys — total, male, and female.
# Left nodes (prod): country (food) / EXIO region (non-food) where labor is
# expended in the food supply chain. Right nodes (cons): country whose food
# consumption demands that labor.
for (sector in c("food", "nonfood")) {
  is_nonfood <- sector == "nonfood"
  l_country <- if (sector == "food") l_food_country else l_nonfood_region

  for (hr_label in c("total", "hr_m", "hr_f")) {
    mat_hr <- switch(hr_label,
      total = l_country$hr_m + l_country$hr_f,
      hr_m  = l_country$hr_m,
      hr_f  = l_country$hr_f
    )
    pcap_label <- switch(hr_label,
      total = "min/cap/day (total)",
      hr_m  = "min/cap/day (male)",
      hr_f  = "min/cap/day (female)"
    )

    res <- if (is_nonfood) agg_to_region_sankey(agg_to_region_matrix(mat_hr)) else agg_to_country_sankey(mat_hr)
    # M.hr matrix: * 1e6 converts to hr, * 60 converts to min → pcap_scale = 6e7
    sk  <- mat_to_sankey(res$mat, scale = 1, pop = res$pop, pcap_label = pcap_label, digits = 2,
                         pcap_scale = 6e7)
    p   <- sankeyNetwork(
      Links = sk$links, Nodes = sk$nodes,
      Source = "source", Target = "target", Value = "value", NodeID = "name",
      NodeGroup = "group", colourScale = sankey_colour_scale,
      sinksRight = FALSE, fontSize = 13, nodeWidth = 20, nodePadding = 10,
      units = "Mhr", iterations = 0
    ) %>% add_pcap_tooltip(sk) %>% add_continent_legend()

    htmlwidgets::saveWidget(p, paste0("results/sankey_hr_", sector, "_", hr_label, ".html"),
                            selfcontained = FALSE)
  }
}

# Sankeys for energy (TJ→PJ, ÷1e3) by food/non-food sector
for (sector in c("food", "nonfood")) {
  is_nonfood = sector == "nonfood"
  l_country = if (sector == "food") l_food_country else l_nonfood_region
  res = if (is_nonfood) agg_to_region_sankey(agg_to_region_matrix(l_country[["en"]])) else agg_to_country_sankey(l_country[["en"]])
  sk  = mat_to_sankey(res$mat, scale = 1e3)
  p   = sankeyNetwork(
    Links = sk$links, Nodes = sk$nodes,
    Source = "source", Target = "target", Value = "value", NodeID = "name",
    NodeGroup = "group", colourScale = sankey_colour_scale,
    sinksRight = FALSE, fontSize = 13, nodeWidth = 20, nodePadding = 10,
    units = "PJ", iterations = 0
  ) %>% add_continent_legend()
  htmlwidgets::saveWidget(p, paste0("results/sankey_", sector, "_en.html"),
                          selfcontained = FALSE)
}

# Country-wise production/consumption/export/import summary (kcal/capita/day)
prod_cty <- data.frame(iso3c = rownames(mat),
                       dom_consump = diag(mat), 
                       export = -(rowSums(mat) - diag(mat)),
                       import = colSums(mat) - diag(mat)) %>%
  mutate(export_perc = export / (dom_consump+import) * 100) %>%
  # join population
  right_join(pop_y, by = "iso3c") %>%
  mutate(prod_pday = dom_consump / pop * 1000 / 365, # Mcal to kcal
         imp_pday = import / pop * 1000 / 365, # Mcal to kcal
         exp_pday = export / pop * 1000 / 365) %>%
  mutate(supp_pday = prod_pday+imp_pday) %>% # Mcal to kcal
  filter(!iso3c %in% c("ROW", "ANT")) %>% # countries with no population data
  filter(!iso3c %in% c("SGP", "MDV"), prod_pday > 1)  # countries not interesting

# Draw a multi-panel grouped bar graph, showing production and export (for top 20 countries with highest export_perc)
prod_cty %>%
  arrange(desc(supp_pday)) %>%
  slice_head(n=25) %>%  
  # order bars by supp_pday (=prod+imp)
  mutate(iso3c = factor(iso3c, levels = iso3c[order(supp_pday, decreasing = TRUE)])) %>%
  
  pivot_longer(cols = c("prod_pday", "imp_pday", "exp_pday"), names_to = "type", values_to = "kcal_pday") %>%
  mutate(#iso3c = factor(iso3c, levels = iso3c[order(export_perc, decreasing = TRUE)]),
    type = factor(type, levels = c("prod_pday", "imp_pday", "exp_pday"), labels = c("Own consumption", "Import", "Export"))) %>%
  ggplot(aes(x=iso3c, y=kcal_pday, fill=type)) +
  geom_bar(stat="identity") +
  # coord_flip() +
  labs(x="Country (ISO3)", y="kcal/capita/day", fill="Type",
       title=paste0("Top 25 countries in total kcal footprint (", year, ")")) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_fill_manual(values=c("Own consumption"="#1f77b4", "Import"="#2ca02c", "Export"="#ff7f0e")) +
  # Tilt x-axis labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(-40000, 30000) 

ggsave(file.path(output_dir, paste0("plots/top25_", year, "_total.png")), width=10, height=6)
# geom_text(aes(label=round(kcal_pday,1)), position=position_stack(vjust=0.5), size=3)

# Make a similar graph for the bottom 20 countries with lowest export_perc (but prod_pday > 5000)
prod_cty %>%
  arrange(supp_pday) %>%
  slice_head(n=25) %>%
  # order bars by supp_pday (=prod+imp)
  mutate(iso3c = factor(iso3c, levels = iso3c[order(supp_pday, decreasing = TRUE)])) %>%
  pivot_longer(cols = c("prod_pday", "imp_pday", "exp_pday"), names_to = "type", values_to = "kcal_pday") %>%
  mutate(#iso3c = factor(iso3c, levels = iso3c[order(export_perc, decreasing = TRUE)]),
    type = factor(type, levels = c("prod_pday", "imp_pday", "exp_pday"), labels = c("Own consumption", "Import", "Export"))) %>%
  ggplot(aes(x=iso3c, y=kcal_pday, fill=type)) +
  geom_bar(stat="identity") +
  # coord_flip() +
  labs(x="Country (ISO3)", y="kcal/capita/day", fill="Type",
       title=paste0("Bottom 25 countries in total kcal footprint (", year, ")")) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_fill_manual(values=c("Own consumption"="#1f77b4", "Import"="#2ca02c", "Export"="#ff7f0e")) +
  # Tilt x-axis labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylim(-40000, 30000) 

ggsave(file.path(output_dir, paste0("plots/bottom25_", year, "_total.png")), width=10, height=6)

prod_cty %>%
  arrange(supp_pday) %>%
  # order bars by supp_pday (=prod+imp)
  mutate(iso3c = factor(iso3c, levels = iso3c[order(supp_pday, decreasing = TRUE)])) %>%
  pivot_longer(cols = c("prod_pday", "imp_pday", "exp_pday"), names_to = "type", values_to = "kcal_pday") %>%
  mutate(#iso3c = factor(iso3c, levels = iso3c[order(export_perc, decreasing = TRUE)]),
    type = factor(type, levels = c("prod_pday", "imp_pday", "exp_pday"), labels = c("Own consumption", "Import", "Export"))) %>%
  ggplot(aes(x=iso3c, y=kcal_pday, fill=type)) +
  geom_bar(stat="identity") +
  # coord_flip() +
  labs(x="Country (ISO3)", y="kcal/capita/day", fill="Type",
       title=paste0("Domestic supply of total kcal footprint (", year, ")")) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_fill_manual(values=c("Own consumption"="#1f77b4", "Import"="#2ca02c", "Export"="#ff7f0e")) +
  # Tilt x-axis labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



######### Total consumption analysis #################

mat_y = readRDS(file.path(output_dir,"calories/calorie_trade_mat_cons.rds"))

cons_cty <- data.frame(iso3c = colnames(mat_y), 
                       dom_consumption = diag(mat_y), 
                       export = -(rowSums(mat_y) - diag(mat_y)),
                       import = colSums(mat_y) - diag(mat_y)) %>%
  mutate(import_perc = import / (dom_consumption+import) * 100) %>%
  # join population
  right_join(pop_y, by = "iso3c") %>%
  mutate(dom_cons_pday = (dom_consumption) / pop * 1000 / 365,
         imp_pday = (import) / pop * 1000 / 365,
         exp_pday = (export) / pop * 1000 / 365)  %>% # Mcal to kcal/day/cap
  mutate(supp_pday = dom_cons_pday+imp_pday) %>%
  filter(!iso3c %in% c("ROW", "ANT")) %>% # countries with no population data
  filter(!iso3c %in% c("SGP", "MDV"), dom_cons_pday+imp_pday > 1000)

# Draw a similar plot for consumption and import (stacked); and export (for top 20 countries with highest import_perc)
cons_cty  %>% # countries not interesting
  arrange(desc(supp_pday)) %>%
  slice_head(n=25) %>%
  mutate(iso3c = factor(iso3c, levels = iso3c[order(supp_pday, decreasing = TRUE)])) %>%
  pivot_longer(cols = c("dom_cons_pday", "imp_pday", "exp_pday"), names_to = "type", values_to = "kcal_pday") %>%
  mutate(type = factor(type, levels = c("dom_cons_pday", "imp_pday", "exp_pday"), 
                       labels = c("Dom. consumption", "Import", "Export"))) %>%
  ggplot(aes(x=iso3c, y=kcal_pday, fill=type)) +
  geom_bar(stat="identity") +
  # coord_flip() +
  labs(x="Country (ISO3)", y="kcal/capita/day", fill="Type",
       title=paste0("Top 20 countries in final kcal consumption (", year, ")")) +
  theme_minimal() +
  theme(legend.position = "top") +
  ylim(-12000, 8000) +
  scale_fill_manual(values=c("Dom. consumption"="#1f77b4", "Import"="#2ca02c", "Export"="#ff7f0e")) +
  # Tilt x-axis labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+
# geom_text(aes(label=round(kcal_pday,1)), position=position_stack(vjust=0.5), size=3)

ggsave(file.path(output_dir, paste0("plots/top25_", year, "_final_cons.png")), width=10, height=6)

# Make a similar graph for the bottom 20 countries with lowest import_perc (but cons_pday > 2000)
cons_cty %>%
  arrange(supp_pday) %>%
  slice_head(n=25) %>%
  mutate(iso3c = factor(iso3c, levels = iso3c[order(supp_pday, decreasing = TRUE)])) %>%
  pivot_longer(cols = c("dom_cons_pday", "imp_pday", "exp_pday"), names_to = "type", values_to = "kcal_pday") %>%
  mutate(type = factor(type, levels = c("dom_cons_pday", "imp_pday", "exp_pday"), 
                       labels = c("Dom. consumption", "Import", "Export"))) %>%
  ggplot(aes(x=iso3c, y=kcal_pday, fill=type)) +
  geom_bar(stat="identity") +
  # coord_flip() +
  labs(x="Country (ISO3)", y="kcal/capita/day", fill="Type",
       title=paste0("Bottom 20 countries in final kcal consumption (", year, ")")) +
  theme_minimal() +
  theme(legend.position = "top") +
  ylim(-12000, 8000) +
  scale_fill_manual(values=c("Dom. consumption"="#1f77b4", "Import"="#2ca02c", "Export"="#ff7f0e")) +
  # Tilt x-axis labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+
# geom_text(aes(label=round(kcal_pday,1)), position=position_stack(vjust=0.5), size=3)

ggsave(file.path(output_dir, paste0("plots/bottom25_", year, "_final_cons.png")), width=10, height=6)


# Societal food loss analysis ####
df_loss = prod_cty %>% left_join(cons_cty, by="iso3c") %>%
  select(iso3c, dom_consump, dom_consumption) %>%
  mutate(loss_perc = (dom_consump - dom_consumption) / dom_consump * 100) %>%
  arrange(desc(loss_perc))


