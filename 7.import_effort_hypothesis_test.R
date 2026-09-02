#### H1/H2/H3 test: why is the "Import effort" energy~time scatter for       ####
#### imported protein so linear across very different importers?             ####
#
# Updated after fixing the RoW-duplication bug in
# 1.1.mrio_convert_indirect.R / 2.analyze_result.R (see chat): l_int_i's
# origin/row axis is now natively EXIO-region-resolved (49 rows: 44 named
# countries + 5 RoW aggregates), not pasted out to 187 FABIO-country rows.
# Consequently fp_impcons$nonfood is now (49 * n_nf) x 187, not
# (187 * n_nf) x 187 -- every computation below that touches non-food data
# had to change shape accordingly. fp_impcons$food is UNCHANGED (still
# 23001 x 187, FABIO-commodity-resolved, no RoW ambiguity of its own -- see
# chat: l_int_d is a per-row diagonal multiply, never had the duplication
# issue).
#
# Three hypotheses tested for the tight linear pattern in the "import effort"
# panel (energy vs. time per 50g of imported-consumed protein):
#  H1 (supply-side). Dominant global exporters share a similar energy/labor
#     ratio (similarly mechanized, capital-intensive production).
#  H2 (demand-side). A country's livestock-vs-plant import mix drives both
#     axes together (feed-conversion inefficiency).
#  H3 (importer's own tech level). An importer's OWN mechanization level
#     (proxied by its own energy/labor ratio as a producer/exporter) predicts
#     where it lands, rather than who it trades with or what it eats.
#
# Per chat discussion: anything combining food+non-food effort is resolved at
# EXIO-region level (49 rows) -- summing e.g. Somalia's and Ethiopia's
# individually-correct food footprints into one "RoW Africa" total is a
# legitimate aggregation, but non-food (and, to match it, protein/kcal here)
# can't honestly be attributed below region granularity for RoW aggregates.
# Food-ONLY figures stay at full 187-country resolution where shown.
#
# Standalone script -- only needs lightweight FABIO metadata plus
# data/fp_impcons_<year>.rds (already used throughout) and FABIO's raw Y.rds.

library(tidyverse)
library(Matrix)
library(data.table)
library(ggrepel)

year = 2020
FABIO_path = "H:/MyDocuments/Data/FABIO/input/"

regions <- fread(file.path(FABIO_path, "regions.csv"))
items   <- fread(file.path(FABIO_path, "items.csv"))
nrreg <- nrow(regions)
nrcom <- nrow(items)
stopifnot(nrreg == 187, nrcom == 123)

reg_map <- readxl::read_xlsx("data/fabio-exiobase.xlsx", sheet = "regions_concordance", col_names = TRUE)
FABIO_reg <- readxl::read_xlsx(file.path(FABIO_path, "fabio_classifications_v2.xlsx"), sheet = "Countries") %>%
  select(-area) %>% rename(ISO = iso3c, FAO_code = area_code) %>% left_join(reg_map) %>%
  mutate(EXIOBASE_code = as.numeric(ifelse(EXIOBASE_code == "NA", 47, EXIOBASE_code)),
         EXIOBASE = ifelse(EXIOBASE == "NA", "RoW Europe", EXIOBASE))

# Data-quality fix: code 41 (TWN's own individually-modeled EXIO region) is
# mislabeled "RoW Asia and Pacific" in the source concordance -- collides
# with the genuine RoW Asia & Pacific aggregate (code 45), breaking anything
# that treats region names as unique and wrongly flagging Taiwan as a RoW
# aggregate member below. Taiwan is one of EXIOBASE's individually-modeled
# economies; correct its label.
stopifnot(identical(sort(unique(FABIO_reg$EXIOBASE[FABIO_reg$EXIOBASE_code == 41])), "RoW Asia and Pacific"))
FABIO_reg$EXIOBASE[FABIO_reg$EXIOBASE_code == 41] <- "Taiwan"

stopifnot(identical(regions$iso3c, FABIO_reg$ISO))
row_countries <- FABIO_reg$ISO[grepl("RoW", FABIO_reg$EXIOBASE)]

n_reg_EXIO <- length(unique(FABIO_reg$EXIOBASE_code))
region_row_of_country <- FABIO_reg$EXIOBASE_code
region_members <- split(seq_len(nrreg), region_row_of_country)  # region index (character) -> FABIO column indices
region_name_of_index <- FABIO_reg$EXIOBASE[match(seq_len(n_reg_EXIO), FABIO_reg$EXIOBASE_code)]
is_row_region <- grepl("RoW", region_name_of_index)

# Item classification: animal-origin food commodities (meat/dairy/eggs/fish +
# the "Livestock" live-animal rows, which carry the labor/energy of raising
# animals even though final protein is recorded under "Livestock products").
# Unaffected by the fix -- food stays FABIO-commodity-resolved throughout.
group_food    <- rep(items$group, times = nrreg)
animal_origin <- group_food %in% c("Livestock", "Livestock products", "Fish")

producer_idx <- rep(seq_len(nrreg), each = nrcom)   # producer-country row index, 1..187, food-sector rows only

# ---- Aggregation helpers -----------------------------------------------

# Food: still 187 FABIO-country row-blocks (unchanged shape).
agg_country_footprint <- function(mat) {
  mat_country = matrix(0, nrow = nrreg, ncol = nrreg)
  nsect = nrow(mat) / nrreg
  for (i in 1:nrreg) mat_country[i, ] = colSums(mat[((i - 1) * nsect + 1):(i * nsect), ])
  rownames(mat_country) = colnames(mat_country) = regions$iso3c
  mat_country
}

# Non-food: now 49 EXIO-region row-blocks (fixed shape) x 187 FABIO-country columns.
agg_exio_region_footprint <- function(mat) {
  mat_region = matrix(0, nrow = n_reg_EXIO, ncol = nrreg)
  nsect = nrow(mat) / n_reg_EXIO
  stopifnot(nsect == round(nsect))
  for (i in 1:n_reg_EXIO) mat_region[i, ] = colSums(mat[((i - 1) * nsect + 1):(i * nsect), ])
  rownames(mat_region) = region_name_of_index
  colnames(mat_region) = regions$iso3c
  mat_region
}

# Collapse a matrix's row and/or column axis from FABIO-country (187) to
# EXIO-region (49) resolution, summing member countries within each region
# block. An axis already region-resolved is left as-is. Mirrors
# 2.analyze_result.R's collapse_axis_to_region()/agg_to_region_matrix().
collapse_axis_to_region <- function(mat, axis) {
  n = if (axis == 1) nrow(mat) else ncol(mat)
  if (n == n_reg_EXIO) return(mat)
  stopifnot(n == nrreg)
  if (axis == 1) {
    out = matrix(0, nrow = n_reg_EXIO, ncol = ncol(mat))
    for (i in seq_len(n_reg_EXIO)) out[i, ] = colSums(mat[region_members[[as.character(i)]], , drop = FALSE])
    rownames(out) = region_name_of_index
    colnames(out) = colnames(mat)
  } else {
    out = matrix(0, nrow = nrow(mat), ncol = n_reg_EXIO)
    for (i in seq_len(n_reg_EXIO)) out[, i] = rowSums(mat[, region_members[[as.character(i)]], drop = FALSE])
    rownames(out) = rownames(mat)
    colnames(out) = region_name_of_index
  }
  out
}
agg_to_region_matrix <- function(mat) {
  m = collapse_axis_to_region(mat, axis = 1)
  m = collapse_axis_to_region(m, axis = 2)
  rownames(m) = colnames(m) = region_name_of_index
  m
}

cat("Loading fp_impcons (~95MB)...\n")
fp_impcons <- readRDS(paste0("data/fp_impcons_", year, ".rds"))

hr_food_country    <- agg_country_footprint(fp_impcons$food$hr_m    + fp_impcons$food$hr_f)      # 187 x 187, M.hr
en_food_country    <- agg_country_footprint(fp_impcons$food$en)                                  # 187 x 187, TJ
hr_nonfood_region  <- agg_exio_region_footprint(fp_impcons$nonfood$hr_m + fp_impcons$nonfood$hr_f)  # 49 x 187, M.hr
en_nonfood_region  <- agg_exio_region_footprint(fp_impcons$nonfood$en)                            # 49 x 187, TJ

# Food, region-region (49x49) for the combined "total" ranking/response below.
hr_food_region_sq <- agg_to_region_matrix(hr_food_country)
en_food_region_sq <- agg_to_region_matrix(en_food_country)

stopifnot(identical(rownames(hr_food_region_sq), rownames(hr_nonfood_region)))

# Non-food's row axis is already region-resolved (49); collapse its column
# axis to region too, for a symmetric 49x49 comparable to the food matrices.
hr_nonfood_region_sq <- collapse_axis_to_region(hr_nonfood_region, axis = 2)
en_nonfood_region_sq <- collapse_axis_to_region(en_nonfood_region, axis = 2)

hr_total_region <- hr_food_region_sq + hr_nonfood_region_sq  # 49 x 49
en_total_region <- en_food_region_sq + en_nonfood_region_sq  # 49 x 49


#### Diagnostic: shape check + food- vs. non-food-sector share of "import effort" ####

cat(sprintf("\nnrow(fp_impcons$nonfood$hr_m) = %d (= n_reg_EXIO(%d) x n_nf) -- confirms the fixed, region-resolved shape.\n",
           nrow(fp_impcons$nonfood$hr_m), n_reg_EXIO))

hr_food_import_by_region <- colSums(hr_food_region_sq) - diag(hr_food_region_sq)
hr_nonfood_import_by_region <- colSums(hr_nonfood_region_sq) - diag(hr_nonfood_region_sq)

df_sector_share <- data.frame(
  exio_region = rownames(hr_food_region_sq),
  hr_nonfood_share = hr_nonfood_import_by_region / (hr_food_import_by_region + hr_nonfood_import_by_region),
  is_row = is_row_region
) %>% filter(!is_row, is.finite(hr_nonfood_share))

cat("\n---- Food- vs. non-food-sector share of 'import effort' (region-resolved) ----\n")
cat(sprintf("  median non-food share of hours: %.1f%%  (range %.1f - %.1f%%)\n",
           median(df_sector_share$hr_nonfood_share) * 100,
           min(df_sector_share$hr_nonfood_share) * 100, max(df_sector_share$hr_nonfood_share) * 100))

row_idx <- which(is_row_region)
row_hr_share_by_region <- colSums(hr_nonfood_region_sq[row_idx, , drop = FALSE]) / (colSums(hr_nonfood_region_sq) - diag(hr_nonfood_region_sq))
cat(sprintf("\n---- RoW-aggregate share of non-food import-effort hours (post-fix; compare to pre-fix 83-93%%) ----\n"))
cat(sprintf("  median: %.1f%%  (range %.1f - %.1f%%)\n",
           median(row_hr_share_by_region[!is_row_region], na.rm = TRUE) * 100,
           min(row_hr_share_by_region[!is_row_region], na.rm = TRUE) * 100,
           max(row_hr_share_by_region[!is_row_region], na.rm = TRUE) * 100))


#### H1: do the biggest suppliers of "import effort" share a similar energy/labor ratio? ####

exported_hr_total <- rowSums(hr_total_region) - diag(hr_total_region)
exported_en_total <- rowSums(en_total_region) - diag(en_total_region)

df_exporters <- data.frame(exio_region = rownames(hr_total_region), exported_hr = exported_hr_total, exported_en = exported_en_total,
                           ratio = exported_en_total / exported_hr_total, is_row = is_row_region) %>%
  filter(is.finite(ratio), exported_hr > 0)

cv <- function(x) sd(x) / mean(x)
top_n <- 20
top_exporters <- df_exporters %>% filter(!is_row) %>% arrange(desc(exported_hr)) %>% slice_head(n = top_n)
cat("\nTop", top_n, "non-RoW suppliers of TOTAL (food+non-food) 'import effort':\n")
print(top_exporters %>% mutate(across(where(is.numeric), ~round(.x, 3))))

pool <- df_exporters %>% filter(!is_row, !exio_region %in% top_exporters$exio_region)
set.seed(1)
rand_cv <- replicate(2000, cv(sample(pool$ratio, top_n)))
cat(sprintf("\nCV(energy/labor ratio): top %d = %.3f | random draws (n=2000): mean %.3f, sd %.3f\n",
           top_n, cv(top_exporters$ratio), mean(rand_cv), sd(rand_cv)))
cat(sprintf("Empirical p-value: %.3f\n", mean(rand_cv <= cv(top_exporters$ratio))))

p_h1 <- ggplot(top_exporters, aes(x = exported_hr, y = exported_en)) +
  geom_point(size = 3, color = "#1f77b4") +
  geom_text_repel(aes(label = exio_region), size = 3.5) +
  geom_smooth(method = "lm", se = FALSE, color = "grey50", linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Total exported effort (M.hr, log)", y = "Total exported energy (TJ, log)",
       title = paste0("H1 -- energy/labor ratio among top ", top_n, " TOTAL import-effort suppliers")) +
  theme_minimal()
print(p_h1)
ggsave("results/h1_exporter_ratio_clustering.pdf", p_h1, width = 8, height = 6)

# H1b: food-sector-only ranking (avoids service-trade contamination).
exported_hr_food <- rowSums(hr_food_region_sq) - diag(hr_food_region_sq)
exported_en_food <- rowSums(en_food_region_sq) - diag(en_food_region_sq)
df_exporters_food <- data.frame(exio_region = rownames(hr_food_region_sq), exported_hr = exported_hr_food, exported_en = exported_en_food,
                                ratio = exported_en_food / exported_hr_food, is_row = is_row_region) %>%
  filter(is.finite(ratio), exported_hr > 0)
top_exporters_food <- df_exporters_food %>% filter(!is_row) %>% arrange(desc(exported_hr)) %>% slice_head(n = top_n)
cat("\nTop", top_n, "non-RoW FOOD-SECTOR-ONLY import-effort suppliers:\n")
print(top_exporters_food %>% mutate(across(where(is.numeric), ~round(.x, 3))))
pool_food <- df_exporters_food %>% filter(!is_row, !exio_region %in% top_exporters_food$exio_region)
set.seed(1); rand_cv_food <- replicate(2000, cv(sample(pool_food$ratio, top_n)))
cat(sprintf("[Food-only] CV: top = %.3f | random: mean %.3f | p = %.3f\n",
           cv(top_exporters_food$ratio), mean(rand_cv_food), mean(rand_cv_food <= cv(top_exporters_food$ratio))))

p_h1b <- ggplot(top_exporters_food, aes(x = exported_hr, y = exported_en)) +
  geom_point(size = 3, color = "#2ca02c") + geom_text_repel(aes(label = exio_region), size = 3.5) +
  geom_smooth(method = "lm", se = FALSE, color = "grey50", linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Food-sector exported effort (M.hr, log)", y = "Food-sector exported energy (TJ, log)",
       title = paste0("H1b -- energy/labor ratio among top ", top_n, " FOOD-sector suppliers")) +
  theme_minimal()
print(p_h1b)
ggsave("results/h1b_exporter_ratio_clustering_food_only.pdf", p_h1b, width = 8, height = 6)

# H1c: non-food-sector-only ranking (already natively region-resolved).
exported_hr_nonfood <- rowSums(hr_nonfood_region_sq) - diag(hr_nonfood_region_sq)
exported_en_nonfood <- rowSums(en_nonfood_region_sq) - diag(en_nonfood_region_sq)
df_exporters_nonfood <- data.frame(exio_region = rownames(hr_nonfood_region_sq), exported_hr = exported_hr_nonfood, exported_en = exported_en_nonfood,
                                   ratio = exported_en_nonfood / exported_hr_nonfood, is_row = is_row_region) %>%
  filter(is.finite(ratio), exported_hr > 0)
top_exporters_nonfood <- df_exporters_nonfood %>% filter(!is_row) %>% arrange(desc(exported_hr)) %>% slice_head(n = top_n)
cat("\nTop", top_n, "non-RoW NON-FOOD-SECTOR-ONLY import-effort suppliers:\n")
print(top_exporters_nonfood %>% mutate(across(where(is.numeric), ~round(.x, 3))))
pool_nonfood <- df_exporters_nonfood %>% filter(!is_row, !exio_region %in% top_exporters_nonfood$exio_region)
set.seed(1); rand_cv_nonfood <- replicate(2000, cv(sample(pool_nonfood$ratio, top_n)))
cat(sprintf("[Non-food-only] CV: top = %.3f | random: mean %.3f | p = %.3f\n",
           cv(top_exporters_nonfood$ratio), mean(rand_cv_nonfood), mean(rand_cv_nonfood <= cv(top_exporters_nonfood$ratio))))

p_h1c <- ggplot(top_exporters_nonfood, aes(x = exported_hr, y = exported_en)) +
  geom_point(size = 3, color = "#d62728") + geom_text_repel(aes(label = exio_region), size = 3.5) +
  geom_smooth(method = "lm", se = FALSE, color = "grey50", linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Non-food-sector exported effort (M.hr, log)", y = "Non-food-sector exported energy (TJ, log)",
       title = paste0("H1c -- energy/labor ratio among top ", top_n, " NON-FOOD-sector suppliers")) +
  theme_minimal()
print(p_h1c)
ggsave("results/h1c_exporter_ratio_clustering_nonfood_only.pdf", p_h1c, width = 8, height = 6)


#### H2: does the livestock/plant mix of each importer's import-effort explain the pattern? ####

# Animal/plant classifier: still FABIO-commodity-resolved (187 countries, no
# RoW ambiguity), restricted to off-diagonal (producer != consumer) rows.
zero_domestic_effort <- function(mat) {
  m <- as(mat, "TsparseMatrix")
  keep <- producer_idx[m@i + 1] != (m@j + 1)
  Matrix::sparseMatrix(i = m@i[keep] + 1, j = m@j[keep] + 1, x = m@x[keep], dims = dim(m))
}
food_hr_imp <- zero_domestic_effort(fp_impcons$food$hr_m + fp_impcons$food$hr_f)
food_en_imp <- zero_domestic_effort(fp_impcons$food$en)
hr_animal <- colSums(food_hr_imp[animal_origin, ]);  hr_plant <- colSums(food_hr_imp[!animal_origin, ])
en_animal <- colSums(food_en_imp[animal_origin, ]);  en_plant <- colSums(food_en_imp[!animal_origin, ])

# Region-level animal share: sum member countries' hr_animal/hr_plant, then
# take the ratio at region level (legitimate aggregation, food has no RoW
# ambiguity to begin with).
animal_share_region <- function(vals_animal, vals_plant) {
  a = sapply(seq_len(n_reg_EXIO), function(i) sum(vals_animal[region_members[[as.character(i)]]]))
  p = sapply(seq_len(n_reg_EXIO), function(i) sum(vals_plant[region_members[[as.character(i)]]]))
  a / (a + p)
}
hr_animal_share_region <- animal_share_region(hr_animal, hr_plant)
en_animal_share_region <- animal_share_region(en_animal, en_plant)

# Protein denominator, rebuilt from FABIO's raw Y (mirrors 1.mrio_convert.R's
# coeff_pro derivation; no need for the full pipeline).
cat("\nLoading FABIO Y.rds...\n")
FABIO_y <- readRDS(file.path(FABIO_path, "Y.rds"))[[as.character(year)]]
FABIO_y_hh <- FABIO_y[, grep("food", colnames(FABIO_y))]
stopifnot(nrow(FABIO_y_hh) == nrreg * nrcom, ncol(FABIO_y_hh) == nrreg)

coeff_pro_vec <- fread(file.path(FABIO_path, "nutrient_coefficients_protein.csv")) %>%
  mutate(protein_per_kg = 10 * protein_per_100g) %>% pull(protein_per_kg) %>% rep(times = nrreg)
coeff_pro_vec[rep(items$comm_group, times = nrreg) == "Live animals"] <- 0

y_hh_t <- as(FABIO_y_hh, "TsparseMatrix")
is_import_flow <- producer_idx[y_hh_t@i + 1] != (y_hh_t@j + 1)
Y_imp <- Matrix::sparseMatrix(i = y_hh_t@i[is_import_flow] + 1, j = y_hh_t@j[is_import_flow] + 1,
                              x = y_hh_t@x[is_import_flow], dims = dim(FABIO_y_hh))
Y_imp_pro <- Matrix::Diagonal(x = coeff_pro_vec) %*% Y_imp * 1000   # g protein

pro_import_g_country <- setNames(colSums(Y_imp_pro), regions$iso3c)
pro_import_g_region <- sapply(seq_len(n_reg_EXIO), function(i) sum(pro_import_g_country[region_members[[as.character(i)]]]))
names(pro_import_g_region) <- region_name_of_index

hr_import_total_region <- colSums(hr_total_region) - diag(hr_total_region)
en_import_total_region <- colSums(en_total_region) - diag(en_total_region)

df_conv <- data.frame(exio_region = rownames(hr_total_region),
                      hr_animal_share = hr_animal_share_region,
                      en_animal_share = en_animal_share_region,
                      hr_import_total = hr_import_total_region,
                      en_import_total = en_import_total_region) %>%
  mutate(pro_import_g = pro_import_g_region[exio_region],
         hr_per_50g_protein = hr_import_total * 1e6 / pro_import_g * 50,
         mj_per_50g_protein = en_import_total * 1e6 / pro_import_g * 50,
         is_row = is_row_region) %>%
  filter(is.finite(hr_per_50g_protein), is.finite(mj_per_50g_protein), pro_import_g > 0, !is_row)

cat("\n---- H2 correlations: animal-origin share vs. conversion factor (region-resolved) ----\n")
cat(sprintf("  hr_animal_share vs hr/50g protein: r = %.3f\n", cor(df_conv$hr_animal_share, df_conv$hr_per_50g_protein)))
cat(sprintf("  en_animal_share vs MJ/50g protein: r = %.3f\n", cor(df_conv$en_animal_share, df_conv$mj_per_50g_protein)))

p_h2 <- ggplot(df_conv, aes(x = mj_per_50g_protein, y = hr_per_50g_protein)) +
  geom_point(aes(color = hr_animal_share, size = pro_import_g), alpha = 0.85) +
  geom_text_repel(aes(label = exio_region), size = 3, max.overlaps = 15) +
  scale_color_viridis_c(option = "C", labels = scales::percent) +
  scale_size_continuous(range = c(1, 8), guide = "none") +
  labs(x = "Energy (MJ / 50 g imported-consumed protein)", y = "Time (hr / 50 g imported-consumed protein)",
       color = "Animal-origin\nshare of\nimport-effort hours",
       title = "H2 -- import-effort line colored by animal-vs-plant mix") +
  theme_minimal() + theme(legend.position = "right")
print(p_h2)
ggsave("results/h2_livestock_share_vs_line_position.pdf", p_h2, width = 10, height = 7)


#### H3: does an importer's OWN mechanization level explain where it sits on the line? ####

own_ratio <- setNames(df_exporters$ratio, df_exporters$exio_region)
df_h3 <- df_conv %>% mutate(own_ratio = own_ratio[exio_region]) %>% filter(is.finite(own_ratio))

cat("\n---- H3 correlations: own (exporter-side) energy/labor ratio vs. import-effort conversion factor ----\n")
cat(sprintf("  own_ratio vs hr/50g protein: r = %.3f\n", cor(df_h3$own_ratio, df_h3$hr_per_50g_protein)))
cat(sprintf("  own_ratio vs MJ/50g protein: r = %.3f\n", cor(df_h3$own_ratio, df_h3$mj_per_50g_protein)))
cat(sprintf("  own_ratio vs the LINE ITSELF (mj/hr): r = %.3f\n",
           cor(df_h3$own_ratio, df_h3$mj_per_50g_protein / df_h3$hr_per_50g_protein)))

p_h3 <- ggplot(df_h3, aes(x = mj_per_50g_protein, y = hr_per_50g_protein)) +
  geom_point(aes(color = own_ratio, size = pro_import_g), alpha = 0.85) +
  geom_text_repel(aes(label = exio_region), size = 3, max.overlaps = 15) +
  scale_color_viridis_c(option = "C", trans = "log10") +
  scale_size_continuous(range = c(1, 8), guide = "none") +
  labs(x = "Energy (MJ / 50 g imported-consumed protein)", y = "Time (hr / 50 g imported-consumed protein)",
       color = "Own (exporter-side)\nenergy/labor\nratio (MJ/hr, log)",
       title = "H3 -- import-effort line colored by importer's OWN mechanization level") +
  theme_minimal() + theme(legend.position = "right")
print(p_h3)
ggsave("results/h3_own_mechanization_vs_line_position.pdf", p_h3, width = 10, height = 7)

fwrite(top_exporters, "output/h1_top_exporters_ratio.csv")
fwrite(top_exporters_food, "output/h1b_top_exporters_ratio_food_only.csv")
fwrite(top_exporters_nonfood, "output/h1c_top_exporters_ratio_nonfood_only.csv")
fwrite(df_conv %>% select(exio_region, hr_animal_share, en_animal_share, pro_import_g, mj_per_50g_protein, hr_per_50g_protein),
      "output/h2_livestock_share_vs_conversion.csv")
fwrite(df_h3 %>% select(exio_region, own_ratio, mj_per_50g_protein, hr_per_50g_protein),
      "output/h3_own_mechanization_vs_conversion.csv")
