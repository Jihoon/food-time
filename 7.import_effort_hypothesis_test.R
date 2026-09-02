#### H1/H2 test: why is the "Import effort" energy~time scatter for imported  ####
#### protein so linear across very different importers?                      ####
#
# Two candidate explanations discussed for the pattern in
# p_tradeoff_protein_totalecon_import_effort (2.analyze_result.R ~L2532):
#
#  H1 (supply-side). The small set of countries that dominate GLOBAL protein
#     exports (grain/oilseed/meat: USA, BRA, RUS, AUS, ARG, FRA, CAN, ...)
#     themselves sit on a similar energy/labor ratio -- similarly mechanized,
#     capital-intensive export agriculture -- so no matter which mix of these
#     an importer draws from, the aggregate ratio comes out similar.
#
#  H2 (demand-side). A country's livestock-vs-plant protein IMPORT MIX is
#     itself a common driver: livestock products cost much more energy and
#     time per gram of protein than cereals/pulses (feed-conversion
#     inefficiency raises both together), so a country's animal-share of
#     imports alone could generate the same line, independent of which
#     specific exporters supply it.
#
# Standalone script -- only needs lightweight FABIO metadata plus two files:
#   data/fp_impcons_<year>.rds  -- energy/labor footprint of X_imp (each
#                                   country's IMPORTED food consumption),
#                                   item- and producer-country-resolved.
#                                   Saved by 2.analyze_result.R section 1.3.
#   FABIO's raw Y.rds (32MB)    -- rebuilds the protein denominator (g
#                                   protein/yr) for the SAME imported-
#                                   consumption bucket, mirroring
#                                   1.mrio_convert.R's coeff_pro derivation,
#                                   without touching the giant L/EXIO matrices.

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

# Flag "Rest of World" aggregate regions (many FABIO countries pasted to one
# EXIO value) so they don't masquerade as distinct exporters in H1.
reg_map <- readxl::read_xlsx("data/fabio-exiobase.xlsx", sheet = "regions_concordance", col_names = TRUE)
FABIO_reg <- readxl::read_xlsx(file.path(FABIO_path, "fabio_classifications_v2.xlsx"), sheet = "Countries") %>%
  select(-area) %>% rename(ISO = iso3c, FAO_code = area_code) %>% left_join(reg_map) %>%
  mutate(EXIOBASE_code = as.numeric(ifelse(EXIOBASE_code == "NA", 47, EXIOBASE_code)),
         EXIOBASE = ifelse(EXIOBASE == "NA", "RoW Europe", EXIOBASE))
stopifnot(identical(regions$iso3c, FABIO_reg$ISO))
row_countries <- FABIO_reg$ISO[grepl("RoW", FABIO_reg$EXIOBASE)]

# Item classification: animal-origin food commodities (meat/dairy/eggs/fish +
# the "Livestock" live-animal rows, which carry the labor/energy of raising
# animals even though final protein is recorded under "Livestock products").
group_food    <- rep(items$group, times = nrreg)
animal_origin <- group_food %in% c("Livestock", "Livestock products", "Fish")

producer_idx <- rep(seq_len(nrreg), each = nrcom)   # producer-country row index, 1..187, food-sector rows only

agg_country_footprint <- function(mat) {
  mat_country = matrix(0, nrow = nrreg, ncol = nrreg)
  nsect = nrow(mat) / nrreg
  for (i in 1:nrreg) mat_country[i, ] = colSums(mat[((i - 1) * nsect + 1):(i * nsect), ])
  rownames(mat_country) = colnames(mat_country) = regions$iso3c
  mat_country
}

cat("Loading fp_impcons (~95MB)...\n")
fp_impcons <- readRDS(paste0("data/fp_impcons_", year, ".rds"))  # list(food=list(hr_m,hr_f,en), nonfood=list(...))

hr_food_country    <- agg_country_footprint(fp_impcons$food$hr_m    + fp_impcons$food$hr_f)
hr_nonfood_country <- agg_country_footprint(fp_impcons$nonfood$hr_m + fp_impcons$nonfood$hr_f)
hr_total_country   <- hr_food_country + hr_nonfood_country                                    # M.hr

# food (23001 = 187x123 items) and non-food (32725 = 187x175 EXIO sectors) rows
# have different dimensions -- must aggregate each to 187x187 BEFORE summing.
en_food_country    <- agg_country_footprint(fp_impcons$food$en)
en_nonfood_country <- agg_country_footprint(fp_impcons$nonfood$en)
en_total_country   <- en_food_country + en_nonfood_country      # TJ


#### Diagnostic: how much of "import effort" is food-sector vs. non-food-sector? ####
# Matters a lot for how to read H1/H2 below: if non-food (freight, packaging,
# wholesale/retail trade, financial & insurance services, ...) dominates the
# total, then (a) H1 should really be tested on non-food suppliers, not food
# exporters, and (b) H2's animal-vs-plant classifier -- which only exists at
# food-item resolution, since fp_nonfood collapses the FABIO-commodity axis
# via l_int_i %*% X -- can only ever explain the (smaller) food-sector slice.

hr_food_import_vec    <- colSums(hr_food_country)    - diag(hr_food_country)
hr_nonfood_import_vec <- colSums(hr_nonfood_country) - diag(hr_nonfood_country)
en_food_import_vec    <- colSums(en_food_country)    - diag(en_food_country)
en_nonfood_import_vec <- colSums(en_nonfood_country) - diag(en_nonfood_country)

df_sector_share <- data.frame(
  country = regions$iso3c,
  hr_nonfood_share = hr_nonfood_import_vec / (hr_food_import_vec + hr_nonfood_import_vec),
  en_nonfood_share  = en_nonfood_import_vec / (en_food_import_vec + en_nonfood_import_vec),
  is_row = regions$iso3c %in% row_countries
) %>% filter(!is_row, is.finite(hr_nonfood_share), is.finite(en_nonfood_share))

cat("\n---- Food- vs. non-food-sector share of 'import effort' ----\n")
cat(sprintf("  median non-food share of hours:  %.1f%%  (range %.1f - %.1f%%)\n",
           median(df_sector_share$hr_nonfood_share) * 100,
           min(df_sector_share$hr_nonfood_share) * 100, max(df_sector_share$hr_nonfood_share) * 100))
cat(sprintf("  median non-food share of energy: %.1f%%  (range %.1f - %.1f%%)\n",
           median(df_sector_share$en_nonfood_share) * 100,
           min(df_sector_share$en_nonfood_share) * 100, max(df_sector_share$en_nonfood_share) * 100))
fwrite(df_sector_share, "output/food_vs_nonfood_share_import_effort.csv")


#### H1: do the biggest suppliers of "import effort" share a similar energy/labor ratio? ####

exported_hr <- rowSums(hr_total_country) - diag(hr_total_country)
exported_en <- rowSums(en_total_country) - diag(en_total_country)

df_exporters <- data.frame(country = regions$iso3c, exported_hr, exported_en,
                           ratio = exported_en / exported_hr,          # TJ/M.hr == MJ/hr
                           is_row = regions$iso3c %in% row_countries) %>%
  filter(is.finite(ratio), exported_hr > 0)

top_n <- 20
top_exporters <- df_exporters %>% filter(!is_row) %>% arrange(desc(exported_hr)) %>% slice_head(n = top_n)
cat("\nTop", top_n, "non-RoW suppliers of 'import effort' (any country's imported-protein bucket):\n")
print(top_exporters %>% mutate(across(where(is.numeric), ~round(.x, 3))))

cv <- function(x) sd(x) / mean(x)
pool <- df_exporters %>% filter(!is_row, !country %in% top_exporters$country)
set.seed(1)
rand_cv <- replicate(2000, cv(sample(pool$ratio, top_n)))

cat(sprintf("\nCV(energy/labor ratio): top %d exporters = %.3f | random %d-country draws (n=2000): mean %.3f, sd %.3f\n",
           top_n, cv(top_exporters$ratio), top_n, mean(rand_cv), sd(rand_cv)))
cat(sprintf("Empirical p-value (share of random draws with CV <= top exporters' CV): %.3f\n",
           mean(rand_cv <= cv(top_exporters$ratio))))

p_h1 <- ggplot(top_exporters, aes(x = exported_hr, y = exported_en)) +
  geom_point(size = 3, color = "#1f77b4") +
  geom_text_repel(aes(label = country), size = 3.5) +
  geom_smooth(method = "lm", se = FALSE, color = "grey50", linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Total exported effort embodied in others' imports (M.hr, log)",
       y = "Total exported energy embodied in others' imports (TJ, log)",
       title = paste0("H1 -- energy/labor ratio among top ", top_n, " non-RoW import-effort suppliers")) +
  theme_minimal()
print(p_h1)
ggsave("results/h1_exporter_ratio_clustering.pdf", p_h1, width = 8, height = 6)

# H1b: same test, but ranked/measured on FOOD-SECTOR-ONLY effort. The total
# (food+non-food) ranking above can be dominated by embodied SERVICES (freight,
# insurance, finance) rather than actual grain/oilseed/meat exports -- e.g. GBR,
# GRC, ZAF showing up as top-20 "import effort suppliers" by total effort is a
# services-trade artifact, not agricultural export dominance, and confounds a
# fair test of the original hypothesis ("dominant global exporters of traded
# protein commodities... cluster on a similar energy/labor ratio").
exported_hr_food <- rowSums(hr_food_country) - diag(hr_food_country)
exported_en_food <- rowSums(en_food_country) - diag(en_food_country)

df_exporters_food <- data.frame(country = regions$iso3c, exported_hr = exported_hr_food, exported_en = exported_en_food,
                                ratio = exported_en_food / exported_hr_food,
                                is_row = regions$iso3c %in% row_countries) %>%
  filter(is.finite(ratio), exported_hr > 0)

top_exporters_food <- df_exporters_food %>% filter(!is_row) %>% arrange(desc(exported_hr)) %>% slice_head(n = top_n)
cat("\nTop", top_n, "non-RoW suppliers of FOOD-SECTOR-ONLY import effort:\n")
print(top_exporters_food %>% mutate(across(where(is.numeric), ~round(.x, 3))))

pool_food <- df_exporters_food %>% filter(!is_row, !country %in% top_exporters_food$country)
set.seed(1)
rand_cv_food <- replicate(2000, cv(sample(pool_food$ratio, top_n)))

cat(sprintf("\n[Food-sector-only] CV(energy/labor ratio): top %d exporters = %.3f | random draws (n=2000): mean %.3f, sd %.3f\n",
           top_n, cv(top_exporters_food$ratio), mean(rand_cv_food), sd(rand_cv_food)))
cat(sprintf("[Food-sector-only] Empirical p-value: %.3f\n",
           mean(rand_cv_food <= cv(top_exporters_food$ratio))))

p_h1b <- ggplot(top_exporters_food, aes(x = exported_hr, y = exported_en)) +
  geom_point(size = 3, color = "#2ca02c") +
  geom_text_repel(aes(label = country), size = 3.5) +
  geom_smooth(method = "lm", se = FALSE, color = "grey50", linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Food-sector exported effort embodied in others' imports (M.hr, log)",
       y = "Food-sector exported energy embodied in others' imports (TJ, log)",
       title = paste0("H1b -- energy/labor ratio among top ", top_n, " non-RoW FOOD-sector suppliers")) +
  theme_minimal()
print(p_h1b)
ggsave("results/h1b_exporter_ratio_clustering_food_only.pdf", p_h1b, width = 8, height = 6)

# H1c: same test again, but for NON-FOOD-SECTOR-ONLY effort. Given the
# diagnostic above, this -- not H1b -- is likely the more relevant supply-side
# test: do the big suppliers of embodied non-food effort (freight, packaging,
# wholesale/retail trade, financial & insurance services, ...) cluster on a
# similar energy/labor ratio, the way H1 originally asked about grain/meat
# exporters specifically?
exported_hr_nonfood <- rowSums(hr_nonfood_country) - diag(hr_nonfood_country)
exported_en_nonfood <- rowSums(en_nonfood_country) - diag(en_nonfood_country)

df_exporters_nonfood <- data.frame(country = regions$iso3c, exported_hr = exported_hr_nonfood, exported_en = exported_en_nonfood,
                                   ratio = exported_en_nonfood / exported_hr_nonfood,
                                   is_row = regions$iso3c %in% row_countries) %>%
  filter(is.finite(ratio), exported_hr > 0)

top_exporters_nonfood <- df_exporters_nonfood %>% filter(!is_row) %>% arrange(desc(exported_hr)) %>% slice_head(n = top_n)
cat("\nTop", top_n, "non-RoW suppliers of NON-FOOD-SECTOR-ONLY import effort:\n")
print(top_exporters_nonfood %>% mutate(across(where(is.numeric), ~round(.x, 3))))

pool_nonfood <- df_exporters_nonfood %>% filter(!is_row, !country %in% top_exporters_nonfood$country)
set.seed(1)
rand_cv_nonfood <- replicate(2000, cv(sample(pool_nonfood$ratio, top_n)))

cat(sprintf("\n[Non-food-only] CV(energy/labor ratio): top %d exporters = %.3f | random draws (n=2000): mean %.3f, sd %.3f\n",
           top_n, cv(top_exporters_nonfood$ratio), mean(rand_cv_nonfood), sd(rand_cv_nonfood)))
cat(sprintf("[Non-food-only] Empirical p-value: %.3f\n",
           mean(rand_cv_nonfood <= cv(top_exporters_nonfood$ratio))))

p_h1c <- ggplot(top_exporters_nonfood, aes(x = exported_hr, y = exported_en)) +
  geom_point(size = 3, color = "#d62728") +
  geom_text_repel(aes(label = country), size = 3.5) +
  geom_smooth(method = "lm", se = FALSE, color = "grey50", linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Non-food-sector exported effort embodied in others' imports (M.hr, log)",
       y = "Non-food-sector exported energy embodied in others' imports (TJ, log)",
       title = paste0("H1c -- energy/labor ratio among top ", top_n, " non-RoW NON-FOOD-sector suppliers")) +
  theme_minimal()
print(p_h1c)
ggsave("results/h1c_exporter_ratio_clustering_nonfood_only.pdf", p_h1c, width = 8, height = 6)


#### H2: does the livestock/plant mix of each importer's import-effort explain the pattern? ####

# Restrict to off-diagonal (producer != consumer) rows before splitting by item
# group, so the classifier matches the "import effort" bucket exactly.
zero_domestic_effort <- function(mat) {
  m <- as(mat, "TsparseMatrix")
  keep <- producer_idx[m@i + 1] != (m@j + 1)
  Matrix::sparseMatrix(i = m@i[keep] + 1, j = m@j[keep] + 1, x = m@x[keep], dims = dim(m))
}

food_hr_imp <- zero_domestic_effort(fp_impcons$food$hr_m + fp_impcons$food$hr_f)
food_en_imp <- zero_domestic_effort(fp_impcons$food$en)

hr_animal <- colSums(food_hr_imp[animal_origin, ]);  hr_plant <- colSums(food_hr_imp[!animal_origin, ])
en_animal <- colSums(food_en_imp[animal_origin, ]);  en_plant <- colSums(food_en_imp[!animal_origin, ])

df_composition <- data.frame(country = regions$iso3c, hr_animal, hr_plant, en_animal, en_plant) %>%
  mutate(hr_animal_share = hr_animal / (hr_animal + hr_plant),
         en_animal_share = en_animal / (en_animal + en_plant))

# Protein denominator for the SAME imported-consumption bucket, rebuilt
# directly from FABIO's raw Y (no need for the full pipeline -- mirrors
# 1.mrio_convert.R's coeff_pro derivation).
cat("\nLoading FABIO Y.rds...\n")
FABIO_y <- readRDS(file.path(FABIO_path, "Y.rds"))[[as.character(year)]]
FABIO_y_hh <- FABIO_y[, grep("food", colnames(FABIO_y))]   # 23001 x 187
stopifnot(nrow(FABIO_y_hh) == nrreg * nrcom, ncol(FABIO_y_hh) == nrreg)

coeff_pro_vec <- fread(file.path(FABIO_path, "nutrient_coefficients_protein.csv")) %>%
  mutate(protein_per_kg = 10 * protein_per_100g) %>% pull(protein_per_kg) %>%
  rep(times = nrreg)
coeff_pro_vec[rep(items$comm_group, times = nrreg) == "Live animals"] <- 0

y_hh_t <- as(FABIO_y_hh, "TsparseMatrix")
is_import_flow <- producer_idx[y_hh_t@i + 1] != (y_hh_t@j + 1)
Y_imp <- Matrix::sparseMatrix(i = y_hh_t@i[is_import_flow] + 1, j = y_hh_t@j[is_import_flow] + 1,
                              x = y_hh_t@x[is_import_flow], dims = dim(FABIO_y_hh))
Y_imp_pro <- Matrix::Diagonal(x = coeff_pro_vec) %*% Y_imp * 1000   # g protein

pro_import_g <- setNames(colSums(Y_imp_pro), regions$iso3c)
hr_import_total_vec <- colSums(hr_total_country) - diag(hr_total_country)   # food + non-food, matches original plot
en_import_total_vec <- colSums(en_total_country) - diag(en_total_country)

# hr_animal_share/en_animal_share (df_composition, above) only exist at
# FOOD-item resolution -- fp_nonfood has already collapsed the FABIO-commodity
# axis via l_int_i %*% X, so non-food effort can't be attributed back to a
# specific food item's animal/plant origin. Given the food-vs-nonfood
# diagnostic above, break the conversion factor into food/non-food pieces too,
# so a weak correlation against the TOTAL can be told apart from "the
# classifier doesn't move even the food-sector piece it's derived from" vs.
# "the food-sector piece moves as expected, but non-food dominates the total
# and isn't explained by food composition at all".
df_conv <- df_composition %>%
  mutate(pro_import_g = pro_import_g[country],
         hr_import_total    = hr_import_total_vec[country],
         en_import_total    = en_import_total_vec[country],
         hr_import_food     = hr_food_import_vec[country],
         hr_import_nonfood  = hr_nonfood_import_vec[country],
         en_import_food     = en_food_import_vec[country],
         en_import_nonfood  = en_nonfood_import_vec[country],
         hr_per_50g_protein         = hr_import_total    * 1e6 / pro_import_g * 50,
         mj_per_50g_protein         = en_import_total     * 1e6 / pro_import_g * 50,
         hr_per_50g_protein_food    = hr_import_food      * 1e6 / pro_import_g * 50,
         hr_per_50g_protein_nonfood = hr_import_nonfood   * 1e6 / pro_import_g * 50,
         mj_per_50g_protein_food    = en_import_food       * 1e6 / pro_import_g * 50,
         mj_per_50g_protein_nonfood = en_import_nonfood    * 1e6 / pro_import_g * 50,
         is_row = country %in% row_countries) %>%
  filter(is.finite(hr_per_50g_protein), is.finite(mj_per_50g_protein), pro_import_g > 0, !is_row)

cat("\n---- H2 correlations: animal-origin share of import-effort vs. conversion factor, by sector ----\n")
cat(sprintf("  hr_animal_share vs TOTAL    hr/50g:  r = %.3f\n", cor(df_conv$hr_animal_share, df_conv$hr_per_50g_protein)))
cat(sprintf("  hr_animal_share vs FOOD     hr/50g:  r = %.3f\n", cor(df_conv$hr_animal_share, df_conv$hr_per_50g_protein_food)))
cat(sprintf("  hr_animal_share vs NONFOOD  hr/50g:  r = %.3f\n", cor(df_conv$hr_animal_share, df_conv$hr_per_50g_protein_nonfood)))
cat(sprintf("  en_animal_share vs TOTAL    MJ/50g:  r = %.3f\n", cor(df_conv$en_animal_share, df_conv$mj_per_50g_protein)))
cat(sprintf("  en_animal_share vs FOOD     MJ/50g:  r = %.3f\n", cor(df_conv$en_animal_share, df_conv$mj_per_50g_protein_food)))
cat(sprintf("  en_animal_share vs NONFOOD  MJ/50g:  r = %.3f\n", cor(df_conv$en_animal_share, df_conv$mj_per_50g_protein_nonfood)))

p_h2 <- ggplot(df_conv, aes(x = mj_per_50g_protein, y = hr_per_50g_protein)) +
  geom_point(aes(color = hr_animal_share, size = pro_import_g), alpha = 0.85) +
  geom_text_repel(aes(label = country), size = 3, max.overlaps = 15) +
  scale_color_viridis_c(option = "C", labels = scales::percent) +
  scale_size_continuous(range = c(1, 8), guide = "none") +
  labs(x = "Energy (MJ / 50 g imported-consumed protein)",
       y = "Time (hr / 50 g imported-consumed protein)",
       color = "Animal-origin\nshare of\nimport-effort hours",
       title = "H2 -- import-effort line colored by animal-vs-plant mix") +
  theme_minimal() + theme(legend.position = "right")
print(p_h2)
ggsave("results/h2_livestock_share_vs_line_position.pdf", p_h2, width = 10, height = 7)


#### H3: does an importer's OWN mechanization level explain where it sits on the line? ####
# H2's animal-origin share came out negatively correlated with hr/50g -- the
# WRONG sign for a feed-conversion story -- which smells like a confound with
# wealth: richer countries tend to both eat more meat AND run lower-labor-hour
# (more mechanized/capital-intensive) supply chains for unrelated reasons, and
# the latter, larger effect could be swamping/reversing the former. Proxy each
# country's own mechanization level with df_exporters$ratio -- the energy/labor
# ratio of what IT exports/produces (already computed above, self-contained,
# no new data) -- and test whether that, not who it trades with or what it
# eats, predicts its position on the import-effort line. If trade partners are
# systematically similarly-developed (gravity/homophily in trade networks), a
# country's own tech level should proxy for its import basket's average tech
# level too.
own_ratio <- setNames(df_exporters$ratio, df_exporters$country)   # own overall (food+non-food) energy/labor ratio

df_h3 <- df_conv %>%
  mutate(own_ratio = own_ratio[country]) %>%
  filter(is.finite(own_ratio))

cat("\n---- H3 correlations: own (exporter-side) energy/labor ratio vs. import-effort conversion factor ----\n")
cat(sprintf("  own_ratio vs hr/50g protein (import effort):  r = %.3f\n", cor(df_h3$own_ratio, df_h3$hr_per_50g_protein)))
cat(sprintf("  own_ratio vs MJ/50g protein (import effort):  r = %.3f\n", cor(df_h3$own_ratio, df_h3$mj_per_50g_protein)))
cat(sprintf("  own_ratio vs the LINE ITSELF (mj_per_50g/hr_per_50g, i.e. each importer's own import-side ratio): r = %.3f\n",
           cor(df_h3$own_ratio, df_h3$mj_per_50g_protein / df_h3$hr_per_50g_protein)))

p_h3 <- ggplot(df_h3, aes(x = mj_per_50g_protein, y = hr_per_50g_protein)) +
  geom_point(aes(color = own_ratio, size = pro_import_g), alpha = 0.85) +
  geom_text_repel(aes(label = country), size = 3, max.overlaps = 15) +
  scale_color_viridis_c(option = "C", trans = "log10") +
  scale_size_continuous(range = c(1, 8), guide = "none") +
  labs(x = "Energy (MJ / 50 g imported-consumed protein)",
       y = "Time (hr / 50 g imported-consumed protein)",
       color = "Own (exporter-side)\nenergy/labor\nratio (MJ/hr, log)",
       title = "H3 -- import-effort line colored by importer's OWN mechanization level") +
  theme_minimal() + theme(legend.position = "right")
print(p_h3)
ggsave("results/h3_own_mechanization_vs_line_position.pdf", p_h3, width = 10, height = 7)


#### Diagnostic: is any of the above actually a RoW-pasting artifact? ####
# All FABIO countries mapped to the same "Rest of World" EXIO aggregate share
# IDENTICAL per-tonne energy/labor intensities (single EXIO value pasted to
# every member -- see reorder_countries_to_FABIO() / CLAUDE.md, and the same
# caveat already flagged in 6.protein_import_effort_decomposition.R). We've
# excluded RoW-mapped countries from being analyzed AS importers/exporters
# throughout (is_row filters above), so the ~40 countries in every plot so far
# are all individually modeled -- but their import BASKETS still draw on the
# ~140 other FABIO countries that ARE RoW-pasted. If a country sources a large
# share of its import effort from RoW aggregates, part of its line position
# reflects just one of 5 pasted technology archetypes, not real supplier
# diversity -- which could manufacture spurious cross-country correlation
# independent of H1/H2/H3.
row_idx <- which(regions$iso3c %in% row_countries)

row_hr_share_vec <- setNames(colSums(hr_total_country[row_idx, ]) / hr_import_total_vec, regions$iso3c)
row_en_share_vec <- setNames(colSums(en_total_country[row_idx, ]) / en_import_total_vec, regions$iso3c)

df_row_share <- df_conv %>%
  transmute(country, hr_per_50g_protein, mj_per_50g_protein,
            row_hr_share = row_hr_share_vec[country],
            row_en_share = row_en_share_vec[country])

cat("\n---- RoW-aggregate share of each analyzed importer's import-effort ----\n")
cat(sprintf("  median RoW share of import-effort hours:  %.1f%% (range %.1f - %.1f%%)\n",
           median(df_row_share$row_hr_share) * 100, min(df_row_share$row_hr_share) * 100, max(df_row_share$row_hr_share) * 100))
cat(sprintf("  median RoW share of import-effort energy: %.1f%% (range %.1f - %.1f%%)\n",
           median(df_row_share$row_en_share) * 100, min(df_row_share$row_en_share) * 100, max(df_row_share$row_en_share) * 100))
cat(sprintf("\nCorrelation(RoW hr-share, hr/50g protein): r = %.3f\n", cor(df_row_share$row_hr_share, df_row_share$hr_per_50g_protein)))
cat(sprintf("Correlation(RoW en-share, MJ/50g protein): r = %.3f\n", cor(df_row_share$row_en_share, df_row_share$mj_per_50g_protein)))

p_row_check <- ggplot(df_row_share, aes(x = mj_per_50g_protein, y = hr_per_50g_protein)) +
  geom_point(aes(color = row_hr_share), size = 3) +
  geom_text_repel(aes(label = country), size = 3, max.overlaps = 15) +
  scale_color_viridis_c(option = "D", labels = scales::percent) +
  labs(x = "Energy (MJ / 50 g imported-consumed protein)",
       y = "Time (hr / 50 g imported-consumed protein)",
       color = "Share of import-effort\nhours from RoW-pasted\nexporters",
       title = "Diagnostic -- import-effort line colored by reliance on RoW-pasted exporters") +
  theme_minimal() + theme(legend.position = "right")
print(p_row_check)
ggsave("results/diagnostic_row_share_vs_line_position.pdf", p_row_check, width = 10, height = 7)
fwrite(df_row_share, "output/diagnostic_row_share_vs_conversion.csv")

# RoW share is huge AND roughly uniform across every importer (no variance to
# correlate against position) -- so every importer's basket is effectively a
# weighted average of the same handful of underlying EXIO-region technology
# archetypes (44 named regions + 5 RoW aggregates = 49 total). A weighted
# average of points that themselves sit near a line will land near that line
# for ANY weighting -- no story about mechanization, diet, or exporter
# identity required. Test the root cause directly: do the 49 EXIO-region
# archetypes themselves fall near a line in energy/labor space? (Using the
# UNFILTERED exported_hr/exported_en from the H1 section -- all 187 FABIO
# countries, RoW-mapped included -- aggregated up to their 49 EXIO regions.)
exio_region_of <- setNames(FABIO_reg$EXIOBASE, FABIO_reg$ISO)

df_region <- data.frame(country = regions$iso3c, exported_hr, exported_en) %>%
  mutate(exio_region = exio_region_of[country]) %>%
  filter(is.finite(exported_hr), is.finite(exported_en), exported_hr > 0, exported_en > 0) %>%
  group_by(exio_region) %>%
  summarise(exported_hr = sum(exported_hr), exported_en = sum(exported_en), n_countries = n(), .groups = "drop") %>%
  mutate(ratio = exported_en / exported_hr, is_row = grepl("RoW", exio_region))

cat("\n---- EXIO-region-level (49 regions) energy/labor archetypes ----\n")
print(df_region %>% arrange(desc(exported_hr)) %>% mutate(across(c(exported_hr, exported_en, ratio), ~round(.x, 3))))

p_region <- ggplot(df_region, aes(x = exported_hr, y = exported_en)) +
  geom_point(aes(color = is_row), size = 3) +
  geom_text_repel(aes(label = exio_region), size = 3, max.overlaps = 20) +
  geom_smooth(method = "lm", se = FALSE, color = "grey50", linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = c("FALSE" = "#1f77b4", "TRUE" = "#d62728")) +
  labs(x = "Total exported effort by EXIO region (M.hr, log)",
       y = "Total exported energy by EXIO region (TJ, log)",
       color = "RoW aggregate",
       title = "Diagnostic -- do the 49 underlying EXIO-region archetypes themselves fall near a line?") +
  theme_minimal()
print(p_region)
ggsave("results/diagnostic_exio_region_ratio.pdf", p_region, width = 9, height = 7)
fwrite(df_region, "output/diagnostic_exio_region_ratio.csv")

fwrite(top_exporters, "output/h1_top_exporters_ratio.csv")
fwrite(top_exporters_food, "output/h1b_top_exporters_ratio_food_only.csv")
fwrite(top_exporters_nonfood, "output/h1c_top_exporters_ratio_nonfood_only.csv")
fwrite(df_conv %>% select(country, hr_animal_share, en_animal_share, pro_import_g,
                          mj_per_50g_protein, hr_per_50g_protein,
                          mj_per_50g_protein_food, hr_per_50g_protein_food,
                          mj_per_50g_protein_nonfood, hr_per_50g_protein_nonfood),
      "output/h2_livestock_share_vs_conversion.csv")
