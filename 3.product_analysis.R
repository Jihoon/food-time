#### 3. Product-level time inequality analysis ####
#
# Analyzes which food products contribute most to labor time inequalities
# across countries, covering:
#   (1) Production-side labor intensity and cross-country inequality by product
#   (2) Consumption-based time footprint attributed to each food product
#   (3) Decomposition of cross-country inequality by product contribution
#   (4) Gender (male/female) labor composition by product
#
# Run this script in the same R session after scripts 0 and 1, which load
# all necessary objects: l_int_d, FABIO_x, FABIO_L, FABIO_y_hh, io, items,
# regions, coeff_cal, coeff_pro, year, yr, FABIO_path, FABIO_reg.

library(tidyverse)
library(data.table)
library(Matrix)
library(ggplot2)
library(ggrepel)
library(gt)

output_dir <- "output"
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

# Load intensity vectors saved by script 1 (M.hr/tonne per product-country row)
l_int_d <- readRDS(paste0("data/FABIO_exio_satellites_food_", year, ".rds"))

# Nutrient coefficients (copied from script 1 for standalone use)
coeff_cal <- fread(file.path(FABIO_path, "nutrient_coefficients.csv")) %>%
  select(kcal_per_kg) %>%
  slice(rep(1:n(), times = nrow(FABIO_reg))) %>%
  mutate(kcal_per_kg = ifelse(
    items$comm_group[match(io$item, items$item)] == "Live animals", 0, kcal_per_kg))

# Production attributed to household food consumption (23001 × 187 mass matrix)
# Expensive: ~30-60 s; skip if already computed in session
if (!exists("FABIO_x_hh")) {
  message("Computing FABIO_x_hh = FABIO_L %*% FABIO_y_hh ...")
  FABIO_x_hh <- FABIO_L %*% FABIO_y_hh
}

# Population (aligned to regions$iso3c)
data(countrypops)
pop_y <- subset(countrypops, year == yr) %>%
  select(iso3c = country_code_3, pop = population) %>%
  right_join(regions %>% select(iso3c), by = "iso3c")
pop_vec <- pop_y$pop  # length 187, matches column order of FABIO matrices

# Countries excluded due to missing population
no_pop_ctys <- c("ROW", "TWN", "ANT")


# ── Helper ──────────────────────────────────────────────────────────────────

# Standard Gini coefficient (size-neutral, bounded [0,1])
gini_coeff <- function(x) {
  x <- sort(x[x > 0 & is.finite(x)])
  if (length(x) < 2) return(NA_real_)
  n <- length(x)
  (2 * sum(seq_along(x) * x) / (n * sum(x))) - (n + 1) / n
}


# ============================================================================
# 1. PRODUCTION-SIDE LABOR INTENSITY BY PRODUCT
#    Measures: for each food product, how labor-intensive is production
#    and how unequal is that intensity across producing countries?
# ============================================================================

# Flat data frame: one row per country × product (23001 rows)
df_int <- data.frame(
  iso3c       = io$iso3c,
  item        = io$item,
  hr_m        = l_int_d$hr_m,   # M.hr per tonne (male)
  hr_f        = l_int_d$hr_f,   # M.hr per tonne (female)
  en          = l_int_d$en,     # TJ per tonne
  production  = FABIO_x,        # tonnes
  kcal_per_kg = coeff_cal$kcal_per_kg
) %>%
  left_join(items %>% select(item, comm_group), by = "item") %>%
  left_join(regions %>% select(iso3c, continent), by = "iso3c") %>%
  filter(comm_group != "Live animals") %>%
  mutate(
    hr_total     = hr_m + hr_f,
    female_share = if_else(hr_total > 0, hr_f / hr_total, NA_real_),
    # Labor per Mcal: (M.hr/tonne) / (Mcal/tonne) = M.hr/Mcal
    # 1 tonne has kcal_per_kg × 1000 kcal = kcal_per_kg × 1e-3 Mcal
    hr_m_per_Mcal    = if_else(kcal_per_kg > 0, hr_m * 1e3 / kcal_per_kg, NA_real_),
    hr_f_per_Mcal    = if_else(kcal_per_kg > 0, hr_f * 1e3 / kcal_per_kg, NA_real_),
    hr_total_per_Mcal = if_else(kcal_per_kg > 0, hr_total * 1e3 / kcal_per_kg, NA_real_)
  )

# Production-weighted product summary
prod_summary <- df_int %>%
  filter(production > 0) %>%
  group_by(item, comm_group) %>%
  summarise(
    n_countries        = n(),
    total_prod_Mt      = sum(production) / 1e6,
    # Global totals (M.hr)
    total_hr_m         = sum(hr_m * production),
    total_hr_f         = sum(hr_f * production),
    total_hr           = total_hr_m + total_hr_f,
    # Production-weighted average intensity
    wtd_hr_per_t       = weighted.mean(hr_total, production),
    wtd_hr_per_Mcal    = weighted.mean(hr_total_per_Mcal, production, na.rm = TRUE),
    wtd_female_share   = weighted.mean(female_share, production, na.rm = TRUE),
    # Cross-country inequality in labor intensity (Gini of hr/tonne)
    gini_hr_intensity  = gini_coeff(hr_total),
    cv_hr_intensity    = sd(hr_total) / mean(hr_total),
    .groups = "drop"
  ) %>%
  mutate(
    global_hr_share_pct = total_hr / sum(total_hr, na.rm = TRUE) * 100
  ) %>%
  arrange(desc(total_hr))

fwrite(prod_summary,
       file.path(output_dir, paste0("product_labor_summary_", year, ".csv")))

message("Top 10 products by global labor footprint:")
print(prod_summary %>% select(item, comm_group, total_hr, global_hr_share_pct, gini_hr_intensity) %>% head(10))


# ============================================================================
# 2. CONSUMPTION-BASED TIME FOOTPRINT BY PRODUCT
#    For each food product, how much labor time is embedded in the quantities
#    consumed across countries?  Inequality measured as Gini of per-capita
#    footprint across countries.
# ============================================================================

# Footprint matrices: fp[i,j] = M.hr in product-origin row i for country j's consumption
fp_hr_m <- Diagonal(x = l_int_d$hr_m) %*% FABIO_x_hh   # 23001 × 187 (M.hr)
fp_hr_f <- Diagonal(x = l_int_d$hr_f) %*% FABIO_x_hh
fp_cal  <- Diagonal(x = coeff_cal$kcal_per_kg * 1e3) %*% FABIO_x_hh  # kcal (×1000 kg/t)

# Aggregate 23001 product-origin rows → 123 product rows (sum over origins)
fp_hr_m_item <- rowsum(as.matrix(fp_hr_m), group = io$item)   # 123 × 187
fp_hr_f_item <- rowsum(as.matrix(fp_hr_f), group = io$item)
fp_cal_item  <- rowsum(as.matrix(fp_cal),  group = io$item)   # kcal

colnames(fp_hr_m_item) <- colnames(fp_hr_f_item) <-
  colnames(fp_cal_item) <- regions$iso3c

# Per-capita footprint:
#   M.hr / person × 1e6 = hr/person  (consistent with script 2 convention)
fp_hr_m_pcap <- t(t(fp_hr_m_item) / pop_vec) * 1e6   # hr/person
fp_hr_f_pcap <- t(t(fp_hr_f_item) / pop_vec) * 1e6
fp_cal_pcap  <- t(t(fp_cal_item)  / pop_vec) / 365   # kcal/person/day

fp_hr_total_pcap <- fp_hr_m_pcap + fp_hr_f_pcap       # hr/person

# Country-level totals (sum over all products)
total_hr_pcap  <- colSums(fp_hr_total_pcap, na.rm = TRUE)
total_cal_pcap <- colSums(fp_cal_pcap, na.rm = TRUE)

valid_ctys <- regions$iso3c[!regions$iso3c %in% no_pop_ctys]

# Per-product consumption summary
cons_summary <- lapply(rownames(fp_hr_m_item), function(prod) {
  hr_t <- fp_hr_total_pcap[prod, valid_ctys]
  cal  <- fp_cal_pcap[prod, valid_ctys]
  hr_t[!is.finite(hr_t)] <- 0
  cal[!is.finite(cal)]   <- 0

  data.frame(
    item                 = prod,
    total_cons_hr_Mhr    = sum(fp_hr_m_item[prod, ] + fp_hr_f_item[prod, ]),
    total_cons_cal_Mcal  = sum(fp_cal_item[prod, ]) / 1e6,
    gini_hr_pcap         = gini_coeff(hr_t),   # inequality of time footprint across countries
    gini_cal_pcap        = gini_coeff(cal),
    # Spearman: do countries with more calories also embed more labor?
    spearman_cal_hr      = if (sum(hr_t > 0) > 5)
                             cor(cal, hr_t, method = "spearman") else NA_real_
  )
}) %>%
  bind_rows() %>%
  left_join(items %>% select(item, comm_group), by = "item") %>%
  filter(comm_group != "Live animals") %>%
  arrange(desc(total_cons_hr_Mhr))

fwrite(cons_summary,
       file.path(output_dir, paste0("product_consumption_footprint_", year, ".csv")))


# ============================================================================
# 3. DECOMPOSITION: WHICH PRODUCTS DRIVE CROSS-COUNTRY TIME INEQUALITY?
#    Uses Shorrocks covariance decomposition:
#      Contribution_k = Cov(h_k, H) / Var(H) × 100%
#    where H_c = total hr/cap in country c, h_{k,c} = hr/cap from product k.
#    Positive contribution → product amplifies inequality.
#    Negative contribution → product dampens inequality.
# ============================================================================

H_vec <- total_hr_pcap[valid_ctys]
H_vec[!is.finite(H_vec)] <- 0
H_valid <- H_vec[H_vec > 0]
var_H   <- var(H_valid)

decomp_df <- lapply(rownames(fp_hr_total_pcap), function(prod) {
  h_k <- fp_hr_total_pcap[prod, valid_ctys]
  h_k[!is.finite(h_k)] <- 0
  # Align to countries with positive total footprint
  h_k_valid <- h_k[names(H_valid)]

  data.frame(
    item              = prod,
    mean_hr_pcap      = mean(h_k_valid),
    cov_with_total    = cov(h_k_valid, H_valid),
    share_of_var_pct  = cov(h_k_valid, H_valid) / var_H * 100,
    total_cons_Mhr    = sum(fp_hr_m_item[prod, ] + fp_hr_f_item[prod, ])
  )
}) %>%
  bind_rows() %>%
  left_join(items %>% select(item, comm_group), by = "item") %>%
  filter(comm_group != "Live animals") %>%
  arrange(desc(share_of_var_pct))

fwrite(decomp_df,
       file.path(output_dir, paste0("product_inequality_decomp_", year, ".csv")))

message("Top 10 products driving cross-country time inequality:")
print(decomp_df %>% select(item, comm_group, share_of_var_pct, total_cons_Mhr) %>% head(10))


# ============================================================================
# 4. GENDER LABOR COMPOSITION BY PRODUCT
#    Which products have the highest female labor share, and how does that
#    vary across continents?
# ============================================================================

gender_by_product <- df_int %>%
  filter(production > 0, is.finite(female_share)) %>%
  group_by(item, comm_group) %>%
  summarise(
    global_female_share = weighted.mean(female_share, production, na.rm = TRUE),
    total_hr_f          = sum(hr_f * production),
    total_hr_m          = sum(hr_m * production),
    total_hr            = total_hr_f + total_hr_m,
    .groups = "drop"
  ) %>%
  arrange(desc(global_female_share))

gender_by_continent <- df_int %>%
  filter(production > 0, is.finite(female_share)) %>%
  group_by(item, comm_group, continent) %>%
  summarise(
    wtd_female_share = weighted.mean(female_share, production, na.rm = TRUE),
    total_hr         = sum(hr_total * production),
    .groups = "drop"
  )

fwrite(gender_by_product,
       file.path(output_dir, paste0("product_gender_summary_", year, ".csv")))


# ============================================================================
# 5. PLOTS
# ============================================================================

# ── 5a. Top 25 products by total global labor footprint ─────────────────────

top25_hr <- prod_summary %>%
  slice_head(n = 25) %>%
  mutate(item = factor(item, levels = rev(item))) %>%
  pivot_longer(c(total_hr_m, total_hr_f), names_to = "gender", values_to = "hr") %>%
  mutate(gender = recode(gender, "total_hr_m" = "Male", "total_hr_f" = "Female"))

ggplot(top25_hr, aes(x = item, y = hr, fill = gender)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Male" = "#4d96d4", "Female" = "#d4704d")) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = NULL, y = "Total labor hours (M.hr)", fill = NULL,
       title = paste0("Top 25 food products by global labor time footprint (", year, ")"),
       subtitle = "Production-side; stacked by gender") +
  theme_minimal() +
  theme(legend.position = "top", axis.text.y = element_text(size = 9))

ggsave(file.path(output_dir, paste0("plots/product_top25_labor_", year, ".png")),
       width = 9, height = 8)


# ── 5b. Total footprint vs. cross-country inequality (Gini) ─────────────────

label_threshold <- prod_summary %>%
  filter(global_hr_share_pct > 1.5 | gini_hr_intensity > 0.65)

ggplot(prod_summary,
       aes(x = global_hr_share_pct, y = gini_hr_intensity,
           size = total_prod_Mt, color = comm_group)) +
  geom_point(alpha = 0.75) +
  geom_text_repel(data = label_threshold, aes(label = item),
                  size = 3, max.overlaps = 20, show.legend = FALSE) +
  scale_size_continuous(range = c(1.5, 12), name = "Production (Mt)") +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Share of global labor footprint (%)",
       y = "Gini of labor intensity across countries",
       color = "Commodity group",
       title = paste0("Food product labor footprint vs. cross-country inequality (", year, ")"),
       subtitle = "Size = production volume; Gini measures how unequally labor-intensive this product is across countries") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(file.path(output_dir, paste0("plots/product_footprint_vs_gini_", year, ".png")),
       width = 11, height = 7)


# ── 5c. Inequality decomposition: contribution to cross-country variance ─────

decomp_plot <- decomp_df %>%
  filter(abs(share_of_var_pct) > 0.5) %>%
  arrange(desc(share_of_var_pct)) %>%
  slice_head(n = 30) %>%
  mutate(item = factor(item, levels = rev(item)),
         direction = if_else(share_of_var_pct > 0, "Amplifies", "Dampens"))

ggplot(decomp_plot, aes(x = item, y = share_of_var_pct, fill = comm_group)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "Contribution to cross-country inequality (%)",
       fill = "Commodity group",
       title = paste0("Food products driving cross-country labor time inequality (", year, ")"),
       subtitle = "Shorrocks decomposition: Cov(product footprint, total footprint) / Var(total footprint)") +
  theme_minimal() +
  theme(legend.position = "top", axis.text.y = element_text(size = 9))

ggsave(file.path(output_dir, paste0("plots/product_inequality_decomp_", year, ".png")),
       width = 10, height = 8)


# ── 5d. Labor intensity per Mcal: most time-intensive products ───────────────

prod_summary %>%
  filter(!is.na(wtd_hr_per_Mcal), wtd_hr_per_Mcal > 0, total_prod_Mt > 0.1) %>%
  arrange(desc(wtd_hr_per_Mcal)) %>%
  slice_head(n = 25) %>%
  mutate(item = factor(item, levels = rev(item))) %>%
  ggplot(aes(x = item, y = wtd_hr_per_Mcal, fill = comm_group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = NULL, y = "Labor intensity (M.hr / Mcal)",
       fill = "Commodity group",
       title = paste0("Most time-intensive food products per unit energy (", year, ")"),
       subtitle = "Production-weighted average across countries") +
  theme_minimal() +
  theme(legend.position = "top", axis.text.y = element_text(size = 9))

ggsave(file.path(output_dir, paste0("plots/product_intensity_per_Mcal_", year, ".png")),
       width = 9, height = 8)


# ── 5e. Gender composition: female share by product (top 30 by total hours) ──

gender_plot_data <- gender_by_product %>%
  filter(total_hr > 0) %>%
  arrange(desc(total_hr)) %>%
  slice_head(n = 30) %>%
  mutate(item = factor(item, levels = item[order(global_female_share)])) %>%
  pivot_longer(c(total_hr_m, total_hr_f), names_to = "gender", values_to = "hr") %>%
  mutate(gender = recode(gender, "total_hr_m" = "Male", "total_hr_f" = "Female"))

ggplot(gender_plot_data, aes(x = item, y = hr, fill = gender)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "white", linewidth = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("Male" = "#4d96d4", "Female" = "#d4704d")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Share of labor hours", fill = NULL,
       title = paste0("Gender composition of labor by food product (", year, ")"),
       subtitle = "Top 30 products by total labor hours; ordered by female share") +
  theme_minimal() +
  theme(legend.position = "top", axis.text.y = element_text(size = 9))

ggsave(file.path(output_dir, paste0("plots/product_gender_share_", year, ".png")),
       width = 9, height = 8)


# ── 5f. Gini of consumption footprint vs. Gini of calorie availability ───────
# Shows whether products with unequal time footprints also have unequal calorie access

merged_gini <- cons_summary %>%
  filter(!is.na(gini_hr_pcap), !is.na(gini_cal_pcap)) %>%
  left_join(prod_summary %>% select(item, global_hr_share_pct, total_prod_Mt), by = "item")

label_thr2 <- merged_gini %>%
  filter(gini_hr_pcap > 0.65 | global_hr_share_pct > 1.5)

ggplot(merged_gini, aes(x = gini_cal_pcap, y = gini_hr_pcap,
                         size = total_cons_cal_Mcal, color = comm_group)) +
  geom_point(alpha = 0.75) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_text_repel(data = label_thr2, aes(label = item),
                  size = 3, max.overlaps = 20, show.legend = FALSE) +
  scale_size_continuous(range = c(1.5, 10), name = "Consumption (Mcal)") +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Gini of calorie availability across countries",
       y = "Gini of labor time footprint across countries",
       color = "Commodity group",
       title = paste0("Calorie inequality vs. labor time inequality by food product (", year, ")"),
       subtitle = "Points above the diagonal → labor more unequally distributed than calories") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(file.path(output_dir, paste0("plots/product_cal_vs_labor_gini_", year, ".png")),
       width = 11, height = 7)

message("Product analysis complete. Outputs saved to '", output_dir, "/'.")
