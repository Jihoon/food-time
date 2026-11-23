#### Domestic embodied energy by sector and commodity: high-energy group vs comparators ####
#
# Run in the session where 0.mrio_prep.R and 2.analyze_result.R are loaded (uses l_int_d,
# l_int_i, FABIO_L, FABIO_y_hh, regions, items, exio_nonfood_sectors). Writes CSVs to results/.
#
# "Domestic" = the diagonal convention of 2.analyze_result.R / 5.nonfood_energy_breakdown.R:
# items produced in country i that feed country i's own food consumption.
#   food-sector energy     = l_int_d$en (per tonne)            x embodied tonnes
#   non-food-sector energy = l_int_i$en (sector x item, per t)  x embodied tonnes
# Per-50 g figures use protein in i's own items reaching i's final food demand, which
# should be close to Figure 3's domestic-effort CF_energy (printed for comparison).

suppressMessages({ library(tidyverse); library(Matrix); library(data.table) })
need <- c("l_int_d", "l_int_i", "FABIO_L", "FABIO_y_hh", "regions", "items", "exio_nonfood_sectors", "FABIO_reg")
miss <- need[!sapply(need, exists)]
if (length(miss)) stop("Missing session objects: ", paste(miss, collapse = ", "))
if (!"en" %in% names(l_int_d) || !"en" %in% names(l_int_i)) stop("l_int_d / l_int_i have no 'en' element: ", paste(names(l_int_d), collapse = ", "))

FABIO_path <- "H:/MyDocuments/Data/FABIO/input/"
nrreg <- nrow(regions); nrcom <- nrow(items); n_nf <- length(exio_nonfood_sectors)
if (!exists("io")) io <- fread(file.path(FABIO_path, "io_labels.csv"))
prot_kg <- fread(file.path(FABIO_path, "nutrient_coefficients_protein.csv"))$protein_per_100g * 10   # g protein per kg, per item
prot_kg[items$comm_group == "Live animals"] <- 0
# l_int_i rows are EXIO-region resolved (49 regions x n_nf sectors, region-major), not per
# FABIO country -- see 1.1.mrio_convert_indirect.R. "Domestic" non-food effort for FABIO
# country i is therefore the row block of i's EXIO region (for a RoW member, its RoW region).
n_exio <- nrow(l_int_i$en) / n_nf
stopifnot(length(prot_kg) == nrcom, n_exio == round(n_exio), length(l_int_d$en) == nrreg * nrcom,
          identical(regions$iso3c, FABIO_reg$ISO))
cat("l_int_i: ", n_exio, "EXIO regions x", n_nf, "non-food sectors\n")
exio_code <- as.integer(FABIO_reg$EXIOBASE_code)                  # 1..49, EXIO region order
if (exists("EXIO_reg")) {                                          # sanity: block 41 should be Taiwan
  nf_reg <- EXIO_reg$region[EXIO_reg$sector %in% exio_nonfood_sectors]
  cat("block 41 region:", unique(nf_reg[(40 * n_nf + 1):(41 * n_nf)]), "| KOR block:",
      unique(nf_reg[((exio_code[regions$iso3c == "KOR"] - 1) * n_nf + 1):(exio_code[regions$iso3c == "KOR"] * n_nf)]), "\n")
}

one_country <- function(i) {
  rows <- ((i - 1) * nrcom + 1):(i * nrcom)
  X_i  <- as.vector(FABIO_L %*% FABIO_y_hh[, i])[rows]           # tonnes of i's own items feeding i
  y_i  <- as.vector(FABIO_y_hh[rows, i])                          # tonnes of i's own items in i's final demand
  food_by_item <- as.numeric(l_int_d$en[rows]) * X_i
  k <- exio_code[i]
  A <- l_int_i$en[((k - 1) * n_nf + 1):(k * n_nf), rows, drop = FALSE]
  nf <- as.matrix(A) * rep(X_i, each = n_nf)                      # sector x item
  list(food_by_item = food_by_item, nf_by_sector = rowSums(nf), nf_by_item = colSums(nf),
       protein_g = sum(y_i * prot_kg) * 1000)
}

breakdown <- function(label, isos) {
  idx <- which(regions$iso3c %in% isos)
  if (!length(idx)) { warning("no FABIO countries for ", label); return(NULL) }
  parts <- lapply(idx, one_country)
  s <- function(k) Reduce(`+`, lapply(parts, `[[`, k))
  food_item <- s("food_by_item"); nf_sector <- s("nf_by_sector"); nf_item <- s("nf_by_item"); P <- s("protein_g")
  total <- sum(food_item) + sum(nf_sector)
  per50 <- function(x) x / P * 50
  list(
    summary = tibble(region = label, n_countries = length(idx), total = total,
                     food_sector_share = sum(food_item) / total, nonfood_share = sum(nf_sector) / total,
                     energy_per_50g_raw = per50(total)),
    sectors = tibble(region = label, sector = exio_nonfood_sectors, energy = nf_sector) %>%
      mutate(share_of_total = energy / total) %>% arrange(desc(energy)) %>% slice_head(n = 8),
    items = tibble(region = label, item = items$item, energy = food_item + nf_item) %>%
      mutate(share_of_total = energy / total) %>% arrange(desc(energy)) %>% slice_head(n = 8))
}

mideast <- FABIO_reg$ISO[FABIO_reg$EXIOBASE == "RoW Middle East"]
targets <- list("South Korea" = "KOR", "Norway" = "NOR", "Greece" = "GRC", "RoW Middle East" = mideast,
                "China" = "CHN", "India" = "IND", "Netherlands" = "NLD", "United States" = "USA")
res <- Filter(Negate(is.null), Map(breakdown, names(targets), targets))

summ <- bind_rows(lapply(res, `[[`, "summary"))
# Units of l_int_*$en are not restated here; scale so the per-50 g value is comparable with
# Figure 3 (MJ). If energy_per_50g_raw is ~1e-6 of the Figure 3 value, the satellite is in TJ.
cat("\n---- Domestic embodied energy: totals, food vs non-food sector share, per 50 g (raw units) ----\n")
print(summ %>% mutate(across(c(food_sector_share, nonfood_share), ~ sprintf("%.0f%%", 100 * .x))), width = Inf)
if (exists("tradeoff_protein_totalecon_domestic_effort")) {
  cat("\nFigure 3 domestic-effort CF_energy (MJ / 50 g), food + non-food, for comparison:\n")
  print(tradeoff_protein_totalecon_domestic_effort %>%
          filter(effort_label == "Domestic effort", exio_region %in% names(targets)) %>%
          group_by(exio_region) %>% summarise(mj_per_50g = first(mj_per_50g_protein), .groups = "drop"))
}
cat("\n---- Top non-food EXIO sectors (share of each region's total domestic energy) ----\n")
print(bind_rows(lapply(res, `[[`, "sectors")) %>% mutate(share = sprintf("%.1f%%", 100 * share_of_total)) %>%
        select(region, sector, share), n = Inf, width = Inf)
cat("\n---- Top FABIO commodities (food + non-food sector energy, share of total) ----\n")
print(bind_rows(lapply(res, `[[`, "items")) %>% mutate(share = sprintf("%.1f%%", 100 * share_of_total)) %>%
        select(region, item, share), n = Inf, width = Inf)

out_dir <- "results"
write.csv(summ, file.path(out_dir, "energy_route_breakdown_summary.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(res, `[[`, "sectors")), file.path(out_dir, "energy_route_breakdown_sectors.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(res, `[[`, "items")), file.path(out_dir, "energy_route_breakdown_items.csv"), row.names = FALSE)
cat("\nWrote energy_route_breakdown_{summary,sectors,items}.csv to results/\n")
