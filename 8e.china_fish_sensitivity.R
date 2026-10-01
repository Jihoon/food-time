#### Does China's Figure 3 position rest on fishing labor? ####
#
# Figure 3 places China far above the time-energy frontier (much more paid time per 50 g of
# domestically produced protein than countries at similar energy), and the Findings read
# that as farming that stopped short of labor-displacing mechanization. About a third of
# China's domestic food hours are fish (9z). EXIOBASE's numbers are taken as given; this
# only asks how much of China's position rests on fish.
#
# For every EXIO-modeled country: domestic paid hours (by gender) and energy embodied in the
# domestically produced and consumed stream (X_dom, own-country rows, direct food-sector
# intensities -- the left panel of Figure 3), per 50 g of the protein in that stream, with
# and without the fish item (its hours, energy and protein all removed). China's position is
# summarized as its residual from a cross-country fit of log time on log energy.
#
# Run in the pipeline session (after 0, 1, 2.analyze_result.R): needs X_dom (or Y_dom),
# FABIO_L, l_int_d, FABIO_y_hh_pro, io, regions, nrcom, row_countries.

library(tidyverse)
library(Matrix)

PROJ <- "h:/MyDocuments/Projects/REMASS/food time"
setwd(PROJ)

need <- c("l_int_d", "FABIO_y_hh_pro", "io", "regions", "nrcom")
miss <- need[!vapply(need, exists, logical(1))]
if (length(miss)) stop("Run 0, 1 and 2.analyze_result.R in this session first; missing: ", paste(miss, collapse = ", "))
if (!exists("X_dom")) {
  if (!exists("Y_dom")) stop("Y_dom / X_dom missing: run 2.analyze_result.R (section 1.3).")
  X_dom <- FABIO_L %*% Y_dom
}
cty <- regions$iso3c
f <- which(io$item[1:nrcom] == "Fish, Seafood"); stopifnot(length(f) == 1)
targets <- which(!cty %in% (if (exists("row_countries")) row_countries else character(0)))

pos <- map_dfr(targets, function(j) {
  rows <- ((j - 1) * nrcom + 1):(j * nrcom)
  t    <- as.numeric(X_dom[rows, j])                                   # tonnes produced in j for j's own consumption
  pro  <- as.numeric(FABIO_y_hh_pro[rows, j])                          # protein of j's domestic final food (g/yr)
  hm <- l_int_d$hr_m[rows] * 1e6 * t; hf <- l_int_d$hr_f[rows] * 1e6 * t; en <- l_int_d$en[rows] * 1e6 * t
  keep <- seq_len(nrcom) != f
  tibble(country = cty[j],
         fish_share_hours   = (hm[f] + hf[f]) / sum(hm + hf),
         fish_share_protein = pro[f] / sum(pro),
         hr_f_all  = sum(hf) / sum(pro) * 50,        hr_m_all  = sum(hm) / sum(pro) * 50,
         mj_all    = sum(en) / sum(pro) * 50,
         hr_f_nofish = sum(hf[keep]) / sum(pro[keep]) * 50, hr_m_nofish = sum(hm[keep]) / sum(pro[keep]) * 50,
         mj_nofish   = sum(en[keep]) / sum(pro[keep]) * 50)
}) %>% filter(is.finite(hr_f_all), hr_f_all > 0, mj_all > 0)

frontier <- function(time, mj) {                    # residual from log time ~ log energy across countries
  fit <- lm(log(time) ~ log(mj)); resid(fit)
}
pos <- pos %>%
  mutate(time_all = hr_f_all + hr_m_all, time_nofish = hr_f_nofish + hr_m_nofish,
         resid_all = frontier(time_all, mj_all), resid_nofish = frontier(time_nofish, mj_nofish),
         rank_time_all = rank(-time_all), rank_time_nofish = rank(-time_nofish),
         rank_resid_all = rank(-resid_all), rank_resid_nofish = rank(-resid_nofish))

cat("\n---- China and the other high-fish countries: Figure 3 (left panel) position with and without fish ----\n")
show <- pos %>% arrange(desc(fish_share_hours)) %>% slice_head(n = 8) %>% pull(country) %>% union("CHN")
print(pos %>% filter(country %in% show) %>% arrange(desc(fish_share_hours)) %>%
        transmute(country, `fish % of hours` = round(100 * fish_share_hours), `fish % of protein` = round(100 * fish_share_protein),
                  `hr/50g all` = round(time_all, 3), `no fish` = round(time_nofish, 3),
                  `MJ/50g all` = round(mj_all, 1), `no fish ` = round(mj_nofish, 1),
                  `above frontier (log) all` = round(resid_all, 2), `no fish  ` = round(resid_nofish, 2),
                  `rank by time all/no fish` = paste(rank_time_all, rank_time_nofish, sep = "/"),
                  `rank above frontier all/no fish` = paste(rank_resid_all, rank_resid_nofish, sep = "/")),
      n = Inf, width = Inf)
cat(sprintf("\nn = %d countries. Ranks: 1 = most time per 50 g / furthest above the frontier.\n", nrow(pos)))
cat("Read: if China stays far above the frontier without fish, the China passage stands on farming;",
    "\n      if it drops toward the cloud, part of 'China stopped short' is fishing labor.\n")
write_csv(pos, "results/china_fish_sensitivity.csv")
cat("Wrote results/china_fish_sensitivity.csv\n")
