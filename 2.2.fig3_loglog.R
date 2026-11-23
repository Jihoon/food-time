#### 2.2.fig3_loglog.R -- Figure 3 on log-log axes ####
#
# Same data, panels, colours and labels as Figure 3 in 2.1.figures.R; only the
# axes change, from linear to log10. On linear axes the rich-country cluster is
# squeezed against zero while India and the RoW regions stretch the scale; on
# log axes the cluster spreads out and equal distances mean equal ratios.
#
# Source AFTER 2.analyze_result.R, in the same R session: the four
# `*_domestic_effort` tables are built there. Nothing is recomputed here, and
# 2.1.figures.R is left untouched (its Figure 3 stays linear).
#
# Output (new files only, no existing figure is overwritten):
#   results/manuscript/fig3_tradeoff_domestic_effort_loglog.pdf
#   results/manuscript/fig3_tradeoff_domestic_effort_loglog.png   (for slides)

library(tidyverse)
library(patchwork)

required_objects = c(
  "year",
  "tradeoff_protein_foodecon_domestic_effort",
  "label_tradeoff_protein_foodecon_domestic_effort",
  "tradeoff_protein_totalecon_domestic_effort",
  "label_tradeoff_protein_totalecon_domestic_effort"
)
missing_objects = required_objects[!sapply(required_objects, exists)]
if (length(missing_objects) > 0) {
  stop("Source 2.analyze_result.R first -- missing objects: ",
       paste(missing_objects, collapse = ", "))
}

# Reuse the manuscript theme when 2.1.figures.R has already been sourced;
# otherwise define the same one so this script also runs on its own.
if (!exists("FIG_BASE")) FIG_BASE = 20
if (!exists("theme_fig")) {
  theme_fig = function(base = FIG_BASE, legend_position = "top", angle_x = 0) {
    theme_minimal(base_size = base) +
      theme(
        plot.title      = element_text(size = base * 1.30, face = "bold", hjust = 0.5,
                                       margin = margin(b = 14)),
        axis.title      = element_text(size = base * 1.20),
        axis.text       = element_text(size = base * 0.85),
        axis.text.x     = element_text(size = base * 0.85, angle = angle_x,
                                       hjust = if (angle_x == 0) 0.5 else 1),
        legend.position = legend_position,
        legend.title    = element_text(size = base * 1.05),
        legend.text     = element_text(size = base * 0.95),
        legend.key.size = unit(1.6, "lines"),
        strip.text      = element_text(size = base * 1.10, face = "bold", color = "black")
      )
  }
}

fig_dir = "results/manuscript"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Same subset as Figure 3: domestic effort only, Taiwan dropped (its non-food
# intensity is an outlier; see 9r.diagnose_twn_nonfood_outlier.R).
fig3_drop_regions = c("Taiwan")
fig3_subset = function(df) {
  df %>% filter(effort_label == "Domestic effort", !exio_region %in% fig3_drop_regions)
}

# A log axis cannot show zero or negative values. Both axes should be strictly
# positive here; if a region has a zero it is dropped with a message instead of
# a silent ggplot warning, so the panel count stays honest.
fig3_positive = function(df, what) {
  keep = is.finite(df$mj_per_50g_protein) & df$mj_per_50g_protein > 0 &
         is.finite(df$hr_per_50g_protein) & df$hr_per_50g_protein > 0
  if (any(!keep)) {
    message(sprintf("fig3 log-log (%s): dropped %d non-positive row(s): %s", what, sum(!keep),
                    paste(unique(df$exio_region[!keep]), collapse = ", ")))
  }
  df[keep, ]
}

fig3_food        = fig3_subset(tradeoff_protein_foodecon_domestic_effort)        %>% fig3_positive("food")
fig3_total       = fig3_subset(tradeoff_protein_totalecon_domestic_effort)       %>% fig3_positive("food + non-food")
fig3_food_label  = fig3_subset(label_tradeoff_protein_foodecon_domestic_effort)  %>% filter(exio_region %in% fig3_food$exio_region)
fig3_total_label = fig3_subset(label_tradeoff_protein_totalecon_domestic_effort) %>% filter(exio_region %in% fig3_total$exio_region)

fig3_size_limits = range(c(fig3_food$g_protein_per_cap_day,
                           fig3_total$g_protein_per_cap_day), na.rm = TRUE)

tradeoff_panel_log = function(df, labels, title) {
  ggplot(df, aes(x = mj_per_50g_protein, y = hr_per_50g_protein)) +
    geom_line(aes(group = exio_region), color = "grey60", linewidth = 0.6) +
    geom_point(aes(color = type, size = g_protein_per_cap_day), alpha = 0.85) +
    ggrepel::geom_text_repel(data = labels, aes(label = exio_region),
                             size = 5, max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = c(hr_f = "#ca2323", hr_m = "#1f77b4"),
                       labels = c(hr_f = "Female", hr_m = "Male")) +
    scale_size_continuous(range = c(1, 8), limits = fig3_size_limits) +
    # Expansion is applied on the log scale, so the same 5% / 12% margins keep
    # the largest points (8 mm radius) clear of the panel border as in Figure 3.
    scale_x_log10(expand = expansion(mult = c(0.05, 0.12)),
                  breaks = scales::breaks_log(n = 6),
                  labels = scales::label_number(drop0trailing = TRUE)) +
    scale_y_log10(expand = expansion(mult = c(0.05, 0.12)),
                  breaks = scales::breaks_log(n = 6),
                  labels = scales::label_number(drop0trailing = TRUE)) +
    annotation_logticks(sides = "bl", color = "grey50", linewidth = 0.3,
                        short = unit(0.08, "cm"), mid = unit(0.14, "cm"),
                        long = unit(0.2, "cm")) +
    labs(x = "Energy used in production (MJ / 50 g protein, log scale)",
         y = "Paid labor time (hr / 50 g protein, log scale)",
         color = "Gender", size = "g protein/cap/day", title = title) +
    theme_fig(legend_position = "right")
}

p_fig3_food_log  = tradeoff_panel_log(fig3_food,  fig3_food_label,  "Food sector")
p_fig3_total_log = tradeoff_panel_log(fig3_total, fig3_total_label, "Food + non-food sectors")

fig3_loglog = (p_fig3_food_log | p_fig3_total_log) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Paid labor time vs. energy per 50 g of domestically-consumed protein (", year,
                   ") — domestic effort, log-log"),
    caption = "Time is paid labor in production; energy is energy used in production. Household (unpaid) time is not included.",
    theme = theme(plot.title = element_text(size = FIG_BASE * 1.5, face = "bold", hjust = 0.5,
                                            margin = margin(b = 16)))) &
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "fig3_tradeoff_domestic_effort_loglog.pdf"),
       fig3_loglog, width = 22, height = 10)
ggsave(file.path(fig_dir, "fig3_tradeoff_domestic_effort_loglog.png"),
       fig3_loglog, width = 22, height = 10, dpi = 200, bg = "white")

cat("Saved log-log Figure 3 to", fig_dir, "\n")


#### Variant: female and male in separate panels ####

# Rows = gender, columns = sector. The connecting line between a region's female
# and male points is dropped (each panel holds one gender), so the two rows are
# read against each other instead. Axis limits are shared down each column, so a
# point sits at the same x in both rows and the same y scale applies to both
# genders.
gender_names = c(hr_f = "Female", hr_m = "Male")

fig3_gender_panel = function(df, g, title, xlim, ylim, show_x) {
  d = df %>% filter(type == g)
  ggplot(d, aes(x = mj_per_50g_protein, y = hr_per_50g_protein)) +
    geom_point(aes(color = type, size = g_protein_per_cap_day), alpha = 0.85) +
    ggrepel::geom_text_repel(aes(label = exio_region), size = 4.5, max.overlaps = 20,
                             show.legend = FALSE) +
    scale_color_manual(values = c(hr_f = "#ca2323", hr_m = "#1f77b4"),
                       labels = gender_names, guide = "none") +
    scale_size_continuous(range = c(1, 8), limits = fig3_size_limits) +
    scale_x_log10(limits = xlim, expand = expansion(mult = c(0.05, 0.12)),
                  breaks = scales::breaks_log(n = 6),
                  labels = scales::label_number(drop0trailing = TRUE)) +
    scale_y_log10(limits = ylim, expand = expansion(mult = c(0.05, 0.12)),
                  breaks = scales::breaks_log(n = 6),
                  labels = scales::label_number(drop0trailing = TRUE)) +
    annotation_logticks(sides = "bl", color = "grey50", linewidth = 0.3,
                        short = unit(0.08, "cm"), mid = unit(0.14, "cm"),
                        long = unit(0.2, "cm")) +
    labs(x = if (show_x) "Energy used in production (MJ / 50 g protein, log scale)" else NULL,
         y = paste0(gender_names[[g]], ": paid labor time (hr / 50 g protein, log scale)"),
         size = "g protein/cap/day", title = title) +
    theme_fig(legend_position = "right")
}

fig3_col = function(df, title) {
  xlim = range(df$mj_per_50g_protein)
  ylim = range(df$hr_per_50g_protein)
  top = fig3_gender_panel(df, "hr_f", title, xlim, ylim, show_x = FALSE) +
    theme(axis.title.x = element_blank())
  bottom = fig3_gender_panel(df, "hr_m", NULL, xlim, ylim, show_x = TRUE)
  top / bottom
}

fig3_loglog_gender = (fig3_col(fig3_food, "Food sector") | fig3_col(fig3_total, "Food + non-food sectors")) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Paid labor time vs. energy per 50 g of domestically-consumed protein (", year,
                   ") — domestic effort, log-log, by gender"),
    caption = "Time is paid labor in production; energy is energy used in production. Household (unpaid) time is not included.",
    theme = theme(plot.title = element_text(size = FIG_BASE * 1.5, face = "bold", hjust = 0.5,
                                            margin = margin(b = 16)))) &
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "fig3_tradeoff_domestic_effort_loglog_by_gender.pdf"),
       fig3_loglog_gender, width = 22, height = 16)
ggsave(file.path(fig_dir, "fig3_tradeoff_domestic_effort_loglog_by_gender.png"),
       fig3_loglog_gender, width = 22, height = 16, dpi = 150, bg = "white")

cat("Saved log-log by-gender Figure 3 to", fig_dir, "\n")


#### Variant: selective labels (frontier, outliers, Figure 2 countries) ####

# With ~50 regions a label on every point hides the pattern. This version labels
# three groups only, once per region (at the midpoint of its female and male
# points):
#   frontier -- no other region has both lower energy and lower time per 50 g
#               protein (bold)
#   outlier  -- the n_outlier regions farthest from the centre of the cloud in
#               standardised log-log space, frontier regions excluded (plain)
#   Figure 2 -- the eight countries in the mosaics, mosaic_isos in
#               2.analyze_result.R (italic)
# A region in several groups takes the first in that order. Every region still
# gets its points and its female-male line.
fig2_isos = if (exists("mosaic_isos")) mosaic_isos else
  c("USA", "FRA", "AUS", "CHN", "IDN", "IND", "BRA", "ZAF")
region_iso = function(r) {
  if (exists("region_to_iso_1to1")) unname(region_to_iso_1to1[as.character(r)])
  else rep(NA_character_, length(r))
}

fig3_select_labels = function(df, n_outlier = 4) {
  reg = df %>%
    group_by(exio_region) %>%
    summarise(mj_per_50g_protein = exp(mean(log(mj_per_50g_protein))),
              hr_per_50g_protein = exp(mean(log(hr_per_50g_protein))),
              .groups = "drop")
  x = reg$mj_per_50g_protein; y = reg$hr_per_50g_protein
  frontier = vapply(seq_along(x), function(i) !any(x < x[i] & y < y[i]), logical(1))
  dist = sqrt(as.numeric(scale(log(x)))^2 + as.numeric(scale(log(y)))^2)
  cand = which(!frontier)
  outlier = seq_along(x) %in% cand[order(-dist[cand])][seq_len(min(n_outlier, length(cand)))]
  fig2 = region_iso(reg$exio_region) %in% fig2_isos
  reg %>%
    mutate(group = case_when(frontier ~ "frontier", outlier ~ "outlier", fig2 ~ "figure2", TRUE ~ NA_character_),
           face  = c(frontier = "bold", outlier = "plain", figure2 = "italic")[group]) %>%
    filter(!is.na(group))
}

tradeoff_panel_log_sel = function(df, labels, title) {
  ggplot(df, aes(x = mj_per_50g_protein, y = hr_per_50g_protein)) +
    geom_line(aes(group = exio_region), color = "grey70", linewidth = 0.5) +
    geom_point(aes(color = type, size = g_protein_per_cap_day), alpha = 0.8) +
    ggrepel::geom_text_repel(data = labels, aes(label = exio_region, fontface = face),
                             size = 5.5, max.overlaps = Inf, box.padding = 0.5,
                             point.padding = 0.3, min.segment.length = 0, segment.color = "grey40",
                             segment.size = 0.3, seed = 1, show.legend = FALSE) +
    scale_color_manual(values = c(hr_f = "#ca2323", hr_m = "#1f77b4"),
                       labels = c(hr_f = "Female", hr_m = "Male")) +
    scale_size_continuous(range = c(1, 8), limits = fig3_size_limits) +
    scale_x_log10(expand = expansion(mult = c(0.05, 0.12)),
                  breaks = scales::breaks_log(n = 6),
                  labels = scales::label_number(drop0trailing = TRUE)) +
    scale_y_log10(expand = expansion(mult = c(0.05, 0.12)),
                  breaks = scales::breaks_log(n = 6),
                  labels = scales::label_number(drop0trailing = TRUE)) +
    annotation_logticks(sides = "bl", color = "grey50", linewidth = 0.3,
                        short = unit(0.08, "cm"), mid = unit(0.14, "cm"),
                        long = unit(0.2, "cm")) +
    labs(x = "Energy used in production (MJ / 50 g protein, log scale)",
         y = "Paid labor time (hr / 50 g protein, log scale)",
         color = "Gender", size = "g protein/cap/day", title = title) +
    theme_fig(legend_position = "right")
}

fig3_labels_food  = fig3_select_labels(fig3_food)
fig3_labels_total = fig3_select_labels(fig3_total)

fig3_loglog_sel = (tradeoff_panel_log_sel(fig3_food, fig3_labels_food, "Food sector") |
                   tradeoff_panel_log_sel(fig3_total, fig3_labels_total, "Food + non-food sectors")) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Paid labor time vs. energy per 50 g of domestically-consumed protein (", year,
                   ") — domestic effort, log-log"),
    caption = "Time is paid labor in production; energy is energy used in production; household (unpaid) time is not included.\nLabelled: frontier regions (bold), the four regions farthest from the centre of the cloud (plain), and the eight Figure 2 countries (italic). Each line joins a region's female and male points.",
    theme = theme(plot.title = element_text(size = FIG_BASE * 1.5, face = "bold", hjust = 0.5,
                                            margin = margin(b = 16)),
                  plot.caption = element_text(size = FIG_BASE * 0.8, hjust = 0.5,
                                              margin = margin(t = 12)))) &
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "fig3_tradeoff_domestic_effort_loglog_labeled.pdf"),
       fig3_loglog_sel, width = 22, height = 10)
ggsave(file.path(fig_dir, "fig3_tradeoff_domestic_effort_loglog_labeled.png"),
       fig3_loglog_sel, width = 22, height = 10, dpi = 200, bg = "white")

cat("Saved selectively labelled log-log Figure 3 to", fig_dir, "\n")
