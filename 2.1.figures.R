#### 2.1.figures.R -- final manuscript figures ####
#
# Re-renders the three manuscript figures at publication font sizes, plus a
# simplified variant of Figure 1. Source AFTER 2.analyze_result.R: every input
# object is built there and listed in `required_objects` below. This script
# computes no new footprints -- presentation only.
#
#   Fig. 1   results/manuscript/fig1_time_footprint.pdf
#            Per-capita food-provisioning time by country and gender, food +
#            non-food, imports split by continent of origin.
#   Fig. 1b  results/manuscript/fig1b_time_footprint_simple.pdf
#            Same bars and same totals, collapsed to unpaid / paid, with the
#            food vs. non-food split kept for domestic, export and import.
#   Fig. 2   results/manuscript/fig2_protein_mosaic.pdf
#            Time and energy conversion factors for domestic vs. imported
#            protein, eight selected countries.
#   Fig. 3   results/manuscript/fig3_tradeoff_domestic_effort.pdf
#            Energy-time tradeoff per 50 g domestically-consumed protein,
#            domestic effort only: food sector (left), food + non-food (right).
#   Fig. S1  results/manuscript/figS1_protein_mosaic_time_all.pdf
#   Fig. S2  results/manuscript/figS2_protein_mosaic_energy_all.pdf
#            Figure 2's two mosaics for every directly-modelled country.
#
# Two .docx hand-offs carry the same figures as vector graphics:
#   manuscript_figures.docx  the four main figures, one per A3 landscape page
#   si_figures.docx          S1 and S2 split across 22 x 22 in pages, since a
#                            single sheet would exceed Word's 22-inch limit

library(tidyverse)
library(patchwork)

required_objects = c(
  "plot_countries", "year",
  # Fig. 1 / 1b
  "b_direct", "direct_ord", "import_levels", "pro_direct", "und_direct",
  "y_max_direct", "y_min_direct", "pro_scale_direct", "und_scale_direct",
  # Fig. 2 and its SI full-sample counterparts
  "df_mosaic", "df_mosaic_energy", "mosaic_fill_colors", "vline_df",
  "protein_max", "protein_max_energy",
  "df_mosaic_nonrow", "df_mosaic_nonrow_energy", "mosaic_nonrow_ord",
  # Fig. 3
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

fig_dir = "results/manuscript"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)


#### Shared publication theme ####

# One knob for the whole figure set. Every other size is a multiple of it, so
# bumping FIG_BASE rescales all four figures consistently.
FIG_BASE = 20

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


#### Figure 1: time footprint by country and gender ####

p_fig1_f = plot_countries(b_direct %>% filter(type == "hr_f"),
                          "Female time (hr/cap/day)", "") +
  geom_point(data = pro_direct, aes(x = country, y = pro_per_cap_day * pro_scale_direct),
             color = "black", size = 3, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(y_min_direct, y_max_direct),
                     sec.axis = sec_axis(~ . / pro_scale_direct, name = "g protein/cap/day"))

p_fig1_m_base = plot_countries(b_direct %>% filter(type == "hr_m"),
                               "Male time (hr/cap/day)", "") +
  geom_point(data = und_direct, aes(x = country, y = undernourishment_pct * und_scale_direct),
             color = "black", size = 3, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(y_min_direct, y_max_direct),
                     sec.axis = sec_axis(~ . / und_scale_direct,
                                         name = "Undernourishment (%)"))

# LUX (male) carries by far the largest import-effort bar and its continent
# makeup isn't readable from color alone, so each import_cont_* segment is
# labelled beside its own block in that block's fill color. With
# position_stack() on negative values the FIRST factor level lands farthest
# from zero, so cumulative-summing in import_levels order reproduces the drawn
# segment bounds exactly.
fig1_fill_scheme = attr(p_fig1_m_base, "fill_scheme")
lux_idx = which(direct_ord == "LUX")

fig1_lux_segments = b_direct %>%
  filter(country == "LUX", type == "hr_m", grepl("^import_cont_", as.character(footprint_type))) %>%
  mutate(footprint_type = factor(as.character(footprint_type), levels = import_levels)) %>%
  arrange(footprint_type) %>%
  mutate(ymax = cumsum(per_capita_value) - sum(per_capita_value),
         ymin = ymax - per_capita_value,
         y_mid = (ymin + ymax) / 2,
         label = paste0(gsub("_nf$", "", gsub("^import_cont_", "", as.character(footprint_type))),
                        ": ", round(per_capita_value, 2), "h"),
         color = fig1_fill_scheme[as.character(footprint_type)]) %>%
  filter(per_capita_value >= 0.15)  # skip slivers too thin for non-overlapping labels

p_fig1_m = p_fig1_m_base +
  { if (length(lux_idx) == 1 && nrow(fig1_lux_segments) > 0)
      geom_text(data = fig1_lux_segments, aes(x = lux_idx + 0.6, y = y_mid, label = label),
                color = fig1_lux_segments$color, hjust = 0, size = 4.5, fontface = "bold",
                inherit.aes = FALSE) }

# ~20 fill categories, so the legend gets its own (smaller) text size and a
# capped row count; everything else follows FIG_BASE.
fig1_legend = theme(legend.text  = element_text(size = FIG_BASE * 0.62),
                    legend.title = element_text(size = FIG_BASE * 0.70),
                    legend.key.size = unit(1.0, "lines"))

fig1 = (p_fig1_f + theme_fig() + fig1_legend) /
       (p_fig1_m + theme_fig(angle_x = 90) + fig1_legend) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

fig1[[1]] <- fig1[[1]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                               axis.title.x = element_blank())
fig1[[2]] <- fig1[[2]] + theme(axis.ticks.x = element_line())

ggsave(file.path(fig_dir, "fig1_time_footprint.pdf"), fig1, width = 20, height = 14)


#### Figure 1b: same bars, collapsed to unpaid / paid x food x trade direction ####

# Collapses Fig. 1's ~20 stacked categories into seven, keeping the food /
# non-food split but dropping the continent-of-origin breakdown of imports and
# the activity breakdown of household time. The seven groups partition
# b_direct exactly (asserted below), so each country's bar totals match Fig. 1.
# "preparation_econ" is paid restaurant/food-service preparation, so it belongs
# with paid food-sector domestic time. "_nf" marks non-food-sector components.
fig1b_groups = c("Unpaid (household)",
                 "Paid: food, domestic", "Paid: non-food, domestic",
                 "Paid: food, export",   "Paid: non-food, export",
                 "Paid: food, import",   "Paid: non-food, import")

# Same hues plot_countries() uses in Fig. 1: warm for unpaid, blues for
# food-sector paid, greens for non-food-sector paid, darkest = domestic.
fig1b_colors = c("Unpaid (household)"       = "#feb24c",
                 "Paid: food, domestic"     = "#08519c",
                 "Paid: non-food, domestic" = "#31a354",
                 "Paid: food, export"       = "#3182bd",
                 "Paid: non-food, export"   = "#74c476",
                 "Paid: food, import"       = "#6baed6",
                 "Paid: non-food, import"   = "#a1d99b")

b_fig1b = b_direct %>%
  filter(type %in% c("hr_m", "hr_f")) %>%
  mutate(ft = as.character(footprint_type),
         sector = ifelse(grepl("_nf$", ft), "non-food", "food"),
         group = case_when(
           grepl("_non\\.econ$", ft)                 ~ "Unpaid (household)",
           ft == "preparation_econ"                  ~ "Paid: food, domestic",
           grepl("^domestic_per_capita", ft)         ~ paste0("Paid: ", sector, ", domestic"),
           grepl("^export_per_capita", ft)           ~ paste0("Paid: ", sector, ", export"),
           grepl("^import_cont_", ft)                ~ paste0("Paid: ", sector, ", import"),
           .default = NA_character_))

# Every b_direct row must land in exactly one group, or the bar totals would
# silently stop matching Fig. 1 -- the whole point of this variant.
stopifnot(!any(is.na(b_fig1b$group)), all(b_fig1b$group %in% fig1b_groups))

b_fig1b = b_fig1b %>%
  group_by(country, type, group) %>%
  summarise(per_capita_value = sum(per_capita_value, na.rm = TRUE), .groups = "drop") %>%
  mutate(group = factor(group, levels = fig1b_groups),
         sign = ifelse(grepl("import$", as.character(group)), -1, 1))

# The seven groups partition b_direct and split positive/negative the same way
# Fig. 1 does, so the stack bounds are identical -- reuse Fig. 1's limits and
# dot scaling rather than recomputing them, and the two figures stay directly
# comparable bar for bar.
plot_fig1b_panel = function(gender, ylabel, dots, dot_col, dot_scale, sec_name) {
  df = b_fig1b %>% filter(type == gender)
  ggplot(df %>% filter(sign == 1),
         aes(x = country, y = per_capita_value, fill = group)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_bar(data = df %>% filter(sign == -1),
             aes(x = country, y = -per_capita_value, fill = group),
             stat = "identity", position = "stack") +
    geom_point(data = dots, aes(x = country, y = .data[[dot_col]] * dot_scale),
               color = "black", size = 3, inherit.aes = FALSE) +
    scale_fill_manual(values = fig1b_colors, name = "Footprint type", drop = FALSE) +
    scale_y_continuous(limits = c(y_min_direct, y_max_direct),
                       sec.axis = sec_axis(~ . / dot_scale, name = sec_name)) +
    labs(x = "Country (ISO3)", y = ylabel)
}

p_fig1b_f = plot_fig1b_panel("hr_f", "Female time (hr/cap/day)",
                             pro_direct, "pro_per_cap_day", pro_scale_direct,
                             "g protein/cap/day")
p_fig1b_m = plot_fig1b_panel("hr_m", "Male time (hr/cap/day)",
                             und_direct, "undernourishment_pct", und_scale_direct,
                             "Undernourishment (%)")

fig1b = (p_fig1b_f + theme_fig()) / (p_fig1b_m + theme_fig(angle_x = 90)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

fig1b[[1]] <- fig1b[[1]] + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                                 axis.title.x = element_blank())
fig1b[[2]] <- fig1b[[2]] + theme(axis.ticks.x = element_line())

ggsave(file.path(fig_dir, "fig1b_time_footprint_simple.pdf"), fig1b, width = 20, height = 13)


#### Figure 2: domestic vs. imported protein mosaic (time and energy) ####

# Labels sit in a gutter beside their own column, not inside or above their
# rectangle. df_mosaic's fits_inside / label_y rule parks every label too small
# to fit at a near-identical height just above its rectangle, so a column of
# thin stacked bands ends up with all four labels on top of each other; simply
# repelling them apart does not help either, because the labels then form one
# vertical stack over a set of anchors that are themselves nearly coincident,
# and no reader can tell which label belongs to which band.
#
# Instead: each column's labels are fanned out to a guaranteed minimum vertical
# spacing in their original stacking order, parked in the margin on that
# column's own side (domestic left, import right), and joined to their band by
# a short leader drawn in the band's own colour. Order, side, and colour give
# three independent cues, and the layout is deterministic rather than the
# outcome of a physics solver, so the figure is identical run to run.
mosaic_text_colors = c("food | Domestic effort"     = "#0d47a1",
                       "food | Import effort"       = "#5b9bd5",
                       "non-food | Domestic effort" = "#b3521a",
                       "non-food | Import effort"   = "#d2854a")

# Push values apart to at least `gap`, preserving order, then recentre on the
# original mean and slide back inside `limits` if the fan overshot an edge.
fan_labels = function(y, gap, limits) {
  ord = order(y)
  ys = y[ord]
  for (i in seq_along(ys)[-1]) ys[i] = max(ys[i], ys[i - 1] + gap)
  ys = ys - mean(ys) + mean(y)
  if (min(ys) < limits[1]) ys = ys + (limits[1] - min(ys))
  if (max(ys) > limits[2]) ys = ys - (max(ys) - limits[2])
  out = numeric(length(y))
  out[ord] = ys
  out
}

mosaic_label_layout = function(df, label_col, x_span, y_span, gap_frac, pad_frac) {
  df %>%
    filter(!is.na(.data[[label_col]])) %>%
    mutate(lab = .data[[label_col]]) %>%
    group_by(country, source) %>%
    mutate(y_rect = (ymin + ymax) / 2,
           y_lab  = fan_labels(y_rect, gap_frac * y_span, c(0, y_span)),
           x_edge = ifelse(source == "Domestic", xmin, xmax),
           x_lab  = ifelse(source == "Domestic", xmin - pad_frac * x_span,
                                                 xmax + pad_frac * x_span),
           h_just = ifelse(source == "Domestic", 1, 0)) %>%
    ungroup()
}

# vline_data = NULL omits the reference line entirely (the SI panels below span
# every country, where a single-country reference adds clutter rather than
# context). gap_frac tracks label_size: smaller type needs less vertical room.
mosaic_panel = function(df, label_col, ylabel, title, vline_data = NULL, xintercept = NA,
                        ncol = 3, label_size = 5, gap_frac = 0.085, pad_frac = 0.05,
                        gutter = c(0.32, 0.34), base = FIG_BASE) {
  lab_df = mosaic_label_layout(df, label_col,
                               x_span = max(df$xmax, na.rm = TRUE),
                               y_span = max(df$ymax, na.rm = TRUE),
                               gap_frac = gap_frac, pad_frac = pad_frac)
  g = ggplot(df) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_grp),
              color = "white", linewidth = 0.4)
  # Reference line before the labels, so a label that happens to sit on it
  # stays readable rather than being struck through.
  if (!is.null(vline_data)) {
    g = g +
      geom_vline(data = vline_data, aes(xintercept = xintercept),
                 linetype = "dashed", color = "black", linewidth = 0.5) +
      geom_text(data = vline_data, aes(x = xintercept, y = Inf, label = "USA"),
                color = "black", size = label_size, fontface = "bold",
                hjust = 1.1, vjust = 1.5)
  }
  g +
    geom_segment(data = lab_df,
                 aes(x = x_lab, y = y_lab, xend = x_edge, yend = y_rect, color = fill_grp),
                 linewidth = 0.4) +
    geom_text(data = lab_df,
              aes(x = x_lab, y = y_lab, label = lab, color = fill_grp, hjust = h_just),
              size = label_size, fontface = "bold") +
    scale_color_manual(values = mosaic_text_colors, guide = "none") +
    # Room for the label gutters on both sides of every facet. The energy
    # mosaic's labels ("5.7 GJ/yr") are the long ones that set this.
    scale_x_continuous(expand = expansion(mult = gutter)) +
    facet_wrap(~country, ncol = ncol) +
    scale_fill_manual(values = mosaic_fill_colors, name = "") +
    labs(x = "Protein supply (g/cap/day)", y = ylabel, title = title) +
    theme_fig(base = base, legend_position = "bottom")
}

p_fig2_time = mosaic_panel(
  df_mosaic, "label_min", "Time conversion factor (min / g protein)",
  "Daily time embodied in domestic/imported protein provisioning",
  vline_data = vline_df, xintercept = protein_max)

p_fig2_energy = mosaic_panel(
  df_mosaic_energy, "label_gj", "Energy conversion factor (MJ / g protein)",
  "Yearly energy embodied in domestic/imported protein provisioning",
  vline_data = vline_df, xintercept = protein_max_energy)

# Both panels share the fill scale, so guides="collect" leaves one centered
# bottom legend. plot_spacer() opens the gap: each mosaic spans 3 facet columns
# within its width unit, so 1/12 = (1/4)/3 is a quarter of one facet's width.
fig2 = (p_fig2_time | plot_spacer() | p_fig2_energy) +
  plot_layout(guides = "collect", widths = c(1, 1/12, 1)) &
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "fig2_protein_mosaic.pdf"), fig2, width = 28, height = 12)


#### SI: the same two mosaics for every directly-modelled country ####

# Figure 2 shows eight selected countries; these are the full-sample versions
# for the SI, one facet per directly-modelled (non-RoW) country. Same builder,
# so the reading conventions carry over exactly -- only the facet grid, type
# size and gutter widths shrink to fit the larger grid. No USA reference line:
# across the whole sample a single-country marker is clutter, not context.
#
# df_mosaic_nonrow carries its stack in hr / g protein, where df_mosaic uses
# min / g protein. Converted here so the SI axis matches the main figure's and
# a reader can move between them without a unit change.
#
# Facet order: 2.analyze_result.R computes mosaic_nonrow_ord (countries by
# domestic time per capita, descending) but its joins leave df_mosaic_nonrow's
# country column as character, so facet_wrap falls back to alphabetical. The
# intended order is restored here.
si_order = function(df) {
  df %>% mutate(country = factor(as.character(country),
                                 levels = intersect(mosaic_nonrow_ord,
                                                    unique(as.character(country)))))
}

df_si_mosaic_time   = df_mosaic_nonrow %>%
  mutate(across(c(ymin, ymax), ~ .x * 60)) %>%
  si_order()
df_si_mosaic_energy = df_mosaic_nonrow_energy %>% si_order()

si_countries = nlevels(df_si_mosaic_time$country)
si_ncol      = 6
si_height    = ceiling(si_countries / si_ncol) * 3.6 + 2

si_mosaic = function(df, label_col, ylabel, title, ncol = si_ncol, label_size = 3.5) {
  mosaic_panel(df, label_col, ylabel, title,
               ncol = ncol, label_size = label_size, gap_frac = 0.07, pad_frac = 0.04,
               gutter = c(0.30, 0.32), base = 16)
}

si_title_time   = "Daily time embodied in domestic/imported protein provisioning"
si_title_energy = "Yearly energy embodied in domestic/imported protein provisioning"
si_ylab_time    = "Time conversion factor (min / g protein)"
si_ylab_energy  = "Energy conversion factor (MJ / g protein)"

p_si_mosaic_time = si_mosaic(
  df_si_mosaic_time, "label_min", si_ylab_time,
  paste0(si_title_time, ", all modelled countries"))

p_si_mosaic_energy = si_mosaic(
  df_si_mosaic_energy, "label_gj", si_ylab_energy,
  paste0(si_title_energy, ", all modelled countries"))

ggsave(file.path(fig_dir, "figS1_protein_mosaic_time_all.pdf"), p_si_mosaic_time,
       width = 26, height = si_height, limitsize = FALSE)
ggsave(file.path(fig_dir, "figS2_protein_mosaic_energy_all.pdf"), p_si_mosaic_energy,
       width = 26, height = si_height, limitsize = FALSE)


#### SI, Word edition: the same mosaics split across page-sized chunks ####

# Word's maximum page dimension is 22 inches, so the single-sheet SI figures
# above (26 x ~30 in) can only go into a .docx shrunk by about a third, which
# costs the label legibility the gutter layout exists to protect. Instead each
# SI mosaic is split into chunks of countries that fit one page at full scale:
# 4 columns x 4 rows per page, 18 x 16.4 in native, inside a 22 x 22 in page.
# Country order is preserved across chunks, so the pages read as one sequence.
SI_DOCX_NCOL = 4
SI_DOCX_ROWS = 4
si_per_page  = SI_DOCX_NCOL * SI_DOCX_ROWS
si_docx_page = 22
si_docx_mar  = 0.4
si_chunk_w   = 18

si_chunks = split(levels(df_si_mosaic_time$country),
                  ceiling(seq_len(si_countries) / si_per_page))

# One page per chunk, for each of the two mosaics. Height follows the rows the
# chunk actually fills, so a short final chunk gets a short figure rather than
# a page of whitespace.
si_docx_pages = unlist(recursive = FALSE, lapply(
  list(list(df = df_si_mosaic_time,   col = "label_min", ylab = si_ylab_time,
            title = si_title_time,    tag = "S1"),
       list(df = df_si_mosaic_energy, col = "label_gj",  ylab = si_ylab_energy,
            title = si_title_energy,  tag = "S2")),
  function(spec) lapply(seq_along(si_chunks), function(i) {
    isos = si_chunks[[i]]
    list(plot = si_mosaic(spec$df %>% filter(country %in% isos) %>%
                            mutate(country = droplevels(country)),
                          spec$col, spec$ylab, spec$title,
                          ncol = SI_DOCX_NCOL, label_size = 3.8),
         w = si_chunk_w,
         h = ceiling(length(isos) / SI_DOCX_NCOL) * 3.6 + 2,
         caption = sprintf("Figure %s (%d of %d). %s. Countries %d-%d of %d, ordered by domestic time per capita.",
                           spec$tag, i, length(si_chunks), spec$title,
                           (i - 1) * si_per_page + 1,
                           (i - 1) * si_per_page + length(isos), si_countries))
  })))


#### Figure 3: energy-time tradeoff, domestic effort ####

# Each region is a vertical segment at its (gender-invariant) energy value,
# running from male to female hours. Domestic effort only -- the effort that
# occurred in the same country that consumed the protein; the import-effort
# facet of the source plots is dropped here. Region-resolved (49 EXIO
# regions), since effort_consumption_df is region-keyed on both sectors.
# Taiwan is excluded: its non-food intensity is an outlier that compresses the
# rest of the scatter (diagnosed in 9r.diagnose_twn_nonfood_outlier.R).
fig3_drop_regions = c("Taiwan")

fig3_subset = function(df) {
  df %>% filter(effort_label == "Domestic effort", !exio_region %in% fig3_drop_regions)
}

fig3_food        = fig3_subset(tradeoff_protein_foodecon_domestic_effort)
fig3_total       = fig3_subset(tradeoff_protein_totalecon_domestic_effort)
fig3_food_label  = fig3_subset(label_tradeoff_protein_foodecon_domestic_effort)
fig3_total_label = fig3_subset(label_tradeoff_protein_totalecon_domestic_effort)

# Shared point-size limits across both panels: same protein supply draws the
# same dot on either side, and patchwork's guides="collect" can then merge the
# two size legends into one instead of printing it twice.
fig3_size_limits = range(c(fig3_food$g_protein_per_cap_day,
                           fig3_total$g_protein_per_cap_day), na.rm = TRUE)

tradeoff_panel = function(df, labels, title) {
  ggplot(df, aes(x = mj_per_50g_protein, y = hr_per_50g_protein)) +
    geom_line(aes(group = exio_region), color = "grey60", linewidth = 0.6) +
    geom_point(aes(color = type, size = g_protein_per_cap_day), alpha = 0.85) +
    ggrepel::geom_text_repel(data = labels, aes(label = exio_region),
                             size = 5, max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = c(hr_f = "#ca2323", hr_m = "#1f77b4"),
                       labels = c(hr_f = "Female", hr_m = "Male")) +
    scale_size_continuous(range = c(1, 8), limits = fig3_size_limits) +
    # Extra upper expansion so the largest points (8 mm radius) aren't clipped
    # by the panel border; the default 5% is sized for point position, not radius.
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(x = "Energy (MJ / 50 g protein)", y = "Time (hr / 50 g protein)",
         color = "Gender", size = "g protein/cap/day", title = title) +
    theme_fig(legend_position = "right")
}

p_fig3_food  = tradeoff_panel(fig3_food,  fig3_food_label,  "Food sector")
p_fig3_total = tradeoff_panel(fig3_total, fig3_total_label, "Food + non-food sectors")

fig3 = (p_fig3_food | p_fig3_total) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Energy vs. time per 50 g of domestically-consumed protein (", year,
                   ") — domestic effort"),
    theme = theme(plot.title = element_text(size = FIG_BASE * 1.5, face = "bold", hjust = 0.5,
                                            margin = margin(b = 16)))) &
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "fig3_tradeoff_domestic_effort.pdf"), fig3, width = 22, height = 10)


#### Figures into .docx files as vector graphics ####

# Word gets EMF+ (Enhanced Metafile), not PNG: the figures stay vector, so they
# scale without pixelation and Word can Ungroup them into editable shapes and
# text. rvg is deliberately not used -- its Word entry point (body_add_vg) is
# defunct as of rvg 0.4.2 and now redirects to a raster path, and dml() objects
# have no body_add() method for .docx (they are the PowerPoint route).
# officer::body_add_gg() would also rasterize. The PDFs above remain the
# archival copies; this .docx is the hand-off format.
# emfPlus = TRUE is required for the alpha blending these figures use.
library(officer)
library(devEMF)

emf_render = function(plot, width, height) {
  path = tempfile(fileext = ".emf")
  emf(path, width = width, height = height, emfPlus = TRUE, emfPlusFont = TRUE)
  print(plot)
  dev.off()
  path
}

# A3 landscape. These canvases are 20-28 inches wide; on A4 they would shrink
# far past legibility, and the whole point of FIG_BASE is readable type.
docx_page_w = 16.54
docx_page_h = 11.69
docx_margin = 0.4
docx_max_w  = docx_page_w - 2 * docx_margin
docx_max_h  = docx_page_h - 2 * docx_margin

# Each figure keeps the aspect ratio of its own ggsave() call above and is
# scaled to whichever of page width / page height binds first.
docx_figures = list(
  list(plot = fig1,  w = 20, h = 14,
       caption = "Figure 1. Per-capita food-provisioning time by country and gender, food and non-food sectors, imports split by continent of origin."),
  list(plot = fig1b, w = 20, h = 13,
       caption = "Figure 1b. As Figure 1, collapsed to unpaid and paid time with the food / non-food split retained. Country totals match Figure 1."),
  list(plot = fig2,  w = 28, h = 12,
       caption = "Figure 2. Time and energy conversion factors for domestically produced versus imported protein, eight selected countries."),
  list(plot = fig3,  w = 22, h = 10,
       caption = "Figure 3. Energy-time tradeoff per 50 g of domestically-consumed protein, domestic effort only: food sector (left), food and non-food (right).")
)

# One figure per page, captioned, scaled to whichever of page width / height
# binds first -- never enlarged past its native size, so type stays at the size
# it was designed at.
write_figure_docx = function(figures, target, page_w, page_h, margin, orient) {
  max_w = page_w - 2 * margin
  max_h = page_h - 2 * margin
  doc = read_docx()
  for (i in seq_along(figures)) {
    f = figures[[i]]
    fit = min(1, max_w / f$w, max_h / f$h)
    if (i > 1) doc = body_add_break(doc)
    doc = doc %>%
      body_add_par(f$caption, style = "Normal") %>%
      body_add_img(src = emf_render(f$plot, f$w, f$h),
                   width = f$w * fit, height = f$h * fit)
  }
  doc = body_set_default_section(doc, prop_section(
    page_size = page_size(width = page_w, height = page_h, orient = orient),
    page_margins = page_mar(top = margin, bottom = margin, left = margin, right = margin,
                            header = 0, footer = 0, gutter = 0)))
  print(doc, target = target)
  invisible(target)
}

write_figure_docx(docx_figures, file.path(fig_dir, "manuscript_figures.docx"),
                  docx_page_w, docx_page_h, docx_margin, "landscape")

# SI mosaics, one chunk of countries per page. 22 x 22 in is Word's largest
# square page, which lets each 18 x 16.4 in chunk embed at full scale.
write_figure_docx(si_docx_pages, file.path(fig_dir, "si_figures.docx"),
                  si_docx_page, si_docx_page, si_docx_mar, "portrait")
