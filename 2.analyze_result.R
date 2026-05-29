#### Result analysis ####


#### 1. Energy and Labor footprint ####

# Load economic (direct food sector) and non-economic (indirect non-food sector) satellite data
l_int_d <- readRDS(file = paste0("data/FABIO_exio_satellites_food_", year, ".rds"))
l_int_i <- readRDS(file = paste0("data/FABIO_exio_satellites_nonfood_", year, ".rds"))

# Footprint summed at the FABIO country level
fp_food <- lapply(l_int_d, function(d) Matrix::Diagonal(x=d) %*% FABIO_x_hh)

n_nf = length(exio_nonfood_sectors) # EXIO non-food sectors
# Footprint for non-food sectors is calculated by multiplying the intensity matrix 
# (in the same dim as satellite (32725x23001)) with the FABIO_x_hh vector (23001x187) to get a matrix of dimension (32725x187), and then summing every 187 rows to get a matrix of dimension (187x187) where rows are origin countries and columns are target countries.

fp_nonfood <- lapply(l_int_i, function(d) {
  fp_list <- vector("list", nrreg)
  for (i in 1:nrreg) {
    print(paste("Processing region", i, "of", nrreg))
    A_i <- d[((i-1)*n_nf + 1):(i*n_nf), ]
    # Build B_i as block diagonal: each column k gets FABIO_x_hh commodities for region k, source i
    B_i <- Matrix::Diagonal(x = as.vector(FABIO_x_hh[, i])) %*% 
      kronecker(Matrix::Diagonal(nrreg), Matrix::Matrix(1, nrcom, 1))
    fp_list[[i]] <- A_i %*% B_i
  }
  do.call(rbind, fp_list) %>% as("CsparseMatrix")
})

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

l_food_country = lapply(fp_food, agg_country_footprint)
l_nonfood_country = lapply(fp_nonfood, agg_country_footprint)

# Validate: sum(l_food_country[[i]]) == sum(fp_food[[i]])

# Validate: find the ten biggest cell from l_food_country[[2]] and check source and destination.
find_top_cells(l_food_country[[2]], matrix_name = names(l_food_country)[2])
find_top_cells(l_food_country[[3]], matrix_name = names(l_food_country)[3])
find_top_cells(l_nonfood_country[[2]], matrix_name = names(l_nonfood_country)[2])
find_top_cells(l_nonfood_country[[3]], matrix_name = names(l_nonfood_country)[3])



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

summary_food = lapply(l_food_country, country_summary)
summary_nonfood = lapply(l_nonfood_country, country_summary)

# Make per capita labor hour footprints daily by dividing by 365
for (i in names(summary_food[2:3])) {
  summary_food[[i]] %>% mutate(across(ends_with("per_capita"), ~ .x / 365)) -> summary_food[[i]]
  summary_nonfood[[i]] %>% mutate(across(ends_with("per_capita"), ~ .x / 365)) -> summary_nonfood[[i]]
}

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
summary_nonfood_df <- summary_nonfood_df %>% mutate(is_row = country %in% row_countries)

# Order of domestic hours (by female)
sum_ord = (summary_food_df %>%   
    filter(type %in% c("hr_f")) %>% 
    group_by(country) %>% 
    summarise(d = sum(domestic_per_capita)) %>% 
    arrange(-d))$country

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
  # filter(type %in% c("hr_m", "hr_f")) %>%
  pivot_longer(cols = c("domestic_per_capita", "export_per_capita", "import_per_capita"), 
               names_to = "footprint_type", values_to = "per_capita_value") %>%
  mutate(country = factor(country, levels = sum_ord),
         footprint_type = factor(footprint_type, levels = c("export_per_capita", "import_per_capita", "domestic_per_capita")))


# Vertically stack df_ghd_gender to summary_food_df_long, and then plot again with the same function plot_countries. This will add the non-economic food time to the existing labor footprint
summary_food_df_long_with_ghd = bind_rows(summary_food_df_long, df_ghd_combined) %>%
  arrange(country, type, footprint_type) %>%
  mutate(is_row = as.character(country) %in% row_countries,
         country = factor(country, levels = sum_ord),
         footprint_type = factor(footprint_type,
                                 levels = c("preparation_non.econ", "processing_non.econ", "growth_collection_non.econ",
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

p_hr_f_nonfood = plot_countries(summary_nonfood_df_long %>% filter(type %in% c("hr_f")),
               "Daily time footprint per capita (hr/cap/day)",
               paste0("Female nonfood time footprint per capita by country (", year, ")"))

p_hr_m_nonfood = plot_countries(summary_nonfood_df_long %>% filter(type %in% c("hr_m")),
               "Daily time footprint per capita (hr/cap/day)",
               paste0("Male nonfood time footprint per capita by country (", year, ")"))

p_en_nonfood = plot_countries(summary_nonfood_df_long %>% filter(type %in% c("en")),
               "Nonfood sector energy footprint per capita (GJ/cap/yr)",
               paste0("Energy footprint per capita by country (", year, ")"))

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
a = summary_food_df_long_with_ghd %>% filter(country %in% partial_cty) %>%
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

p1 = plot_countries(a %>% filter(type %in% c("hr_f")), "Female time footprint per capita (hr/day)", "") +
  geom_point(data = pro_partial, aes(x = country, y = pro_per_cap_day * pro_scale),
             color = "black", size = 2.5, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(-1, 3.5),
                     sec.axis = sec_axis(~ . / pro_scale, name = "g protein/cap/day"))
p2 = plot_countries(a %>% filter(type %in% c("hr_m")), "Male time footprint per capita (hr/day)", "") 
p3 = plot_countries(a %>% filter(type %in% c("en")), "Energy footprint per capita (GJ/yr)", "")

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
ggsave(paste0("results/footprint_all_countries.pdf"), p_combined, width = 18, height = 12)
ggsave(paste0("results/footprint_partial_countries.pdf"), p_combined_partial, width = 18, height = 12)

# Non-food sector: combined female + male + energy, ordered by non-food female time
nonfood_ord = (summary_nonfood_df_long %>%
  filter(type == "hr_m", footprint_type == "domestic_per_capita") %>%
  arrange(-per_capita_value))$country

p_hr_f_nf = plot_countries(
  summary_nonfood_df_long %>% filter(type == "hr_f") %>%
    mutate(country = factor(country, levels = nonfood_ord)),
  "Daily time footprint per capita (hr/cap/day)",
  paste0("Female non-food time footprint per capita by country (", year, ")"))

p_hr_m_nf = plot_countries(
  summary_nonfood_df_long %>% filter(type == "hr_m") %>%
    mutate(country = factor(country, levels = nonfood_ord)),
  "Daily time footprint per capita (hr/cap/day)",
  paste0("Male non-food time footprint per capita by country (", year, ")"))

p_en_nf = plot_countries(
  summary_nonfood_df_long %>% filter(type == "en") %>%
    mutate(country = factor(country, levels = nonfood_ord)),
  "Non-food sector energy footprint per capita (GJ/cap/yr)",
  paste0("Non-food energy footprint per capita by country (", year, ")"))

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
  filter(!is_row) %>%
  mutate(country = factor(country, levels = direct_ord))

p_hr_f_direct = plot_countries(b %>% filter(type == "hr_f"), "Female time footprint per capita (hr/day)", "") + ylim(-1, 3.5)
p_hr_m_direct = plot_countries(b %>% filter(type == "hr_m"), "Male time footprint per capita (hr/day)", "") + ylim(-1, 3.5)

p_combined_direct = p_hr_f_direct / p_hr_m_direct + plot_layout(guides = "collect") & theme(
  legend.position = "top",
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()
)
p_combined_direct[[2]] <- p_combined_direct[[2]] +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.x = element_line())
print(p_combined_direct)

ggsave(paste0("results/footprint_direct_countries.pdf"), p_combined_direct, width = 18, height = 12)





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

p_kcal = plot_countries(summary_kcal_df_long, "Daily kcal supply per capita (kcal/cap/day)", "kcal")
p_protein = plot_countries(summary_pro_df_long, "Daily protein supply per capita (g/cap/day)", "protein")
ggsave("results/kcal_supply.pdf",    p_kcal,    width = 18, height = 6)
ggsave("results/protein_supply.pdf", p_protein, width = 18, height = 6)


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
  mutate(hr_per_100g_protein = per_capita_value_time / per_capita_value_protein * 100)

v_ord_protein_econlabor = (summary_time_protein %>%
  filter(type == "hr_f", cat == "domestic", grepl("per_capita", footprint_type_time)) %>%
  arrange(-hr_per_100g_protein))$country

df_convfac_protein_econlabor    = summary_time_protein %>% filter(grepl("per_capita", footprint_type_time))
df_convfac_protein_nonecon = summary_time_protein %>% filter(country %in% cty_ghd)

v_ord_protein_alltime = (df_convfac_protein_nonecon %>% 
               filter(type %in% c("hr_f"), cat=="domestic") %>% 
               group_by(country) %>% summarise(d = sum(hr_per_100g_protein, na.rm=TRUE)) %>%
               arrange(-d))$country

p_conversion_protein_econlabor = ggplot(
  df_convfac_protein_econlabor %>%
    filter(footprint_type_time == "domestic_per_capita") %>%
    select(country, type, footprint_type_time, hr_per_100g_protein) %>%
    mutate(country = factor(country, levels = v_ord_protein_econ)),
  aes(x = country, y = hr_per_100g_protein, fill = footprint_type_time)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~type, ncol = 1, scales = "fixed") +
  labs(x = "Country (ISO3)", y = "Time per 100 g protein (hr/100g)", fill = "Footprint type",
       title = paste0("Time per 100 g protein by country (", year, ") - Economic time only")) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("domestic_per_capita" = "#1f77b4",
                               "export_per_capita"   = "#2ca02c",
                               "import_per_capita"   = "#ff7f0e")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_conversion_protein_nonecon = ggplot(
  df_convfac_protein_nonecon %>%
    select(country, type, footprint_type_time, hr_per_100g_protein) %>%
    mutate(country = factor(country, levels = v_ord_alltime)),
  aes(x = country, y = hr_per_100g_protein, fill = footprint_type_time)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~type, ncol = 1, scales = "fixed") +
  labs(x = "Country (ISO3)", y = "Time per 100 g protein (hr/100g)", fill = "Footprint type",
       title = paste0("Time per 100 g protein by country (", year, ")")) +
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


#### Country spotlight ####

make_spotlight_data = function(iso) {
  time_energy = bind_rows(
    summary_food_df_long    %>% mutate(sector = "food sector"),
    summary_nonfood_df_long %>% mutate(sector = "non-food sector")
  ) %>%
    filter(as.character(country) == iso) %>%
    mutate(footprint_type = as.character(footprint_type),
           per_capita_value = ifelse(type == "en", per_capita_value / 1e3, per_capita_value))

  nutrition = bind_rows(
    summary_kcal_df_long %>% mutate(type = "kcal",    sector = "—"),
    summary_pro_df_long  %>% mutate(type = "protein", sector = "—")
  ) %>%
    filter(as.character(country) == iso) %>%
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
        labels = c("Female labor\n(hr/cap/day)", "Male labor\n(hr/cap/day)",
                   "Energy\n(GJ/cap/yr)", "Food supply\n(kcal/cap/day)", "Protein supply\n(g/cap/day)"))
    )
}

make_continent_data = function(iso) {
  pop = subset(countrypops, year == yr & country_code_3 == iso)$population
  bind_rows(lapply(list("food sector" = l_food_country, "non-food sector" = l_nonfood_country),
                   function(l_country) {
    bind_rows(lapply(names(l_country), 
      function(metric) {
        mat  = l_country[[metric]]
        cty  = colnames(mat)
        cont = regions$continent[match(cty, regions$iso3c)]
        not_focal = cty != iso

        exp_by_cont = tapply(mat[iso, not_focal],  cont[not_focal], sum, na.rm = TRUE)
        imp_by_cont = tapply(mat[not_focal, iso],  cont[not_focal], sum, na.rm = TRUE)

        denom = if (metric == "en") pop / 1e3 else pop / 1e6 * 365

        bind_rows(
          data.frame(continent = names(exp_by_cont), value = as.numeric(exp_by_cont) / denom,
                    flow = "export"),
          data.frame(continent = names(imp_by_cont), value = as.numeric(imp_by_cont) / denom,
                    flow = "import")
        ) %>% mutate(metric = metric)
      }
      )
    )
  }), .id = "sector") %>%
    mutate(metric = factor(metric,
      levels = c("hr_f", "hr_m", "en"),
      labels = c("Female labor\n(hr/cap/day)", "Male labor\n(hr/cap/day)", "Energy\n(GJ/cap/yr)")))
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
  ggplot(df, aes(x = flow, y = value, fill = continent)) +
    geom_col(width = 0.6, position = "stack") +
    facet_grid(metric ~ sector, scales = "free_y") +
    scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = NULL, fill = "Continent",
         title = paste0(iso, " food trade flows by continent (", year, ")")) +
    theme_minimal() +
    theme(legend.position = "right", strip.text = element_text(size = 10))
}

for (iso in c("USA", "CHN", "IND")) {
  ggsave(paste0("results/", tolower(iso), "_spotlight.pdf"),  plot_spotlight(iso), width = 14, height = 5)
  ggsave(paste0("results/", tolower(iso), "_continent.pdf"),  plot_continent(iso), width = 10, height = 10)
}

# Combined three-country comparison: absolute footprints + domestic/import shares
make_select_country_data = function(isos = c("USA", "AUT", "CHN", "IND")) {
  bind_rows(lapply(isos, function(iso) {
    make_spotlight_data(iso) %>% mutate(country = iso)
  })) %>%
    filter(flow %in% c("domestic", "import")) %>%
    mutate(country = factor(country, levels = isos),
           flow    = factor(flow, levels = c("domestic", "import")))
}

plot_select_country = function(isos = c("USA", "AUT", "CHN", "IND")) {
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
  mutate(mj_per_100g_protein = (per_capita_value_energy / 365) / per_capita_value_protein * 100)

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
  select(country, type, hr_per_100g_protein) %>%
  left_join(summary_energy_protein %>%
              filter(footprint_type == "domestic_per_capita") %>%
              select(country, mj_per_100g_protein, g_protein_per_cap_day = per_capita_value_protein),
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
  summarise(hr_per_100g_protein = sum(hr_per_100g_protein, na.rm = TRUE), .groups = "drop") %>%
  left_join(summary_energy_protein %>%
              filter(footprint_type == "domestic_per_capita") %>%
              select(country, mj_per_100g_protein, g_protein_per_cap_day = per_capita_value_protein),
            by = "country") %>%
  left_join(regions %>% select(iso3c, continent), by = c("country" = "iso3c")) %>%
  drop_na() %>%
  mutate(is_row = as.character(country) %in% row_countries)

# Scatter plot helper: energy (x) vs time (y), sized by nutrition supply per capita
tradeoff_scatter <- function(df, x_col, y_col, size_col, x_lab, y_lab, size_lab, title) {
  has_row <- "is_row" %in% colnames(df)
  base_aes <- aes(x = .data[[x_col]], y = .data[[y_col]], size = .data[[size_col]],
                  color = continent, label = country)
  g <- ggplot(df, base_aes)
  if (has_row) {
    g <- g +
      geom_point(aes(alpha = is_row)) +
      ggrepel::geom_text_repel(aes(alpha = is_row), size = 2.5, max.overlaps = 20, show.legend = FALSE) +
      scale_alpha_manual(values = c("TRUE" = 0.3, "FALSE" = 0.9), guide = "none")
  } else {
    g <- g +
      geom_point(alpha = 0.9) +
      ggrepel::geom_text_repel(size = 2.5, max.overlaps = 20, show.legend = FALSE)
  }
  g +
    scale_size_continuous(range = c(1, 8)) +
    facet_wrap(~type, nrow = 1,
               labeller = labeller(type = c("hr_m" = "Male", "hr_f" = "Female"))) +
    labs(x = x_lab, y = y_lab, size = size_lab,
         color = "Continent", title = title) +
    theme_minimal() +
    theme(legend.position = "right")
}

p_tradeoff_kcal_econlabor = tradeoff_scatter(
  tradeoff_kcal_econlabor, "mj_per_2000kcal", "hr_per_2000kcal", "kcal_per_cap_day",
  "Energy (MJ / 2000 kcal)", "Time (hr / 2000 kcal)", "kcal/cap/day",
  paste0("Energy vs. time to provision 2000 kcal (", year, ") — Econ. labor"))

p_tradeoff_kcal_allwork = tradeoff_scatter(
  tradeoff_kcal_allwork, "mj_per_2000kcal", "hr_per_2000kcal", "kcal_per_cap_day",
  "Energy (MJ / 2000 kcal)", "Time (hr / 2000 kcal)", "kcal/cap/day",
  paste0("Energy vs. time to provision 2000 kcal (", year, ") — Econ. + non-econ labor"))

p_tradeoff_protein_econlabor = tradeoff_scatter(
  tradeoff_protein_econlabor, "mj_per_100g_protein", "hr_per_100g_protein", "g_protein_per_cap_day",
  "Energy (MJ / 100 g protein)", "Time (hr / 100 g protein)", "g protein/cap/day",
  paste0("Energy vs. time to provision 100 g protein (", year, ") — Econ. labor"))

p_tradeoff_protein_allwork = tradeoff_scatter(
  tradeoff_protein_allwork, "mj_per_100g_protein", "hr_per_100g_protein", "g_protein_per_cap_day",
  "Energy (MJ / 100 g protein)", "Time (hr / 100 g protein)", "g protein/cap/day",
  paste0("Energy vs. time to provision 100 g protein (", year, ") — Econ. + non-econ labor"))

p_tradeoff_econlabor = (p_tradeoff_kcal_econlabor | p_tradeoff_protein_econlabor) +
  plot_layout(guides = "collect") & theme(legend.position = "right")

p_tradeoff_allwork = (p_tradeoff_kcal_allwork | p_tradeoff_protein_allwork) +
  plot_layout(guides = "collect") & theme(legend.position = "right")

ggsave(paste0("results/tradeoff_convfac_econlabor.pdf"), p_tradeoff_econlabor, width = 20, height = 8)
ggsave(paste0("results/tradeoff_convfac_allwork.pdf"), p_tradeoff_allwork, width = 20, height = 8)

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
  mutate(hr_per_100g_protein = hr_per_cap_day / g_protein_per_cap_day * 100,
         mj_per_100g_protein = mj_per_cap_day / g_protein_per_cap_day * 100) %>%
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
  mutate(hr_per_100g_protein = hr_per_cap_day / g_protein_per_cap_day * 100,
         mj_per_100g_protein = mj_per_cap_day / g_protein_per_cap_day * 100) %>%
  left_join(regions %>% select(iso3c, continent), by = c("country" = "iso3c")) %>%
  drop_na() %>%
  mutate(is_row = as.character(country) %in% row_countries)

p_tradeoff_kcal_econlabor_consump = tradeoff_scatter(
  tradeoff_kcal_econlabor_consump, "mj_per_2000kcal", "hr_per_2000kcal", "kcal_per_cap_day",
  "Energy (MJ / 2000 kcal)", "Time (hr / 2000 kcal)", "kcal/cap/day",
  paste0("Energy vs. time to provision 2000 kcal (", year, ") — Economic, consumption-based"))

p_tradeoff_kcal_allwork_consump = tradeoff_scatter(
  tradeoff_kcal_allwork_consump, "mj_per_2000kcal", "hr_per_2000kcal", "kcal_per_cap_day",
  "Energy (MJ / 2000 kcal)", "Time (hr / 2000 kcal)", "kcal/cap/day",
  paste0("Energy vs. time to provision 2000 kcal (", year, ") — Total incl. household, consumption-based"))

p_tradeoff_protein_econlabor_consump = tradeoff_scatter(
  tradeoff_protein_econlabor_consump, "mj_per_100g_protein", "hr_per_100g_protein", "g_protein_per_cap_day",
  "Energy (MJ / 100 g protein)", "Time (hr / 100 g protein)", "g protein/cap/day",
  paste0("Energy vs. time to provision 100 g protein (", year, ") — Economic, consumption-based"))

p_tradeoff_protein_allwork_consump = tradeoff_scatter(
  tradeoff_protein_allwork_consump, "mj_per_100g_protein", "hr_per_100g_protein", "g_protein_per_cap_day",
  "Energy (MJ / 100 g protein)", "Time (hr / 100 g protein)", "g protein/cap/day",
  paste0("Energy vs. time to provision 100 g protein (", year, ") — Total incl. household, consumption-based"))

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
  tradeoff_pcap_econlabor, "mj_per_cap_day", "hr_per_cap_day", "g_protein_per_cap_day",
  "Energy (MJ/cap/day)", "Time (hr/cap/day)", "g protein/cap/day",
  paste0("Energy vs. time per capita (", year, ") — Food, Economic"))

p_tradeoff_pcap_allwork = tradeoff_scatter(
  tradeoff_pcap_allwork, "mj_per_cap_day", "hr_per_cap_day", "g_protein_per_cap_day",
  "Energy (MJ/cap/day)", "Time (hr/cap/day)", "g protein/cap/day",
  paste0("Energy vs. time per capita (", year, ") — Food, Total incl. household"))

p_tradeoff_pcap_nonfood_econlabor = tradeoff_scatter(
  tradeoff_pcap_nonfood_econlabor, "mj_per_cap_day", "hr_per_cap_day", "g_protein_per_cap_day",
  "Energy (MJ/cap/day)", "Time (hr/cap/day)", "g protein/cap/day",
  paste0("Energy vs. time per capita (", year, ") — Non-food, Economic"))

ggsave(paste0("results/tradeoff_pcap_econlabor.pdf"),        p_tradeoff_pcap_econlabor,        width = 14, height = 6)
ggsave(paste0("results/tradeoff_pcap_allwork.pdf"),        p_tradeoff_pcap_allwork,        width = 14, height = 6)
ggsave(paste0("results/tradeoff_pcap_nonfood_econlabor.pdf"),  p_tradeoff_pcap_nonfood_econlabor,  width = 14, height = 6)

# Non-food: normalised by kcal and protein (non-food time/energy per unit nutrition)
tradeoff_kcal_nonfood_econlabor = tradeoff_pcap_nonfood_econlabor %>%
  inner_join(kcal_domestic, by = "country") %>%
  mutate(mj_per_2000kcal = mj_per_cap_day / kcal_per_cap_day * 2000,
         hr_per_2000kcal = hr_per_cap_day  / kcal_per_cap_day * 2000)

tradeoff_protein_nonfood_econlabor = tradeoff_pcap_nonfood_econlabor %>%
  mutate(mj_per_100g_protein = mj_per_cap_day / g_protein_per_cap_day * 100,
         hr_per_100g_protein = hr_per_cap_day  / g_protein_per_cap_day * 100)

p_tradeoff_kcal_nonfood_econlabor = tradeoff_scatter(
  tradeoff_kcal_nonfood_econlabor, "mj_per_2000kcal", "hr_per_2000kcal", "kcal_per_cap_day",
  "Energy (MJ / 2000 kcal)", "Time (hr / 2000 kcal)", "kcal/cap/day",
  paste0("Energy vs. time per 2000 kcal (", year, ") — Non-food, Economic"))

p_tradeoff_protein_nonfood_econlabor = tradeoff_scatter(
  tradeoff_protein_nonfood_econlabor, "mj_per_100g_protein", "hr_per_100g_protein", "g_protein_per_cap_day",
  "Energy (MJ / 100 g protein)", "Time (hr / 100 g protein)", "g protein/cap/day",
  paste0("Energy vs. time per 100 g protein (", year, ") — Non-food, Economic"))

p_tradeoff_nonfood_econlabor = (p_tradeoff_kcal_nonfood_econlabor | p_tradeoff_protein_nonfood_econlabor) +
  plot_layout(guides = "collect") & theme(legend.position = "right")

ggsave(paste0("results/tradeoff_convfac_nonfood_econlabor.pdf"), p_tradeoff_nonfood_econlabor, width = 20, height = 8)


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
  cont  <- ifelse(other,
                  sub("^Other ", "", labels),
                  regions$continent[match(labels, regions$iso3c)])
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
# Left nodes (prod): country where labor is expended in the food supply chain.
# Right nodes (cons): country whose food consumption demands that labor.
for (sector in c("food", "nonfood")) {
  l_country <- if (sector == "food") l_food_country else l_nonfood_country

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

    res <- agg_to_country_sankey(mat_hr)
    # M.hr matrix: * 1e6 converts to hr, * 60 converts to min → pcap_scale = 6e7
    sk  <- mat_to_sankey(res$mat, scale = 1e3, pop = res$pop, pcap_label = pcap_label, digits = 2,
                         pcap_scale = 6e7)
    p   <- sankeyNetwork(
      Links = sk$links, Nodes = sk$nodes,
      Source = "source", Target = "target", Value = "value", NodeID = "name",
      NodeGroup = "group", colourScale = sankey_colour_scale,
      sinksRight = FALSE, fontSize = 13, nodeWidth = 20, nodePadding = 10,
      units = "Ghr", iterations = 0
    ) %>% add_pcap_tooltip(sk) %>% add_continent_legend()

    htmlwidgets::saveWidget(p, paste0("results/sankey_hr_", sector, "_", hr_label, ".html"),
                            selfcontained = FALSE)
  }
}

# Sankeys for energy (TJ→EJ, ÷1e6) by food/non-food sector
for (sector in c("food", "nonfood")) {
  l_country = if (sector == "food") l_food_country else l_nonfood_country
  res = agg_to_country_sankey(l_country[["en"]])
  sk  = mat_to_sankey(res$mat, scale = 1e6)
  p   = sankeyNetwork(
    Links = sk$links, Nodes = sk$nodes,
    Source = "source", Target = "target", Value = "value", NodeID = "name",
    NodeGroup = "group", colourScale = sankey_colour_scale,
    sinksRight = FALSE, fontSize = 13, nodeWidth = 20, nodePadding = 10,
    units = "EJ", iterations = 0
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
  filter(!iso3c %in% c("ROW", "TWN", "ANT")) %>% # countries with no population data
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
  filter(!iso3c %in% c("ROW", "TWN", "ANT")) %>% # countries with no population data
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


