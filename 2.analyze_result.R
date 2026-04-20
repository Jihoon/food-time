#### Result analysis ####

# Population
library(gt)
data(countrypops)

#### 1. Energy and Labor footprint ####

# Total consumption-based footprint per country (23001 x 187)
# By taking Diagonal(), we preserve the origin information in the intensity vectors. Target info is lost.
# which can be good enough.


# food-sector footprints
l_int_d = readRDS(file = paste0("data/FABIO_exio_satellites_food_", year, ".rds"))
l_int_i = readRDS(file = paste0("data/FABIO_exio_satellites_nonfood_", year, ".rds"))

# Footprint summed at the FABIO country level
fp_food <- lapply(l_int_d, function(d) Matrix::Diagonal(x=d) %*% FABIO_x_hh)
# fp_nonfood <- lapply(l_int_i, function(d) d %*% FABIO_x_hh)
n_nf = length(exio_nonfood_sectors) # EXIO non-food sectors
# fp_nonfood <- lapply(l_int_i, function(d) {
#   fp <- matrix(0, nrow = n_nf * nrreg, ncol = nrreg)
#   fp <- as(fp, "CsparseMatrix")
#   B_i <- matrix(0, nrcom*nrreg, nrreg)
#   for (i in 1:nrreg) {
#     A_i <- d[((i-1)*n_nf + 1):(i*n_nf), ]
#     for (k in 1:nrreg) {
#       B_i[((k-1)*nrcom + 1):(k*nrcom), k] <- FABIO_x_hh[((k-1)*nrcom + 1):(k*nrcom), i]
#     }
#     fp[((i-1)*n_nf + 1):(i*n_nf), ] <- A_i %*% B_i
#   }
#   fp
# })

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

# Make an index for FABIO region - EXIO nonfood sectors
io_nonfood = expand.grid(sect = exio_nonfood_sectors, iso3c = regions$iso3c) %>% select(iso3c, sect)

# Validate: find the biggest cell from fp_food[[2]] and check if it corresponds to expected source and destination
largest_cell_index = which(fp_food[[2]] == max(fp_food[[2]]), arr.ind = TRUE)
print(paste("Largest cell in fp_food[[2]] is at row", largest_cell_index[1], "and column", largest_cell_index[2]))
print(paste("This corresponds to country", io[largest_cell_index[1],c("iso3c", "item")], "and sector", regions$area[largest_cell_index[2]]))

largest_cell_index = which(fp_nonfood[[2]] == max(fp_nonfood[[2]]), arr.ind = TRUE)
print(paste("Largest cell in fp_food[[2]] is at row", largest_cell_index[1], "and column", largest_cell_index[2]))
print(paste0("This corresponds to country ", io_nonfood[largest_cell_index[1],"iso3c"], "'s '", 
             io_nonfood[largest_cell_index[1],"sect"], 
             "' going to ", regions$area[largest_cell_index[2]]))

# Debug: Find a country block for KIR and get the ten biggest cell values in fp_food[[2]] for that block and show their product and destination countries
kir_block_indices = which(io$iso3c == "KIR")
kir_block_values = fp_food[[2]][kir_block_indices, ]
kir_block_largest_indices = order(kir_block_values, decreasing = TRUE)[1:10]
print("Top 10 largest cells in KIR block of fp_food[[2]]:")
for (index in kir_block_largest_indices) {
  row_col = arrayInd(index, dim(kir_block_values))
  print(paste("Cell at row", row_col[1], "and column", row_col[2], "with value", kir_block_values[index]))
  print(paste("This corresponds to product", io[kir_block_indices[row_col[1]], c("iso3c", "item")], "and destination country", regions$area[row_col[2]]))
}

# Debug: Find a country block for ISL and get the ten biggest cell values in fp_food[[2]] for that block and show their product and destination countries
isl_block_indices = which(io$iso3c == "ISL")
isl_block_values = fp_food[[2]][isl_block_indices, ]
isl_block_largest_indices = order(isl_block_values, decreasing = TRUE)[1:10]
print("Top 10 largest cells in ISL block of fp_food[[2]]:")
for (index in isl_block_largest_indices) {
  row_col = arrayInd(index, dim(isl_block_values))
  print(paste("Cell at row", row_col[1], "and column", row_col[2], "with value", isl_block_values[index]))
  print(paste("This corresponds to product", io[isl_block_indices[row_col[1]], c("iso3c", "item")], "and destination country", regions$area[row_col[2]]))
}

# Debug: Find a country block for KIR and get the ten biggest cell values in fp_nonfood[[2]] for that block and show their product and destination countries
kir_block_indices = which(io_nonfood$iso3c == "KIR")
kir_block_values = fp_nonfood[[2]][kir_block_indices, ]
kir_block_largest_indices = order(kir_block_values, decreasing = TRUE)[1:10]
print("Top 10 largest cells in KIR block of fp_nonfood[[2]]:")
for (index in kir_block_largest_indices) {
  row_col = arrayInd(index, dim(kir_block_values))
  print(paste("Cell at row", row_col[1], "and column", row_col[2], "with value", kir_block_values[index]))
  print(paste("This corresponds to sector", io_nonfood[kir_block_indices[row_col[1]], c("sect")], "and destination country", regions$area[row_col[2]]))
}



# Aggregate all at country level ####
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



# Domestic/import/export summary by country ####

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
    drop_na() 
summary_nonfood_df = bind_rows(lapply(names(summary_nonfood), 
                                      function(x) summary_nonfood[[x]] %>% mutate(type = x)))%>%
    drop_na() 

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
  mutate(country = factor(country, levels = sum_ord),
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
p1 = plot_countries(a %>% filter(type %in% c("hr_f")), "Female time footprint per capita (hr/day)", "")+ ylim(-1, 3.5) 
p2 = plot_countries(a %>% filter(type %in% c("hr_m")), "Male time footprint per capita (hr/day)", "")+ ylim(-1, 3.5) 
p3 = plot_countries(a %>% filter(type %in% c("en")), "Energy footprint per capita (GJ/yr)", "")

p_combined_partial  = p1 / p2 + plot_layout(guides = "collect") & theme(legend.position = "top",
                                                                        axis.text.x = element_blank(),
                                                                        axis.ticks.x = element_blank()
)
p_combined_partial[[2]] <- p_combined_partial[[2]] +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.ticks.x = element_line())
print(p_combined_partial)

# Save it to a PDF
ggsave(paste0("results/footprint_all_countries.pdf"), p_combined, width = 18, height = 12)
ggsave(paste0("results/footprint_partial_countries.pdf"), p_combined_partial, width = 18, height = 12)





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
df_convfac_kcal_econ = summary_time_kcal %>% filter(grepl("per_capita", footprint_type_time)) 

df_convfac_kcal_nonecon = summary_time_kcal %>% filter(country %in% cty_ghd) 

# Plot the distribution of hr/kcal conversion factors by country and by type (domestic, export, import) with different colors for type and different facets for hr_m and hr_f.
# Order countries by domestic_hr_per_2000kcal of hr_f
v_ord_econ = (df_convfac_kcal_econ %>%   
           filter(type %in% c("hr_f"), cat=="domestic") %>% 
           # group_by(country) %>% 
           arrange(-hr_per_2000kcal))$country
v_ord_all = (df_convfac_kcal_nonecon %>% 
               filter(type %in% c("hr_f"), cat=="domestic") %>% 
               group_by(country) %>% summarise(d = sum(hr_per_2000kcal, na.rm=TRUE)) %>%
               arrange(-d))$country

# Plot only the economic conversion factors first, and then add the non-economic ones as a separate plot. This way we can see the difference between the two and also avoid having too many bars in one plot.
p_conversion_kcal_econ = ggplot(df_convfac_kcal_econ %>% 
                                  filter(footprint_type_time=="domestic_per_capita") %>%
                                  select(country, type, footprint_type_time, hr_per_2000kcal) %>%
                                  mutate(country = factor(country, levels = v_ord_econ)),
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
                                mutate(country = factor(country, levels = v_ord_all)),
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

v_ord.tot = (summary_time_kcal %>% drop_na() %>%
           filter(type %in% c("hr_f")) %>% 
           # group_by(country) %>% 
           arrange(-domestic_hr_per_2000kcal.tot))$country
p_conversion_kcal = ggplot(summary_time_kcal %>% drop_na() %>%
                             select(country, type, domestic_hr_per_2000kcal.tot) %>%
                             pivot_longer(cols = starts_with(c("domestic", "export", "import")), 
                                          names_to = "footprint_type", values_to = "hr_per_2000kcal") %>% 
                             mutate(country = factor(country, levels = v_ord.tot)),
                           aes(x=country, y=hr_per_2000kcal, fill=footprint_type)) +
  geom_bar(stat="identity", position="stack") +
  facet_wrap(~type, ncol=1, scales = "fixed") +
  labs(x="Country (ISO3)", y="Time per 2000 kcal (hr/2000kcal)", fill="Footprint type",
       title=paste0("Time per 2000 kcal conversion factors by country (", year, ")")) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_fill_manual(values=c("domestic_hr_per_2000kcal"="#1f77b4", "export_hr_per_2000kcal"="#2ca02c", "import_hr_per_2000kcal"="#ff7f0e")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))      


# Do the same for protein
summary_time_protein = summary_food_df_long %>% 
  filter(type %in% c("hr_m", "hr_f")) %>%
  select(country, type, footprint_type, per_capita_value) %>%
  pivot_wider(names_from = footprint_type, values_from = per_capita_value) %>%
  left_join(summary_pro_df_long %>% select(country, footprint_type, per_capita_value) %>%
              pivot_wider(names_from = footprint_type, values_from = per_capita_value), 
            by = "country", suffix = c("_time", "_protein")) %>%
  mutate(across(ends_with(c("per_capita_time", "per_capita_protein")), ~ .x / 1e6)) %>%
  mutate(domestic_hr_per_g_protein = domestic_per_capita_time / domestic_per_capita_protein,
         export_hr_per_g_protein = export_per_capita_time / export_per_capita_protein,
         import_hr_per_g_protein = import_per_capita_time / import_per_capita_protein)
p_conversion_protein = plot_countries(summary_time_protein %>% 
                                select(country, type, starts_with(c("domestic", "export", "import"))) %>%
                                pivot_longer(cols = starts_with(c("domestic", "export", "import")), 
                                             names_to = "footprint_type", values_to = "hr_per_g_protein") %>%
                                mutate(footprint_type = factor(footprint_type, levels = c("domestic_hr_per_g_protein", "export_hr_per_g_protein", "import_hr_per_g_protein"))),
                              "Time per g protein (hr/g)", "Protein time conversion factors")
         



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
  



# energy_fp = Matrix::Diagonal(x=FABIO_en_int_d) %*% FABIO_x_hh #TJ = EJ/10^6 (sum = 22.5 EJ)
# # Note: "In the United States, food production uses 10.11 quadrillion Btu annually" = 10.7 EJ
# hr_m_fp = Matrix::Diagonal(x=FABIO_hr_m_int_d) %*% FABIO_x_hh #M.hr
# hr_f_fp = Matrix::Diagonal(x=FABIO_hr_f_int_d) %*% FABIO_x_hh #M.hr 
# 
# colnames(energy_fp) = colnames(hr_m_fp) = colnames(hr_f_fp) = regions$iso3c

# Country-wise consumption-based footprint 

consumption = "food" # FABIO FD category
country = "KOR"
extension = "hr_m_int_d" #"en_int_d", "hr_f_int_d"

for (country in regions$iso3c) {
  
  for (extension in names(l_int)) {
    print(paste("Calculating", extension, "footprint for", country))
    
    Y_country <- FABIO_y[, which(fd$iso3c == country)]
    # Y_country <- Yi[, fd$iso3c == country]
    colnames(Y_country) <- fd$fd[fd$iso3c == country]
    
    pop = subset(countrypops, country_code_3 == country & year == yr)$population
    if (country %in% setdiff(regions$iso3c, unique(countrypops$country_code_3))) {pop=NA}
    
    print(paste("Population =", pop))
    
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
    
    int <- l_int[[extension]]
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
    
    cal_consum = sum(country_consump*coeff_cal$kcal_per_kg)/365/pop*1000 #kcal/cap/day
    pro_consum = sum(country_consump*coeff_pro$protein_per_kg)/365/pop*1000 #kcal/cap/day
    
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
    
    data.table::fwrite(data_tot, file=paste0("output/FABIO_", country,"_", year, "_", extension, "_", consumption,".csv"), sep=",")
    
    data_imp_tot = tail(data_tot, 1) %>% rename(importer = item_target) %>% mutate(importer = country)
    
    # Fill mat with data_imp_tot where rownames(mat) match names(data_imp_tot), and colnames(mat) match data_imp_tot$importer
    mat[rownames(mat) %in% names(data_imp_tot), data_imp_tot$importer[1]] <- 
      as.numeric(data_imp_tot[1, names(data_imp_tot) %in% rownames(mat)])
  }
}

# Total Mcal trade between countries
saveRDS(mat, file.path("data/calorie_trade_mat.rds"))



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

pop_y = subset(countrypops, year == yr) %>%
  select(iso3c = country_code_3, pop = population)

pop_y = pop_y %>% right_join(regions, by = "iso3c") %>%
  select(iso3c, pop)

cty = "USA"
imp = round(mat_y[,cty] / pop_y$pop[pop_y$iso3c==cty] / 365, 4)
ex = round(mat_y[cty,] / pop_y$pop[pop_y$iso3c==cty] / 365, 4)

# Domestic consumption
print(paste("Total consumption in", cty, ":", sum(imp, na.rm=TRUE), "kcal/capita/day"))
# Tot imports
print(paste("Total imports in", cty, ":", sum(imp[names(imp)!=cty], na.rm=TRUE), "kcal/capita/day"))
# Tot exports
print(paste("Total exports from", cty, ":", sum(ex[names(ex)!=cty], na.rm=TRUE), "kcal/capita/day"))



######### Total footprint analysis #################

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


