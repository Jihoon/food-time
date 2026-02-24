#### Result analysis ####


#### 1. Energy and Labor footprint ####

# Total consumption-based footprint per country (23001 x 187)
# By taking Diagonal(), we preserve the origin information in the intensity vectors. Target info is lost.
# which can be good enough.

fp_food <- lapply(l_int_d, function(d) Matrix::Diagonal(x=d) %*% FABIO_x_hh)
fp_nonfood <- lapply(l_int_i, function(d) d %*% FABIO_x_hh)

# energy_fp = Matrix::Diagonal(x=FABIO_en_int_d) %*% FABIO_x_hh #TJ = EJ/10^6 (sum = 22.5 EJ)
# # Note: "In the United States, food production uses 10.11 quadrillion Btu annually" = 10.7 EJ
# hr_m_fp = Matrix::Diagonal(x=FABIO_hr_m_int_d) %*% FABIO_x_hh #M.hr
# hr_f_fp = Matrix::Diagonal(x=FABIO_hr_f_int_d) %*% FABIO_x_hh #M.hr 
# 
# colnames(energy_fp) = colnames(hr_m_fp) = colnames(hr_f_fp) = regions$iso3c

# Country-wise consumption-based footprint 
library(gt)

# Direct food-sector footprint
l_int = readRDS(file = paste0("data/FABIO_exio_satellites_food_", year, ".rds"))
# l_int = readRDS(file = paste0("data/FABIO_exio_satellites_nonfood_", year, ".rds"))

consumption = "food"
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

# Population
library(gt)
data(countrypops)

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



