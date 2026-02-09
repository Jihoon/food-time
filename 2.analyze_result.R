#### Result analysis ####

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



