# Re-order the EXIO intensity vectors (9800) to match FABIO region order based on FABIO_reg mapping
reorder_countries_to_FABIO <- function(int_vector, direction=2) {
  # Change just the region ordering from EXIO to FABIO
  # Wx regions' values simply pasted to all FABIO countries within the region 
  # direction: 1 for row-wise, 2 for column-wise (default)
  
  n_sect = length(EXIO_sect) 
  n_countries_EXIO = n_reg_EXIO
  n_countries_FABIO = nrreg
  
  # # Make an index vector by repeating 200 times
  # exio_reg_idx = rep(1:n_countries, each=n_sect)
  
  # Check if any dimensional length of int_vector is multiples of n_countries_EXIO
  if (!any(dim(int_vector) %% n_countries_EXIO == 0)) {
    stop("Input vector length must be a multiple of 49 EXIO regions.")
  }
  
  # Find that multiple to determine the number of sectors (n_sect) in the input vector
  if (direction == 1) {
    n_sect = dim(int_vector)[1] / n_countries_EXIO 
    int_vector = t(int_vector) # Transpose to make it column-wise for easier processing
  } 
  print(paste("Input vector has", n_sect, "sectors per EXIO region."))
  
  # If input is a vector of length 9800
  if (min(dim(int_vector))<2) {
    M = matrix(int_vector, nrow = n_countries_EXIO, ncol = n_sect, byrow = TRUE)
    M_reordered <- M[FABIO_reg$EXIOBASE_code, , drop = FALSE]
    
    # Flatten back to original vector shape
    data_reordered <- as.vector(t(M_reordered)) # Length 37400 = 187 countries * 200 sectors
  }
  # If input is a matrix with 9800 columns, reorder each block of 200 columns
  else {
    n_rows = dim(int_vector)[1]
    data_reordered = matrix(0, nrow = n_rows, ncol = n_countries_FABIO * n_sect)
    
    for (i in 1:n_rows) {
      M = matrix(int_vector[i, ], nrow = n_countries_EXIO, ncol = n_sect, byrow = TRUE)
      M_reordered <- M[FABIO_reg$EXIOBASE_code, , drop = FALSE]
      
      # Flatten back to original vector shape
      data_reordered[i, ] <- as.vector(t(M_reordered)) # Length 37400 = 187 countries * 200 sectors
    }

    data_reordered <- as(data_reordered, "CsparseMatrix") # Convert to sparse matrix format to save memory
  }
  
  if (direction == 1) {
    data_reordered = t(data_reordered) # Transpose back to original orientation
  }
  
  return(data_reordered)
}


# Re-order the FABIO mass X vector (187x200=37400) to match EXIO region order based on FABIO_reg mapping
reorder_countries_to_EXIO <- function(m_vector) {
  n_sect = 200
  n_countries = 187
  
  # # Make an index vector by repeating 200 times
  # exio_reg_idx = rep(1:n_countries, each=n_sect)
  r = FABIO_reg %>% arrange(EXIOBASE_code)
  r_idx = r %>% row_number()
  
  M = matrix(m_vector, nrow = n_countries, ncol = n_sect, byrow = TRUE)
  M_reordered = M[r_idx, , drop = FALSE]
  
  # Aggregate rows by EXIO regions
  M_exio = rowsum(M_reordered, group = r$EXIOBASE_code)
  
  # Flatten back to original vector shape
  data_reordered <- as.vector(t(M_exio)) # mass vector (length 9800) to be used for intensity calculation (J/kg)
  
  return(data_reordered)
}


# Function to convert EXIO intensity vectors to FABIO-based matrices
convert_intensities <- function(sat_EXIO) {
  
  # Intensity by mass (J/kg or hr/kg) by EXIO sec/reg (len = 9800)
  exio_en_int = sat_EXIO$sat_en / exio_mass_x
  exio_hr_m_int = sat_EXIO$sat_hr_m / exio_mass_x
  exio_hr_f_int = sat_EXIO$sat_hr_f / exio_mass_x
  
  # Replace NaN and Inf with 0
  exio_en_int[!is.finite(exio_en_int)] = 0
  exio_hr_m_int[!is.finite(exio_hr_m_int)] = 0
  exio_hr_f_int[!is.finite(exio_hr_f_int)] = 0
  
  # Re-order intensity vectors to match FABIO region order (len = 37400)
  # Wx regions' values simply pasted to all FABIO countries within the region
  exio_en_intm_ordered = reorder_countries_to_FABIO(exio_en_int)
  exio_hr_m_intm_ordered = reorder_countries_to_FABIO(exio_hr_m_int)
  exio_hr_f_intm_ordered = reorder_countries_to_FABIO(exio_hr_f_int)
  
  # Debugging: Check the top intensities of fish sectors in the ordered intensity vectors to see if they correspond to expected EXIO sectors and countries
  fish_exio_hr_m_int = exio_hr_m_intm_ordered[fish_idx_l]
  # top_fish_exio_hr_m_int = sort(fish_exio_hr_m_int, decreasing = TRUE)[1:10]
  idx = order(fish_exio_hr_m_int, decreasing = TRUE)[1:20]
  cat("Top EXIO sectors receiving fish (male) labor intensity (per ton):\n")
  print(data.frame(sector = long_idx[fish_idx_l[idx]], intensity = fish_exio_hr_m_int[idx]))

  fish_exio_hr_f_int = exio_hr_f_intm_ordered[fish_idx_l]
  # top_fish_exio_hr_f_int = sort(fish_exio_hr_f_int, decreasing = TRUE)[1:10]
  idx = order(fish_exio_hr_f_int, decreasing = TRUE)[1:20]
  cat("Top EXIO sectors receiving fish (female) labor intensity (per ton):\n")
  print(data.frame(sector = long_idx[fish_idx_l[idx]], intensity = fish_exio_hr_f_int[idx]))

  # Finally, calculate FABIO-based EXIO energy and labor use matrices for agri-food sects
  FABIO_en = FABIO_x_in_EXIO %*% Diagonal(x=exio_en_intm_ordered)
  FABIO_hr_m = FABIO_x_in_EXIO %*% Diagonal(x=exio_hr_m_intm_ordered)
  FABIO_hr_f = FABIO_x_in_EXIO %*% Diagonal(x=exio_hr_f_intm_ordered)
  # This chunk doesn't work for (indirect) matrix satellites, since Diagonal() expects a vector.
  
  # TODO: FABIO_x_in_EXIO need to be checked again for computing correct order, 
  # but looks ok for now since this will be run after convert_mass_vecs() is called.
  
  # Debugging: Check the top values in the labor hours for FABIO fish sectors.
  fish_fabio_hr_m = rowSums(FABIO_hr_m[fish_idx, ])
  idx = order(fish_fabio_hr_m, decreasing = TRUE)[1:20]
  cat("Top FABIO sectors receiving fish (male) labor hr:\n")
  print(data.frame(iso3c = io$iso3c[fish_idx][idx], sector = io$item[fish_idx][idx], labor_m = fish_fabio_hr_m[idx]))

  fish_fabio_hr_f = rowSums(FABIO_hr_f[fish_idx, ])
  idx = order(fish_fabio_hr_f, decreasing = TRUE)[1:20]
  cat("Top FABIO sectors receiving fish (female) labor hr:\n")
  print(data.frame(iso3c = io$iso3c[fish_idx][idx], sector = io$item[fish_idx][idx], labor_f = fish_fabio_hr_f[idx]))

  # Derive per-capita values for the top fish sectors in FABIO to check if they correspond to expected countries with high fish consumption/production.
  fish_fabio_hr_m_cap = fish_fabio_hr_m / pop_y$pop[match(io$iso3c[fish_idx], pop_y$iso3c)] *2
  idx = order(fish_fabio_hr_m_cap, decreasing = TRUE)[1:20]
  cat("Top FABIO sectors receiving fish (male) labor hr per capita:\n")
  print(data.frame(iso3c = io$iso3c[fish_idx][idx], sector = io$item[fish_idx][idx], labor_m_cap = fish_fabio_hr_m_cap[idx]))

  fish_fabio_hr_f_cap = fish_fabio_hr_f / pop_y$pop[match(io$iso3c[fish_idx], pop_y$iso3c)] *2 # Assume 1:1 gender ratio
  idx = order(fish_fabio_hr_f_cap, decreasing = TRUE)[1:20]
  cat("Top FABIO sectors receiving fish (female) labor hr per capita:\n")
  print(data.frame(iso3c = io$iso3c[fish_idx][idx], sector = io$item[fish_idx][idx], labor_f_cap = fish_fabio_hr_f_cap[idx]))


  return(list(sat_en_FAB = FABIO_en,
              sat_hr_m_FAB = FABIO_hr_m,
              sat_hr_f_FAB = FABIO_hr_f))
}

# Claude made this one below to fix the non-food footprint being too large. 
# But I won't use it since I made the fix myself.

# Inverse of convert_mass_vecs: maps an EXIO mass vector/matrix back to FABIO classification.
# Each EXIO sector's mass is distributed to FABIO products proportionally to their
# contribution recorded in FABIO_x_in_EXIO (built by convert_mass_vecs).
# Total mass is preserved for EXIO sectors that have FABIO mappings (food sectors).
#
# Input:  exio_vec_mass — vector (37400 or 9800) or matrix with 37400 or 9800 columns
# Output: FABIO mass vector (23001) or matrix with 23001 columns
convert_mass_vecs_inv <- function(exio_vec_mass, FABIO_x_in_EXIO_mat = FABIO_x_in_EXIO) {
  n_exio_full = nrreg * length(EXIO_sect)       # 37400 = 187 countries × 200 sectors
  n_exio_agg  = n_reg_EXIO * length(EXIO_sect)  # 9800  =  49 regions  × 200 sectors

  # Coerce vectors to 1-row matrix so the rest of the logic is uniform
  is_vec = is.null(dim(exio_vec_mass))
  if (is_vec) exio_vec_mass = matrix(exio_vec_mass, nrow = 1)

  n_cols = ncol(exio_vec_mass)
  if (n_cols == n_exio_agg) {
    exio_vec_mass = reorder_countries_to_FABIO(exio_vec_mass)  # K×9800 → K×37400
  } else if (n_cols != n_exio_full) {
    stop(paste("Input column count must be", n_exio_full, "(37400) or", n_exio_agg, "(9800)."))
  }

  # Column-normalize FABIO_x_in_EXIO: share of each FABIO product in each EXIO sector
  col_inv = 1 / colSums(FABIO_x_in_EXIO_mat)
  col_inv[!is.finite(col_inv)] = 0
  col_shares = FABIO_x_in_EXIO_mat %*% Diagonal(x = col_inv)  # 23001×37400, columns sum to 1

  result = exio_vec_mass %*% t(col_shares)  # K×37400 %*% 37400×23001 = K×23001

  if (is_vec) as.vector(result) else result
}


# Util functions for summary
get_mat_summary = function(mat) {
  # Summarise the matrix
  prod = rowSums(mat)
  cons = colSums(mat)
  dom = diag(mat)
  exp = prod - dom
  imp = cons - dom
  
  mat_summary = data.frame(iso3c = names(prod), prod, cons, dom, exp, imp) %>%
    mutate(net = exp - imp,
           dom_perc = dom / cons * 100,
           imp_perc = imp / cons * 100,
           exp_perc = exp / prod * 100) %>%
    arrange(desc(prod)) 
  
  subset(countrypops,  year == yr) %>% 
    select(country_code_3, population) %>%
    right_join(mat_summary, by = c("country_code_3" = "iso3c")) %>%
    mutate(prod_pday = prod / population * 1000 / 365, # Mcal to kcal
           imp_pday = imp / population * 1000 / 365,
           exp_pday = exp / population * 1000 / 365,
           cons_pday = (dom+imp) / population * 1000 / 365) %>%
    arrange(desc(prod)) -> mat_summary
  
  return(mat_summary)
}

get_mat_pcap = function(mat) {
  # Divide mat_y columns by population of matching countries
  pop_missing = c("ANT", "TWN", "ROW")
  # mat_Y_cap = mat_y[!rownames(mat_y) %in% pop_missing, !colnames(mat_y) %in% pop_missing] /
  # pop_y$pop[match(colnames(mat_y)[!colnames(mat_y) %in% pop_missing], pop_y$iso3c)]
  mat_cap = t(t(mat) / pop_y$pop[match(colnames(mat), pop_y$iso3c)])* 1000 / 365
  
  return(mat_cap)  
}


get_net_trade <- function(mat) {
  if (!is.matrix(mat) || nrow(mat) != ncol(mat)) {
    stop("Input must be a square matrix")
  }
  # Create masks for upper and lower triangles (excluding diagonal)
  upper <- mat
  upper[lower.tri(mat, diag = TRUE)] <- 0
  lower <- mat
  lower[upper.tri(mat, diag = TRUE)] <- 0
  
  # Subtract upper triangle from lower triangle
  result1 <- lower - t(upper)
  result2 <- upper - t(lower)
  result = result1+result2
  result[result<0] <- 0
  
  return(result)
}



# Make a block diagonal matrix having M as diagonals, repeated for each FABIO region
block_diag_repeat <- function(M, n) {
  Msp <- as(M, "dgCMatrix") 
  kronecker(Diagonal(n), Msp)
}


sparsity <- function(mat) {
  nnz <- sum(mat != 0, na.rm = TRUE)
  1 - nnz / length(mat)
}


# Result / Analysis

plot_countries <- function(df, ylabel, maintitle) {

  # Get domestsic/export/import factor names out of df$footprint_type
  footprint_types = as.character(unique(df$footprint_type))
  name_dom = footprint_types[grepl("domestic", footprint_types)]
  name_exp = footprint_types[grepl("export", footprint_types)]
  name_imp = footprint_types[grepl("import", footprint_types)]
  # print(name_dom, name_exp, name_imp)
  print(footprint_types)

  # Food footprint_types are plain ("domestic_per_capita"); their non-food-sector
  # counterparts carry a "_nf" suffix ("domestic_per_capita_nf") so a food and a
  # non-food component can be stacked/negated side by side in the same bar while
  # staying visually distinguishable (blue family = food, green family = non-food).
  color_for_sector = function(names, food_color, nonfood_color) {
    setNames(ifelse(grepl("_nf$", names), nonfood_color, food_color), names)
  }

  # Import categories disaggregated by continent-of-origin ("import_cont_<Continent>",
  # optionally suffixed "_nf" for non-food) get one shade per continent -- a blues
  # ramp for food, a greens ramp for non-food -- instead of the flat single color
  # used for a plain "import_per_capita" total.
  import_continents = unique(gsub("_nf$", "", gsub("^import_cont_", "", name_imp[grepl("^import_cont_", name_imp)])))
  if (length(import_continents) > 0) {
    pal_food_imp    = setNames(colorRampPalette(RColorBrewer::brewer.pal(9, "Blues")[3:9])(length(import_continents)), import_continents)
    pal_nonfood_imp = setNames(colorRampPalette(RColorBrewer::brewer.pal(9, "Greens")[3:9])(length(import_continents)), import_continents)
    imp_scheme = setNames(sapply(name_imp, function(nm) {
      if (!grepl("^import_cont_", nm)) return(if (grepl("_nf$", nm)) "#a1d99b" else "#6baed6")
      is_nf = grepl("_nf$", nm)
      cont  = gsub("_nf$", "", gsub("^import_cont_", "", nm))
      if (is_nf) pal_nonfood_imp[[cont]] else pal_food_imp[[cont]]
    }), name_imp)
  } else {
    imp_scheme = color_for_sector(name_imp, "#6baed6", "#a1d99b")
  }

  c_scheme = c(color_for_sector(name_dom, "#08519c", "#31a354"),
               color_for_sector(name_exp, "#3182bd", "#74c476"),
               imp_scheme)
  # Check if the first row of df has type starting with "hr_m" or "hr_f" to determine if it's labor or energy footprint
  if (!"type" %in% colnames(df)) { # Nutrient
    part_negative = name_exp
    scale_factor = 1
  } else if (df$type[1] %in% c("hr_m", "hr_f")) { # Labor footprint
    part_negative = name_imp # all import_per_capita[_nf] categories are plotted as negative bars, stacked
    scale_factor = 1
    # Econ (paid) categories stay in the cool blue family (set above); non-econ
    # (household/unpaid) categories get their own warm yellow-orange-red family,
    # so the two groups are distinguishable by hue while individual items are
    # still distinguishable by shade within each family.
    c_scheme = c(c_scheme,
                 "preparation_econ" = "#bdd7e7",
                 "preparation_non.econ" = "#fed976",
                 "processing_non.econ" = "#feb24c",
                 "growth_collection_non.econ" = "#fd8d3c",
                 "energy_non.econ" = "#f03b20",
                 "water_non.econ" = "#bd0026")
    print("Plotting labor footprint with negative import values")
  } else {  # Energy footprint or other types of footprints, default to negative import values
    part_negative = name_exp
    scale_factor = 1e3
  }
  
  # Get first work before "_" of part_negative to determine the type of footprint for labeling
  neg_type = strsplit(as.character(part_negative), "_")[[1]][1]
  pos_type = ifelse(neg_type == "import", "export", "import")
  
  has_row <- "is_row" %in% colnames(df)
  pos_df  <- df %>% filter(!footprint_type %in% part_negative)
  neg_df  <- df %>% filter(footprint_type %in% part_negative)

  if (has_row) {
    g <- ggplot(pos_df, aes(x=country, y=per_capita_value/scale_factor, fill=footprint_type, alpha=is_row))
  } else {
    g <- ggplot(pos_df, aes(x=country, y=per_capita_value/scale_factor, fill=footprint_type))
  }

  g <- g +
    geom_bar(stat="identity", position="stack") +
    labs(x="Country (ISO3)", y=ylabel, fill="Footprint type") +
    theme_minimal() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))

  if (has_row) {
    g <- g + geom_bar(data=neg_df,
                      aes(x=country, y=-per_capita_value/scale_factor, fill=footprint_type, alpha=is_row),
                      stat="identity", position="stack")
  } else {
    g <- g + geom_bar(data=neg_df,
                      aes(x=country, y=-per_capita_value/scale_factor, fill=footprint_type),
                      stat="identity", position="stack")
  }

  g <- g + scale_fill_manual(values=c_scheme)

  if (has_row) {
    g <- g + scale_alpha_manual(values=c("TRUE"=0.3, "FALSE"=0.9), guide="none")
  }

  print(g)
  return(g)
}


# Function to plot the top 10 countries with the highest per capita footprint for each type of footprint

find_top_cells <- function(df, n = 10, matrix_name = "the matrix") {
  if (is.data.frame(df)) {
    mat <- as.matrix(df)
  } else if (is.matrix(df)) {
    mat <- df
  } else {
    stop("Input must be a data frame or matrix")
  }
  
  # Top largest cells
  largest_indices <- order(mat, decreasing = TRUE)[1:min(n, length(mat))]
  print(paste("Top", length(largest_indices), "largest cells in", matrix_name, ":"))
  for (idx in largest_indices) {
    row_col <- arrayInd(idx, dim(mat))
    val <- mat[idx]
    row_name <- if (!is.null(rownames(mat))) rownames(mat)[row_col[1]] else row_col[1]
    col_name <- if (!is.null(colnames(mat))) colnames(mat)[row_col[2]] else row_col[2]
    print(paste("Cell at row", row_col[1], "and column", row_col[2], "with value", val))
    print(paste("This corresponds to", row_name, "and", col_name))
  }
  
  # Top off-diagonal, assuming square
  if (nrow(mat) != ncol(mat)) {
    print("Matrix is not square, skipping off-diagonal analysis")
    return()
  }
  off_diag_vec <- mat[!diag(nrow(mat))]
  largest_offdiag_indices <- order(off_diag_vec, decreasing = TRUE)[1:min(n, length(off_diag_vec))]
  print(paste("Top", length(largest_offdiag_indices), "largest off-diagonal cells in", matrix_name, ":"))
  diag_mask <- diag(nrow(mat))
  flat_idx <- which(!diag_mask, arr.ind = FALSE)
  for (i in 1:length(largest_offdiag_indices)) {
    idx_in_flat <- largest_offdiag_indices[i]
    global_idx <- flat_idx[idx_in_flat]
    row_col <- arrayInd(global_idx, dim(mat))
    val <- mat[global_idx]
    row_name <- if (!is.null(rownames(mat))) rownames(mat)[row_col[1]] else row_col[1]
    col_name <- if (!is.null(colnames(mat))) colnames(mat)[row_col[2]] else row_col[2]
    print(paste("Cell at row", row_col[1], "and column", row_col[2], "with value", val))
    print(paste("This corresponds to", row_name, "and", col_name))
  }
}


