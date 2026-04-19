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
  
  # Finally, calculate FABIO-based EXIO energy and labor use matrices for agri-food sects
  FABIO_en = FABIO_x_in_EXIO %*% Diagonal(x=exio_en_intm_ordered)
  FABIO_hr_m = FABIO_x_in_EXIO %*% Diagonal(x=exio_hr_m_intm_ordered)
  FABIO_hr_f = FABIO_x_in_EXIO %*% Diagonal(x=exio_hr_f_intm_ordered)
  # This chuck doesn't work for (indirect) matrix satellites, since Diagonal() expects a vector.
  
  # TODO: FABIO_x_in_EXIO need to be checked again for computing correct order, 
  # but looks ok for now since this will be run after convert_mass_vecs() is called.
  
  return(list(sat_en_FAB = FABIO_en, 
              sat_hr_m_FAB = FABIO_hr_m, 
              sat_hr_f_FAB = FABIO_hr_f))
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
  footprint_types = unique(df$footprint_type)
  name_dom = footprint_types[grepl("domestic", footprint_types)]
  name_exp = footprint_types[grepl("export", footprint_types)]
  name_imp = footprint_types[grepl("import", footprint_types)]
  # print(name_dom, name_exp, name_imp)
  print(footprint_types)

  c_scheme = c(name_dom="#1f77b4", name_exp="#2ca02c", name_imp="#ff7f0e")
  # Check if the first row of df has type starting with "hr_m" or "hr_f" to determine if it's labor or energy footprint
  if (!"type" %in% colnames(df)) { # Nutrient
    part_negative = name_exp
    scale_factor = 1
  } else if (df$type[1] %in% c("hr_m", "hr_f")) { # Labor footprint
    part_negative = "import_per_capita"
    scale_factor = 1
    c_scheme = c(c_scheme,
                 "preparation_econ" = "#8610ca",
                 "preparation_non.econ" = "#ffa6a6",
                 "processing_non.econ" = "#fc4a4a",
                 "growth_collection_non.econ" = "#ce0303")
    print("Plotting labor footprint with negative import values")
  } else {  # Energy footprint or other types of footprints, default to negative import values
    part_negative = name_exp
    scale_factor = 1e3
  }

  # Get first work before "_" of part_negative to determine the type of footprint for labeling
  neg_type = strsplit(part_negative, "_")[[1]][1]
  pos_type = ifelse(neg_type == "import", "export", "import")

  g = ggplot(df %>%
               filter(footprint_type != part_negative),
             aes(x=country, y=per_capita_value/scale_factor, fill=footprint_type)) +
    # When bars are stacked, make sure that "domestic_per_capita" is at the bottom and what's on top is determined by pos_type.
    geom_bar(stat="identity", position="stack") +
    labs(x="Country (ISO3)", y=ylabel, fill="Footprint type") +
    theme_minimal() +
    theme(legend.position = "top") +
    # Tilt x-axis labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # Add import_per_capita values as negative y axis bars by footprint_type and type
    geom_bar(data = df %>%
               filter(footprint_type == part_negative),
             aes(x=country, y=-per_capita_value/scale_factor, fill=footprint_type),
             stat="identity", position="stack") +
    scale_fill_manual(values=c_scheme) #+
    # labs(fill="Footprint type", title=paste0(maintitle,"\n(positive: domestic+", pos_type, ", negative: ", neg_type, ")"))
    labs(fill="Footprint type")

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


# Reusable stacked bar chart: country × kcal breakdown
plot_country_bars <- function(data, cols, labels, title,
                               sort_desc = TRUE, n = 25,
                               ylim_vals = NULL, filename = NULL,
                               output_dir = NULL) {
  d <- data %>%
    { if (sort_desc) arrange(., desc(supp_pday)) else arrange(., supp_pday) } %>%
    { if (!is.null(n)) slice_head(., n = n) else . } %>%
    mutate(iso3c = factor(iso3c, levels = iso3c[order(supp_pday, decreasing = TRUE)])) %>%
    pivot_longer(cols = all_of(cols), names_to = "type", values_to = "kcal_pday") %>%
    mutate(type = factor(type, levels = cols, labels = labels))

  p <- ggplot(d, aes(x = iso3c, y = kcal_pday, fill = type)) +
    geom_bar(stat = "identity") +
    labs(x = "Country (ISO3)", y = "kcal/capita/day", fill = "Type", title = title) +
    theme_minimal() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = setNames(c("#1f77b4", "#2ca02c", "#ff7f0e"), labels))

  if (!is.null(ylim_vals)) p <- p + ylim(ylim_vals[1], ylim_vals[2])

  print(p)

  if (!is.null(filename) && !is.null(output_dir))
    ggsave(file.path(output_dir, filename), width = 10, height = 6)

  invisible(p)
}
