# Re-order the EXIO intensity vectors (9800) to match FABIO region order based on FABIO_reg mapping
reorder_countries_to_FABIO <- function(int_vector) {
  # Change just the region ordering from EXIO to FABIO
  # Wx regions' values simply pasted to all FABIO countries within the region 
  n_sect = 200
  n_countries_EXIO = 49
  n_countries_FABIO = 187
  
  # # Make an index vector by repeating 200 times
  # exio_reg_idx = rep(1:n_countries, each=n_sect)
  
  if (!(9800 %in% dim(int_vector))) {
    stop("Input vector length must be 9800 (49 EXIO regions x 200 sectors).")
  }
  
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
