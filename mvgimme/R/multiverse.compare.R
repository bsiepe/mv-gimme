#' Compare Results from Multiverse Analysis
#'
#' This function compares the results obtained from a multiverse analysis with a reference model. It performs three levels of comparison: group-level, subgroup-level, and individual-level comparison. The results are combined into a tibble where each row represents one specification, and individual results can be stored in a dataframe inside the tibble.
#'
#' @param l_res A list of results obtained from the multiverse analysis. Each element of the list should be a data object containing the results.
#' @param ref_model The reference model to compare the results with.
#' @return A tibble containing the comparison results at different levels: group-level, subgroup-level, and individual-level comparison along with condition information for each specification.
#' @export
multiverse.compare <- function(l_res,
                               ref_model
){
  
  
  #---- Compare Group-Level
  comp_group <- multiverse.compare.group(
    l_res = l_res, ref_model = ref_model
  )
  
  #---- Compare Subgroup-Level
  comp_subgroup <- multiverse.compare.subgroup(
    l_res = l_res, ref_model = ref_model
  )  
  
  #---- Compare Individual Solutions
  comp_ind <- multiverse.compare.individual(
    l_res = l_res, ref_model = ref_model
  )   
  
  #--- Output
  # Idea: tibble where each row is one specification
  # individual results can be stored in dataframe inside tibble
  comp_res <- dplyr::bind_cols(comp_group, comp_subgroup, comp_ind)
  
  
  # save condition for each specification 
  conds <- do.call(rbind, lapply(l_res, `[[`, "conds"))
  comp_res <- cbind(comp_res, conds)
  
  return(comp_res)
}


multiverse.compare.group <- function(l_res,
                                     ref_model){
  
  #--- Reference model info
  # TODO this will have to account for nonconvergence
  n_ind <- length(ref_model$data)
  n_var <- ref_model$n_vars_total
  
  # indices for temporal and contemporaneous
  temp_ind <- 1:(n_var/2)
  cont_ind <- ((n_var/2)+1):n_var
  
  #--- Adjacency Matrix
  ## Find group effects adjacency matrix
  ref_groupedge <- ifelse(ref_model$path_counts == n_ind, 1, 0)
  # ignore autoregressive coefs
  diag(ref_groupedge[, temp_ind]) <- rep(0, n_var)
  
  # Find group effects adjacency matrix differences
  l_ref_diff <- lapply(l_res, function(x){
    tmp_groupedge <- ifelse(x$path_counts == n_ind, 1, 0)
    diag(tmp_groupedge[, temp_ind]) <- rep(0, n_var)
    diff_groupedge <- ref_groupedge - tmp_groupedge
    return(diff_groupedge)
  }
  )
  
  # Count occurrence of each group effect
  ## list of adjacency matrices
  l_adjacency <- lapply(l_res, function(x){
    tmp_groupedge <- ifelse(x$path_counts == n_ind, 1, 0)
    return(tmp_groupedge)
  }
  )
  
  #--- Heterogeneity
  ## divide no. of group edges by no. of total edges
  l_heterogeneity <- list()
  for(i in 1:length(l_adjacency)){
    # calculate number of estimated edges, group + individual
    l_adjacency_ind <- ifelse(l_res[[i]]$path_counts > 0, 1, 0)
    l_heterogeneity[[i]] <- sum(l_adjacency[[i]])/sum(l_adjacency_ind)
  }
  
  #--- Output
  l_out <- tibble(
    n_ind = n_ind,
    ref_diff_g = l_ref_diff, 
    adjacency_g = l_adjacency,
    hetereogeneity_g = l_heterogeneity
  ) %>% 
    dplyr::mutate(ref_diff_g = as.list(ref_diff_g)) %>% 
    tidyr::unnest(n_ind)
  
  
}


multiverse.compare.subgroup <- function(l_res, 
                                        ref_model){
  
  #--- Reference model info
  n_ind <- length(ref_model$data)
  ref_sim_matrix <- ref_model$sim_matrix
  
  #--- Number of subgroups
  l_n_sub <- lapply(l_res, function(x){
    return(length(x$path_counts_sub))
  }
  )
  

  #--- Size of subgroups
  l_size_sub <- lapply(l_res, function(x){
    return(table(x$fit$sub_membership))
  })
  
  #--- Subgroup similarity matrix & modularity
  ## Use similarity matrix of each specification
  l_pert <- lapply(l_res, function(x){
    l_calc <- list()
    # can use VI and ARI from perturbR package here
    l_calc$vi <- vi.dist(ref_sim_matrix, x$sim_matrix)
    l_calc$ari <- arandi(ref_sim_matrix, x$sim_matrix)
    l_calc$modularity <- as.numeric(x$fit$modularity[1])   
    return(l_calc)
  })

  #--- Subgroup edges
  ## Loop over subgroups
  l_out <- tibble(
    n_sub_g = l_n_sub,
    size_sub_s = l_size_sub,
    l_pert = l_pert
  ) %>% 
    tidyr::unnest_wider(l_pert)
    
  return(l_out)
}



multiverse.compare.individual <- function(l_res,
                                          ref_model){
  
  #--- Reference model 
  n_ind <- length(ref_model$data)
  n_var <- ref_model$n_vars_total
    
  # indices for temporal and contemporaneous
  temp_ind <- 1:(n_var/2)
  cont_ind <- ((n_var/2)+1):n_var
    
  ## Estimates
  ref_path_est_mats <- ref_model$path_est_mats
  
  
  #--- Adjacency matrix
  # ignore autoregressive effects
  ref_adj_mats <- lapply(ref_path_est_mats, function(x){
    tmp <- ifelse(x != 0, 1, 0)
    diag(tmp[, temp_ind]) <- rep(0, n_var)
    return(tmp)
  })
  
  #--- Density
  ref_dens_temp <- lapply(ref_path_est_mats, function(x){
    abs_sum(x[, temp_ind])
  })
  
  ref_dens_cont <- lapply(ref_path_est_mats, function(x){
    abs_sum(x[, cont_ind])
  })
  
  #--- Fit indices
  fit_ind_names <- c("chisq", "df", "npar", "pvalue", "rmsea", "srmr",
                     "nnfi", "cfi", "bic", "aic", "logl")
  ref_fit_ind <- ref_model$fit[names(ref_model$fit) %in% fit_ind_names]
  
  #--- Centrality
  ref_outstrength <- lapply(ref_path_est_mats, function(x){
    colSums(abs(x))
  })

  ########################
  #--- Nonconverg. Checks
  ########################
  # TODO
  # can only add this after I have seen some actual nonconvergence
  
  
  ########################
  #--- Plausibility Checks
  ########################
  
  # Check diagonal of Psi (contemporaneous) matrix
  # Values should be >= 0 & <= 1
  # inspired by https://github.com/aweigard/GIMME_AR_simulations/blob/master/analyze_recovery_Balanced.R
  # 1 means implausible value somewhere in psi matrix
  l_implausible <- lapply(l_res, function(x){
    lapply(x$path_est_mats, function(y){
      sum(ifelse(any(diag(y[,temp_ind]) < 0 | diag(y[,temp_ind]) > 1) , 1, 0))
    })
  })
  sum_implausible <- sapply(l_plausible, sum)
  
  
  ########################
  #--- Compare Edges
  ########################
  #--- Nondirectional recovery
  # Only for contemporaneous effects
  ref_nondir_adj_mats <- lapply(ref_path_est_mats, function(x){
    tmp <- ifelse(x[,cont_ind] != 0, 1, 0)
    tmp <- nondirect_adjacency(tmp)
    return(tmp)
  })
  
  l_diff_nondir_adj <- lapply(l_res, function(x){
    tmp_nondir_adj_mats <- lapply(x$path_est_mats, function(y){
      tmp <- ifelse(x[,cont_ind] != 0, 1, 0)
      tmp <- nondirect_adjacency(tmp)
      return(tmp)
    })
    l_nondir_adj <- Map('-', ref_nondir_adj_mats, tmp_nondir_adj_mats)
    lapply(l_nondir_adj, function(y){as.matrix(y)})
    
  })
  
  #--- Directional recovery
  l_diff_adj <- lapply(l_res, function(x){
    ## adjacency matrices
    # ignore AR coefs
    tmp_adj_mats <- lapply(x$path_est_mats, function(y){
      tmp <- ifelse(x != 0, 1, 0)
      tmp <- diag(tmp[, temp_ind]) <- rep(0, n_var/2)
      return(tmp)
    })
    l_adj <- Map('-', ref_adj_mats, tmp_adj_mats)
    lapply(l_adj, function(y){as.matrix(y)})
  })
  #--- Difference/bias path estimates
  l_diff_ests <- lapply(l_res, function(x){
    l_est <- Map('-', ref_path_est_mats, x$path_est_mats)
    lapply(l_est, function(y){as.matrix(y)})
  })

  
  #--- Density
  l_diff_dens_temp <- lapply(l_res, function(x){
    tmp_denstemp <- lapply(x$path_est_mats, function(y){
      abs_sum(x[, temp_ind])
    })
    Map('-', ref_dens_temp, tmp_denstemp)
  })
  
  l_diff_dens_cont <- lapply(l_res, function(x){
    tmp_denscont <- lapply(x$path_est_mats, function(y){
      abs_sum(x[, cont_ind])
    })
    Map('-', ref_dens_cont, tmp_denstemp)
  })
  
  #--- Compare Fits
  l_diff_fit <- lapply(l_res, function(x){
    return(ref_fit_ind - x$fit[names(x$fit) %in% fit_ind_names])
  })
  
  
  #--- Compare centrality
  # compute outstrength
  l_cent <- lapply(l_res, function(x){
    tmp_outstrength <- lapply(x$path_est_mats, function(y){
      colSums(abs(y))
    })
  })
  
  # compute difference
  l_diff_cent <- lapply(l_res, function(x){
    tmp_outstrength <- lapply(x$path_est_mats, function(y){
      colSums(abs(y))
    })
    Map('-', ref_outstrength, tmp_outstrength)
  })
  
  ## binary: most central node the same?
  # split by temporal and contemporaneous
  central_node_identical <- list()
  
  for(i in seq_along(l_cent)){
    max_temp_ref <- which.max(ref_outstrength[[i]][temp_ind])
    max_temp_mv <- which.max(l_cent[[i]][temp_ind])
    max_cont_ref <- which.max(ref_outstrength[[i]][cont_ind])
    max_cont_mv <- which.max(l_cent[[i]][cont_ind])
    central_node_identical[[i]] <- list()
    central_node_identical[[i]]$temp_identical <- max_temp_ref == max_temp_mv
    central_node_identical[[i]]$cont_identical <- max_cont_ref == max_cont_mv
  }
  
  
  ########################
  #--- Aggregate
  ########################
  # First aggregate by taking mean for each edge across individuals
  # then summarize this again
  # split by averaging over all nonzero differences, or all differences
  
  #--- Nondirected adjacency
  mean_diff_nondir_adj <- lapply(l_diff_nondir_adj, function(x){
    l_tmp <- list()
    l_tmp$nondir_adj_sum_mat_i <- apply(simplify2array(x), 1:2, abs_sum)
    l_tmp$nondir_adj_sum_sum_i <- sum(l_tmp$nondir_adj_sum_mat)
    l_tmp$nondir_adj_sum_mean_i <- mean(l_tmp$nondir_adj_sum_mat)
    return(l_tmp)
  })
  
  #--- Adjacency matrix
  mean_diff_adj <- lapply(l_diff_adj, function(x){
    l_tmp <- list()
    l_tmp$adj_sum_mat_i <- apply(simplify2array(x), 1:2, abs_sum)
    l_tmp$adj_sum_sum_i <- sum(l_tmp$adj_sum_mat)
    l_tmp$adj_sum_mean_i <- mean(l_tmp$adj_sum_mat)
    return(l_tmp)
  })
  
  
  #--- Mean differences of edges
  mean_diff_ests <- lapply(l_diff_ests, function(x){
    l_tmp <- list()
    mean_mat <- apply(simplify2array(x), 1:2, mean)
    l_tmp$mean_nonzero_diff_edge_i <- mean(abs(mean_mat[mean_mat != 0]))
    l_tmp$med_nonzero_diff_edge_i <- stats::median(abs(mean_mat[mean_mat != 0]))
    l_tmp$mean_diff_edge_i <- mean(abs(mean_mat))
    l_tmp$med_diff_edge_i <- stats::median(abs(mean_mat))
    return(l_tmp)
  })
  
  #--- Densities
  mean_diff_dens_temp <- sapply(l_diff_dens_temp, mean)
  mean_diff_dens_cont <- sapply(l_diff_dens_cont, mean)
  
    
  #--- Fits 
  mean_diff_fits <- lapply(l_diff_fit, function(x){
    l_tmp <- list()
    l_tmp$mean_diff_fit_i <- apply(x, 2, abs_mean)
    l_tmp$med_diff_fit_i <- apply(x, 2, abs_med)
    return(l_tmp)
  })
  
  
  #--- Centrality
  # difference across all centrality values
  mean_diff_cent <- lapply(l_diff_cent, function(x){
    l_tmp <- list()
    diff_cent <- rowSums(simplify2array(x))
    l_tmp$mean_diff_cent_i <- mean(abs(diff_cent))
    l_tmp$med_diff_cent_i <- stats::median(abs(diff_cent))
    return(l_tmp)
  })
  
  # central node identical
  sum_temp_central_identical <- sum(sapply(central_node_identical, function(x) x$temp_identical))
  sum_cont_central_identical <- sum(sapply(central_node_identical, function(x) x$cont_identical))

  

  l_out <- tibble(
    l_implausible_i = l_implausible,
    sum_implausible_i = sum_implausible,
    l_diff_nondir_adj_i = l_diff_nondir_adj,
    l_diff_adj_i = l_diff_adj,
    l_diff_ests_i = l_diff_ests,
    l_diff_fit_i = l_diff_fit,
    l_diff_cent_i = l_diff_cent,
    central_node_identical_i = central_node_identical,
    mean_diff_nondir_adj_i = mean_diff_nondir_adj,
    mean_diff_adj_i = mean_diff_adj,
    mean_diff_ests_i = mean_diff_ests,
    mean_diff_cent_i = mean_diff_cent,
    mean_diff_dens_temp_i = mean_diff_dens_temp,
    mean_diff_dens_cont_i = mean_diff_dens_cont,
    sum_temp_central_identical_i = sum_temp_central_identical,
    sum_cont_central_identical_i = sum_cont_central_identical,
    mean_diff_fit_i = mean_diff_fits
  ) %>% 
    tidyr::unnest_wider(c(mean_diff_nondir_adj_i,
                          mean_diff_adj_i,
                          mean_diff_ests_i, 
                          mean_diff_cent_i,
                          mean_diff_fit_i))
  return(l_out)
}








# From perturbr internal functions ----------------------------------------
# https://github.com/cran/perturbR/blob/master/R/vi.dist.R
vi.dist <-
  function(cl1,cl2,parts=FALSE, base=2){ 
    if(length(cl1) != length(cl2)) stop("cl1 and cl2 must have same length")
    
    # entropy 
    ent <- function(cl){
      n <- length(cl)
      p <- table(cl)/n
      -sum(p*log(p, base=base))
    }
    # mutual information
    mi <- function(cl1,cl2){
      p12 <- table(cl1,cl2)/length(cl1)
      p1p2 <- outer(table(cl1)/length(cl1),table(cl2)/length(cl2))
      sum(p12[p12>0]*log(p12[p12>0]/p1p2[p12>0], base=base))
    }
    
    if(!parts) return(ent(cl1) + ent(cl2) -2*mi(cl1,cl2))
    ent1 <- ent(cl1)
    ent2 <- ent(cl2)
    mi12 <- mi(cl1,cl2)
    c("vi"=ent1+ent2-2*mi12, "H(1|2)" =ent1-mi12, "H(2|1)"=ent2 -mi12)
  }

# https://github.com/cran/perturbR/blob/master/R/arandi.R
arandi <-
  function(cl1,cl2, adjust=TRUE){
    if(length(cl1)!=length(cl2)) stop("cl1 and cl2 must have same length")
    tab.1 <- table(cl1)
    tab.2 <- table(cl2)
    tab.12 <- table(cl1,cl2)
    if(adjust){
      correc <- sum(choose(tab.1,2))*sum(choose(tab.2,2))/choose(length(cl2),2)
      return((sum(choose(tab.12,2))-correc)/(0.5*sum(choose(tab.1,2))+0.5*sum(choose(tab.2,2))-correc) )}
    else{ 
      1+(sum(tab.12^2)-0.5*sum(tab.1^2)-0.5*sum(tab.2^2))/choose(length(cl2),2)
    }
  }



# Matrix list summary -----------------------------------------------------
matrix_summary <- function(l_matrix) {
  # Calculate mean, median, sd, min, and max for each matrix element across the list
  arr_list <- simplify2array(l_matrix)
  
  a_mean <- apply(arr_list, 1:2, mean)
  a_median <- apply(arr_list, 1:2, stats::median)
  a_min <- apply(arr_list, 1:2, min)
  a_max <- apply(arr_list, 1:2, max)
  a_sd <- apply(arr_list, 1:2, sd)
  
  # Combine the results into a data frame
  summary_df <- as.data.frame(summary_list)
  
  return(summary_df)
}


# Small helpers -----------------------------------------------------------

abs_mean <- function(x){
  mean(abs(x))
}
abs_med <- function(x){
  stats::median(abs(x))
}

abs_sum <- function(x){
  sum(abs(x))
}


nondirect_adjacency <- function(adj_mat) {
  # Number of ariables
  n_adj_vars <- nrow(adj_mat)
  
  # Initialize symmetrical matrix with 0s
  sym_matrix <- matrix(0, nrow = n_adj_vars, ncol = n_adj_vars)
  
  # Iterate through each cell of the original matrix
  for (i in 1:n_adj_vars) {
    for (j in 1:n_adj_vars) {
      # If there is any effect (1) in either direction, update the symmetrical matrix
      if (n_adj_vars[i, j] == 1 || n_adj_vars[j, i] == 1) {
        sym_matrix[i, j] <- 1
        sym_matrix[j, i] <- 1
      }
    }
  }

  return(sym_matrix)
}




