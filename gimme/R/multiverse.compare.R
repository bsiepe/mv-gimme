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
  
  # reference model info
  n_ind <- length(ref_model$data)
  
  ## Find group effects adjacency matrix
  ref_groupedge <- ifelse(ref_model$path_counts == n_ind, 1, 0)
  
  # Find group effects adjacency matrix
  l_ref_diff <- lapply(l_res, function(x){
    tmp_groupedge <- ifelse(x$path_counts == n_ind, 1, 0)
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

  
  #--- Output
  l_out <- tibble(
    n_ind = n_ind,
    ref_diff_g = l_ref_diff, 
    adjacency_g = l_adjacency
  ) %>% 
    dplyr::mutate(ref_diff_g = as.list(ref_diff_g)) %>% 
    tidyr::unnest(n_ind)
  
  
}


multiverse.compare.subgroup <- function(l_res, 
                                        ref_model){
  
  # reference model info
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
    vi_s = l_pert$vi,
    ari_s = l_pert$ari,
    modularity_s = l_pert$modularity
  ) 
    
  return(l_out)
}



multiverse.compare.individual <- function(l_res,
                                          ref_model){
  
  #--- Reference model 
  n_ind <- length(ref_model$data)
  
  ## Estimates
  ref_path_est_mats <- ref_model$path_est_mats
  
  ## Adjacency matrix
  ref_adj_mats <- lapply(ref_path_est_mats, function(x){
    ifelse(x != 0, 1, 0)
  })
  
  ## Fit indices
  fit_ind_names <- c("chisq", "df", "npar", "pvalue", "rmsea", "srmr",
                     "nnfi", "cfi", "bic", "aic", "logl")
  ref_fit_ind <- ref_model$fit[names(ref_model$fit) %in% fit_ind_names]
  
  ## Centrality
  ref_outstrength <- lapply(ref_path_est_mats, function(x){
    colSums(abs(x))
  })

  
  
  #--- Compare Edges
  l_diff_adj <- lapply(l_res, function(x){
    ## adjacency matrices
    tmp_adj_mats <- lapply(x$path_est_mats, function(y){
      ifelse(y != 0, 1, 0)
    })
    l_adj <- Map('-', ref_adj_mats, tmp_adj_mats)
    lapply(l_adj, function(y){as.matrix(y)})
  })
  ## path estimates
  l_diff_ests <- lapply(l_res, function(x){
    l_est<- Map('-', ref_path_est_mats, x$path_est_mats)
    lapply(l_est, function(y){as.matrix(y)})
  })

  #--- Compare Fits
  l_diff_fit <- lapply(l_res, function(x){
    return(ref_fit_ind - x$fit[names(x$fit) %in% fit_ind_names])
  })
  
  #--- Compare centrality
  l_diff_cent <- lapply(l_res, function(x){
    tmp_outstrength <- lapply(x$path_est_mats, function(y){
      colSums(abs(y))
    })
    Map('-', ref_outstrength, tmp_outstrength)
  })
  
  #--- Aggregate
  # First aggregate by taking mean for each edge across individuals
  # then summarize this again
  # split by averaging over all nonzero differences, or all differences
  
  # Adjacency matrix
  mean_diff_adj <- lapply(l_diff_ests, function(x){
    l_tmp <- list()
    l_tmp$adj_sum <- apply(simplify2array(x), 1:2, abs_sum)
    
    
  })
  
  
  # Mean differences of edges
  mean_diff_ests <- lapply(l_diff_ests, function(x){
    l_tmp <- list()
    mean_mat <- apply(simplify2array(x), 1:2, mean)
    l_tmp$mean_nonzero_diff_edge_i <- mean(abs(mean_mat[mean_mat != 0]))
    l_tmp$med_nonzero_diff_edge_i <- stats::median(abs(mean_mat[mean_mat != 0]))
    l_tmp$mean_diff_edge_i <- mean(abs(mean_mat))
    l_tmp$med_diff_edge_i <- stats::median(abs(mean_mat))
    return(l_tmp)
  })
    
  # Fits 
  mean_diff_fits <- lapply(l_diff_fit, function(x){
    l_tmp <- list()
    l_tmp$mean_diff_fit_i <- apply(x, 2, abs_mean)
    l_tmp$med_diff_fit_i <- apply(x, 2, abs_med)
    return(l_tmp)
  })
  
  
  # Centrality
  mean_diff_cent <- lapply(l_diff_cent, function(x){
    l_tmp <- list()
    diff_cent <- rowSums(simplify2array(x))
    l_tmp$mean_diff_cent_i <- mean(abs(diff_cent))
    l_tmp$med_diff_cent_i <- stats::median(abs(diff_cent))
    return(l_tmp)
  })
  
  
  
  l_out <- tibble(
    n_ind_test = n_ind,
    l_diff_adj_i = l_diff_adj,
    l_diff_ests_i = l_diff_ests,
    l_diff_fit_i = l_diff_fit,
    l_diff_cent_i = l_diff_cent,
    mean_diff_ests_i = mean_diff_ests,
    mean_diff_cent_i = mean_diff_cent,
    mean_diff_fit_i = mean_diff_fits
  ) %>% 
    tidyr::unnest_wider(c(mean_diff_ests_i, 
                          mean_diff_cent_i,
                          mean_diff_fit_i))
  return(l_out)
}







# From perturbr internal functions ----------------------------------------
# https://github.com/cran/perturbR/blob/master/R/vi.dist.R
vi.dist <-
  function(cl1,cl2,parts=FALSE, base=2){ # wenn parts=TRUE, werden die Komponenten der VI ebenfalls berechnet
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
