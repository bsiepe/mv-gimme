#' Multiverse GIMME
#'
#' This function is a wrapper around the `gimme` function. It allows fitting a multiverse of models with different specifications in parallel. 
#'
#' @param data The data for fitting the models.
#' @param groupcutoffs The cutoff values for group-level fit indices. Should be a numeric value or a vector of numeric values.
#' @param subcutoffs The cutoff values for subgroup-level fit indices. Should be a numeric value or a vector of numeric values.
#' @param rmsea.cuts The cutoff values for RMSEA (Root Mean Square Error of Approximation) fit indices.
#' @param srmr.cuts The cutoff values for SRMR (Standardized Root Mean Residual) fit indices.
#' @param nnfi.cuts The cutoff values for NNFI (Non-Normed Fit Index) fit indices.
#' @param cfi.cuts The cutoff values for CFI (Comparative Fit Index) fit indices.
#' @param n.excellent The number of excellent fitting indices to require for model selection.
#' @param n.cores The number of cores to use for parallel fitting. Set to 1 for nonparallel fitting.
#' @param prune_output Logical; if TRUE, the model data will be removed from the output list for reduced memory usage.
#' @param prune_output Logical; if TRUE, every iteration will be saved as .RDS
#' @param ... Additional arguments to be passed to \code{\link[mvgimme::gimme]{mvgimme::gimme}} function.
#' @return A list containing the results of the multiverse analysis for different parameter combinations. Each element of the list corresponds to one specification, and the conditions used for fitting are attached as a data frame to each result.
#' @export
multiverse.gimme <- function(data, 
                             groupcutoffs = .75, 
                             subcutoffs = .51,
                             rmsea.cuts = .05,
                             srmr.cuts = .05,
                             nnfi.cuts = .95,
                             cfi.cuts = .95,
                             n.excellent = 2,
                             n.cores = 1,
                             prune_output = TRUE,
                             save_output = FALSE, 
                             ...){
  
  # Input checks
  
  # Create grid of combinations
  combs <- expand.grid(
    groupcutoffs = groupcutoffs,
    subcutoffs = subcutoffs,
    rmsea.cuts = rmsea.cuts,
    srmr.cuts = srmr.cuts,
    nnfi.cuts = nnfi.cuts,
    cfi.cuts = cfi.cuts,
    n.excellent = n.excellent
  )
  
  # Setup output
  # l_out <- list()
  
  # Nonparallel fitting
  if(n.cores == 1){
    l_out <- lapply(1:nrow(combs), function(i){
      tryCatch({mvgimme::gimme(data = data,
                     groupcutoff = combs[i, "groupcutoffs"],
                     subcutoff = combs[i,"subcutoffs"],
                     rmsea.cut = combs[i,"rmsea.cuts"],
                     srmr.cut = combs[i,"srmr.cuts"],
                     nnfi.cut = combs[i,"nnfi.cuts"],
                     cfi.cut = combs[i,"cfi.cuts"],
                     n.excellent = combs[i,"n.excellent"],
                     ...)}, error = function(e) return(NA))
    })
    
    
  }
  
  # Parallel fitting
  if(n.cores >1){
    
    progressr::handlers(global = TRUE)
    progressr::handlers("txtprogressbar")
    
    # Setup parallelization
    future::plan(future::multisession, workers = n.cores)
    
    # progress bar
    p <- progressr::progressor(along = combs)
    
    # Loop over inputs
    l_out <- future.apply::future_lapply(1:nrow(combs),
                                         future.seed = TRUE,
                                function(i){
      
    res <- tryCatch({mvgimme::gimme(data = data,
                     groupcutoff = combs[i, "groupcutoffs"],
                     subcutoff = combs[i,"subcutoffs"],
                     rmsea.cut = combs[i,"rmsea.cuts"],
                     srmr.cut = combs[i,"srmr.cuts"],
                     nnfi.cut = combs[i,"nnfi.cuts"],
                     cfi.cut = combs[i,"cfi.cuts"],
                     n.excellent = combs[i,"n.excellent"],
                     ...)}, error = function(e) return(NA))
    
    # For progress bar
    Sys.sleep(0.001)
    p()
    
    
    
    return(res)
    
    })
    
    # Stop multisession explicitly
    future::plan(future::sequential)

  }
  
  # Attach conditions to output
  for(i in 1:nrow(combs)){
    l_out[[i]]$conds <- combs[i,]
  }
  
  # Prune output 
  if(isTRUE(prune_output)){
    l_out <- lapply(l_out, function(x){
      x$data <- NULL
      x$psi <- NULL
      x$psi_unst <- NULL
      x$vcov <- NULL
      x$vcovfull <- NULL
      x$path_se_est <- NULL
      x$plots <- NULL
      return(x)
    
    })
  }
  
  
  
  
  # Store results
  return(l_out)
  
  
  
  
  
}
