multiverse.gimme <- function(data, 
                             groupcutoffs = .75, 
                             subcutoffs = .51,
                             rmsea.cuts = .05,
                             srmr.cuts = .05,
                             nnfi.cuts = .95,
                             cfi.cuts = .95,
                             n.excellent = 2,
                             n.cores = 1,
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
    l_out <- lapply(mv_combs, function(x){
      mvgimme::gimme(data = data,
                     ...)
    })
    
    
  }
  
  # Parallel fitting
  if(n.cores >1){
    
    progressr::handlers(global = TRUE)
    progressr::handlers("txtprogressbar")
    
    # Setup parallelization
    future::plan(future::multisession, workers = n.cores)
    
    # progress bar
    progressr::progressor(along = combs)
    
    # Loop over inputs
    l_out <- future.apply::future_lapply(1:nrow(combs), 
                                function(i){
      
    res <- mvgimme::gimme(data = data,
                     groupcutoff = combs[i, "groupcutoffs"],
                     subcutoff = combs[i,"subcutoffs"],
                     rmsea.cut = combs[i,"rmsea.cuts"],
                     srmr.cut = combs[i,"srmr.cuts"],
                     nnfi.cut = combs[i,"nnfi.cuts"],
                     cfi.cut = combs[i,"cfi.cuts"],
                     n.excellent = combs[i,"n.excellent"],
                     plot = FALSE,
                     ...)
    Sys.sleep(0.001)
    p()
    return(res)
    
    })
    
    

  }
  
  
  
  
  
  # Store results
  return(l_out)
  
  
  
  
  
}
