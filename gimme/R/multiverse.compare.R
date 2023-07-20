multiverse.compare <- function(l_out,
                               ref_model
){
  
  
  #---- Compare Group-Level
  
  
  #---- Compare Subgroup-Level
  
  
  #---- Compare Individual Solutions
  
  
  #--- Output
  
}


multiverse.compare.group <- function(l_out,
                                     ref_model){
  
  # reference model info
  n_ind <- length(ref_model$data)
  
  ## Find group effects adjacency matrix
  ref_groupedge <- ifelse(ref_model$path_counts == n_ind, 1, 0)
  
  # Find group effects adjacency matrix
  l_ref_diff <- lapply(l_out, function(x){
    tmp_groupedge <- ifelse(x$path_counts == n_ind, 1, 0)
    diff_groupedge <- ref_groupedge - tmp_groupedge
  }
  )
 
}


multiverse.compare.subgroup <- function(l_out, 
                                        ref_model){
  
  # reference model info
  n_ind <- length(ref_model$data)
  ref_sim_matrix <- ref_model$sim_matrix
  
  # Number of subgroups
  
  # Size of subgroups
  
  # Subgroup membership
  # But which count matrix do we need here?
  
    # can use VI and ARI from perturbR package here
    vi.dist(ref_sim_matrix, x)
    arandi(ref_sim_matrix, x)
    
  # Subgroup edges
  
    
}



multiverse.compare.individual <- function(l_out){
  
  # reference model info
  n_ind <- length(ref_model$data)
  
  ## estimates
  ref_path_est_mats <- ref_model$path_est_mats
  
  ## adjacency matrix
  ref_adj_mats <- lapply(ref_path_est_mats, function(x){
    ifelse(x != 0, 1, 0)
  })
  
  
  # Loop over models and individuals
  l_ind <- list()
  for(i in 1:n_ind){
    l_ind[[i]] <- list()
    
  }
  
  # Aggregate over individuals
  
  
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


