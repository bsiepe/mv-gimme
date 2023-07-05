#' Counts number of excellent fit indices
#' @param indices A vector of fit indices from lavaan.
#' @return The number of fit indices that are excellent.
#' @keywords internal 
count.excellent <- function(indices,
                            rmsea.cut = .05,
                            srmr.cut = .05,
                            nnfi.cut = .95,
                            cfi.cut = .95){
  rmseaE    <- ifelse(indices[4] < rmsea.cut, 1, 0)
  srmrE     <- ifelse(indices[5] < srmr.cut, 1, 0)
  nnfiE     <- ifelse(indices[6] > nnfi.cut, 1, 0)
  cfiE      <- ifelse(indices[7] > cfi.cut, 1, 0)
  excellent <- sum(rmseaE, srmrE, nnfiE, cfiE, na.rm = TRUE)
  return(excellent)
}
