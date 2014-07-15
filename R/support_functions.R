#' Compute log counts per million (cpm).
#' 
#' @param counts expression count matrix.
#' @param lib.size library sizes 
#' (Defaul lib.size = NULL, will use total counts per sample
#' if not specified).
#' @return A list containing the log-cpm maxtrix and library size.
#' @export
log2CPM <- function (counts, lib.size = NULL) {
  
  if (is.null(lib.size)){
    lib.size <- colSums(counts)
  } 
  
  y <- t(log2(t(counts + 0.5) / (lib.size + 1) * 1e+06))
  
  list(y = y, lib.size = lib.size)
}

#' Filter low count genes.
#' 
#' @param counts expression count matrix.
#' @param lib.size library sizes 
#' (Defaul lib.size = NULL, will use total counts per sample
#' if not specified).
#' @param thresh expression threshold (Default thresh = 1 cpm).
#' @param minSamples minimum number of samples required to exceed
#' the threshold (Default minSamples = 2). 
#' @return The filtered count matrix.
#' @export
filterCounts <- function(counts, lib.size = NULL, thresh = 1, minSamples = 2) {
  
  cpms <- 2^log2CPM(counts, lib.size = lib.size)$y
  keep <- rowSums(cpms > thresh) >= minSamples
  counts <- counts[keep, ]
  counts
}

#' Compute residuals.
#' 
#' @param exprs expression matrix 
#' (counts / microarray intensities / methylations measures etc.)
#' @param cond input variable, a vector / factor indicating groups or treatment.
#' @return Residual expression levels.
#' @export
compute.residuals <- function(exprs, cond) {
  
  design <- model.matrix(~0 + as.factor(cond))
  beta.hat <- solve(t(design) %*% design) %*% t(design) %*% t(as.matrix(exprs))
  res <- exprs - t(design %*% beta.hat)
  as.data.frame(res)
} 
