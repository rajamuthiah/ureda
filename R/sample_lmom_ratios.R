#' Weight Coefficient for Sample L-moment. (Internal)
#' 
#' @param r The r-th sample l-moment.
#' @param j The j-th term.
#' @param n The sample size.
#' @return The j-th weight factor of the r-th sample l-moment.
weight.factor <- function (r, j, n) {
  
  upper.limit <- min(j - 1, r - 1)
  i <- 0:upper.limit
  
  factor1 <- (-1)^(r - 1 - i)
  factor2 <- choose(r - 1, i)
  factor3 <- choose(r - 1 + i, i)
  factor4 <- choose(j - 1, i)
  factor5 <- choose(n - 1, i)
  
  sum(factor1 * factor2 * factor3 * factor4 / factor5)
  
}

#' Weight Coefficients for Sample L-moment (Internal).
#' 
#' @param r The r-th sample l-moment.
#' @param n The sample size.
#' @return A vector of weight factors.
weightFactors <- function (r, n) {
  
  sapply(1:n, weight.factor, r = r, n = n)
  
}

#' Weight Coefficients for Sample L-moments.
#' 
#' @param nlmom The number of sample l-moments.
#' @param n The sample size.
#' @return A matrix of weight factors.
#' @export
weightMatrix <- function (nlmom, n) {
  
  if (n < 1) {
    stop("n must be greater than or equal to 1.")
  }
  
  if (nlmom < 1) {
    stop("nlmom must be greater than or equal to 1.")
  }
  
  if (nlmom > n) {
    stop("nlmom must be less than or equal to n.")
  }
  
  res <- sapply(1:nlmom, weightFactors, n = n)
  colnames(res) <- paste0("W", 1:nlmom)
  
  res
  
}

#' Sample L-moments ratios.
#' 
#' This function computes the sample l-moments ratios for a given
#' dataset (vector / matrix / data.frame).
#' @param x Data (vector / matrix / data.frame).
#' @param nlmom The number of l-moments (ratios) to compute (Default = 4).
#' @param ratios Compute sample l-moments ratios ? (logical, Default = TRUE).
#' @return A list holding sample l-ratios, sample l-moments, and weights.
#' @export
shape <- function(x,  nlmom = 4, ratios = TRUE){
  
  if (!is.matrix(x)) {
    x <- t(x) 
  }
  
  n <- ncol(x)
  
  weight.mat <- weightMatrix(nlmom, n)
  x <- t(apply(x, 1, sort))
  lmoms <- (x %*% weight.mat) / n
  colnames(lmoms) <- paste0("l", 1:nlmom)
  
  if (ratios & nlmom > 2) {
    
    lrats <- lmoms
    colnames(lrats) <- c(paste0("l", 1:2), paste0("t", 3:nlmom))
    lrats[, 3:nlmom] <- lrats[, 3:nlmom] / lrats[, 2]
    
  }else{
    
    lrats = NULL
    
  }
  
  list(lrats = lrats, lmoms = lmoms, W = weight.mat)
  
}
