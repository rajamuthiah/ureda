#' Compute the coefficients of the Betas (Internal function).
#' 
#' @param r the number of coefficients to be calculated (a whole number).
#' @return a vector of coefficients of length r.
pstar <- function(r){
  
  k <- 0:r
  
  (-1)^(r - k) * choose(r, k) * choose(r + k, k)
  
}

#' Compute a single term of the Betas (Internal function).
#' 
#' @param x sample (numeric vector).
#' @param k the k-th Beta term.
#' @return a numeric vector of length 1 (the k-th Beta term).
bterm <- function(x, k){
  
  n <- length(x) 
  
  orderd.x <- sort(x, decreasing = FALSE)
  
  j <- (k + 1):n
  
  (1 / n) * (1 / choose(n - 1, k)) * sum(choose(j - 1, k) * orderd.x[j])
  
}

#' Compute all the Beta (Internal function).
#' 
#' Vectorize internal.b
#' @param x numeric vector (ie. sample data).
#' @param r the of Betas terms.
#' @return a numeric vector of length r (Betas).
betas <- function(x, r){
  
  sapply(0:r, bterm, x = x)
  
}  

#' Compute the r-th l-statistic from a given sample (Internal function).
#' 
#' @param smpl sample (numeric vector).
#' @param r the r-th sample l-statistic (a single whole number).
#' @return a numeric vector of length 1 (the r-th l-statistic).
singlelmom <- function(smpl, r){
  
  sum(betas(smpl, r - 1) * pstar(r - 1))
  
}

#' Sample L-moments ratios.
#' 
#' This function computes the l-moments ratios for a given
#' dataset (vector / matrix / data.frame).
#' @param x sample (vector / matrix / data.frame).
#' @param ratios Compute l-moments ratios instead of the raw l-moments? (logical, default = TRUE).
#' @param n the number of l-moments (ratios) to compute. (a whole number, default = 4).
#' @return a numeric vector of length n / matrix with 4 columns containing 
#' sample l-moments / l-moments ratios.
#' @export
shape <- function(x, ratios = TRUE, n = 4){
  
  tol <- .Machine$double.eps^0.5
  is.natural <- n > tol & abs(n - round(n)) < tol
  
  if (!is.natural) {
    
    stop("n must be a natural number (ie. n = 1, 2, 3, ....).")
    
  }
  
  if (is.vector(x)) {
    
    if (n > length(x)) {
      
      stop("n must be less than or equal to sample size (length of x).")
      
    }
    
    res <- sapply(1:n, singlelmom, smpl = x)
    names(res) <- paste0("l", 1:n)
    
    if (ratios) {
      
      res[3:n] <- res[3:n] / res[2]  
      names(res) <- c(paste0("l", 1:2), paste0("t", 3:n))
    }
    
  }else{
    
    if (n > ncol(x)) {
      
      stop("n must be less than or equal to sample size (ncol of x).")
      
    }
    
        
    res <- t(apply(x, 1, function(y) sapply(1:n, singlelmom, smpl = y)))
    colnames(res) <- paste0("l", 1:n)

    if (ratios) {
      
      res[, 3:n] <-  res[, 3:n] / res[, 2]
      colnames(res) <- c(paste0("l", 1:2), paste0("t", 3:n))
      
    }
    
  }
  
  res
  
}
