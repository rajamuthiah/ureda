#' A traditional moments based skewness estimator.
#'
#' @param x sample (numeric vector).  
#' @return the (traditional) moments skewness measure.
#' @export 
skew <- function (x) {
  
  n  <-  length(x)
  
  m3 <- mean((x - mean(x))^3)
  
  m3 <- (n^2 / ((n - 1) * (n - 2))) * m3
  
  m3 / var(x)^1.5
  
}

#' A traditional moments based kurtosis estimator.
#'
#' @param x sample (numeric vector).
#' @return the (traditional) moments skewness measure.
#' @export 
kurt <- function (x) {
  
  n <- length(x)
  
  m2 <- mean((x - mean(x))^2)
  m4 <- mean((x - mean(x))^4)
  
  m4 <- (n^2 / ((n - 2) * (n - 3))) * (((n + 1) / (n - 1)) * m4 - 3 * m2^2)  
  
  m4 / var(x)^2 + 3
  
}