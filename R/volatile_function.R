#' Find volatile genes.
#' 
#' @param t3 a vector containing L-skew estimates
#' for each gene.
#' @param t4 a vector containing L-kurt estimates
#' for each gene.
#' @param plot.fit indicating whether to show intermediary 
#' plots.
#' @return pvalues for each gene.
#' @export
volatile <- function(t3, t4, plot.fit = FALSE){  
  
  # Used in volatile to compute discord measure.
  discord <- function(x, mu, A.inv){
    u <- x - mu    
    t(u) %*% A.inv %*% u
  }
  
  # Fit lowess.
  span <- 0.5
  lowess.fit <- lowess(t3, t4, f=span)
  approx.lowess.fit <- approxfun(lowess.fit, rule=2)
  
  # Adjust t4.
  adj.t4 <- t4 - approx.lowess.fit(t3)
  
  # Compute discordancy pvalues
  shape <- cbind(t3, adj.t4)
  A <- cov(shape)
  mu <- colMeans(shape)
  A.inv <- solve(A)
  D <- apply(shape, 1, discord, mu=mu, A.inv=A.inv)
  pvals <- 1 - pchisq(D, df=2)
  
  if(plot.fit){
    
    oldpar <- par(mfrow=c(2, 2), mar=c(4, 4, 2, 1))
    
    Lab.palette <- colorRampPalette(c("gray75", "gray60", "gray35", "gray20"), space = "Lab")
    data.pts.col <- densCols(t3, t4, colramp=Lab.palette)
    cex <- 0.35
    col <- data.pts.col
    
    xx <- pvals <= 0.0001
    
    # Adjustment plot
    plot(t3, t4,  
         cex = cex, 
         pch = 16,
         xlab = expression(tau[3]),
         ylab = expression(tau[4]),
         main = expression(Kurt-Skew~Trend),
         col = col)
    lines(lowess.fit, col = "red")
    points(t3[xx], t4[xx], pch=4, cex=0.5)
    
    # Adjusted plot
    plot(t3, adj.t4,  
         cex = cex, 
         pch = 16,
         xlab = expression(tau[3]),
         ylab = expression(adjusted~tau[4]),
         main = expression(Kurt-Skew~Adjusted),
         col = col) 
    abline(v=0, col="gray")
    points(t3[xx], adj.t4[xx], pch=4, cex=0.5)
    
    # Plot discordancy scores.
    hist(D, freq=FALSE, breaks=40, 
         main="Discordancy measure (D)")
    curve(dchisq(x, df=2), 0, max(D), col="red", add=TRUE)
    
    # Pvals vs t3
    xx <- !xx
    plot(t3[xx], -log10(pvals)[xx], pch=".",
         xlab = expression(tau[3]),
         main = "P-values", 
         ylim = c(0, 4),
         xlim = range(t3),
         ylab="-log10(pvals)")
    abline(h=-log10(c(0.0001, 0.001, 0.01)), col="red", lty=c(1, 2, 4))
    abline(v=0, col="gray")
    points(t3[!xx], rep(4, sum(!xx)), pch=4, cex=0.5)
    
    par(oldpar)
  }
  
  pvals
}
