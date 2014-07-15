#' The Symmetry Outlier Plot.
#' 
#' @param t3 a vector containing L-skew estimates
#' for each gene.
#' @param t4 a vector containing L-kurt estimates
#' for each gene.
#' @param data.name charcater vector indicating name
#' of dataset.
#' @return The symmetry outlier plot.
#' @export
plotSO <- function(t3, t4,
                   data.name = "Data name (RNA-seq / Microarray)"){

  
  data.pts.cex = 0.3,
  data.pts.pch = 16,
  data.med.lty = 1,
  data.med.col="red", 
  extreme.border.col = "black",
  g = seq(-5, 5, 0.05),
  h = seq(0.0, 0.5, 0.25),
  gh.col="black",
  normal.lty = 1,
  data.pts.col = NULL,
  grid = FALSE,
  labcex = 0.6
  
  
  if(is.null(data.pts.col)){
    Lab.palette = colorRampPalette(c("gray80", "gray65", "gray30", "gray20"), space = "Lab")
    data.pts.col = densCols(t3, t4, colramp=Lab.palette)
  }
  
  
  old.par <- par(mar=c(3.5, 2.5, 1.5, 0.5))
  
  plot(c(-1, 1.18), c(-0.8, 1),
       pch="",
       xaxt="none",
       yaxt="none",
       xlab="",
       ylab="",
       main="Symmetry-Outlier Plot (SO-Plot)",
       cex.main=0.8)
  
  points(t3, t4, 
         col=data.pts.col, 
         cex=data.pts.cex,
         pch=data.pts.pch)
  
  if(grid){
    abline(v=seq(-1, 1, .05), h=seq(-0.8, 1, 0.1), 
           lty=3, col="gray90")
  }
  
  abline(v=c(-1, 1), col="gray60") 
  
  points(seq(-1, 1, 0.01), 0.25*(5*seq(-1, 1, 0.01)^2 - 1),
         type="l", col=extreme.border.col)
  points(c(-1, 1), c(1, 1), type="l", col=extreme.border.col)
  
  boxplot(t4, add=T, horizontal=F, at=1.14, yaxt="none",
          boxwex=0.3, 
          cex=0.8)
  
  boxplot(t3, add=T, horizontal=T, at=-0.75, xaxt="none",
          boxwex=0.3, 
          cex=0.8)
  
  abline(v=c(-.35, -.2, -.05, .05, .2, .35), 
         col="gray40", 
         lty=c(1, 4, 2, 2, 4, 1),
         lwd=0.5)
  
  axis(1, seq(-1, 1, .05), F, tck=0.01)
  tmp1 = c(-1, -.75, -.5, -.35, -.2, -.05, .05, .2, .35, .5, .75, 1) 
  mtext(tmp1, at = tmp1,
        side=1, line=0.2, cex=0.6, las=3)
  mtext(expression(L-skew~(tau[3])), 
        1, line=1.5, cex=labcex)
  mtext(data.name, 1, line=2.2, cex=labcex)
  
  axis(2, seq(-0.8, 1, 0.05), F, tck=0.01)
  tmp2 = seq(-0.8, 1, 0.2)
  mtext(tmp2, at=tmp2, side=2, line=0.2, cex=0.6, las=2)
  mtext(expression(L-kurt~(tau[4])), 2, line=1.5, cex=labcex)
  
  abline(v=median(t3), col=data.med.col, lty=data.med.lty)
  abline(h=median(t4), col=data.med.col, lty=data.med.lty)
  points(median(t3), median(t4), col=data.med.col)
  
  for(k in h){
    tau4 = sapply(g, std.gh.lmom, h=k, lmom=4) / sapply(g, std.gh.lmom, h=k, lmom=2)
    tau3 = sapply(g, std.gh.lmom, h=k, lmom=3) / sapply(g, std.gh.lmom, h=k, lmom=2)
    points(tau3, tau4, type="l", lty=ifelse(k <= .25, 1, 2), col=gh.col)
  }
  
  for(k in seq(-1, 1, 0.25)){
    x = sapply(h, std.gh.lmom, g=k, lmom=3) / sapply(h, std.gh.lmom, g=k, lmom=2)
    y = sapply(h, std.gh.lmom, g=k, lmom=4) / sapply(h, std.gh.lmom, g=k, lmom=2)
    points(x, y, cex=0.5, col=gh.col, pch=16)
  }
  
  abline(h=0.1226, v=0, col="gray20")
  abline(h=0, v=0, col="gray20", lty=2)
  points(0, 0.1226, bg="white", pch=21)
  
  text(c(-.35, -.2, -.05) + .02, rep(.75, 3),
       rev(c("minor", "moderate", "large")), 
       srt=90, cex=0.6, pos=1)
  text(-.6, .75, "extreme/\nvolatile", cex=0.6)
  
  t3.summary <- round(quantile(t3, probs=c(0.25, 0.5, 0.75)), 2)
  t3.summary <- paste0(data.name, 
                       " L-skew: (25%, 50%, 75%) = (",
                       t3.summary[1], ", ",
                       t3.summary[2], ", ",
                       t3.summary[3], ")")
  print(t3.summary)
  
  par(old.par)
}

# Data prep: (filter and normalize) and compute residuals
data.prep <- function(path, data.name, residuals=FALSE){
  
  # Prep data for computation of L-moments ratios.
  print("Loading all datasets. (shapeData.rda)")
  
  load(paste0(path, "data/shapeData.rda"))
  
  if(data.name %in% c("bottomly", "maqc", "pickrell")){
    
    # Load data.
    print(paste0("Grabbing ", data.name, " dataset from shapeData."))
    
    exprs <- data[[data.name]]$exprs
    cond <- data[[data.name]]$cond
    
    # Filter low count genes.
    minSamples <- min(table(cond))
    print(paste0("Minimum samples: ", minSamples))
    exprs <- filterCounts(exprs, thresh=1, minSamples=minSamples) 
    exprs <- exprs + 1
    
    # Normalize using DESeq's method.
    log.exprs <- log2(DESeq.scale(exprs))
    
    print("Data has been filtered and normalized. (Currently on log-scale)")
    
    # Compute residuals.
    if(residuals){
      
      resids.log.exprs <- compute.residuals(log.exprs, cond)
      print("Residuals have been computed. (Currently on log-scale)")
    }
    
  }else{
    
    if(data.name == "hammoud"){
      
      # Load data. (hammoud:RPKM)
      
      print(paste0("Grabbing ", data.name, " dataset from shapeData."))
      
      exprs <- data[[data.name]]$exprs
      cond <- data[[data.name]]$cond
      
      # Filter low RPKM genes.
      minSamples <- min(table(cond))
      print(paste0("Minimum samples: ", minSamples))
      exprs <- exprs[rowSums(exprs > 1) > minSamples,]
      exprs <- exprs + 1
      
      log.exprs <- log2(exprs)
      
      print("Data has been filtered and normalized (RPKM). (Currently on log-scale)")
      
      # Compute residuals.
      if(residuals){
        
        resids.log.exprs <- compute.residuals(log.exprs, cond)
        print("Residuals have been computed. (Currently on log-scale)")
      }
      
    }else{
      
      # Load data (geng:microarray)
      # Data has been normalized and background corrected (see limma guide).
      
      print(paste0("Grabbing ", data.name, " dataset from shapeData."))
      
      log.exprs <- data[[data.name]]$exprs
      cond <- data[[data.name]]$cond
      
      print("Data is already normalized (see limma guide). (Currently on log-scale)")
      
      # Compute residuals.
      if(residuals){
        
        resids.log.exprs <- compute.residuals(log.exprs, cond)
        print("Residuals have been computed. (Currently on log-scale)")
      } 
    } 
  }
  
  
  if(residuals){
    
    res <- list(log.exprs=log.exprs, 
                resids.log.exprs=resids.log.exprs, 
                cond=cond)
    
  }else{
    
    res <- list(log.exprs=log.exprs, cond=cond)
    
  }
  
  print("Removing shapeData. (Done)")
  rm("data")
  res
}


