#' Plot genome-wide BAF, LRR and Hr computed by ROHMOseeker 
#'
#' @param hetseries List obtained with ROHMOseeker containing detected ROH, mCA and plot parameters
#' @param dir Direcory where the plot file shoud be created
#' @export
ROHMOplotter <- function(hetseries, dir){
  
  stopifnot(length(hetseries$plot)>0)
  
  sub <- hetseries$plot$sub
  bafsub <- hetseries$plot$bafsub
  lrrnorm <- hetseries$plot$lrrnorm
  chrs <- hetseries$plot$chrs
  ev <- hetseries$plot$ev
  evrel <- hetseries$plot$evrel
  evlrr <- hetseries$plot$evlrr
  bp <- hetseries$plot$bp
  bprel <- hetseries$plot$bprel
  cl <- hetseries$plot$cl
  cuttot <- hetseries$plot$cuttot
  cutrel <- hetseries$plot$cutrel
  
  
  pos <- 1:length(bafsub)
  xax <- "n"
  xxlab <- paste("Genomic position")
  cols <- rep("red", length(chrs))
  cols[chrs%in%seq(2,length(chrs),2)] <- "orange"
  
  tiff(filename = paste0(dir, '/', sub, ".tif"), height = 1500, width = 2000, res = 200)

  
  plot(pos, lrrnorm, ylab="BAF", xaxt = xax, xlab=xxlab, main= sub,
       col="grey", pch=".",cex=1.3, ylim=c(-1,1))
  
  points(pos, 2*bafsub-1, col=cols, pch=".")
  lines(pos, ev,type="l")
  lines(pos[!is.na(evrel)],evrel[!is.na(evrel)], type="l")
  lines(pos, evlrr,type="l", col="green")
  
  if(length(bp)!=0) {
    for(cl in 1:ncol(bp))
      lines(bp[,cl], c(0,0), col="blue",lwd=3)
  }
  
  if(length(bprel)!=0) {
    for(cl in 1:ncol(bprel))
      lines(bprel[,cl], c(0.5,0.5), col="blue",lwd=3)
  }
  
  abline(h=cuttot,lwd=1.5, lty=2)
  abline(h=cutrel,lwd=1.5, lty=2)
  
  dev.off()
  
}
