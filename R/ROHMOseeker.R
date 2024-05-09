#' Detect ROH and mCA regions from BAF and LRR
#' 
#' @param sub name of the sample to process
#' @param bafsubs Data frame with BAF of the sample to process for each probe. First to third columns should contain probe names, chromosome and position.
#' Column name of the column containing the BAF values should be the same as specified in 'sub'. Extra columns are ignored.
#' @param lrrsubs Data frame with LRR of the sample to process for each probe. First to third columns should contain probe names, chromosome and position.
#' Column name of the column containing the LRR values should be the same as specified in 'sub'. Extra columns are ignored.
#' @param winhet Window size for BAF
#' @param winhetrel Window size for Hr 
#' @param winlrr Window size for LRR
#' @param plot If \code{TRUE}, then plot parameters are returned in the output to use with \code{\link{ROHMOplotter}()}
#' @return A list with ROH and mCA detected coordinates. If specified, also plot parameters
#' 
#' @export
ROHMOseeker <- function(sub, bafsubs, lrrsubs, winhet=1000, winhetrel=500,
                      winlrr=1000, plot=FALSE){

  chrs <- as.numeric(as.character(bafsubs$Chromosome)[!bafsubs$Chromosome%in%c("X","Y")])
  tbchr <- table(chrs)
  
  ## keep autosomes
  chrcuts <- bafsubs[!bafsubs$Chromosome %in% c("X","Y"),2]
  chrcuts <- which(diff(as.numeric(as.character(chrcuts)))==1)
  bafsub <- bafsubs[!bafsubs$Chromosome %in% c("X","Y"),sub]
  
  
  #total loss of heterozygosity
  selindbaf <- 1:length(bafsub)
  het <- -0.5*abs(2*bafsub-1)+0.5
  nn <- !is.na(het)
  het <- het[nn]
  selindbaf <- selindbaf[nn]
  hetav <- as.vector(ma(het, winhet))
  ev <- rep(NA, length(bafsub))
  selend <- ((length(hetav)+1):length(selindbaf))
  ev[selindbaf[-selend]] <- hetav
  
  #remove chromosome ends
  if(length(tbchr)!=1)
  {
    chrends <- rbind(chrcuts-winhet, chrcuts)
    for(cc in 1:(length(tbchr)-1)){ev[chrends[1,cc]:chrends[2,cc]] <- NA}
  }
  cuttot <- 0.02
  part <- ev > cuttot
  part[c(1,chrcuts, length(part))] <- TRUE
  part[chrcuts+1] <- TRUE
  locnna <- which(!is.na(part))
  dp <- diff(part[locnna])
  rb <- locnna[which(dp==1)]
  lb <- locnna[which(dp==-1)]
  
  if(part[!is.na(part)][length(part[!is.na(part)])]==FALSE)
    rb <- c(rb, length(part))
  if(part[!is.na(part)][1]==FALSE)
    lb <- c(1, lb)
  bp <- rbind(lb,rb)
  
  
  #relative loss of heterozygosity (baf split)
  whhet <- which(het > 0.2)
  selindbafmos <- selindbaf[whhet]
  hetavrel <- as.vector(ma(het[whhet], winhetrel))
  evrel <- rep(NA, length(bafsub))
  selend <- ((length(hetavrel)+1):length(selindbafmos))
  evrel[selindbafmos[-selend]] <- hetavrel
  
  #remove chromosome ends
  if(length(tbchr)!=1)
  {
    chrends <- rbind(chrcuts-winhetrel, chrcuts)
    for(cc in 1:(length(tbchr)-1)){evrel[chrends[1,cc]:chrends[2,cc]] <- NA}
  }
  cutrel <- 0.5-0.05
  partrel <- evrel > cutrel
  partrel[c(1,chrcuts, length(partrel))] <- TRUE
  partrel[chrcuts+1] <- TRUE
  locnna <- which(!is.na(partrel))
  dprel <- diff(partrel[locnna])
  rbrel <- locnna[which(dprel==1)]
  lbrel <- locnna[which(dprel==-1)]
  
  if(partrel[!is.na(partrel)][length(partrel[!is.na(partrel)])]==FALSE)
    rbrel <- c(rbrel, length(partrel))

  if(partrel[!is.na(partrel)][1]==FALSE)
    lbrel <- c(1, lbrel)
  
  bprel <- rbind(lbrel,rbrel)
  loh <- cbind(bafsubs[lb,2:3],bafsubs[rb,2:3])
  mos <- cbind(bafsubs[lbrel,2:3],bafsubs[rbrel,2:3])
  
  
  #change in lrr of heterozygosity (baf split)
  lrrsub<- lrrsubs[!lrrsubs$Chromosome %in% c("X","Y"),sub]
  selindlrr <- 1:length(lrrsub)
  llrcent <- (lrrsub-mean(lrrsub))
  lrrnorm <- llrcent/max(abs(llrcent))
  nn <- !is.na(lrrnorm)
  lrrnorm <- lrrnorm[nn]
  selindlrr <- selindlrr[nn]
  lrrav <- as.vector(ma(lrrnorm, winlrr))
  evlrr <- rep(NA, length(lrrsub))
  selend <- ((length(lrrav)+1):length(selindlrr))
  evlrr[selindlrr[-selend]] <- lrrav
  
  #remove chromosome ends
  if(length(tbchr)!=1)
  {
    chrends <- rbind(chrcuts-winlrr, chrcuts)
    for(cc in 1:(length(tbchr)-1)){evlrr[chrends[1,cc]:chrends[2,cc]] <- NA}
  }
  
 
  
  if(length(tbchr)==1) {
    pos <- bafsubs[,3]/10^6
    xax <- "s"
    bp <- rbind(bafsubs[lb,3],bafsubs[rb,3])/10^6
    bprel <- rbind(bafsubs[lbrel,3],bafsubs[rbrel,3])/10^6
    xxlab <- paste("Chr:", " (Mb)", sep=names(tbchr))
  }
  
  if(plot){
    plotlist <- list(sub=sub, bafsub=bafsub, lrrnorm=lrrnorm, chrs=chrs, ev=ev, evrel=evrel,
                     evlrr=evlrr, bp=bp, bprel=bprel, cuttot=cuttot, cutrel=cutrel)
  }
  else{
    plotlist <- list()
  }
 
   list(loh=loh, mos=mos, plot=plotlist)
}


ma <- function(x, n = 5){
  cx <- c(0, cumsum(ifelse(is.na(x), 0, x)))
  cn <- c(0, cumsum(ifelse(is.na(x), 0, 1)))
  rx <- cx[(n+1):length(cx)] - cx[1:(length(cx) - n)]
  rn <- cn[(n+1):length(cx)] - cn[1:(length(cx) - n)]
  rx / rn
}


