# Function to derive pseudo-observations for life years lost with left truncated data
#
# Required to run example code in file "SII_LYL_example.R" for the paper
#
# "A note on the measurement of socioeconomic inequalities in life years lost by cause of death" by
# A. Latouche, PK Andersen, G. Rey, M. Moreno-Betancur
#
# Thanks to Maja Pohar-Perme for the fruitlfull discussion and the R function 
#

mycmp.ps <-  function(start, stop,event,failcode=0,cause,tis,check.times) 
  #cause= the cause of interest. 
{
  
  data <- data.frame(start=start,stop=stop,event=event)
  
  
  if(missing(tis))   tis <- sort(unique(c(data$start,data$stop,check.times)))
  
  k <- length(tis)
  out <- NULL
  out$time <- tis
  
  dnie <- dnit <- yit <- rep(NA,length(tis))
  dnie.mat <- dnit.mat <- yit.mat <- matrix(tis,nrow=nrow(data),ncol=length(tis))
  for(it in 1:length(tis)){ #go through all the time points
    
    yit[it] <- sum(data$start <= tis[it] & data$stop>=tis[it])
    
    #number of events of this cause
    dnie[it] <- sum(data$start <= tis[it] & data$stop==tis[it] & data$event==cause)
    #total number of events of all causes
    dnit[it] <- sum(data$start <= tis[it] & data$stop==tis[it] & data$event!=failcode)
    
    yit.mat[,it] <- (data$start <= tis[it] & data$stop>=tis[it])
    
    #number of events of this cause
    dnie.mat[,it] <- (data$start <= tis[it] & data$stop==tis[it] & data$event==cause)
    #total number of events of all causes
    dnit.mat[,it] <- (data$start <= tis[it] & data$stop==tis[it] & data$event!=failcode)
  }
  
  #yit minus i-th value
  yit.mi <- matrix(yit,nrow=nrow(data),ncol=length(tis),byrow=T) - yit.mat
  dnie.mi <- matrix(dnie,nrow=nrow(data),ncol=length(tis),byrow=T) - dnie.mat
  dnit.mi <- matrix(dnit,nrow=nrow(data),ncol=length(tis),byrow=T) - dnit.mat
  
  dLambdae <- dnie/yit
  dLambdat <- dnit/yit
  dLambdao <- (dnit-dnie)/yit #other cause
  
  dLambdae.mat <- dnie.mi/yit.mi
  dLambdat.mat <- dnit.mi/yit.mi
  dLambdao.mat <- (dnit.mi-dnie.mi)/yit.mi #other cause
  
  
  #if yit equals 0, we assume that the hazard equals 0 - needed for bootstrap
  dLambdae[yit==0] <- 0
  dLambdat[yit==0] <- 0
  dLambdao[yit==0] <- 0
  St <- cumprod(1-dLambdat)  					  #S(t)
  Stprej <- c(1,St[-length(St)])						  #S(t-)
  
  cumince <- cumsum(Stprej*dLambdae)
  cuminco <- cumsum(Stprej*dLambdao)
  
  dolzine <- c(diff(tis),0)
  
  dareae <- dolzine*cumince
  dareat <- dolzine*St
  dareao <- dolzine*cuminco
  
  area.e <- cumsum(dareae)
  
  
  #matrix version - for pseudo-observations:  
  
  #if yit equals 0, we assume that the hazard equals 0 - needed for bootstrap
  dLambdae.mat[,yit==0] <- 0
  dLambdat.mat[,yit==0] <- 0
  dLambdao.mat[,yit==0] <- 0
  
  
  St.mat <- t(apply(1-dLambdat.mat,1,cumprod))      			  #S(t)
  Stprej.mat <- cbind(rep(1,nrow(St.mat)),St.mat[,-ncol(St.mat)])						  #S(t-)
  
  cumince.mat <- t(apply(Stprej.mat*dLambdae.mat,1,cumsum))
  cuminco.mat <- t(apply(Stprej.mat*dLambdao.mat,1,cumsum))
  
  dareae.mat <- matrix(dolzine,nrow=nrow(data),ncol=length(tis),byrow=T)*cumince.mat
  dareat.mat <- matrix(dolzine,nrow=nrow(data),ncol=length(tis),byrow=T)*St.mat
  dareao.mat <- matrix(dolzine,nrow=nrow(data),ncol=length(tis),byrow=T)*cuminco.mat
  
  area.e.mat <- t(apply(dareae.mat,1,cumsum))
  
  
  
  
  #izpis
  
  out$est <- cumince
  
  out$area <- c(0,area.e[-length(area.e)])
  
  
  out$est.mi <- cumince.mat
  out$area.mi <- cbind(0,area.e.mat[,-ncol(area.e.mat)])
  
  
  out$area[yit==0] <- NA
  out$est[yit==0] <- NA
  
  out$area.mi[,yit==0] <- NA
  out$est.mi[,yit==0] <- NA
  
  howmany <- nrow(data)
  ps.area <- howmany*matrix(area.e,nrow=nrow(data),ncol=length(tis),byrow=T) - (howmany-1)*area.e.mat
  
  out$ps.area <- ps.area
  
  out
}
