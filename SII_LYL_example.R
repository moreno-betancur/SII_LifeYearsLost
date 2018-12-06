# Example code for the paper:
#
# "A note on the measurement of socioeconomic inequalities in life years lost by cause of death" by
# A. Latouche, PK Andersen, G. Rey, M. Moreno-Betancur
#
# Thanks to Maja Pohar-Perme for the fruitlfull discussion and the R function 
#
# We use the simulated data set (dataSII.txt) previoulsy considered in Moreno-Betancur et al. Epidemiology (2015), 
# to estimate SII in life years lost, and the decomposition of life years lost.

#### Load pseudo-observation function and dataset ----

source("mycmp.R")

dat<-read.table("dataSII.txt",dec=".",sep=";",header=T)
names(dat)
summary(dat$ageEntry) # Age at entry (more that 30 yrs old)
summary(dat$ageExit)
tau<-80  # Set upper limit. In general, one would plot the survival function to inform the choice of the upper limit tau

#### Life years lost, all-cause mortality ----

# Derive pseudo-observations

ys<-mycmp.ps(start = dat$ageEntry,stop =dat$ageExit,event = dat$status,cause=1,check.times=tau)
# print(dim(ys$ps.area))
# v1<-sort(unique(c(dat$ageEntry,dat$ageExit)))
# length(v1)
# We only need the pseudo observation at time tau

pseu<-ys$ps.area[,ys$time==tau]

# Use the pseudo observation as a response in a linear model with "midpoint" (the socioeconomic rank) as predictor

fit<-lm(pseu~midpoint,data=dat)
allcauseSII<-fit$coefficients[2]
print(fit)
print(confint(fit))

# Which is interpreted as follows: 
# "On average across the socioeconomic scale, the least educated lost 6.69 years of life compared to the most educated."


#### Decomposition by cause of death ----

# Just for illustration, we create artificial causes of death by setting the cause to be 2 if the (usually unobserved) x is 
# greater that the median, and 1 otherwise

dat$status2<-ifelse(dat$x>0.76 & dat$status==1,2,dat$status)

# We can loop over the 2 competing risks to get the cause-specific SII

v<-c()
for (cause in 1:2){
  
  ys<-mycmp.ps(start = dat$ageEntry,stop =dat$ageExit,event = dat$status2,cause=as.numeric(cause),check.times=tau)
  pseu<-ys$ps.area[,ys$time==tau]
 
  fit<-lm(pseu~midpoint,data=dat)
  v[cause]<-fit$coefficients[2]
  print(fit$coefficients[2])
  print(confint(fit))
  
}


# The cause-specific SII  is years lost by cause sum up to the all-cause SII 

sum(v)
allcauseSII


#### Crude estimates of years lost by socioeconomic status ----

#  The years lost up to 80 years for each socioeconomic group is computed as:
#  50 - RestrictedMean(30,80)

library(prodlim)

# We loop over the 4 social groups in the simulated dataset

table(dat$socgroup)
out1<-c()
for (i in 1:4)
{
  dat.tmp<-subset(dat, socgroup==i)
  res<-prodlim(Hist(time=ageExit,event=status,entry=ageEntry)~1,data=dat.tmp)
  fres<-stepfun(res$time[-1],1-res$surv)
  
  intfres<-integrate(fres,lower = min(dat.tmp[,7]),upper=tau,subdivisions=length(res$time[-1]))
  
  # print(paste("Ed  Level", i, sep=":"))
  #  print(intfres$value)
  
  
  out1[i]<-(tau-30)-intfres$value
  #   
  # end edlevel      
}
print(out1)