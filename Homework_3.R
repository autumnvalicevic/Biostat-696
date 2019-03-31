library(maps)
library(maptools)
library(spdep)
library(classInt)
library(RColorBrewer)
library(SpatialEpi)
library(nlme)
library(CARBayes)

# Spatial Analysis of 2016 Election Data

# Download data
election <- read.csv("~/Downloads/election_data_2016.csv")

## Here we are getting the information from the maps package in R on the United States
us <-map("state",fill=T,plot=F)
## This is to get the IDs of the counties in NC that are useful to make the spatial polygon
us.IDs <- sapply(strsplit(us$names,":"),function(x) x[1])
us.poly<- map2SpatialPolygons(us,IDs=us.IDs,proj4string=CRS("+proj=longlat + datum=wgs84"))
plot(us.poly,col="white",axes=TRUE)

row.names(us.poly)
length(row.names(us.poly))

# Remove Alaska and Hawaii
election <- election[c(1,3:11,13:51),]

# Tell R which variable to plot and Quintiles
plotvar <- election$pct.Trump
nclr <- 5

plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=seq(min(plotvar),max(plotvar),length=nclr+1))
class

### To see which are the class breaks (we will need this to make the legend)
colcode <- findColours(class,plotclr)


plot(us.poly,border="black",axes=TRUE,xlim=c(-123,-68))
title(main="Percentage Trump Voters")
plot(us.poly,col=colcode,add=T)

leg.txt<-c("[4.2%,17.38%)","[17.38%,30.56%)","[30.56%,43.74%)","[43.74%,56.92%)",
           "[56.92%,70.1%]")
legend(-125,34,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")


#  Computing Morans I
us.nb <- poly2nb(us.poly)
us.weights <- nb2WB(us.nb)
adj.us <- us.weights$adj
weights.us <- us.weights$weights
num.us <- us.weights$num
us.listw <- nb2listw(us.nb)

## The function moran.test computes Moran I for areal data and returns the observed
## value of I, the expected value under the hypothesis of iid data and 
## also the asymptotic variance (we saw the formula in class)
## If we specify randomisation=FALSE, the function moran.test also
## runs a hypothesis test on I using the asymptotic distribution of I
moran.election <- moran.test(election$pct.Trump,listw=us.listw,randomisation=FALSE)
moran.election


## Linear model to explain the percent of Trump Voters
summary(lm(pct.Trump ~ Gini.coefficient + median.household.income.2016, data = election))
election.lm <- lm(pct.Trump ~ Gini.coefficient + median.household.income.2016, data = election)

## Plot of Gini Coefficient
plotvar <- election$Gini.coefficient
nclr <- 5

plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=seq(min(plotvar),max(plotvar),length=nclr+1))
class

### To see which are the class breaks (we will need this to make the legend)
colcode <- findColours(class,plotclr)


plot(us.poly,border="black",axes=TRUE,xlim=c(-123,-68))
title(main="Gini Coefficient")
plot(us.poly,col=colcode,add=T)

leg.txt<-c("[0.426,0.4492)", "[0.4492,0.4724)", "[0.4724,0.4956)", "[0.4956,0.5188)",  "[0.5188,0.542]")
legend(-125,34,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")


## Plot of median household income in 2016
plotvar <- election$median.household.income.2016
nclr <- 5

plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=seq(min(plotvar),max(plotvar),length=nclr+1))
class

### To see which are the class breaks (we will need this to make the legend)
colcode <- findColours(class,plotclr)


plot(us.poly,border="black",axes=TRUE,xlim=c(-123,-68))
title(main="Median Household Income")
plot(us.poly,col=colcode,add=T)

leg.txt<-c("[41754,49192.2)", "[49192.2,56630.4)", "[56630.4,64068.6)", "[64068.6,71506.8)",   "[71506.8,78945]")
legend(-125,34,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")


## Moran's I for Residuals of the linear model to see if there is still spatial correlation
lm.morantest(lm(election$pct.Trump~election$median.household.income.2016 + election$Gini.coefficient),listw=us.listw)

res.election<- as.numeric(summary(election.lm)$residuals)
moran.test(res.election, us.listw, rank=TRUE)

## Fit Improper CAR Model
N <- length(election$pct.Trump)
rep.us <- rep(1:N,num.us)
W <- matrix(0,N,N)
for(i in 1:N){
  W[i,adj.us[rep.us==i]] <- rep(1,num.us[i])
}

formula <- election$pct.Trump~election$median.household.income.2016+election$Gini.coefficient
model.car <- S.CARleroux(formula=formula, W=W, family="gaussian",
                         rho=1,burnin=20000, n.sample=40000,thin=1, prior.mean.beta=NULL, prior.var.beta=NULL,prior.nu2=NULL,prior.tau2=NULL, verbose=TRUE)

model.car$summary.results

samples.eta <- model.car$samples$phi
post.median.eta <- as.numeric(apply(samples.eta,2,median))
plotvar <- post.median.eta
nclr <- 5

plotclr <- brewer.pal(nclr,"RdBu")[5:1]
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(min(post.median.eta),-0.01359792,-0.01020851,0.01084249,0.01440649,max(post.median.eta)),length=nclr+1)
class
quint <- quantile(post.median.eta, seq(0, 1, 0.2))

colcode <- findColours(class,plotclr)

plot(us.poly,border="black",axes=TRUE,xlim=c(-123,-68))
title(xlab="Longitude",ylab="Latitude",main="Spatial random effects \n CAR model")
plot(us.poly,col=colcode,add=T)

leg.txt<-c("[-0.005875403,-0.001823303)", "[-0.001823303,-0.0001147215)", "[-0.0001147215,0.0007020874)",   "[0.0007020874,0.001529794)","[0.001529794,0.008216111]" )
legend(-125,34,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

### Upper bd 95% CI of spatial effects
upp.bd95ci.eta <- as.numeric(apply(samples.eta,2,quantile,0.975))
plotvar <- upp.bd95ci.eta
nclr <- 5


plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(min(upp.bd95ci.eta),0.2504866,0.2752785,0.3086177,0.3819294,max(upp.bd95ci.eta)),length=nclr+1)

colcode <- findColours(class,plotclr)

plot(us.poly,border="black",axes=TRUE,xlim=c(-123,-68))
title(xlab="Longitude",ylab="Latitude",main="Upper bound of 95% CI for spatial random effects \n CAR model")
plot(us.poly,col=colcode,add=T)

leg.txt<-c("[-18.64418,-5.543735)", "[-5.543735,-1.418444)",  "[-1.418444,1.836683)","[1.836683,17.35417)","[17.35417,7565589]")
legend(-125,34,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")


### Lower bd 95% CI of spatial effects
low.bd95ci.eta <- as.numeric(apply(samples.eta,2,quantile,0.025))
plotvar <- low.bd95ci.eta
nclr <- 5

plotclr <- brewer.pal(nclr,"Blues")[5:1]
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(min(low.bd95ci.eta),-0.3931782,-0.3294695,-0.2837618,-0.2487109,max(low.bd95ci.eta)),length=nclr+1)

colcode <- findColours(class,plotclr)

plot(us.poly,border="black",axes=TRUE,xlim=c(-123,-68))
title(xlab="Longitude",ylab="Latitude",main="Lower bound of 95% CI for spatial random effects \n CAR model")
plot(us.poly,col=colcode,add=T)

leg.txt<-c("[-27.2412,-8.659667)", "[-8.659667,-3.809636)", "[-3.809636,0.04388993)","[0.04388993,4.205935)",    "[4.205935,12.55968]")
legend(-125,34,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

us.w.proper.car <- nb2listw(us.nb,style="B")



## Fit Proper CAR Model
model.proper.car <- spautolm(election$pct.Trump~election$median.household.income.2016+election$Gini.coefficient,listw=us.w.proper.car,family="CAR")
summary(model.proper.car)

## Fit SAR Model
model.sar <- spautolm(election$pct.Trump~election$median.household.income.2016+election$Gini.coefficient,listw=us.listw,family="SAR")
summary(model.sar)



## Using Disease Mapping Techniques to fit model
N <- length(election$Trump.1)
rep.us <- rep(1:N,num.us)
W <- matrix(0,N,N)
for(i in 1:N){
  W[i,adj.us[rep.us==i]] <- rep(1,num.us[i])
}

election$Trump.1 <- round(election$Trump / 1000, digits = 0)

formula <- election$Trump.1~election$median.household.income.2016+election$Gini.coefficient+offset(log(election$No_voters))

model.car.1 <- S.CARleroux(formula=formula, W=W, family="poisson",
                           rho=1,burnin=20000, n.sample=60000,thin=1, prior.mean.beta=NULL, prior.var.beta=NULL,prior.nu2=NULL,prior.tau2=NULL, verbose=TRUE)

model.car$summary.results




## Propensity to vote Republican 

## STILL WORKING ON
samples <- model.car.1$samples$phi
post.median.eta <- as.numeric(apply(samples$psi,2,median))
plotvar <- post.median.eta
nclr <- 5

plotclr <- brewer.pal(nclr,"RdBu")[5:1]
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(min(post.median.eta),-0.02,-0.004,0,0.0,0.006,max(post.median.eta)),length=nclr+1)

colcode <- findColours(class,plotclr)

plot(mi.poly.county,border="black",axes=TRUE,xlim=c(-123,-68))
title(xlab="Longitude",ylab="Latitude",main="Spatial random effects \n CAR model")
plot(mi.poly.county,col=colcode,add=T)

leg.txt<-c("[-0.04;-0.02)","[-0.02;-0.004)","[-0.004;0.0)","[0.0;0.006)","[0.006;0.017]")
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")


### Posterior SD of spatial effects
post.sd.eta <- as.numeric(apply(samples.eta,2,sd))
plotvar <- post.sd.eta
nclr <- 5

plotclr <- brewer.pal(nclr,"YlGnBu")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(0.012,0.013,0.014,0.016,0.017,0.023),length=nclr+1)

colcode <- findColours(class,plotclr)

plot(mi.poly.county,border="black",axes=TRUE,xlim=c(-123,-68))
title(xlab="Longitude",ylab="Latitude",main="SD of spatial random effects \n CAR model")
plot(mi.poly.county,col=colcode,add=T)

leg.txt<-c("[0.012; 0.013)","[0.013; 0.014)","[0.014; 0.016)","[0.016; 0.017)","[0.017; 0.023]")
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")



### Lower bd 95% CI of spatial effects
low.bd95ci.eta <- as.numeric(apply(samples.eta,2,quantile,0.025))
plotvar <- low.bd95ci.eta
nclr <- 5

plotclr <- brewer.pal(nclr,"Blues")[5:1]
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(-0.082,-0.049,-0.032,-0.027,-0.023,-0.011),length=nclr+1)

colcode <- findColours(class,plotclr)

plot(mi.poly.county,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="Longitude",ylab="Latitude",main="Lower bound of 95% CI for spatial random effects \n CAR model")
plot(mi.poly.county,col=colcode,add=T)

leg.txt<-c("[-0.082; -0.049)","[-0.049; -0.032)","[-0.032; -0.027)","[-0.027; -0.023)","[-0.023; -0.011]")
legend(-125,34,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")


### Upper bd 95% CI of spatial effects
upp.bd95ci.eta <- as.numeric(apply(samples.eta,2,quantile,0.975))
plotvar <- upp.bd95ci.eta
nclr <- 5


plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(-0.003,0.0,0.025,0.027,0.035,0.057),length=nclr+1)

colcode <- findColours(class,plotclr)

plot(mi.poly.county,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="Longitude",ylab="Latitude",main="Upper bound of 95% CI for spatial random effects \n CAR model")
plot(mi.poly.county,col=colcode,add=T)

leg.txt<-c("[-0.003; 0.0)","[0.017; 0.025)","[0.025; 0.027)","[0.027; 0.035)","[0.035; 0.057]")
legend(-125,34,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

