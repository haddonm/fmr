



data(ocaa)
data(fish)
data(owaa)
data(param)
pars <- param[1:28,2]
yrs <- as.numeric(fish[,"year"])
ages <- as.numeric(2:10)
onaa <- ocaa[,ages]
owa <- as.matrix(owaa[,ages])
nage <- length(ages)
getssq2(pars,M=0.2,yrs,ages,onaa,owa,cpue=fish[,"obsce"]) #should be 332.5389


# Using aspm--------------------------------------------------

library(fmr)
library(codeutils)
library(hplot)
library(knitr)
library(datalowSA)

dbdir <- getDBdir()
rundir <- pathtopath(dbdir,"A_Code/fmr/data-raw")

fishdat <- readdata(pathtopath(rundir,"westroughy.csv"),verbose=FALSE)
fish <- fishdat$fish
props <- fishdat$props # length-, weight-, maturity- and selectivity-at-age
glb <- fishdat$glb

plotprep(width=9, height=6)
pars <- c(13.9,0.25)
aspmLL(pars,infish=fish,inglb=glb,inprops=props)
fishery <- dynamics(pars,infish=fish,inglb=glb,inprops = props)
plotASPM(fishery)


pars <- c(13.8,0.4)
ans <- fitASPM(pars,minfun=aspmLL,dynfun=dynamicsH,infish=fish,
               inglb=glb,inprops=props)
outfit(ans) # a tidier way of printing the list output from optim
fishery <- dynamicsH(ans$par,infish=fish,inglb=glb,inprops = props)

kable(fishery,digits=c(0,3,3,3,3,3,4,4,3))


ceCI <- getLNCI(fishery[,"PredCE"],ans$par[2],P=95)
plotASPM(fishery,CI=ceCI)

# two parameter model--------------------------------
data(fishdat)
fish <- fishdat$fish
props <- fishdat$props # length-, weight-, maturity- and selectivity-at-age


(glb <- fishdat$glb)

plotprep(width=9, height=6)
pars <- c(13.9,0.25)
aspmLL(pars,infish=fish,inglb=glb,inprops=props)
fishery <- dynamicsH(pars,infish=fish,inglb=glb,inprops = props)
plotASPM(fishery)



pars <- c(13.7,0.25)
ans <- fitASPM(pars,minfun=aspmLL,dynfun=dynamicsH,infish=fish,
               inglb=glb,inprops=props)
outfit(ans) # a tidier way of printing the list output from optim
fishery <- dynamicsH(ans$par,infish=fish,inglb=glb,inprops = props)

kable(fishery,digits=c(0,3,3,3,3,3,4,4,3))


ceCI <- getLNCI(fishery[,"PredCE"],ans$par[2])
plotASPM(fishery,CI=ceCI)

#Three parameter model-------------------------------------------
data(dataspm)
fish <- dataspm$fish
glb <- dataspm$glb
props <- dataspm$props
pars <- c(14,0.19,0.95) # Fit 3 par__aspm__with penalty
# pars <- c(13.2794439,0.1731744,0.4933178) # for a second time through
scalepar <- magnitude(pars)
bestL <- optim(pars,aspmPENLL,method="Nelder-Mead",
               infish=fish,inglb=glb,inprops=props,
               control=list(maxit = 1000,parscale=scalepar))
outfit(bestL)
fisheryPen <- dynamics(bestL$par,infish=fish,inglb=glb,inprops=props)
ceCI <- getLNCI(fisheryPen[,"PredCE"],bestL$par[2])
plotASPM(fisheryPen,CI=ceCI)



#test robustness------------------------------------------------

set.seed(12335)  # to get repeatable results, normally you would not do this
data(fishdat)
fish <- fishdat$fish
glb <- fishdat$glb
props <- fishdat$props
pars <- c(14,0.3)
out <- robustASPM(pars,fish,glb,props,scaler=20,N=15,console=FALSE)
str(out)
print(round(out$results,4))





set.seed(12235)  # to get repeatable results, normally you would not do this
data(dataspm)
fish <- dataspm$fish
glb <- dataspm$glb
props <- dataspm$props
pars <- c(14.0,0.2,0.6)
out <- robustASPM(pars,fish,glb,props,scaler=15,N=10,console=FALSE)
print(round(out$results,3))
print(round(out$range,3))


cor(out$results[,c("LnR0","Depl","-veLL","MSY")])  # correlations between outputs
#plotprep(width=8,height=6)
intensity <- 2   #  how many points overlapping = maximum colour
plotprep(width=8, height=7)
pairs(out$results[,c("LnR0","Depl","-veLL","MSY")],pch=16,
      col=rgb(1,0,0,1/intensity),font=7,font.labels = 7)

# production curve and statistics-------------------------------------

data(dataspm)
fish <- dataspm$fish
glb <- dataspm$glb
props <- dataspm$props
pars <- c(13.75,0.189667,0.6) # Fit 3 par__aspm__with penalty
bestL <- optim(pars,aspmPENLL,method="Nelder-Mead",
               infish=fish,inglb=glb,inprops=props,
               control=list(maxit = 1000, parscale = c(10,1,0.1)))
# two times through
bestL <- optim(bestL$par,aspmPENLL,method="Nelder-Mead",
               infish=fish,inglb=glb,inprops=props,
               control=list(maxit = 1000, parscale = c(10,1,0.1)))
par <- bestL$par
print(par)
prod <- getProduction(exp(par[1]),fish,glb,props,
                      Hrg=c(0.01,0.45,0.005),nyr=50)
head(round(prod,3),6)
tail(round(prod,3),6)
anspen <- prodASPM(prod,target=0.48,console=FALSE,plot=TRUE)
round(anspen,3)


# a phase plot--------------------------------------

fisheryPen <- dynamics(bestL$par,infish=fish,inglb=glb,inprops=props)
outs <- aspmphaseplot(fisheryPen,prod,anspen,Blim=0.2,fnt=7)


# bootstrapping -------------------------------------------

data(dataspm)
fish <- dataspm$fish
glb <- dataspm$glb
props <- dataspm$props
pars <- c(13.5,0.18,0.5)
bestL <- fitASPM(pars,minfun=aspmLL,dynfun=dynamicsH,fish,glb,props)
fishery <- dynamics(bestL$par,fish,glb,props)
kable(fishery,digits=c(0,1,1,3,3,3,3,3,3))


reps <- 100
starttime <- Sys.time()
answer <- bootASPM(fish,glb,props,bestL$par,iter=reps)
Sys.time() - starttime
str(answer,max.level=1)




yrs <- fishery[,"Year"]
nyrs <- length(yrs)
par(mfrow=c(2,2),mai=c(0.45,0.45,0.05,0.05)) 
par(cex=0.85, mgp=c(1.35,0.35,0), font.axis=7,font=7,font.lab=7)  
label <- names(answer$result[1,1,])
label <- label[-3]  # remove CPUE
numvar <- length(label)
bootvar <- answer$result[,nyrs,label[1]]
for (i in 1:numvar) { # i=3
  bootvar <- answer$result[,nyrs,label[i]]
  quantCI <- quantile(bootvar,probs=c(0.05,0.5,0.95),na.rm=TRUE)
  hist(bootvar,breaks=30,main="",xlab=label[i],col="red")
  abline(v=quantCI,col=c(4,4,4),lwd=c(1,2,1))
}


pickvar <- "Deplete"
bootvar <- answer$result[,,pickvar]
yrs <- as.numeric(colnames(bootvar))
nyrs <- length(yrs)
quantCI <- t(apply(bootvar,2,quants))
kable(quantCI,digits=c(3,3,3,3,3,3))




ymax <- getmax(bootvar)
par(mfrow=c(1,1),mai=c(0.45,0.45,0.05,0.05)) 
par(cex=0.85, mgp=c(1.35,0.35,0), font.axis=7,font=7,font.lab=7)
plot(yrs,bootvar[1,],type="n",lwd=1,col=0,ylim=c(0,ymax),
     panel.first = grid(),xlab="",ylab=pickvar)
for (i in 1:reps) lines(yrs,bootvar[i,],lwd=1,col="grey")
lines(yrs,quantCI[,"50%"],lwd=2,col="red")
arrows(x0=yrs,y0=quantCI[,"5%"],y1=quantCI[,"95%"],
       col=2,lwd=1,length=0.035,angle=90,code=3)




# Log-Normal Priors--------------------------------------------------

library(codeutils)
library(hplot)

x <- seq(1:1000)
n <- length(x)
y <- rep(1/n,n)  # rep(0.9,n)#
ycum <- cumsum(y)

plotprep(width=9, height=9)
parset(plots=c(3,1))
plot1(x,y,lwd=2,defpar=FALSE)
plot1(x,ycum,lwd=2,defpar=FALSE)
plot1(log(x),ycum,lwd=2,defpar=FALSE)

# The cumulative distribution of a uniform distribution across a range of values 
# for a variable should be a straight line from the minimum up to one.
# log-transform the variable and the uniformity disappears. 

xL <- seq(1,6.9,length=1000)
yL <- rep(1/n,n)
ycumL <- cumsum(y)

plotprep(width=9, height=9)
parset(plots=c(3,1))
plot1(xL,yL,lwd=2,xlab="Log-scale",defpar=FALSE)
plot1(xL,ycumL,lwd=2,xlab="Log-scale",defpar=FALSE)
plot1(exp(xL),ycumL,lwd=2,xlab="Linear-scale",defpar=FALSE)




# Step 1: Define the log-normal distribution parameters
meanlog <- 0  # Mean of the log of the distribution
sdlog <- 1    # Standard deviation of the log of the distribution

# Generate some log-normal data
set.seed(123)
data <- rlnorm(100, meanlog, sdlog)

# Step 2: Define the uniform prior
# Assuming the prior is uniform between 0 and 10
prior_min <- 0
prior_max <- 10

# Step 3: Combine prior and likelihood
# For simplicity, we use a grid approximation for the posterior

# Define a grid of possible values for the parameter
param_grid <- seq(prior_min, prior_max, length.out = 1000)

# Calculate the likelihood for each value in the grid
likelihood <- sapply(param_grid, function(x) prod(dlnorm(data, meanlog = log(x), sdlog = sdlog)))

plot(param_grid, likelihood, type = "l", main = "Likelihood", xlab = "Parameter", ylab = "Likelihood")

# Calculate the prior for each value in the grid
prior <- dunif(param_grid, min = prior_min, max = prior_max)
plot(param_grid, prior, type = "l", main = "prior", xlab = "Parameter", ylab = "Likelihood")


# Calculate the posterior (unnormalized)
posterior_unnorm <- likelihood * prior

# Normalize the posterior
posterior <- posterior_unnorm / sum(posterior_unnorm)

# Plot the posterior distribution
plot(param_grid, posterior, type = "l", main = "Posterior Distribution", xlab = "Parameter", ylab = "Density")





# HIMI Tags---------------------

rundir <- "C:/Users/Malcolm/Dropbox/A_CodeR/mfmur/data-raw/"
filen <- pathtopath(rundir,"himi-tags.csv")
himi <- read.csv(filen,header=TRUE)

props <- matrix(NA,nrow=11,ncol=11,dimnames=list(2012:2022,2013:2023))
for (i in 1:11) { # i = 1
  vect <- as.numeric(himi[i,"N"]-himi[i,(i+1):11])
  num <- c(himi[i,"N"],vect)
  props[i,] <- as.numeric(himi[i,2:12]/num)
}

yrs <- 2013:2023
plotprep(width=9, height=9)
parset(plots=c(1,1))
plot(yrs[2:11],props[1,2:11],type="l",lwd=2,xlab="Year",xlim=c(2013,2023),
     ylim=c(0,0.07))
for (i in 2:10) lines(yrs[1:(11-i)],props[i,(i+1):11],lwd=2)


yrs <- log(13:23)
propsL <- log(props)
propsL[1,10:11] <- NA
propsL[2,11] <- NA

plotprep(width=9, height=9)
parset(plots=c(1,1))
plot(yrs[2:11],propsL[1,2:11],type="l",lwd=2,xlab="Year",ylim=c(-6,-2),xlim=c(13,21))
for (i in 2:10) lines(yrs[1:(11-i)],propsL[i,(i+1):11],lwd=2)


proportL <- matrix(0,nrow=121,ncol=2)
yrs <- 13:22

count <- 0
for (i in 1:10) { 
  for (j in 2:11) {
    count <- count + 1
    proportL[count,] <- c(yrs[i],propsL[j,i])
  }   
}
proportL

tmp <- proportL[order(proportL[,2]),]

out <- tmp[1:73,]

plotprep(width=9, height=9)
parset(plots=c(1,1))
plot1(log(out[,1]),out[,2],type="p",pch=16)


# find F for a catch------------------------------

library(Rcpp)

fisheryPen

expB <- fisheryPen[,"ExploitB"]
catch <- fisheryPen[,"Catch"]
M <- 0.05


predC <- numeric(length(catch))
yr <- 2

Fyr <- 0.4

library(microbenchmark)

ans2 <- microbenchmark(
  out <- datalowSA::getProduction(inR0=exp(ans$par[1]),infish=fish,
                                  inglb=glb,inprops=props,
                                  Hrg=c(0.02,0.15,0.0025),maxiter=3),
  # out2 <- grsearch(matchC,interval=c(0,1.0),M=M,cyr=catch[2],Byr=expB[2],
  #                  tol = 1e-08),
  out3 <- datalowSA::getProductionC(inR0=exp(ans$par[1]),infish=fish,
                                    inglb=glb,inprops=props,
                                    Hrg=c(0.02,0.15,0.0025)),
  times=10
)
ans2

out <- optimize(matchC,interval=c(0,1.25),M=M,cyr=catch[2],Byr=expB[2],
                maximum=FALSE,tol = .Machine$double.eps^0.25)$minimum



Fvalues <- function(catches,exploitB,M,nyr,maxF=1.0,tol=1e-09) {
  yrF <- numeric(nyr)
  for (yr in 1:nyr)
    yrF[yr] <- optimize(matchC,interval=c(0,maxF),M=M,cyr=catches[yr],
                        Byr=exploitB[yr],maximum=FALSE,tol=tol)$minimum
  return(yrF)
} # end of Fvalues




ans2 <- microbenchmark(
  yrF <- Fvalues(catches=fisheryPen[2:13,"Catch"],
          exploitB=fisheryPen[1:12,"ExploitB"],M=0.05,maxF=1.0),
  yrF <- Fvalues2(catches=fisheryPen[2:13,"Catch"],
                 exploitB=fisheryPen[1:12,"ExploitB"],M=0.05,nyr=12,maxF=1.0),
  times=100
)
ans2




out <- datalowSA::getProduction(inR0=exp(ans$par[1]),infish=fish,
                                inglb=glb,inprops=props,
                                Hrg=c(0.03,0.05,0.0001),maxiter=5)

plot(as.numeric(rownames(out)),out$Yield,type="l",lwd=2,xlim=c(0.03,0.05),
     ylim=c(200,235))



pars <- c(10.25,0.25)
fishery <- dynamicsF(pars,infish=fish,inglb=glb,inprops = props)
fishery

ceCI <- getLNCI(fishery[,"PredCE"],pars[2])
plotASPM(fishery,CI=NA)


bestL <- optim(pars,aspmFLL,method="Nelder-Mead",
               infish=fish,inglb=glb,inprops=props,
               control=list(maxit = 1000, parscale = c(10,0.1)))
bestL
fisheryF <- dynamicsF(bestL$par,infish=fish,inglb=glb,inprops=props)
fisheryF
ceCI <- getLNCI(fisheryF[,"PredCE"],bestL$par[2])
plotASPM(fisheryF,CI=ceCI)

q <- exp(mean(log(fisheryF[,"CPUE"]/fisheryF[,"ExploitB"]),na.rm=TRUE))
B0 <- round(fisheryF[1,"SpawnB"])
exB0 <- round(fisheryF[1,"ExploitB"]) 
msy <- (0.027 * B0)
mcy <- msy * 2/3
{
  cat("q    = ",round(q,4),"\n")
  cat("B0   = ",round(B0,0),"\n")
  cat("exB0 = ",round(exB0,0),"\n")
  cat("MSY  = ",round(msy,0),"\n")
  cat("MCY  = ",round(mcy,0),"\n")
}

library(fmr)
library(codeutils)
library(hplot)

data(westroughy)
fish <- westroughy$fish
glb <- westroughy$glb
props <- westroughy$props
pars <- c(7,0.3)
scalepar <- magnitude(pars)
bestL <- optim(pars,aspmLL,method="Nelder-Mead",dynfun=dynamicsH,infish=fish,
               inglb=glb,inprops=props,
               control=list(maxit = 1000, parscale = scalepar))
outfit(bestL)
fishery <- dynamicsH(bestL$par,fish,glb,props)
fishery

pars <- c(7,0.3)
scalepar <- magnitude(pars)
bestL <- optim(pars,aspmLL,method="Nelder-Mead",dynfun=dynamicsF,infish=fish,
               inglb=glb,inprops=props,
               control=list(maxit = 1000, parscale = scalepar))
outfit(bestL)
fishery <- dynamicsF(bestL$par,fish,glb,props)
fishery


pars <- c(7.0,0.3,-7.7) # logR0, sigCE, depletion
scalepar <- magnitude(pars)
bestL <- optim(pars,aspmLL,method="Nelder-Mead",dynfun=dynamicsH,infish=fish,
               inglb=glb,inprops=props,
               control=list(maxit = 1000, parscale = scalepar))
outfit(bestL)
fishery <- dynamicsH(bestL$par,fish,glb,props)
print(round(fishery,4)) 


pars <- c(7.0,0.3,-7.7) # logR0, sigCE, depletion
scalepar <- magnitude(pars)
bestL <- optim(pars,aspmLL,method="Nelder-Mead",dynfun=dynamicsF,infish=fish,
               inglb=glb,inprops=props,
               control=list(maxit = 1000, parscale = scalepar))
outfit(bestL)
fishery <- dynamicsF(bestL$par,fish,glb,props)
print(round(fishery,4)) 

pars <- c(7.0,0.3,-7.7) # logR0, sigCE, depletion
bestL <- fitASPM(initpar=pars,minfun=aspmLL,dynfun=dynamicsF,infish=fish,
               inglb=glb,inprops=props)
outfit(bestL)
fishery <- dynamicsF(bestL$par,fish,glb,props)
print(round(fishery,4)) 





data("westroughy")
fish <- westroughy$fish
glb <- westroughy$glb
props <- westroughy$props
pars <- c(7,0.3,-7.7)
aspmLL(pars,dynamicsH,fish,glb,props)      # should be -2.277029
bestL <- fitASPM(pars,aspmLL,dynamicsH,infish=fish,inglb=glb,inprops=props)
bestL
fishery <- dynamicsH(bestL$par,fish,glb,props)
round(fishery,4)


# getProduction-------------------------------------

grsearch(f=matchC, c(0.2,0.45), M=0.036, cyr=5117.988, Byr=15760.748, tol = 1e-09)


data("westroughy")
fish <- westroughy$fish
glb <- westroughy$glb
props <- westroughy$props
pars <- c(7.0,0.3)
bestL <- nlminb(start=pars,aspmLL,dynfun=dynamicsH,infish=fish,inglb=glb,
                inprops = props,
              control=list(eval.max=500,iter.max=300,trace=0,rel.tol=1e-08))
prod <- getProduction(exp(bestL$par[1]),infish=fish,inglb=glb,inprops=props,
                      Hrg=c(0.0005,0.07,0.0005),nyr=100)
prod[78:100,]
prod[84:89,]
