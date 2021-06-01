



# test calcnaaC ----------------------------------------------------------------

library(fmr)
library(rutilsMH)
library(Rcpp)


data(ocaa)
data(fish)
data(owaa)
data(param)
pars <- param[1:28,2]

yrs <- as.numeric(fish[,"year"])
ages <- as.numeric(2:10)
onaa <- ocaa[,ages]
owa <- as.matrix(owaa[,ages])
M <- 0.2
out <- calcnaaC(pars,M,yrs,ages,owa)
getssq2(pars,M,yrs,ages,onaa,owa,fish[,"obsce"])


calcnaa2 <- function(pars,M,yrs,ages,owa) {  
  # pars=pars; pnaa=pnaa; M=M; sel=selec ; age=ages 
  sel <- logistic(pars[27],pars[28],ages)
  nyr <- length(yrs)
  nage <- length(ages)
  pnaa <- matrix(0,nrow=nyrs,ncol=nage,dimnames=list(yrs,ages)) 
  enaa <- pnaa
  exB <- numeric(nyr)
  endR <- nyr
  endN <- endR + nage - 1
  f1 <- endN - 1
  pnaa[,1] <- exp(pars[1:endR])
  pnaa[1,2:nage] <- exp(pars[(endR+1):endN])
  for (yr in 2:nyr) { 
    for (age in 2:nage) {  # yr=2; age=2
      sF <- sel[age - 1] * exp(pars[yr + f1])
      pnaa[yr,age] <- pnaa[(yr-1),(age-1)] * exp(-(M + sF))
    }
  }
  for (yr in 1:nyr) {
    enaa[yr,] <- pnaa[yr,] * sel
    exB[yr] <- sum(enaa[yr,] * owa[yr,])/1000 
  }
  return(list(pnaa=pnaa,enaa=enaa,exB=exB,sel=sel))                     
} # end of calcnaa2


out <- calcnaa2(pars,M,yrs,ages,owa)
out2 <- calcnaaC(pars,M,yrs,ages,owa)

round(out$pnaa,5)
round(out2$pnaa,5)


library(microbenchmark)

microbenchmark(
  out <- calcnaa2(pars,M,yrs,ages,owa),
  out2 <- calcnaaC(pars,M,yrs,ages,owa)
)

# end of test calcnaaC ---------------------------------------------------------
