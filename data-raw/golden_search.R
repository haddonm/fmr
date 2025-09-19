


library(Rcpp)
library(microbenchmark)
library(fmr)
library(codeutils)
library(hplot)

cppFunction('double matchCC(double Fy, double M, double cyr, double Byr) {
  return abs((Byr * (1 - exp(-(M + Fy))) * Fy/(M + Fy)) - cyr);
}')

matchC <- function(Fy,M,cyr,Byr) {
  out <- abs((Byr * (1 - exp(-(M + Fy))) * Fy/(M + Fy)) - cyr)
  return(out)
}




data("westroughy")
fish <- westroughy$fish
glb <- westroughy$glb
props <- westroughy$props
pars <- c(7.1,-1,-7.7) # logR0, sigCE, estimate avq
fishery <- dynF(pars,fish,glb,props)
pars <- c(7.1086,-1.148,-7.838)
bestL <- optim(pars,dynF,method="Nelder-Mead",infish=fish,inglb=glb,
               inprops=props,control=list(maxit=1000,parscale = c(10,1,10)))
str(bestL)
outfit(bestL,digits=6)
out <- dynF(bestL$par,fish,glb,props,full=TRUE)
print(round(out$fishery,4))
fishery <- out$fishery
printV(round(fishery[,"catch"] - fishery[,"predC"],3))


expB <- out$fishery$exploitB
M <- glb$M
catch <- out$fishery$catch

ans2 <- microbenchmark(
  result <- grsearch(matchC, interval=c(0,1.25), M=M,cyr=catch[2],Byr=expB[2]),
  result2 <- grsearch(matchCC,interval=c(0,1.0),M=M,cyr=catch[2],Byr=expB[2]),
  #                  tol = 1e-08),
  # out3 <- datalowSA::getProductionC(inR0=exp(bestL$par[1]),infish=fish,
  #                                   inglb=glb,inprops=props,
  #                                   Hrg=c(0.02,0.15,0.0025)),
  times=10
)
ans2














