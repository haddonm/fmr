


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




ans2 <- microbenchmark(
  result <- grsearch(matchC, interval=c(0,1.25), M=M,cyr=catch[2],Byr=expB[2]),
  # out2 <- grsearch(matchC,interval=c(0,1.0),M=M,cyr=catch[2],Byr=expB[2],
  #                  tol = 1e-08),
  out3 <- datalowSA::getProductionC(inR0=exp(ans$par[1]),infish=fish,
                                    inglb=glb,inprops=props,
                                    Hrg=c(0.02,0.15,0.0025)),
  times=10
)
ans2














