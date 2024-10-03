

# reference Deriso et al (1985) and Quinn and Deriso (1999)

suppressPackageStartupMessages({
  require(codeutils)
  library(hplot)

})


prefixdir <- getDBdir()
# datadir <- paste0(prefixdir,"AbaloneData")
# rawdir <- paste0(prefixdir,"A_codeUse/rforcpue_use/")
options("show.signif.stars"=FALSE,"stringsAsFactors"=FALSE,
        "max.print"=50000,"width"=240)

# functions ------------------------

ricker <- function(alpha, beta, Sp) {
  N0 <- alpha * Sp * exp(-beta*Sp)
  return(N0)
} # end of Ricker



alpha <- 0.007643
beta <- 6.336e-10



Sp <- seq(1e7,1e10,1e8)

plotprep(width=9, height=5)
parset()
plot1(Sp/1e6,ricker(alpha,beta,Sp),defpar=FALSE,lwd=2)


data(halibut)
list2env(halibut,envir=environment())

# parameters = nyrs(age8) + (nages-1)(a9-a20) + nyrs(F) + 2(q) + 2 (alpha, beta) + 4(select)


# selectivity ---------------------
ls1 <- c(-1.032,-0.719,-0.494,-0.356,-0.238,-0.087,-0.054)
s1 <- exp(ls1)
ls2 <- c(-1.706,-1.129,-0.749,-0.523,-0.265,-0.205,-0.077)
s2 <- exp(ls2)


plotprep(width=9, height=5)
parset()
plot1(ages,c(s1,rep(1.0,6)),lwd=2,defpar=FALSE)
lines(ages,c(s2,rep(1.0,6)),lwd=2,col=2)


sel <- function(a,b,ages,maxage=0) {
  ans <- (ages^a) * exp(-b*ages)
  maxs <- max(ans)
  ans <- ans/maxs
  if ((maxage > 0) & (maxage %in% ages)) {
    ans[ages>=maxage] <- 1.0
  }
  return(ans)
} # end of sel



a <- 5.66869
b <- 0.3158
ages <- 8:20


sa <- sel(a=a,b=b,ages=ages,maxage=15)

plotprep(width=9, height=5)
parset()
plot1(ages,sa,defpar=FALSE,lwd=2)
lines(ages,c(s2,rep(1.0,6)),lwd=2,col=2)
lines(ages,c(s1,rep(1.0,6)),lwd=2,col=3)



ssq <- function(pars,...) {
  predsa <- sel(pars[1],pars[2],ages)
  ans <- sum((sa - predsa)^2)
  return(ans)
} 

out <- optim(par=c(a,b),ssq,method="Nelder-Mead",ages=c(8:14),sa=s2)

out

cbind(c(s2,rep(1.0,6)),sel(out$par[1],out$par[2],ages=ages,maxage=15))










