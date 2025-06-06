% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aspm.r
\name{matchC}
\alias{matchC}
\title{matchC gives difference between predicted catch and observed catch}
\usage{
matchC(f, M, cyr, Nyr, sel, waa)
}
\arguments{
\item{f}{the trial predicted instantaneous fishing mortality}

\item{M}{the instantaneous natural mortality value}

\item{cyr}{the catch in yr t}

\item{Nyr}{the numbers at age at the end of yr-1}

\item{sel}{the selectivity of the fishing gear in the model}

\item{waa}{the weight -at- age for the species}
}
\value{
the absolute difference between the predicted and observed catches
}
\description{
matchC for a given instantaneous fishing mortality rate
calculates the difference between the predicted and observed catches.
matchC uses the Baranov catch equation to account, when usung the
instantaneous rates, for the fish that died due to fishing that would
have died anyway, due to natural mortality. The equation is:
Ct = (F/(F + M)) x (1 - exp(-(M + F))) x Bexp. Of course, one needs to
account for selectivity aby age and F by year, so check the equation
in the mfmur book. matchC is used by the optimize function which uses
'a combination of a golden section search and successive parabolic
interpolation' to minimize the difference between the absolute
difference between the catch and the predicted catch.
}
\examples{
data(westroughy)
fish <- westroughy$fish
glb <- westroughy$glb
props <- westroughy$props
sel <- props$sela
pars <- c(7.15,-1,-7.7) # logR0, sigCE, estimate avq
out <- dynamicsF(pars,fish,glb,props,full=TRUE)  
matchC(f=0.3,M=0.05,cyr=15340,Nyr=out$Nt[,2],sel=sel,props$waa)
matchC(f=0.35,M=0.05,cyr=15340,Nyr=out$Nt[,2],sel=sel,props$waa)
out2 <- optimize(matchC,interval=c(0,1),M=0.05,cyr=3924,Nyr=out$Nt[,2],
                sel=sel,props$waa)
out2
f <- out2$minimum
(397388.66 * (1 - exp(-(0.05 + f))) * f/(0.05 + f))
# f=0.2679;M=0.0036;cyr=3924;Nyr=Nt[,yr-1];sel=sel;waa=props$waa
}
