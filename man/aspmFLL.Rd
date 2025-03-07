% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aspm.r
\name{aspmFLL}
\alias{aspmFLL}
\title{aspmFLL negative log-likelihood for the ASPM}
\usage{
aspmFLL(par, infish, inglb, inprops)
}
\arguments{
\item{par}{the dynamics relies on many parameters sitting in the global
environment in particular ages, nages, maxage, M, maa, waa, sela, fish,
and nyrs. 'pars' can contain either two or three parameters. 1) is
the log-transformed average unfished recruitment, inR0. 2) is the
variability around the index of relative abundance (cpue) during the
fitting process, and if is present 3) is the initial depletion level
initdepl, which if present will be fitted as well.}

\item{infish}{the fish data.frame from readdata or built in dataset}

\item{inglb}{the glb data.frame from readdata or built in dataset}

\item{inprops}{the props data.frame from readdata or built in dataset}
}
\value{
a scalar that is the negative log-likelihood using Log-Normal
random errors  designed for use with an ASPM using instantaneous
mortality rates rather than harvest rates for the fishing mortality.
}
\description{
aspmFLL is the negative log-likelihood for Log-Normally
distributed data, set-up for use with the ASPM model fitting.
This uses the aspmdynamicsF function to calculate the dynamics, including
the predicted index valuess for calculating the negative log-likelihood.
}
\examples{
\dontrun{
data(westroughy)
fish <- westroughy$fish
glb <- westroughy$glb
props <- westroughy$props
pars <- c(14,0.3)
dynamics(pars,fish,glb,props)
aspmFLL(pars,fish,glb,props)      # should be 5.171
pars <- c(14.0,0.3,0.95) # logR0, sigCE, depletion
aspmdynamicsF(pars,fish,glb,props)    # note the harvest rates of 85\% exploitable biomass
aspmFLL(pars,fish,glb,props)      # should be 114.9547
pars <- c(14,0.3)
bestL <- optim(pars,aspmLL,method="Nelder-Mead",infish=fish,inglb=glb,inprops=props,
               control=list(maxit = 1000, parscale = c(10,0.1)))
outoptim(bestL)
fishery <- aspmdynamicsF(bestL$par,fish,glb,props)
print(round(fishery,4)) 
}
}
