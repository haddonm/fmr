% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aspm.r
\name{plotceASPM}
\alias{plotceASPM}
\title{plotceASPM plots just the fit of the ASPM model to the CPUE data}
\usage{
plotceASPM(infish, CI = NA, defineplot = TRUE)
}
\arguments{
\item{infish}{output from the dynamics function using the optimum parameters}

\item{CI}{the output matrix from getLNCI function. Needs to have the lower CI
in column 1 and the upper CI in the third column.}

\item{defineplot}{defaults to TRUE, determines whether to set up an new
graphics window.}
}
\value{
returns nothing but does plot a graph
}
\description{
plotceASPM plots just the fit of the ASPM model to the CPUE data
and provides a more detailed visual than plotASPM. It has the option of
including lognormal confidence intervals around the predcted CPUE.
}
\examples{
data(westroughy)
fish <- westroughy$fish
glb <- westroughy$glb
props <- westroughy$props
pars <- c(7.1,-1,-7.7)
bestaspm <- fitASPM(pars,minfun=dynamicsH,infish=fish,
                    inglb=glb,inprops=props)
out <- dynamicsH(bestaspm$estimate,fish,glb,props,full=TRUE)
ceCI <- getLNCI(out$fishery[,"predCE"],exp(bestaspm$estimate[2]))
plotceASPM(out$fishery,CI=ceCI)
}
