% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmr_utils.R
\name{getssq2}
\alias{getssq2}
\title{getssq2 calculates the sum of squares for the age-structured model}
\usage{
getssq2(pars, M, yrs, ages, onaa, owa, cpue)
}
\arguments{
\item{pars}{a vector of 28 parameters, 1:9 being the log(rec) for yrs 1929:
1937,10:17 being the log(numbers-at-ages) 10 - 3 for year 1929 (rather
than estimating recruitments for 8 years prior to 1929), 18-26 are log(F)
for years 29:37, and 27:28 are the selectivity parameters sel50, sel95.}

\item{M}{the natural mortality instantaneous rate}

\item{yrs}{a vector of the years of fisheries data, eg 1929:1937}

\item{ages}{a vector of the ages, eg 2:10}

\item{onaa}{the observed numbers-at-age in the catches}

\item{owa}{either the matrix of observed weight-at-age years as rows,
weight-at-age as columns. If a vector of average predicted weight-at-age
is being used it MUST first be converted into a matrix of length(yrs)
rows where the vector is duplicated in each row. Assuming 9 years and 10
ages one could use matrix(data=rep(owa,9),nrow=9,ncol=10,byrow=TRUE)}

\item{cpue}{the vector of cpue for each year}
}
\value{
a scalar holding the combined sum of squared residuals
}
\description{
getssq2 calculates the combined sum of squares residuals for
both the observed numbers-at-age and the CPUE. This uses an example from
Fournier and Archibald, 1982; which was also used as an example in
Haddon, 2011.
}
\examples{
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

getssq2(pars,M=0.2,yrs,ages,onaa,owa,cpue=fish[,"obsce"]) #should be 331.1994
}
