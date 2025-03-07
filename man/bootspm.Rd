% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spm.r
\name{bootspm}
\alias{bootspm}
\title{bootspm conducts a bootstrap analysis on a spm model}
\usage{
bootspm(optpar, fishery, iter = 100, schaefer = TRUE)
}
\arguments{
\item{optpar}{The optimum model parameters from an earlier analysis}

\item{fishery}{the fishery data containing the original observed cpue values}

\item{iter}{the number of bootstrap replicates to be run}

\item{schaefer}{default is TRUE, determines whether a Schaefer or a Fox model
is run}
}
\value{
a list of two matrices. One containing the bootstrap parameters and
the other containing some of the dynamics, including the ModelB, the
bootstrap CPUE sample, the Depletion, and the annual harvest rate.
}
\description{
bootspm conducts a bootstrap analysis on a spm model. It does
this by saving the original fishery data, estimating the cpue residuals,
and multiplying the optimum predicted CPUE by a bootstrap sample of the
log-normal residuals (Haddon, 2011, p311). This bootstrap sample of CPUE
replaces the original fish(,"cpue") (should be square brackets) and
the model is re-fitted. This is repeated iter times and the outputs
reported ready for the derivation of
percentile confidence intervals. The optimum solution is used as the
first bootstrap replicate (it is standard practice to include the
original fit in the bootstrap analysis). If 1000 replicates are run this
procedure can take a couple of minutes on a reasonably fast computer. A
comparison of the mean with the median should provide some notion of any
bias in the mean estimate.
}
\examples{
\dontrun{
 data(dataspm)
 pars <- c(r=0.2,K=6000,Binit=2800)
 ans <- fitSPM(pars,dataspm$fish,schaefer=TRUE,maxiter=1000)
 boots <- bootspm(ans$par,fishery=dataspm$fish,iter=500,schaefer=TRUE)
 dynam <- boots$dynam
 bootpar <- boots$bootpar
 MSY <- bootpar[,"r"]*bootpar[,"K"]/4
 Depl <- dynam[,length(fish[,"year"]),"Depletion"] # pick the last year
 bootpar <- cbind(bootpar,MSY,Depl)
 rows <- colnames(bootpar)
 columns <- c(c(0.025,0.05,0.5,0.95,0.975),"Mean")
 bootCI <- matrix(NA,nrow=length(rows),ncol=length(columns),
                dimnames=list(rows,columns))
 for (i in 1:length(rows)) { 
   tmp <- sort(bootpar[,i])
   qtil <-  quantile(tmp,probs=c(0.025,0.05,0.5,0.95,0.975),na.rm=TRUE)
   bootCI[i,] <- c(qtil,mean(tmp,na.rm=TRUE))
 }
 round(bootCI,3)
}
}
