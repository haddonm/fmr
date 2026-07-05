

#' @title findFs uses an iterative strategy for finding each year's F values
#' 
#' @description findFs is used when searching for the instantaneous F that will
#'     produce the required catch in a given year in an age-structured model
#'     when there are multiple gears, with different selectivity being used. It
#'     uses a method modified from that described in Methot and Wetzel 
#'     (2013, p7). This arrangement is designed solely for annual fisheries.
#'     By including a 'season' term and dividing each estimate of fFyr by
#'     the fraction of a year in each 'season', this adjustment routine could 
#'     be applied to seasonal fisheries.  
#'
#' @param cyr a vector of known catches in a given year for nfleets
#' @param Nyr the numbers at size at the start of the given year or end of the 
#'     year before
#' @param sel a matrix of the selectivity of the fishing gears
#' @param aaw the weight-at-age
#' @param M the instantaneous natural mortality rate
#' @param reps how many internal loops to use finding each F, default = 6
#'
#' @returns the fully selected fishing mortality rates
#' @export
#' 
#' @references Methot, R.D. and C.R. Wetzel (2013) Stock synthesis: A biological 
#'     and statistical framework for fish stock assessment and fishery 
#'     management. Supplementary Material: Appendix A: Technical Description of 
#'     the Stock Synthesis assessment program \emph{Fisherie Research} 142:
#'     86-99. http://dx.doi.org/10.1016/j.fishres.2012.10.012 
#'
#' @examples
#' \dontrun{
#'   require(fmr)
#'  library(codeutils)
#'  library(hplot)
#'  library(knitr)
#'  data("westroughy")
#'  fish <- westroughy$fish
#'  glb <- westroughy$glb
#'  props <- westroughy$props
#'  pars <- c(7,-0.4,-6.7)  
#'  bestFD <- fitASPM(initpar=pars,minfun=dynF,infish=fish,inglb=glb,
#'                    inprops=props,gradtol=1e-05,stepmax=0.1,steptol=1e-07,
#'                    hessian=TRUE,reps=9)
#'  outfit(bestFD,digits=7,title="Instantaneous Rates 2",
#'         parnames=c("LnR0","Ln(sigCE)","Ln(q)"))   
#' }                
findFs <- function(cyr,Nyr,sel,aaw,M,reps=8) {
  sel <- as.matrix(sel)
  nages <- length(sel[,1]) # uses selectivity by age
  nfleet <- ncol(sel)      # separate selectivity by fleet
  predCyr <- numeric(nfleet)
  wata <- aaw/1000
  Byr <- sum((Nyr*sel*wata)) # exploitable biomass by fleet in year yr
  temp1 <- cyr / (Byr + 0.1*cyr)  # next 4 lines = Pope's approximation
  join1 <- 1/(1 + exp(30*(temp1 - 0.95))) # join keeps it differentiable
  tempyr <- (join1 * temp1) + (0.95 * (1 - join1)) #approx mid-yr harvest rate
  fFyr <- -log(1 - tempyr)
  wtNyr <- wata * Nyr
  sF <- matrix(0,nrow=nages,ncol=nfleet)
  for (ft in 1:nfleet) {
    sF[,ft] <- sel[,ft]*fFyr[ft]
    predCyr[ft] <- sum((sF[,ft]/(sF[,ft] + M)) * wtNyr * 
                         (1 - exp(-(M + sF[,ft]))))
  }
  for (i in 1:reps) {
    Zadj <- cyr/(predCyr + 0.0000001)  # limits how precise one can get
    fFyr <- Zadj * fFyr
    for (ft in 1:nfleet) {
      sF[,ft] <- sel[,ft]*fFyr[ft]
      predCyr[ft] <- sum((sF[,ft]/(sF[,ft] + M)) * wtNyr * 
                           (1 - exp(-(M + sF[,ft]))))
    }
  }
  Zadj <- cyr/(predCyr + 0.0000001)  # limits how precise one can get
  fFyr <- Zadj * fFyr
  return(fFyr) 
} # end of findFs

#' @title plottmbprof a wrapper function for the tmbprofile plot
#' 
#' @description plottmbprof provides a simplified interface to the RTMB 
#'     tmbprofile plot. It plots the profile values and adds the optimum and
#'     CI valued confidence intervals.
#'
#' @param inmod the RTMB model from MakeADFun
#' @param parname the character name of the parameter from the model, the
#'     column name of the variable of parameter being profiled
#' @param CI default = 0.95 the probability level of the Confidence intervals
#' @param adjust by how much should the labels be adjusted down and up, 
#'     default = 0.05 = 5% of the y-scale
#'
#' @returns the CI invisibly
#' @export
#'
#' @examples
#' print("wait on examples")
#' # syntax:  plottmbprof(model,"logR0",CI=0.95,adjust=0.05)
#' # inmod=model; parname="logR0"; CI=0.95
plottmbprof <- function(inmod,parname,CI=0.95,adjust=0.05) { 
  prof <- tmbprofile(inmod,parname,trace=F)
  maxy <- getmax(prof$value,mult=1)
  miny <- which.min(prof$value)
  plot(prof,lwd=2,panel.first=grid(),ylim=c(min(prof$value),maxy),yaxs="i")
  optval <- prof[miny,parname]
  abline(v=optval,lwd=3,col=2)
  text(x=optval,y=(1-adjust)*maxy,round(optval,3),cex=1.0,pos=4)
  CI <- confint(prof,level=CI)
  ymin <- prof$value[miny]
  text(x=CI[1],y=(1+adjust)*ymin,round(CI[1],3),cex=1,pos=4)
  text(x=CI[2],y=(1+adjust)*ymin,round(CI[2],3),cex=1,pos=2)
  return(invisible(CI))
} # end of plottmbprof

