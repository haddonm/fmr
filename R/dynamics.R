

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
#' print("wait on example data sets")
#' # cyr=obsC;Nyr=Nt[,yr-1];sel=as.matrix(sel[,pickft]);aaw=aaw;M=M;reps=reps
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
  fFyr <- -log(1 - tempyr)/season
  wtNyr <- wata * Nyr
  sF <- matrix(0,nrow=nages,ncol=nfleet)
  for (ft in 1:nfleet) {
    sF[,ft] <- sel[,ft]*fFyr[ft]
    predCyr[ft] <- sum((sF[,ft]/(sF[,ft] + M)) * wtNyr * 
                         (1 - exp(-(M + sF[,ft]))))
  }
  for (i in 1:reps) {
    Zadj <- cyr/(predCyr + 0.000001)  # limits how precise one can get
    fFyr <- Zadj * fFyr
    for (ft in 1:nfleet) {
      sF[,ft] <- sel[,ft]*fFyr[ft]
      predCyr[ft] <- sum((sF[,ft]/(sF[,ft] + M)) * wtNyr * 
                           (1 - exp(-(M + sF[,ft]))))
    }
  }
  return(fFyr) 
} # end of findFs
