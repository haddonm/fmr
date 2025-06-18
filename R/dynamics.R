

#' @title findFs uses an iterative strategy for finding each year's F values
#' 
#' @description findFs is used when searching for the instantaneous F that will
#'     produce the required catch in a given year in an age-structured model
#'     when there are multiple gears, with different selectivity being used.
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
#' @examples
#' print("wait on example data sets")
findFs <- function(cyr,Nyr,sel,aaw,M,reps=8) {
  # cyr=obsC;Nyr=Nt[,yr-1];sel=as.matrix(sel[,pickft]);aaw=aaw;M=M;reps=reps
  sel <- as.matrix(sel)
  nages <- length(sel[,1])
  nfleet <- ncol(sel)
  predCyr <- numeric(nfleet)
  wata <- aaw/1000
  Byr <- sum((Nyr*sel*wata))
  temp1 <- cyr / Byr
  join1 <- 1/(1 + exp(30*(temp1 - 0.95)))
  tempyr <- (join1 * temp1) + (0.95 * (1 - join1))
  fFyr <- -log(1 - tempyr)
  wtNyr <- wata * Nyr
  sF <- matrix(0,nrow=nages,ncol=nfleet)
  for (ft in 1:nfleet) sF[,ft] <- sel[,ft]*fFyr[ft]
  for (ft in 1:nfleet)
    predCyr[ft] <- sum((sF[,ft]/(sF[,ft] + M)) * wtNyr * 
                         (1 - exp(-(M + sF[,ft]))))
  for (i in 1:reps) {
    Zadj <- cyr/(predCyr)# + 0.0001)
    fFyr <- Zadj * fFyr
    for (ft in 1:nfleet) sF[,ft] <- sel[,ft]*fFyr[ft]
    for (ft in 1:nfleet)
      predCyr[ft] <- sum((sF[,ft]/(sF[,ft] + M)) * wtNyr * 
                           (1 - exp(-(M + sF[,ft]))))
  }
  # cat(Zadj,cyr[1]-predCyr[1],cyr[2]-predCyr[2],sum(cyr),sum(predCyr),"\n")
  return(fFyr) 
} # end of findFs
