

#' @title addnorm - adds a normal distribution to a histogram of a data set.
#'
#' @description  addnorm - adds a normal distribution to a histogram of a data
#'    set. This is generally to be used to illustrate whether log-transformation
#'    normalizes a set of catch or cpue data.
#' @param inhist - is the output from a call to 'hist' (see examples)
#' @param xdata -  is the data that is being plotted in the histogram.
#' @param inc - defaults to a value of 0.01; is the fine grain increment used to
#'    define the normal curve. The histogram will be coarse grained relative to
#'    this.
#' @return a list with a vector of 'x' values and a vector of 'y' values (to be
#'    used to plot the fitted normal probability density function), and a vector
#'    used two called 'stats' containing the mean and sandard deviation of the
#'    input data
#' @export addnorm
#' @examples
#' x <- rnorm(1000,mean=5,sd=1)
#' dev.new(height=6,width=4,noRStudioGD = TRUE)
#' par(mfrow= c(1,1),mai=c(0.5,0.5,0.3,0.05))
#' par(cex=0.85, mgp=c(1.5,0.35,0), font.axis=7)
#' outH <- hist(x,breaks=25,col=3,main="")
#' nline <- addnorm(outH,x)
#' lines(nline$x,nline$y,lwd=3,col=2)
#' print(nline$stats)
addnorm <- function(inhist,xdata,inc=0.01) {
  lower <- inhist$breaks[1]
  upper <- tail(inhist$breaks,1)
  cw <- inhist$breaks[2]-inhist$breaks[1]
  x <- seq(lower,upper, inc) #+ (cw/2)
  avCE <- mean(xdata,na.rm=TRUE)
  sdCE <- sd(xdata,na.rm=TRUE)
  N <- length(xdata)
  ans <- list(x=x,y=(N*cw)*dnorm(x,avCE,sdCE),stats=c(avCE,sdCE,N))
  return(ans)
} # end of addnorm

#' @title calcrmse calcuates the root mean square error for two vectors
#' 
#' @description given a set of observations and predicted values, calcrmse
#'     calculates the root mean square error sqrt(resid^2/n) to provide a 
#'     measure of relative fit. An assumption of normal errors is made so
#'     if the observations and predicted values actually have a log-normal
#'     distribution one should log-transform the input values.
#'
#' @param obs the observed data
#' @param pred the predicted values whose fit to the observations is to 
#'     measured
#'
#' @return a scalar value which is the rmse.
#' @export
#'
#' @examples
#' \dontrun{
#' x <- rep(1,10)
#' y <- c(1.109,1.210,0.947,0.933,0.832,0.864,0.633,0.820,1.004,1.049)
#' calcrmse(log(x),log(y))  # should be 0.1899495
#' }
calcrmse <- function(obs,pred){ # obs=fishery[,"CPUE"]; pred=fishery[,"PredCE"]
  pickCE <- which(!is.na(obs))
  resid <- obs[pickCE] - pred[pickCE]
  rmse <- sqrt(sum(resid^2)/length(resid))
  return(rmse)
}

#' @title getLNCI gets the log-normal confidence intervals
#' 
#' @description getLNCI takes the mean and the standard deviation and produces
#'     the asymmetric log-normal confidence intervals around the mean values 
#'
#' @param av the mean value or a vector of mean values
#' @param se the standard deviation 
#' @param P the percent used for the CI, defaults to 95.
#'
#' @return a vector of three for a single input or a matrix of 3 columns for 
#'     input vectors
#' @export
#'
#' @examples
#' \dontrun{
#'   av <- c(4.0,2.15)
#'   se <- 0.33
#'   getLNCI(av,se,P=95)
#'   se <- c(0.33,0.4)
#'   getLNCI(av,se)
#' }
getLNCI <- function(av,se,P=95){  # av=fissp[,"CPUE"]; se=rmse;P=0.95
  Zmult <- -qnorm((1-(P/100))/2.0)
  lower <- av * exp(-Zmult*se)
  upper <- av * exp(Zmult*se)
  if (length(av) > 1) {
    result <- cbind(lower,av,upper)
    colnames(result) <- c("lower","mean","upper")
  } else {
    result <- c(lower=lower,mean=av,upper=upper)
  }
  return(result)
} # end of getLNCI   


#' @title incol is a utility to determine is a column is present in a matrix
#'
#' @description incol is a utility to determine whether a names columns is
#'     present in a given matrix or data.frame.
#'
#' @param incol the name of the column; defaults to "year" as an example
#' @param inmat the matrix or data.frame within which to search for incol
#'
#' @return TRUE or FALSE
#' @export
#'
#' @examples
#' \dontrun{
#' test <- matrix(c(1,2,3,4),nrow=2,ncol=2,dimnames=list(1:2,c("year","Catch")))
#' print(test)
#' iscol("year",test)
#' iscol("Catch",test)
#' iscol("catch",test)
#' iscol("ages",test)
#' }
iscol <- function(incol="year",inmat) { # incol="ages"; inmat=dat
  if (length(grep(incol,colnames(inmat))) < 1) return(FALSE)
  else return(TRUE)
}


#' @title penalty0 enables the adding of a large penalty as one approaches 0.0
#'
#' @description penalty0 allows for the option of adding a large penalty as
#'     a parameter approaches 0.0 . See negLL1 for example code that 
#'     contains such a parameter. For example, when
#'     fitting an spm sometimes the optimal mathematical model fit can occur
#'     by depressing the r value to 0 or even go negative. Input values 
#'     < 0.006 begin to generate large values as one goes smaller. The
#'     examples below illustrate this.
#'
#' @param x the parameter value that potentially incurs a penalty
#'
#' @return a single value as a penalty to be added to a Log-Likelihood or SSQ
#' @export
#'
#' @examples
#'   penalty0(0.5)
#'   penalty0(0.1)
#'   penalty0(0.01)
#'   penalty0(0.005)
penalty0 <- function(x){
  ans <- 100*exp(-1000*x)
  return(ans)
} # end of penalty0

#' @title penalty1 adds an increasingly large penalty as a value approaches 1.0
#'
#' @description penalty1 allows for the option of adding a large penalty as
#'     a parameter approaches 1.0 and moves to become larger than 1. For 
#'     example, when fitting a surplus production model sometimes the 
#'     optimal mathematical model fit can occur by implying catches greater
#'     than available biomass, implying harvest rates > 1.0. By adding a 
#'     large penalty to such values and adding those to the likelihood 
#'     such strange outcomes can be avoided. This will accept a single
#'     value or a vector.
#'
#' @param x the parameter value that potentially incurs a penalty
#'
#' @return a single value as a penalty to be added to a Log-Likelihood or SSQ
#' @export
#'
#' @examples
#'  x <- c(0.5,0.8,0.88,0.9,0.98,0.99,1.01)
#'  round(cbind(x,penalty1(x)),4)
penalty1 <- function(x){
  nl <- length(x)
  ans <- numeric(nl)
  pick <- which(x > 0.5)
  ans[pick] <- 100.0*(abs((1-abs(x[pick])-0.5))/0.5)^50
  return(ans)
} # end of penalty1