

#' @title calccaa estimates the catch-at-age from teh predicted numbers-at-age
#' 
#' @description calccaa uses the Baranov catch equation to calculate the 
#'     predicted catch-at-age from teh predicted numbers-at-age (such as 
#'     calculates by calcnaa)
#'
#' @param pars the model parameter vector
#' @param pnaa the predicted numbers-at-age from the parameters
#' @param M the instantaneous natural mortality rate
#' @param sel the selectivity-at-age, which can be parameters
#' @param ages a vector of ages 
#'
#' @return a matrix of predicted catch-at-age with yrs as rows and ages as cols
#' @export
#'
#' @examples
#' print("wait on data-sets")
calccaa <- function(pars,pnaa,M,sel,ages) {  
  # pars=pars; pnaa=pnaa; M=M; sel=selec ; age=ages 
  caa <- pnaa
  nyr <- nrow(caa)
  nage <- length(ages)
  f1 <- nyr + nage - 1
  for (yr in 1:nyr) {
    for (age in 1:nage) { # yr=1; age=1
      sF <- sel[age] * exp(pars[yr + f1])
      caa[yr,age] <- (sF/(M + sF)) * pnaa[yr,age] * (1 - exp(-(M + sF)))
    }
  }
  return(caa=caa)                     
} # end of calccaa


#' @title calcnaa estimates the predicted numbers-at-age from the parameters
#' 
#' @description calcnaa estimates the predicted numbers-at-age from the 
#'     parameters, along with the estimate of M, and selectivity-at-age
#'
#' @param pars the model parameter vector
#' @param M the instantaneous natural mortality rate
#' @param sel the selectivity-at-age, which can be parameters
#' @param yrs a vector of yrs, used to label the rows
#' @param ages a vector of ages, used to label the columns 
#'
#' @return  a matrix of predicted numbers-at-age
#' @export
#'
#' @examples
#' print("wait on data-sets")
calcnaa <- function(pars,M,sel,yrs,ages) {  
  # pars=pars; pnaa=pnaa; M=M; sel=selec ; age=ages 
  nyr <- length(yrs)
  nage <- length(ages)
  pnaa <- matrix(0,nrow=nyr,ncol=nage,dimnames=list(yrs,ages))  
  endR <- nyr
  endN <- endR + nage - 1
  f1 <- endN - 1
  pnaa[,1] <- exp(pars[1:endR])
  pnaa[1,2:nage] <- exp(pars[(endR+1):endN])
  for (yr in 2:nyr) { 
    for (age in 2:nage) {  # yr=2; age=2
      sF <- sel[age - 1] * exp(pars[yr + f1])
      pnaa[yr,age] <- pnaa[(yr-1),(age-1)] * exp(-(M + sF))
    }
  }
  return(pnaa=pnaa)                     
} # end of calcnaa

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

#' @title getq calculates catchability from cpue and predicted exploitable biomass
#' 
#' @description getq uses a closed form version of the catchability that uses 
#'     the observed cpue and the predicted exploitable biomass to estimate the
#'     scaling factor between erxploitable biomass and cpue.
#'
#' @param cpue the vector of observed cpue
#' @param exB the output vector of predicted exploitable biomass
#'
#' @return the catchability as a scaler
#' @export
#'
#' @examples
#' cpue <- c(0.369,0.254,0.407,0.357,0.326,0.284,0.403,0.297,0.397)
#' exB <- c(3712,3251,3810,4499,4500,3977,3972,3220,3386) 
#' getq(cpue,exB)  # should be 8.964018e-05 
getq <- function(cpue,exB) {
  if (length(cpue) != length(exB))
    stop("Input vectors of different lengths in function: getq \n")
  nyr <- length(cpue)
  q <- exp(sum(log(cpue/exB))/nyr)
  return(q)
} # end of getq


#' @title getsingle extracts a single number from an input line of characters
#'
#' @description getsingle splits up a text line and translates the first non-
#'     empty character string into a number.
#'
#' @param inline the line of text, usually taken after using readLines
#' @param sep the separator used to divide the numbers from descriptive text.
#'     defaults to a comma.
#'
#' @return a single number
#' @export
#'
#' @examples
#' \dontrun{
#' x <- "12.3 , this is a number"
#' y <- "21.3 # 22.3 # here are two numbers"
#' getsingle(x)
#' getsingle(y,sep="#")
#' }
getsingle <- function(inline,sep=",") {  # inline=dat[41]
  tmp <- unlist(strsplit(inline,sep))
  tmp <- gsub(" ","",tmp)
  tmp <- tmp[nchar(tmp) > 0]
  return(as.numeric(tmp[1]))
}

#' @title getssq2 calculates the sum of squares for the age-structured model
#' 
#' @description getssq2 calculates the combined sum of squares residuals for
#'     both the observed numbers-at-age and the CPUE. This uses an example from
#'     Fournier and Archibald, 1982; which was also used as an example in 
#'     Haddon, 2011.
#'
#' @param pars a vector of 28 parameters, 1:9 being the log(rec) for yrs 1929:
#'     1937,10:17 being the log(numbers-at-ages) 10 - 3 for year 1929 (rather
#'     than estimating recruitments for 8 years prior to 1929), 18-26 are log(F)
#'     for years 29:37, and 27:28 are the selectivity parameters sel50, sel95.
#' @param M the natural mortality instantaneous rate
#' @param yrs a vector of the years of fisheries data, eg 1929:1937
#' @param ages a vector of the ages, eg 2:10
#' @param onaa the observed numbers-at-age in the catches
#' @param owa either the matrix of observed weight-at-age years as rows, 
#'     weight-at-age as columns. If a vector of average predicted weight-at-age 
#'     is being used it MUST first be converted into a matrix of length(yrs)
#'     rows where the vector is duplicated in each row. Assuming 9 years and 10 
#'     ages one could use matrix(data=rep(owa,9),nrow=9,ncol=10,byrow=TRUE)
#' @param cpue the vector of cpue for each year
#'
#' @return a scalar holding the combined sum of squared residuals
#' @export
#'
#' @examples
#' data(ocaa)
#' data(fish)
#' data(owaa)
#' data(param)
#' pars <- param[1:28,2]
#' yrs <- as.numeric(fish[,"year"])
#' ages <- as.numeric(2:10)
#' onaa <- ocaa[,ages]
#' owa <- as.matrix(owaa[,ages])
#' nage <- length(ages)
#' getssq2(pars,M=0.2,yrs,ages,onaa,owa,cpue=fish[,"obsce"]) #should be 332.5389
getssq2 <- function(pars,M,yrs,ages,onaa,owa,cpue) {
  out <- calcnaa(pars,M,yrs,ages,owa)   #calcnaaC(pars,M,yrs,ages,owa)
  pcaa <- calccaa(pars,out$pnaa,M,out$sel,ages)  #calccaaC(pars,out$pnaa,M,out$sel,ages)
  ssq1 <- sum((log(onaa/pcaa)^2),na.rm=TRUE)
  exB <- out$exB
  q <- exp(sum(log(cpue/exB))/length(yrs))
  ssq2 <- sum((log(cpue) - log(q * exB))^2)
  ssq= ssq1 + ssq2
  return(ssq)
} # end of getssq2


#' @title getvector  extracts a vector of numbers from a line of characters
#' 
#' @description getvector when reading in a csv file using readLines, 
#'     getvector extarcts a line of numbers from a specified line within
#'     the readLine object.This function works out how many numbers there
#'     are. If you wish to add a comment at the end of a vector of numbers
#'     it must be separated from tehm by the separator. e.g. a comma
#' @param indat the readLines object
#' @param locate the line number from which to extract the numbers
#' @param sep the separator between numbers, defaults to ","
#'
#' @return a vector of numbers
#' @export
#'
#' @examples
#' \dontrun{
#' x <- "12.3, 15.1, 8.7,10.3,  # this is a vector of numbers"
#' y <- "21.3 # 22.3 # 8.7 # 10.3 # here are four numbers"
#' getvector(x)
#' getvector(y,sep="#")
#' }
getvector <- function(indat,locate,sep=",") { # indat=dat; locate=pick+2;sep=","
  vect <- indat[locate]
  if (length(grep("\t",vect) > 0)) vect <- gsub("\t",sep,vect)
  vect <- unlist(strsplit(vect,sep))
  vect <- gsub(" ","",vect)
  vect1 <- vect
  vect <- suppressWarnings(as.numeric(vect))
  vect <- vect[nchar(vect) > 0]   
  if (!any(vect1 == "NA")) {
    vect <- vect[!is.na(vect)]
  }
  return(vect)
}

#' @title grsearch a golden ratio search algorithm to estimate F in an ASPM
#' 
#' @description grsearch implements a one-dimensional golden-ratio search for
#'     the instantaneous F that gives rise to a required catch level in a
#'     given year. It uses an input function and interval to search and a value 
#'     to match. The algorithm itself is described in many places on the www. 
#'
#' @param f the function whose value needs to match cyr (catch in year)
#' @param interval the interval over which to search
#' @param M the instantaneous natural mortality rate
#' @param cyr the observed catch in a given year
#' @param Byr the exploitable biomass at the start of the year in which the
#'     catches are to be taken
#' @param tol the tolerance of the match between the predicted and observed 
#'     catch
#'
#' @return an instantaneous fishing mortality rate that should generate the 
#'     observed catches from a given exploitable biomass in a year
#' @export
#'
#' @examples
#' matchC <- function(Fy,M,cyr,Byr) {
#'   return(abs((Byr * (1 - exp(-(M + Fy))) * Fy/(M + Fy)) - cyr))
#' }
#' grsearch(f=matchC,interval=c(0.2,0.45),M=0.036,cyr=5117.988,Byr=15760.748,
#'          tol = 1e-09)  # should be 0.4007973
grsearch <- function(f, interval, M, cyr, Byr, tol = 1e-09) {#.Machine$double.eps^0.25) {
  golden_ratio <- (sqrt(5) - 1) / 2  # Define the golden ratio
  lower=interval[1]
  upper=interval[2]
  x1 <- upper - golden_ratio * (upper - lower)  # Calculate initial points
  x2 <- lower + golden_ratio * (upper - lower)
  f1 <- f(x1,M, cyr, Byr)   # Evaluate function at initial points
  f2 <- f(x2,M, cyr, Byr)
  while (abs(upper - lower) > tol) {   # Iteratively narrow the interval
    if (f1 < f2) {
      upper <- x2
      x2 <- x1
      f2 <- f1
      x1 <- upper - golden_ratio * (upper - lower)
      f1 <- f(x1,M, cyr, Byr)
    } else {
      lower <- x1
      x1 <- x2
      f1 <- f2
      x2 <- lower + golden_ratio * (upper - lower)
      f2 <- f(x2,M, cyr, Byr)
    }
  }
  # Return the midpoint of the final interval as the estimated minimum
  return((lower + upper) / 2)
} # end of grsearch



#' @title incol is a utility to determine is a column is present in a matrix
#'
#' @description incol is a utility to determine whether a names columns is
#'     present in a given matrix or data.frame. If the input is neither a 
#'     matrix or a data.frame an error and warning will be thrown.
#'
#' @param incol the name of the column being looked for.
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
iscol <- function(incol,inmat) { # incol="ages"; inmat=pnaa
  if (class(inmat)[1] %in% c("matrix","data.frame")) { 
    if (length(grep(incol,colnames(inmat))) < 1) return(FALSE)
    else return(TRUE) 
  } else {
    stop(cat("input is neither a matrix nor a data.frame \n"))
  }
} # end of iscol

#' @title logistic standard selectivity function
#'
#' @description logistic calculates a Logistic curve that can be used as a
#'     selectivity function, or maturity curve, of wherever a logistic is
#'     required. This version uses the logistic function
#'     1/(1+exp(-log(19.0)*(lens-inL50)/(inL95-inL50))),
#'     which explicitly defines the SM50 and uses SM95 as the second parameter.
#' @param inl50 is the length/age at 50 percent selection/maturity/whatever
#' @param inl95 is the length/age at 5 percent selection/maturity/whatever
#' @param depend a vector of lengths/ages for which the logistic value will be
#'     calculated.
#' @return A vector of length/age(depend) containing predicted logistic values
#' @export
#' 
#' @examples
#' in50 <- 3.398
#' in95 <- 4.386
#' ages <- seq(2,10,1)
#' select <- logistic(inl50=in50,inl95=in95,depend=age)
#' round(cbind(ages,select),5)
logistic <- function(inl50,inl95,depend) {
  L50 <- exp(inl50)
  L95 <- exp(inl95)
  ans <- 1/(1+exp(-log(19.0)*(depend-L50)/(L95 - L50)))
  return(ans)
} # end of logistic

#' @title matchC gives difference between predicted catch and observed catch
#' 
#' @description matchC for a given instantaneous fishing mortality rate
#'     calculates the difference between the predicted and observed catches. 
#'     This is to be used by the optimize function which uses 'a combination of 
#'     a golden section search and successive parabolic interpolation' to 
#'     minimize the difference between the absolute difference between the 
#'     catch and the predicted catch.
#'
#' @param f the trial predicted instantaneous fishing mortality
#' @param M the instantaneous natural mortality value
#' @param cyr the catch in yr t
#' @param Byr the exploitable biomass at the start of yr t or end of t-1
#'
#' @return the absolute difference between the predicted and observed catches
#' @export
#'
#' @examples
#' matchC(f=0.04,M=0.05,cyr=15340,Byr=397388.66)
#' matchC(f=0.04037,M=0.05,cyr=15340,Byr=397388.66)
#' out <- optimize(matchC,interval=c(0,1),M=0.05,cyr=15340,Byr=397388.66)
#' out
#' f <- out$minimum
#' (397388.66 * (1 - exp(-(0.05 + f))) * f/(0.05 + f))
matchC <- function(f,M,cyr,Byr) {
  out <- abs((Byr * (1 - exp(-(M + f))) * f/(M + f)) - cyr)
  return(out)
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

#' @title readdata reads in a standard format data file
#'
#' @description readdata reads in a standard format data file. An example of the
#'     standard format is generated by using the function dataTemplate, which
#'     generates an example datafile which can be used to demonstrate the
#'     methods present in datalowSA, or can be used as a template to edit and
#'     input a new or different dataset.
#'
#' @param filename the filename (including the full path if required) containing
#'     the data in the standard format.
#' @param property does the data include a table of length-at-age, maturity-at-
#'     age, weight-at-age, fecundity-at-age, and that kind of thing (see the
#'     standard format for a list of data.)
#' @param verbose Default = TRUE, which prints out details as data is read in.
#'
#' @return a list of three objects, fish includes the Year, Catch, CPUE, and SE
#'     of the CPUE, glb containing an array of biological properties that can
#'     be global, and props containing age, laa, waa, maa, sela.
#' @export
#'
#' @examples
#' \dontrun{
#' dataTemplate(filename="test.csv", title="A test of the functions")
#' ans <- readdata("test.csv",property=FALSE)
#' str(ans)
#' }             # filename="sardine.csv"; 
readdata <- function(filename,property=FALSE,verbose=TRUE) {  # property=FALSE; filename=paste0(datadir,"eastdeepwatershark.csv")
  dat <- readLines(filename)
  spsname <- scan(file=filename,skip=0,what=character(),nlines = 1,quiet=TRUE)
  spsname <- removeEmpty(gsub(",","",spsname))
  pick <- grep("RESILIENCE",dat)[1]
  resilience <- scan(file=filename,skip=pick,what=character(),nlines = 1,quiet=TRUE)
  resilience <- removeEmpty(gsub(",","",resilience))
  pick <- grep("NYRS",dat)[1] + 1
  nyrs <- getsingle(dat[pick])
  if (verbose) cat(spsname,"\n\n resilience = ",resilience,
                   "  Num Years = ",nyrs,"\n\n")
  # get the fishery data into as many columns as required
  pick <- grep("YEARS",dat)[1]
  columns <- tolower(removeEmpty(unlist(strsplit(dat[pick],","))))
  numcols <- length(columns)
  columns <- c("year",columns[2:numcols])
  skips <- pick:(pick+nyrs-1)
  fish <- as.data.frame(matrix(NA,nrow=nyrs,ncol=length(columns),
                               dimnames=list(1:nyrs,columns)))
  for (i in 1:nyrs) { # i = 1
    fish[i,] <- scan(file=filename,skip=skips[i],sep=",",nlines = 1,quiet=TRUE)[1:numcols]
  }
  if (verbose) {
    cat("fish \n\n")
    print(head(fish,10))
    cat("\n")
  }
  glb=list(resilience=resilience,spsname=spsname)
  pick <- grep("BIOLOGY",dat)[1] + 1
  if (is.na(pick)) {
    props <- NULL
  } else {
    maxage <- getsingle(dat[pick])
    M <- getsingle(dat[pick+1])
    Linf <- getsingle(dat[pick+2])
    K <- getsingle(dat[pick+3])
    t0 <- getsingle(dat[pick+4])
    Waa <- getsingle(dat[pick+5])
    Wab <- getsingle(dat[pick+6])
    M50a <- getsingle(dat[pick+7])
    deltaM <- getsingle(dat[pick+8])
    sela50 <- getsingle(dat[pick+9])
    deltaS <- getsingle(dat[pick+10])
    steep <- getsingle(dat[pick+11])
    R0 <- getsingle(dat[pick+12])
    ages <- 0:maxage
    nages <- length(ages)
    glb <- list(maxage=maxage,M=M,
                Linf=Linf, K=K, t0=t0,
                Waa=Waa, Wab=Wab,
                M50a=M50a, deltaM=deltaM,
                steep=steep,R0=R0,
                sela50=sela50, deltaS=deltaS,
                resilience=resilience,
                nages=nages,ages=ages,nyrs=nyrs,spsname=spsname)
    columns <- c("age","laa","waa","maa","sela","feca")
    props <- as.data.frame(matrix(NA,nrow=nages,ncol=length(columns),
                                  dimnames=list(ages,columns)))   
    if (property) {
      pick <- grep("PROPERTY",dat)[1] + 1
      for (i in 1:nages) {
        props[i,] <- getvector(dat,pick,sep=",") 
        pick <- pick + 1
      }
    } else {
      # now calculate the properties
      laa <- vB(c(Linf,K,t0),ages)
      waa <- (Waa * laa ^ Wab)/1000
      maa <- logist(M50a,deltaM,ages)
      sela <- logist(sela50,deltaS,ages)
      feca <- sela * waa
      props <- as.data.frame(cbind(ages,laa,waa,maa,sela,feca))
    }
  }
  if (verbose) { 
    cat("biology \n")
    print(unlist(glb))
    cat("\n properties \n")
    print(head(props,10))
  }
  #   dat <- readLines(filename)
  pick <- grep("AGE",dat)[1] + 1
  if (is.na(pick)) {
    agedata <- NULL
  } else {
    yrage <- getsingle(dat[pick])
    numsex <- getsingle(dat[(pick+1)])
    ages <- getvector(dat,(pick+2),sep=",")
    agemax <- max(ages)
    nage <- length(ages)
    naa <- matrix(0,nrow=(yrage*numsex),ncol=(nage+2))
    colnames(naa) <- c("year","sex",ages)
    pick <- pick+2  # i = 1
    for (i in 1:(yrage*numsex)) naa[i,] <- getvector(dat,(pick+i),sep=",")
    rownames(naa) <- naa[,1]
    agedata <- list(yrage=yrage,ages=ages,agemax=agemax,nage=nage,naa=naa)
  }
  pick <- grep("LENGTH",dat)[1] + 1
  if (is.na(pick)) {
    lendata <- NULL
  } else {
    yrlen <- getsingle(dat[pick])
    numsexl <- getsingle(dat[(pick+1)])
    lengths <- getvector(dat,(pick+2),sep=",")
    maxlen <- max(lengths)
    nlength <- length(lengths)
    if (verbose) {
      cat("Number of years of Data: ",yrlen,"\n")
      cat("Number of sexes with Data: ",numsexl,"\n")
      cat("Length classes: ",lengths,"\n")
      cat("Number of length classes ",nlength,"\n")
    }
    nal <- matrix(0,nrow=(yrlen*numsexl),ncol=(nlength+2))
    colnames(nal) <- c("year","sex",lengths)
    pick <- pick+2
    for (i in 1:(yrlen*numsexl)) nal[i,] <- getvector(dat,(pick+i),sep=",")
    rownames(nal) <- nal[,1]
    
    lendata <- list(yrlen=yrlen,lengths=lengths,maxlen=maxlen,
                    nlength=nlength,nal=nal)
  }
  ans <- list(fish=fish,glb=glb,props=props,agedata=agedata,lendata=lendata)
  return(ans)
} # end of readdata






#' @title vB calculates the predicted von Bertalanffy length at age
#'
#' @description vB calculates length at age for the von Bertalanffy curve.
#'
#' @param par is a vector the first three cells of which are Linf, K, and t0
#'    for the VB curve; the fourth parameter will be sigma, the standard
#'    deviation of the normal likelihoods used with the residuals
#' @param ages is a vector of ages
#'
#' @return a vector of predicted lengths for the vector of ages in 'ages'
#' @export
#'
#' @examples
#' ages <- seq(0,20,1)
#' pars <- c(Linf=50,K=0.3,t0=-1.0,sigma=1.0) # Linf, K, t0, sigma
#' cbind(ages,vB(pars,ages))
vB <- function(par,ages) {
  return(par[1] * (1 - exp(-par[2]*(ages-par[3]))))
}
