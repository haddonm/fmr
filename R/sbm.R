
#' @title abdynamicsH calculates the size-based dynamics for the population
#' 
#' @description abdynamicsH is the main function for calculating the size-based
#'     population dynamics, which follows the numbers-at-size. It exports 
#'     the full dynamics or a pared down version of only the predicted cpue 
#'     and predicted numbers-at-size in the catch for the calculation of the'
#'     minimum negative log-likelihood. This is the combination of the 
#'     Log-Normal errors on cpue and the Multinomial errors on the size-
#'     composition of the catch. In this case it uses Rcpp and RcppArmadillo
#'     to speed some of the computations
#'
#' @param abstk the stock derived from makeabstock
#' @param pin the log transformed vector of parameters, pindat column 1
#' @param fish the fishery data in two matrices. The first contains the
#'     year, the LML, catches, and the cpue, and the second contains the 
#'     numbers-at-size in the catch, just for those years with observations.
#' @param const contains biological constants, including those for growth
#' @param glb the global object
#' @param full should only the predicted cpue and numbers-at-size in the
#'     catch be output (full =FALSE), or all of the dynamics (full = TRUE)
#' @param initH the constant initial harvest rate that gives rise, approximately,
#'     to the desired initial depletion. default=0, which implies no initial
#'     depletion. default = 0
#'
#' @return only the predicted values or the full dynamics depending on full
#' @export
#'
#' @examples
#' print("wait on example data sets")
#' # abstk=abstk;pin=pindat[,1];fish=fish;const=const;glb=glb;full=TRUE;initH=0 
abdynamicsH <- function(abstk,pin,fish,const,glb,full=TRUE,initH=0) {
  epin <- exp(pin)
  npin <- length(epin)
  R0 <- epin["LnR0"]
  cpue <- fish[,"cpue"] # needed for q estimation 
  pickce <- which(cpue > 0)   
  catch <- fish[,"catch"]
  yrs <- glb$yrs   
  nyr <- length(yrs) #  define matrices, vectors, and constants  
  midpts <- glb$midpts
  N <- glb$Nclass
  recyrs <- glb$recyrs
  pickrec <- match(recyrs,yrs)
  nrec <- length(pickrec)
  recloc <- which(substr(names(pin),1,1) == "d")
  recdevs <- epin[recloc]
  biol <- abstk$biol
  dyn <- abstk$dyn
  B0 <- abstk$B0
  Nt <- abstk$Nt
  NEbeg <- abstk$NEbeg
  WtL <- biol$WtL
  maturity <- biol$maturity
  emergent <- biol$emergence
  MatWt <- WtL * maturity
  steep <- glb$steep
  surv <- exp(-glb$natM)
  surv2 <- exp(-glb$natM/2.0) 
  LMLs <- unique(fish[,"lml"])  # define selectivity by year
  nlml <- length(LMLs)
  sel <- matrix(0,nrow=N,ncol=nyr,dimnames=list(midpts,yrs))
  for (i in 1:nlml) { # i = 1
    pick <- which(fish[,"lml"] == LMLs[i])
    L50 <- fish[pick[1],"lml"] 
    delta <- epin["seldelta"]
    sel[,pick] <- knifelogistic(L50,delta,midpts,knifeedge=LMLs[i])
  }
  gpar <- c(epin["Linf"],const["K"],const["cvG"])
  G <- STMvB(p=gpar,mids=midpts)
  exB0 <- sum(NEbeg[,1] * sel[,1] * WtL)/1e6
  expBU <- numeric(nyr)
  NumN <- matrix(0,nrow=N,ncol=nyr,  # midyr ExploitNAS 
                 dimnames=list(rownames(Nt),colnames(Nt)))
  catchN <- NumN  
  if (initH > 0) {  # deplete the unfished stock if required
    #   for (yr in 2:nyr) {    #  yr=2
    #     yrm1 <- yr - 1
    #     NumN[,yr] <- (surv2 * (G %*% Nt[,yrm1])) #  yrm1 = yr - 1
    #     fishH <- sel[,1] *  initH  
    #     Nt[,yr] <- (surv2 * (1-(fishH)) * NumN[,yr])
    #     matB[yr] <- sum(Nt[,yr] * MatWt)/1e6
    #     recN <- bh(steep,R0,B0,matB[yr],devR=1) 
    #     Nt[1,yr] <- recN
    #   } # year loop
    #   Nt[,1] <- Nt[,nyr]
    #   expB[1] <- sum(Nt[,1] * sel[,1] * WtL)/1e6
    #   expBU[1] <- expB[1]
    #   matB[1] <- sum(Nt[,1] * MatWt)/1e6
    #   depl[1] <- matB[1]/B0
    #   # exdepl[1] <- expB2/exB0
  }
  countrec <- 0         
  for (yr in 2:nyr) { # nyr-1 DYNAMICS  # yr=nyr  # yr=2
    yrm1 <- yr - 1
    NumN[,yr] <- (surv2 * (G %*% Nt[,yrm1])) # NumN=mid-year NAS, before fishing
    expB1 <- sum(NumN[,yr] * sel[,yr] * WtL)/1e6 # prefish exploitB
    expBU[yr] <- expB1
    dyn[yr,"harvest"] <- min(catch[yr]/expB1,0.8)
    fishH <- sel[,yr] * dyn[yr,"harvest"]
    catchN[,yr] <- (fishH * NumN[,yr])
    Nt[,yr] <- (surv2 * ((1-(fishH)) * NumN[,yr]))
    dyn[yr,"predC"] <- sum(catchN[,yr] * WtL)/1e6
    expB2 <- sum(Nt[,yr] * sel[,yr] * WtL)/1e6
    dyn[yr,"exploitB"] <- (expB1 + expB2)/2.0  # half pre- + post fishing expB
    dyn[yr,"matureB"] <- sum(Nt[,yr] * MatWt)/1e6
    devR <- 1.0
    if (yr %in% pickrec) {
      countrec <- countrec + 1
      devR <- recdevs[countrec]
    }
    dyn[yr,"recruit"] <- bh(dyn[yr,"matureB"],R0,B0,steep,devR=devR) 
    # cat(yrs[yr],yr,countrec,devR,recruit[yr],"\n")
    Nt[1,yr] <- dyn[yr,"recruit"]
  } # year loop
  dyn[,"deplet"] <-dyn[,"matureB"]/B0
  # exBdepl <- expB/exB0
  lambda <- glb$lambda   
  qest <- epin["qest"]
  expB <- dyn[2:nyr,"exploitB"]
  dyn[,"predCE"] <- c(NA,qest * (expB ^ lambda))
  if (full) {
    ans <- list(Nt=Nt,catchN=catchN,dyn=dyn,G=G,B0=B0,initH=initH,
                sel=sel)
    
  } else {
    ans <- list(predce=dyn[,"predCE"],predC=dyn[,"predC"],catchN=catchN,
                Nt=Nt,dyn=dyn)
  }
  return(ans)
} # end of abdynamicsH



#' @title getabSBMinit initiates a single sex abalone Size-Based Model
#' 
#' @description getabSBMinit estimates the B0, N0, and NEbeg (initial start of 
#'     year exploitable numbers-at-size), for a single sex abalone stock. It 
#'     requires the unfished recruitment R0, the globals object, glb, the 
#'     biol object of the species being modelled.
#'
#' @param R0 the nominal scale estimated unfished recruitment
#' @param glb the globals object containing natM and the midpts
#' @param biol the biol list containing growth, maturity and weight-at-length
#'
#' @returns a list of B0 the unfished spawning biomass, and N0 the unfished 
#'     numbers-at-size
#' @export
#'
#' @examples
#' data(absaudata)
#' abstk <- makeabstock(fish=absaudata$fish,const=absaudata$const,
#'                      glb=absaudata$glb,pin=absaudata$pindat)
#' initSBM <- getabSBMinit(exp(absaudata$pindat["LnR0",1]),glb=absaudata$glb,
#'                       biol=abstk$biol)
#' print(initSBM$B0)
#' print(round(initSBM$N0))
getabSBMinit <- function(R0,glb,biol) {
  surv <- exp(-glb$natM)
  mids <- glb$midpts
  G <- biol$G
  I <- makeUnit(nrow(G)) # generates a square Unit matrix same size as G
  S <- makeUnit(nrow(G),surv)
  rec <- numeric(glb$Nclass); rec[1] <- 1 # generate mass per recruit
  N0 <- (solve(I - G %*% S)) %*% rec
  A0 <-  sum(biol$maturity * biol$WtL * N0)/1e06
  B0 <- R0 * A0
  initN <- N0 * R0
  return(list(B0=B0,N0=initN))
} # end of getabSBMinit

#' @title getConst extracts 'nb' numbers from a line of text
#'
#' @description getConst parses a line of text and extracts 'nb' pieces of
#'     text as numbers. If there are missing values this will stop the 
#'     function with a warning. If RStudio drops into debug, toggle 'message
#'     only' under Debug/on Error to avoid that, and put NA where the data
#'     are missing
#'
#' @param inline text line to be parsed, usually obtained using readLines
#' @param nb the number of numbers to extract
#' @param index which non-empty object to begin extracting from?
#'
#' @return a vector of length 'nb'
#' @export
#'
#' @examples
#'   inline <- "name , 31,32,33"
#'   getConst(inline=inline,nb=3,index=2)
#'   inline <- "anothername , 31,32,NA"
#'   getConst(inline=inline,nb=3,index=2)
#'   inline <- "finalname , 31,32,"  # missing values not allowed   
#'  # getConst(inline=inline,nb=3,index=2) # unhash and try to see error
getConst <- function(inline,nb,index=2) { # inline=dat[linenum];nb=3;index=2
  ans <- numeric(nb)
  tmp <- removeEmpty(unlist(strsplit(inline,",")))
  if (length(tmp) < (nb+(index-1))) {
    stop(paste("possible problem with data",tmp[1],
               "missing comma or missing NA/data?",sep=" "),"\n")
  } else {
    count <- 0
    for (j in index:(nb+index-1)) {  # j=4
      count <- count + 1
      if (tmp[j] == "NA") {
        ans[count] <- NA
      } else {
        ans[count] <- as.numeric(tmp[j])
      }
    }
  }
  return(ans)
}   # end getConst

#' @title getrmse calculates the rmse of the input 'cpue' series
#'
#' @description getrmse calculates the root mean square error (rmse) of the 
#'     input invar series (defaults to 'cpue') against an input 'year' time 
#'     series. This is primarily designed to generate an alternative estimate 
#'     of the intrinsic variability of a cpue time-series to that which may be 
#'     obtained from a cpue standardization. This enables one to use the Francis
#'     2011 weighting scheme for indices of relative abundance, which he derived 
#'     from Clarke and Hare, 2006.
#'     
#' @references Clark, W.G., and S.R. Hare (2006) Assessment and management of 
#'     Pacific halibut: data, methods, and policy. \emph{Scientific Report 83. 
#'     International Pacific Halibut Commission}, Seattle, Washington. 104p. 
#'     
#' @references Francis, R.I.C.C. (2011) Data weighting in statistical fisheries 
#'     stock assessment models. \emph{Canadian Journal of Fisheries and Aquatic 
#'     Sciences} \bold{68}: 1124-1138.
#'     
#'
#' @param indat the matrix, spmdat, or data.frame containing both a 'year'
#'     column and an invar column (default to 'cpue')
#' @param invar the column name of the variable whose rmse is wanted; defaults 
#'     to 'cpue'
#' @param inyr the column name that points to the 'year' name
#' @param natlog should we estimate the rmse of log-transformed data. default
#'     = TRUE because we are usually dealing with cpue.
#' @return a list of the rmse and the loess predicted values of the invar 
#'     for each year in the time-series
#' @export
#'
#' @examples
#' year <- 1986:1994
#' cpue <- c(1.2006,1.3547,1.0585,1.0846,0.9738,1.0437,0.7759,1.0532,1.284)
#' dat <- as.matrix(cbind(year,cpue))
#' getrmse(dat,invar="cpue")  # should be 0.08265127
#' getrmse(dat,invar="cpue")$rmse
getrmse <- function(indat,invar="cpue",inyr="year",natlog=TRUE){
  if (iscol(inyr,indat) & iscol(invar,indat)) {
    nyr <- dim(indat)[1]
    predictedCE <- rep(NA,nyr)
    varloc <- grep(invar,colnames(indat))
    nvar <- length(varloc)
    if (nvar > 1) {
      obsvar <- rep(NA,nyr)
      for (i in 1:nvar) {
        pick <- which(indat[,varloc[i]] > 0)
        obsvar[pick] <- indat[pick,varloc[i]]
      }
    } else {
      obsvar <- indat[,varloc]
    }
    picky <- which(obsvar > 0)
    if (natlog) {
      model <- loess(log(obsvar[picky]) ~ indat[picky,inyr])
    } else {
      model <- loess(obsvar[picky] ~ indat[picky,inyr])
    }
    predictedCE[picky] <- model$fitted
    rmse <- sqrt(sum(model$residuals^2)/model$n)
    return(list(rmse=rmse,predictedCE=predictedCE,loessmod=model))
  } else {
    warning("Use the invar and inyr arguments to name both 'year' and 'cpue'\n")
  }
} # end of getrmseCE

#' @title getSBMinit initiates a Size-Based Model
#' 
#' @description getSBMinit estimates the B0, N0, and NEbeg (initial start of 
#'     year exploitable numbers-at-size), all by gender. It requires the 
#'     unfished recruitment R0, the globals object, glb, the biol object, and 
#'     the sex being modelled.
#'
#' @param inR0 the nominal scale estimated unfished recruitment
#' @param glb the globals object containing natM and the midpts
#' @param biol the biol list containing growth, maturity and weight-at-length
#' @param sex which gender is being analysed? Either 1, females or 2, males
#'
#' @returns a list of B0 the unfished spawning biomass, and N0 the unfished 
#'     numbers-at-size
#' @export
#'
#' @examples
#' print("Wait on datasets")
getSBMinit <- function(inR0,glb,biol,sex) { # assumes glb inR0 = par["R0"]
  # inR0 =exp(pin["LnR0",1]);biol=biol;sex=1;glb=glb;
  surv <- exp(-glb$natM)
  mids <- glb$midpts
  G <- biol$G[,,sex]
  I <- makeUnit(nrow(G)) # generates a square Unit matrix same size as G
  S <- makeUnit(nrow(G),surv)
  rec <- numeric(glb$Nclass); rec[1] <- 1
  N0 <- (solve(I - G %*% S)) %*% rec
  A0 <-  sum(biol$maturity[,sex] * biol$WtL[,sex] * N0)/1000.0
  B0 <- inR0 * A0
  initN <- N0 * inR0
  return(list(B0=B0,N0=initN))
} # end of getSBMB0



#' @title getsingleNum find a line of text and extracts a single number
#'
#' @description getsingleNum uses grep to find an input line. If the variable
#'     being searched for fails then NULL is returned
#'
#' @param varname the name of the variable to get from intxt
#' @param intxt text to be parsed, usually obtained using readLines
#'
#' @return a single number or, if no value is in the data file a NULL
#' @export
#'
#' @examples
#' \dontrun{
#'  txtlines <- c("replicates, 100","Some_other_text, 52")
#'  getsingleNum("replicates",txtlines)
#'  getsingleNum("eeplicates",txtlines)
#'  getsingleNum("other",txtlines)
#' }
getsingleNum <- function(varname,intxt) {
  begin <- grep(varname,intxt)
  if (length(begin) > 0) {
    return(as.numeric(getConst(intxt[begin],1)))
  } else {
    return(NULL)
  }
}

#' @title getStr obtains a string from an input text line
#'
#' @description  getStr obtains a string from an input text line in
#'     which any parts are separated by ','. Then, after ignoring the
#'     first component, assumed to be a label, it returns the first
#'     nb parts.
#'
#' @param inline input text line with components separated by ','
#' @param nb number of parts to return
#' @param index which part to start returning to the program, default=2
#'
#' @return a vector of character string(s)
#' @export
#'
#' @examples
#'   txt <- "runlabel, development_run, label for this particular run"
#'   getStr(txt,1)
getStr <- function(inline,nb,index=2) { # inline=indat[begin];nb=4;start=1
  tmp <- removeEmpty(unlist(strsplit(inline,",")))
  finish <- nb + (index - 1)
  outconst <- as.character(tmp[index:finish])
  return(outconst)
} # end of getStr

#' @title knifelogistic a Logistic selectivity function with knifeedge option
#'
#' @description knifelogistic a logistic selectivity function with the option
#'     of including a knifeedge to allow for LML. This uses the
#'     logistic function 1/(1+exp(-log(19.0)*(lens-inL50)/(delta))),
#'     where delta = inL95 - inL50. This explicitly defines the SM50
#'     but uses delta (which is SM95-SM50) as the second. This ensures
#'     that when adding variation to parameters, to vary between
#'     populations, when SM95 and SM50 are close together it is not
#'     possible for SM50 to become larger than SM95. Be careful using the
#'     knifeedge option. Strictly knifeedge selectivity would entail the
#'     selectivity values being zero up to the knife-edge and then being 1.0.
#'     This is not what happens here. Instead the knifeedge option literally
#'     sets all values to zero at and below the value of knifeedge but leaves
#'     any curve above that value as it is. Hence this is not strict knife-edge
#'     selectivity. However, it does provide a selectivity curve that reflects
#'     the selectivity of a fishery that uses a Legal Minimum Length or Size.
#'
#' @param L50 is the length at 50 percent selection
#' @param delta is the difference between the 95 percent selection and
#'     the 50 percent selection
#' @param lens a vector of lengths for which the logistic maturity value
#'     will be calculated
#' @param knifeedge defaults to 0. If knifeedge is set to a particular
#'     length then the selectivity less than the value of knifeedge is set
#'     to zero.
#' @param maxLML default = 0. The parameter is to allow for a slot selectivity
#'
#' @return A vector of length(lens) containing the predicted selectivity at
#'    length values
#' @export
#'
#' @examples
#' \dontrun{
#' inL50 <- 100.0
#' delta <- 8.0
#' lens <- seq(2,210,2)
#' select <- knifelogistic(L50=inL50,delta,lens)
#' select <- knifelogistic(L50=inL50,delta,lens,knifeedge=105)
#' select <- knifelogistic(L50=inL50,delta,lens,maxLML=185)
#' }
knifelogistic <- function(L50,delta,lens,knifeedge=0,maxLML=0) {
  ans <- 1/(1+exp(-log(19.0)*(lens-L50)/(delta)))
  if (knifeedge > 0) {
    pick <- which(lens < knifeedge)
    if (length(pick) > 0) ans[pick] <- 0.0
  }
  if (maxLML > 0) {
    pick <- which(lens > maxLML)
    if (length(pick) > 0) ans[pick] <- 0.0
  }
  return(ans)
} # end of knifelogistic

#' @title makeabstock defines an abalone stock's properties prior to assessment
#' 
#' @description makeabstock characterizes an abalone stock and fishery prior to
#'     its use in a stock assessment. The output list contains a 'biol' list 
#'     containing the weight, maturity and emergence from crypsis for each 
#'     size/length or size-class. It also contains the growth transition matrix.
#'     It also contain the unfished biomass B0, then the initial fishery 
#'     dynamics 'dyn', containing catch, cpue, mature and exploitable biomass,
#'     etc. Then 'Nt' the matrix that will contain the numbers-at-time, 
#'     containing the initial numbers-at-size in the first column. The 'NEbeg',
#'     the start of year numbers-at-size of exploitable abalone. The final 
#'     'fishery' object contains each year's selectivity (to allow for LML 
#'     changes through time), the unique LMLs that occur, the years in the data,
#'     the number of years, the survivorship from natural mortality and that 
#'     from half the natural mortality. 
#'
#' @param fish the fishery dependent data including year, LML, catch, and cpue
#' @param const parameters for weight-at, maturity-at, and emergence-at-size,
#'     and the von Bertalanffy parameters Linf, K, and CV.
#' @param glb the global constants object
#' @param pin the vector of initial parameters and whether they are fitted or 
#'     fixed
#' @param listtoenv default=FALSE, should list2env be used to put all list 
#'     objects into the .GlobalEnv.
#'
#'
#' @returns a list of size objects, biol, B0, dyn, Nt, NEbeg, and fishery
#' @export
#'
#' @examples
#' data(absaudata)
#' abstk <- makeabstock(fish=absaudata$fish,const=absaudata$const,
#'                      glb=absaudata$glb,pin=absaudata$pindat,listtoenv=FALSE)
#' str(abstk)
makeabstock <- function(fish,const,glb,pin,listtoenv=FALSE) {
  mids <- glb$midpts
  Nclass <- glb$Nclass
  WtL <- const["Wta"] * mids ^ const["Wtb"]
  maturity <- logistic(L50=const["M50"],delta=const["Mdelta"],depend=mids)
  emergence <- logistic(L50=const["E50"],delta=const["Edelta"],depend=mids)
  G <- STMvB(p=c(const["Linf"],const["K"],const["cvG"]),mids=mids)
  biol <- list(WtL=WtL,maturity=maturity,emergence=emergence,G=G)
  # do structure
  yrs <- fish[,"year"]
  nyr <- length(yrs)
  defmat <- matrix(0,nrow=Nclass,ncol=nyr,dimnames=list(mids,yrs))
  defvect <- numeric(nyr) ; names(defvect) <- yrs
  Nt <- NEbeg <- defmat 
  sel <- matrix(0,nrow=Nclass,ncol=nyr,dimnames=list(mids,yrs))
  colnam <- c("catch","predC","cpue","predCE","matureB","exploitB","midyrexpB",
              "deplet","recruit","harvest")
  dyn <- matrix(0,nrow=nyr,ncol=length(colnam),dimnames=list(yrs,colnam))
  dyn[,"catch"] <- fish[,"catch"]
  dyn[,"cpue"] <- fish[,"cpue"]
  surv <- exp(-glb$natM)
  surv2 <- exp(-glb$natM/2.0) 
  LMLs <- unique(fish[,"lml"])  # define selectivity by year
  nlml <- length(LMLs)
  for (i in 1:nlml) { # i = 1
    pick <- which(fish[,"lml"] == LMLs[i])
    L50 <- fish[pick[1],"lml"] 
    delta <- exp(pin["seldelta",1])
    sel[,pick] <- knifelogistic(L50,delta,mids,knifeedge=LMLs[i])
  }
  fishery <- list(sel=sel,LMLs=LMLs,yrs=yrs,nyr=nyr,surv=surv,surv2=surv2)
  #equilpops
  R0 <- exp(pin["LnR0",1])
  init <- getabSBMinit(R0=R0,glb=glb,biol=biol)
  B0 <- init$B0
  dyn[1,"matureB"] <- B0
  Nt[,1] <- init$N0 
  NEbeg[,1] <- init$N0 * sel[,1] 
  dyn[1,"exploitB"] <- sum(NEbeg[,1] * WtL)/1e06
  dyn[1,"deplet"] <- 1.0
  dyn[1,"recruit"] <- R0
  stock <- list(biol=biol,B0=B0,dyn=dyn,Nt=Nt,NEbeg=NEbeg,fishery=fishery)
  if (listtoenv) list2env(stock,envir=.GlobalEnv)
  return(invisible(stock))
} # end of makeabstock

#' @title makestock set up the size-based model stock 
#' 
#' @description makestock initiates a size-based model of a stock by using the 
#'     biological and fishery parameters to estimate the required maturity, 
#'     selectivity, weight-at-length, unfished numbers-at-size, unfished mature
#'     and exploitable biomass.
#'
#' @param fish the fishery dependent data year, LML, catches cpue
#' @param const the biological constants relating to maturity, growth, and 
#'     weight-at-length
#' @param glb the globals object
#' @param pin the initial starting values of the estimable parameters and 
#'     whether they are estimated or not.
#'
#' @returns a list of the biological properties, B0, matB, expB, Nt, NEbeg,
#'     depl, and fishery details.
#' @export
#'
#' @examples
#' print("wait on data sets")
makestock <- function(fish,const,glb,pin) { 
  # fish=fish; sizecomp=sizecomp; const=const; glb=glb; pin=pin
  # do biology
  gender <- c("F","M")
  mids <- glb$midpts
  Nclass <- glb$Nclass
  WtL <- matrix(0,nrow=Nclass,ncol=2,dimnames=list(mids,gender))
  WtL[,1] <- const["FWta"] * mids ^ const["FWtb"]
  WtL[,2] <- const["MWta"] * mids ^ const["MWtb"]
  maturity <- matrix(0,nrow=Nclass,ncol=2,dimnames=list(mids,gender))
  maturity[,1] <- logistic(L50=const["FM50"],delta=const["Fdelta"],
                           depend=mids)
  maturity[,2] <- logistic(L50=const["MM50"],delta=const["Mdelta"],
                           depend=mids)  
  G <- array(data=0,dim=c(Nclass,Nclass,2),dimnames=list(mids,mids,gender))
  G[,,1] <- STMvB(p=c(const["FLinf"],const["FvbK"],const["FCVvb"]),mids=mids)
  G[,,2] <- STMvB(p=c(const["MLinf"],const["MvbK"],const["MCVvb"]),mids=mids)
  biol <- list(WtL=WtL,maturity=maturity,G=G)
  # do structure
  yrs <- fish[,"year"]
  nyr <- length(yrs)
  yrlab <- c((yrs[1]-1),yrs)
  defarr <- array(0,dim=c(Nclass,(nyr+1),2),dimnames=list(mids,yrlab,gender))
  defvect <- numeric(nyr+1) ; names(defvect) <- yrlab
  Nt <- NEbeg <- defarr 
  sel <- matrix(0,nrow=Nclass,ncol=nyr,dimnames=list(mids,yrs))
  colnam <- c("catch","predC","cpue","predCE","matureB","exploitB","midyrexpB",
              "deplet","recruit","harvest")
  dyn <- matrix(0,nrow=(nyr+1),ncol=length(colnam),dimnames=list(yrlab,colnam))
  dyn[,"catch"] <- c(NA,fish[,"catch"])
  dyn[,"cpue"] <- c(NA,fish[,"cpue"])
  surv <- exp(-glb$natM)
  surv2 <- exp(-glb$natM/2.0) 
  LMLs <- unique(fish[,"lml"])  # define selectivity by year
  nlml <- length(LMLs)
  for (i in 1:nlml) { # i = 1
    pick <- which(fish[,"lml"] == LMLs[i])
    L50 <- exp(pin["selL50",1]); delta <- exp(pin["seldelta",1])
    sel[,pick] <- knifelogistic(L50,delta,mids,knifeedge=LMLs[i])
  }
  fishery <- list(sel=sel,LMLs=LMLs,yrs=yrs,nyr=nyr,surv=surv,surv2=surv2)
  #equilpops
  B0 <- numeric(2); names(B0) <- gender
  R0 <- exp(pin["LnR0",1])
  for (sex in 1:2) {
    init <- getSBMinit(inR0=R0/2.0,glb=glb,biol=biol,sex=sex)
    B0[sex] <- init$B0
    dyn[1,"matureB"] <- dyn[1,"matureB"] + B0[sex]
    Nt[,1,sex] <- init$N0 
    NEbeg[,1,sex] <- init$N0 * sel[,1] 
    dyn[1,"exploitB"] <- dyn[1,"exploitB"] + 
                         sum(NEbeg[,1,sex] * WtL[,sex])/1000.0
  }
  dyn[1,"deplet"] <- 1.0
  dyn[1,"recruit"] <- R0
  stock <- list(biol=biol,B0=B0,dyn=dyn,Nt=Nt,NEbeg=NEbeg,fishery=fishery)
} # end of makestock


#' @title plotSBMdynamics plots the implied dynamics from the assessment
#' 
#' @description plotSBMdynamics generates a 3 x 2 plot of the cpue imposed on
#'     the predicted cpue, the exploitable biomass, the mature biomass,
#'     the harvest rate, the recruitment levels, and the mature biomass
#'     depletion levels. 
#'
#' @param rundir the directory for all ctrl, data, and output files.
#' @param dyn the matrix of the dynamics out of the 'dynamics' function 
#'     that uses the optimum estimate of the parameters
#' @param fish the fishery data used in the assessment
#' @param pars the back-transformed parameters from the assessment
#' @param glb the globals object
#' @param console should the plot go to the console or be saved? Default=TRUE 
#'
#' @return nothing but it does add a plot to the console or to rundir
#' @export
#'
#' @examples
#' \dontrun{
#'     plotSBMdynamics(rundir,dyn,fish,exp(pindat[,1]),glb,console = TRUE)
#' }
#' #  rundir=rundir; console=TRUE; pars=allpin
#' # rundir=rundir;dyn=outdyn$dyn;fish=fish;pars=exp(out$optpar);glb=glb;
#' # fisindex=NULL; predfisindex=NULL; console = TRUE
plotSBMdynamics <- function(rundir,dyn,fish,pars,glb,console=TRUE) {
  # rundir=""; dyn=out$dyn; fish=fish;pars=exp(pindat[,1]);console=TRUE
  if (console) {   filen="" } else {
    filen <- filenametopath(rundir,"optimum_population_dynamics.png")
  }
  numplot <- c(4,2)
  npar <- length(pars)
  years <- glb$yrs
  nyrs <- length(years)
  recyrs <- glb$recyrs
  deviates <- rep(1.0,length(years))
  pickyr <- match(recyrs,years)
  devpars <- npar - length(recyrs) + 1
  deviates[pickyr] <- pars[devpars:npar]
  plotprep(width=9,height=8,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=numplot,byrow=FALSE,margin=c(0.3,0.45,0.05,0.05))
  ymax <- getmax(c(dyn[,"predCE"],fish[,"cpue"]))
  plot1(years,dyn[,"predCE"],lwd=2,col=2,ylab="CPUE",defpar=FALSE,maxy=ymax)
  lines(years,fish[,"cpue"],lwd=2,col=1)
  legend("topright",c("Observed","Predicted"),col=c("black","red"),
         lwd=3,cex=1,bty="n")
  # residuals
  tmp <- dyn[,"cpue"]/dyn[,"predCE"]
  resid <- tmp[which(tmp > 0)]
  residyr <- as.numeric(names(resid))
  nres <- length(resid)  
  maxy <- getmax(resid)
  miny <- getmin(resid,mult=1.2)
  plot(residyr,resid,type="p",pch=16,cex=1,ylab="CPUE Residuals",xlab="",
       ylim=c(miny,maxy),panel.first=grid())
  #  plot1(residyr,resid,type="p",pch=16,cex=1,ylab="CPUE Residuals",defpar=FALSE)
  abline(h=1.0,lwd=1,col="darkgrey")
  for (i in 1:nres) {
    x <- c(residyr[i],residyr[i])
    lines(x,c(resid[i],1.0),lwd=1,col=2)
  }
  plot1(years,fish[,"catch"],lwd=2,ylab="Catch (t)",defpar=FALSE)
  lines(years,dyn[,"predC"],lwd=2,col=2)
  plot1(years,dyn[,"harvest"],lwd=2,ylab="Harvest Rate",defpar=FALSE)
  # plot1(years,dyn[,"exBdepl"],lwd=2,ylab="Exploitable Biomass",defpar=FALSE)
  # label <- paste0("Depl = ",round(dyn[nyrs,"exBdepl"],3))
  #  mtext(text=label,side=3,line=-1.1,cex=1.0)
  plot1(years,dyn[,"deplet"],lwd=2,ylab="Mature Biomass depletion",defpar=FALSE)
  label <- paste0("Depl = ",round(dyn[nyrs,"deplet"],3))
  mtext(text=label,side=3,line=-1.1,cex=1.0)
  recs <- dyn[,"recruit"]
  pickpred <- match(glb$recyrs,as.numeric(names(recs)))  
  plot1(years,dyn[,"recruit"],lwd=2,ylab="Recruitment",defpar=FALSE)
  lines(years[pickpred],recs[pickpred],lwd=2,col="red")
  spb <- dyn[,"matureB"]
  nodev <- bh(spb,pars[1],spb[1],glb$steep,devR=1)
  lines(years,nodev,lwd=2,col=4)
  legend("topright",c("No Deviates","+ Deviates"),col=c("blue","red"),
         lwd=3,cex=1,bty="n")
  plot1(years,deviates,lwd=2,ylab="Recruitment deviates",defpar=FALSE,col="blue")
  lines(years[pickpred],deviates[pickpred],lwd=2,col="red")
  abline(h=1,lwd=1,col="darkgrey")
  if (!console) {
    caption <- paste0("The population dynamics implied by the assessment. ",
                      "Note different years used on residual plot.")
    addplot(filen,rundir=rundir,category="Assessment",caption)
  }
} # end of plotSBMdynamics

#' @title SBMctrltemplate creates a draft control file for Size-Based Models
#' 
#' @description SBMctrltemplate simplifies the generation of the working 
#'     control file for fitting a size-based model to population data. The
#'     file must be a .csv file. We use Me instead of M to denote natural
#'     mortality because M is not unique, which means R can become confused. 
#'
#' @param rundir the full path to the working directory being used
#' @param infile what name will be given to the control .csv file
#' 
#' @seealso{
#'  \link{readSBMctrl}
#' }
#'
#' @return nothing but it does generate a csv file into the rundir
#' @export
#'
#' @examples
#' \dontrun{
#'  rundir <- tempdir()
#'  SBMctrltemplate(rundir,infile="draftctrlLBM.csv")
#'  setup <- readSBMctrl(rundir,infile="draftctrlLBM.csv")
#'  # indir=rundir; filename="control2.csv"  #
#' }
SBMctrltemplate <- function(rundir,infile="controlLBM.csv") { 
  filename <- pathtopath(rundir,infile)
  cat("DESCRIPTION \n",
      file=filename,append=FALSE)
  cat("Control file containing details of a particular run. Modify the  \n",
      file=filename,append=TRUE)
  cat("contents to suit your own situation. In particular. modify the  \n",
      file=filename,append=TRUE)
  cat("contents of this description to suit the population being fitted \n",
      file=filename,append=TRUE)
  cat("\n\n",file=filename,append=TRUE)
  cat("START \n",file=filename,append=TRUE)
  cat("runlabel, An_Example, label for particular run \n",
      file=filename,append=TRUE)
  cat("datafile, dataSBM.csv, name of input data file \n",
      file=filename,append=TRUE)
  cat("pinfile, pinSBM.csv, name of the input parameter file  \n",
      file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("CONSTANTS, 31, \n",file=filename,append=TRUE)
  cat("minsc, 5, centre of minimum size class \n",file=filename,append=TRUE)
  cat("cw, 0.5, class width mm \n",file=filename,append=TRUE)
  cat("Nclass, 31, number of size classes \n",file=filename,append=TRUE)
  cat("FWta, 0.00156937, weight-at-length intercept \n",file=filename,append=TRUE)
  cat("FWtb, 2.76902, weight-at-length gradient \n",file=filename,append=TRUE)
  cat("FM50, 10, SAM50% \n",file=filename,append=TRUE)
  cat("Fdelta, 1.5, SAM95% \n",file=filename,append=TRUE)
  cat("FLinf, 16.0, maximum vb length \n",file=filename,append=TRUE)
  cat("FvbK, 0.4, female vonB K parameter  \n",
      file=filename,append=TRUE)
  cat("FCVvb, 0.1, maximum StDev for variation of female growth \n",
      file=filename,append=TRUE)
  cat("MWta, 0.00156937, weight-at-length intercept \n",
      file=filename,append=TRUE)
  cat("MWtb, 2.76902, weight-at-length gradient \n",file=filename,append=TRUE)
  cat("MM50, 10, SAM50% \n",file=filename,append=TRUE)
  cat("Mdelta, 1.5, SAM95% \n",file=filename,append=TRUE)
  cat("MLinf, 17.0, maximum vb length \n",file=filename,append=TRUE)
  cat("MvbK, 0.375, female vonB K parameter  \n",
      file=filename,append=TRUE)
  cat("MCVvb, 0.1, maximum StDev for variation of male growth \n",
      file=filename,append=TRUE)
  cat("recyr1, 1995, 1st year of rec devs  \n",file=filename,append=TRUE)
  cat("recyr2, 2017, last year of rec devs  \n",file=filename,append=TRUE)
  cat("natM, 0.69, natural mortality \n",file=filename,append=TRUE)
  cat("sigmaR, 0.6, permissible variation in recruitment deviates  \n",
      file=filename,append=TRUE)
  cat("steep, 0.5, recruitment steepness  \n",file=filename,append=TRUE)
  cat("wtsc, 0, weighting given to size-composition data  \n",
      file=filename,append=TRUE)  
  cat("sigce, 0, weighting given to cpue likelihood 0=use rmse from loess  \n",
      file=filename,append=TRUE)  
  cat("lambda, 1, exponent for cpue vs expB relationship  \n",
      file=filename,append=TRUE)  
  # cat("Radjy1, 1986, rec bias adjustment year 1 \n",file=filename,append=TRUE)
  # cat("Radjy2, 2003, rec bias adjustment year 2 \n",file=filename,append=TRUE)
  # cat("Radjy3, 2007, rec bias adjustment year 3 \n",file=filename,append=TRUE)
  # cat("Radjy4, 2016, rec bias adjustment year 4 \n",file=filename,append=TRUE)
  # cat("adjmax, 0.775, maximum bias adjustment \n",file=filename,append=TRUE)
  # cat("initdepl, 1.0, maximum MaxDL constraint \n",file=filename,append=TRUE)
  cat("omega, 1,1,0,0,emphasis on each data stream, ce, sizecomp, fisindex, fiscomp \n",
      file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
} # end of SBMctrltemplate



#' @title SBMdatatemplate creates a draft data file for LBM
#' 
#' @description SBMdatatemplate simplifies the generation of the working 
#'     data file for fitting a size-based model to population data. The
#'     file must be a .csv file. The data input includes the catches, the
#'     cpue, any FIS data, and size-composition data. 
#'
#' @param rundir the full path to the working directory being used
#' @param outfile what name will be given to the data .csv file
#' 
#' @seealso{
#'  \link{readSBMdata}
#' }
#'
#' @return nothing but it does generate a csv file into the rundir
#' @export
#'
#' @examples
#' \dontrun{
#'  rundir <- tempdir()
#'  SBMdatatemplate(rundir,outfile="dataSBM.csv")
#'  setup <- readSBMdata(rundir,outfile="dataSBM.csv")
#' }
SBMdatatemplate <- function(rundir,outfile="dataSBM.csv") { 
  # rundir=rundir; filename="dataSBM.csv"
  filename <- pathtopath(rundir,outfile)
  cat("A_rock_Lobster_example, label for particular run \n\n\n",
      file=filename,append=FALSE)
  cat("FISHERY, 34, 4, \n",file=filename,append=TRUE)
  cat("year, lml, catch, cpue \n",file=filename,append=TRUE)
  cat("1988,  9, 0.21, NA,, \n",file=filename,append=TRUE)
  cat("1989,  9, 0.22, NA,, \n",file=filename,append=TRUE)  
  cat("1990,  9, 0.08, NA,, \n",file=filename,append=TRUE)  
  cat("1991,  9, 0.81, NA,, \n",file=filename,append=TRUE)  
  cat("1992,  9, 5.06, NA,, \n",file=filename,append=TRUE)  
  cat("1993,  9, 0.29, NA,, \n",file=filename,append=TRUE)  
  cat("1994,  9, 3.18, NA,, \n",file=filename,append=TRUE)  
  cat("1995,  9, 37.83, NA ,, \n",file=filename,append=TRUE) 
  cat("1996,  9, 56.37, NA,, \n",file=filename,append=TRUE)  
  cat("1997,  9, 49.87, NA,, \n",file=filename,append=TRUE)  
  cat("1998,  9, 77.64, NA,, \n",file=filename,append=TRUE)  
  cat("1999,  9, 129.3, NA,, \n",file=filename,append=TRUE)  
  cat("2000,  9, 152.45, NA,, \n",file=filename,append=TRUE)  
  cat("2001,  9, 206.48, NA,, \n",file=filename,append=TRUE)  
  cat("2002,  9, 116.67, NA,, \n",file=filename,append=TRUE)  
  cat("2003,  9, 103.84, 1.0926,,, \n",file=filename,append=TRUE)  
  cat("2004,  9, 186.59, 1.1498,,, \n",file=filename,append=TRUE)  
  cat("2005,  9, 140.83, 1.2342,,, \n",file=filename,append=TRUE)  
  cat("2006,  9, 203.26, 1.0603,,, \n",file=filename,append=TRUE)  
  cat("2007,  9, 238.51, 1.1737,,, \n",file=filename,append=TRUE)  
  cat("2008,  9, 241.87, 1.1877,,, \n",file=filename,append=TRUE)   
  cat("2009,  9, 198.05, 1.3129, \n",file=filename,append=TRUE) 
  cat("2010,  9.5, 117.93, 0.9446, \n",file=filename,append=TRUE)  
  cat("2011,  9.5, 141.16, 0.8929, \n",file=filename,append=TRUE)  
  cat("2012,  9.5, 153.95, 1.0906, \n",file=filename,append=TRUE) 
  cat("2013,  9.5, 177.58, 1.1302, \n",file=filename,append=TRUE)  
  cat("2014,  9.5, 176.78, 1.0182, \n",file=filename,append=TRUE) 
  cat("2015,  9.5, 123.38, 0.6780, \n",file=filename,append=TRUE)  
  cat("2016,  9.5, 195.36, 0.8337, \n",file=filename,append=TRUE)   
  cat("2017,  9.5, 185.29, 1.3069, \n",file=filename,append=TRUE)  
  cat("2018,  9.5, 162.85, 0.8025, \n",file=filename,append=TRUE)   
  cat("2019,  9.5, 104.57, 0.6291, \n",file=filename,append=TRUE) 
  cat("2020,  9.5, 115.0, 0.8932, \n",file=filename,append=TRUE) 
  cat("2021,  9.5, 129.1, 0.9424, \n",file=filename,append=TRUE)
  cat("\n\n",file=filename,append=TRUE)
  cat("SIZECOMP, 31, 17, 2,  \n",file=filename,append=TRUE)
  cat("2001,2002,2003,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019  \n",file=filename,append=TRUE)
  cat("5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20.0,  \n",
      file=filename,append=TRUE)
  cat("0,0,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("1,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("4,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("4,4,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("11,15,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("39,30,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("48,42,18,1,0,0,0,0,1,0,0,1,0,0,0,0,1,  \n",file=filename,append=TRUE)
  cat("60,66,55,14,52,31,8,11,15,12,16,35,49,13,13,26,53,  \n",
      file=filename,append=TRUE)
  cat("78,103,70,132,194,124,59,90,67,77,84,87,67,65,66,59,123,  \n",
      file=filename,append=TRUE)
  cat("72,153,106,163,250,97,51,72,108,114,123,48,85,49,52,64,100,  \n",
      file=filename,append=TRUE)
  cat("72,142,86,191,203,85,45,63,112,94,83,57,61,44,58,54,96,  \n",
      file=filename,append=TRUE)
  cat("88,253,76,237,175,77,75,58,107,51,106,45,50,47,68,76,70,  \n",
      file=filename,append=TRUE)
  cat("103,196,53,214,147,71,102,60,72,36,67,32,54,49,97,66,47,  \n",
      file=filename,append=TRUE)
  cat("110,121,36,254,119,63,92,74,64,43,49,35,37,30,57,70,34,  \n",
      file=filename,append=TRUE)
  cat("56,85,11,259,106,40,81,63,46,21,25,21,22,19,51,48,36,  \n",
      file=filename,append=TRUE)
  cat("51,53,4,188,87,34,84,42,29,16,13,16,20,10,43,43,19,  \n",
      file=filename,append=TRUE)
  cat("16,35,1,115,44,21,67,30,14,12,11,8,11,4,15,25,24,  \n",
      file=filename,append=TRUE)
  cat("5,20,0,59,19,23,16,27,8,10,8,2,3,1,4,13,12,  \n",
      file=filename,append=TRUE)
  cat("1,8,0,17,16,19,9,11,2,4,4,2,3,1,2,9,10,  \n",file=filename,append=TRUE)
  cat("6,5,0,16,8,9,7,1,2,1,0,1,2,0,1,2,5,  \n",file=filename,append=TRUE)
  cat("0,0,0,2,1,4,1,1,2,0,1,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("0,2,0,1,3,0,0,0,0,0,0,0,0,0,0,1,0,  \n",file=filename,append=TRUE)
  cat("0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("0,0,3,0,0,0,0,0,0,0,0,0,1,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("3,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("2,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("7,18,2,0,0,0,0,0,0,0,0,1,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("25,35,3,1,0,0,0,0,0,0,0,0,0,1,0,0,0,  \n",file=filename,append=TRUE)
  cat("44,44,7,0,1,0,0,0,1,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("53,51,50,42,77,43,13,39,23,26,34,30,45,30,13,20,83,  \n",
      file=filename,append=TRUE)
  cat("62,64,95,158,238,109,55,110,76,109,118,78,93,73,50,64,138,  \n",
      file=filename,append=TRUE)
  cat("93,115,89,202,275,107,52,94,78,121,139,81,96,58,56,51,149,  \n",
      file=filename,append=TRUE)
  cat("84,122,80,201,241,95,57,93,123,81,124,66,71,55,55,63,144,  \n",
      file=filename,append=TRUE)
  cat("92,173,59,245,257,103,51,84,116,74,136,64,77,43,74,80,119,  \n",
      file=filename,append=TRUE)
  cat("76,147,43,237,198,70,65,68,89,60,75,42,54,56,62,73,76,  \n",
      file=filename,append=TRUE)
  cat("87,196,57,250,154,80,73,66,110,33,62,40,63,39,82,77,81,  \n",
      file=filename,append=TRUE)
  cat("64,153,52,237,152,59,75,76,84,52,74,29,49,65,80,75,71,  \n",
      file=filename,append=TRUE)
  cat("70,166,47,273,138,66,94,81,58,34,62,50,50,50,102,74,81,  \n",
      file=filename,append=TRUE)
  cat("83,109,60,276,118,58,113,63,43,43,49,38,42,36,111,77,86,  \n",
      file=filename,append=TRUE)
  cat("85,96,25,310,117,46,91,47,44,25,50,23,31,30,81,61,44,  \n",
      file=filename,append=TRUE)
  cat("74,65,13,255,80,27,60,67,40,24,37,31,22,13,59,48,42,  \n",
      file=filename,append=TRUE)
  cat("62,44,2,223,70,20,52,41,39,16,33,17,23,5,36,46,25,  \n",
      file=filename,append=TRUE)
  cat("20,36,0,122,34,7,35,40,17,7,12,12,7,10,9,19,15,  \n",
      file=filename,append=TRUE)
  cat("14,24,0,58,16,7,13,13,8,1,3,3,7,3,2,10,11,  \n",file=filename,append=TRUE)
  cat("5,4,0,12,6,1,4,10,2,2,2,1,3,1,1,3,3,  \n",file=filename,append=TRUE)
  cat("0,4,0,2,2,0,0,4,0,0,0,3,2,0,0,2,1,  \n",file=filename,append=TRUE)
  cat("1,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  \n",file=filename,append=TRUE)
  cat("  \n",file=filename,append=TRUE)
} # end of SBMdatatemplate

#' @title SBMpintemplate generates an example of a pin file's format
#' 
#' @description SBMpintemplate is used to illustrate the required format of
#'     the pin file. It may not hold examples of all parameters but it does
#'     provide a template in the correct format that can be edited. 
#'
#' @param rundir the full path to the working directory being used
#' @param infile the name to be given to the pin file
#' @param title text to be placed at the head of the pin file as a label
#' 
#' @seealso{
#'  \link{readSBMpin}
#' }
#'
#' @return nothing but it does generate a file named 'infile' within 'rundir'
#' @export
#'
#' @examples
#' print("wait on data")
SBMpintemplate <- function(rundir,infile,title="A_rock_Lobster") {
  filename <- pathtopath(rundir,infile)
  cat(title,"\n\n",file=filename,append=FALSE)
  cat("PARAMNUM, 27 ,the number of parameters \n",file=filename,append=TRUE)
  cat("LnR0 , 14.1804 , 1 , log init recruitment \n",file=filename,append=TRUE)
  cat("qest , -0.45638 , 1 , log catchability \n",file=filename,append=TRUE)
  cat("selL50 , 2.302585 , 1 , log sel L50 \n",file=filename,append=TRUE)  
  cat("seldelta, 0.4054651 , 1 , log seldelta \n",file=filename,append=TRUE)
  cat("d1995 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d1996 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d1997 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d1998 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d1999 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2000 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2001 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2002 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2003 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2004 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2005 , -0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2006 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2007 , -0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2008 , -0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2009 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2010 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2011 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2012 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2013 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2014 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2015 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2016 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("d2017 , 0.001 , 0 , log rec dev \n",file=filename,append=TRUE) 
  cat("\n\n",file=filename,append=TRUE)
} # end of LBMpintemplate

#' @title readSBMctrl reads the control file and generates working R object
#' 
#' @description readSBMctrl reads in the control file and outputs the ctrl,
#'     constants, and glb objects
#'
#' @param rundir the full path to the working directory being used
#' @param infile the name of the control file found in rundir
#' @param verbose should progress be reported to the console, default=TRUE
#' 
#' @seealso{
#'  \link{SBMctrltemplate}
#' }
#'
#' @return a list containing the ctrl, constants, and glb objects 
#' @export
#'
#' @examples
#' \dontrun{
#'  rundir <- tempdir()
#'  SBMctrltemplate(rundir,infile="ctrlSBM.csv")
#'  setup <- readSBMctrl(rundir,infile="ctrlSBM.csv") 
#' }
readSBMctrl <- function(rundir,infile="controlLBM.csv",verbose=TRUE) {
  # rundir=rundir; infile="ctrlSBM.csv"; verbose=TRUE
  filenames <- dir(rundir)
  if (length(grep(infile,filenames)) != 1)
    stop(cat(infile," not found in ",rundir," \n"))
  filename <- pathtopath(rundir,infile)
  indat <- readLines(filename)   # reads the whole file as character strings
  begin <- grep("START",indat) + 1
  runlabel <- getStr(indat[begin],1)
  datafile <- getStr(indat[begin+1],1)
  pinfile <- getStr(indat[begin+2],1)
  filedata <- TRUE
  if (length(grep(datafile,filenames)) != 1) {
    warning("population data file not found \n")
    filedata <- FALSE
  }
  if (length(grep(pinfile,filenames)) != 1) {
    warning("initial parameter file not found \n")
    filedata <- FALSE
  }
  if (verbose) {
    if (filedata) cat("All required files appear to be present \n")
    else cat("Not all files appear to be present \n")
  }
  pindat <- readSBMpin(rundir=rundir,infile=pinfile)
  recloc <- which(substr(rownames(pindat),1,1) == "d")
  recnames <- rownames(pindat)[recloc]
  recyrs <- as.numeric(substr(recnames,2,5))
  if(verbose) cat("Now reading biological constants from control file \n")
  minsc <-  getsingleNum("minsc",indat) # minimum size class
  cw    <- getsingleNum("cw",indat) # class width
  Nclass <- getsingleNum("Nclass",indat) # number of classes
  midpts <- seq(minsc,minsc+((Nclass-1)*cw),cw)
  FWta <- getsingleNum("FWta",indat)
  FWtb <- getsingleNum("FWtb",indat)
  FM50 <- getsingleNum("FM50",indat)
  Fdelta <- getsingleNum("Fdelta",indat)
  FLinf <- getsingleNum("FLinf",indat)
  FvbK <- getsingleNum("FvbK",indat)
  FCVvb <- getsingleNum("FCVvb",indat)
  MWta <- getsingleNum("MWta",indat)
  MWtb <- getsingleNum("MWtb",indat)
  MM50 <- getsingleNum("MM50",indat)
  Mdelta <- getsingleNum("Mdelta",indat)
  MLinf <- getsingleNum("MLinf",indat)
  MvbK <- getsingleNum("MvbK",indat)
  MCVvb <- getsingleNum("MCVvb",indat)
  # Now fishery constants
  recyr1 <- recyrs[1]
  recyr2 <- recyrs[length(recyrs)]
  natM <- getsingleNum("natM",indat)  
  sigR <- getsingleNum("sigmaR",indat)  
  steep <- getsingleNum("steep",indat)  
  wtsc <- getsingleNum("wtsc",indat)
  sigce <- getsingleNum("sigce",indat) # set = 0 if starting
  lambda <- getsingleNum("lambda",indat)
  # initdepl <- getsingleNum("initdepl",indat)
  # if (initdepl == 1) initdepl <- NULL
  pickL <- grep("omega",indat)
  omega <- getConst(indat[pickL],4,2) # need four values
  #adjbias <- calcadjbias(y1=y1,y2=y2,y3=y3,y4=y4,bmax=bmax,recyrs=recyrs)
  if (verbose) cat("All data appears to have been read successfully \n")
  constants <- c(FWta,FWtb,FM50,Fdelta,FLinf,FvbK,FCVvb,MWta,MWtb,MM50,
                 Mdelta,MLinf,MvbK,MCVvb)
  names(constants) <- c("FWta","FWtb","FM50","Fdelta","FLinf","FvbK","FCVvb",
                        "MWta","MWtb","MM50","Mdelta","MLinf","MvbK","MCVvb")
  id <- c(natM,sigR,steep,wtsc,sigce,lambda)
  names(id) <- c("M","sigmaR","steepness","wtsc","sigce","lambda")
  ctrl <- list(runlabel=runlabel,ctrlfile=infile,datafile=datafile,
               pinfile=pinfile,rundir=rundir,id=id)
  glb <- list(midpts=midpts,Nclass=Nclass,recyrs=recyrs,recdev1=recloc[1],
              natM=natM,sigR=sigR,steep=steep,wtsc=wtsc,sigce=sigce,
              lambda=lambda,#initdepl=initdepl,
              omega=omega,recpar=NULL) #,phase=1)
              #adjbias=adjbias,y1=y1,y2=y2,y3=y3,y4=y4,
  return(list(ctrl=ctrl,constants=constants,glb=glb))
} # end of readSBMctrl


#' @title readSBMdata reads the fishery and sizecomp data for a SBM model
#' 
#' @description readSBMdata is used to read in the fishery data - year, catch,
#'     cpue, and size composition data. The format for the data file can be 
#'     obtained by using the SBMdatatemplate function.
#'
#' @param rundir the complete path to the directory containing the scenario 
#'     files
#' @param infile the name of the csv file containing the data
#'
#' @returns a list containing a matrix of the fishery data, fish, and an array 
#'     of the size-composition data size-class x year x gender, sizecomp.
#' @export
#'
#' @examples
#' print("wait on available datasets") #
readSBMdata <- function(rundir,infile) { 
  #   rundir=rundir; infile="dataSBM.csv";
  filename <- pathtopath(rundir,infile)
  dat <- readLines(filename)
  linenum <- grep("TITLE",dat)
  title <- dat[linenum+1]
  linenum <- grep("FISHERY",dat)
  dimen <- getConst(dat[linenum],nb=2,index=2)
  columns <- getStr(dat[linenum+1],nb=dimen[2],index=1)
  fish <- matrix(0,nrow=dimen[1],ncol=dimen[2]); colnames(fish) <- columns
  numline <- linenum + 2
  for (i in 1:dimen[1]) { # i = 1
    fish[i,] <- getConst(dat[numline],nb=dimen[2],index=1)
    numline <- numline+1
  }
  rownames(fish) <- fish[,1]
  linenum <- grep("SIZECOMP",dat)
  dimen <- getConst(dat[linenum],nb=3,index=2)
  sizecomp <- array(data=0,dim=dimen)
  yrs <- getConst(dat[linenum+1],nb=dimen[2],index=1)
  sizeclass <- getConst(dat[(linenum+2)],nb=dimen[1],index=1)
  if (dimen[3]==1) sexes <- c("X") else sexes <- c("F","M")
  dimnames(sizecomp) <- list(sizeclass,yrs,sexes)
  numline <- linenum + 3
  for (sex in 1:dimen[3]) {
    for (size in 1:dimen[1]) {
      sizecomp[size,,sex] <- getConst(dat[numline],nb=dimen[2],index=1)
      numline <- numline + 1
    }
  }
  return(list(fish=fish,sizecomp=sizecomp))
} # end of readSBMdata

#' @title readLBMpin reads in an initial parameters file
#' 
#' @description readSBMpin reads in the initial parameter values from the 
#'     pin file named in the control file and expected to be found in rundir.
#'     Its format is very specific and the SBMpintemplate function should be 
#'     used to provide an example.
#'
#' @param rundir the full path to the working directory being used
#' @param infile the name of the pin file
#' 
#' @seealso{
#'  \link{SBMpintemplate}
#' }
#'
#' @return a vector of parameters in a fixed order
#' @export
#'
#' @examples
#' print("wait on suitable internal data")
readSBMpin <- function(rundir, infile) { # rundir=rundir; infile="pinSBM.csv"
  filenames <- dir(rundir)
  if (length(grep(infile,filenames)) != 1)
    stop(cat(infile," not found in ",rundir," \n"))
  filename <- pathtopath(rundir,infile)
  indat <- readLines(filename)   # reads the whole file as character strings
  numpar <- getsingleNum("PARAMNUM",indat)
  columns <- c("param","phase")
  rows <- NULL
  param <- matrix(0,nrow=numpar,ncol=length(columns))
  colnames(param) <- columns
  begin <- grep("PARAMNUM",indat)
  for (i in 1:numpar) { # i=1
    tmp <- removeEmpty(unlist(strsplit(indat[begin+i],",")))
    param[i,] <- as.numeric(tmp[2:3])
    rows <- c(rows,tmp[1])
  }
  rownames(param) <- rows
  return(param)
} # end of readLBMpin

#' @title STMvB Generates a Size Transition Matrix for von Bertalanffy growth
#'
#' @description STMvB With the input of the three parameters inside a vector,
#'     and a vector of initial lengths or mid-points of size classes STMvB
#'     generates a square transition matrix with the probabilities of
#'     growing from each initial size into the same or larger sizes.
#'     Negative growth is disallowed so the maximum size class = Linf. All 
#'     columns in the matrix sum to one.The von Bertalanffy growth curve 
#'     implies that the growth increments decline linearly with initial size.
#'     This, in turn implies that when an animal reaches the Linf it will stay
#'     at that size. This uses the Fabens version of the vB curve, which has
#'     different implications to the age-based vB curve. The implementation 
#'     here does not have a plus-group in the final size. When the curve 
#'     implies there would be negative growth the cut-off means that initially
#'     some of the later columns do not sum to 1.0. Rather than lumping all the 
#'     'missing' fish into the largest size class the total for each column is
#'     used to rescale the proportion of fish in each of the size classes into
#'     which it grows.
#'     
#' @param p a vector of three parameters in the following order
#'     Linf the maximum size of all animals, K the growth rate parameter that
#'     determines how quickly each individual approaches the maximum size, and
#'     growCV, the coefficient of variation of the growth such that it 
#'     describes the standard deviation of any Normal distribution used to 
#'     describe the distribution or spread of sizes resulting from a given 
#'     initial size. It does this through the sd = expected mean size (initial 
#'     size + mean increment) multiplied by the growCV (sd = x_bar x growCV). 
#' @param mids a vector of initial lengths which also define the width of
#'     each size class thus, mids from 1 - 60 would imply size classes 1, 2, 3,
#'     ... 60. With the units matching the original measurements.
#'     
#' @return A square matrix with dimension = the length of the mids vector
#' @export
#'
#' @references Haddon, M. (2011) \emph{Modelling and Quantitative Methods in 
#'     Fisheries} 2nd edition. Boca Raton: Chapman & Hall/CRC. 449 pp. 
#'     ISBN: 978-1-58488-561-0
#'     
#' @references Fabens, A.J. (1965) Properties and fitting of the von Bertalanffy 
#'     growth curve. \emph{Growth}, \emph{29}: 265-289.
#' 
#' @examples
#'  pars <- c(Linf=75.0,K=0.275,CVvb=0.08)
#'  lens <- seq(3,73,5)
#'  M <- 0.25; surv <- exp(-M)
#'  G <- STMvB(p=pars,mids=lens)
#'  print(round(G,4))
#'  nyr <- 20
#'  Nt <- matrix(0,nrow=15,ncol=nyr,dimnames=list(lens,2000:(1999+nyr)))
#'  rec <- numeric(15); rec[1] <- 100
#'  Nt[,1] <- rec
#'  for (i in 1:(nyr-1)) Nt[,i+1] <- G %*% (surv * Nt[,i]) + rec
#'  round(Nt,1)
STMvB <- function(p,mids) { #    # p <- pars; mids <- sizes
  n <- length(mids)
  G <- matrix(0,nrow=n,ncol=n,dimnames=list(mids,mids))
  cw <- mids[2]-mids[1]
  Linf <- p[1]  
  K <- p[2]  
  CV <- p[3]
  inc <- (Linf - mids) * (1 - exp(-K))
  meanL <- mids + inc
  sd <- meanL * CV 
  for (j in 1:n) {
    for (i in 1:n) {
      Prob <- (1-pnorm(mids[i],meanL[j],sd[j],FALSE))
      if (i < j)  { G[i,j] <- 0.0 }
      if (i == j) { G[i,j] <- Prob }
      if (i > j)  { G[i,j] <- Prob - (1-pnorm(mids[i-1],meanL[j],
                                              sd[j],FALSE)) }
    }
  }
  csum <- colSums(G)
  for (i in 1:n) G[,i] <- G[,i]*(1/csum[i])
  class(G) <- "STM"
  return(G)
} # end of STMvB


#' @title vbfabens calculates predicted mean growth increments for vB Fabens
#' 
#' @description vbfabens calculate the predicted mean growth increment 
#'     expected from the Fabens version of the von Bertalanffy growth curve.
#'     This is designed to be fitted to tagging data which is used to estimate
#'     the mean increment of different sized animals. Assuming vB growth
#'     implies a linear decrease in the growth increment with increasing size.
#'
#' @param p a vector of two vB parameters Linfinity, the maximum size, and K 
#'     the rate of growth towards the maximum size.
#' @param mids the center of each size class for which the mean expected 
#'     groeth increment will be calculated.
#'
#' @returns a vector of growth increments the samelength as the vector of mid-
#'     points.
#' @export
#'
#' @examples
#'  mids <- seq(3,75,5)
#'  p <- c(75,0.26)
#'  Linc <- vbfabens(p,mids)
#'  print(cbind(mids,Linc))
vbfabens <- function(p,mids) {
  const <- (1 - exp(-p[2]))
  predG <- (p[1] - mids) * const
  return(predG)
} # end of vbfabens




