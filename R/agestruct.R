
#' @title bhsim - calculates number of bh from input spawning biomass
#'
#' @description bhsim - calculates number of recruits from input spawning
#'     biomass assumes steep, R0, B0, and sigmaR are available. It
#'     also assumes the use of the Beverton - Holt stock recruitment curve. 
#'     Only used in simdynF
#'
#' @param inSB input spawning biomass in tonnes.
#' @param R0 unfished recruitment levels
#' @param B0 unfished spawning biomass
#' @param steep the steepness of the Beverton-Holt SR relationship
#' @param sigmaR recruitment variability
#' 
#' @return a scaler representing the number of recruits
#' @export
#' 
#' @examples
#' \dontrun{
#'  library(fmr)
#'  data(const)
#'  globals <- initiateglobals(const)
#'  pop <- makestock(glb=globals)
#'  bh(2500,pop$R0,pop$B0,sigmaR=1e-07,globals$steep) # = 1763871
#' } 
bhsim <- function(inSB,R0,B0,steep,sigmaR) {
  epsilon <- exp(rnorm(1,mean=0,sd=sigmaR) - (sigmaR * sigmaR)/2)
  recs <- ((4*steep*R0*inSB)/(((1-steep)*B0)+(5*steep-1)*inSB)) * epsilon
  return(recs)
} # end of bhsim

#' @title calcsel estimates selectivity
#'
#' @param ageorlen a vector of ages or lengths against which to estimate the 
#'     selectivity
#' @param fishprop the fishery properties for a particular fleet from the
#'     fishery list object out of getconstants 
#' @param typeselect what sort of selectivity curve to use. Currently, the only
#'     options are 'logistic' and 'domed' selectivity
#' @param zerocomps a vector of ages or sizes that are never selected. Identify
#'     their indices in the outsel vector here. eg. if age 0 is never selected
#'     then set zerocomps = 1, if ages 1-2 then zerocomps = c(1,2).
#'     
#'
#' @returns a vector of selectivity for each ageorlen
#' @export
#'
#' @examples
#' ages <- 0:25
#' pars <- c(0, 3, 1.0) # catchability + logistic parameters q, sel50, deltasel
#' cbind(ages,calcsel(ageorlen=ages,fishprop=pars,typeselect="logistic"))
#' pars <- c(0,5,9,16,23,-5,-2)
#' cbind(ages,calcsel(ageorlen=ages,fishprop=pars,typeselect="domed",
#'       zerocomp=1))
calcsel <- function(ageorlen,fishprop,typeselect,zerocomps=NULL) {
  propfish <- fishprop[2:length(fishprop)] # pick only selectivity terms
  npar <- ifelse(typeselect == "logistic",2,6)
  pickp <- propfish[1:npar]
  outsel <- switch(typeselect,
                   logistic = 1/(1+exp(-(ageorlen - pickp[1])/pickp[2])),
                   domed = domedsel(pickp,ageorlen),
                   stop("Unknown selectivity chosen. "))
  if (length(zerocomps) > 0) outsel[1:zerocomps] <- 0.0
  return(outsel)
} # end of calcsel

#' @title domedsel calculates domed selectivity curves
#' 
#' @description domedsel uses 6 parameters and a set of mean size or age classes 
#'     to calculate a domed selectivity curve with a maximum of 1.0 (rescaling 
#'     can be done outside the function), but has parameters for the selectivity 
#'     of the initial and final size/age classes. There is an ascending limb and 
#'     a descending limb with the potential of a plateau in between. The six 
#'     parameters are:
#'     
#'     1. the age/size where selectivity first becomes 1.0, peak1
#'     
#'     2. the size/age where selectivity first begins to decline, peak2
#'     
#'     3. the steepness of the ascending limb, asc ln(width)
#'     
#'     4. the steepness of the descending limb, dsc  ln(width)
#'     
#'     5. the selectivity of the first age/size class, selmin, and 
#'     
#'     6. the selectivity of the last age/size class, selmax 
#'     
#'     The selectivity of the first and last composition classes, selmin and 
#'     selmax, are the inverse logit transformation of the value used in the 
#'     calculations.
#'     The descending limb of any dome shaped selectivity curves imply that the 
#'     fishing gear used is unable to collect all representatives of the larger 
#'     or older classes. The predicted numbers of smaller or younger animals, 
#'     that are only partially selected, are inflated because of the partial 
#'     selection. If any larger or older animals are, in fact, caught, then the 
#'     same inflation can happen to those animals as a result of the partial 
#'     selection implied by the dome shape. Small and young animals weight 
#'     very little, the same cannot be said for the larger or older animals. 
#'     Some people refer to the extra biomass this phenomenon can imply as 
#'     'ghost biomass', even though it might be real. Whatever the case, when 
#'     using dome shaped selectivity it is best to be aware of this issue and 
#'     to be cautious about how this is interpreted. The 20* terms in the J1 
#'     and J2 factors are required to force the joins to be as effective as
#'     required (see Methot and Wetzel). 
#'
#' @param p a vector of six parameters.
#' @param L a vector of the mean of nL age/size classes
#'
#' @return a vector of selectivities
#' @export
#' 
#' @references Methot, R.D. and C.R, Wetzel (2013) Stock synthesis: A biological 
#'     and statistical framework for fish stock assessment and fishery management. 
#'     Supplementary material, Appendix A. Equs A1.30 onwards. 
#'     \emph{Fisheries Research} 142:86-99.
#'     
#'     Hurtado-Ferro, F., Punt, A.E., and K.T. Hill (2014) Use of multiple 
#'     selectivity patterns as a proxy for spatial structure. 
#'     \emph{Fisheries Research} 158:102-115.
#'
#' @examples
#'   L <- seq(1,60,1)
#'   p <- c(25,32,16,33,-5,-2)
#'   sel <- domedsel(p,L)
#'   plot(L,sel,type="l",xlab="Age",ylab="Selectivity",lwd=2)
domedsel <- function(p,L) {
  nL <- length(L)
  J1 <- 1/(1 + exp(-20*((L - p[1])/(1 + abs(L - p[1])))))
  J2 <- 1/(1 + exp(-20*((L - p[2])/(1 + abs(L - p[2])))))   
  comp1 <- 1/(1 + exp(-p[5])) # forced to be 0 - 1
  comp2 <- exp((-(L - p[1])^2)/p[3])
  comp3 <- exp((-(L[1] - p[1])^2)/p[3])
  asc <- comp1 + (1 - comp1) * ((comp2 - comp3)/(1 - comp3))
  comp4 <- 1/(1 + exp(-p[6]))
  comp5 <- exp((-(L - p[2])^2)/p[4])
  comp6 <- exp((-(L[nL] - p[2])^2)/p[4])
  dsc <- 1 + (comp4 - 1) * ((comp5 - 1)/(comp6 - 1))
  sel <- (asc * (1 - J1)) + J1 * (1 - J2 + dsc * J2)
  sel <- sel/max(sel) # to ensure a maximum = 1.0
  return(sel)
} # end of domedsel

#' @title readASPMdata - reads the input constants that condition the model
#'
#' @description readASPMdata - reads the input constants that condition the
#'     model. The constants include all those produced by the function
#'     template2F1S: M, surv, Linf, K, t0, growCV, Wta, Wtb, steep,
#'     age50M, deltaM, B0, sigmaR, etc. The data are read in from 'infile' and
#'     must contain sections entitled: BIOLOGY, FISHERY, STRUCTURE, and
#'     HISTORICALCATCH
#' @param infile - defaults to 'constants.csv'. Needs to be a csv file
#' 
#' @return a list of vectors 'biology' and 'fish' plus a lists of 'fishbiol',
#'     'fishery', and 'glb'
#' @export
#' 
#' @examples
#' \dontrun{
#'   rundir <- tempdir()
#'   template2F1S(rundir,filename="Flt2Stock1.csv")
#'   const2 <- readASPMdata(infile=paste0(rundir,"//Flt2Stock1.csv"))
#'   str(const2)
#' }
readASPMdata <- function(infile="constants.csv") { #  infile=filen
  #  infile=pathtopath(rundir,"F21S.csv")
  # infile=pathtopath(rundir,"test2F1S.csv")
  datain <- readLines(con = infile)
  # structure------
  randseed <- getsingleNum("randseed",datain)
  set.seed(randseed)
  #lens <- getvect("LFstruct",datain,n=3)
  age <- getvect("agestruct",datain,n=3)
 # lengths <- seq(lens[1],lens[2],lens[3])
  ages <- seq(age[1],age[2],age[3])  
  nfleet <- getsingleNum("nfleet",datain)
  linenum <- grep("fleets",datain) # get fleetnames
  tmp <- removeEmpty(unlist(strsplit(datain[linenum],",",fixed=TRUE)))
  fleets <- tmp[2:length(tmp)]
  linenum <- grep("selecttype",datain) # get selectivity type by fleet
  tmp <- removeEmpty(unlist(strsplit(datain[linenum],",",fixed=TRUE)))
  seltype <- tmp[2:(nfleet+1)]   
  # fishery------
  selfleet <- paste0(fleets,"S")
  fishery <- makelist(selfleet)
  for (i in 1:nfleet) {  # i = 1
    if (seltype[i] == "logistic") {
      label <- c("qc","sel50","deltasel")      
      fishery[[i]] <- getvect(selfleet[i],datain,n=length(label))
      names(fishery[[i]]) <- c("qc","sel50","deltasel")     
    }
    if (seltype[i] == "domed") {
      fishery[[i]] <- getvect(selfleet[i],datain,n=7)
      names(fishery[[i]]) <- c("qc","peak1","peak2","asc","dsc","selmin",
                               "selmax") 
    }
  }
  # if (nsex == 1) {
  #   sexes <- c("c")
  # } else {
  #   sexes <- c("f","m")
  # }
  # if (nregion == 1) {
  #   sexreg <- paste(regions,sexes,sep="_")
  #   fltreg <- paste(fleets,regions,sep="_")
  # } else {
  #   sexreg <- c(paste(regions[1],sexes,sep="_"),
  #               paste(regions[2],sexes,sep="_"))
  #   fltreg <- c(paste(fleets,regions[1],sep="_"),
  #               paste(fleets,regions[2],sep="_"))
  # }
  # biology------
  rows <- c("Linf","K","t0","growCV","Wta","Wtb","age50M","deltaM")
  nrows <- length(rows)
#  ncols <- length(sexreg)
  ncols=1
  biology <- matrix(0,nrow=nrows,ncol=ncols,dimnames=list(rows,"mixed"))
  biology["Linf",] <- getvect("Linf",datain,ncols)
  biology["K",] <- getvect("K",datain,ncols)
  biology["t0",] <- getvect("t0",datain,ncols)
  biology["growCV",] <- getvect("growCV",datain,ncols)
  biology["Wta",] <- getvect("Wta",datain,ncols)
  biology["Wtb",] <- getvect("Wtb",datain,ncols)
  biology["age50M",] <- getvect("Age50M",datain,ncols)
  biology["deltaM",] <- getvect("deltaM",datain,ncols)
  M <- getsingleNum("M",datain)
  steep <- getsingleNum("steep",datain)
  setup <- grep("HISTORICALCATCH", datain)  # histcatch------
  ncat <- getsingleNum("HISTORICALCATCH",datain)
  columns <- c("year",fleets,paste0(fleets,"CE"))
  numcol <- length(columns)
  fish <- matrix(0,nrow=ncat,ncol=numcol)
  for (i in 1:ncat) 
    fish[i,] <- getConst(datain[(setup + i)],nb=numcol,index=1)
  colnames(fish) <- columns
  rownames(fish) <- fish[,"year"]
  for (i in 1:numcol) fish[which(fish[,i] == 0),i] <- NA
  cenames <- paste0(fleets,"CE")
  outrmse <- makelist(cenames)
  for (flt in 1:nfleet) 
    outrmse[[cenames[flt]]] <- getrmse(fish,invar=cenames[flt])
  sigCE <- sapply(outrmse,"[[","rmse")
  
    
  nyrs <- ncat
  initdepl <- getvect("initdepl",datain,1)
  startyr <- fish[1,"year"]
  endyr <- fish[ncat,"year"]
  glb <- list(ages,length(ages),nyrs,startyr,endyr,fleets,nfleet,
              seltype,initdepl,M,steep,initdepl,sigCE)
  names(glb) <- c("ages","nages","nyrs","startyr","endyr","fleets","nfleet",
                  "selecttype","initdepl","M","steep","initdepl","sigCE")
  ans <- list(biology,fishery,glb,fish,outrmse)
  names(ans) <- c("biology","fishery","glb","fish","outrmse")
  return(ans)
} # end of readASPMdata

#' @title getconstants - reads the input constants that condition the model
#'
#' @description getconstants - reads the input constants that condition the
#'     model. The constants include all those produced by the function
#'     constantfileTemplate: M, surv, Linf, K, t0, growCV, Wta, Wtb, steep,
#'     age50M, deltaM, B0, sigmaR, etc. The data are read in from 'infile' and
#'     must contain sections entitled: BIOLOGY, FISHERY, STRUCTURE, and
#'     HISTORICALCATCH
#' @param infile - defaults to 'constants.csv'. Needs to be a csv file
#' 
#' @return a list of two vectors 'biology' and 'fishery' plus a list of two
#'     vectors of lengths and ages
#' @export
#' 
#' @examples
#' require(fmr)
#' data(const)
#' str(const)
#' # infile=pathtopath(datraw,"F2-A1-S1.csv")
getconstants <- function(infile="constants.csv") { #  infile=filen
  datain <- readLines(con = infile)
  # structure------
  randseed <- getsingleNum("randseed",datain)
  set.seed(randseed)
  nregion <- getsingleNum("nregion",datain)
  linenum <- grep("regname",datain) # get fleetnames
  tmp <- removeEmpty(unlist(strsplit(datain[linenum],",",fixed=TRUE))) 
  regions <- tmp[2:length(tmp)]
  nsex <- getsingleNum("nsex",datain)
  lens <- getvect("LFstruct",datain,n=3)
  age <- getvect("agestruct",datain,n=3)
  lengths <- seq(lens[1],lens[2],lens[3])
  ages <- seq(age[1],age[2],age[3])  
  nfleet <- getsingleNum("fleets",datain)
  linenum <- grep("fleetname",datain) # get fleetnames
  tmp <- removeEmpty(unlist(strsplit(datain[linenum],",",fixed=TRUE)))
  fleets <- tmp[2:length(tmp)]
  linenum <- grep("selecttype",datain) # get selectivity type by fleet
  tmp <- removeEmpty(unlist(strsplit(datain[linenum],",",fixed=TRUE)))
  seltype <- tmp[2:(nfleet+1)]   
  # fishery------
  selfleet <- paste0(fleets,"S")
  fishery <- makelist(selfleet)
  for (i in 1:nfleet) {  # i = 1
    if (seltype[i] == "logistic") {
      label <- c("qc","sel50","deltasel")      
      fishery[[i]] <- getvect(selfleet[i],datain,n=length(label))
      names(fishery[[i]]) <- c("qc","sel50","deltasel")     
    }
    if (seltype[i] == "domed") {
      fishery[[i]] <- getvect(selfleet[i],datain,n=7)
      names(fishery[[i]]) <- c("qc","peak1","peak2","asc","dsc","selmin",
                               "selmax") 
    }
  }
  if (nsex == 1) {
    sexes <- c("c")
  } else {
    sexes <- c("f","m")
  }
  if (nregion == 1) {
    sexreg <- paste(regions,sexes,sep="_")
    fltreg <- paste(fleets,regions,sep="_")
  } else {
    sexreg <- c(paste(regions[1],sexes,sep="_"),
                paste(regions[2],sexes,sep="_"))
    fltreg <- c(paste(fleets,regions[1],sep="_"),
                paste(fleets,regions[2],sep="_"))
  }
  # biology------
  rows <- c("M","surv","Linf","K","t0","growCV","Wta","Wtb","age50M","deltaM")
  nrows <- length(rows)
  ncols <- length(sexreg)
  biology <- matrix(0,nrow=nrows,ncol=ncols,dimnames=list(rows,sexreg))
  biology["M",] <- getvect(varname="natM",intxt=datain,n=ncols)
  biology["surv",] <- exp(-biology["M",])
  biology["Linf",] <- getvect("Linf",datain,ncols)
  biology["K",] <- getvect("K",datain,ncols)
  biology["t0",] <- getvect("t0",datain,ncols)
  biology["growCV",] <- getvect("growCV",datain,ncols)
  biology["Wta",] <- getvect("Wta",datain,ncols)
  biology["Wtb",] <- getvect("Wtb",datain,ncols)
  biology["age50M",] <- getvect("Age50M",datain,ncols)
  biology["deltaM",] <- getvect("deltaM",datain,ncols)
  # rec
  B0    <- getsingleNum("B0",datain)
  sigmaR <- getsingleNum("sigmaR",datain)
  origsigR <- sigmaR
  steep <- getsingleNum("steep",datain)
  rec=c(B0=B0,sigmaR=sigmaR,origsigR=origsigR,steep=steep)
  splitR <- getvect("R0split",datain,nregion)
  setup <- grep("HISTORICALCATCH", datain)  # histcatch------
  ncat <- getsingleNum("HISTORICALCATCH",datain)
  columns <- c(fleets,"year")
  numcol <- length(columns)
  histcat <- matrix(0,nrow=ncat,ncol=numcol)
  for (i in 1:ncat) 
    histcat[i,] <- getConst(datain[(setup + i)],nb=numcol,index=1)
  colnames(histcat) <- columns
  rownames(histcat) <- histcat[,"year"]
  histcat[which(histcat == 0)] <- NA
  nyrs <- ncat
  # Rsplit <- rep(1,nyrs)   # this allows for proces error in the recruitment
  # if (nregion > 1) {      # split between regions
  #   Rsplit <- rnorm(nyrs,mean=splitR[1],sd=splitR[2])
  #   Rsplit[1] <- splitR[1]
  # }
  sigmaCE <- getvect("sigmaCE",datain,nregion)
  initdepl <- getvect("initdepl",datain,nregion)
  startyr <- histcat[1,"year"]
  endyr <- histcat[ncat,"year"]
  fishbiol <- list(rec=rec,Rsplit=splitR,initdepl=initdepl,sigmaCE=sigmaCE)
  struct <- list(lengths,length(lengths),ages,length(ages),nyrs,startyr,endyr,
                 fleets,nfleet,seltype,nregion,regions,nsex,sexes,sexreg,fltreg)
  names(struct) <- c("sizes","nsizes","ages","nages","nyrs","startyr",
                     "endyr","fleets","nfleet","selecttype","nregion","regions",
                     "nsex","sexes","sexreg","fltreg")
  ans <- list(biology,fishbiol,fishery,struct,histcat,randseed)
  names(ans) <- c("biology","fishbiol","fishery","structure","histC",
                  "randseed")
  return(ans)
} # end of getconstants


#' @title getvect extracts 'n' numbers from an identified line of text
#'
#' @description getvect parses a line of text and extracts 'n' pieces of
#'     text as numbers. The specific line within the intxt, is identified by
#'     the text at the start of the line found using varname
#'
#' @param varname the name of the vector of numbers to be parsed
#' @param intxt the text data file read in using readLines
#' @param n the number of values to extract from the line
#'
#' @returns a vector of n values
#' @export
#'
#' @examples
#' txt1 <- "trawl,   1.4E-05,    2,    0.75,"  
#' txt2 <- "autoline,9.0E-06,    5,    1.25,"
#' tmp <- rbind(txt1,txt2)
#' print(getvect("autoline",tmp,3))
#' print(getvect("trawl",tmp,3))
getvect <- function(varname,intxt,n) { # varname=fleets[i];intxt=datain;n=length(label)
  begin <- grep(varname,intxt)
  nvar <- length(begin)
  if (nvar > 0) {
    if (nvar > 1) {
      tmp <- NULL
      numC <- nchar(varname)      
      for (j in 1:nvar) { # j=1
        start <- substr(intxt[begin[j]],1,numC)
        if (start == varname) tmp <- begin[j]
      }
      begin <- tmp
    }   
    return(as.vector(getConst(intxt[begin],nb=n,index=2),mode="numeric"))
  } else {
    return(NULL)
  }
} # end of getvect

#' @title initiateglobals - calculates fishery constants as vectors
#'
#' @description initiateglobals - calculates fishery constants as vectors,
#'     including length-at-age, weight-at-age, selectivity-a-age, length-at
#'     50% selection, maturity-at-age, and a growth transition matrix. It now
#'     includes sex, fleets, and regions as partitions of the model.
#'
#' @param consts the output from 'getconstants'
#' @param agesorsizes is selectivity by 'ages' or 'sizes', default = 'ages'
#' @param selabove are there ages or sizes that are never selected. Identify
#'     the upper index in the outsel vector here. eg. if age 0 is never selected
#'     then set zerocomps = 1, if ages 0-3 then zerocomps = 4
#' 
#' @return a list containing LaA, WaA, Sel, MatA, growtran, and L50, and more
#' @export
#' 
#' @examples
#' \dontrun{
#' data(simconst2)
#' globals <- initiateglobals(simconst2,agesorsizes="ages",selabove=1)
#' str(globals,max.level=1)
#' }
initiateglobals <- function(consts,agesorsizes,selabove=0) { 
  #  consts=simconst2; agesorsizes="ages";selabove=1
  biology <- consts$biology
  fishbiol <- consts$fishbiol
  fishery <- consts$fishery
  struct <- consts$structure
  histC <- consts$histC
  randseed <- consts$randseed
  sizes <- struct$sizes
  nsizes <- struct$nsizes  
  cw <- sizes[2] - sizes[1]
  ages <- struct$ages
  nages <- struct$nages
  maxage <- max(ages)
  nyrs <- struct$nyrs
  startyr <- struct$startyr
  endyr <- struct$endyr
  nfleet <- struct$nfleet
  fleets <- struct$fleets
  seltype <- struct$selecttype
  nregion <- struct$nregion
  regions <- struct$regions
  nsex <- struct$nsex
  sexes <- struct$sexes
  sexreg <- struct$sexreg
  nsexreg <- length(sexreg)
  fltreg <- struct$fltreg
  nfltreg <- length(fltreg)
  LaA <- matrix(0,nrow=nages,ncol=nsexreg,dimnames=list(ages,sexreg))
  WaA <- matrix(0,nrow=nages,ncol=nsexreg,dimnames=list(ages,sexreg))
  MatA <- matrix(0,nrow=nages,ncol=nsexreg,dimnames=list(ages,sexreg))
  WaL <- matrix(0,nrow=nsizes,ncol=nsexreg,dimnames=list(sizes,sexreg))
  growtran <- array(data=0,dim=c(nsizes,nages,nsex,nregion),
                    dimnames=list(sizes,ages,sexes,regions)) 
  if (agesorsizes == "ages") {
    Sel <- matrix(0,nrow=nages,ncol=nfleet,dimnames=list(ages,fleets))
  } else {
    Sel <- matrix(0,nrow=nsizes,ncol=nfleet,dimnames=list(sizes,fleets))
  }
  columns <- c("meanL","sdL","L95","U95")
  numcol <- length(columns)
  gpar <- as.matrix(biology[c("Linf","K","t0","growCV"),])
  grow <- array(data=0,dim=c(nages,numcol,nsex,nregion),
                dimnames=list(ages,columns,sexes,regions))
  for (i in 1:nsex) {
    for (j in 1:nregion) {
      g <- ifelse(j == 1,i,i+j)
      grow[,"meanL",i,j] <- gpar[1,g] * (1 - exp(-gpar[2,g]*(ages - gpar[3,g])))
      grow[,"sdL",i,j] <- gpar[4,g]*grow[,"meanL",i,j]
      grow[,"L95",i,j] <- qnorm(0.025,mean=grow[,"meanL",i,j],
                                sd=grow[,"sdL",i,j])
      grow[,"U95",i,j] <- qnorm(0.975,mean=grow[,"meanL",i,j],
                                sd=grow[,"sdL",i,j])
      LaA[,g] <- grow[,"meanL",i,j]
      WaA[,g] <- biology["Wta",g] * LaA[,g] ^ biology["Wtb",g]
      WaL[,g] <- biology["Wta",g] * (sizes+(cw/2.0)) ^ biology["Wtb",g]
      MatA[,g] <- 1/(1 + exp(-(ages -biology["age50M",g])/biology["deltaM",g]))
      sdage <- grow[,"sdL",i,j]
      for (age in 0:maxage) {
        usea <- age+1
        x = dnorm(sizes,mean=LaA[usea,g],sd=sdage[usea])
        growtran[,usea,i,j] <- x/sum(x,na.rm=TRUE) 
      }
    }
  }
  for (j in 1:nfleet) { # j = 2 # fleet includes fleet and sex and area
    propfish <- fishery[[j]]
    if (agesorsizes == "ages") {
       Sel[,j] <- calcsel(ageorlen=ages,fishprop=propfish,
                          typeselect=seltype[j],zerocomps=selabove)
       } else {
       Sel[,j] <- calcsel(ageorlen=sizes,fishprop=propfish,
                            typeselect=seltype[j],zerocomps=selabove)  
    }
  }   
  # growtran <- makeSTM(pars,glb$sizes,funL=vBfabens)
  qc <- unlist(lapply(fishery,"[[","qc"))
  ans <- list(ages,maxage,nages,sizes,nsizes,LaA,WaA,Sel,MatA,growtran,
              grow,biology,fishbiol,fishbiol$rec,fishbiol$initdepl,WaL,gpar,
              qc,nyrs,startyr:endyr,nfleet,fleets,seltype,nregion,regions,histC,
              nsex,sexes,sexreg,nsexreg,fltreg,nfltreg,biology["M",1],randseed)
  names(ans) <- c("ages","maxage","nages","sizes","nsizes","LaA","WaA","Sel",
                  "MatA","growtran","grow","biology","fishbiol","rec",
                  "initdepl","wal","gpar","qc","nyrs","years","nfleet",
                  "fleets","seltype","nregion","regions","histC",
                  "nsex","sexes","sexreg","nsexreg","fltreg","nfltreg","natM",
                  "randseed")
  return(ans)
}  # end of initiateglobals

#' @title makedataset generates a csv file from siumulated data
#' 
#' @description makedataset generates a csv file from simulated data produced
#'     by getconstants, initiateglobals, and simdynF.
#'
#' @param const the output object from getconstants
#' @param glb the output object from initiateglobals
#' @param out the output object from simdynF
#' @param filename the name of the datafile, MUST be a '.csv' file
#' @param rundir the directory path to the subdirectory in which all analysis 
#'     data and results are stored, default = ''
#' @param datatitle line 1 of data file decribing its contents, default = ''
#' @param catchSD standard deviation of normal variation imposed on reported
#'     catches. default = 1e-08, effectively no variation
#' @param cpueSD standard deviation of normal variation imposed on reported
#'     CPUE. default = 1e-08, effectively no variation
#'     
#' @seealso{
#'    \link{getconstants}, \link{initiateglobals}. \link{simdynF}
#' }     
#'
#' @returns nothing but it does write a data file to rundir/filename
#' @export
#'
#' @examples
#' \dontrun{
#'  const2 <- getconstants(pathtopath(datraw,"F2-A1-S1.csv"))
#'  glb <- initiateglobals(consts=const2,selabove=1)
#'  stock <- makestock(glb=glb)
#'  out <- simdynF(glb=glb,stk=stock,reps=5,full=TRUE)
#'  catchSD=1e-08
#'  cpueSD=1e-08
#'  sigR=1e-08
#'  filename="F21S.csv"
#'  datatitle <- "A two-fleet one stock/region fishery data set"
#'  makedataset(const2,glb,out,filename,rundir='',datatitle="test",
#'              catchSD=1e-08,cpueSD=1e-08,sigR=1e-08)
#'  # now open 'F21S.csv
#' }
makedataset <- function(const,glb,out,filename,rundir='',datatitle="",
                        catchSD=1e-08,cpueSD=1e-08) {
  
  # const=const2;glb=glb;out=out;stock=stock;filename="F2S1new.csv";
  # rundir=rundir;datatitle="";
  # catchSD=1e-08;cpueSD=1e-08
  
  histC <- samplefishery(out,glb,
                         errors=c(catchSD=catchSD,cpueSD=cpueSD))
  writedatafile(const,glb,histC,rundir,filename,datatitle)
} # end of makedataset


#' @title makeprops generates a matrix of biological and selectivity properties
#' 
#' @description makeprops is a wrapper function that uses a constants file 
#'     to generate a matrix of the constant biological properties for a fishery
#'     including length-at-age, weight-at-age, maturity-at-age, and the 
#'     selectivity of each fleet
#'
#' @param const a structured data file in the format defined by template2F1S
#' @param selabove a vector of ages or sizes that are never selected. Identify 
#'     their indices in the selecitivyt vector here. eg. if age 0 is never 
#'     selected then set zerocomps = 1, if ages 1-2 then zerocomps = c(1,2).
#'     
#' @seealso{
#'   \link{template2F1S}, \link{calcsel},
#' }      
#'
#' @returns a matrix of biological and selectivity propoerties
#' @export
#'
#' @examples
#' \dontrun{
#'   rundir <- tempdir()
#'   dirExists(rundir,verbose=TRUE)
#'   template2F1S(rundir,filename="test2F1S.csv")
#'   const2 <- getconstants(pathtopath(rundir,"test2F1S.csv"))
#'   fish <- const2$histC
#'   glb <- const2$structure
#'   props <- makeprops(const2,selabove=1)
#'   print(props)
#' }
makeprops <- function(const,selabove=1) {  # const = const2; selabove=1
  biol <- const$biology
  glb <- const$glb
  fishery <- const$fishery
  ages <- glb$ages
  nfleet <- glb$nfleet
  fleets <- glb$fleets
  columns <- c("age","laa","waa","maa",fleets)
  props <- as.data.frame(matrix(0,nrow=glb$nages,ncol=length(columns),
                                dimnames=list(ages,columns)))
  props[,"age"] <- ages 
  vbpars <- c(linf=biol["Linf",],K=biol["K",],t0=biol["t0",])
  props[,"laa"] <- vB(vbpars,glb$ages)
  props[,"waa"] <- biol["Wta",] * (props[,"laa"] ^ biol["Wtb",])
  props[,"maa"] <- logist(inL50=biol["age50M",],delta=biol["deltaM",],
                          depend=ages)
  for (flt in 1:nfleet) {
    propfish <- fishery[[flt]]
    props[,fleets[flt]] <- calcsel(ageorlen=ages,fishprop=propfish,
                                   typeselect=glb$selecttype[flt],
                                   zerocomps=selabove)
  }
  return(props=props)
} # end of makeprops

#' @title makesimstock sets up the containers for the stock dynamics
#'
#' @description makestock sets up the containers for the stock dynamics. The
#'     objective is to have an object that can be passed to substantive
#'     functions so the dynamics can be followed easily without using the
#'     global environment and remembering all the particular variables.
#'
#' @param glb the constants obtained from constants.csv including those
#'     relating to the biology, the fishery, and the structural details.
#'
#' @return a list containg the more important placeholders for the stock and
#'     fishery dynamics.
#' @export
#'
#' @examples
#' \dontrun{
#'  library(fmr)
#'  data(const)
#'  globals <- initiateglobals(const)
#'  stock <- makesimstock(globals)
#'  str(stock)
#' }
makesimstock <- function(glb) {  # glb=glb
  ages <- glb$ages
  years <- glb$years
  sizes <- glb$sizes
  nages <- glb$nages
  maxage <- max(glb$ages)
  natM <- glb$natM
  rec <- glb$fishbiol$rec
  nyrs <- length(glb$years)
  yearlab <- c((years[1]-1),years)
  nsizes <- glb$nsizes
  nfleet <- glb$nfleet
  fleets <- glb$fleets
  nregion <- glb$nregion
  regions <- glb$regions
  sexreg <- glb$sexreg
  fltreg <- glb$fltreg
  Sel <- glb$Sel[,1:nfleet]
  
  NaA <- matrix(0,nrow=nages,ncol=(nyrs+1),  # currently A1-S1
                dimnames=list(seq(0,maxage,1),yearlab))
  Bsp <- matrix(0,nrow=(nyrs+1),ncol=1,dimnames=list(yearlab,"SpawnB"))
  columns <- paste0(fleets,"_ExB")
  ExpB <- matrix(0,nrow=(nyrs+1),ncol=nfleet,dimnames=list(yearlab,columns))
  columns <- paste0(fleets,"_catchB")
  catchB <- matrix(0,nrow=(nyrs+1),ncol=nfleet,dimnames=list(yearlab,columns))
  catchB[1,] <- NA
  catchB[2:(nyrs+1),] <- glb$histC[,1:nfleet]
  Harvest <- matrix(0,nrow=(nyrs+1),ncol=1,dimnames=list(yearlab,"Harvest"))
  CatchN <- matrix(0,nrow=nages,ncol=(nyrs+1),
                   dimnames=list(seq(0,maxage,1),yearlab))
  popsizeD <- matrix(0,nrow=nsizes,ncol=(nyrs+1),
                     dimnames=list(sizes,yearlab))
  columns <- paste0(fleets,"_CE")   
  cpue <- matrix(0,nrow=(nyrs+1),ncol=nfleet,dimnames=list(yearlab,columns))
  # Now calculate the unfished state
  ans <- unfished(glb)   # Nt, R0, A0,SpB0,expB0,N0mat
  NaA[,1] <- ans$Nt
  unfishedNaA <- ans$Nt
  A0 <- ans$A0
  R0 <- ans$R0
  B0 <- ans$B0
  ExB0 <- ans$ExB0
  Bsp[1] <- B0
  #ExB0 <- ans$expB0
  if (nfleet == 1) {
    ExpB[1] <- ExB0
  } else {
    ExpB[1,] <- ExB0
  }
  unfishedNaL <- glb$growtran[,,1,1] %*% NaA[,1]
  popsizeD[,1] <- unfishedNaL
  stk <- list(NaA=NaA,Bsp=Bsp,ExpB=ExpB,catchB=catchB,Harvest=Harvest,
              CatchN=CatchN,popsizeD=popsizeD,cpue=cpue,
              unfishedNaA=unfishedNaA,unfishedNaL=unfishedNaL,
              A0=A0,R0=R0,B0=B0,ExB0=ExB0,natM=natM,steep=rec["steep"],
              sigmaR=rec["sigmaR"])
  return(stk)
}  # end of makesimstock

#' @title plotdynfish generates a plot of a fisheries dynamics
#' 
#' @description plotdynfish generates a plot of a fisheries dynamics. This 
#'     includes plots of predicted CPUE vs observed CPUE, the residuals for the
#'     same (separate plots for different gears), the spawning biomass 
#'     depletion, the recruitment levels, the catches, and the instantaneous
#'     F levels
#'
#' @param outfish a data.frame of the fishery dynamics, as a minimum it must 
#'     obviously include observed and predicted catches and CPUE, spawning
#'     biomass depletion, recruitment, and instantaneous F values. The column
#'     headings of which can be identified in teh columns argument, in which
#'     the order is important.
#' @param console should the plot go to the console, default = TRUE
#' @param prepplot should the plotprep function be used? default = TRUE
#' @param addtitle The default filename if 'fishery_dynamics.png' use addtitle
#'     to add text in front of that. eg addtitle='Speciesname_'
#' @param rundir the directory into which to save the file if console=FALSE,
#'     default = ''
#' @param width the width of the plot within plotprep if used, default=8
#' @param height the height of the plot within plotprep if used, default=7
#' @param nfleet default = 2, determines how many gears are expected can only
#'     be 1 or 2. If nfleet = 1 consider reducing height
#' @param columns column headings identifying the components to be plotted. The
#'     default is: columns=c('year','twlC','aulnC','twlCE','aulnCE', 'twlPCE',
#'     'aulnPCE','deplsB','recruit','twlPF','aulnPF','Trawl','Autoline'), which
#'     identify the year column, the catch columns, the observed CPUE columns, 
#'     the predicted CE cols, the depletion and recruitment cols, then the
#'     instantaneous F estimates by gear, and finally the gear names. If nfleet
#'     = 1 then, obviously one only lists the single gear.
#'
#' @returns nothing but it does generate a plot
#' @export
#'
#' @examples
#' print("Wait on example data - again")
plotdynfish <- function(outfish,console=TRUE,addtitle="",prepplot=TRUE,
                        rundir="",width=8,height=7,nfleet=2,
                        columns=c("year","twlC","aulnC","twlCE","aulnCE",
                                  "twlPCE","aulnPCE","deplsB","recruit",
                                  "twlPF","aulnPF","Trawl","Autoline")) {
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  numcol <- length(columns)
  year <- columns[1]
  catch <- columns[2:(2+nfleet-1)]
  obsCE <- columns[(nfleet+2):(2*nfleet+1)]
  predCE <- columns[(2*nfleet+2):(3*nfleet+1)]
  depl <- columns[(3*nfleet+2)]
  recruit <- columns[(3*nfleet+3)]
  instF <- columns[(3*nfleet+4):(4*nfleet+3)]
  gears <- tail(columns,nfleet)
  fishery <- replacezeros(outfish)
  yrs <- fishery[,"year"]
  if (console & prepplot) {
    plotprep(width=width,height=height,cex=1.0,filename="",verbose=FALSE)
  } else {
    filen <- pathtopath(rundir,paste0(addtitle,"fishery_dynamics.png"))
    plotprep(width=width,height=height,cex=1.0,newdev=TRUE,filename=filen,
             verbose=FALSE)
  }
  if (nfleet == 1) {
    parset(plots=c(3,2),margin=c(0.3,0.4,0.05,0.05),byrow=FALSE)
  } else {
    parset(plots=c(4,2),margin=c(0.3,0.4,0.05,0.05),byrow=FALSE)
  }
  
  maxy <- getmax(fishery[,c(obsCE[1],predCE[1])])
  plot(yrs,fishery[,predCE[1]],type="l",lwd=2,col=1,
       ylab=paste0(gears[1]," CPUE"),ylim=c(0,maxy),yaxs="i",xlab="",
       panel.first=grid())
  points(yrs,fishery[,obsCE[1]],pch=16,cex=1.25,col=2)
  lines(yrs,fishery[,obsCE[1]],lwd=1,col=2,lty=3)
  if (nfleet > 1) {
    legend("topright",c("Predicted","Observed"),col=c(1:nfleet),lwd=3,bty="n",
           cex=1.2)
  }
  if (nfleet == 2) {
    maxy <- getmax(fishery[,c(obsCE[2],predCE[2])])
    plot(yrs,fishery[,predCE[2]],type="l",lwd=2,col=1,
         ylab=paste0(gears[2]," CPUE"),
         ylim=c(0,maxy),yaxs="i",xlab="",panel.first=grid())
    points(yrs,fishery[,obsCE[2]],pch=16,cex=1.25,col=2)
    lines(yrs,fishery[,obsCE[2]],lwd=1,col=2,lty=3)
  }
  maxy <- getmax(fishery[,depl])
  plot(yrs,fishery[,depl],type="l",lwd=2,col=1,ylab="Spawning Depletion",
       ylim=c(0,maxy),yaxs="i",xlab="",panel.first=grid())
  abline(h=c(0.4,0.2),lwd=c(1,1),col=c(3,2))
  maxy <- getmax(fishery[,recruit])
  plot(yrs,fishery[,recruit],type="l",lwd=2,col=1,ylab="Recruitment",
       ylim=c(0,maxy),yaxs="i",xlab="",panel.first=grid())
  twlresid <- fishery[,obsCE[1]]/fishery[,predCE[1]]
  maxy <- getmax(twlresid); miny <- getmin(twlresid)
  plot(yrs,twlresid,type="p",pch=16,col=1,cex=1,
       ylab=paste0(gears[1]," CPUE Residuals"),
       ylim=c(miny,maxy),yaxs="i",xlab="",panel.first=grid())
  lines(yrs,twlresid,lwd=1,col=2,lty=2)
  abline(h=1,lwd=1,col=1)
  if (nfleet == 2) {
    aulnresid <- fishery[,obsCE[2]]/fishery[,predCE[2]]
    maxy <- getmax(aulnresid); miny <- getmin(aulnresid)
    plot(yrs,aulnresid,type="p",pch=16,col=1,cex=1,
         ylab=paste0(gears[2]," CPUE Residuals"),
         ylim=c(miny,maxy),yaxs="i",xlab="",panel.first=grid())
    lines(yrs,aulnresid,lwd=1,col=2,lty=2)
    abline(h=1,lwd=1,col=1)
  }
  if (nfleet == 2) {
    totC <- rowSums(fishery[,c(catch)],na.rm=TRUE)
    totC[which(totC == 0)] <- NA
    maxy <- getmax(totC)
  } else {
    maxy <- getmax(fishery[,catch[1]])
  }
  plot(yrs,fishery[,catch[1]],type="l",lwd=2,col=4,ylab="Catches (t)",
       ylim=c(0,maxy),yaxs="i",xlab="",panel.first=grid())
  if (nfleet == 2) {
    lines(yrs,fishery[,catch[2]],lwd=2,col=2)
    lines(yrs,totC,lwd=2,col=4)
    legend("topleft",c(gears,"Total"),col=c(1,2,4),lwd=3,bty="n",
           cex=1.1)
  }
  maxy <- getmax(fishery[,c(instF)])
  plot(yrs,fishery[,instF[1]],type="l",lwd=2,col=1,ylab="Instantaneous F",
       ylim=c(0,maxy),yaxs="i",xlab="",panel.first=grid())
  if (nfleet == 2) {
    lines(yrs,fishery[,instF[2]],lwd=2,col=2)
    legend("topleft",c(gears),col=c(1:nfleet),lwd=3,bty="n",cex=1.1)
  }
} # end of plotdynfish


#' @title samplefishery generates catcges and CPUE with-without errors
#' 
#' @description samplefishery extracts the catches and the cpue from the
#'     simulated fishery and can add error to both, as desired.
#'
#' @param out the output of simdynF, which generates teh dynamics of the 
#'     simullated fishery
#' @param glb the globals object in teh simulation, from initiateglobals
#' @param errors the variarion imposed on the catches and the cpue, named
#'     catchSD and cpueSD, both with default values = 1e-08, which do not
#'     change the simulated values.
#'
#' @returns histcatch, a matrix of year, catch-by-fleet, and cpue-by-fleet
#' @export
#'
#' @examples
#' \dontrun{
#'   library(codeutils)
#'   library(hplot)
#'   library(agestruct) # obviously setup your own direcdtory structure
#'   data(const2)
#'   glb <- initiateglobals(consts=const2,selabove=1)
#'   stock <- makestock(glb=glb)
#'   out <- simdynF(glb=glb,stk=stock,reps=8,full=TRUE)
#'   histC <- samplefishery(out,glb)
#'   print(histC)
#' }
samplefishery <- function(out,glb,
                          errors=c(catchSD=1e-08,cpueSD=1e-08)) {
  fleets <- glb$fleets
  nflet <- glb$nfleet
  fishery <- out$fishery
  catchnames <- paste0(fleets,"C")
  catches <- fishery[,catchnames]
  nobs <- nrow(catches) * ncol(catches)
  catches <- round(catches * rnorm(nobs,1,sd=errors["catchSD"]))
  cpuenames <- paste0(fleets,"PCE")
  cpue <- fishery[,cpuenames]
  nobs <- nrow(cpue) * ncol(cpue)
  cpue <- cpue * rnorm(nobs,1,sd=errors["cpueSD"])
  histcatch <- as.matrix(cbind(year=fishery[,"year"],catches,cpue))
  return(invisible(histcatch))
} # end of samplefishery

#' @title simdynF describe the ASPM dynamics
#'
#' @description simdynF summarizes the dynamics of an Age-Structured
#'     Integrated Model (ASM). Fishing mortality is implemented as instantaneous 
#'     rates rather than as annual harvest rates. Nt is the numbers-at-age at 
#'     the start (or end) of a year, catchN is the predicted numbers-at-age in 
#'     the catch of each fleet, NumC is the total numbers-at-age in the catch, 
#'     Lt is the predicted numbers-at-size at the start of each year, Gtran*Nt, 
#'     LCflt is the numbers-at-size in the catch for each fleet, Gtran*catchN,
#'     at the end of each year, and LC is the total numbers-at-size in the 
#'     catch, which is Gtran * NumC. These numbers-at-size ignore the fact that 
#'     larger fish are differentially taken by size due to selectivity (a weak 
#'     assumption in purely age-structured models.
#'
#' @param glb the glb data.frame from readdata or built in dataset const2.
#' @param stk the stock object from makestock
#' @param reps how many iterations to use in findF, default = 6
#' @param full should all outputs from dynamics be given. When fitting the 
#'     model, set this to FALSE. Once fitted, change this to TRUE to get all
#'     the required outputs.
#' 
#' @return if !full then a data.frame containing the fishery dynamics according 
#'     to the input arguments glb and stk. Includes Year, Catch, PredC, SpawnB, 
#'     ExploitB, FullH, CPUE, PredCE, Deplete, Recruit, FullF. if full, then a 
#'     list of the fishery dynamics, plus Nt, NumC, catchN, Lt, LC, and LCflt.
#' @export
#'
#' @examples
#' data(simconst2)
#' glb <- initiateglobals(simconst2,agesorsizes="ages",selabove=1)
#' stock <- makesimstock(glb=glb)
#' out <- simdynF(glb=glb,stk=stock,reps=6,full=TRUE)
#' str(out)
#' print(round(out$fishery,3)) 
simdynF <- function(glb,stk,reps=6,full=FALSE) { 
  # glb=glb; stk=stock; reps=6; full=TRUE; 
  aaw <- glb$WaA
  aam <- glb$MatA
  sel <- as.matrix(glb$Sel)
  R0 <- stk$R0
  B0 <- stk$B0
  sigR <- stk$sigmaR
  M <- glb$natM
  nfleet <- glb$nfleet
  fleets <- glb$fleets
  years <- glb$years
  steep <- glb$fishbiol$rec["steep"]
  initdepl <- glb$initdepl
  allyrs <- c((years[1]-1),years)
  if (initdepl < 1.0) {
    dep <- doDepletion(indepl=initdepl,stk,glb,inc=0.02)
    spb <- SpB(dep$Ndepl,aam,aaw)
    Rinit <- bhsim(spb,R0,B0,steep,sigmaR=1e-07)
  } else {
    Rinit <- R0
  }
  nyrs <- glb$nyrs
  nyr1 <- nyrs + 1
  ages <- glb$ages
  nages <- glb$nages
  maxage <- glb$maxage
  catchcol <- paste0(fleets,"C")
  predCcol <- paste0(fleets,"PC")
  exBcols <- paste0(fleets,"eB")
  CEcols <- paste0(fleets,"CE")
  predCEcols <- paste0(fleets,"PCE")
  yrFcols <- paste0(fleets,"PF") 
  Nt <- matrix(0,nrow=nages,ncol=nyr1,dimnames=list(ages,allyrs))
  NumC <- Nt
  sizes <- glb$sizes
  nsizes <- glb$nsizes
  Lt <- LC <- matrix(0,nrow=nsizes,ncol=nyr1,dimnames=list(sizes,allyrs))
  LCflt <- array(0,dim=c(nsizes,nyrs,nfleet),dimnames=list(sizes,years,fleets))
  Gtran <- glb$growtran[,,1,1]
  columns <- c("year",catchcol,predCcol,CEcols,predCEcols,exBcols,
               "spawnB","deplsB","recruit",yrFcols)
  fishery <- matrix(NA,nrow=nyr1,ncol=length(columns),
                    dimnames=list(0:nyrs,columns))
  fishery[,"year"] <- as.numeric(rownames(stk$catchB))
  fishery[,catchcol] <- stk$catchB
  fishery[1,exBcols] <- stk$ExpB[1,]
  fishery[1,"spawnB"] <- stk$Bsp[1]
  fishery[1,"deplsB"] <- glb$initdepl
  fishery[1,"recruit"] <- stk$R0
  hS <- exp(-M/2)   # for midyear CPUE
  surv <- exp(-M)
  Nt[,1] <- stk$NaA[,1]
  Lt[,1] <- Gtran %*% Nt[,1]
  obscatch <- as.matrix(fishery[,catchcol])  # retain original catch data
  obscatch[which(is.na(obscatch))] <- 0
  predCN <- matrix(NA,nrow=nages,ncol=nfleet,dimnames=list(ages,fleets))
  catchN <- array(NA,dim=c(nages,nyrs,nfleet),
                  dimnames=list(ages,years,fleets))
  set.seed(glb$randseed)
  #for (yr in 1:18) { #nyr1) {  # yr=40
  for (yr in 1:nyrs) {  # yr=1
    obsC <- obscatch[(yr+1),]
    spb <- fishery[yr,"spawnB"]  #SpB(Nt[,(yr-1)],aam,aaw)
    Nt[1,(yr+1)] <- bhsim(spb,R0,B0,steep,sigmaR=sigR)
    fishery[(yr+1),"recruit"] <- Nt[1,(yr+1)]
    nextEN <- nextNT <- numeric(nages)  # exploitable and spawning NaA
    if (nfleet == 1) {  # Single fleet
      yrF <- findFs(obsC,Nyr=Nt[,yr],sel=sel,aaw=aaw,M=M,reps=reps)
      sF <- sel * yrF
      psF <- sF[2:(nages-1)]
      msF <- sF[nages]
      catchN[,yr,1] <- (sF/(M + sF)) * Nt[,yr] * (1 - exp(-(M + sF)))
      fishery[(yr+1),predCcol] <- sum(catchN[,yr,1] * aaw/1000.0)
      fishery[(yr+1),yrFcols] <- yrF
      Nt[1,(yr+1)] <- Nt[1,(yr+1)] - catchN[1,yr,1]
      Nt[2:(nages-1),(yr+1)] <- (Nt[1:(nages-2),yr] * exp(-(M + psF)))
      Nt[nages,(yr+1)] <- (Nt[nages,yr] + Nt[(nages-1),yr]) *
                                             exp(-(M + msF))
      nextEN[2:(nages-1)] <- (Nt[1:(nages-2),yr] * exp(-(M + psF*(yrF/2.0))))
      nextEN[nages] <- (Nt[nages,yr] + Nt[(nages-1),yr]) * 
                                          exp(-(M + msF * (yrF/2.0)))
    } else {  # Multi-fleet one positive catch solution
      if (countgtzero(obsC) < nfleet) {
        pickft <- which(obsC > 0)
        yrF <- findFs(obsC[pickft],Nyr=Nt[,yr],sel=sel[,pickft],
                      aaw=aaw,M=M,reps=reps)
        sF <- sel[,pickft] * yrF
        psF <- sF[2:(nages-1)] # sel of ages 1 - (maxage-1)
        msF <- sF[nages]       # sel of maxage = plusgroup  
        predCN[,pickft] <- (sF/(M + sF)) * Nt[,yr] * (1 - exp(-(M + sF)))
        fishery[(yr+1),predCcol[pickft]] <- sum(predCN[,pickft] * aaw/1000.0)
        fishery[(yr+1),yrFcols[pickft]] <- yrF
        catchN[,yr,pickft] <- predCN[,pickft]
        Nt[1,(yr+1)] <- Nt[1,(yr+1)] - catchN[1,yr,pickft]
        # main dynamics
        Nt[2:(nages-1),(yr+1)] <- (Nt[1:(nages-2),yr] * exp(-(M + psF)))
        Nt[nages,(yr+1)] <- (Nt[nages,yr] + Nt[(nages-1),yr] * 
                                               exp(-(M + msF)))
        nextEN[2:(nages-1)] <- (Nt[1:(nages-2),yr] * 
                                  exp(-(M + psF*(yrF/2.0))))
        nextEN[nages] <- (Nt[nages,yr] + Nt[(nages-1),yr]) * 
                                            exp(-(M + msF * (yrF/2.0)))
      } else {  # Multi-fleet all with catches 
        yrF <- numeric(nfleet)
        for (i in 1:nfleet) { # i = 1
          yrF[i] <- findFs(cyr=obsC[i],Nyr=Nt[,yr],sel=sel[,i],aaw=aaw,M=M,
                           reps=reps)
          predCN[,i] <- ((sel[,i] * yrF[i])/(M + (sel[,i] * yrF[i]))) * 
                          Nt[,yr] * (1 - exp(-(M + sel[,i] * yrF[i]))) 
          fishery[(yr+1),predCcol[i]] <- sum(predCN[,i] * aaw/1000.0)
          fishery[(yr+1),yrFcols[i]] <- yrF[i]
        }
        catchN[,yr,] <- predCN
        Nt[1,(yr+1)] <- Nt[1,(yr+1)] - sum(catchN[1,yr,],na.rm=TRUE)
        nextEN <- nextNt <- numeric(nages)
        multe1 <- multe2 <- mult2 <- mult1 <- exp(-M)
        for (ft in 1:nfleet) {
          mult1 <- mult1 * exp(-sel[2:(nages-1),ft] * yrF[ft])
          mult2 <- mult2 * exp(-sel[nages,ft] * yrF[ft])
          multe1 <- multe1 * exp(-sel[2:(nages-1),ft] * (yrF[ft]/2.0))
          multe2 <- multe2 * exp(-sel[nages,ft] * (yrF[ft]/2.0))
        }
        nextNt[2:(nages-1)] <- (Nt[1:(nages-2),yr] * mult1)
        nextNt[nages] <- (Nt[nages,yr] + Nt[(nages-1),yr]) * mult2 
        Nt[2:nages,(yr+1)] <- nextNt[2:nages]
        nextEN[2:(nages-1)] <- (Nt[1:(nages-2),yr] * multe1)
        nextEN[nages] <- (Nt[nages,yr] + Nt[(nages-1),yr]) * multe2
      } # end of all-fleets loop 
    } # end of multifleet loop
    fishery[(yr+1),"spawnB"] <- SpB(Nt[,(yr+1)],aam,aaw) 
    if (nfleet == 1) { 
      fishery[(yr+1),exBcols] <- ExB(nextEN,sel[,1],glb$WaA)
      NumC[,(yr+1)] <- catchN[,yr,] 
      } else {
        for (ft in 1:nfleet) 
          fishery[(yr+1),exBcols[ft]] <- ExB(nextEN,sel[,ft],glb$WaA)
      NumC[,(yr+1)] <- rowSums(catchN[,yr,],na.rm=TRUE)
    }    
  } # end of yr loop
  fishery[,"deplsB"] <- fishery[,"spawnB"]/B0
  ExpB <- as.matrix(fishery[,exBcols])
  catchnames <- catchcol
  for (ft in 1:nfleet) { # ft =1
    pickC <- which(fishery[,catchnames[ft]] > 0)
    tmp <- ExpB[pickC,ft] * glb$qc[ft]
    fishery[pickC,predCEcols[ft]] <- tmp/(mean(tmp)) # mean = 1.0 
  }
  if (full) {
    for (yr in 2:nyr1) {
      Lt[,yr] <- Gtran %*% Nt[,yr]
      LC[,yr] <- Gtran %*% NumC[,yr]
      for (ft in 1:nfleet) LCflt[,(yr-1),ft] <- Gtran %*% catchN[,(yr-1),ft]
    }    
    out <- list(fishery=as.data.frame(fishery),Nt=Nt,NumC=NumC,catchN=catchN,
                Lt=Lt,LC=LC,LCflt=LCflt)
    return(out)
  } else {
    return(as.data.frame(fishery))
  }
} # end of simdynF


#' @title template2F1S data template for a 2 fleet 1 region/stock fishery
#'
#' @param rundir directory in which to find the data file and run the analysis
#' @param filename the name of the data file to be produced, default=
#'     'F21S.csv'.
#'
#' @return the function write a data file to rundir and returns the filename
#' @export
#'
#' @examples
#' \dontrun{
#'   rundir <- tempdir()
#'   dirExists(rundir,verbose=TRUE)
#'   template2F1S(rundir,filename="Flt2Stock1.csv")
#'   dir(rundir)
#' }
template2F1S <- function(rundir,filename="F2S1.csv") {
  filename <- pathtopath(rundir,filename)
  cat("Data for a 2 Fleet 1 Region 1 Stock model  \n\n",
      file=filename,append=FALSE)
  cat("#STRUCTURE,,, \n",file=filename,append=FALSE)
  cat("randseed, 6924062, for repeatability \n",file=filename,append=TRUE)
  # cat("nregion, 1,,, number of regions, imples 1 stock  \n",
  #     file=filename,append=TRUE)
  # cat("regname, east,,, labels for region  \n",file=filename,append=TRUE)
  # cat("nsex, 1,,, number of sexes \n",file=filename,append=TRUE)
 #  cat("LFstruct,0,72,1, sequence for lengths \n",file=filename,append=TRUE)
  cat("agestruct,0,30,1, sequence for ages \n",file=filename,append=TRUE)
  cat("nfleet,2,,, fleetnumbers \n",file=filename,append=TRUE)
  fleets <- c("twl", "auln")
  cat("fleets, twl, auln  \n",file=filename,append=TRUE)
  cat("selecttype, logistic, domed,  # currently only logistic or domed \n",
      file=filename,append=TRUE)
  cat("initdepl, 1.0,, the initial depletion level \n",
      file=filename,append=TRUE)
  cat("\n\n",file=filename, append=TRUE)
  selnames <- paste0(fleets,"S")
  cat("#FISHERY,    qc,     sel50,   deltas,  \n",file=filename, append=TRUE)
  cat(selnames[1],",  1.4E-04,    4.5,    0.75, \n",file=filename, append=TRUE)
  cat(selnames[2],",  9.0E-05,    10,  15,    10,  30,   -7,   0.5,  \n",
      file=filename, append=TRUE)
  cat("#         qc         peak1 peak2  asc  dsc  selmin selmax \n",
      file=filename, append=TRUE)
  cat("\n\n",file=filename, append=TRUE)
  cat("#BIOLOGY,,, \n",
      file=filename, append=TRUE)
  cat("M,         0.21 ,,, \n",file=filename, append=TRUE)
  cat("Linf,	    56,,, \n",file=filename, append=TRUE)
  cat("K,	       0.20,,, \n",file=filename, append=TRUE)
  cat("t0,	      -0.1,,, \n",file=filename, append=TRUE)
  cat("growCV,	   0.075,,, \n",file=filename, append=TRUE)
  cat("Wta,	      5.88E-06,,, \n",file=filename, append=TRUE)
  cat("Wtb,	      3.31,,, \n",file=filename, append=TRUE)
  cat("steep,      0.7,,, \n",file=filename, append=TRUE)
  cat("Age50M,	   3,,, \n",file=filename, append=TRUE)
  cat("deltaM,	   0.75,,, \n",file=filename, append=TRUE)
  cat("\n\n",file=filename, append=TRUE)
  cat("#HISTORICALCATCH,46,,, \n",file=filename, append=TRUE)
  cat("1975,NA,NA,NA,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1976,9,NA,2.2648,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1977,20,NA,1.7314,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1978,29,NA,1.7318,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1979,58,NA,1.4174,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1980,75,NA,1.3511,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1981,82,NA,1.8418,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1982,108,NA,1.9532,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1983,133,NA,1.8696,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1984,163,NA,1.9487,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1985,167,NA,1.5997,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1986,208,NA,1.8995,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1987,271,NA,1.429,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1988,271,NA,1.2055,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1989,285,NA,1.4616,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1990,289,NA,1.649,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1991,297,NA,1.8187,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1992,291,NA,1.5823,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1993,280,NA,1.8183,NA,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1994,298,13,1.6311,3.9495,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1995,308,28,1.3231,2.2962,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1996,316,46,1.2011,2.1636,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1997,330,60,0.7488,2.251,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1998,363,89,1.3769,2.1577,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("1999,430,111,1.3956,2.1396,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2000,445,145,1.1556,1.4165,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2001,510,156,0.4726,1.7317,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2002,525,156,0.7709,1.0261,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2003,500,155,0.5209,1.2184,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2004,460,175,0.7361,0.8835,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2005,410,168,0.465,0.5055,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2006,394,145,0.4892,0.7799,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2007,355,128,0.4421,0.4345,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2008,320,124,0.3189,0.5281,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2009,304,123,0.3249,0.3872,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2010,275,120,0.3362,0.2824,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2011,271,118,0.1487,0.3024,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2012,262,110,0.1654,0.1923,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2013,246,99,0.2291,0.1788,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2014,146,54,0.2392,0.2359,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2015,146,54,0.2898,0.2093,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2016,146,54,0.3255,0.2273,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2017,146,54,0.2036,0.3024,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2018,146,54,0.2552,0.3559,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2019,146,54,0.3104,0.316,yr_tc_ac, \n",file=filename,append=TRUE)
  cat("2020,146,54,0.2797,0.301,yr_tc_ac, \n",file=filename,append=TRUE)
  return(invisible(filename))
} # end of template2F1S

#' @title unfished - generates the numbers at age for an unfished population
#'
#' @description unfished - generates the numbers at age for an unfished
#'    population, and determines the recruitment dynamics
#' @param glb - the global constants object containing biology, fishery and
#'     structure
#' @return a list containing N0, R0, A0, SpB0, and expB0
#' @export
#' @examples
#' \dontrun{
#' library(fmr)
#' datafileTemplate(filename="C:/A_CSIRO/Rcode/agestructrun/constants.csv")
#' const <- getconstants("C:/A_CSIRO/Rcode/agestructrun/constants.csv")
#' globals <- initiateglobals(const)
#' unfish <- unfished(globals)
#' }
unfished <- function(glb) {   # glb <- globals
  maxage <- glb$maxage
  nages <- glb$nages
  M <- glb$biology["M",]
  surv = exp(-M) 
  Nt <- numeric(nages)
  Nt[1] <- 1        # sets R0 = 1 for initiation NaA[0,0]
  for (age in 2:(nages-1)) Nt[age] <- Nt[age-1] * surv 
  Nt[nages] <- (Nt[nages-1] * surv)/(1-surv)
  # Estimate the biomass A0 that would generate a recruitment of 1.0
  A0 <- sum(Nt * glb$MatA * glb$WaA)/1000.0
  B0 <- glb$rec["B0"]
  R0 <- B0 / A0
  Nt <- R0*Nt # Generate the initial age distribution given R0
  nfleet <- glb$nfleet
  ExB0 <- numeric(nfleet)
  for (flt in 1:nfleet) ExB0[flt] <- sum(Nt * glb$Sel[,flt] * glb$WaA)/1000.0
 # SpB0 <- sum(Nt * glb$MatA * glb$WaA)/1000.0
  res <- list(Nt, R0, A0, B0,ExB0)
  names(res) <- c("Nt","R0","A0","B0","ExB0")
  return(res)
}  # end of unfished


#' @title writedatafile writes out a simulated data file to rundir/filename
#' 
#' @description writedatafile writes out a simulated data file to 
#'     rundir/filename that is suitable for use with dynF2, ie only catches and
#'     cpue. The intent is to add the facility to include age- and size-
#'     composition data later.
#'
#' @param const the object made by getconstants aimed at simulating dynamics
#' @param glb the globals object obtained from initiategeglobals
#' @param histC an object inside the output from getconstants
#' @param rundir the directory path to the subdirectory in which all analysis 
#'     data and results are stored.
#' @param filename  the name of the data file, MUST be a '.csv' file
#' @param datatitle line 1 of data file decribing its contents, default = ''
#'
#' @returns nothing but it does write a data file to rundir/filename
#' @export
#'
#' @examples
#' print("see help for makedatafile, which uses writedatafile")
writedatafile <- function(const,glb,histC,rundir,filename,datatitle) {
  filename <- pathtopath(rundir,filename)
  cat(datatitle,"  \n\n",file=filename,append=FALSE)
  newseed <- glb$randseed
  cat("#STRUCTURE,,, \n",file=filename,append=TRUE)
  cat("randseed,",newseed,", for repeatability \n",file=filename,append=TRUE)
  struct <- const$structure
  ages <- struct$ages
  strage <- paste0(c(ages[1],tail(ages,1),ages[2]-ages[1]),collapse=",")
  cat("agestruct,",strage,", sequence for ages \n",file=filename,append=TRUE)
  cat("nfleet,",glb$nfleet,",,, fleetnumbers \n",file=filename,append=TRUE)
  flts <- paste0(glb$fleets,collapse=",")
  cat("fleets,",flts,",  \n",file=filename,append=TRUE)
  sels <- paste0(glb$seltype,collapse=",")
  cat("selecttype,",sels,",# currently only logistic or domed \n",
      file=filename,append=TRUE)
  cat("initdepl,",glb$initdepl,",, the initial depletion level \n",
      file=filename,append=TRUE)
  cat("\n\n",file=filename, append=TRUE)
  fishsel <- const$fishery
  label <- paste0("#FISHERY, ",paste0(names(fishsel[[1]]),collapse=","))
  cat(label,",  \n",file=filename, append=TRUE)
  for (flt in 1:glb$nfleet) { # flt=1
    cat(paste0(glb$fleets[flt],"S"),",",paste0(fishsel[[flt]],collapse=","),
        ", \n",file=filename, append=TRUE)
  }
  cat("#         qc         peak1 peak2  asc  dsc  selmin selmax \n",
      file=filename, append=TRUE)
  cat("\n\n",file=filename, append=TRUE)
  biol <- const$biology
  steep <- const$fishbiol$rec["steep"]
  cat("#BIOLOGY,,, \n",
      file=filename, append=TRUE)
  cat("M,        ",biol["M",],",,, \n",file=filename, append=TRUE)
  cat("Linf,	   ",biol["Linf",],",, \n",file=filename, append=TRUE)
  cat("K,	       ",biol["K",],",, \n",file=filename, append=TRUE)
  cat("t0,	     ",biol["t0",],",, \n",file=filename, append=TRUE)
  cat("growCV,	 ",biol["growCV",],",, \n",file=filename, append=TRUE)
  cat("Wta,	     ",biol["Wta",],",, \n",file=filename, append=TRUE)
  cat("Wtb,	     ",biol["Wtb",],",, \n",file=filename, append=TRUE)
  cat("steep,    ",steep,",, \n",file=filename, append=TRUE)
  cat("Age50M,	 ",biol["age50M",],",, \n",file=filename, append=TRUE)
  cat("deltaM,	 ",biol["deltaM",],",, \n",file=filename, append=TRUE)
  cat("\n\n",file=filename, append=TRUE)
  cat("#HISTORICALCATCH,45,,, \n",file=filename, append=TRUE)
  numrow <- nrow(histC)
  for (i in 1:numrow) {
    dat <- paste0(round(histC[i,],4),collapse=",")
    label <- paste0(dat,", fishdat, \n")
    cat(label,file=filename, append=TRUE)
  }
  return(filename)
} # end of writedatafile

