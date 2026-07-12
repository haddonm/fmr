
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
  #  infile=pathtopath(rundir,"test2F1S.csv")
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
  fishery <- makelist(fleets)
  for (i in 1:nfleet) {  # i = 1
    fltlabel <- paste0("fleet",i)     
    if (seltype[i] == "logistic") {
      label <- c("qc","sel50","deltasel")      
      fishery[[i]] <- getvect(fltlabel,datain,n=length(label))
      names(fishery[[i]]) <- c("qc","sel50","deltasel")     
    }
    if (seltype[i] == "domed") {
      fishery[[i]] <- getvect(fltlabel,datain,n=7)
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
  biology["Wta",] <- getvect("WaLa",datain,ncols)
  biology["Wtb",] <- getvect("WaLb",datain,ncols)
  biology["age50M",] <- getvect("Age50M",datain,ncols)
  biology["deltaM",] <- getvect("deltaM",datain,ncols)
  M <- getsingleNum("M",datain)
  steep <- getsingleNum("steepness",datain)
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
getvect <- function(varname,intxt,n) {
  begin <- grep(varname,intxt)
  if (length(begin) > 0) {
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
#' @param selabove are there ages or sizes that are never selected. Identify
#'     the upper index in the outsel vector here. eg. if age 0 is never selected
#'     then set zerocomps = 1, if ages 0-3 then zerocomps = 4
#' 
#' @return a list containing LaA, WaA, Sel, MatA, growtran, and L50, and more
#' @export
#' 
#' @examples
#' data(const)
#' globals <- initiateglobals(const)
#' str(globals,max.level=1)
initiateglobals <- function(consts,selabove=0) { #  consts=const; selabove=1
  biology <- consts$biology
  fishbiol <- consts$fishbiol
  fishery <- consts$fishery
  struct <- consts$structure
  histC <- consts$histC
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
  Sel <- matrix(0,nrow=nsizes,ncol=nfleet,dimnames=list(sizes,fleets))
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
    Sel[,j] <- calcsel(ageorlen=sizes,fishprop=propfish,
                       typeselect=seltype[j],zerocomps=selabove)
  }   
  # growtran <- makeSTM(pars,glb$sizes,funL=vBfabens)
  qc <- unlist(lapply(fishery,"[[","qc"))
  ans <- list(ages,maxage,nages,sizes,nsizes,LaA,WaA,Sel,MatA,growtran,
              grow,biology,fishbiol,fishbiol$rec,fishbiol$initdepl,WaL,gpar,
              qc,nyrs,startyr:endyr,nfleet,fleets,seltype,nregion,regions,histC,
              nsex,sexes,sexreg,nsexreg,fltreg,nfltreg)
  names(ans) <- c("ages","maxage","nages","sizes","nsizes","LaA","WaA","Sel",
                  "MatA","growtran","grow","biology","fishbiol","rec",
                  "initdepl","wal","gpar","qc","nyrs","years","nfleet",
                  "fleets","seltype","nregion","regions","histC",
                  "nsex","sexes","sexreg","nsexreg","fltreg","nfltreg")
  return(ans)
}  # end of initiateglobals

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
  props[,"maa"] <- logist(inL50=biol["age50M",],delta=biol["deltaM",],depend=ages)
  for (flt in 1:nfleet) {
    propfish <- fishery[[flt]]
    props[,fleets[flt]] <- calcsel(ageorlen=ages,fishprop=propfish,
                                   typeselect=glb$selecttype[flt],zerocomps=selabove)
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
#' library(fmr)
#' data(const)
#' globals <- initiateglobals(const)
#' stock <- makesimstock(globals)
#' str(stock)
makesimstock <- function(glb) { # glb=glb
  ages <- glb$ages
  years <- glb$years
  sizes <- glb$sizes
  nages <- glb$nages
  maxage <- max(glb$ages)
  nyrs <- length(glb$years)
  yearlab <- c((years[1]-1),years)
  nsizes <- glb$nsizes
  nfleet <- glb$nfleet
  fleets <- glb$fleets
  nregion <- glb$nregion
  regions <- glb$regions
  sexreg <- glb$sexreg
  fltreg <- glb$fltreg
  NaA <- matrix(0,nrow=nages,ncol=(nyrs+1),  # currently A1-S1
                dimnames=list(seq(0,maxage,1),yearlab))
  Bsp <- matrix(0,nrow=(nyrs+1),ncol=1,dimnames=list(yearlab,"SpawnB"))
  ans <- unfished(glb)
  Bsp[1] <- ans$B0
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
   # Nt, R0, A0,SpB0,expB0,N0mat
  NaA[,1] <- ans$Nt
  unfishedNaA <- ans$Nt
  A0 <- ans$A0
  R0 <- ans$R0
  B0 <- ans$B0
  # ExB0 <- ans$expB0
  # if (nfleet == 1) {
  #   ExpB[1] <- ExB0
  # } else {
  #   ExpB[1,] <- ExB0
  # }
  unfishedNaL <- glb$growtran %*% NaA[,1]
  popsizeD[,1] <- unfishedNaL
  stk <- list(NaA=NaA,Bsp=Bsp,ExpB=ExpB,catchB=catchB,Harvest=Harvest,
              CatchN=CatchN,popsizeD=popsizeD,cpue=cpue,
              unfishedNaA=unfishedNaA, unfishedNaL=unfishedNaL,
              A0=A0,R0=R0,B0=B0) #,ExB0=ExB0)
  return(stk)
}  # end of makesimstock



#' @title template2F1S data template for a 2 fleet 1 region/stock fishery
#'
#' @param rundir directory in which to find the data file and run the analysis
#' @param filename the name of the data file to be produced, default=
#'     'Fleet2Region1.csv'.
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
template2F1S <- function(rundir,filename="Fleet2Region1.csv") {
  filename <- pathtopath(rundir,filename)
  cat("Data for a 2 Fleet 1 Region 1 Stock model  \n\n",
      file=filename,append=FALSE)
  cat("#STRUCTURE,,, \n",file=filename,append=FALSE)
  cat("randseed, 8684569, for repeatability \n",file=filename,append=TRUE)
  # cat("nregion, 1,,, number of regions, imples 1 stock  \n",
  #     file=filename,append=TRUE)
  # cat("regname, east,,, labels for region  \n",file=filename,append=TRUE)
  # cat("nsex, 1,,, number of sexes \n",file=filename,append=TRUE)
 #  cat("LFstruct,0,72,1, sequence for lengths \n",file=filename,append=TRUE)
  cat("agestruct,0,30,1, sequence for ages \n",file=filename,append=TRUE)
  cat("nfleet,2,,, fleetnumbers \n",file=filename,append=TRUE)
  cat("fleets, twl, auln,  \n",file=filename,append=TRUE)
  cat("selecttype, logistic, domed,  # currently only logistic or domed \n",
      file=filename,append=TRUE)
  cat("initdepl, 1.0,, the initial depletion level \n",
      file=filename,append=TRUE)
  cat("\n\n",file=filename, append=TRUE)
  cat("#FISHERY,    qc,     sel50,   deltas,  \n",file=filename, append=TRUE)
  cat("fleet1,   1.4E-04,    4.5,    0.75, \n",file=filename, append=TRUE)
  cat("fleet2,   9.0E-05,    10,  15,    10,  30,   -7,   0.5,  \n",
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
  cat("WaLa,	     5.88E-06,,, \n",file=filename, append=TRUE)
  cat("WaLb,	     3.31,,, \n",file=filename, append=TRUE)
  cat("steepness, 0.7,,, \n",file=filename, append=TRUE)
  cat("Age50M,	   3,,, \n",file=filename, append=TRUE)
  cat("deltaM,	   0.75,,, \n",file=filename, append=TRUE)
  cat("\n\n",file=filename, append=TRUE)
  cat("#HISTORICALCATCH,45,,, \n",file=filename, append=TRUE)
  cat("1975,0,0,0,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1976,8,0,1.817950618,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1977,18,0,1.879496675,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1978,26,0,1.59365282,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1979,52,0,1.800166222,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1980,68,0,1.897771455,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1981,74,0,1.751619141,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1982,97,0,1.818011726,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1983,120,0,1.72671127,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1984,147,0,1.585619174,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1985,150,0,1.822948399,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1986,187,0,1.502407019,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1987,244,0,1.53118141,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1988,244,0,1.681680949,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1989,256,0,1.372333195,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1990,260,0,1.487280263,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1991,267,0,1.427421894,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1992,262,0,1.260075269,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1993,252,0,1.38153378,0,yr-tc-ac, \n",file=filename,append=TRUE)
  cat("1994,268,12,1.240763546,2.48431198,abov, \n",file=filename,append=TRUE)
  cat("1995,277,25,1.231613821,2.197692164,abov, \n",file=filename,append=TRUE)
  cat("1996,284,41,1.282217138,2.639122034,abov, \n",file=filename,append=TRUE)
  cat("1997,297,54,1.286378001,2.230568499,abov, \n",file=filename,append=TRUE)
  cat("1998,327,80,1.085279131,2.154980718,abov, \n",file=filename,append=TRUE)
  cat("1999,387,100,1.047047386,1.946129626,abov, \n",file=filename,append=TRUE)
  cat("2000,400,130,0.832557651,2.009846266,abov, \n",file=filename,append=TRUE)
  cat("2001,459,140,0.760019702,1.414459913,abov, \n",file=filename,append=TRUE)
  cat("2002,472,140,0.483588224,1.210468368,abov, \n",file=filename,append=TRUE)
  cat("2003,450,140,0.641863502,1.135852524,abov, \n",file=filename,append=TRUE)
  cat("2004,414,158,0.55298688,0.952771825,abov, \n",file=filename,append=TRUE)
  cat("2005,369,151,0.576119575,0.670068941,abov, \n",file=filename,append=TRUE)
  cat("2006,355,130,0.468798891,0.342674742,abov, \n",file=filename,append=TRUE)
  cat("2007,320,115,0.342799012,0.641914919,abov, \n",file=filename,append=TRUE)
  cat("2008,288,112,0.364003382,0.397247141,abov, \n",file=filename,append=TRUE)
  cat("2009,274,111,0.425092835,0.181578642,abov, \n",file=filename,append=TRUE)
  cat("2010,248,108,0.314242866,0.663728398,abov, \n",file=filename,append=TRUE)
  cat("2011,244,106,0.306935687,0.250832979,abov, \n",file=filename,append=TRUE)
  cat("2012,236,99,0.315813651,0.460257868,abov, \n",file=filename,append=TRUE)
  cat("2013,221,89,0.243150639,0.363334106,abov, \n",file=filename,append=TRUE)
  cat("2014,110,44,0.147133082,0.228131383,abov, \n",file=filename,append=TRUE)
  cat("2015,110,44,0.215242012,0.535650776,abov, \n",file=filename,append=TRUE)
  cat("2016,110,44,0.263228871,0.38243359,abov, \n",file=filename,append=TRUE)
  cat("2017,110,44,0.122266773,0.323630679,abov, \n",file=filename,append=TRUE)
  cat("2018,110,44,0.453325032,0.275539674,abov, \n",file=filename,append=TRUE)
  cat("2019,110,44,0.260213502,0.747151256,abov, \n",file=filename,append=TRUE)
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
#' library(agestruct)
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
  R0 <- glb$rec["R0"]
  B0 <- R0 * A0
  # Generate the initial age distribution given R0
  Nt <- R0*Nt
  SpB0 <- sum(Nt * glb$MatA * glb$WaA)/1000.0
  res <- list(Nt, R0, A0, B0)
  names(res) <- c("Nt","R0","A0","B0")
  return(res)
}  # end of unfished