

#' @title aspmindexfit plots an aspm fit to CPUE or Index data
#' 
#' @description aspmindexfit plots an aspm fit to CPUE or Index data with the
#'     option of including 95% Log-Normal CI around the model fit.
#'
#' @param infish the matrix of fishery dynamics from the dynamicsH, dynamicsF,
#'     of dynF functions
#' @param glb the globals object for the aspm analysis used
#' @param CI the output of the getLNCI function that defines Log-Normal CI
#' @param rundir default = "", give a full path if saving a png file
#' @param CIlwd width of CI lines default = 1
#' @param CIcol colour of CI lines, default = 4 = blue
#' @param console should the plot be sent to the console, default = TRUE. If
#'     set = FALSE, then rundir can be set to identify where the png of the 
#'     plot will be saved
#'     
#' @seealso{
#'    \link{getLNCI}, \link{dynamicsH}, \link{dynamicsF}, \link{dynF}
#' }     
#'
#' @return nothing but it does generate a plot
#' @export
#'
#' @examples
#' data("westroughy")
#' fish <- westroughy$fish; glb <- westroughy$glb
#' props <- westroughy$props
#' pars <- c(7.1,-1,-7.7)
#' ans <- fitASPM(pars,dynamicsH,infish=fish,inglb=glb,inprops =props)
#' out <- dynamicsH(ans$estimate,infish=fish,inglb=glb,inprops=props,
#'                  full=TRUE)
#' ceCI <- getLNCI(out$fishery[,"predCE"],exp(ans$estimate[2]))
#' aspmindexfit(infish=out$fishery,glb=glb,CI=ceCI,console=TRUE)
aspmindexfit <- function(infish,glb,CI=NULL,rundir="",CIlwd=1,CIcol=4,
                         console=TRUE) {
  colnames(infish) <- tolower(colnames(infish))
  yrs <- infish[,"year"]
  cpue <- infish[,"cpue"]
  predce <- infish[,"predce"]
  ymax <- getmax(c(cpue,predce))
  if (inherits(CI,"matrix")) ymax <- getmax(CI[,"upper"]) 
  filen=""
  if (!console) filen <- pathtopath(rundir,paste0(glb$spsname,"_aspm.png"))
  if (console) plotprep(width=9,height=5,cex=1.0,filename=filen)
  parset()
  plot(yrs,cpue,type="p",pch=16,col=2,cex=1.0,ylim=c(0,ymax),
       yaxs="i",xlab="",panel.first=grid(),ylab="Relative Abundance Index")
  lines(yrs,predce,lwd=2,col=1)
  if (inherits(CI,"matrix"))  {
    segments(x0=yrs,y0=CI[,1],x1=yrs,y1=CI[,3],lwd=CIlwd,col=CIcol)
  }
  if (!console) dev.off()
} # end of aspmindexfit

#' @title aspmphaseplot - plots the phase plot of harvest rate vs biomass
#' 
#' @description aspmphaseplot uses the output from displayModel to plot up 
#'     the phase plot of harvest rate vs Biomass, marked with the limit and
#'     default targets. It identifies the start and end years (green and red
#'     dots) and permits the stock status to be determined visually. It also 
#'     plots out the catch time-series and harvest rate time-series to aid in
#'     interpretation of the phase plot.
#'
#' @param fishery the object output by the function dynamics, containing the 
#'     fishery dynamics (Year, Catch, PredC, SpawnB, ExploitB, FullH, CPUE,
#'     PredCE, and Deplete).
#' @param prod the matrix containing the production data from the function
#'     getProductionC
#' @param ans the vector of results from the function prodASPM
#' @param Blim the limit reference point, defaults to 0.2 so that 0.2B0 is used.
#' @param filename default is empty. If a filename is put here a .png file
#'     with that name will be put into the working directory. 
#' @param resol the resolution of the png file, defaults to 200 dpi
#' @param fnt the font used in the plot and axes. Default=7, bold Times. Using
#'     6 gives Times, 1 will give SansSerif, 2 = bold Sans
#'
#' @return an invisible list of B0, Bmsy, Hmsy, and Hlim.
#' @export
#'
#' @examples
#' print("wait on a better example")
aspmphaseplot <- function(fishery,prod,ans,Blim=0.2,filename="",resol=200,
                          fnt=7) {
   lenfile <- nchar(filename)
   if (lenfile > 3) {
      end <- substr(filename,(lenfile-3),lenfile)
      if (end != ".png") filename <- paste0(filename,".png")
      png(filename=filename,width=5.5,height=5.0,units="in",res=resol)
   } 
   B0 <- ans["B0"] 
   Bmsy <- ans["Bmsy"]
   Hmsy <- ans["Hmsy"]
   Btarg <- ans["Btarg"]
   Htarg=ans["Htarg"]
   pickL <- which.closest(Blim,prod[,"Depletion"])
   Hlim <- prod[pickL,"Harvest"]
   Hmax <- getmax(c(fishery[,"FullH"],Hlim))
   pickyr <- which(fishery[,"FullH"] > 0)
   numval <- dim(fishery)[1]
   par(mai=c(0.45,0.45,0.05,0.05),oma=c(0.0,0,0.0,0.0)) 
   par(cex=0.85, mgp=c(1.35,0.35,0), font.axis=fnt,font=fnt,font.lab=fnt)  
   layout(matrix(c(1,2)),heights=c(3,1))
   plot(fishery[,"SpawnB"],fishery[,"FullH"],type="l",lwd=2,col=1,xlab="Biomass",
        ylab="Annual Harvest Rate",ylim=c(0,Hmax),yaxs="i",xlim=c(0,B0))
   points(fishery[,"SpawnB"],fishery[,"FullH"],pch=16,cex=1.0,col=4)
   points(fishery[2,"SpawnB"],fishery[2,"FullH"],pch=16,cex=1.5,col=3)
   points(fishery[numval,"SpawnB"],fishery[numval,"FullH"],pch=16,cex=1.5,col=2)
   abline(v=c(Blim*B0,Btarg,B0),col=c(2,3,3),lty=2)
   abline(h=c(Htarg,Hlim),col=c(3,2),lty=2)
   text(Btarg,0.05*Hmax,"Btarg",cex=1.0,font=fnt,pos=4)
   text(Blim*B0,0.05*Hmax,"Blim",cex=1.0,font=fnt,pos=4)
   text(0,0.95*Htarg,"Htarg",cex=1.0,font=fnt,pos=4)
   text(0,0.95*Hlim,"Hlim",cex=1.0,font=fnt,pos=4)
   yrs <- fishery[pickyr,"Year"]
   catch <- fishery[pickyr,"Catch"]
   harvest <- fishery[pickyr,"FullH"]
   par(mai=c(0.3,0.45,0.05,0.45)) 
   cmax <- getmax(catch)
   plot(yrs,catch,type="l",lwd=2,col=2,ylab="",xlab="",
        ylim=c(0,cmax),yaxs="i",panel.first=grid(ny=0))
   par(new=TRUE)
   plot(yrs,harvest,type="l",lwd=2,col=4,ylim=c(0,Hmax),yaxt="n",ylab="",
        yaxs="i",xlab="")
   points(yrs[1],harvest[1],pch=16,cex=1.5,col=3)
   points(yrs[numval],harvest[numval],pch=16,cex=1.5,col=2)
   abline(h=c(Htarg),col=c(4),lty=2)
   ym2 <- round(Hmax,2)
   axis(side=4,at=seq(0,ym2,length=3),labels = seq(0,ym2,length=3))
   mtext("Catch (t)",side=2,outer=F,line=1.2,font=fnt,cex=1.0,col=2) 
   mtext("Harvest Rate",side=4,outer=F,line=1.1,font=fnt,cex=1.0,col=4) 
   if (lenfile > 0) {
      outfile <- paste0(getwd(),"/",filename)
      print(outfile)
      dev.off()
   }
   result <- list(B0=B0,Btarg=Btarg,Bmsy=Bmsy,Hmsy=Hmsy,Htarg=Htarg,Hlim=Hlim)
   return(invisible(result))
} # end of aspmphaseplot

#' @title bh calculates the expected Beverton-Holt recruitment
#'
#' @description bh calculate the expected Beverton-Holt stock recruitment level
#'     from the available spawning biomass, the steepness, R0 and B0. This would
#'     be used when fitting a model to data. 
#'
#' @param spb the current spawning or mature biomass
#' @param h the steepness of the Beverton-Holt stock recruitment curve
#' @param R0 the unfished average recruitment level
#' @param B0 the unfished spawning biomass.
#'
#' @return the expected Beverton-Holt recruitment level, a real number in the 
#'     linear scale
#' @export
#'
#' @examples
#' rec <- bh(10000,0.75,1500000,30000)
#' print(rec)   # should be 1285714
#' bh(spb=30000,h=0.75,R0=1500000,B0=30000)  # should be 1500000
bh <- function(spb,h,R0,B0) {
   a <- (4 * h * R0)/(5 * h - 1)
   b <- (B0 * (1 - h))/(5 * h -1)
   return((a * spb)/(b + spb))
} # end of bh


#' @title doDepletion - depletes the stock to the declared initdep level
#'
#' @description doDepletion - depletes the stock to the declared initdep level
#'     There is a printout to the screen of the final outcome. The procedure
#'     assumes you have the unfished NaA in column 1 of the stock$NaA object.
#'     The depletion is arrived at by literally searching for the first
#'     constant harvest rate that leads to the required depletion.
#'
#' @param inR0 the nominal scale unfished recruitment
#' @param indepl the nominal scale initial depletion level to be searched for
#' @param inprops the props object from the data object
#' @param inglb the globals object from the data object
#' @param inc the starting value and step in harvest rate used when finding the 
#'     selected depletion level; defaults to 0.02
#' @param Numyrs the number of years of fishing while searching for the
#'     desired depletion. Some number between 40 - 50 seems to work well.
#'     
#' @return A list containing the vector of numbers-at-age at the required 
#'     depletion level, the Fvalue required to achieve that depletion, and
#'     the actual depletion level achieved
#' @export
#' @examples
#' \dontrun{
#' data(westroughy)
#' fish <- westroughy$fish
#' glb <- westroughy$glb
#' props <- westroughy$props
#' dep <- doDepletion(glb$R0,indepl=0.4,props,glb,inc=0.02)
#' print(dep)
#' }  
doDepletion <- function(inR0,indepl,inprops,inglb,inc=0.02,Numyrs=50) {
#  inR0=glb$R0;indepl=0.4;inprops=props;inglb=glb;inc=0.02;Numyrs=50
   maxage <- inglb$maxage
   Nages <- inglb$nages
   M <- inglb$M
   steep <- inglb$steep
   maa <- inprops$maa
   waa <- inprops$waa
   sela <- inprops$sela
   Nt <- matrix(0,nrow=(maxage+1),ncol=Numyrs,
                dimnames=list(seq(0,maxage,1),seq(0,(Numyrs-1),1)))
   Frange <- seq(0.02,2*M,inc)
   NF <- length(Frange)
   unfish <- unfished(inglb,inprops,inR0)
   B0 <- unfish$B0
   naa <- unfish$N0
   for (hr in 1:NF) {  # hr=1
      dHarv <- Frange[hr]  # apply continually increasing H
      # Fy <- -log(1-dHarv)
      SpawnB <- B0
      Nt[,1] <- naa
      year <- 2     # instead of zero to allow for indexing matrices
      repeat {
         Nt[1,year] <- bh(SpawnB,steep,inR0,B0)  #  Age 0
         Nt[2:(Nages-1),year] <- Nt[1:(Nages-2),(year-1)] * exp(-M/2)
         Nt[Nages,year] <- (Nt[(Nages-1),(year-1)]*exp(-M/2)) + (Nt[(Nages),(year-1)]*exp(-M/2))
         # Now do the fishing
         Nt[,year] <- (Nt[,year] * (1-(sela * dHarv))) * exp(-M/2)
         SpawnB <- SpB(Nt[,year],maa,waa)
         depletion <- SpawnB/B0
         year <- year + 1
         if ((depletion <= indepl) | (year > Numyrs)) break
      }
      if (depletion <= indepl) break
   }
   ans <- list(Ndepl=Nt[,(year-1)],Fval=Frange[hr],depl=depletion,SpawnB=SpawnB)
   return(ans)
} # End of DoDepletion

#' @title dynF describe the ASPM dynamics
#'
#' @description dynF summarizes the dynamics of an Age-Structured
#'     Production Model (ASPM). Fitting the ASPM entails estimating the unfished
#'     recruitment level (R0), which is input as a parameter. In this case 
#'     fishing mortality is implemented as instantaneous rates rather that
#'     as annual harvest rates. The maximize F = 4.0 
#'
#' @param pars the dynamics relies on many parameters sitting in the global
#'     environment in particular ages, nages, maxage, M, fish, nyrs, and
#'     maa, waa, sela, which are contained in props. 'pars' can contain either 
#'     two or three parameters. 1) is the log-transformed average unfished 
#'     recruitment, inR0. 2) is the variability around the index of relative 
#'     abundance (cpue or survey index) during the fitting process, and if is 
#'     present 3) if is present, is the catchability 'q',
#'     which alternatively can be estimated using the closed form. 4) If 
#'     present this would be the log of the initail depletion.
#' @param infish the fish data.frame from readdata or built in dataset
#' @param inglb the glb data.frame from readdata or built in dataset
#' @param inprops the props data.frame from readdata or built in dataset
#' @param waa the character name of the weight-at-age
#' @param maa the character name of the maturity-at-age
#' @param sela the character name of the selectivity-at-age
#' @param full should all outputs from dynamics be given. When fitting the 
#'     model, set this to FALSE. Once fitted, change this to TRUE to get all
#'     the required outputs.
#' @param reps how many loops within the findF function to find each F
#' 
#' @seealso{
#'  \link{dynamicsF}, \link{dynamicsF}, \link{findF}
#' } 
#' 
#' @return a data.frame containing the fishery dynamics according to the input
#'     parameter inR0. Includes Year, Catch, PredC, SpawnB, ExploitB, FullH,
#'     CPUE, PredCE, Deplete, Recruit, FullF.
#' @export
#'
#' @examples
#' \dontrun{
#' data(westroughy)
#' fish <- westroughy$fish
#' glb <- westroughy$glb
#' props <- westroughy$props
#' pars <- c(7.1,-1,-7.7) # logR0, sigCE, estimate avq
#' fishery <- dynF(pars,fish,glb,props) 
#' pars <- c(7.1,-1,-7.7)
#' bestL <- optim(pars,dynF,method="Nelder-Mead",infish=fish,inglb=glb,
#'                inprops=props,control=list(maxit=1000,parscale = c(10,1,10)))
#' str(bestL)
#' codeutils::outfit(bestL,digits=6)
#' out <- dynF(bestL$par,fish,glb,props,full=TRUE)
#' print(round(out$fishery,4)) 
#' }
#' # pars=bestL$par;infish=fish;inglb=glb;inprops=props
#' # waa="waa";maa="maa";sela="sela"
dynF <- function(pars,infish,inglb,inprops,
                 waa="waa",maa="maa",sela="sela",full=FALSE,reps=6) { 
  aaw <- inprops[,waa]
  aam <- inprops[,maa]
  sel <- inprops[,sela]
  epars <- exp(pars)
  R0 <- epars[1]
  sigCE <- epars[2]
  avq <- epars[3]
  B0 <- getB0(R0,inglb,inprops)   
  nyrs <- length(infish[,"year"])
  nyrs1 <- nyrs + 1
  nages <- inglb$nages
  maxage <- inglb$maxage
  Nt <- matrix(0,nrow=nages,ncol=nyrs1,dimnames=list(0:maxage,0:nyrs))
  columns <- c("year","catch","predC","spawnB","exploitB","cpue",
               "predCE","deplete","recruit","fullF","fullH")
  fishery <- matrix(NA,nrow=nyrs1,ncol=length(columns),
                    dimnames=list(0:nyrs,columns))
  fishery[,"year"] <- c((infish$year[1]-1),infish$year)
  fishery[,"catch"] <- c(NA,infish$catch)
  fishery[,"cpue"] <- c(NA,infish$cpue)
  fishery[1,"recruit"] <- R0  
  catch <- fishery[,"catch"]
  M <- inglb$M
  surv <- exp(-M)
  Nt[1,1] <- R0   # get unfished (yr=1) Numbers-at-age index 0 - maxage
  for (age in 1:(maxage-1)) Nt[age+1,1] <- Nt[age,1] * surv
  Nt[maxage+1,1] <- (Nt[maxage,1] * surv)/(1-surv)
  for (yr in 2:nyrs1) {  # yr=2
    spb <- SpB(Nt[,(yr-1)],aam,aaw)
    exb <- ExB(Nt[,(yr-1)],sel,aaw)
    Nt[1,yr] <- bh(spb,inglb$steep,R0,B0)
    fishery[yr,"recruit"] <- Nt[1,yr]
    yrF <- findFs(catch[yr],Nt[,yr-1],sel,aaw,M,reps=reps)
    # yrF <- optimize(matchC,interval=c(0,4.0),M=M,cyr=catch[yr],
    #                 Nyr=Nt[,yr-1],sel=sel,waa=aaw,maximum=FALSE,
    #                 tol = 1e-09)$minimum
    sF <- sel * yrF
    psF <- sF[2:nages]
    msF <- sF[nages]
    catchN <- (sF/(M + sF)) * (Nt[,yr-1] * (1 - exp(-(M + sF))))  
    fishery[yr,"predC"] <- sum(catchN * aaw)/1000.0
    fishery[yr,"fullF"] <- yrF
    Nt[2:nages,yr] <- (Nt[1:(nages-1),(yr-1)] * exp(-(M + psF)))
    Nt[nages,yr] <- Nt[nages,yr] + (Nt[nages,yr-1] * exp(-(M + msF)))
    fishery[(yr-1),c("spawnB","exploitB")] <- c(spb,exb)
    fishery[yr,"fullH"] <- 1-exp(-yrF)
  }
  spb <- SpB(Nt[,yr],aam,aaw)   # to complete final year
  exb <- ExB(Nt[,yr],sel,aaw)
  fishery[yr,c("spawnB","exploitB")] <- c(spb,exb)
  fishery[,"deplete"] <- fishery[,"spawnB"]/B0
  ExpB <- fishery[1:nyrs,"exploitB"]
  fishery[2:nyrs1,"predCE"] <- ExpB * avq
  pick <- which(fishery[,"cpue"] > 0)
  f <- -sum(dnorm(log(fishery[pick,"cpue"]),log(fishery[pick,"predCE"]),
                  sigCE,log=TRUE),na.rm=TRUE)
  if (full) {
    absdiffC <- sum(abs(fishery[,"catch"] - fishery[,"predC"]),na.rm=TRUE)
    out <- list(fishery=as.data.frame(fishery),Nt=Nt,B0=B0,R0=R0,avq=avq,
                LL=f,diffC=absdiffC)
    return(out)
  } else {
    return(f)
  }
} # end of dynF

#' @title dynamicsF describe the ASPM dynamics using instantaneous F
#'
#' @description dynamicsF summarizes the dynamics of an Age-Structured
#'     Production Model (ASPM). Fitting the ASPM entails estimating the unfished
#'     recruitment level (R0), which is input as a parameter. In this case 
#'     fishing mortality is implemented as instantaneous rates rather that
#'     as annual harvest rates. The maximize F = 4.0 
#'
#' @param pars the dynamics relies on many parameters sitting in the global
#'     environment in particular ages, nages, maxage, M, fish, nyrs, and
#'     maa, waa, sela, which are contained in props. 'pars' can contain either 
#'     two or three parameters. 1) is the log-transformed average unfished 
#'     recruitment, inR0. 2) is the variability around the index of relative 
#'     abundance (cpue or survey index) during the fitting process, and if is 
#'     present 3) if is present, is the catchability 'q',
#'     which alternatively can be estimated using the closed form. 4) If 
#'     present this would be the log of the initail depletion.
#' @param infish the fish data.frame from readdata or built in dataset
#' @param inglb the glb data.frame from readdata or built in dataset
#' @param inprops the props data.frame from readdata or built in dataset
#' @param waa the character name of the weight-at-age
#' @param maa the character name of the maturity-at-age
#' @param sela the character name of the selectivity-at-age
#' @param full should all outputs from dynamics be given. When fitting the 
#'     model, set this to FALSE. Once fitted, change this to TRUE to get all
#'     the required outputs.
#' 
#' @seealso{
#'  \link{dynamicsH}
#' } 
#' 
#' @return a data.frame containing the fishery dynamics according to the input
#'     parameter inR0. Includes Year, Catch, PredC, SpawnB, ExploitB, FullH,
#'     CPUE, PredCE, Deplete, Recruit, FullF.
#' @export
#'
#' @examples
#' \dontrun{
#' data(westroughy)
#' fish <- westroughy$fish
#' glb <- westroughy$glb
#' props <- westroughy$props
#' pars <- c(7.0,0.3,-7.7) # logR0, sigCE, estimate avq
#' fishery <- dynamicsF(pars,fish,glb,props,full=TRUE) 
#' bestL <- optim(pars,dynamicsF,method="Nelder-Mead",infish=fish,inglb=glb,
#'                inprops=props,control=list(maxit = 1000,parscale = c(10,0.1)))
#' str(bestL)
#' fishery <- dynamicsF(bestL$par,fish,glb,props)
#' print(round(fishery,4)) 
#' }
dynamicsF <- function(pars,infish,inglb,inprops,
                      waa="waa",maa="maa",sela="sela",full=FALSE) { 
  # pars=pars;infish=fish;inglb=glb;inprops=props;waa="waa";maa="maa";sela="sela";maxF=1.0
  aaw <- inprops[,waa]
  aam <- inprops[,maa]
  sel <- inprops[,sela]
  epars <- exp(pars)
  R0 <- epars[1]
  sigCE <- epars[2]
  avq <- epars[3]
  B0 <- getB0(R0,inglb,inprops)   
  nyrs <- length(infish[,"year"])
  nages <- inglb$nages
  maxage <- inglb$maxage
  Nt <- matrix(0,nrow=nages,ncol=(nyrs+1),dimnames=list(0:(nages-1),0:nyrs))
  columns <- c("year","catch","predC","spawnB","exploitB","cpue",
               "predCE","deplete","recruit","fullF","fullH")
  fishery <- matrix(NA,nrow=(nyrs+1),ncol=length(columns),
                    dimnames=list(0:nyrs,columns))
  fishery[,"year"] <- c((infish$year[1]-1),infish$year)
  fishery[,"catch"] <- c(NA,infish$catch)
  fishery[,"cpue"] <- c(NA,infish$cpue)
  fishery[1,"recruit"] <- R0  
  catch <- fishery[,"catch"]
  M <- inglb$M
  surv <- exp(-M)
  Nt[,1] <- R0   # get unfished (yr=1) Numbers-at-age
  for (age in 1:(maxage-1)) Nt[age+1,1] <- Nt[age,1] * surv
  Nt[maxage+1,1] <- (Nt[maxage,1] * surv)/(1-surv)
  for (yr in 2:(nyrs+1)) {  # yr=2
    spb <- SpB(Nt[,(yr-1)],aam,aaw)
    exb <- ExB(Nt[,(yr-1)],sel,aaw)
    Nt[1,yr] <- bh(spb,inglb$steep,R0,B0)
    fishery[yr,"recruit"] <- Nt[1,yr]
    yrF <- optimize(matchC,interval=c(0,4.0),M=M,cyr=catch[yr],
                    Nyr=Nt[,yr-1],sel=sel,waa=aaw,maximum=FALSE,
                    tol = 1e-09)$minimum
    sF <- sel * yrF
    psF <- sF[2:nages]
    msF <- sF[nages]
    catchN <- (sF/(M + sF)) * (Nt[,yr-1] * (1 - exp(-(M + sF))))  
    fishery[yr,"predC"] <- sum(catchN * aaw)/1000.0
    fishery[yr,"fullF"] <- yrF
    Nt[2:nages,yr] <- (Nt[1:(nages-1),(yr-1)] * exp(-(M + psF)))
    Nt[nages,yr] <- Nt[nages,yr] + (Nt[nages,yr-1] * exp(-(M + msF)))
    fishery[(yr-1),c("spawnB","exploitB")] <- c(spb,exb)
    fishery[yr,"fullH"] <- 1-exp(-yrF)
  }
  spb <- SpB(Nt[,yr],aam,aaw)   # to complete final year
  exb <- ExB(Nt[,yr],sel,aaw)
  fishery[yr,c("spawnB","exploitB")] <- c(spb,exb)
  fishery[,"deplete"] <- fishery[,"spawnB"]/B0
  ExpB <- fishery[1:nyrs,"exploitB"]
  fishery[2:(nyrs+1),"predCE"] <- ExpB * avq
  pick <- which(fishery[,"cpue"] > 0)
  f <- -sum(dnorm(log(fishery[pick,"cpue"]),log(fishery[pick,"predCE"]),
                  sigCE,log=TRUE),na.rm=TRUE)
  if (full) {
    absdiffC <- sum(abs(fishery[,"catch"] - fishery[,"predC"]),na.rm=TRUE)
    out <- list(fishery=as.data.frame(fishery),Nt=Nt,B0=B0,R0=R0,avq=avq,
                LL=f,diffC=absdiffC)
    return(out)
  } else {
    return(f)
  }
} # end of dynamicsF


#' @title dynamicsH describe the ASPM dynamics using annual harvest rates
#'
#' @description dynamicsH summarizes the dynamics of the Age-Structured
#'     Production Model (ASPM) in which the catches are represented as annual
#'     harvest rates. Fitting the ASPM entails estimating the unfished
#'     recruitment level (R0), which is input as a parameter. There may be 
#'     other parameters as described in the pars section. 
#'
#' @param pars the dynamics relies on many parameters sitting in the function's
#'     environment, these are ages, nages, maxage, M, maa, waa, sela, fish,
#'     and nyrs. 'pars' can contain either two to three log-transformed 
#'     parameters: 1) is the log of average unfished recruitment, inR0. 2) is 
#'     the variability around the index of relative abundance (cpue), not used 
#'     directly in the dynamics but rather in the estimation of the likelihoods 
#'     during the fitting process, 3) if is present, is the catchability 'q',
#'     which alternatively can be estimated using the closed form. 4) If 
#'     present this would be the log of the initail depletion.
#' @param infish the fish data.frame from readdata or an internal dataset
#' @param inglb the glb data.frame from readdata or an internal dataset
#' @param inprops the props data.frame from readdata or an internal dataset
#' @param waa the character name of the weight-at-age
#' @param maa the character name of the maturity-at-age
#' @param sela the character name of the selectivity-at-age
#' @param full should all outputs from dynamics be given. When fitting the 
#'     model, set this to FALSE. Once fitted, change this to TRUE to get all
#'     the required outputs.
#'     
#' @seealso{
#'  \link{dynamicsF}
#' } 
#' 
#' @return a data.frame containing the fishery dynamics according to the input
#'     parameter inR0. In particular it includes the Catch and PredC, and the
#'     CPUE and PredCE, which can be used in a maximum likelihood context.
#' @export
#'
#' @examples
#' data(westroughy)
#' fish <- westroughy$fish
#' glb <- westroughy$glb
#' props <- westroughy$props
#' pars=c(7.064324,-1.257487,-7.694248)
#' bestL <- optim(pars,dynamicsH,method="Nelder-Mead",infish=fish,inglb=glb,
#'                inprops=props,control=list(maxit=1000,parscale = c(10,1,10)))
#' str(bestL)
#' codeutils::outfit(bestL,digits=6)
#' out <- dynamicsH(bestL$par,infish=fish,inglb=glb,inprops=props,full=TRUE)
#' print(round(out$fishery,4))
#' # pars=c(7.064324,-1.257487,-7.694248);infish=fish;inglb=glb;
#' # inprops=props;full=TRUE; waa="waa";maa="maa";sela="sela"
dynamicsH <- function(pars,infish,inglb,inprops,
                      waa="waa",maa="maa",sela="sela",full=FALSE) {  
  aaw <- inprops[,waa]  # setup the model structure 
  aam <- inprops[,maa]
  sel <- inprops[,sela]
  epars <- exp(pars) # back-transform the parameters
  R0 <- epars[1]
  sigCE <- epars[2] # used in the estimation of the -ve Log-Likelihood
  B0 <- getB0(R0,inglb,inprops)  
  nyrs <- length(infish[,"year"])
  nages <- inglb$nages
  maxage <- inglb$maxage
  Nt <- matrix(0,nrow=nages,ncol=(nyrs+1),dimnames=list(0:maxage,0:nyrs))
  columns <- c("year","catch","predC","spawnB","exploitB","cpue",
               "predCE","deplete","recruit","fullF","fullH")
  fishery <- matrix(NA,nrow=(nyrs+1),ncol=length(columns),
                    dimnames=list(0:nyrs,columns))
  fishery[,"year"] <- c((infish$year[1]-1),infish$year) # include year 0
  fishery[,"catch"] <- c(NA,infish$catch)
  fishery[,"cpue"] <- c(NA,infish$cpue)
  hS <- exp(-inglb$M/2)  # survivorship from half natural mortality
  surv <- exp(-inglb$M)
  # calculate unfished numbers-at-age given inR0 over next 3 lines
  Nt[,1] <- R0   # Nt column index 1 = year 0
  for (age in 1:(maxage-1)) Nt[age+1,1] <- Nt[age,1] * surv
  Nt[maxage+1,1] <- (Nt[maxage,1] * surv)/(1-surv)
  for (yr in 2:(nyrs+1)) {  # yr=2# do dynamics for yr 2 - nyrs+1 columns
    spb <- SpB(Nt[,(yr-1)],aam,aaw)
    exb <- ExB(Nt[,(yr-1)]*hS,sel,aaw)
    Nt[1,yr] <- bh(spb,inglb$steep,R0,B0)
    fishery[yr,"recruit"] <- Nt[1,yr]
    harvest <- min((fishery[yr,"catch"]/exb),0.975) # sets upper limit to H
    hrate <- sel * harvest
    fishery[yr,"fullH"] <- harvest
    Ct <- (Nt[,(yr-1)] * hS) * hrate
    Nt[2:nages,yr] <- ((Nt[1:(nages-1),(yr-1)] * hS) - Ct[1:(nages-1)]) * hS
    Nt[nages,yr] <- Nt[nages,yr] + ((Nt[nages,yr-1] * hS) - Ct[nages]) * hS
    fishery[(yr-1),c("spawnB","exploitB")] <- c(spb,exb) # save in 'fishery'
    fishery[yr,"predC"] <- sum(Ct * aaw)/1000
  }
  fishery[,"fullF"] <- -log(1 - fishery[,"fullH"]) # translate annual H to F
  spb <- SpB(Nt[,yr],aam,aaw)   # complete final year
  exb <- ExB(Nt[,yr]*hS,sel,aaw)
  fishery[yr,c("spawnB","exploitB")] <- c(spb,exb)
  fishery[,"deplete"] <- fishery[,"spawnB"]/B0
  ExpB <- fishery[1:nyrs,"exploitB"]  # calculate predicted CPUE
  avq <- ifelse(length(epars) > 2, epars[3], # estimate or closed form
              exp(mean(log(infish$cpue/fishery[1:nyrs,"exploitB"]),na.rm=TRUE)))
  fishery[2:(nyrs+1),"predCE"] <- ExpB * avq
  pick <- which(fishery[,"cpue"] > 0)  # calculate -ve log-likelihood
  f <- -sum(dnorm(log(fishery[pick,"cpue"]),log(fishery[pick,"predCE"]),
                  sigCE,log=TRUE),na.rm=TRUE)  # uses log-Normal errors
  if (full) {  # output everything or only -ve log-likelihood
    absdiffC <- sum(abs(fishery[,"catch"] - fishery[,"predC"]),na.rm=TRUE)
    out <- list(fishery=as.data.frame(fishery),Nt=Nt,B0=B0,R0=R0,avq=avq,
                LL=f,diffC=absdiffC)
    return(out)
  } else {
    return(f)
  }
} # end of dynamicsH

#' @title ExB calculate exploitable biomass from numbers-at-age
#'
#' @description ExB calculates the spawning biomass from a vector of
#'     numbers-at-age, selectivity-at-age, and Weight-at-age.
#'
#' @param invect the numbers-at-age as a vector
#' @param SelA selectivity-at-age vector
#' @param WeightA weight-at-age vector as kilograms
#' 
#' @return ExB a scalar as tonnes.
#' @export
#' @examples
#' \dontrun{
#' data(fishdat)
#' str(fishdat)
#' glb <- fishdat$glb
#' fish <- fishdat$fish
#' props <- fishdat$props
#' unfish <- unfished(glb,props,glb$R0)
#' N1 <- unfish$N0 * exp(-glb$M/2)  # no growth
#' ExB(N1,props$sela,props$waa)  # should be 17999.6
#' }
ExB <- function(invect, SelA, WeightA) {
  ans <- sum(SelA * WeightA * invect)/1000.0
  return(ans)
}


#' @title findF uses an iterative strategy for finding each year's F value
#' 
#' @description findF is used when searching for the instantaneous F that will
#'     produce the required catch in a given year in an age-structured model. It
#'     needs generalizing to operate with multiple gears.
#'
#' @param cyr the known catch in a given year
#' @param Nyr the numbers at size at the start of the given year or end of the 
#'     year before
#' @param sel the selectivity of the fishing gear
#' @param aaw the weight-at-age
#' @param M the instantaneous natural mortality rate
#' @param Fmax the limit on themaximum F allowed, default = 3.0 ~H = 0.95
#' @param reps how many internal loops to use finding each F, default = 6
#'
#' @returns the fully selected fishing mortality rate
#' @export
#'
#' @examples
#' print("wait on example data sets")
findF <- function(cyr,Nyr,sel,aaw,M,Fmax=3.0,reps=8) {
  # cyr=catch[yr]; Nyr <- Nt[,yr-1]; sel=sel; aaw=aaw  
  wata <- aaw/1000
  Byr <- sum((Nyr*sel*wata))
  temp1 <- cyr / Byr
  join1 <- 1/(1 + exp(30*(temp1 - 0.95)))
  tempyr <- (join1 * temp1) + (0.95 * (1 - join1))
  fFyr <- -log(1 - tempyr)
  sF <- sel * fFyr
  Zyr <- sF + M
  predCyr <- sum((sF/Zyr) * (wata * Nyr) * (1 - exp(-Zyr)))
  for (i in 1:reps) {
    Zadj <- cyr/(predCyr + 0.00001)
    cat(Zadj,"\n")
    fFyr <- Zadj * fFyr
    sF <- sel * fFyr
    Zyr <- sF + M
    predCyr <- sum((sF/Zyr) * (wata * Nyr) * (1 - exp(-Zyr)))
  }
  cat("\n")
  return(fFyr) 
} # end of findF

#' @title fitASPM fits an age-structured production model
#'
#' @description fitASPM fits an age-structured production model that can have 
#'     up to three parameters, R0 the unfished recruitment level, se the 
#'     variation around the estimated CPUE, and avq the catchability 
#'     coefficient. The inverse Hessian, when using log-likelihoods, provides
#'     access to the variance-covariance matrix, which in turn can be used
#'     to generate standard errors on the model parameters. This is used when
#'     using asypmtotic errors to caharacterize uncertainty. See section 6.5
#'     in Haddon (2021).
#'     
#' @references Haddon, M. (2021) Using R for Modelling and Quantitative 
#'     Methods in Fisheries, CRC Press / Chapman & Hall/ Boca Raton 337p.
#'     ISBN: 9780367469894
#'
#' @param initpar a vector of 2 or 3 numbers that are the initial parameter
#'     values given to the estimate of logR0, and the estimate of the variation
#'     around the CPUE data that the model is to be fitted to, and finally, the
#'     catchability coefficient.
#' @param minfun the dynamics function outputning a -ve log-likelihood 
#' @param infish the fish data.frame from readdata or built in dataset
#' @param inglb the glb data.frame from readdata or built in dataset
#' @param inprops the props data.frame from readdata or built in dataset
#' @param gradtol default = 1e-04 a positive scalar giving the tolerance at 
#'     which the scaled gradient is considered close enough to zero to 
#'     terminate the algorithm. The scaled gradient is a measure of the 
#'     relative change in f in each direction p[i] divided by the relative 
#'     change in p[i] (copied frm nlm help).
#' @param stepmax default = 0.1 a positive scalar which gives the maximum 
#'     allowable scaled step length. stepmax is used to prevent steps which 
#'     would cause the optimization function to overflow, to prevent the 
#'     algorithm from leaving the area of interest in parameter space, or to 
#'     detect divergence in the algorithm. stepmax would be chosen small enough 
#'     to prevent the first two of these occurrences, but should be larger than 
#'     any anticipated reasonable step. (copied frm nlm help).
#' @param steptol default = 1e-08 A positive scalar providing the minimum 
#'     allowable relative step length. (copied frm nlm help).
#' @param hessian should hessian be estimated at the optimum, default = FALSE
#' 
#' @return a list containing the optimal output
#' @export
#'
#' @examples
#' data("westroughy")
#' fish <- westroughy$fish
#' glb <- westroughy$glb
#' props <- westroughy$props
#' pars <- c(7,0.3,-7.7)
#' bestL <- fitASPM(pars,dynamicsH,infish=fish,inglb=glb,inprops=props,
#'                  hessian=TRUE)
#' bestL
#' out <- dynamicsH(bestL$estimate,fish,glb,props,full=TRUE)
#' round(out$fishery,4)
fitASPM <- function(initpar,minfun,infish,inglb,inprops,gradtol=1e-04,
                    stepmax=0.1,steptol=1e-08,hessian=FALSE) { 
  paramscale = magnitude(initpar)
  bestL <- optim(par=initpar,fn=minfun,method="Nelder-Mead",
                 infish=infish,inglb=inglb,inprops=inprops,
                 control=list(maxit = 1000, parscale = paramscale))
  bestL <- nlm(f=minfun,p=bestL$par,hessian=hessian,
               infish=infish,inglb=inglb,inprops=inprops,
               gradtol=gradtol,
               stepmax=stepmax,steptol=steptol,
               iterlim=1000)
  return(bestL)
} # end of fitASPM

#' @title getB0 calculates the B0 from biological properties and R0
#'
#' @description getB0 calculates the B0 from biological properties of M,
#'     maxage, maa, and waa, plus the input of R0 on the nominal scale (the 
#'     hypothetical unfished recruitment level). This is used in the 'dynamics'
#'     function.
#'
#' @param inR0 the estimate of unfished recruitment
#' @param inglb the glb data.frame from readdata or built in dataset
#' @param inprops the props data.frame from readdata or built in dataset
#' 
#' @return a single number that is the estimate of B0
#' @export
#'
#' @examples
#' \dontrun{
#' data(fishdat)
#' glb <- fishdat$glb
#' props <- fishdat$props
#' getB0(1275000,glb,props) # shoud give 19429.76
#' getB0(1000000,glb,props) # should give 15239.03
#' }
getB0 <- function(inR0,inglb,inprops) { # assumes glb inR0 = par["R0"]
   maxage <- inglb$maxage
   surv <- exp(-inglb$M)
   Nt <- numeric(maxage+1)
   Nt[1] <- 1  # calculate numbers-at-age per recruit
   for (age in 1:(maxage-1)) Nt[age+1] <- Nt[age] * surv
   Nt[maxage+1] <- (Nt[maxage] * surv)/(1-surv)
   A0 <-  sum(inprops$maa * inprops$waa * Nt)/1000.0
   B0 <- inR0 * A0
   return(B0)
}  # end of getB0

#' @title getProduction estimates MSY for ASPM using annual harvest rates
#'
#' @description getProduction takes the optimum estimate of R0 from ASPM and
#'     estimates the production curve, from which it gains the MSY, Bmsy, Hmsy,
#'     and Dmsy (depletion at MSY). It generate the production curve by stepping
#'     through the dynamics of the fishery over a series of constant harvest
#'     rates for however many iterations of nyr required.
#'
#' @param inR0 the optimum estimate of R0 from ASPM on nominal scale
#' @param infish the fish data.frame from readdata or built in dataset
#' @param inglb the glb data.frame from readdata or built in dataset
#' @param inprops the props data.frame from readdata or built in dataset#'
#' @param Hrg the desired sequence of harvest rates to be used to define the
#'     production curve. The default is c(0.025,0.4,0.025), but it is a good
#'     idea to plot teh harvest rate against yield and manually adjust the
#'     sequence values to focus attention and detail on where the MSY falls.
#' @param nyr The number of years making up a single iteration while searching
#'     for the equilibrium yield. default = 50.
#' @param maxiter defaults to 3. Determines how many runs through the nyr
#'     steps are conducted to reach equilibrium. One can trial different values
#'     but 3 is a reasonable compromise. It is best to test different maxiter
#'     values to ensure equilibria are reached in each yield estimate.
#'
#' @return a data.frame containing the sequence of harvest rates, the related
#'     spawning biomass, the exploitable biomass, the yield, and depletion.
#' @export
#'
#' @examples
#' \dontrun{
#'   data("westroughy")
#'   fish <- westroughy$fish
#'   glb <- westroughy$glb
#'   props <- westroughy$props
#'   pars <- c(7.0,0.3)
#'   bestL <- nlminb(start=pars,objective=dynamicsH,infish=fish,inglb=glb,
#'                   inprops = props,
#'                control=list(eval.max=500,iter.max=300,trace=0,rel.tol=1e-08))
#'   prod <- getProduction(exp(bestL$par[1]),infish=fish,inglb=glb,
#'                         inprops=props,Hrg=c(0.0005,0.07,0.0005),nyr=100)
#'   prod[78:100,]
#'   prod[84:89,]  # check out H = 0.0425 giving the MSY.
#' }
getProduction <- function(inR0,infish,inglb,inprops,Hrg=c(0.025,0.4,0.025),
                          nyr=50,maxiter=3) {
  maxage <- inglb$maxage; M <- inglb$M;  steep <- inglb$steep
  nages <- inglb$nages; aam <- inprops$maa; aaw <- inprops$waa
  sel <- inprops$sela
  surv <- exp(-M)
  hS <- exp(-M/2)
  B0 <- getB0(inR0,inglb,inprops)   
  getyield <- function(NAA,harvest) {
    for (iter in 1:maxiter) {  # iterate until equilibrium yield reached
      for (yr in 2:nyr) {  # yr=nyrs+1
        spb <- SpB(NAA[,(yr-1)],aam,aaw)
        NAA[1,yr] <- bh(spb,steep,inR0,B0)
        Ct <- (NAA[,(yr-1)] * hS) * (harvest * sel)
        NAA[2:nages,yr] <- ((NAA[1:(nages-1),(yr-1)] * hS) - Ct[1:(nages-1)]) * hS
        NAA[nages,yr] <- NAA[nages,yr] + ((NAA[nages,yr-1] * hS) - Ct[nages]) * hS
      }
      newcatch <- sum(Ct * aaw)/1000
      NAA[,1] <- NAA[,nyr]
    }
    spb <- SpB(NAA[,nyr],aam,aaw)
    exb <- ExB(NAA[,nyr]*hS,sel,aaw)
    return(c(spb=spb,exb=exb,yield=newcatch))
  } # end of getyield
  Nt <- matrix(0,nrow=(maxage+1),ncol=nyr,dimnames=list(0:maxage,1:nyr))
  hrange <- seq(Hrg[1],Hrg[2],Hrg[3])
  nH <- length(hrange)
  columns <- c("Harvest","SpawnB","ExploitB","Yield","Depletion")
  production <- matrix(NA,nrow=(nH+1),ncol=length(columns),
                       dimnames=list(c(0,hrange),columns))
  production[,"Harvest"] <- c(0,hrange)
  Nt[1,1] <- inR0
  for (age in 1:(maxage-1)) Nt[age+1,1] <- Nt[age,1] * surv
  Nt[maxage+1,1] <- (Nt[maxage,1] * surv)/(1-surv)
  production[1,"SpawnB"] <- SpB(Nt,aam,aaw)
  production[1,"ExploitB"] <- ExB(Nt*hS,sel,aaw)
  for (hnum in 1:nH)
    production[(hnum+1),2:4] <- getyield(Nt,hrange[hnum])
  production[,"Depletion"] <- production[,"SpawnB"]/production[1,"SpawnB"]
  return(as.data.frame(production))
} # end of getProduction


#' @title logist Logistic selectivity function with knifeedge option
#'
#' @description logist calculates a Logistic curve that can be used as a
#'     selectivity function, or maturity curve, of wherever a logistic is
#'     required. This version uses the logistic function
#'     1/(1+exp(-log(19.0)*(lens-inL50)/(inL95-inL50))),
#'     which explicitly defines the SM50 and uses SM95 as the second parameter.
#' @param inL50 is the length at 50 percent selection/maturity/whatever
#' @param delta is the difference in selection/maturity/whatever between
#'     inL50 and inL95
#' @param depend a vector of lengths/ages for which the logistic value will be
#'     calculated.
#' @param knifeedge defaults to 0. If knifeedge is set to a particular length or
#'     age then the logistic value <= the value of knifeedge is set to
#'     zero, which is essentially knife-edge. Allows for knife-edge selectivity
#' @return A vector of length(depend) containing the predicted logistic values
#' @export
#' @examples
#' \dontrun{
#' in50 <- 100.0
#' deltaS <- 8.0
#' lens <- seq(2,210,2)
#' select <- logist(inL50=in50,delta=deltaS,depend=lens)
#' selectk <- logist(in50,deltaS,lens,knifeedge=105)
#' round(cbind(lens[35:70],select[35:70],selectk[35:70]),5)
#' } 
logist <- function(inL50,delta,depend,knifeedge=0) {
   ans <- 1/(1+exp(-log(19.0)*(depend-inL50)/(delta)))
   if (knifeedge > 0) {
      pick <- which(depend <= knifeedge)
      if (length(pick) > 0) ans[pick] <- 0.0
   }
   return(ans)
}

#' @title MaA an alternative logistic function commonly used for maturity
#'
#' @description MaA - the logistic function exp(a+bxdepend)/(1+exp(a+bxdepend)),
#'     which can also be expressed as 1/(1+(1/exp(a + b x depend))). This has
#'     the property that the SM50 = -a/b and the interquartile distance is
#'     2.Ln(3)/b.
#' @param ina is the intercept of the exponential function
#' @param inb is the gradient of the exponential function
#' @param depend is a vector of lengths/ages for which the logistic maturity
#'     value will be calculated
#' @return A vector of length(depend) containing the predicted maturity values
#'
#' @export
#' @examples
#' a <- -14.383
#' b <- 0.146017
#' lens <- seq(2,210,2)
#' round(MaA(a,b,depend=lens),5) # length based
#' round(MaA(-2.5,0.95,0:25),5)   # age based
MaA <- function(ina,inb,depend) {
   ans <- exp(ina+inb*depend)/(1+exp(ina+inb*depend))
   return(ans)
}

#' @title matchC gives difference between predicted catch and observed catch
#' 
#' @description matchC for a given instantaneous fishing mortality rate
#'     calculates the difference between the predicted and observed catches. 
#'     matchC uses the Baranov catch equation to account, when usung the
#'     instantaneous rates, for the fish that died due to fishing that would
#'     have died anyway, due to natural mortality. The equation is:
#'     Ct = (F/(F + M)) x (1 - exp(-(M + F))) x Bexp. Of course, one needs to
#'     account for selectivity aby age and F by year, so check the equation
#'     in the mfmur book. matchC is used by the optimize function which uses 
#'     'a combination of a golden section search and successive parabolic 
#'     interpolation' to minimize the difference between the absolute 
#'     difference between the catch and the predicted catch.
#'
#' @param f the trial predicted instantaneous fishing mortality
#' @param M the instantaneous natural mortality value
#' @param cyr the catch in yr t
#' @param Nyr the numbers at age at the end of yr-1
#' @param sel the selectivity of the fishing gear in the model
#' @param waa the weight -at- age for the species
#'
#' @return the absolute difference between the predicted and observed catches
#' @export
#'
#' @examples
#' data(westroughy)
#' fish <- westroughy$fish
#' glb <- westroughy$glb
#' props <- westroughy$props
#' sel <- props$sela
#' pars <- c(7.15,-1,-7.7) # logR0, sigCE, estimate avq
#' out <- dynamicsF(pars,fish,glb,props,full=TRUE)  
#' matchC(f=0.3,M=0.05,cyr=15340,Nyr=out$Nt[,2],sel=sel,props$waa)
#' matchC(f=0.35,M=0.05,cyr=15340,Nyr=out$Nt[,2],sel=sel,props$waa)
#' out2 <- optimize(matchC,interval=c(0,1),M=0.05,cyr=3924,Nyr=out$Nt[,2],
#'                 sel=sel,props$waa)
#' out2
#' f <- out2$minimum
#' (397388.66 * (1 - exp(-(0.05 + f))) * f/(0.05 + f))
#' # f=0.2679;M=0.0036;cyr=3924;Nyr=Nt[,yr-1];sel=sel;waa=props$waa
matchC <- function(f,M,cyr,Nyr,sel,waa) {
  sF <- sel * f 
  catchN <- (sF/(M + sF)) * (Nyr * (1 - exp(-(M + sF))))
  out <- abs((sum(catchN * waa)/1000.0) - cyr)
  return(out)
} # en do fmatchC


#' @title plotASPM plots catch, CPUE, Spawning Biomass and Harvest Rate
#' 
#' @description plotASPM after running fitASPM the optimum parameters can be 
#'     put through the dynamics function to generate a dataframe containing
#'     the optimum dynamics. These can be plotted using plotASPM, which plots
#'     out the catches, the Spawning Biomass, the relative CPUE and its fit to
#'     the observed CPUE, and the harvest rate. This routine is still under 
#'     development to include more options.
#'
#' @param infish an object generated by the dynamics function
#' @param CI defaults to NA, if confidence intervals around the cpue have been
#'     obtained using getLNCI, then the resulting matrix will generate 95pc CIs 
#' @param defineplot define the plot size and character outside the plot or
#'     automatically inside. Defaults to TRUE
#' @param target target depletion level. Defaults to 0.48
#' @param usef defines the font to use usef(ont),default = 7 bold times
#' @param png save a png file with the name in 'png', default = "", which
#'     means no file produced
#'
#' @return Nothing, but it does plot six graphs in a single plot.
#' @export
#'
#' @examples
#' \dontrun{
#' data(westroughy)
#' fish <- westroughy$fish
#' glb <- westroughy$glb
#' props <- westroughy$props
#' pars <- c(7.1,-1,-7.7)  # logR0 logceSE logavq
#' bestaspm <- fitASPM(pars,minfun=dynamicsH,infish=fish,
#'                    inglb=glb,inprops=props)
#' out <- dynamicsH(bestaspm$estimate,fish,glb,props,full=TRUE)
#' plotASPM(out$fishery,defineplot=TRUE)
#' ceCI <- getLNCI(out$fishery[,"predCE"],exp(bestspm$estimate[2]))
#' plotASPM(out$fishery,CI=ceCI)
#' }  # infish=fisheryPen; CI=ceCI; defineplot=TRUE; target=0.48; usef=7;png=""
#' # infish=outH$fishery; CI=ceCI;defineplot=TRUE; target=0.48; usef=7;png=""
plotASPM <- function(infish,CI=NA,defineplot=TRUE, target=0.48,usef=7,png="") { 
   if (nchar(png) > 0) defineplot=FALSE
   if (defineplot) { 
      if (names(dev.cur()) %in% c("null device", "RStudioGD"))
         dev.new(width = 7, height = 5.5, noRStudioGD = TRUE)
   }
   par(mfrow=c(3,2),mai=c(0.25,0.4,0.1,0.05),oma=c(0.0,0,0.0,0.0),tck=-0.02) 
   par(cex=0.85,mgp=c(1.35,0.35,0),font.axis=usef,font=usef,font.lab=usef)  
   if (nchar(png) > 0) {
      dev.off()
      graphics.off()
      graphfile <- png
      if (file.exists(graphfile)) file.remove(graphfile)
      png(filename=graphfile,width=210,height=160,units="mm",res=200) 
      par(mfrow=c(3,2),mai=c(0.25,0.4,0.1,0.05),oma=c(0.0,0,0.0,0.0),tck=-0.02) 
      par(cex=0.85,mgp=c(1.35,0.35,0),font.axis=usef,font=usef,font.lab=usef)
   }
   yrs <- infish$year
   nyrs <- length(yrs)
   # plot catches
   ymax <- getmax(infish$catch)
   plot(yrs,infish$catch,type="l",lwd=2,ylim=c(0,ymax),yaxs="i",xlab="",
        panel.first=grid(),ylab="Catch (t)")
   # plot Spawning Biomass
   ymax <- getmax(infish$spawnB)
   plot(yrs,infish$spawnB,type="l",lwd=2,ylim=c(0,ymax),yaxs="i",xlab="",
        panel.first=grid(),ylab="Spawning Biomass (t)")
   # plot CPUE
   ymax <- getmax(c(infish$cpue,infish$predCE))
   if ("matrix" %in% class(CI)) ymax <- getmax(CI[,"upper"]) 
   plot(yrs,infish$cpue,type="p",pch=16,col=2,cex=1.0,ylim=c(0,ymax),yaxs="i",
        xlab="",panel.first=grid(),ylab="Relative Abundance Index")
   lines(yrs,infish$predCE,lwd=2,col=1)
   if ("matrix" %in% class(CI)) {
      segments(x0=yrs,y0=CI[,1],x1=yrs,y1=CI[,3],lwd=1,col=4)
   }
   # plot harvest rate
   ymax <- getmax(infish$fullH)
   plot(yrs,infish$fullH,type="l",lwd=2,ylim=c(0,ymax),yaxs="i",xlab="",
        panel.first=grid(),ylab="Annual Harvest Rate")
   # plot the residuals
   pickCE <- which(infish[,"cpue"] > 0)
   resid <- infish[pickCE,"cpue"]/infish[pickCE,"predCE"]
   nresid <- length(resid)
   ymax <- getmax(resid);    ymin <- getmin(resid,mult=1.1)
   plot(yrs[pickCE],resid,"n",ylim=c(ymin,ymax),ylab="LogN Residuals",xlab="")
   grid()
   abline(h=1.0,col=1)
   segments(x0=yrs[pickCE],y0=rep(1.0,nresid),x1=yrs[pickCE],y1=resid,lwd=2,col=2)
   points(yrs[pickCE],resid,pch=16,col=1,cex=1.0)
   rmseresid <- sqrt(sum(resid^2)/nresid)
   text(min(yrs[pickCE]),ymin*1.05,paste("rmse = ",round(rmseresid,3),sep=""),
        font=7,cex=1.0,pos=4)
   # plot the depletion level
   ymax <- getmax(infish$deplete)
   plot(yrs,infish$deplete,type="l",lwd=2,ylim=c(0,ymax),yaxs="i",xlab="",
        panel.first=grid(),ylab="Depletion")
   abline(h=c(0.2,target),col=c(2,3),lwd=1)
   text(yrs[nyrs-1],0.9,round(infish$deplete[nrow(infish)],3),cex=1.0)
   if (nchar(png) > 0) dev.off()
} # end of plotASPM

#' @title plotceASPM plots just the fit of the ASPM model to the CPUE data
#' 
#' @description plotceASPM plots just the fit of the ASPM model to the CPUE data
#'     and provides a more detailed visual than plotASPM. It has the option of 
#'     including lognormal confidence intervals around the predcted CPUE.
#'
#' @param infish output from the dynamics function using the optimum parameters
#' @param CI the output matrix from getLNCI function. Needs to have the lower CI
#'    in column 1 and the upper CI in the third column.
#' @param defineplot defaults to TRUE, determines whether to set up an new 
#'     graphics window.
#'
#' @return returns nothing but does plot a graph
#' @export
#'
#' @examples
#' data(westroughy)
#' fish <- westroughy$fish
#' glb <- westroughy$glb
#' props <- westroughy$props
#' pars <- c(7.1,-1,-7.7)
#' bestaspm <- fitASPM(pars,minfun=dynamicsH,infish=fish,
#'                     inglb=glb,inprops=props)
#' out <- dynamicsH(bestaspm$estimate,fish,glb,props,full=TRUE)
#' ceCI <- getLNCI(out$fishery[,"predCE"],exp(bestaspm$estimate[2]))
#' plotceASPM(out$fishery,CI=ceCI)
plotceASPM <- function(infish,CI=NA,defineplot=TRUE) { 
# infish=out$fishery; CI=ceCI; defineplot=TRUE
   if (defineplot) { 
      if (names(dev.cur()) %in% c("null device", "RStudioGD"))
         dev.new(width = 7, height = 4.5, noRStudioGD = TRUE)
   } 
   par(mfrow=c(1,1),mai=c(0.4,0.5,0.1,0.05),oma=c(0.0,0,0.0,0.0)) 
   par(cex=1.0, mgp=c(1.35,0.35,0), font.axis=7,font=7,font.lab=7,tck=-0.02)  
   yrs <- infish$year
   if (inherits(CI,"matrix")) { 
      confint <- CI[,3]
    } else {   
      confint <- NA  
   }
   ymax <- getmax(c(infish$cpue,infish$predCE,confint))
   plot(yrs,infish$cpue,type="p",ylim=c(0,ymax),yaxs="i",
        xlab="",panel.first=grid(),ylab="Relative CPUE")
   if (inherits(CI,"matrix")) {
      arrows(x0=yrs,y0=CI[,1],x1=yrs,y1=CI[,3],code=3,length=0.03,angle=90,
             lwd=1,col=4)
   }   
   points(yrs,infish$cpue,pch=16,col=2,cex=1.0)     
   lines(yrs,infish$predCE,lwd=2,col=1)
} # end of plotceASPM

#' @title plotprops generates a 2x2 plot of the fishery properties
#' 
#' @description plotprops generates a 2 x 2 plot of the fishery properties of
#'     the length-at-age, maturity-at-age, weight-at-age, and selectivity-at-age
#'
#' @param rundir the directory in which the analysis is being run
#' @param props the matrix of fishery properties including laa, waa, maa, and 
#'     sela in columns 2 - 5 
#' @param console should the plot go to the console or be saved as a png file
#'     into rundir? default=TRUE ie plot to console
#'
#' @return nothing but it does generate a plot
#' @export
#'
#' @examples
#' data("westroughy")
#' plotprops(rundir="",westroughy$props,console=TRUE)
plotprops <- function(rundir,props,console=TRUE) {
  if (console) {
    filen <- "" 
  } else {
    filen <- pathtopath(rundir,"fishery_properties.png")
  }
  label <- c("Length-at-Age","Weight-at-Age","Maturity-at-Age",
             "Selectivity-at-Age")
  ages <- props[,1]
  if (console) plotprep(width=9, height=7,filename=filen)
  parset(plots=c(2,2),margin=c(0.25,0.5,0.1,0.1),outmargin=c(1,0,0,0))
  for (i in 2:5) {
    plot(ages,props[,i],type="l",lwd=3,xlab="",ylab=label[i-1],
         panel.first=grid())
  }
  mtext("Age Years",side=1,line=-0.1,outer=TRUE,cex=1.1)
} # end of plotprops

#' @title prodASPM summarizes ASPM statistics and plots the productivity
#'
#' @description prodASPM summarizes ASPM statistics and plots the productivity.
#'     The statistics it returns are the MSY, Bmsy, Hmsy (harvest rate that
#'     leads at equilibrium to MSY), Dmsy (the depletion level at which MSY
#'     occurs, and the predicted B0). These are printed to the console if the
#'     parameter 'console' is set to TRUE. In addition, a plot of the production
#'     curve, (yield vs Spawning Biomass), yield vs Harvest Rate, and yield vs
#'     Depletion level is produced if the parameter plot is set to TRUE
#'
#' @param inprod The production matrix generated by getProductionC
#' @param target in addition to the MSY what biomass target, in terms of
#'     target x B0 is wanted. Default = 0.48
#' @param console print the results directly to the console; default = TRUE
#' @param plot plot the yield vs Spawning biomass, harvest rate, and depletion.
#'     defaults to TRUE.
#'
#' @return a vector containing seven output statistics. It also prints to the
#'     console and plots the results if those parameters are set TRUE
#' @export
#'
#' @examples
#' \dontrun{
#' data(fishdat)
#' fish <- fishdat$fish
#' glb <- fishdat$glb
#' props <- fishdat$props
#' pars <- c(14,-1,-5)
#' bestL <- optim(pars,dynamicsH,method="Nelder-Mead",infish=fish,inglb=glb,
#'                inprops=props,control=list(maxit=1000,parscale=c(10,0.1)))
#' prod <- getProduction(exp(bestL$par[1]),infish=fish,inglb=glb,inprops=props,
#'                       Hrg=c(0.0005,0.07,0.0005),nyr=50)
#'  plotprep(width=7,height=6)
#' spsprod <- prodASPM(prod, console=TRUE, plot=TRUE)
#' }
prodASPM <- function(inprod, target=0.48, console=TRUE, plot=TRUE) {
   pickMSY <- which.max(inprod[,"Yield"])
   MSY <- inprod[pickMSY,"Yield"]
   Hmsy <- inprod[pickMSY,"Harvest"]
   Bmsy <-  inprod[pickMSY,"SpawnB"]
   B0 <- inprod[1,"SpawnB"]
   Dmsy <- inprod[pickMSY,"SpawnB"]/B0
   pickTarg <- which.closest(target,inprod[,"Depletion"])
   targC <- inprod[pickTarg,"Yield"]
   Htarg <- inprod[pickTarg,"Harvest"]
   Btarg <- inprod[pickTarg,"SpawnB"]
   if (console) {
      cat("MSY   = ", MSY,"\n")
      cat("Bmsy  = ", Bmsy,"\n")
      cat("Hmsy  = ", Hmsy,"\n")
      cat("Dmsy  = ", Dmsy,"\n")
      cat("B0    = ", B0, "\n")
      cat("targC = ", targC,"\n")
      cat("Htarg = ", Htarg,"\n")
      cat("Btarg = ", Btarg,"\n")
   }
   if (plot) {
      par(mfrow=c(3,1),mai=c(0.45,0.45,0.05,0.05),oma=c(0.1,0,0.0,0.0))
      par(cex=1.0, mgp=c(1.35,0.35,0), font.axis=7,font=7,font.lab=7)
      ymax <- getmax(inprod[,"Yield"],mult=1.1)
      plot(inprod[,"SpawnB"],inprod[,"Yield"],type="l",lwd=2,col=1,ylim=c(0,ymax),yaxs="i",
           xlab="Spawning Biomass (t)",ylab="Surplus Production (t)")
      grid()
      text(0.8*B0,0.9*MSY,paste0("MSY = ",round(MSY,3)),cex=1,pos=4)
      text(0.8*B0,0.75*MSY,paste0("targC = ",round(targC,3)),cex=1,pos=4)
      text(0.75*Bmsy,0.075*ymax,paste0("Bmsy = ",round(Bmsy,3)),cex=1,pos=4)
      abline(h=c(0,MSY),col=c(1,2))
      abline(v=c(Bmsy,Btarg,inprod[1,"SpawnB"]),col=c(2,3,"grey"))
      plot(rev(inprod[,"Harvest"]),rev(inprod[,"Yield"]),type="l",lwd=2,col=1,ylim=c(0,ymax),yaxs="i",
           xlim=c(max(inprod[,"Harvest"]),0),xlab="Annual Harvest Rate",ylab="Surplus Production (t)")
      grid()
      text(0.95*Hmsy,0.075*ymax,paste0("Hmsy = ",round(Hmsy,4)),cex=1,pos=4)
      abline(h=c(0,MSY),col=c(1,2))
      abline(v=c(Hmsy,Htarg),col=c(2,3))
      plot(inprod[,"Depletion"],inprod[,"Yield"],type="l",lwd=2,col=1,ylim=c(0,ymax),yaxs="i",
           xlab="Stock Depletion Level",ylab="Surplus Production (t)")
      grid()
      text(0.75*Dmsy,0.075*ymax,paste0("Dmsy = ",round(Dmsy,4)),cex=1,pos=4)
      abline(h=c(0,MSY),col=c(1,2))
      abline(v=c(Dmsy,target),col=c(2,3))
   }
   ans <- c(MSY=MSY,Bmsy=Bmsy,Hmsy=Hmsy,Dmsy=Dmsy,B0=B0,
            targC=targC,Htarg=Htarg,Btarg=Btarg)
   return(ans)
} # end of prodASPM

#' @title robustASPM conducts a robustness test on the quality of fit of an ASPM
#' 
#' @description robustASPM conducts a robustness test on the quality of fit of 
#'     an ASPM. This is done by using the original optimal model parameters or 
#'     the original guessed parameter values, add random variation to each of 
#'     them, and re-fit the model. This process needs to be repeated multiple 
#'     times. This should enable an analysis of the stability of the modelling 
#'     outcomes. If the optimum parameters are used then add more variation, if
#'     initial guesses are used you may need to select different starting points
#'     so that the random variation covers the parameter space reasonably well.
#'
#' @param inpar the parameter set to begin the trials with, this must be a 
#'    named vector of the optimum parameter values = LogR0, logsigCE, logavq
#' @param dynfun the dynamics function giving rise to -veLL, usually dynF, or
#'     dynamicsF, or dynamicsH
#' @param fish the fisheries data: at least year, catch, and cpue
#' @param glb the global variables containing the biological information
#' @param props the properties calculated from the globals
#' @param N the number of random trials to run; defaults to 10, which is too few
#' @param scaler the divisor that sets the degree of normal random variation to 
#'     add to the parameter values; default = 15 the smaller the value the more
#'     variable the outcome
#' @param Hrange a vector of three numbers denoting the range of harvest rates
#'     to use when characterizing the productivity implied by each fitted 
#'     parameter set. defaults to c(0.01, 0.45, 0.005)
#' @param numyrs the number of years used to drive each harvest rate to an
#'     equilibrium population structure
#' @param verbose print progress count to the screen? default = TRUE
#'
#' @return a list of results from each run, the range of values across runs, and
#'     the median values.
#' @export
#'
#' @examples
#' print("wait on data sets")
robustASPM <- function(inpar,dynfun,fish,glb,props,N=10,scaler=15,
                       Hrange=c(0.01,0.45,0.005),numyrs=50,verbose=TRUE) {
  origpar <- (inpar) # return Ln(R0) to nominal scale
  pars <- cbind(rnorm(N,mean=origpar[1],sd=abs(origpar[1])/scaler),
                rnorm(N,mean=origpar[2],sd=abs(origpar[2])/scaler),
                rnorm(N,mean=origpar[3],sd=abs(origpar[3])/scaler))
  columns <- c("iLnR0","isigmaCE","iavq","-veLL","LnR0","LsigCE","Lavq",
               "R0","sigCE","avq","MSY","B0","pardist","Iters") # prefix i implies input
  results <- matrix(0,nrow=N,ncol=length(columns),dimnames=list(1:N,columns))
  for (i in 1:N) { # i = 1   # 
    dist <- 0
    bestSP <- fitASPM(pars[i,],dynfun,fish,glb,props)
    opar <- bestSP$estimate
    for (pickp in 1:3) dist <- dist + (pars[i,pickp] - inpar[pickp])^2
    dist <- sqrt(dist)
    prod <- getProduction(exp(opar[1]),fish,glb,props,
                          Hrg=Hrange,nyr=numyrs)
    anspen <- prodASPM(prod,console=FALSE,plot=FALSE)
    results[i,] <- c(pars[i,],bestSP$minimum,bestSP$estimate,
                     exp(bestSP$estimate),anspen["MSY"],anspen["B0"],dist,i)
    if (verbose) cat(i," ")
  }
  ordres <- results[order(results[,"-veLL"]),] # see best and worst fit
  rownames(ordres) <- 1:N
  bounds <- apply(results,2,range)
  medvalues <- apply(results,2,median)
  return(list(results=ordres,range=bounds,medians=medvalues))
} # end of robustASMP

#' @title SpB - calculate spawning biomass from a vector of numbers-at-age
#'
#' @description SpB - calculates the spawning biomass from a vector of
#'     numbers-at-age, Maturity-at-age, and Weight-at-age.
#'
#' @param invect - the numbers-at-age as a vector
#' @param MatureA maturity at age vector
#' @param WeightA weight at age vector as kilograms
#' @return SpB - a scalar as tonnes.
#' @export
#' @examples
#' \dontrun{
#' data(fishdat)
#' str(fishdat)
#' glb <- fishdat$glb
#' fish <- fishdat$fish
#' props <- fishdat$props
#' unfish <- unfished(glb,props,glb$R0)
#' N0 <- unfish$N0
#' SpB(N0,props$maa,props$waa)  # should be 18326.52
#' unfish$B0
#' }
SpB <- function(invect, MatureA, WeightA) {
   ans <- sum(MatureA * WeightA * invect)/1000.0
   return(ans)
}

#' @title unfished generates the numbers at age for an unfished population
#'
#' @description unfished generates the numbers at age for an unfished
#'     population, and determines the recruitment dynamics. It requires the
#'     input of R0. The output includes the unfished numbers-at-age N0 plus 
#'     B0, A0, and R0.
#'     
#' @param glob the global constants object containing biology and structure
#' @param props the data.frame containing laa, maa, waa, maa, and sela
#' @param inR0 the log of unfished recruitment from B0
#' 
#' @return a list containing N0, R0, A0, and B0
#' @export
#' 
#' @examples
#' data("westroughy")
#' glb <- westroughy$glb  # contains a guess at log(R0)
#' props <- westroughy$props
#' unfish <- unfished(glob=glb,props=props,inR0=glb$R0)
#' print(unfish)
unfished <- function(glob,props,inR0) {
   R0 <- exp(inR0)
   maxage <- glob$maxage
   hsurv <- exp(-glob$M/2)
   surv <- exp(-glob$M)
   Nt <- numeric(maxage+1)
   Nt2 <- numeric(maxage+1)  # Needed for calculation of Exploitable biomass
   # now calculate numbers-at-age per recruit
   Nt[1] <- 1    # sets R0 = 1 for initiation NaA[0,0]
   for (age in 1:(maxage-1)) Nt[age+1] <- Nt[age] * surv
   Nt[maxage+1] <- (Nt[maxage] * surv)/(1-surv)
   # Estimate the biomass A0 generated by a recruitment of 1.0
   A0 <- SpB(Nt,props$maa,props$waa) # to get tonnes
   B0 <- R0 * A0
   # Now Generate the initial age distribution using R0
   Nages <- length(glob$ages)
   N0 <- R0*Nt
   res <- list(N0=N0, B0=B0, R0=R0, A0=A0)
   return(res)
}  # end of unfished




