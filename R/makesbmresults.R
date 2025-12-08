




#' @title plotbiology generates a pot of the maturity, emergence, WtL, and MatWt
#' 
#' @description plotbiology only generates a 2 x 2 plot of maturity-at-length,
#'     emergence-at-length, weight-at-length, and a combination of maturity 
#'     and weight-at-length (to simplify a few repeated calculations). If a 
#'     filen is given then the plotfile is automatically stored in rundir
#'     and prepared for inclusion in the results webpage.
#'
#' @param rundir the directory for all ctrl, data, and output files.
#' @param biol the matrix of biological properties (except growth).
#' @param console should the plot go to the console or be saved? Default=TRUE 
#' 
#' @seealso{
#'  \link{makebiol}, \link{makehtml}
#' }
#'
#' @return nothing but it does generate a plot to the console or to rundir
#' @export
#'
#' @examples
#' print("wait on data sets")
#' #   plotbiology(rundir=resdir,biol=biol,console=TRUE)
plotbiology <- function(rundir,biol,console=TRUE) {
  if (console) {   filen="" } else {
    filen <- pathtopath(rundir,"biological_properties.png")
  }
  mids <- as.numeric(rownames(biol))
  plotprep(width=8,height=7,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=c(2,2))
  plot1(mids,biol[,"maturity"],lwd=2,defpar = FALSE,ylab="Maturity")
  plot1(mids,biol[,"emergence"],lwd=2,defpar = FALSE,ylab="Emergence")
  plot1(mids,biol[,"WtL"],lwd=2,defpar = FALSE,ylab="Weight-at-Length")
  plot1(mids,biol[,"matwt"],lwd=3,defpar = FALSE,ylab="Maturity x Weight")
  lines(mids,biol[,"mature"]*max(biol[,"WtL"]),lwd=2,col=2)
  lines(mids,biol[,"WtL"],lwd=2,col=3)
  legend("topleft",c("MatWt","Maturity","Weight"),col=c(1,2,3),lwd=3,bty="n",cex=1.5)
  if (!console) {
    caption <- "Biological properties used in the assessment."
    addplot(filen,rundir=rundir,category="biology",caption)
  }
  if (!console) addtable(round(biol,6),"biological_properties.csv",rundir,
                         category="biology",caption = "Biological properties")
} # end of plotbiology

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
    filen <- pathtopath(rundir,"optimum_population_dynamics.png")
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
  plot1(years,dyn[,"deplexB"],lwd=2,ylab="Exploitable Biomass",defpar=FALSE)
  label <- paste0("Depl = ",round(dyn[nyrs,"deplexB"],3))
  mtext(text=label,side=3,line=-1.1,cex=1.0)
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

#' @title plotsizecomp generates a plot of available size composition data
#' 
#' @description plotsizecomp generates a plot of available size composition 
#'     data with the option of plotting the predicted size-composition of
#'     catch on top of it; this latter option is only possible when plotting 
#'     the distributions as proportions.
#'
#' @param rundir the directory for all ctrl, data, and output files.
#' @param incomp filters size-composition data from the fishery
#' @param lml a vector of the lml in each year represented by incomp
#' @param catchN the predicted size-composition of the catch
#' @param start which size-class to start from, default=NA, which means all 
#'     size-classes in the samples will be used.
#' @param proportion should the plots be as proportions or counts,default=TRUE
#' @param console should the plot be sent ot the console or saved to rundir.
#'     default=TRUE
#' @param width the width of the plot, default = 10
#' @param height the height of the plot, default = 9
#' @param fnt what font to use, default font = 7, bold times
#' @param outcategory which tab to put plot, default=NAS, but could be FIS
#' @param filename the final file name
#'
#' @return nothing but it does generate a plot
#' @export
#'
#' @examples
#' \dontrun{
#'   plotsizecomp(rundir,sizecomp,lml=LML,catchN=out$catchN,proportion=TRUE,
#'                console=FALSE)
#' }
plotsizecomp <- function(rundir,incomp,lml=NULL,catchN=NULL,start=NA,
                         proportion=TRUE,console=TRUE,width=10,height=9,fnt=7,
                         outcategory="NAS",filename="sizecomp_nas_by_year.png") { 
  # incomp=sizecomp; Nt=out$Nt; proportion=TRUE
  if (console) { filen <- "" } else {
    filen <- pathtopath(rundir,filename)
  }
  mids <- as.numeric(rownames(incomp))
  nmids <- length(mids)
  if (is.numeric(start)) {
    mids <- mids[start:nmids]
    nmids <- length(mids)
  }
  yrsize <- as.numeric(colnames(incomp))
  nyrs <- length(yrsize)
  nobs <- colSums(incomp,na.rm=TRUE)
  if (proportion) {
    incomp <- prop.table(incomp,margin=2)
    if (!is.null(catchN)) {
      years <- as.numeric(colnames(catchN))
      midpts <- as.numeric(rownames(catchN))
      picky <- match(yrsize,years)
      picks <- match(mids,midpts)
      pNt <- prop.table(catchN[picks,picky],margin=2)
    }
  }
  plotprep(width=width,height=height,newdev=FALSE,filename=filen,
           verbose=FALSE,usefont=fnt)
  parset(plots=pickbound(nyrs),margin=c(0.225,0.225,0.05,0.05),
         outmargin=c(1,1,0,0),byrow=FALSE)
  for (yr in 1:nyrs) {
    y <- incomp[,yr]
    ymax <- getmax(y)
    if ((proportion) & (!is.null(catchN))) ymax <- getmax(c(y,pNt[,yr]))
    plot(mids,y,type="l",lwd=2,ylab="",ylim=c(0,ymax),yaxs="i",
         panel.first = grid(),xlab="")
    mtext(text=yrsize[yr],side=4,line=-1.1,cex=1.0,font=7)
    text(mids[trunc(nmids*0.65)],0.9*ymax,nobs[yr],cex=1.25,pos=4)
    if ((proportion) & (!is.null(catchN))) {
      lines(mids,pNt[,yr],lwd=2,col=2)      
    }
    if (!is.null(lml)) abline(v=lml[yr],lwd=1,col="blue")
  }
  mtext("Shell Length (mm)",side=1,line=-0.2,outer=TRUE,cex=1.2)
  mtext("Year of Fishery",side=2,line=-0.2,outer=TRUE,cex=1.2)  
  if (!console) {
    caption <- "Catch size-composition data used with SBM"
    addplot(filen,rundir,category=outcategory,caption=caption)
  }
} # end of plotsizecomp
