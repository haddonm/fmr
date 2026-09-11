

#' @title comparevars plots a given variables changes between two scenarios 
#' 
#' @description comparevars can be used when comparing the dynamics of two
#'     alternative scenarios when fitting an integrated assessment model. It 
#'     requires the fishery dynamics of each scenario and plots whichever 
#'     variables are being considered (see example). 
#'
#' @param yrs the years of the fishery dynamics
#' @param var1 a variable from the the dynamics of the first scenario
#' @param var2 the same variable from teh second scenario
#' @param varname The name of the variable being compared
#' @param scenarios a character vector of the names of the two scenarios
#' @param console should the plot be printed to the console or a file? 
#'     default = TRUE
#' @param rundir the full path of the rundir where analyses are occurring.
#'     defauklt = ''
#' @param prepplot default = FALSe, implying the plot will be part of a 
#'     composite whose structure is defined elsewhere. If TRUE, then a single
#'     comparison will be produced.
#' @param legpos default = 'topright' the position of the legend, which is
#'     simply the contents of the scenarios argument.
#' @param cex default = 1, font size for single plots, ie prepplot=TRUE
#'
#' @returns if console = FALSe it returns a filename and a file is generated,
#'     if console = TRUE, nothing is returned but it does generate a plot
#' @export
#'
#' @examples
#' \dontrun{  # illustrate typical use
#'   fish49 <- opt49$fishery
#'   fish45 <- opt45$fishery
#'   yrs <- fish49[,"year"]
#'   plotprep(width=10,height=8)
#'   parset(plots=c(2,2),cex=1)
#'   comparevars(yrs=yrs,var1=fish49[,"deplete"],var2=fish45[,"deplete"],
#'               varname="Spawning Biomass Depletion",
#'               scenarios=c("fish49","fish45"),console=TRUE,prepplot=FALSE)
#'   comparevars(yrs=yrs,var1=fish49[,"twlPCE"],var2=fish45[,"twlPCE"],
#'               varname="Predicted CPUE",
#'               scenarios=c("fish49","fish45"),console=TRUE,prepplot=FALSE)
#'   comparevars(yrs=yrs,var1=fish49[,"spawnB"],var2=fish45[,"spawnB"],
#'               varname="Predicted Spawning Biomass",
#'               scenarios=c("fish49","fish45"),console=TRUE,prepplot=FALSE)
#'   comparevars(yrs=yrs,var1=fish49[,"recruit"],var2=fish45[,"recruit"],
#'               varname="Predicted Recruitment",scenarios=c("fish49","fish45"),
#'               console=TRUE,prepplot=FALSE)
#' }
comparevars <- function(yrs,var1,var2,varname,scenarios,console=TRUE,rundir="",
                        prepplot=FALSE,legpos="topright",cex=1.0) {
  if (console) {
    filen <- ""
  } else {
    filen <- pathtopath(rundir,paste0("compare_IA",varname,".png"))
  }
  if (prepplot) {
    plotprep(width=8,height=4.5,filename=filen)
    parset(cex=cex)
  }
  maxy <- getmax(c(var1,var2))
  plot(yrs,var1,type="l",lwd=2,col=1,ylab=varname,xlab="",ylim=c(0,maxy),
       panel.first=grid())
  lines(yrs,var2,lwd=2,col=2)
  legend(legpos,legend=scenarios,lwd=3,col=c(1,2),bty="n",cex=1.2)
  if (!console) {
    dev.off()
    return (filen)
  }
} # end of comparevars

#' @title initialdynamics plots dynamics from initial parameter guesses
#' 
#' @description initialdynamics generates a plot of the predicted stock 
#'     cpue given a set off initial parameters. This can be used to find,
#'     using trial and error, a set of scaling parameters (R0 and q, and 
#'     sometimes selectivity), that keep the predicted stock cpue off the
#'     zero line and that intersect with the observed CPUE. Only plots
#'     one fleet at a time.
#' 
#' @param x the input data, at least a matrix of 'year' and 'cpue'
#' @param year character name of the year variable, default = 'year'
#' @param cpue character name of the CPUE or index variable, default = 'cpue'
#' @param predCE character name of the predicted CPUE,  default = 'predCE'
#' @param width default = 9, the width of the plot
#' @param height default = 6, the height of the plot
#' @param result the BO, total likelihood and cpueLL as a column matrix. 
#'     default = NULL, which means nothing added to plot
#' @param legcex default = 1.25 the font size for the legend if used
#' @param legloc thge location of the legend if used. default = 'topright'
#' @param ... other potential inputs, plotting parameters, etc.
#' 
#' @returns nothing but it does produce a plot 
#' @export 
#' 
#' @examples 
#' # think of something 
initialdynamics <- function(x,year='year',cpue='cpue',predCE='predCE',width=9,
                            height=6,result=NULL,legcex=1.25,legloc="topright",...) { 
  
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  fishery <- replacezeros(x)
  yrs <- fishery[,year]
  plotprep(width=width,height=height,cex=1.0,verbose=FALSE) 
  parset(plots=c(1,1),margin=c(0.3,0.4,0.05,0.05),byrow=FALSE)
  maxy <- getmax(fishery[,c(cpue,predCE)])
  plot(yrs,fishery[,cpue],type="p",pch=16,cex=1.0,col=1,
       ylab="CPUE",ylim=c(0,maxy),yaxs="i",xlab="",
       panel.first=grid())
  lines(yrs,fishery[,predCE],lwd=2,col=2)
  if (!is.null(result)) {
    result <- round(result,2)
    label <- c(paste0("B0     = ",result[1,]),paste0("totalLL = ",result[2,]),
               paste0("cpueLL = ",result[3,]))
    legend(legloc,label,col=0,lwd=0,bty="n",cex=legcex)
  }
} # end of initialdynamics

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
#' @param target target depletion level. Defaults to 0.48
#' @param usef defines the font to use usef(ont),default = 7 bold times
#' @param rundir default '', otherwise give a full path if saving a file
#' @param console default=TRUE, which implies the plot will go to console. If
#'     set TRUE it will save pathtopath(rundir,'scenario.png')
#' @param scenario default = '', but if a file is to be saved then this should\
#'     be given a name for the file, if left as '', then it will be called 
#'     ASPM.png
#' @param maxy default = 0, which means use the data to set ymax. If maxy > 0
#'     the ymax will be set to that value. Used to allow details to be seen
#'
#' @return Nothing, but it does plot six graphs in a single plot.
#' @export
#'
#' @examples
#' require(fmr)
#' data(westroughy)
#' fish <- westroughy$fish
#' glb <- westroughy$glb
#' props <- westroughy$props
#' pars <- c(7.1,-1,-7.7)  # logR0 logceSE logavq
#' bestaspm <- fitASPM(pars,minfun=dynamicsH,infish=fish,
#'                    inglb=glb,inprops=props)
#' out <- dynamicsH(bestaspm$estimate,fish,glb,props,full=TRUE)
#' plotASPM(out$fishery, console=TRUE)
#' ceCI <- getLNCI(out$fishery[,"predCE"],exp(bestaspm$estimate[2]))
#' plotASPM(out$fishery,CI=ceCI, console=TRUE)
#' # infish=fishery; CI=ceCI;defineplot=TRUE; target=0.48; usef=7;rundir=""
#' # console=TRUE; scenario=""; maxy=10
plotASPM <- function(infish,CI=NA,target=0.48,usef=7,
                     rundir="",console=TRUE,scenario="",maxy=0) { 
  filen <- ""
  if (!console) {
    if (nchar(scenario) == 0) scenario <- "ASPM" 
    filen <- pathtopath(rundir,paste0(scenario,".png"))
  }
  plotprep(width=7, height=5.5,filename=filen,verbose=FALSE)
  parset(plots=c(3,2),margin=c(0.25,0.4,0.1,0.05),cex=0.85,font=usef)
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
  if (maxy > 0) ymax <- maxy
  plot(yrs,infish$cpue,type="p",pch=16,col=2,cex=1.0,ylim=c(0,ymax),yaxs="i",
       xlab="",panel.first=grid(),ylab="Relative Abundance Index")
  lines(yrs,infish$predCE,lwd=2,col=1)
  if ("matrix" %in% class(CI)) {
    segments(x0=yrs,y0=CI[,1],x1=yrs,y1=CI[,3],lwd=1,col=4)
  }
  # plot harvest rate
  ymax <- getmax(infish$fullH)
  if (!is.null(infish$fullF)) ymax <- getmax(c(infish$fullH,infish$fullF))   
  plot(yrs,infish$fullH,type="l",lwd=2,ylim=c(0,ymax),yaxs="i",xlab="",
       panel.first=grid(),ylab="Annual Rate")
  if (!is.null(infish$fullF)) {
    lines(yrs,infish$fullF,lwd=2,col=2)
    legend("topleft",c("Harvest","InstantF"),col=c(1,2),lwd=3,bty="n",cex=1.0)
  }   
  # plot the residuals
  pickCE <- which(infish[,"cpue"] > 0)
  resid <- log(infish[pickCE,"cpue"]/infish[pickCE,"predCE"])
  nresid <- length(resid)
  ymax <- getmax(resid);    ymin <- getmin(resid,mult=1.1)
  plot(yrs[pickCE],resid,"n",ylim=c(ymin,ymax),ylab="LogN Residuals",xlab="")
  grid()
  abline(h=0.0,col=1)
  segments(x0=yrs[pickCE],y0=rep(0.0,nresid),x1=yrs[pickCE],y1=resid,lwd=2,col=2)
  points(yrs[pickCE],resid,pch=16,col=1,cex=1.0)
  rmseresid <- sqrt(sum(resid^2)/nresid)
  text(min(yrs[pickCE]),ymin*0.85,paste("rmse = ",round(rmseresid,3),sep=""),
       font=7,cex=1.0,pos=4)
  # plot the depletion level
  ymax <- getmax(infish$deplete)
  plot(yrs,infish$deplete,type="l",lwd=2,ylim=c(0,ymax),yaxs="i",xlab="",
       panel.first=grid(),ylab="Depletion")
  abline(h=c(0.2,target),col=c(2,3),lwd=1)
  text(yrs[nyrs-1],0.9,round(infish$deplete[nrow(infish)],3),cex=1.0)
  if (!console) dev.off()
  return(invisible(filen))
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


#' @title plotcompfit plots observed vs predicted composition proportions
#' 
#' @description plotcompfit generates a plot of the predicted proportional 
#'     composition in the catches versus the observed proportional composition
#'     data to illustrate the fit of the model to the observed data. The xlabel 
#'     of all plots is left blank, but the range of each x-axis is printed at
#'     the bottom of all. 
#'
#' @param obscomp the observed composition data from the fishery as a matrix 
#'     with ages or sizes as rownames and years as column names
#' @param predcomp the predicted composition data from the model as a matrix 
#'     with ages or sizes as rownames and years as column names. Only the years that match the observed years will be plotted.
#' @param analysis a character name for the comparison being made, 
#'     default='composition_fit'. Whatever is used should have no spaces as it is 
#'     also used in teh filename, if used
#' @param ylabel a character label for the y-axis, this is added to the 
#'     contents of the analysis argument, default='opt49_ages'
#' @param console should the plot got to the console or file, default = TRUE
#' @param outdir the directory into which to send the file is console = FALSE,
#'     default =''
#' @param predcol the plotted colour of the predicted line, default = 'red'
#' @param cex the generic size of the font used in the plots
#' @param topcex the font size used for the top labels = year and sample size
#' @param maxcat the maximum age or size to be plotted. Default = 0, which 
#'     implies all catcegories are plotted. This can be useful when each plot 
#'     has a long tail of very few observations.
#'
#' @returns It generates a plot and, invisibly, a list of the filename and 
#'     caption in case it is saved as a file, as well as the final obscomp
#'     and predcomp, if differences are to be plotted. 
#' @export
#'
#' @examples
#' \dontrun{
#'   predcomp <- opt49$catchN[,2:46]
#'   obscomp <- const$agecomp$twl
#'   plotcompfit(obscomp,predcomp,analysis="Model Fit",ylabel="opt49",
#'               console=TRUE,outdir="",predcol="red",xlabel="",cex=1.0,
#'               topcex=0.9,maxcat=15)
#' }
plotcompfit <- function(obscomp,predcomp,analysis="Composition_Fit",
                        ylabel="opt49_ages",console=TRUE,outdir="",
                        predcol="red",cex=1.0,topcex=0.9,maxcat=0) {
  sampsize <- round(colSums(obscomp,na.rm=TRUE),1)  
  fstsamp <- which(sampsize > 0)
  picksamp <- min(fstsamp):max(fstsamp)
  if ((length(picksamp > 1)) & (ncol(obscomp) > 1)) {
    obscomp <- obscomp[,picksamp]
    predcomp <- predcomp[,picksamp]
  }
  sampsize <- round(colSums(obscomp,na.rm=TRUE),1)   
  Nsamp <- ncol(obscomp)   
  obscomp <- prop.table(obscomp,margin=2)
  predcomp <- prop.table(predcomp,margin=2)  
  compcl <- as.numeric(rownames(obscomp))  # expects size or age classes 
  if ((maxcat > 1) & (maxcat < max(compcl))) {
    obscomp <- obscomp[1:maxcat,]
    predcomp <- predcomp[1:maxcat,]
  }
  compcl <- as.numeric(rownames(obscomp))  # expects size or age classes    
  label <- as.numeric(colnames(obscomp))   # expects years    
  if (Nsamp > 45) {
    warning(cat(ylabel," Composition data limited to maximum 45 years \n"))
    obscomp <- obscomp[,1:45]
    predcomp <- predcomp[,1:45]
    Nsamp <- 45
    sampsize <- sampsize[1:45]
    label <- label[1:45]
  }
  addyrs <- paste0(label[1],"_",label[length(label)])
  filen <- ""
  if (!console) {
    filen <- paste0(outdir,"/comp_fit_for_",analysis,"_",ylabel,"_",addyrs,".png")
  }
  caption <- paste0("Comp Model deviations for ",ylabel,"-composition data")
  if (Nsamp <= 25) {
    nr <- 5
    nc <- ceiling(Nsamp/5)
    hgt <- 6
  }
  if (Nsamp <= 45) {
    nr <- ceiling(Nsamp/5)
    nc <- 5
    hgt <- 12
  }
  plotprep(width=10,height=hgt,newdev=!console,filename=filen,cex=cex,
           verbose=FALSE)
  parset(outmargin=c(1,2,1,1),margin=c(0.2,0.2,0,0))
  matfor <- matrix(c(1:(nr*nc)),nr,nc,byrow=TRUE)
  layout(matfor,heights=rep(1,(nr*nc)),TRUE)
  maxy <- getmax(c(obscomp,predcomp))
  if (sampsize[1] > 0) {
    plot(compcl,obscomp[,1],type="l",lwd=3,col=1,ylim=c(0,maxy),xaxt="n",
         ylab="")
    lines(compcl,predcomp[,1],lwd=3,col=2)
  } else {  
    plotnull(xvals=as.numeric(rownames(obscomp))) 
  }
  mtext(label[1],side=3,outer=FALSE,cex=topcex,line=-1)
  mtext(trunc(sampsize[1]),side=3,outer=FALSE,line=-1,cex=topcex,adj=1)
  if (Nsamp > 1) {
    for (i in 2:Nsamp) {
      if (sampsize[i] > 0) {
        plot(compcl,obscomp[,i],type="l",lwd=3,col=1,ylim=c(0,maxy),xaxt="n",
             ylab="")
        lines(compcl,predcomp[,i],lwd=3,col=2)  
      } else {  plotnull() }    
      mtext(label[i],side=1,outer=FALSE,line=-0.75,cex=topcex)
      mtext(label[i],side=3,outer=FALSE,cex=topcex,line=-1)
      mtext(trunc(sampsize[i]),side=3,outer=FALSE,line=-1,cex=topcex,adj=1)
    }
  } 
  # else {
  #     for (i in 2:20) {
  #       if (sampsize[i] > 0) {
  #         plot(compcl,obscomp[,i],type="l",lwd=3,col=1,ylim=c(0,maxy),xaxt="n",
  #              ylab="")
  #         lines(compcl,predcomp[,i],lwd=3,col=2)
  #       } else {  plotnull()  }    
  #       mtext(label[i],side=3,outer=FALSE,cex=topcex,line=-1)
  #       mtext(trunc(sampsize[i]),side=3,outer=FALSE,line=-1,cex=topcex,adj=1)
  #     }
  #     if (sampsize[21] > 0) {
  #       plot(compcl,obscomp[,21],type="l",lwd=3,col=1,ylim=c(0,maxy),xaxt="n",
  #            ylab="")
  #       lines(compcl,predcomp[,21],lwd=3,col=2)
  #     } else {  
  #       plotnull(xvals=as.numeric(rownames(obscomp))) 
  #     }
  #     mtext(label[21],side=3,outer=FALSE,cex=topcex,line=-1)
  #     mtext(trunc(sampsize[21]),side=3,outer=FALSE,line=-1,cex=topcex,adj=1)
  #     if (Nsamp > 21) {
  #       for (i in 22:min(40,Nsamp)) {
  #         if (sampsize[i] > 0) {
  #           plot(compcl,obscomp[,i],type="l",lwd=3,col=1,ylim=c(0,maxy),
  #                xaxt="n",ylab="")
  #           lines(compcl,predcomp[,i],lwd=3,col=2)
  #         } else {  plotnull()  }    
  #         mtext(label[i],side=3,outer=FALSE,cex=topcex,line=-1)
  #         mtext(trunc(sampsize[i]),side=3,outer=FALSE,line=-1,cex=topcex,adj=1)
  #       }
  #     }
  #     if (Nsamp > 40) {
  #       if (sampsize[41] > 0) {
  #         plot(compcl,obscomp[,41],type="l",lwd=3,col=1,ylim=c(0,maxy),
  #              xaxt="n",ylab="")
  #         lines(compcl,predcomp[,41],lwd=3,col=2)
  #       } else {  
  #         plotnull(xvals=as.numeric(rownames(obscomp))) 
  #       }
  #       mtext(label[21],side=3,outer=FALSE,cex=topcex,line=-1)
  #       mtext(trunc(sampsize[21]),side=3,outer=FALSE,line=-1,cex=topcex,adj=1)
  #       if (Nsamp > 41) {
  #         for (i in 42:Nsamp) {
  #           if (sampsize[i] > 0) {
  #             plot(compcl,obscomp[,i],type="l",lwd=3,col=1,ylim=c(0,maxy),
  #                  xaxt="n",ylab="")
  #             lines(compcl,predcomp[,i],lwd=3,col=2)
  #           } else {  plotnull()  }    
  #           mtext(label[i],side=3,outer=FALSE,cex=topcex,line=-1)
  #           mtext(trunc(sampsize[i]),side=3,outer=FALSE,line=-1,cex=topcex,
  #                 adj=1)
  #         }
  #       }
  #     }
  #   }  
  # }
  txtlabel <- paste0(ylabel,"  ",analysis)
  mtext(text=txtlabel,side=2,outer=TRUE,cex=1.1,line=0.2)
  xtxtlabel <- paste0("Categories ",min(compcl)," - ",max(compcl))
  xtxtlabel <- paste0(xtxtlabel,"  Black lines = Observed Proportions")
  mtext(text=xtxtlabel,side=1,outer=TRUE,cex=1.1,line=-0.25)
  if (!console) dev.off()
  return(invisible(list(filename=filen,caption=caption,obscomp=obscomp,
                        predcomp=predcomp)))
} # end of plotcompfit

#' @title plotcompfit plots observed vs predicted composition proportions
#' 
#' @description plotcompfit generates a plot of the predicted proportional 
#'     composition in the catches versus the observed proportional composition
#'     data to illustrate the fit of the model to the observed data. The xlabel 
#'     of all plots is left blank, but the range of each x-axis is printed at
#'     the bottom of all. 
#'
#' @param diffcomp the difference between the observed composition data and the
#'     predicted. Dervied from the plotcompfit function
#'     with ages or sizes as rownames and years as column names
#' @param analysis a character name for the comparison being made, 
#'     default='deviations'. Whatever is used should have no spaces as it is 
#'     also used in the filename, if used
#' @param ylabel a character label for the y-axis, this is added to the 
#'     contents of the analysis argument, default='opt49_ages'
#' @param console should the plot got to the console or file, default = TRUE
#' @param outdir the directory into which to send the file is console = FALSE,
#'     default =''
#' @param cex the generic size of the font used in the plots
#' @param topcex the font size used for the top labels = year and sample size
#' @param relmin should the absolute difference, topright text, be relative to
#'     the minimum relative difference or teh actual relative difference,
#'     default = FALSE, giving the real absolute difference
#'
#' @returns It generates a plot and, invisibly, a list of the filename and 
#'     caption in case it is saved as a file. 
#' @export
#'
#' @examples
#' \dontrun{
#'   diffcomp <- compfit$obscomp - compfit$predcomp  # from plotcompfit
#'   diffout <- plotdiffcomp(diffcomp,analysis="composition_deviations",
#'                           ylabel="opt49_age",console=TRUE,outdir="",
#'                           cex=1.0,topcex=0.9) 
#' }
#' #  diffcomp;analysis="composition_deviations";ylabel="opt49_age"
#' #  console=TRUE;outdir="";cex=1.0;topcex=0.9
plotdiffcomp <- function(diffcomp,analysis="deviations",ylabel="opt49_ages",
                         console=TRUE,outdir="",cex=1.0,topcex=0.9,
                         relmin=FALSE) {
  Nsamp <- ncol(diffcomp) 
  sampos <- apply(diffcomp,2,function(x){sum(abs(x))})
  if (Nsamp <= 25) {
    nr <- 5
    nc <- ceiling(Nsamp/5)
    hgt <- 6
  }
  if (Nsamp <= 45) {
    nr <- ceiling(Nsamp/5)
    nc <- 5
    hgt <- 12
  }
  compcl <- as.numeric(rownames(diffcomp))  # expects size or age classes 
  ncompcl <- length(compcl)
  if (Nsamp > 45) {
    warning(cat(ylabel," Composition data limited to maximum 45 years \n"))
    diffcomp <- diffcomp[,1:45]
  }
  label <- as.numeric(colnames(diffcomp))   # expects years    
  addyrs <- paste0(label[1],"_",label[length(label)])
  filen <- ""
  if (!console) {
    filen <- paste0(outdir,"/agecomp_fit_for_",analysis,"_",addyrs,".png")
  }
  caption <- paste0("Agecomp Model fit ",ylabel,"-composition data for ",
                    analysis)
  reldiff <- round(sampos*1000)
  if (relmin) reldiff <- reldiff - min(reldiff)
  plotprep(width=10,height=hgt,newdev=!console,filename=filen,cex=cex,
           verbose=FALSE)
  parset(outmargin=c(1,2,1,1),margin=c(0.2,0.2,0,0))
  matfor <- matrix(c(1:(nr*nc)),nr,nc,byrow=TRUE)
  layout(matfor,heights=rep(1,(nr*nc)),TRUE)
  maxy <- getmax(diffcomp)
  miny <- getmin(diffcomp)
  if (sampos[1] > 0) {
    plot(compcl,diffcomp[,1],type="p",pch=16,col=1,cex=1.0,ylim=c(miny,maxy),
         xaxt="n",ylab="")
    abline(h=0,lwd=1,col=1)
    for (j in 1:ncompcl) lines(c(compcl[j],compcl[j]),c(0,diffcomp[j,1]),
                               lwd=2,col=2)
    mtext(reldiff[1],side=3,outer=FALSE,line=-1,cex=topcex,adj=1)
  } else {  
    plotnull(xvals=as.numeric(rownames(diffcomp))) 
  }
  mtext(label[1],side=3,outer=FALSE,cex=topcex,line=-1)  
  if (Nsamp > 1) {
    for (i in 2:Nsamp) {
      if (sampos[i] > 0) {
        plot(compcl,diffcomp[,i],type="p",pch=16,col=1,cex=1.0,
             ylim=c(miny,maxy),xaxt="n",ylab="")
        abline(h=0,lwd=1,col=1)
        for (j in 1:ncompcl) lines(c(compcl[j],compcl[j]),c(0,diffcomp[j,i]),
                                   lwd=2,col=2)
        mtext(reldiff[i],side=3,outer=FALSE,line=-1,cex=topcex,adj=1)
      } else {  plotnull() }    
      mtext(label[i],side=3,outer=FALSE,cex=topcex,line=-1)
    }
  }
  txtlabel <- paste0(ylabel,"  ",analysis)
  mtext(text=txtlabel,side=2,outer=TRUE,cex=1.1,line=0.2)
  xtxtlabel <- paste0("Categories ",min(compcl)," - ",max(compcl))
  xtxtlabel <- paste0(xtxtlabel,"  Black Dots = Deviations")
  mtext(text=xtxtlabel,side=1,outer=TRUE,cex=1.1,line=-0.25)
  if (!console) dev.off()
  return(invisible(list(filename=filen,caption=caption)))
} # end of plotdiffcomp

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
#' @param obsdata default = TRUE. To plot CPUE residuals the observed CPUE is 
#'     required. When plotting simulated data there are no 'observed' so set
#'     obsdata to FALSE, which omits the residual plots. 
#' @param year character name of the year column, default ='year'
#' @param recruit character name of the recruitment column, the predrec column
#'     contains the predicted recruitment without the deviate value. default =
#'     'recruit' and 'predrec'. If no predrec values leave as default
#' @param depl character name of the spawning biomass depletion column default 
#'     = 'deplsB'
#' @param spawnB character name of the spawning biomass column default 
#'     = 'spawnB'
#' @param gears character names of the fishing gear names, not columns.
#'     default = 'Trawl' and 'Autoline'
#' @param catch character names of the catches columns. default = 'twl',"auln 
#' @param instF character name of the ycolumns containing the estimates of 
#'     instantaneous F for each gear through time. default='twlPF', 'aulnPF"
#' @param cecols character name of the column of cpue information, predicted
#'     values for each gear followed by the observed values for that gear (if
#'     present). default = 'twlPCE','twlCE','aulnPCE','aulnCE'
#'
#' @returns nothing but it does generate a plot
#' @export
#'
#' @examples
#' print("Wait on example data - again")
#' print("Can now plot simulated dynamics")
plotdynfish <- function(outfish,console=TRUE,addtitle="",prepplot=TRUE,
                        rundir="",width=8,height=7,nfleet=2,obsdata=TRUE,
                        year="year",recruit=c("recruit","predrec"),
                        depl="deplsB",spawnB="spawnB",
                        gears=c("Trawl","Autoline"),
                        catch=c("twl","auln"),instF=c("twlPF","aulnPF"),
                        cecols=c("twlPCE","twlCE","aulnPCE","aulnCE")) {

  # outfish=outIA$fishery;console=TRUE;addtitle="";prepplot=TRUE;rundir=""
  #              width=8;height=7;nfleet=1;obsdata=TRUE
  #              year="year";recruit=c("recruit","predrec")
  #              depl="deplete";gears=c("Trawl");catch=c("twl");instF=c("fullF") 
  #              cecols=c("twlPCE","twlCE"); spawnB="spawnB"
  
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  fishery <- replacezeros(outfish)
  yrs <- fishery[,year]
  if (console & prepplot) {
    plotprep(width=width,height=height,cex=1.0,filename="",verbose=FALSE)
  } else {
    filen <- pathtopath(rundir,paste0(addtitle,"fishery_dynamics.png"))
    plotprep(width=width,height=height,cex=1.0,newdev=TRUE,filename=filen,
             verbose=FALSE)
  }
  if ((nfleet == 1) || (!obsdata)) {
    parset(plots=c(4,2),margin=c(0.3,0.4,0.05,0.05),byrow=FALSE)
  } else {
    parset(plots=c(5,2),margin=c(0.3,0.4,0.05,0.05),byrow=FALSE)
  }
  if (obsdata) {
    maxy <- getmax(fishery[,cecols[1:2]])
  } else { 
    maxy <- getmax(fishery[,cecols[1]]) 
  }
  plot(yrs,fishery[,cecols[1]],type="l",lwd=2,col=1,
       ylab=paste0(gears[1]," CPUE"),ylim=c(0,maxy),yaxs="i",xlab="",
       panel.first=grid())
  if (obsdata) {
    points(yrs,fishery[,cecols[2]],pch=16,cex=1.1,col=2)
    lines(yrs,fishery[,cecols[2]],lwd=1,col=2,lty=3)
  }
  if ((nfleet > 1) & (obsdata)) {
    legend("topright",c("Predicted","Observed"),col=c(1:nfleet),lwd=3,bty="n",
           cex=1.2)
  }
  if (nfleet == 2) {  # could generalize to nfleet > 2 with loop here
    if (obsdata) { 
      usece <- fishery[,cecols[3:4]]
      maxy <- getmax(usece)
    } else { 
      usece <- fishery[,cecols[2]]
      maxy <- getmax(usece) 
    }
    plot(yrs,usece[,1],type="l",lwd=2,col=1,
         ylab=paste0(gears[2]," CPUE"),
         ylim=c(0,maxy),yaxs="i",xlab="",panel.first=grid())
    if (obsdata) {
      points(yrs,usece[,2],pch=16,cex=1.1,col=2)
      lines(yrs,usece[,2],lwd=1,col=2,lty=3)
    }
  }
  # spawning depletion----------------
  maxy <- getmax(fishery[,depl])
  plot(yrs,fishery[,depl],type="l",lwd=2,col=1,ylab="Spawning Depletion",
       ylim=c(0,maxy),yaxs="i",xlab="",panel.first=grid())
  abline(h=c(0.4,0.2),lwd=c(1,1),col=c(3,2))
  # catches-----------------
  if (nfleet == 2) {
    totC <- rowSums(fishery[,c(catch)],na.rm=TRUE)
    totC[which(totC == 0)] <- NA
    maxy <- getmax(totC)
  } else {
    maxy <- getmax(fishery[,catch[1]])
  }
  plot(yrs,fishery[,catch[1]],type="l",lwd=2,col=1,ylab="Catches (t)",
       ylim=c(0,maxy),yaxs="i",xlab="",panel.first=grid())
  if (nfleet == 2) {
    lines(yrs,fishery[,catch[2]],lwd=2,col=2)
    lines(yrs,totC,lwd=2,col=4)
    legend("topleft",c(gears,"Total"),col=c(1,2,4),lwd=3,bty="n",
           cex=1.1)
  }
  # instantaneous F-------------------
  maxy <- getmax(fishery[,c(instF)])
  plot(yrs,fishery[,instF[1]],type="l",lwd=2,col=1,ylab="Instantaneous F",
       ylim=c(0,maxy),yaxs="i",xlab="",panel.first=grid())
  if (nfleet == 2) {
    lines(yrs,fishery[,instF[2]],lwd=2,col=2)
    legend("topleft",c(gears),col=c(1:nfleet),lwd=3,bty="n",cex=1.1)
  }
  # cpue residuals--------------------
  if (obsdata) {
    twlresid <- fishery[,cecols[2]]/fishery[,cecols[1]]    
    maxy <- getmax(twlresid); miny <- getmin(twlresid)
    plot(yrs,twlresid,type="p",pch=16,col=1,cex=1,
         ylab=paste0(gears[1]," CPUE Residuals"),
         ylim=c(miny,maxy),yaxs="i",xlab="",panel.first=grid())
    lines(yrs,twlresid,lwd=1,col=2,lty=2)
    abline(h=1,lwd=1,col=1)
  }
  if ((nfleet == 2) & (obsdata)) {
    usece <- fishery[,cecols[3:4]]
    aulnresid <- usece[,2]/usece[,1]
    maxy <- getmax(aulnresid); miny <- getmin(aulnresid)
    plot(yrs,aulnresid,type="p",pch=16,col=1,cex=1,
         ylab=paste0(gears[2]," CPUE Residuals"),
         ylim=c(miny,maxy),yaxs="i",xlab="",panel.first=grid())
    lines(yrs,aulnresid,lwd=1,col=2,lty=2)
    abline(h=1,lwd=1,col=1)
  }
  # spawning biomass-----------------
  maxy <- getmax(fishery[,spawnB])
  plot(yrs,fishery[,spawnB],type="l",lwd=2,col=1,ylab="Spawning Biomass",
       ylim=c(0,maxy),yaxs="i",xlab="",panel.first=grid())
  # recruitment---------------------
  maxy <- getmax(fishery[,recruit[1]])
  plot(yrs,fishery[,recruit[1]],type="l",lwd=2,col=1,ylab="Recruitment",
       ylim=c(0,maxy),yaxs="i",xlab="",panel.first=grid())
  if (length(fishery[,recruit[2]]) > 0) {
    lines(yrs,fishery[,recruit[2]],lwd=2,col=2)
  }
  if (length(fishery[,recruit[2]]) > 0) {
    recdevs <- fishery[,recruit[1]]/fishery[,recruit[2]]
    maxy <- getmax(recdevs); miny <- getmin(recdevs)
    plot(yrs,recdevs,type="p",pch=16,cex=1,ylim=c(miny,maxy),xlab="",
         ylab="Recruitment Deviates",panel.first=grid())
    lines(yrs,recdevs,lwd=1,col="grey")
    abline(h=1.0,lwd=1,col=1)
    pick1 <- which(recdevs == 1.0)
    if (length(pick1 > 0)) points(yrs[pick1],recdevs[pick1],pch=16,cex=1,col=2)
  }
} # end of plotdynfish

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
#' @return thefilen used, which defaults '', but generates a plot
#' @export
#'
#' @examples
#' data("westroughy")
#' plotprops(rundir="",westroughy$props,console=TRUE)
#' # rundir=rundir; props=props; console=TRUE
plotprops <- function(rundir,props,console=TRUE) {
  filen <- "" 
  if (!console)  filen <- pathtopath(rundir,"fishery_properties.png")
  label <- c("Length-at-Age","Weight-at-Age","Maturity-at-Age",
             "Selectivity-at-Age")  
  numcol <- ncol(props)
  ages <- props[,"age"]
  plotprep(width=9, height=7,filename=filen,verbose=FALSE)
  parset(plots=c(2,2),margin=c(0.25,0.5,0.1,0.1),outmargin=c(1,0,0,0))
  for (i in 2:4) {
    maxy <- getmax(props[,i])
    plot(ages,props[,i],type="l",lwd=3,xlab="",ylab=label[i-1],ylim=c(0,maxy),
         yaxs="i",panel.first=grid())
  }
  mtext("Age Years",side=1,line=-0.1,outer=TRUE,cex=1.1)
  if (numcol == 5) {
    plot(ages,props[,numcol],type="l",lwd=3,xlab="",ylab=label[4],
         panel.first=grid())
  } else {
    fleets <- colnames(props)
    maxy <- getmax(props[,5])
    plot(ages,props[,5],type="l",lwd=3,xlab="",ylab=label[4],
         ylim=c(0,maxy),yaxs="i",panel.first=grid())
    for (flt in 6:numcol) lines(ages,props[,flt],lwd=3,col=(flt-4))
    legend("bottomright",fleets[5:numcol],col=c(1:(numcol-4)),lwd=3,
           cex=1.2,bty="n")
  }
  return(invisible(filen))
} # end of plotprops
