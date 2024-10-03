

library(fmr)
library(codeutils)

dbdir <- getDBdir()
datadir <- pathtopath(dbdir,"A_Code/fmr/data-raw")

# Fournier & Archibald 1982---------------------------
source(pathtopath(datadir,"util_functions.R"))

# observed catch at age from Fournier and Archibald
ocaa <- read.csv(file=pathtopath(datadir,"ocaa.csv"),header=TRUE)
colnames(ocaa) <- c("year",seq(2,10,1))
ocaa

# observed weight-at-age from Fournier and Archibald
owaa <- read.csv(file=pathtopath(datadir,"owaa.csv"),header=TRUE)
colnames(owaa) <- c("year",seq(2,10,1))
owaa

# fishery data from Fournier and Archibald
fish <- read.csv(file=pathtopath(datadir,"fishery.csv"),header=TRUE)
fish[,"year"] <- 29:37
fish

# initial parameters from Fournier and Archibald 
param <- read.csv(file=pathtopath(datadir,"param.csv"),header=TRUE)
param


# western roughy 2018-------------------------------------
westroughy <- readdata(pathtopath(datadir,"westroughy.csv"),verbose=FALSE)



save(param,file=pathtopath(datadir,"param.RData"))
save(ocaa,file=pathtopath(datadir,"ocaa.RData"))
save(owaa,file=pathtopath(datadir,"owaa.RData"))
save(fish,file=pathtopath(datadir,"fish.RData"))
save(westroughy,file=pathtopath(datadir,"westroughy.RData"))

# check and transfer -------------------------------------------------------

tools::checkRdaFiles(paths=datadir)
tools::resaveRdaFiles(paths=datadir,compress="auto")
tools::checkRdaFiles(paths=datadir)



# Deriso et al 1989- Halibut-----------------------------------------------
# halibut from Deriso et al 1999
dbdir <- getDBdir()
wdir <- pathtopath(dbdir,"/A_Code/fmr/data-raw/")

caatab <- read.csv(file=pathtopath(wdir,"table8_16a.csv"),header=TRUE)

mfectab <- read.csv(file=pathtopath(wdir,"table8_16b.csv"),header=TRUE)

waatab <- read.csv(file=pathtopath(wdir,"table8_16c.csv"),header=TRUE)

years <- caatab[,"year"]
nyrs <- length(years)
ages <- as.numeric(gsub("a","",colnames(mfectab[,2:ncol(mfectab)])))
nages <- length(ages)
caacols <- ncol(caatab)
effort <- caatab[,caacols]; names(effort) = years
caa <- caatab[,2:(caacols-1)]; rownames(caa) <- years
mfeccols <- ncol(mfectab)
natM <- mfectab[1,2:mfeccols]
fec <- mfectab[2,2:mfeccols]
waacols <- ncol(waatab)
waa <- waatab[,2:waacols]; rownames(waa) <- years

halibut <- list(years=years,nyrs=nyrs,ages=ages,nages=nages,caa=caa,natM=natM,
                fec=fec,waa=waa,effort=effort)


datadir <- pathtopath(dbdir,"/A_Code/fmr/data/")

save(halibut,file=pathtopath(datadir,"halibut.RData"))


# make data for Francis 1992-----------------------------
library(codeutils)
library(hplot)
library(fmr)
ddir <- getDBdir()
rundir <- pathtopath(ddir,"/A_CodeR/mfmur/data-raw/")
source(pathtopath(rundir,"source_fmr.R"))


fish <- read.csv(file=pathtopath(rundir,"francis_1992.csv"))
colnames(fish) <- c("year","catch","index","cv")

data(fishdat)
glb <- fishdat$glb
glb$maxage <- 100
glb$M <- 0.05
glb$Linf <- 42.5
glb$K <- 0.059
glb$t0 <- -0.346
glb$Waa <- 0.0963
glb$Wab <- 2.68
glb$M50a <- 23
glb$deltaM <- 2
glb$sela50 <- 23
glb$deltaS <- 1
glb$steep <- 0.95
glb$R0 <- 20.0
glb$nages <- 101
glb$ages <- 0:100
glb$nyrs <- 12
glb$spsname <- "Orange_Roughy"


columns <- c("age","laa","waa","maa","sela")
props <- as.data.frame(matrix(0,nrow=glb$nages,ncol=length(columns),
                              dimnames=list(glb$ages,columns)))
props[,"age"] <- glb$ages
vbpars <- c(glb$Linf,glb$K,glb$t0)
props[,"laa"] <- vB(vbpars,glb$ages)

wtatage <- function(wtpar,laa) {
  wtaa <- wtpar[1] * (laa ^ wtpar[2])
  return(wtaa)
} # end of wtatage

wtpar <- c(glb$Waa,glb$Wab)
props[,"waa"] <- wtatage(wtpar=wtpar,laa=props[,"laa"]) 
props[,"maa"] <- logist(inL50=glb$M50a,delta=glb$deltaM,depend=glb$ages)
props[,"sela"] <- logist(inL50=glb$sela50,delta=glb$deltaS,depend=glb$ages)

francis92 <- list(fish=fish,glb=glb,props=props)

dbdir <- getDBdir()
datadir <- pathtopath(dbdir,"/A_Code/fmr/data/")

save(francis92,file=pathtopath(datadir,"francis92.RData"))

# check and transfer -------------------------------------------------------

tools::checkRdaFiles(paths=datadir)
tools::resaveRdaFiles(paths=datadir,compress="auto")
tools::checkRdaFiles(paths=datadir)





