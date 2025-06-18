

library(fmr)
library(codeutils)

dbdir <- getDBdir()
datadir <- pathtopath(dbdir,"A_CodeR/fmr/notpublic")

# Beverton_Holt 1957---------------------------
source(pathtopath(datadir,"util_functions.R"))

# observed catch at age from Beverton_Holt 57 
ocaa <- read.csv(file=pathtopath(datadir,"ocaa.csv"),header=TRUE)
colnames(ocaa) <- c("year",seq(2,10,1))
ocaa

# observed weight-at-age from Beverton_Holt 57 
owaa <- read.csv(file=pathtopath(datadir,"owaa.csv"),header=TRUE)
colnames(owaa) <- c("year",seq(2,10,1))
owaa

# fishery data from Beverton_Holt 57 
fish <- read.csv(file=pathtopath(datadir,"fishery.csv"),header=TRUE)
fish[,"year"] <- 29:37
fish

# initial parameters from Beverton_Holt 57 
param <- read.csv(file=pathtopath(datadir,"param.csv"),header=TRUE)
param


# western roughy 2018-------------------------------------
westroughy <- readdata(pathtopath(datadir,"westroughy.csv"),verbose=FALSE)



save(param,file=pathtopath(datadir,"param.RData"))
save(ocaa,file=pathtopath(datadir,"ocaa.RData"))
save(owaa,file=pathtopath(datadir,"owaa.RData"))
save(fish,file=pathtopath(datadir,"fish.RData"))


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


# western roughy 2018  make data------------------------------------
library(fmr)
library(codeutils)

dbdir <- getDBdir()
datadir <- pathtopath(dbdir,"A_CodeR/fmr/data-raw")
source(pathtopath(datadir,"util_functions.R"))

westroughy <- readdata(pathtopath(datadir,"westroughy.csv"),verbose=FALSE)
glb <- westroughy$glb
fish <- westroughy$fish

glb$R0 <-  8.0

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

westroughy <- list(fish=fish,glb=glb,props=props)

save(westroughy,file=pathtopath(datadir,"westroughy.RData"))

# check and transfer -------------------------------------------------------

tools::checkRdaFiles(paths=datadir)
tools::resaveRdaFiles(paths=datadir,compress="auto")
tools::checkRdaFiles(paths=datadir)


# Francis 1992-----------------------------
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

# Fournier - Archibald 1982---------------------
library(codeutils)
ddir <- getDBdir()
datadir <- pathtopath(ddir,"/A_CodeR/fmr/data/")

dat <- c(400,1000,600,200,300,50,10,20,5,2,
         1000,236,556,304,91,125,19,4,7,2,
         700,595,135,297,151,43,55,8,2,3,
         300,416,339,72,148,70,19,24,3,1,
         67,172,212,143,25,42,17,4,5,1,
         97,39,89,92,51,7,11,4,1,1,
         118,54,17,29,21,9,1,1,0,0,
         241,70,31,9,14,10,4,0,1,0,
         1007,144,40,17,5,7,5,2,0,0,
         301,596,81,21,8,2,3,2,1,0,
         478,169,281,29,6,2,0,0,0,0,
         470,273,84,113,9,1,0,0,0,0,
         378,280,157,46,58,4,1,0,0,0,
         917,207,122,48,9,8,0,0,0,0,
         357,512,96,43,12,2,1,0,0,0,
         861,211,289,50,20,5,1,1,0,0,
         738,512,120,153,25,9,2,0,0,0,
         629,429,271,55,60,8,3,1,0,0,
         465,364,225,122,21,19,2,1,0,0,
         463,276,206,118,59,9,8,1,0,0)

fournarch82 <- matrix(dat,nrow=20,ncol=10,byrow=TRUE,dimnames=list(1:20,4:13))

save(fournarch82,file=pathtopath(datadir,"fournarch82.RData"))



tools::checkRdaFiles(paths=datadir)
tools::resaveRdaFiles(paths=datadir,compress="auto")
tools::checkRdaFiles(paths=datadir)


