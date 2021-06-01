

library(fmr)
library(rutilsMH)


datadir <- "C:/Users/User/Dropbox/A_Code/fmr/data-raw/"


ocaa <- read.csv(file=paste0(datadir,"ocaa.csv"),header=TRUE)

colnames(ocaa) <- c("year",seq(2,10,1))
ocaa

owaa <- read.csv(file=paste0(datadir,"owaa.csv"),header=TRUE)
colnames(owaa) <- c("year",seq(2,10,1))
owaa

fish <- read.csv(file=paste0(datadir,"fishery.csv"),header=TRUE)
fish[,"year"] <- 29:37
fish

param <- read.csv(file=paste0(datadir,"param.csv"),header=TRUE)
param


save(param,file=paste0(datadir,"param.RData"))
save(ocaa,file=paste0(datadir,"ocaa.RData"))
save(owaa,file=paste0(datadir,"owaa.RData"))
save(fish,file=paste0(datadir,"fish.RData"))

# check and transfer -------------------------------------------------------

tools::checkRdaFiles(paths=datadir)
tools::resaveRdaFiles(paths=datadir,compress="auto")
tools::checkRdaFiles(paths=datadir)





