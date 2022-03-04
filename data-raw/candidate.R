



data(ocaa)
data(fish)
data(owaa)
data(param)
pars <- param[1:28,2]
yrs <- as.numeric(fish[,"year"])
ages <- as.numeric(2:10)
onaa <- ocaa[,ages]
owa <- as.matrix(owaa[,ages])
nage <- length(ages)
getssq2(pars,M=0.2,yrs,ages,onaa,owa,cpue=fish[,"obsce"]) #should be 332.5389





