

library(codeutils)
library(hplot)

datadir <- "X:/csiro/A_CSIRO/Rcode/SESSF/ss3/fld2016/data/"

load(paste0(datadir,"FLDCAF.RDA"))

str1(spage)
properties(spage)

head(spage)

columns <- c(1,3,4,9,10,11,12,13,14,15,18,21,22,23,25)

head(spage[,columns],20)

table(spage$ZoneOK,spage$ValidL)

pick <- which((spage$ValidL == TRUE) & (spage$Gear == "OT") &
                (spage$Length > -99))
length(pick)
dim(spage)

age <- droplevels(spage[pick,])

properties(age)

pick2 <- which(age$Age2 > -99)

plotprep(width=9, height=5,newdev=FALSE)
parset()
plot1(age$Age1[pick2],age$Age2[pick2],type="p",xlab="Age1",ylab="Age2")

model <- lm(age$Age2[pick2] ~ age$Age1[pick2] - 1)
abline(model,lwd=2,col=2)

abyy <- tapply(age$Age1,list(age$Length,age$Year),length)
