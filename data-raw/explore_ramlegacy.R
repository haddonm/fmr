
library(codeutils)
library(ramlegacy)
load_ramlegacy()


str1(metadata)


datadir <- "C:/Users/Malcolm/Dropbox/A_Code/fmr/data-raw/"

load(paste0(datadir,"DBdata.RData"))

pick <- which(metadata$stockid == "NZSNAPNZ7")

metadata[pick,]




str1(cpue.data)


cpue <- cpue.data[,"NZSNAPNZ7"]


pick <- which(timeseries$stockid == "NZSNAPNZ7")
length(pick)

tsdat <- timeseries[pick,]

unique(tsdat$tsid)

pick1 <- which(tsdat$tsid == "TC-MT")
tsdat[pick1,]

tsdat[(tsdat$tsid == "TC-MT"),c("tsyear","tsvalue")]

tsdat[(tsdat$tsid == "CPUE-index"),c("tsyear","tsvalue")]


tsdat[(tsdat$tsid == "TL-MT"),]



