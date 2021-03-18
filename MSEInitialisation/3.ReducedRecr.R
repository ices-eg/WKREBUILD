#poor recruitment scenarios

#this script superceded by code appended to 2.FLStocks.R

rm(list=ls())
gc()

library(FLCore)

Base <- "WGWIDE19"

mainDir <- getwd()

load(file=file.path(mainDir,Base,"RData","WHOM_SS19_FLS_V3.RData"))

#assessment terminal year
tyr <- range(FLSs.1k[,,,,,1])["maxyear"]

for(iter in seq(2,1000)){stock.n(FLSs.1k[1:5,as.character(tyr),,,,iter]) <- 0.5*as.numeric(stock.n(FLSs.1k[1:5,as.character(tyr),,,,iter]))}

save(FLSs.1k,file=file.path(getwd(),Base,"RData","WHOM_SS19_FLS_V3_RR.RData"))

#multipliers for ages 0-4 based on an assumption that the recruitment over the last 5 years is 
#equivalent to the mean of the 5 lowest recorded recruitments
#see excel sheet RecruitmentScenarios.xlsx

FLSs.1k.RR_5lowest <- FLSs.1k

plot(FLSs.1k.RR_5lowest)

for(iter in seq(1,1000)){
  cat(iter,"\n")
  #age 0
  stock.n(FLSs.1k.RR_5lowest[1,as.character(tyr),,,,iter]) <- 0.35*as.numeric(stock.n(FLSs.1k[1,as.character(tyr),,,,iter]))
  #age 1
  stock.n(FLSs.1k.RR_5lowest[2,as.character(tyr),,,,iter]) <- 0.20*as.numeric(stock.n(FLSs.1k[2,as.character(tyr),,,,iter]))
  #age 2
  stock.n(FLSs.1k.RR_5lowest[3,as.character(tyr),,,,iter]) <- 0.31*as.numeric(stock.n(FLSs.1k[3,as.character(tyr),,,,iter]))
  #age 3
  stock.n(FLSs.1k.RR_5lowest[4,as.character(tyr),,,,iter]) <- 0.35*as.numeric(stock.n(FLSs.1k[4,as.character(tyr),,,,iter]))
  #age 4
  stock.n(FLSs.1k.RR_5lowest[5,as.character(tyr),,,,iter]) <- 0.25*as.numeric(stock.n(FLSs.1k[5,as.character(tyr),,,,iter]))
  
}

plot(FLSs.1k.RR_5lowest)

save(FLSs.1k.RR_5lowest,file=file.path(getwd(),Base,"RData","WHOM_SS19_FLS_V3_RR_5lowest.RData"))
