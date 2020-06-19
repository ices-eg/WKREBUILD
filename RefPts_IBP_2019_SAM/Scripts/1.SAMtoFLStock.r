# ====================================================================================
# 1.SamtoFLStock.r
# ====================================================================================

library(FLCore)
library(FLSAM)
library(stockassessment)

source("../SAM2FLR/SAM2FLSAM.r")
source("../SAM2FLR/SAM2FLSTOCK.r")

WG18SAM   <- SAM2FLSTOCK(sao.name="WHOM_2018", temp="D:/temp")
save(WG18SAM, file="RefPts_IBP_2019_SAM/RData/WGWIDE2018SAM.RData")

WG19SAM   <- SAM2FLSTOCK(sao.name="WHOM_2019", temp="D:/temp")
save(WG19SAM, file="RefPts_IBP_2019_SAM/RData/WGWIDE2019SAM.RData")
