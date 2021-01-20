# ====================================================================================
# 1.SamtoFLStock.r
# ====================================================================================

# install.packages("FLCore", repos="http://flr-project.org/R")
library(FLCore)

# devtools::install_github('fishfollower/SAM/stockassessment', ref='components')
library(stockassessment)

# install.packages("FLSAM", repos="http://flr-project.org/R")
library(FLSAM)

source("../SAM2FLR/SAM2FLSAM.r")
source("../SAM2FLR/SAM2FLSTOCK.r")
source("../SAM2FLR/writing.R")

resdir <- "D:/iWGWIDE/2020/06. Data/hom.27.2a4a5b6a7a-ce-k8/output/SAM"

WG2018SAM   <- SAM2FLSTOCK(sao.name="WHOM_2018", temp="D:/temp")
WG2019SAM   <- SAM2FLSTOCK(sao.name="WHOM_2019", temp="D:/temp")
WG2020SAM   <- SAM2FLSTOCK(sao.name="WHOM_2019", temp="D:/temp")

save(WG2018SAM, file=file.path(resdir,"WGWIDE2018SAM.RData"))
save(WG2019SAM, file=file.path(resdir,"WGWIDE2019SAM.RData"))
save(WG2020SAM, file=file.path(resdir,"WGWIDE2020SAM.RData"))

