#tester

rm(list=ls())
gc()

library(FLCore)

fGetFLStockSelection <- function(stk){
  sel <- matrix(FLCore::harvest(stk), ncol = dim(stk)[2])
  Fbar <- matrix(FLCore::fbar(stk), ncol = dim(stk)[2])
  sel <- sweep(sel, 2, Fbar, "/")
  sel <- sel/max(sel[,seq(dim(sel)[2]-10,dim(sel)[2])])  #last 10 years to avoid noise at start
  sel
}

load(file=file.path("C:","Stocks","hom_27_2a4a5b6a7a-ce-k8","MP_MSE","wk_WKREBUILD","EqSimWHM","RData","lMSE_WGWIDE19_FLStocks_1k15PG.RData"))
load(file=file.path(getwd(),"WGWIDE19","WHOM_SS19_FLS_orig.RData"))

for (i in 1:1000){if (any(!FLCore::iterMeans(ssb(lWHM[[i]])/ssb(FLSs.1k[,,,,,i]))==1)) {cat(i,"\n")}}
for (i in 1:1000){if (any(!FLCore::iterMeans(fbar(lWHM[[i]])/fbar(FLSs.1k[,,,,,i]))==1)) {cat(i,"\n")}}
for (i in 1:1000){if (any(!FLCore::iterMeans(rec(lWHM[[i]])/rec(FLSs.1k[,,,,,i]))==1)) {cat(i,"\n")}}

oStk2 <- lWHM[[2]]
nStk2 <- FLSs.1k[,,,,,2]

stock.wt(oStk2)/stock.wt(nStk2)
catch.wt(oStk2)/catch.wt(nStk2)

oSel <- fGetFLStockSelection(oStk2)
nSel <- fGetFLStockSelection(nStk2)
oSel/nSel
