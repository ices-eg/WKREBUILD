#setup

rm(list=ls())
gc()
try(dev.off(),silent=FALSE)

#library(devtools)
#install_github("ices-tools-prod/msy")
#install.packages("ggplotFL", repos="http://flr-project.org/R")

library(msy)
library(ggplot2)
library(ggpubr)
library(knitr)
library(icesAdvice)
library(tidyverse)
library(r4ss)
library(Cairo)
library(FLCore)
library(knitr)
library(plyr)

ac <- function(str){as.character(str)}
fmt <- function(RP,dgt=3,fmt="f"){ac(formatC(RP,digits=dgt,format=fmt,big.mark = ","))}
fPrettyOP <- function(dfAll,Ass,Scenario,Type){
  op <- ""
  dfSub <- dfAll[dfAll$Assessment==Ass & dfAll$Scenario==Scenario & dfAll$Type==Type,]
  op <- paste0(op,"Blim = ",fmt(dfSub[dfSub$RP=="Blim",]$Val  )," t","\n")
  op <- paste0(op,"Bpa = ",fmt(dfSub[dfSub$RP=="Bpa",]$Val,dgt=0,fmt="d")," t","\n")
  cat("Flim =",fmt(dfSub[dfSub$RP=="Flim",]$Val,dgt=3,fmt="f"),"\n")
  cat("Fpa =",fmt(dfSub[dfSub$RP=="Fpa",]$Val,dgt=3,fmt="f"),"\n")
  cat("Initial Fmsy =",fmt(dfSub[dfSub$RP=="FMSY",]$Val,dgt=3,fmt="f"),"\n")
  Fmsy <- dfSub[dfSub$RP=="FMSY",]$Val
  note <- ""
  if (dfSub[dfSub$RP=="Fpa",]$Val<dfSub[dfSub$RP=="FMSY",]$Val){
    Fmsy <- dfSub[dfSub$RP=="Fpa",]$Val
    note <- "(Fpa)"
  }
  cat("MSY Btrigger =",fmt(dfSub[dfSub$RP=="MSYBtrigger",]$Val,dgt=0,fmt="d"),"t","\n")
  cat("Fp05 =",fmt(dfSub[dfSub$RP=="FP05",]$Val,dgt=3,fmt="f"),"\n")
  if (dfSub[dfSub$RP=="FP05",]$Val < Fmsy){
    Fmsy <- dfSub[dfSub$RP=="FP05",]$Val
    note <- "(Fp05)"
  }
  cat("Final Fmsy =",fmt(Fmsy,dgt=3,fmt="f"),note,"\n")
  op
}

roundUp <- function(x,to=1000)
{
  to*(x%/%to + as.logical(x%%to))
}

fileFormat <- "png"
plot.h <- 5
plot.w <- 6

RData.dir <- file.path(getwd(),"RData")
graphics.dir <- file.path(getwd(),"Graphics")
logs.dir <- file.path(getwd(),"Logs")
source.dir <- file.path(getwd(),"Source")
ass.dir <- file.path(getwd(),"..","..","..","Assessment")
stf.dir <- file.path(getwd(),"STF")

#load historical assessments
#WKWIDE2017
#source(file=file.path(source.dir,"SS3toFLStock_function.R"))

#WK17 <- fSS3toFLStock(dirSS3output = file.path(ass.dir,"2017 WKWIDE","FINAL_RUN","NO_PELGAS_RecDev_FixedFec"),
#                           stockname = "WKWIDE17", fbar_amin = 1, fbar_amax = 10)

#save off to avoid having to create FLStock object in future analyses
#WHOM.WKWIDE2017 <- WK17
#save(WHOM.WKWIDE2017, file=file.path(RData.dir,"WKWIDE2017.RData"))

#WKWIDE2017
load(file=file.path(RData.dir,"WKWIDE2017.RData"))
WK17 <- WHOM.WKWIDE2017
name(WK17) <- "WKWIDE17"

#WGWIDE2017
load(file=file.path(RData.dir,"WGWIDE2017.RData"))
WG17 <- WHOM.WGWIDE2017
name(WG17) <- "WGWIDE17"

#WGWIDE2018
load(file=file.path(RData.dir,"WGWIDE2018.RData"))
WG18 <- WHOM.WGWIDE2018
name(WG18) <- "WGWIDE18"

#WGWIDE2019
load(file=file.path(RData.dir,"WGWIDE2019.RData"))
WG19 <- WHOM.WGWIDE2019
name(WG19) <- "WGWIDE19"

rm(WHOM.WKWIDE2017,WHOM.WGWIDE2017,WHOM.WGWIDE2018,WHOM.WGWIDE2019)

#create a data frame to store all assessment results
dfWHM <- data.frame(
  Assessment = rep("WK17",length=dim(rec(WK17))[2]),
  Year = seq(min(as.numeric(dimnames(stock.n(WK17))[[2]])),max(as.numeric(dimnames(stock.n(WK17))[[2]]))),
  Rec = as.numeric(rec(WK17)),
  Rec.rel = as.numeric(rec(WK17))/mean(as.numeric(rec(WK17))),
  SSB = as.numeric(ssb(WK17)),
  SSB.rel = as.numeric(ssb(WK17))/mean(as.numeric(ssb(WK17))),
  FBar = as.numeric(fbar(WK17)),
  FBar.rel = as.numeric(fbar(WK17))/mean(as.numeric(fbar(WK17))),
  stringsAsFactors = FALSE
)

dfWHM <- dplyr::bind_rows(
  dfWHM,
  data.frame(
    Assessment = rep("WG17",length=dim(rec(WG17))[2]),
    Year = seq(min(as.numeric(dimnames(stock.n(WG17))[[2]])),max(as.numeric(dimnames(stock.n(WG17))[[2]]))),
    Rec = as.numeric(rec(WG17)),
    Rec.rel = as.numeric(rec(WG17))/mean(as.numeric(rec(WG17))),
    SSB = as.numeric(ssb(WG17)),
    SSB.rel = as.numeric(ssb(WG17))/mean(as.numeric(ssb(WG17))),
    FBar = as.numeric(fbar(WG17)),
    FBar.rel = as.numeric(fbar(WG17))/mean(as.numeric(fbar(WG17))),
    stringsAsFactors = FALSE
  )
)


dfWHM <- dplyr::bind_rows(
  dfWHM,
  data.frame(
    Assessment = rep("WG18",length=dim(rec(WG18))[2]),
    Year = seq(min(as.numeric(dimnames(stock.n(WG18))[[2]])),max(as.numeric(dimnames(stock.n(WG18))[[2]]))),
    Rec = as.numeric(rec(WG18)),
    Rec.rel = as.numeric(rec(WG18))/mean(as.numeric(rec(WG18))),
    SSB = as.numeric(ssb(WG18)),
    SSB.rel = as.numeric(ssb(WG18))/mean(as.numeric(ssb(WG18))),
    FBar = as.numeric(fbar(WG18)),
    FBar.rel = as.numeric(fbar(WG18))/mean(as.numeric(fbar(WG18))),
    stringsAsFactors = FALSE
  )
)


dfWHM <- dplyr::bind_rows(
  dfWHM,
  data.frame(
    Assessment = rep("WG19",length=dim(rec(WG19))[2]),
    Year = seq(min(as.numeric(dimnames(stock.n(WG19))[[2]])),max(as.numeric(dimnames(stock.n(WG19))[[2]]))),
    Rec = as.numeric(rec(WG19)),
    Rec.rel = as.numeric(rec(WG19))/mean(as.numeric(rec(WG19))),
    SSB = as.numeric(ssb(WG19)),
    SSB.rel = as.numeric(ssb(WG19))/mean(as.numeric(ssb(WG19))),
    FBar = as.numeric(fbar(WG19)),
    FBar.rel = as.numeric(fbar(WG19))/mean(as.numeric(fbar(WG19))),
    stringsAsFactors = FALSE
  )
)


#some basic plots

#stock summary
Cairo(file = file.path(graphics.dir,paste0("WK17_Summary",".",fileFormat)),
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
p<-plot(WK17)
print(p)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WG17_Summary",".",fileFormat)),
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
p<-plot(WG17)
print(p)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WG18_Summary",".",fileFormat)),
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
p<-plot(WG18)
print(p)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WG19_Summary",".",fileFormat)),
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
p<-plot(WG19)
print(p)
dev.off()


#SSB Plots
Cairo(file = file.path(graphics.dir,paste0("SSB_Abs_Comparison",".",fileFormat)),
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)

#bottom, left, top, right
#par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mar=c(5.1, 4.1, 4.1, 4.1))

plot(dfWHM$Year[dfWHM$Assessment=="WK17"],
     dfWHM$SSB[dfWHM$Assessment=="WK17"]/1e6,
     type="l", xlab="Year", ylab="SSB (Mt)", xlim=c(1982,2018), ylim=c(0,5),
     col="blue", lty=2)

lines(dfWHM$Year[dfWHM$Assessment=="WG17"],
      dfWHM$SSB[dfWHM$Assessment=="WG17"]/1e6, lty=2, col="red")

lines(dfWHM$Year[dfWHM$Assessment=="WG18"],
      dfWHM$SSB[dfWHM$Assessment=="WG18"]/1e6, lty=2, col="green")

lines(dfWHM$Year[dfWHM$Assessment=="WG19"],
      dfWHM$SSB[dfWHM$Assessment=="WG19"]/1e6, lty=2, col="black")

#candidate biomass reference points (Bloss, B1982, B2001, B2003)
points(x=c(1982,2001,2003,dfWHM$Year[which(dfWHM$SSB==min(dfWHM$SSB[dfWHM$Assessment=="WK17"]))]),
       y=c(dfWHM$SSB[dfWHM$Assessment=="WK17" & dfWHM$Year==1982],
           dfWHM$SSB[dfWHM$Assessment=="WK17" & dfWHM$Year==2001],
           dfWHM$SSB[dfWHM$Assessment=="WK17" & dfWHM$Year==2003],
           min(dfWHM$SSB[dfWHM$Assessment=="WK17"]))/1e6,
       pch=20, col="blue",cex=1.5)

points(x=c(1982,2001,2003,dfWHM$Year[which(dfWHM$SSB==min(dfWHM$SSB[dfWHM$Assessment=="WG17"]))]),
       y=c(dfWHM$SSB[dfWHM$Assessment=="WG17" & dfWHM$Year==1982],
           dfWHM$SSB[dfWHM$Assessment=="WG17" & dfWHM$Year==2001],
           dfWHM$SSB[dfWHM$Assessment=="WG17" & dfWHM$Year==2003],
           min(dfWHM$SSB[dfWHM$Assessment=="WG17"]))/1e6,
       pch=20, col="red",cex=1.5)

points(x=c(1982,2001,2003,dfWHM$Year[which(dfWHM$SSB==min(dfWHM$SSB[dfWHM$Assessment=="WG18"]))]),
       y=c(dfWHM$SSB[dfWHM$Assessment=="WG18" & dfWHM$Year==1982],
           dfWHM$SSB[dfWHM$Assessment=="WG18" & dfWHM$Year==2001],
           dfWHM$SSB[dfWHM$Assessment=="WG18" & dfWHM$Year==2003],
           min(dfWHM$SSB[dfWHM$Assessment=="WG18"]))/1e6,
       pch=20, col="green",cex=1.5)

points(x=c(1982,2001,2003,dfWHM$Year[which(dfWHM$SSB==min(dfWHM$SSB[dfWHM$Assessment=="WG19"]))]),
       y=c(dfWHM$SSB[dfWHM$Assessment=="WG19" & dfWHM$Year==1982],
           dfWHM$SSB[dfWHM$Assessment=="WG19" & dfWHM$Year==2001],
           dfWHM$SSB[dfWHM$Assessment=="WG19" & dfWHM$Year==2003],
           min(dfWHM$SSB[dfWHM$Assessment=="WG19"]))/1e6,
       pch=20, col="black",cex=1.5)

legend("topright", legend=c("WKWIDE 2017", "WGWIDE 2017", "WGWIDE 2018", "WGWIDE 2019"), lty=2, col=c("blue","red","green","black"),bty="n")

par(new=TRUE)

#relative change in absolute SSB between WKWIDE2017 and WGWIDE2018
plot(seq(1982,2015),
     100*(dfWHM$SSB[dfWHM$Assessment=="WG17" & dfWHM$Year %in% seq(1982,2015)]/
            dfWHM$SSB[dfWHM$Assessment=="WK17"]-1), type="l", col="red",
     axes=FALSE, ylim=c(-100,100), xlab="", ylab="", xlim=c(1982,2017))
lines(seq(1982,2016),
      100*(dfWHM$SSB[dfWHM$Assessment=="WG18" & dfWHM$Year %in% seq(1982,2016)]/
             dfWHM$SSB[dfWHM$Assessment=="WG17"]-1), type="l", col="green")
lines(seq(1982,2017),
      100*(dfWHM$SSB[dfWHM$Assessment=="WG19" & dfWHM$Year %in% seq(1982,2017)]/
             dfWHM$SSB[dfWHM$Assessment=="WG18"]-1), type="l", col="black")
abline(h=0,lty=2)
axis(side=4)
mtext(text="Absolute SSB revision (%)",side=4,line=2)

dev.off()


#relative SSB
Cairo(file = file.path(graphics.dir,paste0("SSB_Rel_Comparison",".",fileFormat)),
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)

par(mar=c(5.1, 4.1, 4.1, 4.1))

plot(dfWHM$Year[dfWHM$Assessment=="WK17"],
     dfWHM$SSB.rel[dfWHM$Assessment=="WK17"],
     type="l", xlab="Year", ylab="Relative SSB", xlim=c(1982,2018), ylim=c(0,2.5),
     col="blue", lty=2)

abline(h=1)
lines(dfWHM$Year[dfWHM$Assessment=="WG17"],
      dfWHM$SSB.rel[dfWHM$Assessment=="WG17"], col="red", lty=2)

lines(dfWHM$Year[dfWHM$Assessment=="WG18"],
      dfWHM$SSB.rel[dfWHM$Assessment=="WG18"], col="green", lty=2)

lines(dfWHM$Year[dfWHM$Assessment=="WG19"],
      dfWHM$SSB.rel[dfWHM$Assessment=="WG19"], col="black", lty=2)

legend("topright", legend=c("WKWIDE 2017", "WGWIDE 2017","WGWIDE 2018","WGWIDE 2019"), lty=2, col=c("blue","red","green","black"),bty="n")

#important candidate biomass reference points (Bloss, B1982, B2001)
points(x=c(1982,2001,dfWHM$Year[which(dfWHM$SSB.rel==min(dfWHM$SSB.rel[dfWHM$Assessment=="WK17"]))]),
       y=c(dfWHM$SSB.rel[dfWHM$Assessment=="WK17" & dfWHM$Year==1982],
           dfWHM$SSB.rel[dfWHM$Assessment=="WK17" & dfWHM$Year==2001],
           min(dfWHM$SSB.rel[dfWHM$Assessment=="WK17"])),
       pch=20, col="blue",cex=1.5)

points(x=c(1982,2001,dfWHM$Year[which(dfWHM$SSB.rel==min(dfWHM$SSB.rel[dfWHM$Assessment=="WG17"]))]),
       y=c(dfWHM$SSB.rel[dfWHM$Assessment=="WG17" & dfWHM$Year==1982],
           dfWHM$SSB.rel[dfWHM$Assessment=="WG17" & dfWHM$Year==2001],
           min(dfWHM$SSB.rel[dfWHM$Assessment=="WG17"])),
       pch=20, col="red",cex=1.5)

points(x=c(1982,2001,dfWHM$Year[which(dfWHM$SSB.rel==min(dfWHM$SSB.rel[dfWHM$Assessment=="WG18"]))]),
       y=c(dfWHM$SSB.rel[dfWHM$Assessment=="WG18" & dfWHM$Year==1982],
           dfWHM$SSB.rel[dfWHM$Assessment=="WG18" & dfWHM$Year==2001],
           min(dfWHM$SSB.rel[dfWHM$Assessment=="WG18"])),
       pch=20, col="green",cex=1.5)

points(x=c(1982,2001,dfWHM$Year[which(dfWHM$SSB.rel==min(dfWHM$SSB.rel[dfWHM$Assessment=="WG19"]))]),
       y=c(dfWHM$SSB.rel[dfWHM$Assessment=="WG19" & dfWHM$Year==1982],
           dfWHM$SSB.rel[dfWHM$Assessment=="WG19" & dfWHM$Year==2001],
           min(dfWHM$SSB.rel[dfWHM$Assessment=="WG19"])),
       pch=20, col="black",cex=1.5)

par(new=TRUE)
#relative change in relative SSB between WKWIDE2017 and WGWIDE2018
plot(seq(1982,2015),
     100*(dfWHM$SSB.rel[dfWHM$Assessment=="WG17" & dfWHM$Year %in% seq(1982,2015)]/
            dfWHM$SSB.rel[dfWHM$Assessment=="WK17"]-1), type="l", col="red",
     axes=FALSE, ylim=c(-100,100), xlab="", ylab="", xlim=c(1982,2017))
lines(seq(1982,2016),
      100*(dfWHM$SSB.rel[dfWHM$Assessment=="WG18" & dfWHM$Year %in% seq(1982,2016)]/
             dfWHM$SSB.rel[dfWHM$Assessment=="WG17"]-1), col="green")
lines(seq(1982,2017),
      100*(dfWHM$SSB.rel[dfWHM$Assessment=="WG19" & dfWHM$Year %in% seq(1982,2017)]/
             dfWHM$SSB.rel[dfWHM$Assessment=="WG18"]-1), col="black")
abline(h=0,lty=2)
axis(side=4)
mtext(text="Relative SSB revision (%)",side=4,line=2)

dev.off()



#FBar Plots
Cairo(file = file.path(graphics.dir,paste0("FBar_Abs_Comparison",".",fileFormat)),
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)

#bottom, left, top, right
#par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mar=c(5.1, 4.1, 4.1, 4.1))

plot(dfWHM$Year[dfWHM$Assessment=="WK17"],
     dfWHM$FBar[dfWHM$Assessment=="WK17"],
     type="l", xlab="Year", ylab="FBar (1-10)", xlim=c(1982,2018), ylim=c(0,0.25),
     col="blue", lty=2)

lines(dfWHM$Year[dfWHM$Assessment=="WG17"],
      dfWHM$FBar[dfWHM$Assessment=="WG17"], lty=2, col="red")

lines(dfWHM$Year[dfWHM$Assessment=="WG18"],
      dfWHM$FBar[dfWHM$Assessment=="WG18"], lty=2, col="green")

lines(dfWHM$Year[dfWHM$Assessment=="WG19"],
      dfWHM$FBar[dfWHM$Assessment=="WG19"], lty=2, col="black")

#abline(h=0.1079,col="red")

legend("topright", legend=c("WKWIDE 2017", "WGWIDE 2017", "WGWIDE 2018", "WGWIDE 2019"), 
       lty=2, col=c("blue","red","green","black"),bty="n")

par(new=TRUE)
#relative change in absolute FBar between assessments
plot(seq(1982,2015),
     100*(dfWHM$FBar[dfWHM$Assessment=="WG17" & dfWHM$Year %in% seq(1982,2015)]/
            dfWHM$FBar[dfWHM$Assessment=="WK17"]-1), type="l", col="red",
     axes=FALSE, ylim=c(-100,100), xlab="", ylab="", xlim=c(1982,2017))

lines(seq(1982,2016),
      100*(dfWHM$FBar[dfWHM$Assessment=="WG18" & dfWHM$Year %in% seq(1982,2016)]/
             dfWHM$FBar[dfWHM$Assessment=="WG17"]-1), col="green")

lines(seq(1982,2017),
      100*(dfWHM$FBar[dfWHM$Assessment=="WG19" & dfWHM$Year %in% seq(1982,2017)]/
             dfWHM$FBar[dfWHM$Assessment=="WG18"]-1), col="black")

abline(h=0,lty=2)
axis(side=4)
mtext(text="Absolute FBar revision (%)",side=4,line=2)

dev.off()


#relative FBar
#FBar Plots
Cairo(file = file.path(graphics.dir,paste0("FBar_Rel_Comparison",".",fileFormat)),
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)

#bottom, left, top, right
#par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mar=c(5.1, 4.1, 4.1, 4.1))

plot(dfWHM$Year[dfWHM$Assessment=="WK17"],
     dfWHM$FBar.rel[dfWHM$Assessment=="WK17"],
     type="l", xlab="Year", ylab="FBar (1-10)", xlim=c(1982,2018), ylim=c(0,2),
     col="blue", lty=2)

lines(dfWHM$Year[dfWHM$Assessment=="WG17"],
      dfWHM$FBar.rel[dfWHM$Assessment=="WG17"], lty=2, col="red")

lines(dfWHM$Year[dfWHM$Assessment=="WG18"],
      dfWHM$FBar.rel[dfWHM$Assessment=="WG18"], lty=2, col="green")

lines(dfWHM$Year[dfWHM$Assessment=="WG19"],
      dfWHM$FBar.rel[dfWHM$Assessment=="WG19"], lty=2, col="black")

legend("topright", legend=c("WKWIDE 2017", "WGWIDE 2017", "WGWIDE 2018", "WGWIDE 2019"), 
       lty=2, col=c("blue","red","green"),bty="n")

par(new=TRUE)
#relative change in absolute FBar between assessments
plot(seq(1982,2015),
     100*(dfWHM$FBar.rel[dfWHM$Assessment=="WG17" & dfWHM$Year %in% seq(1982,2015)]/
            dfWHM$FBar.rel[dfWHM$Assessment=="WK17"]-1), type="l", col="red",
     axes=FALSE, ylim=c(-100,100), xlab="", ylab="", xlim=c(1982,2017))

lines(seq(1982,2016),
      100*(dfWHM$FBar.rel[dfWHM$Assessment=="WG18" & dfWHM$Year %in% seq(1982,2016)]/
             dfWHM$FBar.rel[dfWHM$Assessment=="WG17"]-1), col="green")

lines(seq(1982,2017),
      100*(dfWHM$FBar.rel[dfWHM$Assessment=="WG19" & dfWHM$Year %in% seq(1982,2017)]/
             dfWHM$FBar.rel[dfWHM$Assessment=="WG18"]-1), col="black")

abline(h=0,lty=2)
axis(side=4)
mtext(text="Relative FBar revision (%)",side=4,line=2)

dev.off()

#Recruitment Plots
Cairo(file = file.path(graphics.dir,paste0("Recr_Abs_Comparison",".",fileFormat)),
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)

#bottom, left, top, right
#par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mar=c(5.1, 4.1, 4.1, 4.1))

plot(dfWHM$Year[dfWHM$Assessment=="WK17"],
     dfWHM$Rec[dfWHM$Assessment=="WK17"],
     type="l", xlab="Year", ylab="Recruitment", xlim=c(1982,2018), ylim=c(0,20e6),
     col="blue", lty=2)

lines(dfWHM$Year[dfWHM$Assessment=="WG17"],
      dfWHM$Rec[dfWHM$Assessment=="WG17"], lty=2, col="red")

lines(dfWHM$Year[dfWHM$Assessment=="WG18"],
      dfWHM$Rec[dfWHM$Assessment=="WG18"], lty=2, col="green")

lines(dfWHM$Year[dfWHM$Assessment=="WG19"],
      dfWHM$Rec[dfWHM$Assessment=="WG19"], lty=2, col="black")

#abline(h=0.1079,col="red")

legend("topright", legend=c("WKWIDE 2017", "WGWIDE 2017", "WGWIDE 2018", "WGWIDE 2019"), 
       lty=2, col=c("blue","red","green", "black"),bty="n")

par(new=TRUE)

#relative change in absolute FBar between assessments
plot(seq(1982,2015),
     100*(dfWHM$Rec[dfWHM$Assessment=="WG17" & dfWHM$Year %in% seq(1982,2015)]/
            dfWHM$Rec[dfWHM$Assessment=="WK17"]-1), type="l", col="red",
     axes=FALSE, ylim=c(-100,100), xlab="", ylab="", xlim=c(1982,2017))

lines(seq(1982,2016),
      100*(dfWHM$Rec[dfWHM$Assessment=="WG18" & dfWHM$Year %in% seq(1982,2016)]/
             dfWHM$Rec[dfWHM$Assessment=="WG17"]-1), col="green")

lines(seq(1982,2017),
      100*(dfWHM$Rec[dfWHM$Assessment=="WG19" & dfWHM$Year %in% seq(1982,2017)]/
             dfWHM$Rec[dfWHM$Assessment=="WG18"]-1), col="black")

abline(h=0,lty=2)
axis(side=4)
mtext(text="Absolute Recr revision (%)",side=4,line=2)

dev.off()
