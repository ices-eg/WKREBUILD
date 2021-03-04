#Report plots

#code to generate plots for the evaluation report

# ================================================================================================================
# 03 EqSim summarize
# 
# Summarize results of EqSim simulations
#
# 06/07/2020 tested on 1000 iters of SAM assessment
# ================================================================================================================

rm(list=ls())
gc()

library(tidyverse)

#computer specific locations
#Drive    <- "D:"
#Base.dir <- file.path(Drive,"GIT")
Drive    <- "C:"
Base.dir <- file.path(Drive,"Stocks","hom_27_2a4a5b6a7a-ce-k8","MP_MSE")
MSE.dir <- file.path(Base.dir,"wk_WKREBUILD","EqSimWHM")
Data.dir <- file.path(MSE.dir,"Data")   
RData.dir <- file.path(MSE.dir,"RData")
Res.dir <- file.path(MSE.dir, "Results")
stats.dir <- file.path(Res.dir, "Stats")
Scripts.dir <- file.path(MSE.dir, "Scripts")
Source.dir <- file.path(MSE.dir,"R")

source(file.path(MSE.dir,"R","utilities.R"))
sapply(list.files(path=file.path(Source.dir), pattern=".R", full.names=TRUE), source)

#recruitment modelling
#historic data (assessments)
WG19 <- loadRData(file.path(RData.dir,"WGWIDE19.RData")) %>% FLCore::setPlusGroup(., 15)
WG20 <- loadRData(file.path(RData.dir,"WGWIDE20.RData")) %>% FLCore::setPlusGroup(., 15)

dfAssessments <- data.frame(WG = rep("WG19",dim(WG19)[2]),
                            Mod = rep("SS3",dim(WG19)[2]),
                            Yr = seq(range(WG19)["minyear"],range(WG19)["maxyear"]),
                            SSB = as.numeric(FLCore::ssb(WG19)),
                            FBar = as.numeric(fbar(WG19)),
                            Rec = as.numeric(rec(WG19)))

dfAssessments <- dplyr::bind_rows(dfAssessments,
                            data.frame(WG = rep("WG20",dim(WG20)[2]),
                            Mod = rep("SS3",dim(WG20)[2]),
                            Yr = seq(range(WG20)["minyear"],range(WG20)["maxyear"]),
                            SSB = as.numeric(FLCore::ssb(WG20)),
                            FBar = as.numeric(fbar(WG20)),
                            Rec = as.numeric(rec(WG20))))


source(file = file.path(Scripts.dir,"OMs.R"))
OM_WGWIDE19 <- OM2.2
OM_WGWIDE20 <- OM2.3

Blimloss <- OM_WGWIDE19$refPts$Blim
SegregBlim  <- function(ab, ssb) log(ifelse(ssb >= Blimloss, ab$a * Blimloss, ab$a * ssb))
set.seed(1)
SRR <- eqsr_fit(window(WG19,1995,2018), remove.years = c(2018), nsamp=1000, models = c("SegregBlim"))

x.mult <- 1.1
y.mult <- 1.4
minSSB <- min(dfAssessments$SSB, max(dfAssessments$SSB) * 0.0125)
maxSSB <- max(dfAssessments$SSB) * x.mult
maxrec <- max(dfAssessments$SSB) * y.mult

ssb_eval <- seq(minSSB, maxSSB, length.out = 100)
n_mods <- 2000

sample_rec <- function(i) {
  FUN <-  match.fun(modset$model[i])
  exp(FUN(modset[i,], ssb_eval) + stats::rnorm(length(ssb_eval), sd = modset $ cv[i]) )
}

modset <- SRR$sr.sto
ids <- sample(1:nrow(modset), n_mods, replace = TRUE)
rec_sim <- sapply(ids, sample_rec)

out <- data.frame(grp = rep(1:length(ssb_eval), n_mods),
                  mid.grp = rep(ssb_eval, n_mods),
                  ssb = jitter(rep(ssb_eval, n_mods), 2), # jitter for nices plotting
                  rec = c(rec_sim),
                  model = rep(modset[ids,"model"], each = length(ssb_eval)))

#model vs obs
ggRec <- ggplot(data = out, mapping = aes(x=rec)) +
  stat_ecdf(geom = "line", pad=FALSE) + xlim(0,1.5e7) +
  stat_ecdf(data = filter(dfAssessments,WG=="WG19" & Yr>=1995 & Yr<=2018), mapping = aes(x=Rec), geom="point", pad=FALSE)


#some line plots of recruitment

if (OM$desc=="WGWIDE20"){SRR <- eqsr_fit(window(FLS,1995,2019), remove.years = c(2019), nsamp=niters, models = c("SegregBlim"))}


simRuns <- loadRData(file=file.path(Res.dir))

ggSSB <- ggplot(data = dfAssessments, mapping = aes(x=Yr, y=SSB/1e6, group=WG)) +
  geom_line() + ylim(0,6) + ylab("SSB (Mt)")

ggRecWG19 <- ggplot(data = filter(dfAssessments,WG=="WG19" & Yr>=1995 & Yr<=2018), mapping = aes(x=Rec)) +
  stat_ecdf(geom="point", pad=FALSE) + xlim(0,1.5e7)

ggRecWG20 <- ggplot(data = filter(dfAssessments,WG=="WG20" & Yr>=1995 & Yr<=2019), mapping = aes(x=Rec)) +
  stat_ecdf(geom="point", pad=FALSE) + xlim(0,1.5e7)



#catch opportunity histograms

dfCatchOps <- data.frame(Year=rep(seq(2016,2020),each=100),
                          catch = rnorm(500))

ggCatchOps <- ggplot(data = dfCatchOps, mapping = aes(catch, fill = cut(catch,10))) +
  geom_histogram(binwidth=0.25, show.legend = FALSE) + facet_wrap(Year,ncol=1,strip.position="right") +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggCatchOps

#use the selected HCR (5.23 - update this to some 'baseline' run?)
runRef <- "WHOM_SS3_OM2.2_MP5.23_1000_23"
Sim.dir <- file.path(Res.dir,runRef)

load(file=file.path(Sim.dir,paste0(runRef,"_SimRuns.Rdata")))

#target F = 0.075
sim <- SimRuns[["0.075"]]

ages <- seq(0,15)
nages <- length(ages)
simYears <- seq(2018,length=23)
nsimYears = length(simYears)
iters <- seq(1,1000)
niters <- 1000

cw <- array(sim$catW,
            dim = c(nages,nsimYears,niters),
            dimnames = list(age=ages, simYear=simYears, iter=iters))

sw <- array(sim$stkW,
            dim = c(nages,nsimYears,niters),
            dimnames = list(age=ages, simYear=simYears, iter=iters))

N <- array(sim$N,
            dim = c(nages,nsimYears,niters),
            dimnames = list(age=ages, simYear=simYears, iter=iters))

rec <- N[c("0"),,]

#get the historic data
dfSW <- readr::read_delim(file = file.path(Data.dir,"StockWeights.dat"),delim=",")
dfCW <- readr::read_delim(file = file.path(Data.dir,"CatchWeights.dat"),delim=",")
#single df
dfWeights <- dplyr::left_join(
  dfSW %>% pivot_longer(cols = paste0("Age",seq(0,15)), names_to = "Age", values_to = "SW", names_prefix = "Age"),
  dfCW %>% pivot_longer(cols = paste0("Age",seq(0,15)), names_to = "Age", values_to = "CW", names_prefix = "Age"),
  by=c("Year","Age"))

#dfWeights <- within(dfWeights, Age <- factor(Age, levels = ac(ages)))

histYears <- seq(range(dfWeights$Year)[1],dfWeights$Year[2])
allYears <- seq(min(histYears),max(simYears))

#Stock objects - point estimates and the 1000 iterations
FLStockfile    <- "WGWIDE19.RData"
FLS <- loadRData(file.path(RData.dir,FLStockfile)) %>% FLCore::setPlusGroup(., 15)
dfSAWeights = data.frame(Age = ages, iter=0, CW = rowMeans(catch.wt(FLS)), SW = rowMeans(stock.wt(FLS)), stringsAsFactors = FALSE)
dfSAWeights <- within(dfSAWeights, Age <- factor(Age, levels = ac(ages)))

FLStockSimfile <- "MSE_WGWIDE19_FLStocks_1k15PG.RData"
FLSs <- loadRData(file.path(RData.dir,FLStockSimfile))

#extract catch, stock weight info from each iteration
FLSs.CW <- lapply(lapply(FLSs,FLCore::catch.wt),rowMeans)
FLSs.SW <- lapply(lapply(FLSs,FLCore::stock.wt),rowMeans)
iCW <- matrix(unlist(FLSs.CW),nrow=nages,ncol=length(FLSs),dimnames=list(age=ac(ages),iter=seq(1,length(FLSs))))
iSW <- matrix(unlist(FLSs.SW),nrow=nages,ncol=length(FLSs),dimnames=list(age=ac(ages),iter=seq(1,length(FLSs))))
dfSAWeights <- dplyr::bind_rows(dfSAWeights,data.frame(Age = ac(ages), iter=rep(seq(1,length(FLSs)),each=nages), CW = c(iCW), SW = c(iSW)))
#make Age a factor so plot order appropriate
dfSAWeights <- within(dfSAWeights, Age <- factor(Age, levels = ac(ages)))


#expand df to store simulated values
dfWeights$Upper <- NA
dfWeights$Med <- NA
dfWeights$Lower <- NA

for (a in ages){
  lower <- apply(sw[ac(a),,],MARGIN=1,FUN=stats::quantile,probs=c(0.05))
  med <- apply(sw[ac(a),,],MARGIN=1,FUN=stats::quantile,probs=c(0.5))
  upper <- apply(sw[ac(a),,],MARGIN=1,FUN=stats::quantile,probs=c(0.95))
  
  dfWeights <- dplyr::bind_rows(dfWeights,
                                 data.frame(Year=simYears,Age=ac(rep(a,length(simYears))),Lower=lower,Med=med,Upper=upper))
}

dfWeights <- within(dfWeights, Age <- factor(Age, levels = ac(ages)))

#swap geom_segment for geom_line to clip horizontal extents
#eg geom_segment(aes(x=2,xend=4,y=20,yend=20))

gSW <- ggplot(data = dfWeights, mapping = aes(x=Year,y=SW)) +
  geom_line(aes(group=Age)) +
  geom_hline(data = filter(dfSAWeights,iter==0), mapping=aes(yintercept=SW), col="red") +
  geom_hline(data = dfSAWeights %>% group_by(Age) %>% summarise(mn = min(SW), mx = max(SW)), mapping=aes(yintercept=mx), col="red", lty=2) +
  geom_hline(data = dfSAWeights %>% group_by(Age) %>% summarise(mn = min(SW), mx = max(SW)), mapping=aes(yintercept=mn), col="red", lty=2) +
  facet_wrap(~Age) + ylab("Stock Weight")


simYears

#age plots
a<-7
lower <- apply(sw[ac(a),,],MARGIN=1,FUN=stats::quantile,probs=c(0.05))
med <- apply(sw[ac(a),,],MARGIN=1,FUN=stats::quantile,probs=c(0.5))
upper <- apply(sw[ac(a),,],MARGIN=1,FUN=stats::quantile,probs=c(0.95))

plot(simYears,lower,type="l",ylim=c(0,0.5),xlim=c(min(allYears),max(allYears)))
lines(simYears,upper,type="l")
lines(simYears,med)
for (ii in seq(99,length=10,by=100)){lines(simYears,sw[ac(a),,ii],col="grey")}
lines(dfWeights$Year[dfWeights$Age==a],dfWeights$SW[dfWeights$Age==a],col="red")

sw %>% 
#percentiles of sw
quantile(sw[6,'2018',],probs=c(0.05,0.5,0.95))
