library(tidyverse)

rm(list=ls())
gc()
# try(dev.off(),silent=TRUE)
# try(sink(),silent=TRUE)

Drive <- "D:"

#Base.dir <- file.path(Drive,"Stocks","hom_27_2a4a5b6a7a-ce-k8")
Base.dir <- file.path(Drive,"GIT")

#historic assessments here (not required)
#Assessment.Dir <- file.path(Base.dir,"Assessment")

#moved code to wk_WKREBUILD location 08/06/2020
#MSE.dir <- file.path(Base.dir,"MP_MSE","MSE 2019","whm.MSE2019.WKREBUILD")
#MSE.dir <- file.path(Base.dir,"MP_MSE","wk_WKREBUILD","EqSimWHM")
MSE.dir <- file.path(Base.dir,"wk_WKREBUILD","EqSimWHM")

#this is where the results of the 1000 assessment runs for initialisation of the MSE are saved
Data.dir <- file.path(MSE.dir,"Data")              
#any useful stuff contained in RData files
RData.dir <- file.path(MSE.dir,"RData")
#debug/verbose output
Log.dir <- file.path(MSE.dir,"Logs")              
#Simulation and statistical outputs
Res.dir <- file.path(MSE.dir, "Results")

Source.dir <- file.path(MSE.dir,"R")              #R functions
Scripts.dir <- file.path(MSE.dir, "Scripts")      #R scripts

#OMs, MPs
source(file = file.path(Scripts.dir,"OMs.R"))
source(file = file.path(Scripts.dir,"MPs.R"))

#source functions
source(file = file.path(Source.dir,"eqsim_Compare_runs.R"))

runs2Compare <- c("OM2.2_MP1.0","OM2.2_MP2.1")

#debug("fCompare_runs")
for (stat in c("Risk3")){
  fCompare_runs(runs2Compare = runs2Compare, 
                Res.dir = Res.dir, 
                Plot.dir = Res.dir,
                PerfStat = stat, 
                TargetFs = c(0, 0.05, 0.074, 0.1, 0.108, 0.2, 0.3),
                simYears = simYears, 
                lStatPer = lStatPer, 
                Blim = OM$refPts$Blim)}
#undebug("fCompare_runs")

load(file="SSB.RData")
ggplot(data = dfAll, 
       aes(x=factor(RunRef), y=Val, fill=factor(Period, levels = c("ST","MT","LT")))) + 
  geom_boxplot(outlier.size = 0.2) + 
  facet_wrap(vars(paste("Ftgt = ",Ftgt))) +
  # ggtitle(paste(StatName,StatUnit)) + 
  theme(legend.position="none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=12, face="bold"),
        axis.text.x = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="red", fill="#CCCCFF"))

load(file="RISK3.RData")
load(file="dfAll.RData")

dfAll %>% 
  mutate(t = paste("Ftgt =", Ftgt)) %>% 
  ggplot(aes(x=factor(RunRef), y=Val, fill=factor(Period, levels = c("ST","MT","LT")))) + 
  geom_col(position="dodge") + 
  geom_hline(yintercept = 0.05, col="black", linetype=2) +  
  facet_wrap(vars(t)) +
  theme(legend.position="none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=12, face="bold"),
        axis.text.x = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="red", fill="#CCCCFF"))

ggplot(data = dfAll,
       aes(x=factor(RunRef), y=Val, fill=factor(Period, levels = c("ST","MT","LT")))) +
  geom_col(position="dodge") + 
  geom_hline(yintercept = 0.05, col="black", linetype=2) +  
  facet_wrap(vars(paste("Ftgt = ",Ftgt))) +
  # ggtitle(paste(StatName,StatUnit)) + 
  theme(legend.position="none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=12, face="bold"),
        axis.text.x = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="red", fill="#CCCCFF"))
