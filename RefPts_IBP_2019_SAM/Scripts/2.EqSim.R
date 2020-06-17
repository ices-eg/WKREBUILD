# ==================================================================================
# 2. EqSim.R
#
# Run the EqSim scenario 5.3 that was selected in IBPHOM 2019
# ==================================================================================

library(msy)
library(ggplot2)
#library(ggpubr)
library(knitr)
#library(icesAdvice)
library(tidyverse)
#library(r4ss)
library(Cairo)
library(FLCore)
library(knitr)
library(plyr)

# each scenario was run without the 1982 data point
# source(file=file.path(getwd(),"Scripts","0.Setup.R"))

RPs <- c("Blim","Bpa","Flim","Fpa","MSYBtrigger","FMSY","FP05","FMSY_final")
nRPs <- length(RPs)

SRR.models <- c("SegregBlim")
rhologRec <- FALSE
recruitment.trim = c(3, -3)
Fcv <- 0.212; 
Fphi <- 0.423
nsim <- 1000

SRR.desc <- paste0(SRR.models,collapse="_")

dfResults <- data.frame("Assessment" = character(),
                        "Scenario" = character(),
                        "SRR" = character(),
                        "BioYrs" = character(),
                        "SelYrs" = character(),
                        "Fcv" = numeric(),
                        "Fphi" = numeric(),
                        "RP" = character(),
                        "Type" = character(),
                        "Val" = numeric(),
                        "RunTime" = character(),
                        stringsAsFactors = FALSE)

rt <- format(Sys.time(), "%d-%b-%Y %H.%M")


#fs to scan
#Fscan <- seq(0.01,0.2,len=20) #low res
Fscan <- seq(0.01,0.2,len=40) #high res

#Base assessments
ass <- c("WG18")
ass <- c("WG18","WG19","WG19SAM")

a <- "WG18"

for (a in ass) {
  
  cat("EqSim analysis based on output from",a,"assessment","\n**************************************************\n")
  
  stk <- get(a)
  stk.name <- a
  
  #terminal year
  ymax <- as.integer(max(dimnames(stk)$year))
  ymin <- as.integer(min(dimnames(stk)$year))
  
  #trim terminal year to ignore
  stk <- window(stk,ymin,ymax-1)
  ymax <- ymax - 1
  
  #biological years (last 10)
  bio.years <- c(ymax-9,ymax)

  #selection years (last 10)
  sel.years <- c(ymax-9,ymax)

  cat("Biology years",bio.years,"\n")
  cat("Selection years",sel.years,"\n")
  cat("Mean SSB",mean(ssb(stk)),"\n")
  cat("Mean FBar",mean(fbar(stk)),"\n")

  ###############################################CASE 5###########################################
  #Bpa = SSB for 2003
  #including 1982
  
  Bpa <- as.integer(ssb(stk)[,"2003"])
  Blim <- as.integer(Bpa/1.4)
  SegregBlim  <- function(ab, ssb) log(ifelse(ssb >= Blim, ab$a * Blim, ab$a * ssb))
  
  # NOTE: THIS HAS NOT ACTUALLY BEEN APPLIED (MP 17/6/2020)
  # check that the Blim is not lower than Bloss, if so, set the breakpoint to Bloss
  # this should only be the case with the WG18 assessment
  # Bloss5<-as.numeric(min(ssb(window(stk,1983,ymax))))
  # if (Blim<Bloss5) {
  #   cat("modifying SegregBlim as Bloss>Blim for ",a,"\n")
  #   cat("Blim = ", Blim, "Bloss5 = ", Bloss5,"\n")
  #   SegregBlim  <- function(ab, ssb) log(ifelse(ssb >= Bloss5, ab$a * Bloss5, ab$a * ssb))
  # }
  
  ###############################################CASE 5.3###########################################
  #Bpa = SSB for 2003, 1995 on
  
  cat("Fitting SRR, 1995 on...\n")
  
  SRR <- msy::eqsr_fit(window(stk,1995,ymax),
                     remove.years = c(),
                     nsamp=nsim, models = SRR.models)
  
  # cat("Breakpoint,Mean=",Blim,SRR$sr.det$a*Blim,"\n")
  
  #Flim
  cat("Flim calculation...\n")
  SIM.Flim <- eqsim_run(SRR,
                           bio.years = bio.years,
                           bio.const = FALSE,
                           sel.years = sel.years,
                           sel.const = FALSE,
                           Fscan = Fscan,
                           Fcv = 0, Fphi = 0,
                           rhologRec = rhologRec,
                           Blim = Blim, Bpa = Bpa,
                           recruitment.trim = recruitment.trim,
                           Btrigger = 0,
                           verbose = TRUE)
  
  Flim <- SIM.Flim$Refs2["catF","F50"]
  Fpa <- Flim/1.4
  
  #step 1 - no Btrigger but including stochasticity in 
  #population and fishery and assessment/advice error
  #Fmsy - 1st run, include assessment error
  #for initial FMSY candidate
  cat("Initial Fmsy calculation...\n")
  SIM1.Fmsy <- eqsim_run(SRR,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim, Bpa = Bpa,
                            recruitment.trim = recruitment.trim,
                            Btrigger = 0,
                            verbose = TRUE)
  
  Fmsy <- SIM1.Fmsy$Refs2["lanF","medianMSY"]
  
  if (Fmsy>Fpa) {Fmsy_final <- Fpa} else {Fmsy_final <- Fmsy}
  
  #Step 2 select MSYBtrigger
  #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
  MSYBtrigger <- Bpa
  
  #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
  cat("Precautionary check FP05...\n")
  SIM2.Fmsy <- eqsim_run(SRR,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim, Bpa = Bpa,
                            recruitment.trim = recruitment.trim,
                            Btrigger = MSYBtrigger5,
                            verbose = TRUE)
  
  # eqsim_plot(SIM2.Fmsy,catch=TRUE)
  # eqsim_plot_range(SIM2.Fmsy, type="median")
  
  Fp05 <- SIM2.Fmsy$Refs2["catF","F05"]
  
  #Final Fmsy
  if (Fmsy_final>Fp05) {Fmsy_final <- Fp05}
  

  cat("**********************************************\n")
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("5.3",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Abs",nRPs),
                                           "Val" = c(Blim,Bpa,Flim,Fpa,MSYBtrigger,Fmsy,Fp05,Fmsy_final),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))


} #loop over assessments

save(dfResults, file = file.path(RData.dir, "dfResults.RData"))

# overview

dfResults %>% 
  filter(RP %in% c("Blim","MSYBtrigger","Flim","Fpa", "FMSY","FMSY_final","FP05")) %>% 
  dplyr::select(Assessment, RP, Val) %>% 
  pivot_wider(values_from = Val, names_from =Assessment) %>% 
  pander::pandoc.table()

