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

ac <- function(str){as.character(str)}
fmt <- function(RP,dgt=3,fmt="f"){ac(formatC(RP,digits=dgt,format=fmt,big.mark = ","))}
roundUp <- function(x,to=1000) {to*(x%/%to + as.logical(x%%to))}

subfolder <- "RefPts_IBP_2019_SAM"
RData.dir <- file.path(getwd(),subfolder, "RData")
graphics.dir <- file.path(getwd(),subfolder, "Graphics")
logs.dir <- file.path(getwd(),subfolder, "Logs")
source.dir <- file.path(getwd(),subfolder, "Source")

#WGWIDE2018
WG18 <- get(load(file=file.path(RData.dir,"WGWIDE2018.RData")))
name(WG18) <- "WGWIDE18"
WG18df <-
  as.data.frame(WG18) %>% 
  bind_rows(as.data.frame(fbar(WG18)) %>% mutate(slot="fbar")) %>% 
  bind_rows(as.data.frame(ssb(WG18)) %>% mutate(slot="ssb")) %>% 
  mutate(run="WG18")

#WGWIDE2019
WG19 <- get(load(file=file.path(RData.dir,"WGWIDE2019.RData")))
name(WG19) <- "WGWIDE19"
WG19df <-
  as.data.frame(WG19) %>% 
  bind_rows(as.data.frame(fbar(WG19)) %>% mutate(slot="fbar")) %>% 
  bind_rows(as.data.frame(ssb(WG19)) %>% mutate(slot="ssb")) %>% 
  mutate(run="WG19")

#WGWIDE2019SAM

WG19SAM <- get(load(file=file.path(RData.dir,"WGWIDE2019SAM.RData")))
name(WG19SAM) <- "WGWIDE19SAM"

WG19SAMdf <-
  as.data.frame(WG19SAM) %>% 
  bind_rows(as.data.frame(fbar(WG19SAM)) %>% mutate(slot="fbar")) %>% 
  bind_rows(as.data.frame(ssb(WG19SAM)) %>% mutate(slot="ssb")) %>% 
  mutate(run="WG19SAM")


#plot
# bind_rows(WG18df, WG19df, WG19SAMdf) %>% 
#   filter(slot %in% c("ssb","fbar")) %>% 
#   
#   ggplot(aes(x=year,y=data)) +
#   theme_bw() +
#   geom_line(aes(colour=run)) +
#   facet_wrap(~slot, scales="free_y")

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

#fs to scan
#Fscan <- seq(0.01,0.2,len=20) #low res
Fscan <- seq(0.01,0.2,len=40) #high res

#Base assessments
ass <- c("WG18")
ass <- c("WG19","WG19SAM")
# ass <- c("WG18","WG19","WG19SAM")

# a <- "WG18"

for (a in ass) {

  rt <- format(Sys.time(), "%d-%b-%Y %H.%M")
  
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
                            Btrigger = MSYBtrigger,
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

# save(dfResults, file = file.path(RData.dir, "dfResults.RData"))

# overview

dfResults <- bind_rows(dfResults_wg18, dfResults)

dfResults %>% 
  group_by(Assessment, RP) %>% 
  filter(row_number() == 1) %>% 
  filter(RP %in% c("Blim","MSYBtrigger","Flim","Fpa", "FMSY","FMSY_final","FP05")) %>% 
  dplyr::select(Assessment, RP, Val) %>% 
  pivot_wider(values_from = Val, names_from =Assessment) %>% 
  pander::pandoc.table()


# dfResults_wg18 <- dfResults
