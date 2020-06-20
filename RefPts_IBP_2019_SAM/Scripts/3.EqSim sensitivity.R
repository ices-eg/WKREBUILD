# ==================================================================================
# 3. EqSim sensitivity.R
#
# Run the EqSim scenario 5.3 that was selected in IBPHOM 2019
# Test sensitivity to different values of Fcv and Fphi
#
# 18/06/2020 MP
# 19/06/2020 MP
# ==================================================================================

library(msy)
library(FLCore)
library(pander)
library(tidyverse)

root.dir <- "D:/GIT/wk_WKREBUILD"
subfolder <- "RefPts_IBP_2019_SAM"
RData.dir <- file.path(root.dir,subfolder, "RData")
graphics.dir <- file.path(root.dir,subfolder, "Graphics")
logs.dir <- file.path(root.dir,subfolder, "Logs")
source.dir <- file.path(root.dir,subfolder, "Source")

#WGWIDE2018
WG18 <- get(load(file=file.path(RData.dir,"WGWIDE2018.RData")))
name(WG18) <- "WGWIDE18"
WG18df <-
  as.data.frame(WG18) %>% 
  bind_rows(as.data.frame(fbar(WG18)) %>% mutate(slot="fbar")) %>% 
  bind_rows(as.data.frame(ssb(WG18)) %>% mutate(slot="ssb")) %>% 
  mutate(Assessment="WG18")

#WGWIDE2019
# WG19 <- get(load(file=file.path(RData.dir,"WGWIDE2019.RData")))
# name(WG19) <- "WGWIDE19"
# WG19df <-
#   as.data.frame(WG19) %>% 
#   bind_rows(as.data.frame(fbar(WG19)) %>% mutate(slot="fbar")) %>% 
#   bind_rows(as.data.frame(ssb(WG19)) %>% mutate(slot="ssb")) %>% 
#   mutate(Assessment="WG19")

#WGWIDE2018SAM
# WG18SAM <- get(load(file=file.path(RData.dir,"WGWIDE2018SAM.RData")))
# name(WG18SAM) <- "WGWIDE18SAM"
# WG18SAMdf <-
#   as.data.frame(WG18SAM) %>% 
#   bind_rows(as.data.frame(fbar(WG18SAM)) %>% mutate(slot="fbar")) %>% 
#   bind_rows(as.data.frame(ssb(WG18SAM)) %>% mutate(slot="ssb")) %>% 
#   mutate(Assessment="WG18SAM")

#WGWIDE2019SAM
# WG19SAM <- get(load(file=file.path(RData.dir,"WGWIDE2019SAM.RData")))
# name(WG19SAM) <- "WGWIDE19SAM"
# WG19SAMdf <-
#   as.data.frame(WG19SAM) %>% 
#   bind_rows(as.data.frame(fbar(WG19SAM)) %>% mutate(slot="fbar")) %>% 
#   bind_rows(as.data.frame(ssb(WG19SAM)) %>% mutate(slot="ssb")) %>% 
#   mutate(Assessment="WG19SAM")


#plot
# bind_rows(WG18df, WG19df, WG18SAMdf, WG19SAMdf) %>%
#   filter(slot %in% c("ssb","fbar")) %>%
# 
#   ggplot(aes(x=year,y=data)) +
#   theme_bw() +
#   geom_line(aes(colour=Assessment), size=1) +
#   facet_wrap(~slot, scales="free_y")

# each scenario was run without the 1982 data point
# source(file=file.path(getwd(),"Scripts","0.Setup.R"))

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

RPs <- c("Blim","Bpa","Flim","Fpa","MSYBtrigger","FMSY_init","FMSY_interim", "FP05","FMSY_final")
nRPs <- length(RPs)

SRR.models <- c("SegregBlim")
rhologRec <- FALSE
recruitment.trim = c(3, -3)
nsim <- 100

SRR.desc <- paste0(SRR.models,collapse="_")


#fs to scan
#Fscan <- seq(0.01,0.2,len=20) #low res
Fscan <- seq(0.01,0.2,len=40) #high res
Fcvs  <- c(0.1, 0.2, 0.3, 0.4, 0.5)  
Fphis <- c(0.1, 0.2, 0.3, 0.4, 0.5)
# Fcvs  <- c(0.3)  
# Fphis <- c(0.5)

#Base assessments
ass <- c("WG18")
# ass <- c("WG19","WG19SAM")
# ass <- c("WG18", "WG18SAM", "WG19", "WG19SAM")

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
  
  cat("Flim= ",Flim,", Fpa= ",Fpa, "\n")
  
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
  
  Fmsy_init <- SIM1.Fmsy$Refs2["lanF","medianMSY"]
  
  if (Fmsy_init>Fpa) {Fmsy_interim <- Fpa} else {Fmsy_interim <- Fmsy_init}
  
  cat("Fmsy_init= ",Fmsy_init,", Fmsy_interim= ",Fmsy_interim, "\n")
  
  #Step 2 select MSYBtrigger
  #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
  MSYBtrigger <- Bpa
  
  #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
  for (Fphi in Fphis) {
    for (Fcv in Fcvs) {
      cat("Precautionary check FP05. Fcv=",Fcv,", Fphi=",Fphi,"\n")
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
      if (Fmsy_interim>Fp05) {Fmsy_final <- Fp05} else {Fmsy_final <- Fmsy_interim}
      
      cat("Fp05= ",Fp05,", Fmsy_final= ",Fmsy_final, "\n")
      
      cat("**********************************************\n")
      
      dfResults <- 
        dplyr::bind_rows(dfResults,
                        data.frame("Assessment" = rep(a,nRPs),
                                   "Scenario" = rep("5.3",nRPs),
                                   "SRR" = rep(SRR.desc,nRPs),
                                   "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                   "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                   "Fcv" = rep(Fcv,nRPs),
                                   "Fphi" = rep(Fphi,nRPs),
                                   "RP" = RPs,
                                   "Type" = rep("Abs",nRPs),
                                   "Val"  = c(Blim,Bpa,Flim,Fpa,MSYBtrigger,
                                              Fmsy_init,Fmsy_interim, Fp05,Fmsy_final),
                                   "RunTime" = rep(as.character(rt),nRPs),
                                   stringsAsFactors = FALSE))
      
    } # end of Fcvs
  } # end of Fphi 


} #loop over assessments

save(dfResults, file = file.path(RData.dir, "dfResults_sens.RData"))
# load(file = file.path(RData.dir, "dfResults.RData"))

# overview

# dfResults <- bind_rows(dfResults_wg18, dfResults)

# Table of results
# dfResults %>% 
#   group_by(Assessment, Fcv, Fphi, RP) %>% 
#   filter(row_number() == 1) %>% 
#   filter(RP %in% c("Blim","MSYBtrigger","FMSY","FMSY_final","FP05")) %>% 
#   dplyr::select(Assessment, Fcv, Fphi, RP, Val) %>% 
#   pivot_wider(values_from = Val, names_from =Assessment) %>% 
#   pander::pandoc.table()

Fpa <- dfResults %>% filter(RP=="Fpa") %>% distinct(Val) %>% as.numeric()

dfResults %>% 
  filter(RP %in% c("FMSY_final","FMSY_init", "FMSY_interim",  "FP05","Fpa")) %>% 
  mutate(RP = factor(RP, 
                     levels=c("FMSY_init", "Fpa", "FMSY_interim",
                              "FP05","FMSY_final"))) %>% 
  ggplot(aes(x=Fcv, y=Val)) +
  theme_bw() +
  geom_line(aes(colour=as.character(Fphi))) +
  geom_point(aes(colour=as.character(Fphi))) +
  geom_hline(aes(yintercept=Fpa), linetype="dashed") +
  facet_wrap(~RP, ncol=5)



dfResults %>% 
  filter(Fcv == 0.1, Fphi==0.5) %>% 
  dplyr::select(Assessment, Fcv, Fphi, RP, Val) %>% 
  pandoc.table()

if (Fmsy_interim>Fp05) {Fmsy_final <- Fp05} else {Fmsy_final <- Fmsy_interim}
Fmsy_init <- dfResults %>% filter(RP=="FMSY_init", Fcv == 0.1, Fphi==0.5) %>% 
  distinct(Val) %>% as.numeric()
Fmsy_interim <- dfResults %>% filter(RP=="FMSY_interim", Fcv == 0.1, Fphi==0.5) %>% 
  distinct(Val) %>% as.numeric()
Fp05 <- dfResults %>% filter(RP=="FP05", Fcv == 0.1, Fphi==0.5) %>% 
  distinct(Val) %>% as.numeric()

if (Fmsy_interim>Fp05) {Fmsy_final <- Fp05} else {Fmsy_final <- Fmsy_interim}


