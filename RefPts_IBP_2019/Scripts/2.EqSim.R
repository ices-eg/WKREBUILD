#repeat of WKWIDE 2017 RP analysis with extension to WGWIDE 2017 and 2018 assessments for comparison
#plus additional scenarios

#each scenario was run with and without the 1982 data point
source(file=file.path(getwd(),"Scripts","0.Setup.R"))

RPs <- c("Blim","Bpa","Flim","Fpa","MSYBtrigger","FMSY","FP05","FMSY_final")
nRPs <- length(RPs)

#runName, SRR

#Baseline run, segmented regression at Blim/Bloss
#no recruitment autocorrelation
#default recruitment trimming, selection and biology years (last 10) 
#and assessment/advice #uncertainty/error
runName <- "Baseline"
SRR.models <- c("SegregBlim")
rhologRec <- FALSE
recruitment.trim = c(3, -3)
Fcv <- 0.212; Fphi <- 0.423

#updated baseline wiyh assessment error as supplied by Martin,
#values similar (smaller) than the defaults
#runName <- "AdviceError"
#SRR.models <- c("SegregBlim")
#rhologRec <- FALSE
#recruitment.trim = c(3, -3)
#Fcv <- 0.2; Fphi <- 0.4

#sensitivity test - 3 years for selection
#runName <- "Sens3yrSelection"
#SRR.models <- c("SegregBlim")
#rhologRec <- FALSE
#recruitment.trim = c(3, -3)
#Fcv <- 0.212; Fphi <- 0.423

#sensitivity test - 3 years for biology
#runName <- "Sens3yrBiology"
#SRR.models <- c("SegregBlim")
#rhologRec <- FALSE
#recruitment.trim = c(3, -3)
#Fcv <- 0.212; Fphi <- 0.423

#sentivity test - extreme trimming for high recruitments, leave the low
#runName <- "SensExremeHighRecr"
#SRR.models <- c("SegregBlim")
#rhologRec <- FALSE
#recruitment.trim = c(2, -3)
#Fcv <- 0.212; Fphi <- 0.423

#sensitivity test - recruitment autocorrelation
#runName <- "SensRecrAuto"
#SRR.models <- c("SegregBlim")
#rhologRec <- TRUE
#recruitment.trim = c(3, -3)
#Fcv <- 0.212; Fphi <- 0.423

#Ricker only
#runName <- "RickerSRR"
#SRR.models <- c("Ricker")
#rhologRec <- FALSE
#recruitment.trim = c(3, -3)
#Fcv <- 0.212; Fphi <- 0.423

#Segmented regression only (free breakpoint)
#runName <- "SegRegSRR"
#SRR.models <- c("Segreg")
#rhologRec <- FALSE
#recruitment.trim = c(3, -3)
#Fcv <- 0.212; Fphi <- 0.423

#Beverton Holt only
#runName <- "BevHoltSRR"
#SRR.models <- c("Bevholt")
#rhologRec <- FALSE
#recruitment.trim = c(3, -3)
#Fcv <- 0.212; Fphi <- 0.423

#3 model mix (default with EqSim)
#runName <- "3ModelSRR"
#SRR.models <- c("Ricker","Bevholt","Segreg")
#rhologRec <- FALSE
#recruitment.trim = c(3, -3)
#Fcv <- 0.212; Fphi <- 0.423

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

#cases (1 = Blim set to the lowest SSB that produced a high recruitment, 
#       2 = Bpa set to Bloss,
#       3 = Blim = Bloss,
#       4 = Blim set to the lowest SSB that produced a high recruitment)

#cases <- c(1,2,3,4)
cases <- c(1,2,3,4,5,6)
#cases <- c(3,4)
#cases <- c(5)
cases.desc <- paste0(cases,collapse="_")

#subcases (1 = all years, 2 = ex 1982, 3 = 1985 on)
subcases <- c(1,2,3)
#subcases <- c(2,3)
subcases.desc <- paste0(subcases,collapse="_")

#fs to scan
#Fscan <- seq(0.01,0.2,len=20) #low res
Fscan <- seq(0.01,0.2,len=40) #high res

#Base assessments
ass <- c("WK17","WG17","WG18")

sink(file = file.path(logs.dir,
                      paste0(runName,"_","cases_",cases.desc,"_subcases_",subcases.desc,"_SRR_",SRR.desc,"_",
                 as.character(format(Sys.time(), "%d-%b-%Y %H.%M")),".txt")), 
     split = TRUE)

sigmalnSSB <- 0.194

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
#bio.years <- c(ymax-2,ymax)

#selection years (last 10)
sel.years <- c(ymax-9,ymax)
#selection years (last 3)
#sel.years <- c(ymax-2,ymax)

cat("Biology years",bio.years,"\n")
cat("Selection years",sel.years,"\n")
cat("Mean SSB",mean(ssb(stk)),"\n")
cat("Mean FBar",mean(fbar(stk)),"\n")

#2017 STF 
dfSTF2017 <- read.delim(file = file.path(stf.dir,"STF_WGWIDE2017.csv"),
                        sep=",")
#2018 STF
dfSTF2018 <- read.delim(file = file.path(stf.dir,"STF_WGWIDE2018_IMY1000.csv"),
                        sep=",")
dfSTF2018$ImY <- 1000

for (ImY in seq(2,200)){
  tdfSTF2018 <- read.delim(file = file.path(stf.dir,paste0("STF_WGWIDE2018_IMY",ImY,"000.csv")),
                          sep=",")
  tdfSTF2018$ImY <- 1000*ImY
  dfSTF2018 <- dplyr::bind_rows(dfSTF2018,tdfSTF2018)
}

if (1 %in% cases) {
  
  Blim1 <- min(dfWHM$SSB[dfWHM$Year %in% c(1982,2001) & dfWHM$Assessment==stk.name])
  Bpa1 <- Blim1*exp(1.645*sigmalnSSB)
  #Bpa1 <- 1.4*Blim1
  SegregBlim  <- function(ab, ssb) log(ifelse(ssb >= Blim1, ab$a * Blim1, ab$a * ssb))
  
  ###############################################CASE 1.1###########################################
  #assuming spasmodic recruitment with Blim set to the lowest SSB that produced a high recruitment
  #including 1982
  
  if (1 %in% subcases) {
    
  cat("Case1.1 (all years, including 1982)\n")
  
  #using all SR data 
  cat("Fitting SRR...\n")
  SRR1.1 <- eqsr_fit(stk,
                  remove.years = c(),
                  nsamp=1000, models = SRR.models)
  
  cat("Breakpoint,Mean=",Blim1,SRR1.1$sr.det$a*Blim1,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen1_1_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsr_plot(SRR1.1)
  dev.off()
  
  #Flim
  cat("Flim calculation...\n")
  SIM.Flim1.1 <- eqsim_run(SRR1.1,
                           bio.years = bio.years,
                           bio.const = FALSE,
                           sel.years = sel.years,
                           sel.const = FALSE,
                           Fscan = Fscan,
                           Fcv = 0, Fphi = 0,
                           rhologRec = rhologRec,
                           Blim = Blim1, Bpa = Bpa1,
                           recruitment.trim = recruitment.trim,
                           Btrigger=0,
                           verbose = FALSE)
  
  eqsim_plot(SIM.Flim1.1,catch=TRUE)
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen1_1_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM.Flim1.1,catch=TRUE)
  dev.off()
  
  Flim1.1 <- SIM.Flim1.1$Refs2["catF","F50"]
  #Fpa <- Flim*exp(-1.645*sdlogFBar)
  Fpa1.1 <- Flim1.1/1.4
  #cat("Flim,Fpa=",Flim1.1,Fpa1.1,"\n")
  
  #step 1 - no Btrigger but including stochasticity in 
  #population and fishery and assessment/advice error
  #Fmsy - 1st run, include assessment error
  #for initial FMSY candidate
  
  cat("Initial Fmsy calculation...\n")
  SIM1.Fmsy1.1 <- eqsim_run(SRR1.1,
                         bio.years = bio.years,
                         bio.const = FALSE,
                         sel.years = sel.years,
                         sel.const = FALSE,
                         Fscan = Fscan,
                         Fcv = Fcv, Fphi = Fphi,
                         rhologRec = rhologRec,
                         Blim = Blim1, Bpa = Bpa1,
                         recruitment.trim = recruitment.trim,
                         Btrigger=0,
                         verbose = FALSE)
  
  Fmsy1.1 <- SIM1.Fmsy1.1$Refs2["lanF","medianMSY"]
  #eqsim_plot_range(SIM1.Fmsy1.1, type="median")

  if (Fmsy1.1>Fpa1.1) {Fmsy1.1_final <- Fpa1.1} else {Fmsy1.1_final <- Fmsy1.1}
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen1_1_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM1.Fmsy1.1, type="median")
  dev.off()
  
  #Step 2 select MSYBtrigger
  #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
  MSYBtrigger1 <- Bpa1
  
  #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
  cat("Precautionary check FP05...\n")
  SIM2.Fmsy1.1 <- eqsim_run(SRR1.1,
                            bio.years = bio.years,
                         bio.const = FALSE,
                         sel.years = sel.years,
                         sel.const = FALSE,
                         Fscan = Fscan,
                         Fcv = Fcv, Fphi = Fphi,
                         rhologRec = rhologRec,
                         Blim = Blim1, Bpa = Bpa1,
                         recruitment.trim = recruitment.trim,
                         Btrigger = MSYBtrigger1,
                         verbose = FALSE)
  
  eqsim_plot(SIM2.Fmsy1.1,catch=TRUE)
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen1_1_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM2.Fmsy1.1,catch=TRUE)
  dev.off()
  
  eqsim_plot_range(SIM2.Fmsy1.1, type="median")
  
  Cairo(file = file.path(graphics.dir,runName,paste0("_Scen1_1_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM2.Fmsy1.1, type="median")
  dev.off()
  
  Fp051.1 <- SIM2.Fmsy1.1$Refs2["catF","F05"]
  #cat(Fmsy1.1,Fp051.1,"\n")

  #Final Fmsy
  if (Fmsy1.1_final>Fp051.1) {Fmsy1.1_final <- Fp051.1}
  
  #summary of results
  cat("**********************************************\n")
  cat("Case 1.1, including 1982\n")
  cat("Blim (abs,rel) = ",Blim1,Blim1/mean(ssb(stk)),"\n")
  cat("Bpa (abs,rel) = ",Bpa1,Bpa1/mean(ssb(stk)),"\n")
  cat("Flim (abs,rel) = ",Flim1.1,Flim1.1/mean(fbar(stk)),"\n")
  cat("Fpa (abs,rel) = ",Fpa1.1,Fpa1.1/mean(fbar(stk)),"\n")
  cat("MSYBtrigger (abs,rel) = ",MSYBtrigger1,MSYBtrigger1/mean(ssb(stk)),"\n")
  cat("FMSY (abs,rel) = ",Fmsy1.1,Fmsy1.1/mean(fbar(stk)),"\n")
  cat("FP05 (abs,rel) = ",Fp051.1,Fp051.1/mean(fbar(stk)),"\n")
  cat("FMSY (abs,rel) = ",Fmsy1.1_final,Fmsy1.1_final/mean(fbar(stk)),"\n")
  
  if(a=="WK17"){    
    #2018 F/catch advice
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy1.1_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger1, "\n")
    AdviceF_2018 <- Fmsy1.1_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger1)
    #linear interpolate to get the associated catch
    AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
    
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy1.1_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger1, "\n")
    
    AdviceF_2019 <- Fmsy1.1_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger1){
      AdviceF_2019 <- Fmsy1.1_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger1)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                           y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                           xout=AdviceF_2019,method="linear")$y
    
    cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
    cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
  }
  
  if (a=="WG17") {
    
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy1.1_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger1, "\n")
    
    #2018 F/catch advice
    AdviceF_2018_1.1_WG17 <- Fmsy1.1_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger1)
    #linear interpolate to get the associated catch
    AdviceC_2018_1.1_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_1.1_WG17,method="linear")$y
    cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_1.1_WG17,AdviceC_2018_1.1_WG17,"\n")
  }
  
  if (a=="WG18") {
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018_1.1_WG17/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy1.1_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger1, "\n")
    
    AdviceF_2019_1.1_WG18 <- Fmsy1.1_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger1){
      AdviceF_2019_1.1_WG18 <- Fmsy1.1_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger1)
    }

    #linear interpolate to get the associated catch
    AdviceC_2019_1.1_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                    y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                    xout=AdviceF_2019_1.1_WG18,method="linear")$y
    
    cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_1.1_WG18,AdviceC_2019_1.1_WG18,"\n")
  }
  
  cat("**********************************************\n")
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("1.1",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Abs",nRPs),
                                           "Val" = c(Blim1,Bpa1,Flim1.1,Fpa1.1,MSYBtrigger1,Fmsy1.1,Fp051.1,Fmsy1.1_final),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))

  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("1.1",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Rel",nRPs),
                                           "Val" = c(Blim1/mean(ssb(stk)),
                                                     Bpa1/mean(ssb(stk)),
                                                     Flim1.1/mean(fbar(stk)),
                                                     Fpa1.1/mean(fbar(stk)),
                                                     MSYBtrigger1/mean(ssb(stk)),
                                                     Fmsy1.1/mean(fbar(stk)),
                                                     Fp051.1/mean(fbar(stk)),
                                                     Fmsy1.1_final/mean(fbar(stk))),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  
  }
  
  if (2 %in% subcases) {
    
  ###############################################CASE 1.2###########################################
  #assuming spasmodic recruitment with Blim set to the lowest SSB that produced a high recruitment
  #all years, excluding 1982
  
  cat("Case1.2, all years, excluding 1982\n")
  
  #clean plots
  try(dev.off(),silent=FALSE)
  cat("Fitting SRR...\n")
  SRR1.2 <- eqsr_fit(stk,
                     remove.years = c(1982),
                     nsamp=1000, models = SRR.models)
  
  cat("Breakpoint,Mean=",Blim1,SRR1.2$sr.det$a*Blim1,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen1_2_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsr_plot(SRR1.2)
  dev.off()
  
  #Flim
  cat("Flim calculation...\n")
  SIM.Flim1.2 <- eqsim_run(SRR1.2,
                           bio.years = bio.years,
                           bio.const = FALSE,
                           sel.years = sel.years,
                           sel.const = FALSE,
                           Fscan = Fscan,
                           Fcv = 0, Fphi = 0,
                           rhologRec = rhologRec,
                           Blim = Blim1, Bpa = Bpa1,
                           recruitment.trim = recruitment.trim,
                           Btrigger=0,
                           verbose = FALSE)
  
  #eqsim_plot(SIM.Flim1.2,catch=TRUE)
  Flim1.2 <- SIM.Flim1.2$Refs2["catF","F50"]
  #Fpa <- Flim*exp(-1.645*sdlogFBar)
  Fpa1.2 <- Flim1.2/1.4
  # cat("Flim,Fpa=",Flim1.2,Fpa1.2,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen1_2_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM.Flim1.2,catch=TRUE)
  dev.off()
  
  
  #step 1 - no Btrigger but including stochasticity in 
  #population and fishery and assessment/advice error
  #Fmsy - 1st run, include assessment error
  #for initial FMSY candidate
  cat("Initial Fmsy calculation...\n")
  SIM1.Fmsy1.2 <- eqsim_run(SRR1.2,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim1, Bpa = Bpa1,
                            recruitment.trim = recruitment.trim,
                            Btrigger=0,
                            verbose = FALSE)
  
  Fmsy1.2 <- SIM1.Fmsy1.2$Refs2["lanF","medianMSY"]
  #eqsim_plot_range(SIM1.Fmsy1.2, type="median")

  if (Fmsy1.2>Fpa1.2) {Fmsy1.2_final <- Fpa1.2} else {Fmsy1.2_final <- Fmsy1.2}
  
    
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen1_2_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM1.Fmsy1.2, type="median")
  dev.off()
  
  #Step 2 select MSYBtrigger
  #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
  MSYBtrigger1 <- Bpa1
  
  #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
  cat("Precautionary check FP05...\n")
  SIM2.Fmsy1.2 <- eqsim_run(SRR1.2,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim1, Bpa = Bpa1,
                            recruitment.trim = recruitment.trim,
                            Btrigger = MSYBtrigger1,
                            verbose = FALSE)
  
  #eqsim_plot(SIM2.Fmsy1.2,catch=TRUE)
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen1_2_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM2.Fmsy1.2,catch=TRUE)
  dev.off()
  
  eqsim_plot_range(SIM2.Fmsy1.2, type="median")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen1_2_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM2.Fmsy1.2, type="median")
  dev.off()
  
  Fp051.2 <- SIM2.Fmsy1.2$Refs2["catF","F05"]
  #cat(Fmsy1.2,Fp051.2,"\n")

  #Final Fmsy
  if (Fmsy1.2_final>Fp051.2) {Fmsy1.2_final <- Fp051.2}
  
  #summary of results
  cat("**********************************************\n")
  cat("Case 1.2, excluding 1982\n")
  cat("Blim (abs,rel) = ",Blim1,Blim1/mean(ssb(window(stk,1983,ymax))),"\n")
  cat("Bpa (abs,rel) = ",Bpa1,Bpa1/mean(ssb(window(stk,1983,ymax))),"\n")
  cat("Flim (abs,rel) = ",Flim1.2,Flim1.2/mean(fbar(window(stk,1983,ymax))),"\n")
  cat("Fpa (abs,rel) = ",Fpa1.2,Fpa1.2/mean(fbar(window(stk,1983,ymax))),"\n")
  cat("MSYBtrigger (abs,rel) = ",MSYBtrigger1,MSYBtrigger1/mean(ssb(window(stk,1983,ymax))),"\n")
  cat("FMSY (abs,rel) = ",Fmsy1.2,Fmsy1.2/mean(fbar(window(stk,1983,ymax))),"\n")
  cat("FP05 (abs,rel) = ",Fp051.2,Fp051.2/mean(fbar(window(stk,1983,ymax))),"\n")
  cat("FMSY (abs,rel) = ",Fmsy1.2_final,Fmsy1.2_final/mean(fbar(window(stk,1983,ymax))),"\n")
  
  if(a=="WK17"){    
    #2018 F/catch advice
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy1.2_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger1, "\n")
    AdviceF_2018 <- Fmsy1.2_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger1)
    #linear interpolate to get the associated catch
    AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
    
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy1.2_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger1, "\n")
    
    AdviceF_2019 <- Fmsy1.2_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger1){
      AdviceF_2019 <- Fmsy1.2_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger1)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                           y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                           xout=AdviceF_2019,method="linear")$y
    cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
    cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
  }
  
  if (a=="WG17") {
    
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy1.2_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger1, "\n")
    
    #2018 F/catch advice
    AdviceF_2018_1.2_WG17 <- Fmsy1.2_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger1)
    #linear interpolate to get the associated catch
    AdviceC_2018_1.2_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_1.2_WG17,method="linear")$y
    cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_1.2_WG17,AdviceC_2018_1.2_WG17,"\n")
  }
  
  if (a=="WG18") {
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018_1.2_WG17/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy1.2_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger1, "\n")
    
    AdviceF_2019_1.2_WG18 <- Fmsy1.2_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger1){
      AdviceF_2019_1.2_WG18 <- Fmsy1.2_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger1)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019_1.2_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                    y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                    xout=AdviceF_2019_1.2_WG18,method="linear")$y
    
    cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_1.2_WG18,AdviceC_2019_1.2_WG18,"\n")
  }
  
  
  cat("**********************************************\n")

  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("1.2",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Abs",nRPs),
                                           "Val" = c(Blim1,Bpa1,Flim1.2,Fpa1.2,MSYBtrigger1,Fmsy1.2,Fp051.2,Fmsy1.2_final),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))

  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("1.2",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Rel",nRPs),
                                           "Val" = c(Blim1/mean(ssb(window(stk,1983,ymax))),
                                                     Bpa1/mean(ssb(window(stk,1983,ymax))),
                                                     Flim1.2/mean(fbar(window(stk,1983,ymax))),
                                                     Fpa1.2/mean(fbar(window(stk,1983,ymax))),
                                                     MSYBtrigger1/mean(ssb(window(stk,1983,ymax))),
                                                     Fmsy1.2/mean(fbar(window(stk,1983,ymax))),
                                                     Fp051.2/mean(fbar(window(stk,1983,ymax))),
                                                     Fmsy1.2_final/mean(fbar(window(stk,1983,ymax)))),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))

  }
  
  if (3 %in% subcases) {
    
  ###############################################CASE 1.3###########################################
  #assuming spasmodic recruitment with Blim set to the lowest SSB that produced a high recruitment
  #1995 onwards
  
  cat("Case1.3, 1995 on\n")
  
  #clean plots
  try(dev.off(),silent=TRUE)
  cat("Fitting SRR...\n")
  
  SRR1.3 <- eqsr_fit(window(stk,1995,ymax),
                     remove.years = c(),
                     nsamp=1000, models = SRR.models)
  
  cat("Breakpoint,Mean=",Blim1,SRR1.3$sr.det$a*Blim1,"\n")

  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen1_3_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsr_plot(SRR1.3)
  dev.off()  
  
  #Flim
  cat("Flim calculation...\n")
  SIM.Flim1.3 <- eqsim_run(SRR1.3,
                           bio.years = bio.years,
                           bio.const = FALSE,
                           sel.years = sel.years,
                           sel.const = FALSE,
                           Fscan = Fscan,
                           Fcv = 0, Fphi = 0,
                           rhologRec = rhologRec,
                           Blim = Blim1, Bpa = Bpa1,
                           recruitment.trim = recruitment.trim,
                           Btrigger=0,
                           verbose = FALSE)
  
  #eqsim_plot(SIM.Flim1.3,catch=TRUE)
  Flim1.3 <- SIM.Flim1.3$Refs2["catF","F50"]
  #Fpa <- Flim*exp(-1.645*sdlogFBar)
  Fpa1.3 <- Flim1.3/1.4
  # cat("Flim,Fpa=",Flim1.2,Fpa1.2,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen1_3_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM.Flim1.3,catch=TRUE)
  dev.off()
  
  
  #step 1 - no Btrigger but including stochasticity in 
  #population and fishery and assessment/advice error
  #Fmsy - 1st run, include assessment error
  #for initial FMSY candidate
  cat("Initial Fmsy calculation...\n")
  SIM1.Fmsy1.3 <- eqsim_run(SRR1.3,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim1, Bpa = Bpa1,
                            recruitment.trim = recruitment.trim,
                            Btrigger=0,
                            verbose = FALSE)
  
  Fmsy1.3 <- SIM1.Fmsy1.3$Refs2["lanF","medianMSY"]
  #eqsim_plot_range(SIM1.Fmsy1.3, type="median")
  
  if (Fmsy1.3>Fpa1.3) {Fmsy1.3_final <- Fpa1.3} else {Fmsy1.3_final <- Fmsy1.3}
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen1_3_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM1.Fmsy1.3, type="median")
  dev.off()
  
  #Step 2 select MSYBtrigger
  #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
  MSYBtrigger1 <- Bpa1
  
  #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
  cat("Precautionary check FP05...\n")
  SIM2.Fmsy1.3 <- eqsim_run(SRR1.3,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim1, Bpa = Bpa1,
                            recruitment.trim = recruitment.trim,
                            Btrigger = MSYBtrigger1,
                            verbose = FALSE)
  
  #eqsim_plot(SIM2.Fmsy1.3,catch=TRUE)
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen1_3_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM2.Fmsy1.3,catch=TRUE)
  dev.off()
  
  #eqsim_plot_range(SIM2.Fmsy1.3, type="median")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen1_3_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM2.Fmsy1.3, type="median")
  dev.off()
  
  Fp051.3 <- SIM2.Fmsy1.3$Refs2["catF","F05"]
  #cat(Fmsy1.3,Fp051.3,"\n")
  
  #Final Fmsy
  if (Fmsy1.3_final>Fp051.3) {Fmsy1.3_final <- Fp051.3}
  
  #summary of results
  cat("**********************************************\n")
  cat("Case 1.3, 1995 on\n")
  cat("Blim (abs,rel) = ",Blim1,Blim1/mean(ssb(window(stk,1995,ymax))),"\n")
  cat("Bpa (abs,rel) = ",Bpa1,Bpa1/mean(ssb(window(stk,1995,ymax))),"\n")
  cat("Flim (abs,rel) = ",Flim1.3,Flim1.3/mean(fbar(window(stk,1995,ymax))),"\n")
  cat("Fpa (abs,rel) = ",Fpa1.3,Fpa1.3/mean(fbar(window(stk,1995,ymax))),"\n")
  cat("MSYBtrigger (abs,rel) = ",MSYBtrigger1,MSYBtrigger1/mean(ssb(window(stk,1995,ymax))),"\n")
  cat("FMSY (abs,rel) = ",Fmsy1.3,Fmsy1.3/mean(fbar(window(stk,1995,ymax))),"\n")
  cat("FP05 (abs,rel) = ",Fp051.3,Fp051.3/mean(fbar(window(stk,1995,ymax))),"\n")
  cat("FMSY (abs,rel) = ",Fmsy1.3_final,Fmsy1.3_final/mean(fbar(window(stk,1995,ymax))),"\n")
  
  if(a=="WK17"){    
    #2018 F/catch advice
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy1.3_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger1, "\n")
    AdviceF_2018 <- Fmsy1.3_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger1)
    #linear interpolate to get the associated catch
    AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
    
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy1.3_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger1, "\n")
    
    AdviceF_2019 <- Fmsy1.3_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger1){
      AdviceF_2019 <- Fmsy1.3_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger1)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                           y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                           xout=AdviceF_2019,method="linear")$y
    cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
    cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
  }
  
  if (a=="WG17") {
    
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy1.3_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger1, "\n")
    
    #2018 F/catch advice
    AdviceF_2018_1.3_WG17 <- Fmsy1.3_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger1)
    #linear interpolate to get the associated catch
    AdviceC_2018_1.3_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_1.3_WG17,method="linear")$y
    cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_1.3_WG17,AdviceC_2018_1.3_WG17,"\n")
  }
  
  if (a=="WG18") {
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018_1.3_WG17/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy1.3_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger1, "\n")
    
    AdviceF_2019_1.3_WG18 <- Fmsy1.3_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger1){
      AdviceF_2019_1.3_WG18 <- Fmsy1.3_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger1)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019_1.3_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                    y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                    xout=AdviceF_2019_1.3_WG18,method="linear")$y
    
    cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_1.3_WG18,AdviceC_2019_1.3_WG18,"\n")
  }
  
  cat("**********************************************\n")
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("1.3",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Abs",nRPs),
                                           "Val" = c(Blim1,Bpa1,Flim1.3,Fpa1.3,MSYBtrigger1,Fmsy1.3,Fp051.3,Fmsy1.3_final),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("1.3",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Rel",nRPs),
                                           "Val" = c(Blim1/mean(ssb(window(stk,1995,ymax))),
                                                     Bpa1/mean(ssb(window(stk,1995,ymax))),
                                                     Flim1.3/mean(fbar(window(stk,1995,ymax))),
                                                     Fpa1.3/mean(fbar(window(stk,1995,ymax))),
                                                     MSYBtrigger1/mean(ssb(window(stk,1995,ymax))),
                                                     Fmsy1.3/mean(fbar(window(stk,1995,ymax))),
                                                     Fp051.3/mean(fbar(window(stk,1995,ymax))),
                                                     Fmsy1.3_final/mean(fbar(window(stk,1995,ymax)))),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  }
}


if (2 %in% cases) {
  
  Bpa2 <- min(dfWHM$SSB[dfWHM$Assessment==stk.name])
  Blim2 <- Bpa2/1.4
  #for this scenario it was not considered precautionary to have the segmented regression
  #breakpoint set lower than Bloss (Bpa here)
  SegregBlim  <- function(ab, ssb) log(ifelse(ssb >= Bpa2, ab$a * Bpa2, ab$a * ssb))
  
  if (1 %in% subcases) {
    
  ###############################################CASE 2.1###########################################
  #assuming spasmodic recruitment with Bpa set to Bloss
  #including 1982
  
  cat("Case2.1\n")
  
  cat("Fitting SRR...\n")
  SRR2.1 <- eqsr_fit(stk,
                     remove.years = c(),
                     nsamp=1000, models = SRR.models)

  cat("Breakpoint,Mean=",Blim2,SRR2.1$sr.det$a*Blim2,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen2_1_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsr_plot(SRR2.1)
  dev.off()  
  
  #Flim
  cat("Flim\n")
  SIM.Flim2.1 <- eqsim_run(SRR2.1,
                           bio.years = bio.years,
                           bio.const = FALSE,
                           sel.years = sel.years,
                           sel.const = FALSE,
                           Fscan = Fscan,
                           Fcv = 0, Fphi = 0,
                           rhologRec = rhologRec,
                           Blim = Blim2, Bpa = Bpa2,
                           recruitment.trim = recruitment.trim,
                           Btrigger = 0,
                           verbose = FALSE)
  
  eqsim_plot(SIM.Flim2.1,catch=TRUE)
  Flim2.1 <- SIM.Flim2.1$Refs2["catF","F50"]
  #Fpa <- Flim*exp(-1.645*sdlogFBar)
  Fpa2.1 <- Flim2.1/1.4
  #cat("Flim,Fpa=",Flim2.1,Fpa2.1,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_vWKWIDE2017_Scen2_1_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM.Flim2.1,catch=TRUE)
  dev.off()
  
  
  #step 1 - no Btrigger but including stochasticity in 
  #population and fishery and assessment/advice error
  #Fmsy - 1st run, include assessment error
  #for initial FMSY candidate
  cat("Fmsy\n")
  SIM1.Fmsy2.1 <- eqsim_run(SRR2.1,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim2, Bpa = Bpa2,
                            recruitment.trim = recruitment.trim,
                            Btrigger=0,
                            verbose = FALSE)
  
  Fmsy2.1 <- SIM1.Fmsy2.1$Refs2["lanF","medianMSY"]
  
  #eqsim_plot_range(SIM1.Fmsy2.1, type="median")
  
  if (Fmsy2.1>Fpa2.1) {Fmsy2.1_final <- Fpa2.1} else {Fmsy2.1_final <- Fmsy2.1}
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen2_1_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM1.Fmsy2.1, type="median")
  dev.off()
  
  #Step 2 select MSYBtrigger
  #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
  MSYBtrigger2 <- Bpa2
  
  #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
  cat("FP05\n")
  SIM2.Fmsy2.1 <- eqsim_run(SRR2.1,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim2, Bpa = Bpa2,
                            recruitment.trim = recruitment.trim,
                            Btrigger = MSYBtrigger2,
                            verbose = FALSE)
  
  eqsim_plot(SIM2.Fmsy2.1,catch=TRUE)
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen2_1_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM2.Fmsy2.1,catch=TRUE)
  dev.off()
  
  eqsim_plot_range(SIM2.Fmsy2.1, type="median")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen2_1_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM2.Fmsy2.1, type="median")
  dev.off()
  
  Fp052.1 <- SIM2.Fmsy2.1$Refs2["catF","F05"]
  cat(Fmsy2.1,Fp052.1,"\n")
  
  #Final Fmsy
  if (Fmsy2.1_final>Fp052.1) {Fmsy2.1_final <- Fp052.1}
  
  #summary of results
  cat("**********************************************\n")
  cat("Case 2.1, including 1982\n")
  cat("Blim (abs,rel) = ",Blim2,Blim2/mean(ssb(stk)),"\n")
  cat("Bpa (abs,rel) = ",Bpa2,Bpa2/mean(ssb(stk)),"\n")
  cat("Flim (abs,rel) = ",Flim2.1,Flim2.1/mean(fbar(stk)),"\n")
  cat("Fpa (abs,rel) = ",Fpa2.1,Fpa2.1/mean(fbar(stk)),"\n")
  cat("MSYBtrigger (abs,rel) = ",MSYBtrigger2,MSYBtrigger2/mean(ssb(stk)),"\n")
  cat("FMSY (abs,rel) = ",Fmsy2.1,Fmsy2.1/mean(fbar(stk)),"\n")
  cat("FP05 (abs,rel) = ",Fp052.1,Fp052.1/mean(fbar(stk)),"\n")
  cat("FMSY (abs,rel) = ",Fmsy2.1_final,Fmsy2.1_final/mean(fbar(stk)),"\n")
  
  if(a=="WK17"){    
    #2018 F/catch advice
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy2.1_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger2, "\n")
    AdviceF_2018 <- Fmsy2.1_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger2)
    #linear interpolate to get the associated catch
    AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
    
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy2.1_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger2, "\n")
    
    AdviceF_2019 <- Fmsy2.1_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger2){
      AdviceF_2019 <- Fmsy2.1_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger2)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                           y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                           xout=AdviceF_2019,method="linear")$y
    cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
    cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
  }
  
  if (a=="WG17") {
    
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy2.1_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger2, "\n")
    
    #2018 F/catch advice
    AdviceF_2018_2.1_WG17 <- Fmsy2.1_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger2)
    #linear interpolate to get the associated catch
    AdviceC_2018_2.1_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_2.1_WG17,method="linear")$y
    cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_2.1_WG17,AdviceC_2018_2.1_WG17,"\n")
  }
  
  if (a=="WG18") {
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018_2.1_WG17/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy2.1_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger2, "\n")
    
    AdviceF_2019_2.1_WG18 <- Fmsy2.1_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger2){
      AdviceF_2019_2.1_WG18 <- Fmsy2.1_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger2)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019_2.1_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                    y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                    xout=AdviceF_2019_2.1_WG18,method="linear")$y
    
    cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_2.1_WG18,AdviceC_2019_2.1_WG18,"\n")
  }
  
  cat("**********************************************\n")
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("2.1",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Abs",nRPs),
                                           "Val" = c(Blim2,Bpa2,Flim2.1,Fpa2.1,MSYBtrigger2,Fmsy2.1,Fp052.1,Fmsy2.1_final),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("2.1",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Rel",nRPs),
                                           "Val" = c(Blim2/mean(ssb(stk)),
                                                     Bpa2/mean(ssb(stk)),
                                                     Flim2.1/mean(fbar(stk)),
                                                     Fpa2.1/mean(fbar(stk)),
                                                     MSYBtrigger2/mean(ssb(stk)),
                                                     Fmsy2.1/mean(fbar(stk)),
                                                     Fp052.1/mean(fbar(stk)),
                                                     Fmsy2.1_final/mean(fbar(stk))),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  }
  
  if (2 %in% subcases) {
    
  ###############################################CASE 2.2###########################################
  #assuming spasmodic recruitment with Bpa set to Bloss
  #excluding 1982
  
  cat("Fitting SRR...\n")
  SRR2.2 <- eqsr_fit(stk,
                     remove.years = c(1982),
                     nsamp=1000, models = SRR.models)
  
  #cat("Breakpoint,Mean=",Blim2,SRR2.2$sr.det$a*Blim2,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen2_2_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsr_plot(SRR2.2)
  dev.off()  
  
  #Flim
  cat("Flim\n")
  SIM.Flim2.2 <- eqsim_run(SRR2.2,
                           bio.years = bio.years,
                           bio.const = FALSE,
                           sel.years = sel.years,
                           sel.const = FALSE,
                           Fscan = Fscan,
                           Fcv = 0, Fphi = 0,
                           rhologRec = rhologRec,
                           Blim = Blim2, Bpa = Bpa2,
                           recruitment.trim = recruitment.trim,
                           Btrigger = 0,
                           verbose = FALSE)
  
  eqsim_plot(SIM.Flim2.2,catch=TRUE)
  Flim2.2 <- SIM.Flim2.2$Refs2["catF","F50"]
  #Fpa <- Flim*exp(-1.645*sdlogFBar)
  Fpa2.2 <- Flim2.2/1.4
  #cat("Flim,Fpa=",Flim2.2,Fpa2.2,"\n")
  
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen2_2_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM.Flim2.2,catch=TRUE)
  dev.off()
  
  #step 1 - no Btrigger but including stochasticity in 
  #population and fishery and assessment/advice error
  #Fmsy - 1st run, include assessment error
  #for initial FMSY candidate
  cat("Fmsy\n")
  SIM1.Fmsy2.2 <- eqsim_run(SRR2.2,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim2, Bpa = Bpa2,
                            recruitment.trim = recruitment.trim,
                            Btrigger = 0,
                            verbose = FALSE)
  
  Fmsy2.2 <- SIM1.Fmsy2.2$Refs2["lanF","medianMSY"]
  #eqsim_plot_range(SIM1.Fmsy2.2, type="median")
  
  if (Fmsy2.2>Fpa2.2) {Fmsy2.2_final <- Fpa2.2} else {Fmsy2.2_final <- Fmsy2.2}
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen2_2_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM1.Fmsy2.2, type="median")
  dev.off()
  
  #Step 2 select MSYBtrigger
  #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
  MSYBtrigger2 <- Bpa2
  
  #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
  cat("FP05\n")
  SIM2.Fmsy2.2 <- eqsim_run(SRR2.2,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim2, Bpa = Bpa2,
                            recruitment.trim = recruitment.trim,
                            Btrigger = MSYBtrigger2,
                            verbose = FALSE)
  
  eqsim_plot(SIM2.Fmsy2.2,catch=TRUE)
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen2_2_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM2.Fmsy2.2,catch=TRUE)
  dev.off()
  
  eqsim_plot_range(SIM2.Fmsy2.2, type="median")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen2_2_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM2.Fmsy2.2, type="median")
  dev.off()
  
  Fp052.2 <- SIM2.Fmsy2.2$Refs2["catF","F05"]
  cat(Fmsy2.2,Fp052.2,"\n")
  
  #Final Fmsy
  if (Fmsy2.2_final>Fp052.2) {Fmsy2.2_final <- Fp052.2}
  
  #summary of results
  cat("**********************************************\n")
  cat("Case 2.2, excluding 1982\n")
  cat("Blim (abs,rel) = ",Blim2,Blim2/mean(ssb(window(stk,1983,ymax))),"\n")
  cat("Bpa (abs,rel) = ",Bpa2,Bpa2/mean(ssb(window(stk,1983,ymax))),"\n")
  cat("Flim (abs,rel) = ",Flim2.2,Flim2.2/mean(fbar(window(stk,1983,ymax))),"\n")
  cat("Fpa (abs,rel) = ",Fpa2.2,Fpa2.2/mean(fbar(window(stk,1983,ymax))),"\n")
  cat("MSYBtrigger (abs,rel) = ",MSYBtrigger2,MSYBtrigger2/mean(ssb(window(stk,1983,ymax))),"\n")
  cat("FMSY (abs,rel) = ",Fmsy2.2,Fmsy2.2/mean(fbar(window(stk,1983,ymax))),"\n")
  cat("FP05 (abs,rel) = ",Fp052.2,Fp052.2/mean(fbar(window(stk,1983,ymax))),"\n")
  cat("FMSY (abs,rel) = ",Fmsy2.2_final,Fmsy2.2_final/mean(fbar(window(stk,1983,ymax))),"\n")
  
  if(a=="WK17"){    
    #2018 F/catch advice
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy2.2_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger2, "\n")
    AdviceF_2018 <- Fmsy2.2_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger2)
    #linear interpolate to get the associated catch
    AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
    
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy2.2_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger2, "\n")
    
    AdviceF_2019 <- Fmsy2.2_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger2){
      AdviceF_2019 <- Fmsy2.2_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger2)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                           y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                           xout=AdviceF_2019,method="linear")$y
    cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
    cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
  }
  
  if (a=="WG17") {
    
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy2.2_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger2, "\n")
    
    #2018 F/catch advice
    AdviceF_2018_2.2_WG17 <- Fmsy2.2_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger2)
    #linear interpolate to get the associated catch
    AdviceC_2018_2.2_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_2.2_WG17,method="linear")$y
    cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_2.2_WG17,AdviceC_2018_2.2_WG17,"\n")
  }
  
  if (a=="WG18") {
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018_2.2_WG17/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy2.2_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger2, "\n")
    
    AdviceF_2019_2.2_WG18 <- Fmsy2.2_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger2){
      AdviceF_2019_2.2_WG18 <- Fmsy2.2_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger2)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019_2.2_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                    y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                    xout=AdviceF_2019_2.2_WG18,method="linear")$y
    
    cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_2.2_WG18,AdviceC_2019_2.2_WG18,"\n")
  }
  
  cat("**********************************************\n")
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("2.2",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Abs",nRPs),
                                           "Val" = c(Blim2,Bpa2,Flim2.2,Fpa2.2,MSYBtrigger2,Fmsy2.2,Fp052.2,Fmsy2.2_final),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("2.2",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Rel",nRPs),
                                           "Val" = c(Blim2/mean(ssb(window(stk,1983,ymax))),
                                                     Bpa2/mean(ssb(window(stk,1983,ymax))),
                                                     Flim2.2/mean(fbar(window(stk,1983,ymax))),
                                                     Fpa2.2/mean(fbar(window(stk,1983,ymax))),
                                                     MSYBtrigger2/mean(ssb(window(stk,1983,ymax))),
                                                     Fmsy2.2/mean(fbar(window(stk,1983,ymax))),
                                                     Fp052.2/mean(fbar(window(stk,1983,ymax))),
                                                     Fmsy2.2_final/mean(fbar(window(stk,1983,ymax)))),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  }
  
  if (3 %in% subcases) {
    
  ###############################################CASE 2.3###########################################
  #assuming spasmodic recruitment with Bpa set to Bloss
  #from 1995 onwards
  
  cat("Fitting SRR, 1995 on\n")
  SRR2.3 <- eqsr_fit(window(stk,1995,ymax),
                     remove.years = c(),
                     nsamp=1000, models = SRR.models)
  
  cat("Breakpoint,Mean=",Blim2,SRR2.3$sr.det$a*Blim2,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen2_3_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsr_plot(SRR2.3)
  dev.off()  
  
  #Flim
  cat("Flim calculation...\n")
  SIM.Flim2.3 <- eqsim_run(SRR2.3,
                           bio.years = bio.years,
                           bio.const = FALSE,
                           sel.years = sel.years,
                           sel.const = FALSE,
                           Fscan = Fscan,
                           Fcv = 0, Fphi = 0,
                           rhologRec = rhologRec,
                           Blim = Blim2, Bpa = Bpa2,
                           recruitment.trim = recruitment.trim,
                           Btrigger = 0,
                           verbose = FALSE)
  
  eqsim_plot(SIM.Flim2.3,catch=TRUE)
  Flim2.3 <- SIM.Flim2.3$Refs2["catF","F50"]
  #Fpa <- Flim*exp(-1.645*sdlogFBar)
  Fpa2.3 <- Flim2.3/1.4
  #cat("Flim,Fpa=",Flim2.2,Fpa2.2,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen2_3_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM.Flim2.3,catch=TRUE)
  dev.off()
  
  #step 1 - no Btrigger but including stochasticity in 
  #population and fishery and assessment/advice error
  #Fmsy - 1st run, include assessment error
  #for initial FMSY candidate
  cat("Initial Fmsy calculation...\n")
  SIM1.Fmsy2.3 <- eqsim_run(SRR2.3,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim2, Bpa = Bpa2,
                            recruitment.trim = recruitment.trim,
                            Btrigger = 0,
                            verbose = FALSE)
  
  Fmsy2.3 <- SIM1.Fmsy2.3$Refs2["lanF","medianMSY"]
  #eqsim_plot_range(SIM1.Fmsy2.3, type="median")
  
  if (Fmsy2.3>Fpa2.3) {Fmsy2.3_final <- Fpa2.3} else {Fmsy2.3_final <- Fmsy2.3}
  
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen2_3_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM1.Fmsy2.3, type="median")
  dev.off()
  
  #Step 2 select MSYBtrigger
  #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
  MSYBtrigger2 <- Bpa2
  
  #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
  cat("Precautionary check FP05...\n")
  SIM2.Fmsy2.3 <- eqsim_run(SRR2.3,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim2, Bpa = Bpa2,
                            recruitment.trim = recruitment.trim,
                            Btrigger = MSYBtrigger2,
                            verbose = FALSE)
  
  eqsim_plot(SIM2.Fmsy2.3,catch=TRUE)
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen2_3_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM2.Fmsy2.3,catch=TRUE)
  dev.off()
  
  eqsim_plot_range(SIM2.Fmsy2.3, type="median")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen2_3_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM2.Fmsy2.3, type="median")
  dev.off()
  
  Fp052.3 <- SIM2.Fmsy2.3$Refs2["catF","F05"]
  cat(Fmsy2.3,Fp052.3,"\n")
  
  #Final Fmsy
  if (Fmsy2.3_final>Fp052.3) {Fmsy2.3_final <- Fp052.3}
  
  #summary of results
  cat("**********************************************\n")
  cat("Case 2.3, from 1995\n")
  cat("Blim (abs,rel) = ",Blim2,Blim2/mean(ssb(window(stk,1995,ymax))),"\n")
  cat("Bpa (abs,rel) = ",Bpa2,Bpa2/mean(ssb(window(stk,1995,ymax))),"\n")
  cat("Flim (abs,rel) = ",Flim2.3,Flim2.3/mean(fbar(window(stk,1995,ymax))),"\n")
  cat("Fpa (abs,rel) = ",Fpa2.3,Fpa2.3/mean(fbar(window(stk,1995,ymax))),"\n")
  cat("MSYBtrigger (abs,rel) = ",MSYBtrigger2,MSYBtrigger2/mean(ssb(window(stk,1995,ymax))),"\n")
  cat("FMSY (abs,rel) = ",Fmsy2.3,Fmsy2.3/mean(fbar(window(stk,1995,ymax))),"\n")
  cat("FP05 (abs,rel) = ",Fp052.3,Fp052.3/mean(fbar(window(stk,1995,ymax))),"\n")
  cat("FMSY (abs,rel) = ",Fmsy2.3_final,Fmsy2.3_final/mean(fbar(window(stk,1995,ymax))),"\n")
  
  if(a=="WK17"){    
    #2018 F/catch advice
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy2.3_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger2, "\n")
    AdviceF_2018 <- Fmsy2.3_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger2)
    #linear interpolate to get the associated catch
    AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
    
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy2.3_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger2, "\n")
    
    AdviceF_2019 <- Fmsy2.3_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger2){
      AdviceF_2019 <- Fmsy2.3_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger2)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                           y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                           xout=AdviceF_2019,method="linear")$y
    cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
    cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
  }
  
  if (a=="WG17") {
    
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy2.3_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger2, "\n")
    
    #2018 F/catch advice
    AdviceF_2018_2.3_WG17 <- Fmsy2.3_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger2)
    #linear interpolate to get the associated catch
    AdviceC_2018_2.3_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_2.3_WG17,method="linear")$y
    cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_2.3_WG17,AdviceC_2018_2.3_WG17,"\n")
  }
  
  if (a=="WG18") {
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018_2.3_WG17/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy2.3_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger2, "\n")
    
    AdviceF_2019_2.3_WG18 <- Fmsy2.3_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger2){
      AdviceF_2019_2.3_WG18 <- Fmsy2.3_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger2)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019_2.3_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                    y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                    xout=AdviceF_2019_2.3_WG18,method="linear")$y
    
    cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_2.3_WG18,AdviceC_2019_2.3_WG18,"\n")
  }
  
  cat("**********************************************\n")
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("2.3",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Abs",nRPs),
                                           "Val" = c(Blim2,Bpa2,Flim2.3,Fpa2.3,MSYBtrigger2,Fmsy2.3,Fp052.3,Fmsy2.3_final),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("2.3",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Rel",nRPs),
                                           "Val" = c(Blim2/mean(ssb(window(stk,1995,ymax))),
                                                     Bpa2/mean(ssb(window(stk,1995,ymax))),
                                                     Flim2.3/mean(fbar(window(stk,1995,ymax))),
                                                     Fpa2.3/mean(fbar(window(stk,1995,ymax))),
                                                     MSYBtrigger2/mean(ssb(window(stk,1995,ymax))),
                                                     Fmsy2.3/mean(fbar(window(stk,1995,ymax))),
                                                     Fp052.3/mean(fbar(window(stk,1995,ymax))),
                                                     Fmsy2.3_final/mean(fbar(window(stk,1995,ymax)))),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  }
}

if (3 %in% cases) {
  
  Blim3 <- min(dfWHM$SSB[dfWHM$Assessment==stk.name])
  Bpa3 <- Blim3*exp(1.645*sigmalnSSB)
  #Bpa3 <- 1.4*Blim3
  SegregBlim  <- function(ab, ssb) log(ifelse(ssb >= Blim3, ab$a * Blim3, ab$a * ssb))
  
  if (1 %in% subcases) {
    
  ###############################################CASE 3.1###########################################
  #no SRR with Blim=Bloss, including 1982
  
  cat("Case3.1\n")
  
  #using all SR data 
  cat("Fitting SRR...\n")
  SRR3.1 <- eqsr_fit(stk,
                     remove.years = c(),
                     nsamp=1000, models = SRR.models)
  
  cat("Breakpoint,Mean=",Blim3,SRR3.1$sr.det$a*Blim3,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen3_1_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsr_plot(SRR3.1)
  dev.off()  
  
  #Flim
  cat("Flim calculation...\n")
  SIM.Flim3.1 <- eqsim_run(SRR3.1,
                           bio.years = bio.years,
                           bio.const = FALSE,
                           sel.years = sel.years,
                           sel.const = FALSE,
                           Fscan = Fscan,
                           Fcv = 0, Fphi = 0,
                           rhologRec = rhologRec,
                           Blim = Blim3, Bpa = Bpa3,
                           recruitment.trim = recruitment.trim,
                           Btrigger = 0,
                           verbose = FALSE)
  
  #eqsim_plot(SIM.Flim3.1,catch=TRUE)
  Flim3.1 <- SIM.Flim3.1$Refs2["catF","F50"]
  #Fpa <- Flim*exp(-1.645*sdlogFBar)
  Fpa3.1 <- Flim3.1/1.4
  #cat("Flim,Fpa=",Flim3.1,Fpa3.1,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen3_1_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM.Flim3.1,catch=TRUE)
  dev.off()
  
  #step 1 - no Btrigger but including stochasticity in 
  #population and fishery and assessment/advice error
  #Fmsy - 1st run, include assessment error
  #for initial FMSY candidate
  
  cat("Initial Fmsy calculation...\n")
  SIM1.Fmsy3.1 <- eqsim_run(SRR3.1,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim3, Bpa = Bpa3,
                            recruitment.trim = recruitment.trim,
                            Btrigger=0,
                            verbose = FALSE)
  
  Fmsy3.1 <- SIM1.Fmsy3.1$Refs2["lanF","medianMSY"]
  #eqsim_plot_range(SIM1.Fmsy3.1, type="median")
  
  if (Fmsy3.1>Fpa3.1) {Fmsy3.1_final <- Fpa3.1} else {Fmsy3.1_final <- Fmsy3.1}
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen3_1_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM1.Fmsy3.1, type="median")
  dev.off()
  
  #Step 2 select MSYBtrigger
  #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
  MSYBtrigger3 <- Bpa3
  
  cat("FP05\n")
  #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
  SIM2.Fmsy3.1 <- eqsim_run(SRR3.1,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim3, Bpa = Bpa3,
                            recruitment.trim = recruitment.trim,
                            Btrigger = MSYBtrigger3,
                            verbose = FALSE)
  
  #eqsim_plot(SIM2.Fmsy3.1,catch=TRUE)
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen3_1_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM2.Fmsy3.1,catch=TRUE)
  dev.off()
  
  #eqsim_plot_range(SIM2.Fmsy3.1, type="median")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen3_1_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM2.Fmsy3.1, type="median")
  dev.off()
  
  Fp053.1 <- SIM2.Fmsy3.1$Refs2["catF","F05"]
  #cat(Fmsy3.1,Fp053.1,"\n")
  
  #Final Fmsy
  if (Fmsy3.1_final>Fp053.1) {Fmsy3.1_final <- Fp053.1}
  
  #summary of results
  cat("**********************************************\n")
  cat("Case 3.1, including 1982\n")
  cat("Blim (abs,rel) = ",Blim3,Blim3/mean(ssb(stk)),"\n")
  cat("Bpa (abs,rel) = ",Bpa3,Bpa3/mean(ssb(stk)),"\n")
  cat("Flim (abs,rel) = ",Flim3.1,Flim3.1/mean(fbar(stk)),"\n")
  cat("Fpa (abs,rel) = ",Fpa3.1,Fpa3.1/mean(fbar(stk)),"\n")
  cat("MSYBtrigger (abs,rel) = ",MSYBtrigger3,MSYBtrigger3/mean(ssb(stk)),"\n")
  cat("FMSY (abs,rel) = ",Fmsy3.1,Fmsy3.1/mean(fbar(stk)),"\n")
  cat("FP05 (abs,rel) = ",Fp053.1,Fp053.1/mean(fbar(stk)),"\n")
  cat("FMSY (abs,rel) = ",Fmsy3.1_final,Fmsy3.1_final/mean(fbar(stk)),"\n")
  
  if(a=="WK17"){    
    #2018 F/catch advice
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy3.1_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger3, "\n")
    AdviceF_2018 <- Fmsy3.1_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger3)
    #linear interpolate to get the associated catch
    AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
    
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy3.1_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger3, "\n")
    
    AdviceF_2019 <- Fmsy3.1_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger3){
      AdviceF_2019 <- Fmsy3.1_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger3)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                           y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                           xout=AdviceF_2019,method="linear")$y
    cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
    cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
  }
  
  if (a=="WG17") {
    
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy3.1_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger3, "\n")
    
    #2018 F/catch advice
    AdviceF_2018_3.1_WG17 <- Fmsy3.1_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger3)
    #linear interpolate to get the associated catch
    AdviceC_2018_3.1_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_3.1_WG17,method="linear")$y
    cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_3.1_WG17,AdviceC_2018_3.1_WG17,"\n")
  }
  
  if (a=="WG18") {
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018_3.1_WG17/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy3.1_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger3, "\n")
    
    AdviceF_2019_3.1_WG18 <- Fmsy3.1_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger3){
      AdviceF_2019_3.1_WG18 <- Fmsy3.1_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger3)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019_3.1_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                    y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                    xout=AdviceF_2019_3.1_WG18,method="linear")$y
    
    cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_3.1_WG18,AdviceC_2019_3.1_WG18,"\n")
  }
  
  cat("**********************************************\n")
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("3.1",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Abs",nRPs),
                                           "Val" = c(Blim3,Bpa3,Flim3.1,Fpa3.1,MSYBtrigger3,Fmsy3.1,Fp053.1,Fmsy3.1_final),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("3.1",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Rel",nRPs),
                                           "Val" = c(Blim3/mean(ssb(stk)),
                                                     Bpa3/mean(ssb(stk)),
                                                     Flim3.1/mean(fbar(stk)),
                                                     Fpa3.1/mean(fbar(stk)),
                                                     MSYBtrigger3/mean(ssb(stk)),
                                                     Fmsy3.1/mean(fbar(stk)),
                                                     Fp053.1/mean(fbar(stk)),
                                                     Fmsy3.1_final/mean(fbar(stk))),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  }
  
  if (2 %in% subcases) {
    
  ###############################################CASE 3.2###########################################
  #no SRR with Blim=Bloss, excluding 1982
  
  #using all SR data, except 1982
  cat("Fitting SRR, excluding 1982...\n")
  SRR3.2 <- eqsr_fit(stk,
                     remove.years = c(1982),
                     nsamp=1000, models = SRR.models)
  
  cat("Breakpoint,Mean=",Blim3,SRR3.2$sr.det$a*Blim3,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen3_2_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsr_plot(SRR3.2)
  dev.off()  
  
  #Flim
  cat("Flim calculation...\n")
  SIM.Flim3.2 <- eqsim_run(SRR3.2,
                           bio.years = bio.years,
                           bio.const = FALSE,
                           sel.years = sel.years,
                           sel.const = FALSE,
                           Fscan = Fscan,
                           Fcv = 0, Fphi = 0,
                           rhologRec = rhologRec,
                           Blim = Blim3, Bpa = Bpa3,
                           recruitment.trim = recruitment.trim,
                           Btrigger = 0,
                           verbose = FALSE)
  
  #eqsim_plot(SIM.Flim3.2,catch=TRUE)
  Flim3.2 <- SIM.Flim3.2$Refs2["catF","F50"]
  #Fpa <- Flim*exp(-1.645*sdlogFBar)
  Fpa3.2 <- Flim3.2/1.4
  #cat("Flim,Fpa=",Flim3.2,Fpa3.2,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen3_2_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM.Flim3.2,catch=TRUE)
  dev.off()
  
  
  #step 1 - no Btrigger but including stochasticity in 
  #population and fishery and assessment/advice error
  #Fmsy - 1st run, include assessment error
  #for initial FMSY candidate
  
  cat("Initial Fmsy calculation...\n")
  SIM1.Fmsy3.2 <- eqsim_run(SRR3.2,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim3, Bpa = Bpa3,
                            recruitment.trim = recruitment.trim,
                            Btrigger=0,
                            verbose = FALSE)
  
  Fmsy3.2 <- SIM1.Fmsy3.2$Refs2["lanF","medianMSY"]
  
  if (Fmsy3.2>Fpa3.2) {Fmsy3.2_final <- Fpa3.2} else {Fmsy3.2_final <- Fmsy3.2}
  
  #eqsim_plot(SIM1.Fmsy3.2)
  #eqsim_plot_range(SIM1.Fmsy3.2, type="median")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen3_2_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM1.Fmsy3.2, type="median")
  dev.off()
  
  #Step 2 select MSYBtrigger
  #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
  MSYBtrigger3 <- Bpa3
  
  cat("Precautionary check Fp05...\n")
  #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
  SIM2.Fmsy3.2 <- eqsim_run(SRR3.2,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim3, Bpa = Bpa3,
                            recruitment.trim = recruitment.trim,
                            Btrigger = MSYBtrigger3,
                            verbose = FALSE)
  
  #eqsim_plot(SIM2.Fmsy3.2,catch=TRUE)
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen3_2_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM2.Fmsy3.2,catch=TRUE)
  dev.off()
  
  #eqsim_plot_range(SIM2.Fmsy3.2, type="median")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen3_2_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM2.Fmsy3.2, type="median")
  dev.off()
  
  Fp053.2 <- SIM2.Fmsy3.2$Refs2["catF","F05"]
  #cat(Fmsy3.2,Fp053.2,"\n")

  #Final Fmsy
  if (Fmsy3.2_final>Fp053.2) {Fmsy3.2_final <- Fp053.2}
  
  cat("**********************************************\n")
  cat("Case 3.2, excluding 1982\n")
  cat("Blim (abs,rel) = ",Blim3,Blim3/mean(ssb(window(stk,1983,ymax))),"\n")
  cat("Bpa (abs,rel) = ",Bpa3,Bpa3/mean(ssb(window(stk,1983,ymax))),"\n")
  cat("Flim (abs,rel) = ",Flim3.2,Flim3.2/mean(fbar(window(stk,1983,ymax))),"\n")
  cat("Fpa (abs,rel) = ",Fpa3.2,Fpa3.2/mean(fbar(window(stk,1983,ymax))),"\n")
  cat("MSYBtrigger (abs,rel) = ",MSYBtrigger3,MSYBtrigger3/mean(ssb(window(stk,1983,ymax))),"\n")
  cat("FMSY (abs,rel) = ",Fmsy3.2,Fmsy3.2/mean(fbar(window(stk,1983,ymax))),"\n")
  cat("FP05 (abs,rel) = ",Fp053.2,Fp053.2/mean(fbar(window(stk,1983,ymax))),"\n")
  cat("FMSY (abs,rel) = ",Fmsy3.2_final,Fmsy3.2_final/mean(fbar(window(stk,1983,ymax))),"\n")
  
  if(a=="WK17"){    
    #2018 F/catch advice
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy3.2_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger3, "\n")
    AdviceF_2018 <- Fmsy3.2_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger3)
    #linear interpolate to get the associated catch
    AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
    
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy3.2_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger3, "\n")
    
    AdviceF_2019 <- Fmsy3.2_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger3){
      AdviceF_2019 <- Fmsy3.2_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger3)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                           y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                           xout=AdviceF_2019,method="linear")$y
    cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
    cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
  }
  
  if (a=="WG17") {
    
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy3.2_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger3, "\n")
    
    #2018 F/catch advice
    AdviceF_2018_3.2_WG17 <- Fmsy3.2_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger3)
    #linear interpolate to get the associated catch
    AdviceC_2018_3.2_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_3.2_WG17,method="linear")$y
    cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_3.2_WG17,AdviceC_2018_3.2_WG17,"\n")
  }
  
  if (a=="WG18") {
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018_3.2_WG17/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy3.2_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger3, "\n")
    
    AdviceF_2019_3.2_WG18 <- Fmsy3.2_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger3){
      AdviceF_2019_3.2_WG18 <- Fmsy3.2_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger3)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019_3.2_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                             y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                             xout=AdviceF_2019_3.2_WG18,method="linear")$y
    
    cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_3.2_WG18,AdviceC_2019_3.2_WG18,"\n")
  }
  
  cat("**********************************************\n")
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("3.2",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Abs",nRPs),
                                           "Val" = c(Blim3,Bpa3,Flim3.2,Fpa3.2,MSYBtrigger3,Fmsy3.2,Fp053.2,Fmsy3.2_final),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("3.2",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Rel",nRPs),
                                           "Val" = c(Blim3/mean(ssb(window(stk,1983,ymax))),
                                                     Bpa3/mean(ssb(window(stk,1983,ymax))),
                                                     Flim3.2/mean(fbar(window(stk,1983,ymax))),
                                                     Fpa3.2/mean(fbar(window(stk,1983,ymax))),
                                                     MSYBtrigger3/mean(ssb(window(stk,1983,ymax))),
                                                     Fmsy3.2/mean(fbar(window(stk,1983,ymax))),
                                                     Fp053.2/mean(fbar(window(stk,1983,ymax))),
                                                     Fmsy3.2_final/mean(fbar(window(stk,1983,ymax)))),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
 
  } 
  
  if (3 %in% subcases) {
    
  ###############################################CASE 3.3###########################################
  #no SRR with Blim=Bloss, 1995 on
  
  #using all SR data, except 1982
  cat("Fitting SRR, 1995 on...\n")
  SRR3.3 <- eqsr_fit(window(stk,1995,ymax),
                     remove.years = c(),
                     nsamp=1000, models = SRR.models)
  
  cat("Breakpoint,Mean=",Blim3,SRR3.3$sr.det$a*Blim3,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen3_3_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsr_plot(SRR3.3)
  dev.off()  
  
  #Flim
  cat("Flim calculation...\n")
  SIM.Flim3.3 <- eqsim_run(SRR3.3,
                           bio.years = bio.years,
                           bio.const = FALSE,
                           sel.years = sel.years,
                           sel.const = FALSE,
                           Fscan = Fscan,
                           Fcv = 0, Fphi = 0,
                           rhologRec = rhologRec,
                           Blim = Blim3, Bpa = Bpa3,
                           recruitment.trim = recruitment.trim,
                           Btrigger = 0,
                           verbose = FALSE)
  
  #eqsim_plot(SIM.Flim3.3,catch=TRUE)
  Flim3.3 <- SIM.Flim3.3$Refs2["catF","F50"]
  #Fpa <- Flim*exp(-1.645*sdlogFBar)
  Fpa3.3 <- Flim3.3/1.4
  #cat("Flim,Fpa=",Flim3.3,Fpa3.3,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen3_3_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM.Flim3.3,catch=TRUE)
  dev.off()
  
  
  #step 1 - no Btrigger but including stochasticity in 
  #population and fishery and assessment/advice error
  #Fmsy - 1st run, include assessment error
  #for initial FMSY candidate
  
  cat("Initial Fmsy calculation...\n")
  SIM1.Fmsy3.3 <- eqsim_run(SRR3.3,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim3, Bpa = Bpa3,
                            recruitment.trim = recruitment.trim,
                            Btrigger=0,
                            verbose = FALSE)
  
  Fmsy3.3 <- SIM1.Fmsy3.3$Refs2["lanF","medianMSY"]
  #eqsim_plot(SIM1.Fmsy3.3)
  #eqsim_plot_range(SIM1.Fmsy3.3, type="median")
  
  if (Fmsy3.3>Fpa3.3) {Fmsy3.3_final <- Fpa3.3} else {Fmsy3.3_final <- Fmsy3.3}
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen3_3_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM1.Fmsy3.3, type="median")
  dev.off()
  
  #Step 2 select MSYBtrigger
  #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
  MSYBtrigger3 <- Bpa3
  
  cat("Precautionary check Fp05...\n")
  #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
  SIM2.Fmsy3.3 <- eqsim_run(SRR3.3,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim3, Bpa = Bpa3,
                            recruitment.trim = recruitment.trim,
                            Btrigger = MSYBtrigger3,
                            verbose = FALSE)
  
  #eqsim_plot(SIM2.Fmsy3.2,catch=TRUE)
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen3_3_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM2.Fmsy3.3,catch=TRUE)
  dev.off()
  
  #eqsim_plot_range(SIM2.Fmsy3.3, type="median")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen3_3_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM2.Fmsy3.3, type="median")
  dev.off()
  
  Fp053.3 <- SIM2.Fmsy3.3$Refs2["catF","F05"]
  #cat(Fmsy3.3,Fp053.3,"\n")

  #Final Fmsy
  if (Fmsy3.3_final>Fp053.3) {Fmsy3.3_final <- Fp053.3}
  
  cat("**********************************************\n")
  cat("Case 3.3, 1995 on\n")
  cat("Blim (abs,rel) = ",Blim3,Blim3/mean(ssb(window(stk,1995,ymax))),"\n")
  cat("Bpa (abs,rel) = ",Bpa3,Bpa3/mean(ssb(window(stk,1995,ymax))),"\n")
  cat("Flim (abs,rel) = ",Flim3.3,Flim3.3/mean(fbar(window(stk,1995,ymax))),"\n")
  cat("Fpa (abs,rel) = ",Fpa3.3,Fpa3.3/mean(fbar(window(stk,1995,ymax))),"\n")
  cat("MSYBtrigger (abs,rel) = ",MSYBtrigger3,MSYBtrigger3/mean(ssb(window(stk,1995,ymax))),"\n")
  cat("FMSY (abs,rel) = ",Fmsy3.3,Fmsy3.3/mean(fbar(window(stk,1995,ymax))),"\n")
  cat("FP05 (abs,rel) = ",Fp053.3,Fp053.3/mean(fbar(window(stk,1995,ymax))),"\n")
  cat("FMSY (abs,rel) = ",Fmsy3.3_final,Fmsy3.3_final/mean(fbar(window(stk,1995,ymax))),"\n")
  
  if(a=="WK17"){    
    #2018 F/catch advice
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy3.3_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger3, "\n")
    AdviceF_2018 <- Fmsy3.3_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger3)
    #linear interpolate to get the associated catch
    AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
    
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy3.3_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger3, "\n")
    
    AdviceF_2019 <- Fmsy3.3_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger3){
      AdviceF_2019 <- Fmsy3.3_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger3)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                           y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                           xout=AdviceF_2019,method="linear")$y
    cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
    cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
  }
  
  if (a=="WG17") {
    
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy3.3_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger3, "\n")
    
    #2018 F/catch advice
    AdviceF_2018_3.3_WG17 <- Fmsy3.3_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger3)
    #linear interpolate to get the associated catch
    AdviceC_2018_3.3_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_3.3_WG17,method="linear")$y
    cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_3.3_WG17,AdviceC_2018_3.3_WG17,"\n")
  }
  
  if (a=="WG18") {
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018_3.3_WG17/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy3.3_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger3, "\n")
    
    AdviceF_2019_3.3_WG18 <- Fmsy3.3_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger3){
      AdviceF_2019_3.3_WG18 <- Fmsy3.3_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger3)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019_3.3_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                    y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                    xout=AdviceF_2019_3.3_WG18,method="linear")$y
    
    cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_3.3_WG18,AdviceC_2019_3.3_WG18,"\n")
  }
  
  cat("**********************************************\n")
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("3.3",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Abs",nRPs),
                                           "Val" = c(Blim3,Bpa3,Flim3.3,Fpa3.3,MSYBtrigger3,Fmsy3.3,Fp053.3,Fmsy3.3_final),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("3.3",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Rel",nRPs),
                                           "Val" = c(Blim3/mean(ssb(window(stk,1995,ymax))),
                                                     Bpa3/mean(ssb(window(stk,1995,ymax))),
                                                     Flim3.3/mean(fbar(window(stk,1995,ymax))),
                                                     Fpa3.3/mean(fbar(window(stk,1995,ymax))),
                                                     MSYBtrigger3/mean(ssb(window(stk,1995,ymax))),
                                                     Fmsy3.3/mean(fbar(window(stk,1995,ymax))),
                                                     Fp053.3/mean(fbar(window(stk,1995,ymax))),
                                                     Fmsy3.3_final/mean(fbar(window(stk,1995,ymax)))),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  }
  
}

if (4 %in% cases) {

  ###############################################CASE 4###########################################
  #Bpa = SSB for 2001
  #including 1982
  
  cat("Case4.1\n")
  
  Bpa4 <- dfWHM$SSB[dfWHM$Year==2001 & dfWHM$Assessment==stk.name]
  Blim4 <- Bpa4/1.4
  SegregBlim  <- function(ab, ssb) log(ifelse(ssb >= Blim4, ab$a * Blim4, ab$a * ssb))
  
  if (1 %in% subcases) {
    
  #using all SR data 
  cat("Fitting SRR, all data...\n")
  SRR4.1 <- eqsr_fit(stk,
                     remove.years = c(),
                     nsamp=1000, models = SRR.models)
  
  cat("Breakpoint,Mean=",Blim4,SRR4.1$sr.det$a*Blim4,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen4_1_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsr_plot(SRR4.1)
  dev.off()  
  
  #Flim
  cat("Flim calculation...\n")
  SIM.Flim4.1 <- eqsim_run(SRR4.1,
                           bio.years = bio.years,
                           bio.const = FALSE,
                           sel.years = sel.years,
                           sel.const = FALSE,
                           Fscan = Fscan,
                           Fcv = 0, Fphi = 0,
                           rhologRec = rhologRec,
                           Blim = Blim4, Bpa = Bpa4,
                           recruitment.trim = recruitment.trim,
                           Btrigger = 0,
                           verbose = FALSE)
  
  eqsim_plot(SIM.Flim4.1,catch=TRUE)
  Flim4.1 <- SIM.Flim4.1$Refs2["catF","F50"]
  #Fpa <- Flim*exp(-1.645*sdlogFBar)
  Fpa4.1 <- Flim4.1/1.4
  #cat("Flim,Fpa=",Flim4.1,Fpa4.1,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen4_1_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM.Flim4.1,catch=TRUE)
  dev.off()
  
  #step 1 - no Btrigger but including stochasticity in 
  #population and fishery and assessment/advice error
  #Fmsy - 1st run, include assessment error
  #for initial FMSY candidate
  cat("Initial Fmsy calculation...\n")
  SIM1.Fmsy4.1 <- eqsim_run(SRR4.1,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim4, Bpa = Bpa4,
                            recruitment.trim = recruitment.trim,
                            Btrigger = 0,
                            verbose = FALSE)
  
  Fmsy4.1 <- SIM1.Fmsy4.1$Refs2["lanF","medianMSY"]
  #eqsim_plot_range(SIM1.Fmsy4.1, type="median")
  
  if (Fmsy4.1>Fpa4.1) {Fmsy4.1_final <- Fpa4.1} else {Fmsy4.1_final <- Fmsy4.1}
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen4_1_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM1.Fmsy4.1, type="median")
  dev.off()
  
  #Step 2 select MSYBtrigger
  #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
  MSYBtrigger4 <- Bpa4
  
  #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
  cat("Precautionary check FP05...\n")
  SIM2.Fmsy4.1 <- eqsim_run(SRR4.1,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim4, Bpa = Bpa4,
                            recruitment.trim = recruitment.trim,
                            Btrigger = MSYBtrigger4,
                            verbose = FALSE)
  
  eqsim_plot(SIM2.Fmsy4.1,catch=TRUE)
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen4_1_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM2.Fmsy4.1,catch=TRUE)
  dev.off()
  
  eqsim_plot_range(SIM2.Fmsy4.1, type="median")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen4_1_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM2.Fmsy4.1, type="median")
  dev.off()
  
  Fp054.1 <- SIM2.Fmsy4.1$Refs2["catF","F05"]
  #cat(Fmsy4.1,Fp054.1,"\n")
  
  #Final Fmsy
  if (Fmsy4.1_final>Fp054.1) {Fmsy4.1_final <- Fp054.1}
  
  #summary of results
  cat("**********************************************\n")
  cat("Case 4.1, including 1982\n")
  cat("Blim (abs,rel) = ",Blim4,Blim4/mean(ssb(stk)),"\n")
  cat("Bpa (abs,rel) = ",Bpa4,Bpa4/mean(ssb(stk)),"\n")
  cat("Flim (abs,rel) = ",Flim4.1,Flim4.1/mean(fbar(stk)),"\n")
  cat("Fpa (abs,rel) = ",Fpa4.1,Fpa4.1/mean(fbar(stk)),"\n")
  cat("MSYBtrigger (abs,rel) = ",MSYBtrigger4,MSYBtrigger4/mean(ssb(stk)),"\n")
  cat("FMSY (abs,rel) = ",Fmsy4.1,Fmsy4.1/mean(fbar(stk)),"\n")
  cat("FP05 (abs,rel) = ",Fp054.1,Fp054.1/mean(fbar(stk)),"\n")
  cat("FMSY (abs,rel) = ",Fmsy4.1_final,Fmsy4.1_final/mean(fbar(stk)),"\n")
  
  if(a=="WK17"){    
    #2018 F/catch advice
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy4.1_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger4, "\n")
    AdviceF_2018 <- Fmsy4.1_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger4)
    #linear interpolate to get the associated catch
    AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
    
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy4.1_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger4, "\n")
    
    AdviceF_2019 <- Fmsy4.1_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger4){
      AdviceF_2019 <- Fmsy4.1_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger4)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                           y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                           xout=AdviceF_2019,method="linear")$y
    cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
    cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
  }
  
  if (a=="WG17") {
    
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy4.1_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger4, "\n")
    
    #2018 F/catch advice
    AdviceF_2018_4.1_WG17 <- Fmsy4.1_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger4)
    #linear interpolate to get the associated catch
    AdviceC_2018_4.1_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_4.1_WG17,method="linear")$y
    cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_4.1_WG17,AdviceC_2018_4.1_WG17,"\n")
  }
  
  if (a=="WG18") {
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018_4.1_WG17/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy4.1_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger4, "\n")
    
    AdviceF_2019_4.1_WG18 <- Fmsy4.1_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger4){
      AdviceF_2019_4.1_WG18 <- Fmsy4.1_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger4)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019_4.1_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                    y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                    xout=AdviceF_2019_4.1_WG18,method="linear")$y
    
    cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_4.1_WG18,AdviceC_2019_4.1_WG18,"\n")
  }
  
  cat("**********************************************\n")
  
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("4.1",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Abs",nRPs),
                                           "Val" = c(Blim4,Bpa4,Flim4.1,Fpa4.1,MSYBtrigger4,Fmsy4.1,Fp054.1,Fmsy4.1_final),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("4.1",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Rel",nRPs),
                                           "Val" = c(Blim4/mean(ssb(stk)),
                                                     Bpa4/mean(ssb(stk)),
                                                     Flim4.1/mean(fbar(stk)),
                                                     Fpa4.1/mean(fbar(stk)),
                                                     MSYBtrigger4/mean(ssb(stk)),
                                                     Fmsy4.1/mean(fbar(stk)),
                                                     Fp054.1/mean(fbar(stk)),
                                                     Fmsy4.1_final/mean(fbar(stk))),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  }
  
  if (2 %in% subcases) {
    
  ###############################################CASE 4.2###########################################
  #Bpa = SSB for 2001
  #excluding 1982
  cat("Case 4.2\n")
  cat("Fitting SRR, excluding 1982...\n")
  SRR4.2 <- eqsr_fit(stk,
                     remove.years = c(1982),
                     nsamp=1000, models = SRR.models)
  
  cat("Breakpoint,Mean=",Blim4,SRR4.2$sr.det$a*Blim4,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen4_2_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsr_plot(SRR4.2)
  dev.off()  
  
  #Flim
  cat("Flim calculation...\n")
  SIM.Flim4.2 <- eqsim_run(SRR4.2,
                           bio.years = bio.years,
                           bio.const = FALSE,
                           sel.years = sel.years,
                           sel.const = FALSE,
                           Fscan = Fscan,
                           Fcv = 0, Fphi = 0,
                           rhologRec = rhologRec,
                           Blim = Blim4, Bpa = Bpa4,
                           recruitment.trim = recruitment.trim,
                           Btrigger = 0,
                           verbose = FALSE)
  
  eqsim_plot(SIM.Flim4.2,catch=TRUE)
  Flim4.2 <- SIM.Flim4.2$Refs2["catF","F50"]
  #Fpa <- Flim*exp(-1.645*sdlogFBar)
  Fpa4.2 <- Flim4.2/1.4
  #cat("Flim,Fpa=",Flim4.2,Fpa4.2,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen4_2_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM.Flim4.2,catch=TRUE)
  dev.off()
  
  #step 1 - no Btrigger but including stochasticity in 
  #population and fishery and assessment/advice error
  #Fmsy - 1st run, include assessment error
  #for initial FMSY candidate
  cat("Initial Fmsy calculation...\n")
  SIM1.Fmsy4.2 <- eqsim_run(SRR4.2,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim4, Bpa = Bpa4,
                            recruitment.trim = recruitment.trim,
                            Btrigger = 0,
                            verbose = FALSE)
  
  Fmsy4.2 <- SIM1.Fmsy4.2$Refs2["lanF","medianMSY"]
  #eqsim_plot_range(SIM1.Fmsy4.2, type="median")
  
  if (Fmsy4.2>Fpa4.2) {Fmsy4.2_final <- Fpa4.2} else {Fmsy4.2_final <- Fmsy4.2}
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen4_2_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM1.Fmsy4.2, type="median")
  dev.off()
  
  #Step 2 select MSYBtrigger
  #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
  MSYBtrigger4 <- Bpa4
  
  #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
  cat("Precautionary check FP05...\n")
  SIM2.Fmsy4.2 <- eqsim_run(SRR4.2,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim4, Bpa = Bpa4,
                            recruitment.trim = recruitment.trim,
                            Btrigger = MSYBtrigger4,
                            verbose = FALSE)
  
  eqsim_plot(SIM2.Fmsy4.2,catch=TRUE)
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen4_2_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM2.Fmsy4.2,catch=TRUE)
  dev.off()
  
  eqsim_plot_range(SIM2.Fmsy4.2, type="median")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen4_2_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM2.Fmsy4.2, type="median")
  dev.off()
  
  Fp054.2 <- SIM2.Fmsy4.2$Refs2["catF","F05"]
  #cat(Fmsy4.2,Fp054.2,"\n")

  #Final Fmsy
  if (Fmsy4.2_final>Fp054.2) {Fmsy4.2_final <- Fp054.2}
  
  #summary of results
  cat("**********************************************\n")
  cat("Case 4.2, excluding 1982\n")
  cat("Blim (abs,rel) = ",Blim4,Blim4/mean(ssb(window(stk,1983,ymax))),"\n")
  cat("Bpa (abs,rel) = ",Bpa4,Bpa4/mean(ssb(window(stk,1983,ymax))),"\n")
  cat("Flim (abs,rel) = ",Flim4.2,Flim4.2/mean(fbar(window(stk,1983,ymax))),"\n")
  cat("Fpa (abs,rel) = ",Fpa4.2,Fpa4.2/mean(fbar(window(stk,1983,ymax))),"\n")
  cat("MSYBtrigger (abs,rel) = ",MSYBtrigger4,MSYBtrigger4/mean(ssb(window(stk,1983,ymax))),"\n")
  cat("FMSY (abs,rel) = ",Fmsy4.2,Fmsy4.2/mean(fbar(window(stk,1983,ymax))),"\n")
  cat("FP05 (abs,rel) = ",Fp054.2,Fp054.2/mean(fbar(window(stk,1983,ymax))),"\n")
  cat("FMSY (abs,rel) = ",Fmsy4.2_final,Fmsy4.2_final/mean(fbar(window(stk,1983,ymax))),"\n")
  
  if(a=="WK17"){    
    #2018 F/catch advice
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy4.2_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger4, "\n")
    AdviceF_2018 <- Fmsy4.2_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger4)
    #linear interpolate to get the associated catch
    AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
    
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy4.2_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger4, "\n")
    
    AdviceF_2019 <- Fmsy4.2_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger4){
      AdviceF_2019 <- Fmsy4.2_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger4)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                           y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                           xout=AdviceF_2019,method="linear")$y
    cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
    cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
  }
  
  if (a=="WG17") {
    
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy4.2_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger4, "\n")
    
    #2018 F/catch advice
    AdviceF_2018_4.2_WG17 <- Fmsy4.2_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger4)
    #linear interpolate to get the associated catch
    AdviceC_2018_4.2_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_4.2_WG17,method="linear")$y
    cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_4.2_WG17,AdviceC_2018_4.2_WG17,"\n")
  }
  
  if (a=="WG18") {
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018_4.2_WG17/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy4.2_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger4, "\n")
    
    AdviceF_2019_4.2_WG18 <- Fmsy4.2_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger4){
      AdviceF_2019_4.2_WG18 <- Fmsy4.2_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger4)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019_4.2_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                    y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                    xout=AdviceF_2019_4.2_WG18,method="linear")$y
    
    cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_4.2_WG18,AdviceC_2019_4.2_WG18,"\n")
  }
  
  cat("**********************************************\n")
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("4.2",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Abs",nRPs),
                                           "Val" = c(Blim4,Bpa4,Flim4.2,Fpa4.2,MSYBtrigger4,Fmsy4.2,Fp054.2,Fmsy4.2_final),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("4.2",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Rel",nRPs),
                                           "Val" = c(Blim4/mean(ssb(window(stk,1983,ymax))),
                                                     Bpa4/mean(ssb(window(stk,1983,ymax))),
                                                     Flim4.2/mean(fbar(window(stk,1983,ymax))),
                                                     Fpa4.2/mean(fbar(window(stk,1983,ymax))),
                                                     MSYBtrigger4/mean(ssb(window(stk,1983,ymax))),
                                                     Fmsy4.2/mean(fbar(window(stk,1983,ymax))),
                                                     Fp054.2/mean(fbar(window(stk,1983,ymax))),
                                                     Fmsy4.2_final/mean(fbar(window(stk,1983,ymax)))),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  
  }
  
  if (3 %in% subcases) {
    
  ###############################################CASE 4.3###########################################
  #Bpa = SSB for 2001, 1995 on
  
  cat("Case 4.3\n")
  cat("Fitting SRR, 1995 on...\n")
  
  SRR4.3 <- eqsr_fit(window(stk,1995,ymax),
                     remove.years = c(),
                     nsamp=1000, models = SRR.models)
  
  cat("Breakpoint,Mean=",Blim4,SRR4.3$sr.det$a*Blim4,"\n")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen4_3_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsr_plot(SRR4.3)
  dev.off()  
  
  #Flim
  cat("Flim calculation...\n")
  SIM.Flim4.3 <- eqsim_run(SRR4.3,
                           bio.years = bio.years,
                           bio.const = FALSE,
                           sel.years = sel.years,
                           sel.const = FALSE,
                           Fscan = Fscan,
                           Fcv = 0, Fphi = 0,
                           rhologRec = rhologRec,
                           Blim = Blim4, Bpa = Bpa4,
                           recruitment.trim = recruitment.trim,
                           Btrigger = 0,
                           verbose = FALSE)
  
  eqsim_plot(SIM.Flim4.3,catch=TRUE)
  Flim4.3 <- SIM.Flim4.3$Refs2["catF","F50"]
  #Fpa <- Flim*exp(-1.645*sdlogFBar)
  Fpa4.3 <- Flim4.3/1.4
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen4_3_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM.Flim4.3,catch=TRUE)
  dev.off()
  
  #step 1 - no Btrigger but including stochasticity in 
  #population and fishery and assessment/advice error
  #Fmsy - 1st run, include assessment error
  #for initial FMSY candidate
  cat("Initial Fmsy calculation...\n")
  SIM1.Fmsy4.3 <- eqsim_run(SRR4.3,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim4, Bpa = Bpa4,
                            recruitment.trim = recruitment.trim,
                            Btrigger = 0,
                            verbose = FALSE)
  
  Fmsy4.3 <- SIM1.Fmsy4.3$Refs2["lanF","medianMSY"]
  #eqsim_plot_range(SIM1.Fmsy4.3, type="median")
  
  if (Fmsy4.3>Fpa4.3) {Fmsy4.3_final <- Fpa4.3} else {Fmsy4.3_final <- Fmsy4.3}
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen4_3_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM1.Fmsy4.3, type="median")
  dev.off()
  
  #Step 2 select MSYBtrigger
  #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
  MSYBtrigger4 <- Bpa4
  
  #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
  cat("Precautionary check FP05...\n")
  SIM2.Fmsy4.3 <- eqsim_run(SRR4.3,
                            bio.years = bio.years,
                            bio.const = FALSE,
                            sel.years = sel.years,
                            sel.const = FALSE,
                            Fscan = Fscan,
                            Fcv = Fcv, Fphi = Fphi,
                            rhologRec = rhologRec,
                            Blim = Blim4, Bpa = Bpa4,
                            recruitment.trim = recruitment.trim,
                            Btrigger = MSYBtrigger4,
                            verbose = FALSE)
  
  eqsim_plot(SIM2.Fmsy4.3,catch=TRUE)
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen4_3_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot(SIM2.Fmsy4.3,catch=TRUE)
  dev.off()
  
  eqsim_plot_range(SIM2.Fmsy4.3, type="median")
  
  Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen4_3_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
  eqsim_plot_range(SIM2.Fmsy4.3, type="median")
  dev.off()
  
  Fp054.3 <- SIM2.Fmsy4.3$Refs2["catF","F05"]
  #cat(Fmsy4.2,Fp054.2,"\n")

  #Final Fmsy
  if (Fmsy4.3_final>Fp054.3) {Fmsy4.3_final <- Fp054.3}
  
  #summary of results
  cat("**********************************************\n")
  cat("Case 4.3, 1995 on\n")
  cat("Blim (abs,rel) = ",Blim4,Blim4/mean(ssb(window(stk,1995,ymax))),"\n")
  cat("Bpa (abs,rel) = ",Bpa4,Bpa4/mean(ssb(window(stk,1995,ymax))),"\n")
  cat("Flim (abs,rel) = ",Flim4.3,Flim4.3/mean(fbar(window(stk,1995,ymax))),"\n")
  cat("Fpa (abs,rel) = ",Fpa4.3,Fpa4.3/mean(fbar(window(stk,1995,ymax))),"\n")
  cat("MSYBtrigger (abs,rel) = ",MSYBtrigger4,MSYBtrigger4/mean(ssb(window(stk,1995,ymax))),"\n")
  cat("FMSY (abs,rel) = ",Fmsy4.3,Fmsy4.3/mean(fbar(window(stk,1995,ymax))),"\n")
  cat("FP05 (abs,rel) = ",Fp054.3,Fp054.3/mean(fbar(window(stk,1995,ymax))),"\n")
  cat("FMSY (abs,rel) = ",Fmsy4.3_final,Fmsy4.3_final/mean(fbar(window(stk,1995,ymax))),"\n")
  
  if(a=="WK17"){    
    #2018 F/catch advice
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy4.3_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger4, "\n")
    AdviceF_2018 <- Fmsy4.3_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger4)
    #linear interpolate to get the associated catch
    AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
    
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy4.3_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger4, "\n")
    
    AdviceF_2019 <- Fmsy4.3_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger4){
      AdviceF_2019 <- Fmsy4.3_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger4)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                           y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                           xout=AdviceF_2019,method="linear")$y
    cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
    cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
  }
  
  if (a=="WG17") {
    
    cat("2018 adviceF calculation\n")
    cat("Fmsy = ",Fmsy4.3_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger4, "\n")
    
    #2018 F/catch advice
    AdviceF_2018_4.3_WG17 <- Fmsy4.3_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger4)
    #linear interpolate to get the associated catch
    AdviceC_2018_4.3_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_4.3_WG17,method="linear")$y
    cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_4.3_WG17,AdviceC_2018_4.3_WG17,"\n")
  }
  
  if (a=="WG18") {
    #2019 F/catch advice
    #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
    ImYC <- 1000*round(AdviceC_2018_4.3_WG17/1000)
    #2019 SSB
    dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
    
    cat("2019 adviceF calculation\n")
    cat("Fmsy = ",Fmsy4.3_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger4, "\n")
    
    AdviceF_2019_4.3_WG18 <- Fmsy4.3_final
    if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger4){
      AdviceF_2019_4.3_WG18 <- Fmsy4.3_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger4)
    }
    
    #linear interpolate to get the associated catch
    AdviceC_2019_4.3_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                    y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                    xout=AdviceF_2019_4.3_WG18,method="linear")$y
    
    cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_4.3_WG18,AdviceC_2019_4.3_WG18,"\n")
  }
  
  cat("**********************************************\n")
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("4.3",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Abs",nRPs),
                                           "Val" = c(Blim4,Bpa4,Flim4.3,Fpa4.3,MSYBtrigger4,Fmsy4.3,Fp054.3,Fmsy4.3_final),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  
  dfResults <- dplyr::bind_rows(dfResults,
                                data.frame("Assessment" = rep(a,nRPs),
                                           "Scenario" = rep("4.3",nRPs),
                                           "SRR" = rep(SRR.desc,nRPs),
                                           "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                           "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                           "Fcv" = rep(Fcv,nRPs),
                                           "Fphi" = rep(Fphi,nRPs),
                                           "RP" = RPs,
                                           "Type" = rep("Rel",nRPs),
                                           "Val" = c(Blim4/mean(ssb(window(stk,1995,ymax))),
                                                     Bpa4/mean(ssb(window(stk,1995,ymax))),
                                                     Flim4.3/mean(fbar(window(stk,1995,ymax))),
                                                     Fpa4.3/mean(fbar(window(stk,1995,ymax))),
                                                     MSYBtrigger4/mean(ssb(window(stk,1995,ymax))),
                                                     Fmsy4.3/mean(fbar(window(stk,1995,ymax))),
                                                     Fp054.3/mean(fbar(window(stk,1995,ymax))),
                                                     Fmsy4.3_final/mean(fbar(window(stk,1995,ymax)))),
                                           "RunTime" = rep(as.character(rt),nRPs),
                                           stringsAsFactors = FALSE))
  
  }
  
  }


if (5 %in% cases) {
  
  ###############################################CASE 5###########################################
  #Bpa = SSB for 2003
  #including 1982
  
  cat("Case5.1\n")
  
  Bpa5 <- dfWHM$SSB[dfWHM$Year==2003 & dfWHM$Assessment==stk.name]
  Blim5 <- Bpa5/1.4
  SegregBlim  <- function(ab, ssb) log(ifelse(ssb >= Blim5, ab$a * Blim5, ab$a * ssb))
  
  #check that the Blim is not lower than Bloss, if so, set the breakpoint to Bloss
  #this should only be the case with the WG18 assessment
  Bloss5<-as.numeric(min(ssb(window(stk,1983,ymax))))
  
  if (Blim5<Bloss5) {
    cat("modifying SegregBlim as Bloss>Blim for ",a,"\n")
    cat("Blim5 = ", Blim5, "Bloss5 = ", Bloss5,"\n")
    SegregBlim  <- function(ab, ssb) log(ifelse(ssb >= Bloss5, ab$a * Bloss5, ab$a * ssb))
  }
  
  if (1 %in% subcases) {
    
    #using all SR data 
    cat("Fitting SRR, all data...\n")
    SRR5.1 <- eqsr_fit(stk,
                       remove.years = c(),
                       nsamp=1000, models = SRR.models)
    
    cat("Breakpoint,Mean=",Blim5,SRR5.1$sr.det$a*Blim5,"\n")
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen5_1_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsr_plot(SRR5.1)
    dev.off()  
    
    #Flim
    cat("Flim calculation...\n")
    SIM.Flim5.1 <- eqsim_run(SRR5.1,
                             bio.years = bio.years,
                             bio.const = FALSE,
                             sel.years = sel.years,
                             sel.const = FALSE,
                             Fscan = Fscan,
                             Fcv = 0, Fphi = 0,
                             rhologRec = rhologRec,
                             Blim = Blim5, Bpa = Bpa5,
                             recruitment.trim = recruitment.trim,
                             Btrigger = 0,
                             verbose = FALSE)
    
    eqsim_plot(SIM.Flim5.1,catch=TRUE)
    Flim5.1 <- SIM.Flim5.1$Refs2["catF","F50"]
    Fpa5.1 <- Flim5.1/1.4

    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen5_1_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot(SIM.Flim5.1,catch=TRUE)
    dev.off()
    
    #step 1 - no Btrigger but including stochasticity in 
    #population and fishery and assessment/advice error
    #Fmsy - 1st run, include assessment error
    #for initial FMSY candidate
    cat("Initial Fmsy calculation...\n")
    SIM1.Fmsy5.1 <- eqsim_run(SRR5.1,
                              bio.years = bio.years,
                              bio.const = FALSE,
                              sel.years = sel.years,
                              sel.const = FALSE,
                              Fscan = Fscan,
                              Fcv = Fcv, Fphi = Fphi,
                              rhologRec = rhologRec,
                              Blim = Blim5, Bpa = Bpa5,
                              recruitment.trim = recruitment.trim,
                              Btrigger = 0,
                              verbose = FALSE)
    
    Fmsy5.1 <- SIM1.Fmsy5.1$Refs2["lanF","medianMSY"]
    
    if (Fmsy5.1>Fpa5.1) {Fmsy5.1_final <- Fpa5.1} else {Fmsy5.1_final <- Fmsy5.1}
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen5_1_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot_range(SIM1.Fmsy5.1, type="median")
    dev.off()
    
    #Step 2 select MSYBtrigger
    #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
    MSYBtrigger5 <- Bpa5
    
    #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
    cat("Precautionary check FP05...\n")
    SIM2.Fmsy5.1 <- eqsim_run(SRR5.1,
                              bio.years = bio.years,
                              bio.const = FALSE,
                              sel.years = sel.years,
                              sel.const = FALSE,
                              Fscan = Fscan,
                              Fcv = Fcv, Fphi = Fphi,
                              rhologRec = rhologRec,
                              Blim = Blim5, Bpa = Bpa5,
                              recruitment.trim = recruitment.trim,
                              Btrigger = MSYBtrigger5,
                              verbose = FALSE)
    
    eqsim_plot(SIM2.Fmsy5.1,catch=TRUE)
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen5_1_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot(SIM2.Fmsy5.1,catch=TRUE)
    dev.off()
    
    eqsim_plot_range(SIM2.Fmsy5.1, type="median")
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen5_1_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot_range(SIM2.Fmsy5.1, type="median")
    dev.off()
    
    Fp055.1 <- SIM2.Fmsy5.1$Refs2["catF","F05"]
    
    #Final Fmsy
    if (Fmsy5.1_final>Fp055.1) {Fmsy5.1_final <- Fp055.1}
    
    #summary of results
    cat("**********************************************\n")
    cat("Case 5.1, including 1982\n")
    cat("Blim (abs,rel) = ",Blim5,Blim5/mean(ssb(stk)),"\n")
    cat("Bpa (abs,rel) = ",Bpa5,Bpa5/mean(ssb(stk)),"\n")
    cat("Flim (abs,rel) = ",Flim5.1,Flim5.1/mean(fbar(stk)),"\n")
    cat("Fpa (abs,rel) = ",Fpa5.1,Fpa5.1/mean(fbar(stk)),"\n")
    cat("MSYBtrigger (abs,rel) = ",MSYBtrigger5,MSYBtrigger5/mean(ssb(stk)),"\n")
    cat("FMSY (abs,rel) = ",Fmsy5.1,Fmsy5.1/mean(fbar(stk)),"\n")
    cat("FP05 (abs,rel) = ",Fp055.1,Fp055.1/mean(fbar(stk)),"\n")
    cat("FMSY (abs,rel) = ",Fmsy5.1_final,Fmsy5.1_final/mean(fbar(stk)),"\n")
    
    if(a=="WK17"){    
      #2018 F/catch advice
      cat("2018 adviceF calculation\n")
      cat("Fmsy = ",Fmsy5.1_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger5, "\n")
      AdviceF_2018 <- Fmsy5.1_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger5)
      #linear interpolate to get the associated catch
      AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
      
      #2019 F/catch advice
      #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
      ImYC <- 1000*round(AdviceC_2018/1000)
      #2019 SSB
      dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
      
      cat("2019 adviceF calculation\n")
      cat("Fmsy = ",Fmsy5.1_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger5, "\n")
      
      AdviceF_2019 <- Fmsy5.1_final
      if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger5){
        AdviceF_2019 <- Fmsy5.1_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger5)
      }
      
      #linear interpolate to get the associated catch
      AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                             y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                             xout=AdviceF_2019,method="linear")$y
      cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
      cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
    }
    
    if (a=="WG17") {
      
      cat("2018 adviceF calculation\n")
      cat("Fmsy = ",Fmsy5.1_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger5, "\n")
      
      #2018 F/catch advice
      AdviceF_2018_5.1_WG17 <- Fmsy5.1_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger5)
      #linear interpolate to get the associated catch
      AdviceC_2018_5.1_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_5.1_WG17,method="linear")$y
      cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_5.1_WG17,AdviceC_2018_5.1_WG17,"\n")
    }
    
    if (a=="WG18") {
      #2019 F/catch advice
      #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
      ImYC <- 1000*round(AdviceC_2018_5.1_WG17/1000)
      #2019 SSB
      dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
      
      cat("2019 adviceF calculation\n")
      cat("Fmsy = ",Fmsy5.1_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger5, "\n")
      
      AdviceF_2019_5.1_WG18 <- Fmsy5.1_final
      if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger5){
        AdviceF_2019_5.1_WG18 <- Fmsy5.1_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger5)
      }
      
      #linear interpolate to get the associated catch
      AdviceC_2019_5.1_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                      y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                      xout=AdviceF_2019_5.1_WG18,method="linear")$y
      
      cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_5.1_WG18,AdviceC_2019_5.1_WG18,"\n")
    }
    
    cat("**********************************************\n")
    
    
    dfResults <- dplyr::bind_rows(dfResults,
                                  data.frame("Assessment" = rep(a,nRPs),
                                             "Scenario" = rep("5.1",nRPs),
                                             "SRR" = rep(SRR.desc,nRPs),
                                             "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                             "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                             "Fcv" = rep(Fcv,nRPs),
                                             "Fphi" = rep(Fphi,nRPs),
                                             "RP" = RPs,
                                             "Type" = rep("Abs",nRPs),
                                             "Val" = c(Blim5,Bpa5,Flim5.1,Fpa5.1,MSYBtrigger5,Fmsy5.1,Fp055.1,Fmsy5.1_final),
                                             "RunTime" = rep(as.character(rt),nRPs),
                                             stringsAsFactors = FALSE))
    
    dfResults <- dplyr::bind_rows(dfResults,
                                  data.frame("Assessment" = rep(a,nRPs),
                                             "Scenario" = rep("5.1",nRPs),
                                             "SRR" = rep(SRR.desc,nRPs),
                                             "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                             "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                             "Fcv" = rep(Fcv,nRPs),
                                             "Fphi" = rep(Fphi,nRPs),
                                             "RP" = RPs,
                                             "Type" = rep("Rel",nRPs),
                                             "Val" = c(Blim5/mean(ssb(stk)),
                                                       Bpa5/mean(ssb(stk)),
                                                       Flim5.1/mean(fbar(stk)),
                                                       Fpa5.1/mean(fbar(stk)),
                                                       MSYBtrigger5/mean(ssb(stk)),
                                                       Fmsy5.1/mean(fbar(stk)),
                                                       Fp055.1/mean(fbar(stk)),
                                                       Fmsy5.1_final/mean(fbar(stk))),
                                             "RunTime" = rep(as.character(rt),nRPs),
                                             stringsAsFactors = FALSE))
  }
  
  if (2 %in% subcases) {
    
    ###############################################CASE 5.2###########################################
    #Bpa = SSB for 2003
    #excluding 1982
    cat("Case 5.2\n")
    cat("Fitting SRR, excluding 1982...\n")
    SRR5.2 <- eqsr_fit(stk,
                       remove.years = c(1982),
                       nsamp=1000, models = SRR.models)
    
    cat("Breakpoint,Mean=",Blim5,SRR5.2$sr.det$a*Blim5,"\n")
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen5_2_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsr_plot(SRR5.2)
    dev.off()  
    
    #Flim
    cat("Flim calculation...\n")
    SIM.Flim5.2 <- eqsim_run(SRR5.2,
                             bio.years = bio.years,
                             bio.const = FALSE,
                             sel.years = sel.years,
                             sel.const = FALSE,
                             Fscan = Fscan,
                             Fcv = 0, Fphi = 0,
                             rhologRec = rhologRec,
                             Blim = Blim5, Bpa = Bpa5,
                             recruitment.trim = recruitment.trim,
                             Btrigger = 0,
                             verbose = FALSE)
    
    eqsim_plot(SIM.Flim5.2,catch=TRUE)
    Flim5.2 <- SIM.Flim5.2$Refs2["catF","F50"]
    Fpa5.2 <- Flim5.2/1.4
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen5_2_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot(SIM.Flim5.2,catch=TRUE)
    dev.off()
    
    #step 1 - no Btrigger but including stochasticity in 
    #population and fishery and assessment/advice error
    #Fmsy - 1st run, include assessment error
    #for initial FMSY candidate
    cat("Initial Fmsy calculation...\n")
    SIM1.Fmsy5.2 <- eqsim_run(SRR5.2,
                              bio.years = bio.years,
                              bio.const = FALSE,
                              sel.years = sel.years,
                              sel.const = FALSE,
                              Fscan = Fscan,
                              Fcv = Fcv, Fphi = Fphi,
                              rhologRec = rhologRec,
                              Blim = Blim5, Bpa = Bpa5,
                              recruitment.trim = recruitment.trim,
                              Btrigger = 0,
                              verbose = FALSE)
    
    Fmsy5.2 <- SIM1.Fmsy5.2$Refs2["lanF","medianMSY"]
    
    if (Fmsy5.2>Fpa5.2) {Fmsy5.2_final <- Fpa5.2} else {Fmsy5.2_final <- Fmsy5.2}
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen5_2_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot_range(SIM1.Fmsy5.2, type="median")
    dev.off()
    
    #Step 2 select MSYBtrigger
    #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
    MSYBtrigger5 <- Bpa5
    
    #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
    cat("Precautionary check FP05...\n")
    SIM2.Fmsy5.2 <- eqsim_run(SRR5.2,
                              bio.years = bio.years,
                              bio.const = FALSE,
                              sel.years = sel.years,
                              sel.const = FALSE,
                              Fscan = Fscan,
                              Fcv = Fcv, Fphi = Fphi,
                              rhologRec = rhologRec,
                              Blim = Blim5, Bpa = Bpa5,
                              recruitment.trim = recruitment.trim,
                              Btrigger = MSYBtrigger5,
                              verbose = FALSE)
    
    eqsim_plot(SIM2.Fmsy5.2,catch=TRUE)
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen5_2_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot(SIM2.Fmsy5.2,catch=TRUE)
    dev.off()
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen5_2_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot_range(SIM2.Fmsy5.2, type="median")
    dev.off()
    
    Fp055.2 <- SIM2.Fmsy5.2$Refs2["catF","F05"]
    
    #Final Fmsy
    if (Fmsy5.2_final>Fp055.2) {Fmsy5.2_final <- Fp055.2}
    
    #summary of results
    cat("**********************************************\n")
    cat("Case 5.2, excluding 1982\n")
    cat("Blim (abs,rel) = ",Blim5,Blim5/mean(ssb(window(stk,1983,ymax))),"\n")
    cat("Bpa (abs,rel) = ",Bpa5,Bpa5/mean(ssb(window(stk,1983,ymax))),"\n")
    cat("Flim (abs,rel) = ",Flim5.2,Flim5.2/mean(fbar(window(stk,1983,ymax))),"\n")
    cat("Fpa (abs,rel) = ",Fpa5.2,Fpa5.2/mean(fbar(window(stk,1983,ymax))),"\n")
    cat("MSYBtrigger (abs,rel) = ",MSYBtrigger5,MSYBtrigger5/mean(ssb(window(stk,1983,ymax))),"\n")
    cat("FMSY (abs,rel) = ",Fmsy5.2,Fmsy5.2/mean(fbar(window(stk,1983,ymax))),"\n")
    cat("FP05 (abs,rel) = ",Fp055.2,Fp055.2/mean(fbar(window(stk,1983,ymax))),"\n")
    cat("FMSY (abs,rel) = ",Fmsy5.2_final,Fmsy5.2_final/mean(fbar(window(stk,1983,ymax))),"\n")
    
    if(a=="WK17"){    
      #2018 F/catch advice
      cat("2018 adviceF calculation\n")
      cat("Fmsy = ",Fmsy5.2_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger5, "\n")
      AdviceF_2018 <- Fmsy5.2_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger5)
      #linear interpolate to get the associated catch
      AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
      
      #2019 F/catch advice
      #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
      ImYC <- 1000*round(AdviceC_2018/1000)
      #2019 SSB
      dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
      
      cat("2019 adviceF calculation\n")
      cat("Fmsy = ",Fmsy5.2_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger5, "\n")
      
      AdviceF_2019 <- Fmsy5.2_final
      if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger5){
        AdviceF_2019 <- Fmsy5.2_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger5)
      }
      
      #linear interpolate to get the associated catch
      AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                             y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                             xout=AdviceF_2019,method="linear")$y
      cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
      cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
    }
    
    if (a=="WG17") {
      
      cat("2018 adviceF calculation\n")
      cat("Fmsy = ",Fmsy5.2_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger5, "\n")
      
      #2018 F/catch advice
      AdviceF_2018_5.2_WG17 <- Fmsy5.2_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger5)
      #linear interpolate to get the associated catch
      AdviceC_2018_5.2_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_5.2_WG17,method="linear")$y
      cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_5.2_WG17,AdviceC_2018_5.2_WG17,"\n")
    }
    
    if (a=="WG18") {
      #2019 F/catch advice
      #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
      ImYC <- 1000*round(AdviceC_2018_5.2_WG17/1000)
      #2019 SSB
      dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
      
      cat("2019 adviceF calculation\n")
      cat("Fmsy = ",Fmsy5.2_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger5, "\n")
      
      AdviceF_2019_5.2_WG18 <- Fmsy5.2_final
      if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger5){
        AdviceF_2019_5.2_WG18 <- Fmsy5.2_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger5)
      }
      
      #linear interpolate to get the associated catch
      AdviceC_2019_5.2_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                      y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                      xout=AdviceF_2019_5.2_WG18,method="linear")$y
      
      cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_5.2_WG18,AdviceC_2019_5.2_WG18,"\n")
    }
    
    cat("**********************************************\n")
    
    dfResults <- dplyr::bind_rows(dfResults,
                                  data.frame("Assessment" = rep(a,nRPs),
                                             "Scenario" = rep("5.2",nRPs),
                                             "SRR" = rep(SRR.desc,nRPs),
                                             "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                             "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                             "Fcv" = rep(Fcv,nRPs),
                                             "Fphi" = rep(Fphi,nRPs),
                                             "RP" = RPs,
                                             "Type" = rep("Abs",nRPs),
                                             "Val" = c(Blim5,Bpa5,Flim5.2,Fpa5.2,MSYBtrigger5,Fmsy5.2,Fp055.2,Fmsy5.2_final),
                                             "RunTime" = rep(as.character(rt),nRPs),
                                             stringsAsFactors = FALSE))
    
    dfResults <- dplyr::bind_rows(dfResults,
                                  data.frame("Assessment" = rep(a,nRPs),
                                             "Scenario" = rep("5.2",nRPs),
                                             "SRR" = rep(SRR.desc,nRPs),
                                             "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                             "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                             "Fcv" = rep(Fcv,nRPs),
                                             "Fphi" = rep(Fphi,nRPs),
                                             "RP" = RPs,
                                             "Type" = rep("Rel",nRPs),
                                             "Val" = c(Blim5/mean(ssb(window(stk,1983,ymax))),
                                                       Bpa5/mean(ssb(window(stk,1983,ymax))),
                                                       Flim5.2/mean(fbar(window(stk,1983,ymax))),
                                                       Fpa5.2/mean(fbar(window(stk,1983,ymax))),
                                                       MSYBtrigger5/mean(ssb(window(stk,1983,ymax))),
                                                       Fmsy5.2/mean(fbar(window(stk,1983,ymax))),
                                                       Fp055.2/mean(fbar(window(stk,1983,ymax))),
                                                       Fmsy5.2_final/mean(fbar(window(stk,1983,ymax)))),
                                             "RunTime" = rep(as.character(rt),nRPs),
                                             stringsAsFactors = FALSE))
    
  }
  
  if (3 %in% subcases) {
    
    ###############################################CASE 5.3###########################################
    #Bpa = SSB for 2003, 1995 on
    
    cat("Case 5.3\n")
    cat("Fitting SRR, 1995 on...\n")
    
    SRR5.3 <- eqsr_fit(window(stk,1995,ymax),
                       remove.years = c(),
                       nsamp=1000, models = SRR.models)
    
    #cat("Breakpoint,Mean=",Blim5,SRR5.3$sr.det$a*Blim5,"\n")
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen5_3_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsr_plot(SRR5.3)
    dev.off()  
    
    #Flim
    cat("Flim calculation...\n")
    SIM.Flim5.3 <- eqsim_run(SRR5.3,
                             bio.years = bio.years,
                             bio.const = FALSE,
                             sel.years = sel.years,
                             sel.const = FALSE,
                             Fscan = Fscan,
                             Fcv = 0, Fphi = 0,
                             rhologRec = rhologRec,
                             Blim = Blim5, Bpa = Bpa5,
                             recruitment.trim = recruitment.trim,
                             Btrigger = 0,
                             verbose = FALSE)
    
    Flim5.3 <- SIM.Flim5.3$Refs2["catF","F50"]
    Fpa5.3 <- Flim5.3/1.4
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen5_3_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot(SIM.Flim5.3,catch=TRUE)
    dev.off()
    
    #step 1 - no Btrigger but including stochasticity in 
    #population and fishery and assessment/advice error
    #Fmsy - 1st run, include assessment error
    #for initial FMSY candidate
    cat("Initial Fmsy calculation...\n")
    SIM1.Fmsy5.3 <- eqsim_run(SRR5.3,
                              bio.years = bio.years,
                              bio.const = FALSE,
                              sel.years = sel.years,
                              sel.const = FALSE,
                              Fscan = Fscan,
                              Fcv = Fcv, Fphi = Fphi,
                              rhologRec = rhologRec,
                              Blim = Blim5, Bpa = Bpa5,
                              recruitment.trim = recruitment.trim,
                              Btrigger = 0,
                              verbose = FALSE)
    
    Fmsy5.3 <- SIM1.Fmsy5.3$Refs2["lanF","medianMSY"]
    
    if (Fmsy5.3>Fpa5.3) {Fmsy5.3_final <- Fpa5.3} else {Fmsy5.3_final <- Fmsy5.3}
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen5_3_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot_range(SIM1.Fmsy5.3, type="median")
    dev.off()
    
    #Step 2 select MSYBtrigger
    #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
    MSYBtrigger5 <- Bpa5
    
    #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
    cat("Precautionary check FP05...\n")
    SIM2.Fmsy5.3 <- eqsim_run(SRR5.3,
                              bio.years = bio.years,
                              bio.const = FALSE,
                              sel.years = sel.years,
                              sel.const = FALSE,
                              Fscan = Fscan,
                              Fcv = Fcv, Fphi = Fphi,
                              rhologRec = rhologRec,
                              Blim = Blim5, Bpa = Bpa5,
                              recruitment.trim = recruitment.trim,
                              Btrigger = MSYBtrigger5,
                              verbose = FALSE)
    
    eqsim_plot(SIM2.Fmsy5.3,catch=TRUE)
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen5_3_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot(SIM2.Fmsy5.3,catch=TRUE)
    dev.off()
    
    eqsim_plot_range(SIM2.Fmsy5.3, type="median")
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen5_3_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot_range(SIM2.Fmsy5.3, type="median")
    dev.off()
    
    Fp055.3 <- SIM2.Fmsy5.3$Refs2["catF","F05"]
    #Final Fmsy
    if (Fmsy5.3_final>Fp055.3) {Fmsy5.3_final <- Fp055.3}
    
    #summary of results
    cat("**********************************************\n")
    cat("Case 5.3, 1995 on\n")
    cat("Blim (abs,rel) = ",Blim5,Blim5/mean(ssb(window(stk,1995,ymax))),"\n")
    cat("Bpa (abs,rel) = ",Bpa5,Bpa5/mean(ssb(window(stk,1995,ymax))),"\n")
    cat("Flim (abs,rel) = ",Flim5.3,Flim5.3/mean(fbar(window(stk,1995,ymax))),"\n")
    cat("Fpa (abs,rel) = ",Fpa5.3,Fpa5.3/mean(fbar(window(stk,1995,ymax))),"\n")
    cat("MSYBtrigger (abs,rel) = ",MSYBtrigger5,MSYBtrigger5/mean(ssb(window(stk,1995,ymax))),"\n")
    cat("FMSY (abs,rel) = ",Fmsy5.3,Fmsy5.3/mean(fbar(window(stk,1995,ymax))),"\n")
    cat("FP05 (abs,rel) = ",Fp055.3,Fp055.3/mean(fbar(window(stk,1995,ymax))),"\n")
    cat("FMSY (abs,rel) = ",Fmsy5.3_final,Fmsy5.3_final/mean(fbar(window(stk,1995,ymax))),"\n")
    
    if(a=="WK17"){    
      #2018 F/catch advice
      cat("2018 adviceF calculation\n")
      cat("Fmsy = ",Fmsy5.3_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger5, "\n")
      AdviceF_2018 <- Fmsy5.3_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger5)
      #linear interpolate to get the associated catch
      AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
      
      #2019 F/catch advice
      #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
      ImYC <- 1000*round(AdviceC_2018/1000)
      #2019 SSB
      dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
      
      cat("2019 adviceF calculation\n")
      cat("Fmsy = ",Fmsy5.3_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger5, "\n")
      
      AdviceF_2019 <- Fmsy5.3_final
      if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger5){
        AdviceF_2019 <- Fmsy5.3_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger5)
      }
      
      #linear interpolate to get the associated catch
      AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                             y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                             xout=AdviceF_2019,method="linear")$y
      cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
      cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
    }
    
    if (a=="WG17") {
      
      cat("2018 adviceF calculation\n")
      cat("Fmsy = ",Fmsy5.3_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger5, "\n")
      
      #2018 F/catch advice
      AdviceF_2018_5.3_WG17 <- Fmsy5.3_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger5)
      #linear interpolate to get the associated catch
      AdviceC_2018_5.3_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_5.3_WG17,method="linear")$y
      cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_5.3_WG17,AdviceC_2018_5.3_WG17,"\n")
    }
    
    if (a=="WG18") {
      #2019 F/catch advice
      #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
      ImYC <- 1000*round(AdviceC_2018_5.3_WG17/1000)
      #2019 SSB
      dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
      
      cat("2019 adviceF calculation\n")
      cat("Fmsy = ",Fmsy5.3_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger5, "\n")
      
      AdviceF_2019_5.3_WG18 <- Fmsy5.3_final
      if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger5){
        AdviceF_2019_5.3_WG18 <- Fmsy5.3_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger5)
      }
      
      #linear interpolate to get the associated catch
      AdviceC_2019_5.3_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                      y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                      xout=AdviceF_2019_5.3_WG18,method="linear")$y
      
      cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_5.3_WG18,AdviceC_2019_5.3_WG18,"\n")
    }
    
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
                                             "Val" = c(Blim5,Bpa5,Flim5.3,Fpa5.3,MSYBtrigger5,Fmsy5.3,Fp055.3,Fmsy5.3_final),
                                             "RunTime" = rep(as.character(rt),nRPs),
                                             stringsAsFactors = FALSE))
    
    dfResults <- dplyr::bind_rows(dfResults,
                                  data.frame("Assessment" = rep(a,nRPs),
                                             "Scenario" = rep("5.3",nRPs),
                                             "SRR" = rep(SRR.desc,nRPs),
                                             "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                             "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                             "Fcv" = rep(Fcv,nRPs),
                                             "Fphi" = rep(Fphi,nRPs),
                                             "RP" = RPs,
                                             "Type" = rep("Rel",nRPs),
                                             "Val" = c(Blim5/mean(ssb(window(stk,1995,ymax))),
                                                       Bpa5/mean(ssb(window(stk,1995,ymax))),
                                                       Flim5.3/mean(fbar(window(stk,1995,ymax))),
                                                       Fpa5.3/mean(fbar(window(stk,1995,ymax))),
                                                       MSYBtrigger5/mean(ssb(window(stk,1995,ymax))),
                                                       Fmsy5.3/mean(fbar(window(stk,1995,ymax))),
                                                       Fp055.3/mean(fbar(window(stk,1995,ymax))),
                                                       Fmsy5.3_final/mean(fbar(window(stk,1995,ymax)))),
                                             "RunTime" = rep(as.character(rt),nRPs),
                                             stringsAsFactors = FALSE))
    
  }
  
}


if (6 %in% cases) {
  
  ###############################################CASE 6###########################################
  #Blim = SSB for 2003
  #including 1982
  
  cat("Case.1\n")
  
  Blim6 <- dfWHM$SSB[dfWHM$Year==2003 & dfWHM$Assessment==stk.name]
  Bpa6 <- 1.4*Blim6
  SegregBlim  <- function(ab, ssb) log(ifelse(ssb >= Blim6, ab$a * Blim6, ab$a * ssb))
  
  if (1 %in% subcases) {
    
    #using all SR data 
    cat("Fitting SRR, all data...\n")
    SRR6.1 <- eqsr_fit(stk,
                       remove.years = c(),
                       nsamp=1000, models = SRR.models)
    
    cat("Breakpoint,Mean=",Blim6,SRR6.1$sr.det$a*Blim6,"\n")
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen6_1_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsr_plot(SRR6.1)
    dev.off()  
    
    #Flim
    cat("Flim calculation...\n")
    SIM.Flim6.1 <- eqsim_run(SRR6.1,
                             bio.years = bio.years,
                             bio.const = FALSE,
                             sel.years = sel.years,
                             sel.const = FALSE,
                             Fscan = Fscan,
                             Fcv = 0, Fphi = 0,
                             rhologRec = rhologRec,
                             Blim = Blim6, Bpa = Bpa6,
                             recruitment.trim = recruitment.trim,
                             Btrigger = 0,
                             verbose = FALSE)
    
    eqsim_plot(SIM.Flim6.1,catch=TRUE)
    Flim6.1 <- SIM.Flim6.1$Refs2["catF","F50"]
    Fpa6.1 <- Flim6.1/1.4
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen6_1_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot(SIM.Flim6.1,catch=TRUE)
    dev.off()
    
    #step 1 - no Btrigger but including stochasticity in 
    #population and fishery and assessment/advice error
    #Fmsy - 1st run, include assessment error
    #for initial FMSY candidate
    cat("Initial Fmsy calculation...\n")
    SIM1.Fmsy6.1 <- eqsim_run(SRR6.1,
                              bio.years = bio.years,
                              bio.const = FALSE,
                              sel.years = sel.years,
                              sel.const = FALSE,
                              Fscan = Fscan,
                              Fcv = Fcv, Fphi = Fphi,
                              rhologRec = rhologRec,
                              Blim = Blim6, Bpa = Bpa6,
                              recruitment.trim = recruitment.trim,
                              Btrigger = 0,
                              verbose = FALSE)
    
    Fmsy6.1 <- SIM1.Fmsy6.1$Refs2["lanF","medianMSY"]
    
    if (Fmsy6.1>Fpa6.1) {Fmsy6.1_final <- Fpa6.1} else {Fmsy6.1_final <- Fmsy6.1}
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen6_1_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot_range(SIM1.Fmsy6.1, type="median")
    dev.off()
    
    #Step 2 select MSYBtrigger
    #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
    MSYBtrigger6 <- Bpa6
    
    #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
    cat("Precautionary check FP05...\n")
    SIM2.Fmsy6.1 <- eqsim_run(SRR6.1,
                              bio.years = bio.years,
                              bio.const = FALSE,
                              sel.years = sel.years,
                              sel.const = FALSE,
                              Fscan = Fscan,
                              Fcv = Fcv, Fphi = Fphi,
                              rhologRec = rhologRec,
                              Blim = Blim6, Bpa = Bpa6,
                              recruitment.trim = recruitment.trim,
                              Btrigger = MSYBtrigger6,
                              verbose = FALSE)
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen6_1_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot(SIM2.Fmsy6.1,catch=TRUE)
    dev.off()
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen6_1_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot_range(SIM2.Fmsy6.1, type="median")
    dev.off()
    
    Fp056.1 <- SIM2.Fmsy6.1$Refs2["catF","F05"]
    
    #Final Fmsy
    if (Fmsy6.1_final>Fp056.1) {Fmsy6.1_final <- Fp056.1}
    
    #summary of results
    cat("**********************************************\n")
    cat("Case 6.1, including 1982\n")
    cat("Blim (abs,rel) = ",Blim6,Blim6/mean(ssb(stk)),"\n")
    cat("Bpa (abs,rel) = ",Bpa6,Bpa6/mean(ssb(stk)),"\n")
    cat("Flim (abs,rel) = ",Flim6.1,Flim6.1/mean(fbar(stk)),"\n")
    cat("Fpa (abs,rel) = ",Fpa6.1,Fpa6.1/mean(fbar(stk)),"\n")
    cat("MSYBtrigger (abs,rel) = ",MSYBtrigger6,MSYBtrigger6/mean(ssb(stk)),"\n")
    cat("FMSY (abs,rel) = ",Fmsy6.1,Fmsy6.1/mean(fbar(stk)),"\n")
    cat("FP05 (abs,rel) = ",Fp056.1,Fp056.1/mean(fbar(stk)),"\n")
    cat("FMSY (abs,rel) = ",Fmsy6.1_final,Fmsy6.1_final/mean(fbar(stk)),"\n")
    
    if(a=="WK17"){    
      #2018 F/catch advice
      cat("2018 adviceF calculation\n")
      cat("Fmsy = ",Fmsy6.1_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger6, "\n")
      AdviceF_2018 <- Fmsy6.1_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger6)
      #linear interpolate to get the associated catch
      AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
      
      #2019 F/catch advice
      #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
      ImYC <- 1000*round(AdviceC_2018/1000)
      #2019 SSB
      dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
      
      cat("2019 adviceF calculation\n")
      cat("Fmsy = ",Fmsy6.1_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger6, "\n")
      
      AdviceF_2019 <- Fmsy6.1_final
      if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger6){
        AdviceF_2019 <- Fmsy6.1_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger6)
      }
      
      #linear interpolate to get the associated catch
      AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                             y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                             xout=AdviceF_2019,method="linear")$y
      cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
      cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
    }
    
    if (a=="WG17") {
      
      cat("2018 adviceF calculation\n")
      cat("Fmsy = ",Fmsy6.1_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger6, "\n")
      
      #2018 F/catch advice
      AdviceF_2018_6.1_WG17 <- Fmsy6.1_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger6)
      #linear interpolate to get the associated catch
      AdviceC_2018_6.1_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_6.1_WG17,method="linear")$y
      cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_6.1_WG17,AdviceC_2018_6.1_WG17,"\n")
    }
    
    if (a=="WG18") {
      #2019 F/catch advice
      #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
      ImYC <- 1000*round(AdviceC_2018_6.1_WG17/1000)
      #2019 SSB
      dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
      
      cat("2019 adviceF calculation\n")
      cat("Fmsy = ",Fmsy6.1_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger6, "\n")
      
      AdviceF_2019_6.1_WG18 <- Fmsy6.1_final
      if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger6){
        AdviceF_2019_6.1_WG18 <- Fmsy6.1_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger6)
      }
      
      #linear interpolate to get the associated catch
      AdviceC_2019_6.1_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                      y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                      xout=AdviceF_2019_6.1_WG18,method="linear")$y
      
      cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_6.1_WG18,AdviceC_2019_6.1_WG18,"\n")
    }
    
    cat("**********************************************\n")
    
    
    dfResults <- dplyr::bind_rows(dfResults,
                                  data.frame("Assessment" = rep(a,nRPs),
                                             "Scenario" = rep("6.1",nRPs),
                                             "SRR" = rep(SRR.desc,nRPs),
                                             "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                             "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                             "Fcv" = rep(Fcv,nRPs),
                                             "Fphi" = rep(Fphi,nRPs),
                                             "RP" = RPs,
                                             "Type" = rep("Abs",nRPs),
                                             "Val" = c(Blim6,Bpa6,Flim6.1,Fpa6.1,MSYBtrigger6,Fmsy6.1,Fp056.1,Fmsy6.1_final),
                                             "RunTime" = rep(as.character(rt),nRPs),
                                             stringsAsFactors = FALSE))
    
    dfResults <- dplyr::bind_rows(dfResults,
                                  data.frame("Assessment" = rep(a,nRPs),
                                             "Scenario" = rep("6.1",nRPs),
                                             "SRR" = rep(SRR.desc,nRPs),
                                             "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                             "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                             "Fcv" = rep(Fcv,nRPs),
                                             "Fphi" = rep(Fphi,nRPs),
                                             "RP" = RPs,
                                             "Type" = rep("Rel",nRPs),
                                             "Val" = c(Blim6/mean(ssb(stk)),
                                                       Bpa6/mean(ssb(stk)),
                                                       Flim6.1/mean(fbar(stk)),
                                                       Fpa6.1/mean(fbar(stk)),
                                                       MSYBtrigger6/mean(ssb(stk)),
                                                       Fmsy6.1/mean(fbar(stk)),
                                                       Fp056.1/mean(fbar(stk)),
                                                       Fmsy6.1_final/mean(fbar(stk))),
                                             "RunTime" = rep(as.character(rt),nRPs),
                                             stringsAsFactors = FALSE))
  }
  
  if (2 %in% subcases) {
    
    ###############################################CASE 6.2###########################################
    #Blim = SSB for 2003, excluding 1982
    
    cat("Case 6.2\n")
    cat("Fitting SRR, excluding 1982...\n")
    SRR6.2 <- eqsr_fit(stk,
                       remove.years = c(1982),
                       nsamp=1000, models = SRR.models)
    
    #cat("Breakpoint,Mean=",Blim6,SRR6.2$sr.det$a*Blim6,"\n")
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen6_2_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsr_plot(SRR6.2)
    dev.off()  
    
    #Flim
    cat("Flim calculation...\n")
    SIM.Flim6.2 <- eqsim_run(SRR6.2,
                             bio.years = bio.years,
                             bio.const = FALSE,
                             sel.years = sel.years,
                             sel.const = FALSE,
                             Fscan = Fscan,
                             Fcv = 0, Fphi = 0,
                             rhologRec = rhologRec,
                             Blim = Blim6, Bpa = Bpa6,
                             recruitment.trim = recruitment.trim,
                             Btrigger = 0,
                             verbose = FALSE)
    
    Flim6.2 <- SIM.Flim6.2$Refs2["catF","F50"]
    Fpa6.2 <- Flim6.2/1.4
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen6_2_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot(SIM.Flim6.2,catch=TRUE)
    dev.off()
    
    #step 1 - no Btrigger but including stochasticity in 
    #population and fishery and assessment/advice error
    #Fmsy - 1st run, include assessment error
    #for initial FMSY candidate
    cat("Initial Fmsy calculation...\n")
    SIM1.Fmsy6.2 <- eqsim_run(SRR6.2,
                              bio.years = bio.years,
                              bio.const = FALSE,
                              sel.years = sel.years,
                              sel.const = FALSE,
                              Fscan = Fscan,
                              Fcv = Fcv, Fphi = Fphi,
                              rhologRec = rhologRec,
                              Blim = Blim6, Bpa = Bpa6,
                              recruitment.trim = recruitment.trim,
                              Btrigger = 0,
                              verbose = FALSE)
    
    Fmsy6.2 <- SIM1.Fmsy6.2$Refs2["lanF","medianMSY"]
    
    if (Fmsy6.2>Fpa6.2) {Fmsy6.2_final <- Fpa6.2} else {Fmsy6.2_final <- Fmsy6.2}
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen6_2_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot_range(SIM1.Fmsy6.2, type="median")
    dev.off()
    
    #Step 2 select MSYBtrigger
    #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
    MSYBtrigger6 <- Bpa6
    
    #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
    cat("Precautionary check FP05...\n")
    SIM2.Fmsy6.2 <- eqsim_run(SRR6.2,
                              bio.years = bio.years,
                              bio.const = FALSE,
                              sel.years = sel.years,
                              sel.const = FALSE,
                              Fscan = Fscan,
                              Fcv = Fcv, Fphi = Fphi,
                              rhologRec = rhologRec,
                              Blim = Blim6, Bpa = Bpa6,
                              recruitment.trim = recruitment.trim,
                              Btrigger = MSYBtrigger6,
                              verbose = FALSE)
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen6_2_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot(SIM2.Fmsy6.2,catch=TRUE)
    dev.off()
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen6_2_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot_range(SIM2.Fmsy6.2, type="median")
    dev.off()
    
    Fp056.2 <- SIM2.Fmsy6.2$Refs2["catF","F05"]
    
    #Final Fmsy
    if (Fmsy6.2_final>Fp056.2) {Fmsy6.2_final <- Fp056.2}
    
    #summary of results
    cat("**********************************************\n")
    cat("Case 6.2, excluding 1982\n")
    cat("Blim (abs,rel) = ",Blim6,Blim6/mean(ssb(window(stk,1983,ymax))),"\n")
    cat("Bpa (abs,rel) = ",Bpa6,Bpa6/mean(ssb(window(stk,1983,ymax))),"\n")
    cat("Flim (abs,rel) = ",Flim6.2,Flim6.2/mean(fbar(window(stk,1983,ymax))),"\n")
    cat("Fpa (abs,rel) = ",Fpa6.2,Fpa6.2/mean(fbar(window(stk,1983,ymax))),"\n")
    cat("MSYBtrigger (abs,rel) = ",MSYBtrigger6,MSYBtrigger6/mean(ssb(window(stk,1983,ymax))),"\n")
    cat("FMSY (abs,rel) = ",Fmsy6.2,Fmsy6.2/mean(fbar(window(stk,1983,ymax))),"\n")
    cat("FP05 (abs,rel) = ",Fp056.2,Fp056.2/mean(fbar(window(stk,1983,ymax))),"\n")
    cat("FMSY (abs,rel) = ",Fmsy6.2_final,Fmsy6.2_final/mean(fbar(window(stk,1983,ymax))),"\n")
    
    if(a=="WK17"){    
      #2018 F/catch advice
      cat("2018 adviceF calculation\n")
      cat("Fmsy = ",Fmsy6.2_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger6, "\n")
      AdviceF_2018 <- Fmsy6.2_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger6)
      #linear interpolate to get the associated catch
      AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
      
      #2019 F/catch advice
      #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
      ImYC <- 1000*round(AdviceC_2018/1000)
      #2019 SSB
      dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
      
      cat("2019 adviceF calculation\n")
      cat("Fmsy = ",Fmsy6.2_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger6, "\n")
      
      AdviceF_2019 <- Fmsy6.2_final
      if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger6){
        AdviceF_2019 <- Fmsy6.2_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger6)
      }
      
      #linear interpolate to get the associated catch
      AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                             y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                             xout=AdviceF_2019,method="linear")$y
      cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
      cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
    }
    
    if (a=="WG17") {
      
      cat("2018 adviceF calculation\n")
      cat("Fmsy = ",Fmsy6.2_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger6, "\n")
      
      #2018 F/catch advice
      AdviceF_2018_6.2_WG17 <- Fmsy6.2_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger6)
      #linear interpolate to get the associated catch
      AdviceC_2018_6.2_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_6.2_WG17,method="linear")$y
      cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_6.2_WG17,AdviceC_2018_6.2_WG17,"\n")
    }
    
    if (a=="WG18") {
      #2019 F/catch advice
      #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
      ImYC <- 1000*round(AdviceC_2018_6.2_WG17/1000)
      #2019 SSB
      dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
      
      cat("2019 adviceF calculation\n")
      cat("Fmsy = ",Fmsy6.2_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger6, "\n")
      
      AdviceF_2019_6.2_WG18 <- Fmsy6.2_final
      if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger6){
        AdviceF_2019_6.2_WG18 <- Fmsy6.2_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger6)
      }
      
      #linear interpolate to get the associated catch
      AdviceC_2019_6.2_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                      y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                      xout=AdviceF_2019_6.2_WG18,method="linear")$y
      
      cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_6.2_WG18,AdviceC_2019_6.2_WG18,"\n")
    }
    
    cat("**********************************************\n")
    
    dfResults <- dplyr::bind_rows(dfResults,
                                  data.frame("Assessment" = rep(a,nRPs),
                                             "Scenario" = rep("6.2",nRPs),
                                             "SRR" = rep(SRR.desc,nRPs),
                                             "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                             "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                             "Fcv" = rep(Fcv,nRPs),
                                             "Fphi" = rep(Fphi,nRPs),
                                             "RP" = RPs,
                                             "Type" = rep("Abs",nRPs),
                                             "Val" = c(Blim6,Bpa6,Flim6.2,Fpa6.2,MSYBtrigger6,Fmsy6.2,Fp056.2,Fmsy6.2_final),
                                             "RunTime" = rep(as.character(rt),nRPs),
                                             stringsAsFactors = FALSE))
    
    dfResults <- dplyr::bind_rows(dfResults,
                                  data.frame("Assessment" = rep(a,nRPs),
                                             "Scenario" = rep("6.2",nRPs),
                                             "SRR" = rep(SRR.desc,nRPs),
                                             "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                             "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                             "Fcv" = rep(Fcv,nRPs),
                                             "Fphi" = rep(Fphi,nRPs),
                                             "RP" = RPs,
                                             "Type" = rep("Rel",nRPs),
                                             "Val" = c(Blim6/mean(ssb(window(stk,1983,ymax))),
                                                       Bpa6/mean(ssb(window(stk,1983,ymax))),
                                                       Flim6.2/mean(fbar(window(stk,1983,ymax))),
                                                       Fpa6.2/mean(fbar(window(stk,1983,ymax))),
                                                       MSYBtrigger6/mean(ssb(window(stk,1983,ymax))),
                                                       Fmsy6.2/mean(fbar(window(stk,1983,ymax))),
                                                       Fp056.2/mean(fbar(window(stk,1983,ymax))),
                                                       Fmsy6.2_final/mean(fbar(window(stk,1983,ymax)))),
                                             "RunTime" = rep(as.character(rt),nRPs),
                                             stringsAsFactors = FALSE))
    
  }
  
  if (3 %in% subcases) {
    
    ###############################################CASE 6.3###########################################
    #Blim = SSB for 2003, 1995 on
    
    cat("Case 6.3\n")
    cat("Fitting SRR, 1995 on...\n")
    
    SRR6.3 <- eqsr_fit(window(stk,1995,ymax),
                       remove.years = c(),
                       nsamp=1000, models = SRR.models)
    
    #cat("Breakpoint,Mean=",Blim5,SRR5.3$sr.det$a*Blim5,"\n")
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen6_3_SRR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsr_plot(SRR6.3)
    dev.off()  
    
    #Flim
    cat("Flim calculation...\n")
    SIM.Flim6.3 <- eqsim_run(SRR6.3,
                             bio.years = bio.years,
                             bio.const = FALSE,
                             sel.years = sel.years,
                             sel.const = FALSE,
                             Fscan = Fscan,
                             Fcv = 0, Fphi = 0,
                             rhologRec = rhologRec,
                             Blim = Blim6, Bpa = Bpa6,
                             recruitment.trim = recruitment.trim,
                             Btrigger = 0,
                             verbose = FALSE)
    
    Flim6.3 <- SIM.Flim6.3$Refs2["catF","F50"]
    Fpa6.3 <- Flim6.3/1.4
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen6_3_Flim",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot(SIM.Flim6.3,catch=TRUE)
    dev.off()
    
    #step 1 - no Btrigger but including stochasticity in 
    #population and fishery and assessment/advice error
    #Fmsy - 1st run, include assessment error
    #for initial FMSY candidate
    cat("Initial Fmsy calculation...\n")
    SIM1.Fmsy6.3 <- eqsim_run(SRR6.3,
                              bio.years = bio.years,
                              bio.const = FALSE,
                              sel.years = sel.years,
                              sel.const = FALSE,
                              Fscan = Fscan,
                              Fcv = Fcv, Fphi = Fphi,
                              rhologRec = rhologRec,
                              Blim = Blim6, Bpa = Bpa6,
                              recruitment.trim = recruitment.trim,
                              Btrigger = 0,
                              verbose = FALSE)
    
    Fmsy6.3 <- SIM1.Fmsy6.3$Refs2["lanF","medianMSY"]
    
    if (Fmsy6.3>Fpa6.3) {Fmsy6.3_final <- Fpa6.3} else {Fmsy6.3_final <- Fmsy6.3}
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen6_3_Fmsy_Range_noAR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot_range(SIM1.Fmsy6.3, type="median")
    dev.off()
    
    #Step 2 select MSYBtrigger
    #Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
    MSYBtrigger6 <- Bpa6
    
    #Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
    cat("Precautionary check FP05...\n")
    SIM2.Fmsy6.3 <- eqsim_run(SRR6.3,
                              bio.years = bio.years,
                              bio.const = FALSE,
                              sel.years = sel.years,
                              sel.const = FALSE,
                              Fscan = Fscan,
                              Fcv = Fcv, Fphi = Fphi,
                              rhologRec = rhologRec,
                              Blim = Blim6, Bpa = Bpa6,
                              recruitment.trim = recruitment.trim,
                              Btrigger = MSYBtrigger6,
                              verbose = FALSE)
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen6_3_Fmsy_Summary_AR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot(SIM2.Fmsy6.3,catch=TRUE)
    dev.off()
    
    Cairo(file = file.path(graphics.dir,runName,paste0(stk.name,"_Scen6_3_Fmsy_Range_AR",".",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    eqsim_plot_range(SIM2.Fmsy6.3, type="median")
    dev.off()
    
    Fp056.3 <- SIM2.Fmsy6.3$Refs2["catF","F05"]
    #Final Fmsy
    if (Fmsy6.3_final>Fp056.3) {Fmsy6.3_final <- Fp056.3}
    
    #summary of results
    cat("**********************************************\n")
    cat("Case 6.3, 1995 on\n")
    cat("Blim (abs,rel) = ",Blim6,Blim6/mean(ssb(window(stk,1995,ymax))),"\n")
    cat("Bpa (abs,rel) = ",Bpa6,Bpa6/mean(ssb(window(stk,1995,ymax))),"\n")
    cat("Flim (abs,rel) = ",Flim6.3,Flim6.3/mean(fbar(window(stk,1995,ymax))),"\n")
    cat("Fpa (abs,rel) = ",Fpa6.3,Fpa6.3/mean(fbar(window(stk,1995,ymax))),"\n")
    cat("MSYBtrigger (abs,rel) = ",MSYBtrigger6,MSYBtrigger6/mean(ssb(window(stk,1995,ymax))),"\n")
    cat("FMSY (abs,rel) = ",Fmsy6.3,Fmsy6.3/mean(fbar(window(stk,1995,ymax))),"\n")
    cat("FP05 (abs,rel) = ",Fp056.3,Fp056.3/mean(fbar(window(stk,1995,ymax))),"\n")
    cat("FMSY (abs,rel) = ",Fmsy6.3_final,Fmsy6.3_final/mean(fbar(window(stk,1995,ymax))),"\n")
    
    if(a=="WK17"){    
      #2018 F/catch advice
      cat("2018 adviceF calculation\n")
      cat("Fmsy = ",Fmsy6.3_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger6, "\n")
      AdviceF_2018 <- Fmsy6.3_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger6)
      #linear interpolate to get the associated catch
      AdviceC_2018 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018,method="linear")$y
      
      #2019 F/catch advice
      #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
      ImYC <- 1000*round(AdviceC_2018/1000)
      #2019 SSB
      dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
      
      cat("2019 adviceF calculation\n")
      cat("Fmsy = ",Fmsy6.3_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger6, "\n")
      
      AdviceF_2019 <- Fmsy6.3_final
      if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger6){
        AdviceF_2019 <- Fmsy6.3_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger6)
      }
      
      #linear interpolate to get the associated catch
      AdviceC_2019 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                             y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                             xout=AdviceF_2019,method="linear")$y
      cat("Revised 2018 advice F/Catch (WK17 settings)= ",AdviceF_2018,AdviceC_2018,"\n")
      cat("Revised 2019 advice F/Catch (WK17 settings)= ",AdviceF_2019,AdviceC_2019,"\n")
    }
    
    if (a=="WG17") {
      
      cat("2018 adviceF calculation\n")
      cat("Fmsy = ",Fmsy6.3_final, "SSB Jan 1 = ",dfSTF2017$SSB_2018[1], "MSYBtrigger = ", MSYBtrigger6, "\n")
      
      #2018 F/catch advice
      AdviceF_2018_6.3_WG17 <- Fmsy6.3_final*(dfSTF2017$SSB_2018[1]/MSYBtrigger6)
      #linear interpolate to get the associated catch
      AdviceC_2018_6.3_WG17 <- approx(x=dfSTF2017$Fbar,y=dfSTF2017$Catch_2018,xout=AdviceF_2018_6.3_WG17,method="linear")$y
      cat("Revised 2018 advice F/Catch (WG17 settings)= ",AdviceF_2018_6.3_WG17,AdviceC_2018_6.3_WG17,"\n")
    }
    
    if (a=="WG18") {
      #2019 F/catch advice
      #IMY catch is assumed to be the AdviceC_2018, rounded to the nearest 1000
      ImYC <- 1000*round(AdviceC_2018_6.3_WG17/1000)
      #2019 SSB
      dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]
      
      cat("2019 adviceF calculation\n")
      cat("Fmsy = ",Fmsy6.3_final, "SSB Jan 1 = ",dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1], "MSYBtrigger = ", MSYBtrigger6, "\n")
      
      AdviceF_2019_6.3_WG18 <- Fmsy6.3_final
      if (dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]<MSYBtrigger6){
        AdviceF_2019_6.3_WG18 <- Fmsy6.3_final*(dfSTF2018$SSB_2019[dfSTF2018$ImY==ImYC][1]/MSYBtrigger6)
      }
      
      #linear interpolate to get the associated catch
      AdviceC_2019_6.3_WG18 <- approx(x=dfSTF2018[dfSTF2018$ImY==ImYC,]$Fbar,
                                      y=dfSTF2018[dfSTF2018$ImY==ImYC,]$Catch_2019,
                                      xout=AdviceF_2019_6.3_WG18,method="linear")$y
      
      cat("Revised 2019 advice F/Catch (WG18 settings)= ",AdviceF_2019_6.3_WG18,AdviceC_2019_6.3_WG18,"\n")
    }
    
    cat("**********************************************\n")
    
    dfResults <- dplyr::bind_rows(dfResults,
                                  data.frame("Assessment" = rep(a,nRPs),
                                             "Scenario" = rep("6.3",nRPs),
                                             "SRR" = rep(SRR.desc,nRPs),
                                             "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                             "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                             "Fcv" = rep(Fcv,nRPs),
                                             "Fphi" = rep(Fphi,nRPs),
                                             "RP" = RPs,
                                             "Type" = rep("Abs",nRPs),
                                             "Val" = c(Blim6,Bpa6,Flim6.3,Fpa6.3,MSYBtrigger6,Fmsy6.3,Fp056.3,Fmsy6.3_final),
                                             "RunTime" = rep(as.character(rt),nRPs),
                                             stringsAsFactors = FALSE))
    
    dfResults <- dplyr::bind_rows(dfResults,
                                  data.frame("Assessment" = rep(a,nRPs),
                                             "Scenario" = rep("6.3",nRPs),
                                             "SRR" = rep(SRR.desc,nRPs),
                                             "BioYrs" = rep(paste(bio.years,collapse=","),nRPs),
                                             "SelYrs" = rep(paste(sel.years,collapse=","),nRPs),
                                             "Fcv" = rep(Fcv,nRPs),
                                             "Fphi" = rep(Fphi,nRPs),
                                             "RP" = RPs,
                                             "Type" = rep("Rel",nRPs),
                                             "Val" = c(Blim6/mean(ssb(window(stk,1995,ymax))),
                                                       Bpa6/mean(ssb(window(stk,1995,ymax))),
                                                       Flim6.3/mean(fbar(window(stk,1995,ymax))),
                                                       Fpa6.3/mean(fbar(window(stk,1995,ymax))),
                                                       MSYBtrigger6/mean(ssb(window(stk,1995,ymax))),
                                                       Fmsy6.3/mean(fbar(window(stk,1995,ymax))),
                                                       Fp056.3/mean(fbar(window(stk,1995,ymax))),
                                                       Fmsy6.3_final/mean(fbar(window(stk,1995,ymax)))),
                                             "RunTime" = rep(as.character(rt),nRPs),
                                             stringsAsFactors = FALSE))
    
  }
  
}

} #loop over assessments

dfResults$runName <- runName

save(dfResults, file = file.path(RData.dir,paste0(runName,"_","cases",cases.desc,"_subcases_",
                            subcases.desc,"_SRR_",SRR.desc,"_",
                            as.character(format(Sys.time(), "%d-%b-%Y %H.%M")),".RData")))

sink()

