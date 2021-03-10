# ================================================================================================================
# EqSim HCR simulator
# 
# 02_EqSim_simulate.R
# 
# The EqSim simulator was developed by David Miller and Colin Millar
# Further enhanced by Andy Campbell and Martin Pastoors and applied to Western Horse mackerel in June 2020
#
# 24/06/2020 generic option; code now independent of fish stock
# 25/06/2020 tested on mackerel stock
# 27/06/2020 tested on 1000 iters of SAM assessment
# 29/06/2020 modified and bugchecked by Andy Campbell
# 01/07/2020 included additional features by Martin Pastoors
# ================================================================================================================

# source(file.path(getwd(),"Scripts","01_EqSim_setup.R"))

#Note: niters and nyr could be included in the OM or MP definitions

#basic simulation settings
#niters <- 10000
#niters <- 1000
niters <- 1000
nyr <- 23

# simulation periods
per1 <- 5
per2 <- 5
# per3 is simply the remainder

# set up the OM =========================================================================================================

#OM <- OM2.1   #WGWIDE 2019, const weights, selection

# WHOM SS 2019
# stock          <- "WHOM"
# assess         <- "SS3"
# assessyear     <- "2019"
# FLStockfile    <- "WGWIDE19.RData"
# FLStockSimfile <- "WHOM_SS19_FLS_V1.RData"    #V1 iterations as single FLStock
# FLStockSimfile <- "WHOM_SS19_FLS_V2.RData"    #V2 new draw, contains variability in selection and weights
# FLStockSimfile <- "WHOM_SS19_FLS_V3.RData"    #V3 improved variability in selection and weights
# OM             <- OM2.2                       #WGWIDE SS 2019, stochastic weights, selection

# WHOM SS 2020
# stock          <- "WHOM"
# assess         <- "SS3"
# assessyear     <- "2020"
# FLStockfile    <- "WGWIDE20.RData"
# FLStockSimfile <- "WHOM_SS20_FLS_V2.RData"
# OM             <- OM2.3                       #WGWIDE SS 2020, stochastic weights, selection

# WHOM SAM 2019
stock          <- "WHOM"
assess         <- "SAM"
assessyear     <- "2019"
FLStockfile    <- "WHOM_SAM19_FLS_WGWIDE.RData"
FLStockSimfile <- "WHOM_SAM19_FLS_converged.RData"
OM             <- OM2.4                         #WGWIDE SAM 2019, stochastic weights, selection

# WHOM SAM 2019
# stock          <- "WHOM"
# assess         <- "SAM"
# assessyear     <- "2020"
# FLStockfile    <- "WHOM_SAM20_FLS_WGWIDE.RData"
# FLStockSimfile <- "WHOM_SAM20_FLS_converged.RData"
# OM             <- OM2.5                       #WGWIDE SAM 2020, stochastic weights, selection

#assessment FLStock
FLS <- loadRData(file.path(RData.dir,FLStockfile)) %>% FLCore::setPlusGroup(., 15)

#Blim <- min(ssb(FLS))
#The IBP in 2019 selected SSB in 2003 as a proxy for Bpa and derived Blim from this (Bpa/1.4)
#given the lack of any clear SRR and sensitivity of the proportions of mixed models
#to individual data points, a segmented regression with breakpoint at Blim is the default SRR

# MP 24/6/2020: this seems wrong; segmented regression should go through Blim, not through Bloss
# Blimloss <- min(OM$refPts$Blim,min(ssb(FLS)))
Blimloss <- OM$refPts$Blim
SegregBlim  <- function(ab, ssb) log(ifelse(ssb >= Blimloss, ab$a * Blimloss, ab$a * ssb))

set.seed(1)

#segmented regression with breakpoint at Blim, from 1995 excluding terminal
if (grepl("WGWIDE19",OM$desc)) {SRR <- eqsr_fit(window(FLS,1995,2018), remove.years = c(2018), nsamp=niters, models = c("SegregBlim"))}
if (grepl("WGWIDE20",OM$desc)) {SRR <- eqsr_fit(window(FLS,1995,2019), remove.years = c(2019), nsamp=niters, models = c("SegregBlim"))}

#emf(file = file.path(Res.dir,paste0(OM$desc,"_SRR.emf")), width = 7, height = 7)
#eqsr_plot(SRR)
#dev.off()

#reassign FLStock object with updated weights into stk slot of SRR 
SRR$stk <- FLS

#start with simulated initial populations
FLSs        <- loadRData(file.path(RData.dir,FLStockSimfile))
FLSs@stock  <- ssb(FLSs)

#add required number of stochastic FLStocks to FIT object
#SRR$stks <- FLSs[(length(FLSs)-niters+1):(length(FLSs))]
SRR$stks <- FLSs

#start,end,vectors of observation and simulation years
#simulation starts in assessment terminal year
minObsYear <- range(SRR$stk)["minyear"]
maxObsYear <- range(SRR$stk)["maxyear"]
obsYears <- ac(seq(minObsYear,maxObsYear))
yStart <- as.numeric(maxObsYear)
yEnd <- yStart + nyr - 1
simYears <- ac(seq(yStart,yEnd))
numWorm <- 5

# assessment dataframe
runName   <- paste(stock,assess,assessyear, OM$code,niters,sep="_")
dfassess  <- fassess_df(runName=runName, FLSs=FLSs, OM = OM, numWorm = numWorm) 

# Loop over management procedures

# mp <- c("MP5.11")
# for (mp in c("MP5.00","MP5.01","MP5.10","MP5.11","MP5.20","MP5.21")) {
# for (mp in c("MP5.00","MP5.23")) {
for (mp in c("MP5.00","MP5.01","MP5.03",
            "MP5.10","MP5.11","MP5.13",
            "MP5.20","MP5.21","MP5.23")) {
  
  MP <- get(mp)
  
  invisible(gc())
  
  runName <- paste(stock,assess,assessyear, OM$code,niters, MP$code,nyr,sep="_")
  

  #exploitation constraints
  if (grepl("WGWIDE19",OM$desc)) {
    #2018 catch known, 2019 as assumed during WGWIDE 2019, 2020 as advised
    dfExplConstraints <- data.frame("Type" = c("Catch","Catch","Catch"), 
                                    "YearNum" = c("1","2","3"),
                                    "Val" = c(101682,110381,83954), 
                                    stringsAsFactors = FALSE)
  }
    
  if (grepl("WGWIDE20",OM$desc)) {
    #2019 catch known, 2020 as assumed during WGWIDE 2020, 2021 as advised
    dfExplConstraints <- data.frame("Type" = c("Catch","Catch","Catch"), 
                                    "YearNum" = c("1","2","3"),
                                    "Val" = c(124947,69527,81376), 
                                    stringsAsFactors = FALSE)
  }
  
  #test for recruitment failure, keep exploration constant at 80kt (regardless of HCR) for first 10 years during which recruitment
  #failure is simulated (1/10 of normal). Then, all HCRs should start from same point
  #dfExplConstraints <- data.frame("Type" = rep("Catch",10), 
  #                                "YearNum" = as.character(seq(1,10)),
  #                                "Val" = c(101682,110381,83954,rep(80000,7)), 
  #                                stringsAsFactors = FALSE)
  
  
  #min/max TAC
  if (!is.na(MP$minTAC)) {
    dfExplConstraints <- dplyr::bind_rows(dfExplConstraints,
                                          data.frame("Type" = "MinTAC",
                                                     "YearNum" = "all",
                                                     "Val" = MP$minTAC,
                                                     stringsAsFactors = FALSE))
  }
  
  #min/max TAC
  if (!is.na(MP$maxTAC)) {
    dfExplConstraints <- dplyr::bind_rows(dfExplConstraints,
                                          data.frame("Type" = "MaxTAC",
                                                     "YearNum" = "all",
                                                     "Val" = MP$maxTAC,
                                                     stringsAsFactors = FALSE))
  }
  
  #IAV
  if (!any(is.na(MP$TAC_IAV))) {
    if (length(MP$TAC_IAV)==2){
      dfExplConstraints <- dplyr::bind_rows(dfExplConstraints,
                                            data.frame("Type" = c("IAVInc","IAVDec"),
                                                       "YearNum" = c("all","all"),
                                                       "Val" = c(MP$TAC_IAV[1],MP$TAC_IAV[2]),
                                                       stringsAsFactors = FALSE))
    } else {
      stop("IAV needs to be vector of length 2 (limit on increase, limit on decrease)")
      #assume same for both
      dfExplConstraints <- dplyr::bind_rows(dfExplConstraints,
                                            data.frame("Type" = c("IAVInc","IAVDec"),
                                                       "YearNum" = c("all","all"),
                                                       "Val" = c(MP$TAC_IAV[1],MP$TAC_IAV[1]),
                                                       stringsAsFactors = FALSE))
    }
  }
  
  #statistical periods for reporting
  lStatPer <- lStatPer2 <- list()
  #create a list for the output statistical periods
  #annual statistics, for each historic and simulated year
  for (y in seq(minObsYear,yEnd)){lStatPer[[ac(y)]]<-c(y,y)}
  for (y in seq(maxObsYear+1,yEnd)){lStatPer2[[ac(y)]]<-c(y,y)}
  
  #Short (first 5 after constraints), Medium (next 5) and Long Term (next 20)
  yConstraints <- 3
  lStatPer[['CU']] <- lStatPer2[['CU']] <- c(yStart, yStart+yConstraints-1)
  lStatPer[['ST']] <- c(yStart+yConstraints,yStart+yConstraints+(per1-1))
  lStatPer2[['ST']] <- c(yStart+yConstraints+1,yStart+yConstraints+(per1-1))
  lStatPer[['MT']] <- lStatPer2[['MT']] <- c(yStart+yConstraints+per1,yStart+yConstraints+(per1+per2-1))
  lStatPer[['LT']] <- lStatPer2[['LT']] <- c(yStart+yConstraints+per1+per2,yStart+nyr-1)
  lStatPer[['HI']] <- c(range(FLSs)[["minyear"]], lStatPer[['CU']][1]-1)
  
  set.seed(1)      # DO YOU NEED TO SET THIS IN THE LOOP?
   
  sim <- eqsim_run(fit = SRR,
                  bio.years = OM$BioYrs,
                  bio.const = OM$BioConst,
                  sel.years = OM$SelYrs,
                  sel.const = OM$SelConst,
                  Fscan = fGetValsScan(MP$F_target,OM$refPts),
                  Fcv = MP$Obs$cvF,
                  Fphi = MP$Obs$phiF,
                  SSBcv = MP$Obs$cvSSB,
                  SSBphi = MP$Obs$phiSSB,
                  Blim = OM$refPts$Blim,
                  Nrun = nyr,
                  calc.RPs = FALSE,
                  dfExplConstraints = dfExplConstraints,
                  Btrigger = fGetValsScan(MP$B_trigger,OM$refPts),
                  HCRName = paste0("fHCR_",MP$HCRName))

  SimRuns <- sim$simStks
  
  #save output
  
  #create a folder for the output and save simRuns data
  dir.create(path = file.path(Res.dir,runName), showWarnings = TRUE, recursive = TRUE)
  save(SimRuns,file = file.path(Res.dir,runName,paste0(runName,"_SimRuns.RData")))
  
  #Write the output to dropbox dir (necessary to save entire image?)
  #save.image(file = file.path(dropbox.dir,paste0(runName,"_Workspace.Rdata")))
  #save.image(file = file.path(Res.dir,runName,paste0(runName,"_Workspace.Rdata")))
  
  #Percentiles to report, number of worm lines for plots
  percentiles = c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975)
  
  #Create list objects to store stocks, stats
  Stocks <- AllStats <- list()
  #create a template by extending the time frame out to the final year of the simulation
  stockTemplate <- window(SRR$stk, end=yEnd)
  
  #set the maturity, natural mortality and f/m proportions
  mat(stockTemplate)[,simYears] <- mat(stockTemplate)[,ac(yStart-1)]
  m(stockTemplate)[,simYears] <- m(stockTemplate)[,ac(yStart-1)]
  m.spwn(stockTemplate)[,simYears] <- m.spwn(stockTemplate)[,ac(yStart-1)]
  harvest.spwn(stockTemplate)[,simYears] <- harvest.spwn(stockTemplate)[,ac(yStart-1)]
  #quick look
  #plot(stockTemplate)
  
  #extend object to store neessary number of iterations
  stockTemplate <- FLCore::propagate(stockTemplate,niters)
  
  #populate each stock object
  #year dimension is trimmed to the actual simulation period (may be longer depending on HCR implemented)
  for (ii in names(SimRuns)) {
    
    cat("Calculating statistics for run with f =",ii,"\n")
    
    Stocks[[ii]] <- stockTemplate
    stock.n(Stocks[[ii]])[,simYears,,,,] <- SimRuns[[ii]]$N[,1:(yEnd-yStart+1),]
    harvest(Stocks[[ii]])[,simYears,,,,] <- SimRuns[[ii]]$F[,1:(yEnd-yStart+1),]
    catch.n(Stocks[[ii]])[,simYears,,,,] <- SimRuns[[ii]]$C[,1:(yEnd-yStart+1),]
    catch.wt(Stocks[[ii]])[,simYears,,,,] <- SimRuns[[ii]]$catW[,1:(yEnd-yStart+1),]
    #landings.n(Stocks[[ii]])[,simYears,,,,] <- SimRuns[[ii]]$L[,1:(yEnd-yStart+1),]
    landings.n(Stocks[[ii]])[,simYears,,,,] <- SimRuns[[ii]]$C[,1:(yEnd-yStart+1),]
    stock.wt(Stocks[[ii]])[,simYears,,,,] <- SimRuns[[ii]]$stkW[,1:(yEnd-yStart+1),]
    mat(Stocks[[ii]])[,simYears,,,,] <- SimRuns[[ii]]$Mat[,1:(yEnd-yStart+1),]
    discards.wt(Stocks[[ii]]) <- catch.wt(Stocks[[ii]])
    landings.wt(Stocks[[ii]]) <- catch.wt(Stocks[[ii]])
    discards.n(Stocks[[ii]])[,simYears,,,,] <- 0
    catch(Stocks[[ii]])[,simYears] <- apply(catch.n(Stocks[[ii]])[,simYears]*catch.wt(Stocks[[ii]])[,simYears],2:6,sum)
    landings(Stocks[[ii]])[,simYears] <- apply(landings.n(Stocks[[ii]])[,simYears]*landings.wt(Stocks[[ii]])[,simYears],2:6,sum)
    discards(Stocks[[ii]])[,simYears] <- apply(discards.n(Stocks[[ii]])[,simYears]*discards.wt(Stocks[[ii]])[,simYears],2:6,sum)
    
    #statistics
    Stats <- list()
    
    #store the operating model and harvest rules
    Stats[["OM"]] <- OM
    Stats[["MP"]] <- MP
    Stats[["simYears"]] <- as.integer(simYears)
    Stats[["obsYears"]] <- as.integer(obsYears)
    
    #SSB
    SSB.true <- ssb(Stocks[[ii]])
    Stats[["SSB"]][["val"]] <- fStatPercs(SSB.true, lStatPer=lStatPer)
    Stats[["SSB"]][["worm"]] <- FLCore::iter(SSB.true,1:numWorm)
    
    #time to recovery after falling below Blim
    firstBelow <- recTimeBlim <- recTimeBpa <- rep(NA,dim(SSB.true)[6])
    names(firstBelow) <- names(recTimeBlim) <- names(recTimeBpa) <- as.character(seq(1,dim(SSB.true)[6]))
    
    #drop the unused dimensions
    SSB1 <- drop(as.array(SSB.true))
    
    #iterations during which SSB fell below Blim
    anyBelow <- apply(SSB1,MARGIN=2,FUN = function(x){any(x<OM$refPts$Blim)})
    #exclude those that never fell below Blim
    SSB1[,!anyBelow]<-NA
    recTimeBlim[names(anyBelow[!anyBelow])] <- NA
    recTimeBpa[names(anyBelow[!anyBelow])] <- NA
    
    #year in which SSB fell below Blim
    firstBelow[names(anyBelow[anyBelow])] <- apply(SSB1[,anyBelow], MARGIN=2, FUN = function(x){min(which(x<OM$refPts$Blim))})
    #remove SSB records prior to the year when it fell below Blim
    for (i in 1:length(firstBelow)){if (!is.na(firstBelow[i])) {SSB1[1:firstBelow[i],i]<-NA}}
    
    #Blim############
    #did any iterations recover above Blim?
    anyBackAboveBlim <- apply(SSB1[,names(anyBelow[anyBelow])], MARGIN=2, FUN = function(x){any(x>OM$refPts$Blim, na.rm=TRUE)})
    #set the recovery time for those that did not recover to 0
    recTimeBlim[names(anyBackAboveBlim[!anyBackAboveBlim])] <- 0
    
    #iterations which fell below and subsequently recovered
    recoveredBlim <- intersect(names(anyBelow[anyBelow]),names(anyBackAboveBlim[anyBackAboveBlim]))
    #get recovery time
    recTimeBlim[recoveredBlim] <- apply(SSB1[,recoveredBlim], MARGIN=2, FUN = function(x){min(which(x>OM$refPts$Blim))})
    
    #recovery time
    rBlim <- recTimeBlim-firstBelow
    
    Stats[["recBlim"]][["val"]] <- quantile(rBlim[rBlim>0],probs=percentiles,na.rm=T)
    Stats[["recBlim"]][["nobelow"]] <- sum(!anyBelow)
    Stats[["recBlim"]][["norecover"]] <- sum(!anyBackAboveBlim)
    
    #Bpa############
    #did any iterations recover above Bpa?
    anyBackAboveBpa <- apply(SSB1[,names(anyBelow[anyBelow])], MARGIN=2, FUN = function(x){any(x>OM$refPts$Bpa, na.rm=TRUE)})
    #set the recovery time for those that did not recover to 0
    recTimeBpa[names(anyBackAboveBpa[!anyBackAboveBpa])] <- 0
    
    #iterations which fell below and subsequently recovered
    recoveredBpa <- intersect(names(anyBelow[anyBelow]),names(anyBackAboveBpa[anyBackAboveBpa]))
    #get recovery time
    recTimeBpa[recoveredBpa] <- apply(SSB1[,recoveredBpa], MARGIN=2, FUN = function(x){min(which(x>OM$refPts$Bpa))})
    
    #recovery time  
    rBpa <- recTimeBpa-firstBelow
    
    Stats[["recBpa"]][["val"]] <- quantile(rBpa[rBpa>0],probs=percentiles,na.rm=T)
    Stats[["recBpa"]][["nobelow"]] <- sum(!anyBelow)
    Stats[["recBpa"]][["norecover"]] <- sum(!anyBackAboveBpa)
    
    #SSB error
    tSSB <- FLQuant(SimRuns[[ii]]$SSBratio[1:(yEnd-yStart+1),],
                    dim = c(1,yEnd-yStart+1,1,1,1,niters),
                    dimnames = list(age="all",year=ac(seq(yStart,yEnd)),unit="unique",season="all",
                                    area="unique",iter=ac(seq(1,niters))))
    
    Stats[["SSBratio"]][["val"]] <- fStatPercs(tSSB, lStatPer = lStatPer)
    Stats[["SSBratio"]][["worm"]] <- FLCore::iter(tSSB,1:numWorm)
    
    SSB.dev <- SimRuns[[ii]][["SSBdev"]]
    Stats[["SSB.dev"]] <- SSB.dev
    
    #FBar - realised F
    FBar <- fbar(Stocks[[ii]])
    Stats[["FBar"]][["val"]] <- fStatPercs(FBar, lStatPer=lStatPer)
    Stats[["FBar"]][["worm"]] <- FLCore::iter(FBar,1:numWorm)
    
    #FBar error
    tFBar <- FLQuant(SimRuns[[ii]]$Fratio[1:(yEnd-yStart+1),],
                     dim = c(1,yEnd-yStart+1,1,1,1,niters),
                     dimnames = list(age="all",year=ac(seq(yStart,yEnd)),unit="unique",season="all",
                                     area="unique",iter=ac(seq(1,niters))))
    
    #browser()
    Stats[["Fratio"]][["val"]] <- fStatPercs(tFBar, lStatPer = lStatPer)
    Stats[["Fratio"]][["worm"]] <- FLCore::iter(tFBar,1:numWorm)
    
    Fdev <- SimRuns[[ii]][["Fdev"]]
    Stats[["Fdev"]] <- Fdev
    
    #yield
    Catch <- catch(Stocks[[ii]])
    Stats[["Catch"]][["val"]] <- fStatPercs(Catch, lStatPer=lStatPer)
    Stats[["Catch"]][["worm"]] <- FLCore::iter(Catch,1:numWorm)
    
    #TAC
    tTAC <- FLQuant(SimRuns[[ii]]$TAC[1:(yEnd-yStart+1),],
                    dim = c(1,yEnd-yStart+1,1,1,1,niters),
                    dimnames = list(age="all",year=ac(seq(yStart,yEnd)),unit="unique",season="all",
                                    area="unique",iter=ac(seq(1,niters))))
    
    Stats[["TAC"]][["val"]] <- fStatPercs(tTAC, lStatPer = lStatPer)
    Stats[["TAC"]][["worm"]] <- FLCore::iter(tTAC,1:numWorm)
    
    #IAV - multiplied by 100 this stat is the absolute percentage change (no indication of up or down)
    IAV <- abs(1-Catch[,as.character(seq(yStart+1,yEnd))]/Catch[,as.character(seq(yStart,yEnd-1))])
    #replace Inf with NA (NA results from comparing with zero catch)
    IAV <- ifelse(is.finite(IAV),IAV,NA)
    
    Stats[["IAV"]][["val"]] <- fStatPercs(IAV, lStatPer = lStatPer2)
    Stats[["IAV"]][["worm"]] <- FLCore::iter(IAV,1:numWorm)
    
    #IAV increases/decreases
    IAVup <- IAVdown <- IAVupdown <- Catch[,as.character(seq(yStart+1,yEnd))]/Catch[,as.character(seq(yStart,yEnd-1))] - 1
    IAVup[IAVup<0] <- NA
    IAVdown[IAVdown>0] <- NA
    
    Stats[["IAVupdown"]][["worm"]] <- FLCore::iter(IAVupdown,1:numWorm)
    Stats[["IAVup"]][["val"]] <- fStatPercs(IAVup, lStatPer = lStatPer2)
    Stats[["IAVdown"]][["val"]] <- fStatPercs(IAVdown, lStatPer = lStatPer2)
    
    #Recruitment
    Rec <- rec(Stocks[[ii]])
    Stats[["Rec"]][["val"]] <- fStatPercs(Rec, lStatPer=lStatPer)
    Stats[["Rec"]][["worm"]] <- FLCore::iter(Rec,1:numWorm)
    
    #probability that SSB is below RP (should this be true or observed SSB?)
    Stats[["pBlim"]][["val"]] <- fStatRisk(SSB = SSB.true, RP = OM$refPts$Blim, lStatPer = lStatPer)
    Stats[["pBpa"]][["val"]] <- fStatRisk(SSB = SSB.true, RP = OM$refPts$Bpa, lStatPer = lStatPer)
    Stats[["pExt"]][["val"]] <- fStatExtinct(SSB = SSB.true, depletion=0.01, firstYear = maxObsYear)

    AllStats[[ac(ii)]] <- Stats
    
  }
  
  settings <- fGetSettings(
    stats = AllStats, lStatPer=lStatPer, SimRuns=SimRuns, 
    FLStockfile=FLStockfile, FLStockSimfile=FLStockSimfile,
    OM=OM, MP=MP, niters=niters, nyr=nyr)
  # save(settings,file = file.path(Res.dir,runName,paste0(runName,"_eqSim_Settings.Rdata")))
  
  df  <- fsummary_df(
    runName=runName, simRuns = simRuns, FLSs=FLSs,
    Res.dir = Res.dir, Plot.dir = Plot.dir,
    lStatPer = lStatPer, simYears = simYears, xlab = MP$xlab,
    OM = OM, MP = MP,
    Fbarrange=c(range(FLS)[["minfbar"]], range(FLS)[["maxfbar"]]),
    numWorm = numWorm, dfassess=dfassess) 
  
  ## Save data
  lStats <- list(stats = AllStats, runName = runName, lStatPer = lStatPer, 
                 OM = OM, MP = MP,
                 settings=settings, df=df)
  dir.create(path = file.path(Res.dir,"Stats"), showWarnings = TRUE, recursive = TRUE)
  
  cat("Saving eqSim Stats RData file", "\n")
  save(lStats,file = file.path(Res.dir,"Stats",paste0(runName,"_eqSim_Stats.Rdata")))
  
  # Save settings
  # settings <- fGetSettings(lStats, SimRuns, FLStockfile = FLStockfile,FLStockSimfile = FLStockSimfile,OM=OM, MP=MP, niters=niters, nyr=nyr)
  # save(settings,file = file.path(Res.dir,runName,paste0(runName,"_eqSim_Settings.Rdata")))
  
  #generate the stock/stat trajectories
  cat("Saving fPlotTraj files", "\n")
  fPlotTraj(sim = lStats, plot.dir = file.path(Res.dir,runName), lStatPer = lStatPer)
  
  cat("Saving fPlotSummary files", "\n")
  suppressWarnings(fPlotSummary(sim = lStats, plot.dir = Res.dir, lStatPer = lStatPer, FtoPlot=fGetValsScan(MP$F_target,OM$refPts)))
  
  cat("Saving fTabulateStats", "\n")
  fTabulateStats(sim = lStats, setting= settings, plot.dir = Res.dir)
  
} # end of for loop MPs


