#eqsim_summary_df
# Explort a single run into a df

# runName <- c("WHOM_SS_OM2.2_MP2.1_1000_20")
# lStatPer = lStatPer
# Blim <- OM$refPts$Blim
# Fbarrange=c(1,10)
# stats <- AllStats

# run=runName
# simRuns = SimRuns
# Res.dir = Res.dir
# Plot.dir = Plot.dir
# lStatPer = lStatPer
# simYears = simYears
# xlab = MP$xlab
# OM = OM
# Fbarrange=c(range(FLS)[["minfbar"]], range(FLS)[["maxfbar"]])

fsummary_df <- function(runName, ftgt, simRuns,
                        Res.dir, Plot.dir, 
                        lStatPer, OM,  
                        simYears, xlab, 
                        Fbarrange=c(1,10)){

  #produces a dataframe
  #comparing the supplied performance statistic for ST,MT and LT for each runName
  #grouped by runName
  
  require(tidyverse)
  options(dplyr.summarise.inform = FALSE)
  
  dfAll <- data.frame(stringsAsFactors = FALSE)

  Blim = OM$refPts$Blim
  Bpa  = OM$refPts$Bpa
  
  #get the data
  # run <- "OM2.2_MP2.1_1000_20"
  # cat(run,"\n")
  
  #load the output of the simulation and the summary statistics
  # load(file = file.path(Res.dir,run,paste0(run,"_SimRuns.RData")))
  # load(file = file.path(Res.dir,run,paste0(run,"_eqSim_Stats.RData")))
  
  years <-
    data.frame(period="CU", year=seq(an(lStatPer$CU[1]),
                                     an(lStatPer$CU[2])),
               stringsAsFactors = FALSE) %>% 
    bind_rows(data.frame(period="ST", year=seq(an(lStatPer$ST[1]),
                                               an(lStatPer$ST[2])),
                         stringsAsFactors = FALSE)) %>% 
    bind_rows(data.frame(period="MT", year=seq(an(lStatPer$MT[1]),
                                               an(lStatPer$MT[2])),
                         stringsAsFactors = FALSE)) %>% 
    bind_rows(data.frame(period="LT", year=seq(an(lStatPer$LT[1]),
                                               an(lStatPer$LT[2])),
                         stringsAsFactors = FALSE)) %>% 
    mutate(year = an(year))
  
  #units
  units <- data.frame(
    PerfStat = c("TAC"   ,"SSB"   ,"CW"    ,"Harvest","Rec"     ,"IAV", "IAVUp","IAVDown", "pblim", "pbpa", "catW"),
    unit     = c("tonnes","tonnes","tonnes","1/year" ,"Millions","perc","perc" ,"perc"   , "prob" , "prob", "kg"),
    stringsAsFactors = FALSE
  )
  
  stock  <- stringr::word(runName,1,sep="_")
  assess <- stringr::word(runName,2,sep="_")
  OMname <- stringr::word(runName,3,sep="_")
  MPname <- stringr::word(runName,4,sep="_")
  niters <- stringr::word(runName,5,sep="_")
  nyrs   <- stringr::word(runName,6,sep="_")
  
  
  # ftgt <- 0.1
  for (ftgt in an(names(SimRuns))){
    
    invisible(gc())
    
    cat("Creating dataframe for run with f = ", ftgt, "\n")
    
    #simulation op
    t <- SimRuns[[ac(ftgt)]]

    #simulation stats
    # t2 <- lOp$stats[[ac(ftgt)]]
    # t2 <- stats[[ac(ftgt)]]

    # dimnames
    slots <- names(t)[!grepl("SimYears", names(t))]
    for (n in slots) {dimnames(t[[n]])$year <- simYears}
    
    #catch numbers
    Cnum <- t[["C"]]
    
    #catch weights
    Cwgt <- t[["catW"]]
    
    #catch weight (tons)
    CW <- apply(Cnum*Cwgt,2:3,sum) 
    # dimnames(CW)$year <- simYears
    t$CW <- CW

    # IAV
    IAV <- abs(1-CW[-1,]/CW[-nrow(CW),])
    IAV <- ifelse(is.finite(IAV),IAV,NA)
    t$IAV <- IAV
    
    # IAVup and IAVdown
    IAVUp <- IAVDown <- IAVUpDown <- (CW[-1,]/CW[-nrow(CW),]) - 1
    
    IAVUp[!is.finite(IAVUp)]<-NA
    IAVDown[!is.finite(IAVDown)]<-NA
    
    IAVUp[IAVUp<0] <- NA
    IAVDown[IAVDown>0] <- NA
    
    t$IAVUp <- IAVUp
    t$IAVDown <- IAVDown
    
    #F management
    # dimnames(t$Fmgmt)$year <- simYears
    
    #harvest
    Fbar <- t[["F"]]
    # dimnames(Harvest)$year <- simYears
    Fbar <- apply(Fbar[ac(Fbarrange[1]:Fbarrange[2]),,],2:3, mean)
    t$Fbar <- Fbar
    # dimnames(Harvest)
    
    #abundance
    Abd <- t[["N"]]
    
    #stock weights
    SW <- t[["stkW"]]
    
    #maturity
    Mat <- t[["Mat"]]
    
    #recruitment
    Rec <- t[["N"]][1,,]
    dimnames(Rec)$year <- simYears
    t$Rec <- Rec
    
    #SSB (Mt)
    SSB <- apply(Abd*SW*Mat,2:3,sum)
    dimnames(SSB)$year <- simYears
    t$SSB <- SSB

    t$pblim <- apply(SSB[]<Blim,1,sum)/dim(SSB)[2]
    t$pbpa  <- apply(SSB[]<Bpa,1,sum)/dim(SSB)[2]

    recov <- 
      as.data.frame(SSB) %>% 
      rownames_to_column(var="year") %>% 
      pivot_longer(cols=2:(ncol(SSB)+1), 
                   names_to = "iter", 
                   values_to = "v") %>% 
      group_by(iter) %>% 
      arrange(iter, year) %>% 
      mutate(v1 = lag(v, n=1),
             v2 = lag(v, n=2)) %>% 
      mutate(recovblim = ifelse(v2 >= Blim & v1 >= Blim & v >= Blim, TRUE, FALSE),
             recovbpa  = ifelse(v2 >= Bpa & v1 >= Bpa & v >= Bpa, TRUE, FALSE)) %>% 
      filter(!is.na(v1), !is.na(v2)) %>%
      
      group_by(year) %>% 
      mutate(n=max(as.numeric(iter))) %>% 
      group_by(year) %>% 
      summarize(
        n         = max(n, na.rm=TRUE),
        recovblim = sum(recovblim, na.rm=TRUE),
        recovbpa  = sum(recovbpa, na.rm=TRUE)
      ) %>% 
      mutate(
        recovblim = recovblim / n,
        recovbpa  = recovbpa / n)
    
    t$recovblim <- recov %>% dplyr::select(year, value=recovblim)
    t$recovbpa <- recov %>% dplyr::select(year, value=recovbpa)
    
    # item <- "C"
    # item <- "SSB"
    # item <- "catW"
    # item <- "IAV"
    # item <- "pblim"
    # item <- "recovblim" 
    for (item in c("stkW", "catW", "Sel", "N", "C", "F",
                   "SSB","CW", "Harvest", "Fbar", "Fmgmt","Rec",
                   "IAV", "IAVUp", "IAVDown", "pblim", "pbpa",
                   "recovblim", "recovbpa")) {
      
      # cat(item, "\n")
      
      if (item %in% c("N", "C", "F", "stkW", "catW", "Sel")) {
        
        x <- 
          as.data.frame.table(t[[(item)]], responseName = "value", stringsAsFactors = FALSE) %>% 
          # rownames_to_column(var="iter") %>% 
          # pivot_longer(cols=2:(nrow(t[[(item)]])+1), 
          #              names_to = "year", 
          #              values_to = "value") %>% 
          mutate(
            PerfStat = item,
            year     = an(year),
            RunRef   = runName,
            stock    = stock,
            assess   = assess,
            OM       = OMname,
            MP       = MPname, 
            niters   = niters,
            nyrs     = nyrs,
            Label = xlab,
            Ftgt = ftgt) %>%
          
          left_join(years, by="year") %>% 
          left_join(units, by="PerfStat") %>% 
          mutate(PerfStat = case_when(
            PerfStat == "C"     ~ "catch.n",
            PerfStat == "N"     ~ "stock.n",
            PerfStat == "F"     ~ "harvest",
            PerfStat == "stkW"  ~ "stock.wt", 
            PerfStat == "catW"  ~ "catch.wt", 
            PerfStat == "Sel"   ~ "sel",
            TRUE    ~ PerfStat))
        
        # print(head(x))
        
        dfAll <- dplyr::bind_rows(dfAll,x)
        
      
        } else if (item %in% c("SSB","CW",
                      "Fbar", "Fmgmt","Rec",
                      "IAV", "IAVUp", "IAVDown")) {
        
        x <- 
          as.data.frame(t(t[[(item)]])) %>% 
          rownames_to_column(var="iter") %>% 
          pivot_longer(cols=2:(nrow(t[[(item)]])+1), 
                       names_to = "year", 
                       values_to = "value") %>% 
          mutate(
            PerfStat = item,
            year     = an(year),
            RunRef = runName,
            stock    = stock,
            assess   = assess,
            OM       = OMname,
            MP       = MPname, 
            niters   = niters,
            nyrs     = nyrs,
            Label = xlab,
            Ftgt = ftgt) %>%
          
          left_join(years, by="year") %>% 
          left_join(units, by="PerfStat") %>% 
          mutate(PerfStat = case_when(
            PerfStat == "SSB"   ~ "stock",
            PerfStat == "CW"    ~ "catch",
            TRUE    ~ PerfStat))
        
        # print(head(x))
        
        dfAll <- dplyr::bind_rows(dfAll,x)
        
      } else if (item %in% c("pblim","pba")) {
        
        x <- 
          as.data.frame(t[[(item)]]) %>% 
          setNames("value") %>% 
          mutate(
            PerfStat = item,
            year     = years$year,
            RunRef   = runName,
            stock    = stock,
            assess   = assess,
            OM       = OMname,
            MP       = MPname, 
            niters   = niters,
            nyrs     = nyrs,
            Label    = xlab,
            Ftgt     = ftgt) %>%
          
          left_join(years, by="year") %>% 
          left_join(units, by="PerfStat") 
        
        # print(head(x))
        
        dfAll <- dplyr::bind_rows(dfAll,x)
        
      } else if (item %in% c("recovblim","recovbpa")){
        
        x <- 
          t[[(item)]] %>% 
          mutate(year = an(year)) %>% 
          mutate(
            PerfStat = item,
            RunRef = runName,
            stock    = stock,
            assess   = assess,
            OM       = OMname,
            MP       = MPname, 
            niters   = niters,
            nyrs     = nyrs,
            Label = xlab,
            Ftgt = ftgt) %>%
          
          left_join(years, by="year") %>% 
          left_join(units, by="PerfStat") 
        
        # print(head(x))
        
        dfAll <- dplyr::bind_rows(dfAll,x)
        
      }
        # end of if statement

    } # end of for loop (item)

  } # end of for loop: ftgt
    
  # dfAll <-
  #   dfAll %>% 
  #   tidyr::separate(RunRef, 
  #                   into=c("stock","assess", "OM","MP","niters","nyrs"), 
  #                   sep="_",
  #                   remove = FALSE) 
    
  dfAll
  
} # end of function

# dfAll %>% 
#   filter(PerfStat == "Yield") %>% 
#   mutate(Ftgt = as.numeric(Ftgt)) %>% 
#   mutate(Period = factor(Period, levels=c("ST","MT","LT"))) %>% 
#   ggplot(aes(x=Ftgt, y=Val, group=Period)) +
#   theme_bw() +
#   geom_bar(stat="identity") +
#   facet_grid(RunRef ~ Period)





