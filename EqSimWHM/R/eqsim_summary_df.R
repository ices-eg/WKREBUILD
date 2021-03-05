#eqsim_summary_df
# Explort a single run into a df

# run=runName
# FLSs = FLSs
# OM = OM
# numWorm = numWorm

fassess_df <- function(runName, FLSs, OM, numWorm=5){
  
  #produces a dataframe

  require(tidyverse)
  options(dplyr.summarise.inform = FALSE)
  
  stock      <- stringr::word(runName,1,sep="_")
  assess     <- stringr::word(runName,2,sep="_")
  assessyear <- stringr::word(runName,3,sep="_")
  OMname     <- stringr::word(runName,4,sep="_")
  niters     <- an(stringr::word(runName,6,sep="_"))
  nyrs       <- an(stringr::word(runName,7,sep="_"))
  Blim       <- OM$refPts$Blim
  Bpa        <- OM$refPts$Bpa
  

  # historical assessment
  dfhistfbar <- 
    as.data.frame(fbar(FLSs)) %>% 
    mutate(slot="fbar") %>% 
    mutate_if(is.factor, as.character) 
  
  dfhistrec  <- 
    as.data.frame(rec(FLSs)) %>% 
    mutate(slot="rec", age=as.character(age)) %>% 
    mutate_if(is.factor, as.character) 
  
  
  dfhistpblim <- 
    as.data.frame(ssb(FLSs)) %>% 
    mutate(perfstat = "pblim") %>% 
    mutate(period = "HI") %>% 
    mutate_if(is.factor, as.character) %>% 
    dplyr::select(perfstat, age, year, period, value=data) %>%
    group_by(perfstat, age, year, period ) %>%
    summarize(mean   = sum((value< OM[["refPts"]][["Blim"]])) / niters) 
  
  x <-
    as.data.frame(FLSs) %>% 
    bind_rows(dfhistfbar) %>% 
    bind_rows(dfhistrec) %>% 
    mutate(period = "HI") 
  
  wormiters <-
    x %>% 
    distinct(iter) %>% 
    # arrange(iter) %>% 
    filter(row_number() <= 5)
  
  dfhist <- 
    x %>% 
    dplyr::select(perfstat=slot, age, year, period, value=data) %>%
    
    # filter(perfstat=="stock") %>% View()
    
    group_by(perfstat, age, year, period ) %>%
    summarize(mean   = mean(value, na.rm=TRUE),
              median = median(value, na.rm=TRUE),
              upper  = quantile(value, probs=0.975, na.rm=TRUE),
              lower  = quantile(value, probs=0.025, na.rm=TRUE)) %>%
    
    bind_rows(dfhistpblim) %>% 
    
    mutate(method = "EqSim") %>%
    mutate(
      runname   = runName,
      stock    = stock,
      assess   = assess,
      assessyear = assessyear,
      om       = OMname,
      niters   = niters,
      nyrs     = nyrs,
      blim     = Blim,
      bpa      = Bpa) %>%
    
    ungroup()
  
  wormhist <- 
    x %>% 
    filter(iter %in% wormiters$iter) %>% 
    dplyr::select(perfstat=slot, age, year, period, iter, value=data)  %>%
    mutate(method = "EqSim") %>%
    mutate(
      runname   = runName,
      stock    = stock,
      assess   = assess,
      assessyear = assessyear,
      om       = OMname,
      niters   = niters,
      nyrs     = nyrs,
      blim     = Blim,
      bpa      = Bpa) %>%
    
    ungroup()
  
  bind_rows(dfhist, wormhist)
  
} # end of fassess_df function




# ==========================================================================================================


runName=runName
simRuns = SimRuns
FLSs = FLSs
Res.dir = Res.dir
Plot.dir = Plot.dir
lStatPer = lStatPer
simYears = simYears
xlab = MP$xlab
OM = OM
Fbarrange=c(range(FLS)[["minfbar"]], range(FLS)[["maxfbar"]])
numWorm = numWorm
dfassess = dfassess

fsummary_df <- function(runName, simRuns, FLSs, 
                        Res.dir, Plot.dir, 
                        lStatPer, OM,  
                        simYears, xlab, 
                        Fbarrange=c(1,10),
                        numWorm=5, dfassess){

  #produces a dataframe
  #comparing the supplied performance statistic for ST,MT and LT for each runName
  #grouped by runName
  
  require(tidyverse)
  options(dplyr.summarise.inform = FALSE)
  
  stock  <- stringr::word(runName,1,sep="_")
  assess <- stringr::word(runName,2,sep="_")
  assessyear <- stringr::word(runName,3,sep="_")
  OMname <- stringr::word(runName,4,sep="_")
  MPname <- stringr::word(runName,5,sep="_")
  niters <- an(stringr::word(runName,6,sep="_"))
  nyrs   <- an(stringr::word(runName,7,sep="_"))
  Blim = OM$refPts$Blim
  Bpa  = OM$refPts$Bpa
  
  years <-
    data.frame(period="HI", year=seq(an(lStatPer$HI[1]),
                                     an(lStatPer$HI[2])),
               stringsAsFactors = FALSE) %>% 
    bind_rows(data.frame(period="CU", year=seq(an(lStatPer$CU[1]),
                                               an(lStatPer$CU[2])),
                         stringsAsFactors = FALSE)) %>% 
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
    perfstat = c("TAC"   ,"SSB"   ,"CW"    ,"Harvest","Rec"     ,"IAV", "IAVUp","IAVDown", "pblim", "pbpa", "catW"),
    unit     = c("tonnes","tonnes","tonnes","1/year" ,"Millions","perc","perc" ,"perc"   , "prob" , "prob", "kg"),
    stringsAsFactors = FALSE
  )
  
  dfall <- df <- worms <- data.frame(stringsAsFactors = FALSE)
  
  # Find out how to use the stock units
  # l <-
  #   units(FLSs) %>% 
  #   data.frame(t(matrix(unlist(l), ncol=length(l), byrow=FALSE)))
  
  
  # ftgt <- 0.0
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

    wormiters <-
      as.data.frame.table(t[[("C")]], responseName = "value", stringsAsFactors = FALSE) %>% 
      distinct(iter) %>% 
      mutate(iter=as.numeric(iter)) %>% 
      arrange(iter) %>% 
      filter(row_number() <= 5) %>% 
      mutate(iter=as.character(iter))
    
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
            perfstat   = item,
            year       = an(year),
            runname    = runName,
            stock      = stock,
            assess     = assess,
            assessyear = assessyear,
            method     = "EqSim",
            om         = OMname,
            mp         = MPname, 
            niters     = niters,
            nyrs       = nyrs,
            label      = xlab,
            ftgt       = ftgt,
            blim       = Blim,
            bpa        = Bpa) %>%

          left_join(years, by="year") %>% 
          left_join(units, by="perfstat") %>% 
          mutate(perfstat = case_when(
            perfstat == "C"     ~ "catch.n",
            perfstat == "N"     ~ "stock.n",
            perfstat == "F"     ~ "harvest",
            perfstat == "stkW"  ~ "stock.wt", 
            perfstat == "catW"  ~ "catch.wt", 
            perfstat == "Sel"   ~ "sel",
            TRUE    ~ perfstat))
        
        # print(head(x))
        
        df <- dplyr::bind_rows(df,x)
        
      
        } else if (item %in% c("SSB","CW","Fbar", "Fmgmt","Rec","IAV", "IAVUp", "IAVDown")) {
        # item <- "SSB"        
        x <- 
          as.data.frame(t(t[[(item)]])) %>% 
          rownames_to_column(var="iter") %>% 
          pivot_longer(cols=2:(nrow(t[[(item)]])+1), 
                       names_to = "year", 
                       values_to = "value") %>% 
          mutate(
            perfstat   = item,
            year       = an(year),
            runname    = runName,
            stock      = stock,
            assess     = assess,
            assessyear = assessyear,
            method     = "EqSim",
            om         = OMname,
            mp         = MPname, 
            niters     = niters,
            nyrs       = nyrs,
            label      = xlab,
            ftgt       = ftgt,
            blim       = Blim,
            bpa        = Bpa) %>%
          
          left_join(years, by="year") %>% 
          left_join(units, by="perfstat") %>% 
          mutate(perfstat = case_when(
            perfstat == "SSB"   ~ "stock",
            perfstat == "CW"    ~ "catch",
            perfstat == "Fbar"  ~ "fbar",
            perfstat == "Fmgmt" ~ "fmgmt",
            perfstat == "Rec"  ~ "rec",
            perfstat == "IAV"  ~ "iav", 
            perfstat == "IAVUp"  ~ "iavup", 
            perfstat == "IAVDown"  ~ "iavdown",
            TRUE    ~ perfstat))
        
        # print(head(x))
        
        df <- dplyr::bind_rows(df,x)
        
      } else if (item %in% c("pblim","pba")) {
        
        # item <- "pblim"        
        
        x <- 
          as.data.frame(t[[(item)]]) %>% 
          setNames("value") %>% 
          mutate(
            perfstat   = item,
            year       = lStatPer[["CU"]][1]:lStatPer[["LT"]][2],
            runname    = runName,
            stock      = stock,
            assess     = assess,
            assessyear = assessyear,
            method     = "EqSim",
            om         = OMname,
            mp         = MPname, 
            niters     = niters,
            nyrs       = nyrs,
            label      = xlab,
            ftgt       = ftgt,
            blim       = Blim,
            bpa        = Bpa) %>%
          
          left_join(years, by="year") %>% 
          left_join(units, by="perfstat") 
        
        # print(head(x))
        
        df <- dplyr::bind_rows(df,x)
        
      } else if (item %in% c("recovblim","recovbpa")){
        
        # item <- "recovblim"
        
        x <- 
          t[[(item)]] %>% 
          mutate(
            perfstat   = item,
            year       = an(year),
            runname    = runName,
            stock      = stock,
            assess     = assess,
            assessyear = assessyear,
            method     = "EqSim",
            om         = OMname,
            mp         = MPname, 
            niters     = niters,
            nyrs       = nyrs,
            label      = xlab,
            ftgt       = ftgt,
            blim       = Blim,
            bpa        = Bpa) %>%
          
          left_join(years, by="year") %>% 
          left_join(units, by="perfstat") 
        
        # print(head(x))
        
        df <- dplyr::bind_rows(df,x)
        
      } # end of if statement

    } # end of for loop (item)

    worms <-
      bind_rows(filter(df, iter %in% wormiters$iter)) 

    df <-
      df %>% 
      group_by(perfstat, age, year, period ) %>%
      summarize(mean   = mean(value, na.rm=TRUE),
                median = median(value, na.rm=TRUE),
                upper  = quantile(value, probs=0.975, na.rm=TRUE),
                lower  = quantile(value, probs=0.025, na.rm=TRUE)) %>%
      mutate(method = "EqSim") %>%
      mutate(
        runname   = runName,
        stock    = stock,
        assess   = assess,
        assessyear = assessyear,
        om       = OMname,
        mp       = MPname, 
        niters   = niters,
        nyrs     = nyrs,
        label    = xlab,
        ftgt     = ftgt,
        blim     = Blim,
        bpa      = Bpa) 
      
    dfall <-
      dfall %>% 
      
      # add history
      bind_rows(mutate(dfassess, ftgt=ftgt)) %>% 
      
      # add worms
      bind_rows(worms) %>% 
      
      # add summaries
      bind_rows(df)
    
  } # end of for loop: ftgt
    
  # dfAll <-
  #   dfAll %>% 
  #   tidyr::separate(runname, 
  #                   into=c("stock","assess", "OM","MP","niters","nyrs"), 
  #                   sep="_",
  #                   remove = FALSE) 
  
  
  dfall
  
} # end of function

# dfAll %>% 
#   filter(perfstat == "Yield") %>% 
#   mutate(ftgt = as.numeric(ftgt)) %>% 
#   mutate(Period = factor(Period, levels=c("ST","MT","LT"))) %>% 
#   ggplot(aes(x=ftgt, y=Val, group=Period)) +
#   theme_bw() +
#   geom_bar(stat="identity") +
#   facet_grid(runname ~ Period)





