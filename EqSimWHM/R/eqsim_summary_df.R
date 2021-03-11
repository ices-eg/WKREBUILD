#eqsim_summary_df

# runName=runName
# OM=OM
# MP=MP
# ftgt=as.numeric(ii)
# Stats=Stats
# mystock=Stocks[[ii]]
# lStatPer = lStatPer
# Fbarrange=c(range(FLS)[["minfbar"]], range(FLS)[["maxfbar"]])

fsummary_df <- function(runName, OM, MP, ftgt,
                        Stats, lStatPer, 
                        Fbarrange=c(1,10)){
  # runName=runName
  # OM=OM
  # MP=MP
  # ftgt=as.numeric(ii)
  # Stats=Stats
  # mystock=Stocks[[ii]]
  # lStatPer = lStatPer
  
  #produces a dataframe
  #comparing the supplied performance statistic for ST,MT and LT for each runName
  #grouped by runName
  
  require(tidyverse)
  options(dplyr.summarise.inform = FALSE)
  
  invisible(gc())
  
  stock      <- stringr::word(runName,1,sep="_")
  assess     <- stringr::word(runName,2,sep="_")
  assessyear <- stringr::word(runName,3,sep="_")
  niters     <- an(stringr::word(runName,5,sep="_"))
  nyrs       <- an(stringr::word(runName,7,sep="_"))
  blim       <- OM$refPts$Blim
  bpa        <- OM$refPts$Bpa
  
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
  
  #units CAN WE TAKE THOSE DIRECTLY FROM THE FLSTOCK OBJECT????
  units <- data.frame(
    perfstat = c("TAC"   ,"SSB"   ,"CW"    ,"Harvest","Rec"     ,"IAV", "IAVUp","IAVDown", "pblim", "pbpa", "catW"),
    unit     = c("tonnes","tonnes","tonnes","1/year" ,"Millions","perc","perc" ,"perc"   , "prob" , "prob", "kg"),
    stringsAsFactors = FALSE
  )
  
  # cat("Creating dataframe for run with f = ", ftgt, "\n")
  
  df <- 
    data.frame(stringsAsFactors = FALSE) %>% 
    
    # adding summary data
    bind_rows(
         as.data.frame(Stats$SSB$val)  %>%     
           mutate_if(is.factor, as.character) %>% 
           mutate(perfstat = "stock", iter="all") %>% 
           rename(metric=age),
         as.data.frame(Stats$Rec$val)  %>% 
           mutate_if(is.factor, as.character) %>% 
           mutate(perfstat = "rec",   iter="all") %>% 
           rename(metric=age),
         as.data.frame(Stats$FBar$val) %>% 
           mutate_if(is.factor, as.character) %>% 
           mutate(perfstat = "fbar",  iter="all") %>% 
           rename(metric=age),
         as.data.frame(Stats$Catch$val)%>% 
           mutate_if(is.factor, as.character) %>% 
           mutate(perfstat = "catch",  iter="all") %>% 
           rename(metric=age),
    ) %>% 
    mutate(metric = case_when(
        metric == "50%" ~ "median",
        metric == "5%" ~ "lower",
        metric == "95%" ~ "upper",
        TRUE           ~ metric)
    ) %>% 
    filter(metric %in% c("mean","median","lower","upper")) %>% 
    filter(!is.na(year)) %>% 
    tidyr::pivot_wider(names_from = metric, values_from=data) %>% 
    
    # now bind the worms
    bind_rows(
      as.data.frame(Stats$SSB$worm)  %>%     
        mutate_if(is.factor, as.character) %>% 
        mutate(perfstat = "stock", metric="worm") ,
      as.data.frame(Stats$Rec$worm)  %>% 
        mutate_if(is.factor, as.character) %>% 
        mutate(perfstat = "rec",   metric="worm", age=ac(age)),
      as.data.frame(Stats$FBar$val) %>% 
        mutate_if(is.factor, as.character) %>% 
        mutate(perfstat = "fbar",  metric="worm"),
      as.data.frame(Stats$Catch$val)%>% 
        mutate_if(is.factor, as.character) %>% 
        mutate(perfstat = "catch",  metric="worm"),
    ) %>% 
    filter(!is.na(year)) %>% 
    
    # add descriptors
    mutate(
      stock      = stock,
      assess     = assess,
      assessyear = assessyear, 
      ftgt       = an(ftgt),
      om         = OM$code,
      mp         = MP$code,
      niters     = niters,
      nyrs       = nyrs,
      blim       = blim,
      bpa        = bpa,
    )

} # end of function

# glimpse(df)
  
