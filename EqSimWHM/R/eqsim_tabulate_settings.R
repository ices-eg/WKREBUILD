fGetSettings <- function(dfStats, dfSimRuns){
  
  tData <- 
    # OM
    data.frame(unlist(dfStats[[1]][[1]][["OM"]], use.names=TRUE), stringsAsFactors = FALSE) %>% 
    tibble::rownames_to_column() %>% 
    setNames(c("desc","value")) %>%
    mutate(class="OM") %>% 
    
    # MP
    bind_rows(
      data.frame(unlist(dfStats[[1]][[1]][["MP"]], use.names=TRUE), stringsAsFactors = FALSE) %>% 
        tibble::rownames_to_column() %>% 
        setNames(c("desc","value")) %>% 
        mutate(class="MP")
    ) %>% 
    
    # Other
    bind_rows(
      data.frame(desc=c("niters","nyr"),
                 value=c(ac(dim(dfSimRuns[[1]][["N"]])[[3]]), ac(dim(dfSimRuns[[1]][["N"]])[[2]])),
                 class="OTHER",
                 stringsAsFactors = FALSE)) %>% 
    
    dplyr::select(class, desc, value)
  
  return(tData)
}

fMakeSettingsTable <- function(dfSettings){
  
  myft <- 
    flextable::flextable(data = dfSettings) %>% 
    flextable::theme_vanilla() %>%
    bold(part = "header") %>%
    autofit()
  
  return(myft)
}

