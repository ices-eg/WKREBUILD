# ================================================================================================================
# 03 EqSim summarize
# 
# Summarize results of EqSim simulations
#
# 06/07/2020 tested on 1000 iters of SAM assessment
# ================================================================================================================

# rm(list=ls())
gc()

# library(tidyverse)

#basic display setting
# niters <- 1000
# maxyr <- 2040

stats.dir <- file.path(Res.dir, "Stats")

file.list <- list.files(path=stats.dir, pattern="eqSim_Stats", full.names=TRUE)

df <- data.frame(stringsAsFactors = FALSE)
for (i in 1:length(file.list)) {
  print(file.list[i])
  df <- bind_rows(df, loadRData(file.list[i])[["df"]])
}
df <- df %>%  mutate(simulator = "EQSIM")
df <- df %>%  mutate(period=ifelse(is.na(period) & year < (an(assessyear)-1), "HI",period)) 

rebuiltThreshold <- 0.5


# summary of simulations
tmp <-
  df %>% 
  dplyr::select(stock, assess, assessyear, om, mp, niters, nyrs, ftgt) %>% 
  distinct() %>% 
  group_by(stock, assess, assessyear, om, mp, niters, nyrs) %>% 
  summarise(
    ftgt     = paste(ftgt, collapse=",")
  ) 

sim_meta <-
  df %>% 
  dplyr::select(stock, assess, assessyear, om, mp, niters, nyrs, perfstat) %>% 
  distinct() %>% 
  group_by(stock, assess, assessyear, om, mp, niters, nyrs) %>% 
  summarise(
    perfstat     = paste(perfstat, collapse=",")
  ) %>% 
  left_join(tmp)


# mystock      =  "WHOM"
# myassess     = "SAM"
# myassessyear = c("2019","2020")
# myom         = c("OM2.4","OM2.5")
# myniters     = "1000"
# mymp         = c("MP5.23")
# mynyrs       = "23"
# myftgt       = c(0.0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15)
# myperfstat   = "stock"
# mycolour     = "blue"
# myyintercept = 834000
# myvalue      = "median"
# myfacets     = c("assessyear","ftgt")
# myperfstat   = "pblim"
# myvalue      = "mean"

# mystock      =  "WHOM"
# myassess     = c("SS3")
# myassessyear = c("2019")
# mysimulator  = "EQSIM"
# myom         = c("OM2.2")
# myniters     = "1000"
# mymp         = c("MP5.23")
# mynyrs       = "23"
# myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15)
# myperfstat   = "stock"
# mycolour     = "blue"
# myyintercept =  611814
# myvalue      = "median"
# myfacets     = c("mp","ftgt")
# myfirstyear  = 2000
# mylastyear   = as.numeric(NA)

# plot function
plotvar <- function(mystock      =  "WHOM",
                    myassess     = "SAM",
                    myassessyear = c("2019","2020"),
                    mysimulator  = "EQSIM",
                    myom         = c("OM2.4","OM2.5"),
                    myniters     = "1000",
                    mymp         = c("MP5.23"),
                    mynyrs       = "23",
                    myftgt       = c(0.0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15),
                    myperfstat   = "stock",
                    mycolour     = "blue",
                    myyintercept = 834000,
                    myvalue      = "median",
                    myfacets     = c("assessyear","ftgt"),
                    myfirstyear  = as.numeric(NA),
                    mylastyear   = as.numeric(NA)) {
  
  
  # history
  h <-
    df %>% 
    filter(period     %in% c("HI")) %>% 
    filter(perfstat   %in% myperfstat) %>% 
    filter(stock      %in% mystock) %>% 
    filter(assess     %in% myassess) %>% 
    filter(assessyear %in% myassessyear) %>% 
    filter(simulator  %in% mysimulator) %>% 
    filter(om         %in% myom) %>% 
    filter(mp         %in% mymp) %>% 
    filter(nyrs       %in% mynyrs) %>% 
    filter(ftgt       %in% myftgt) %>% 
    
    {if(!is.na(myfirstyear)) filter(., year >= myfirstyear) else (.)} %>%    
    {if(!is.na(mylastyear))  filter(., year <= mylastyear) else (.)} %>%   
    
    mutate(age = as.numeric(age)) 

  # future
  d <- 
    df %>%
    filter(period     %in% c("CU","ST","MT",'LT')) %>% 
    filter(perfstat   %in% myperfstat) %>% 
    filter(stock      %in% mystock) %>% 
    filter(assess     %in% myassess) %>% 
    filter(assessyear %in% myassessyear) %>% 
    filter(simulator  %in% mysimulator) %>% 
    filter(om         %in% myom) %>% 
    filter(mp         %in% mymp) %>%
    filter(nyrs       %in% mynyrs) %>% 
    filter(ftgt       %in% myftgt) %>% 
    
    {if(!is.na(myfirstyear)) filter(., year >= myfirstyear) else (.)} %>%    
    {if(!is.na(mylastyear))  filter(., year <= mylastyear) else (.)} %>%   
    
    mutate(age  = as.numeric(age)) 
    
  # periods
  p <-
    bind_rows(h, d) %>% 
    distinct(assessyear, period, year) %>% 
    
    {if(!is.na(myfirstyear)) filter(., year >= myfirstyear) else (.)} %>%    
    {if(!is.na(mylastyear))  filter(., year <= mylastyear) else (.)} %>%   
    
    group_by(assessyear, period) %>% 
    summarise(
      minyear = min(year),
      maxyear = max(year)
    )
  
  # title
  t <- paste(
    toupper(mystock),
    paste(unique(myassess), collapse="-"),
    paste(unique(myassessyear), collapse="-"),
    paste(unique(mysimulator), collapse="-"),
    paste(unique(myom), collapse="-"),
    paste(unique(myniters), collapse="-"),
    paste(unique(mymp), collapse="-"),
    paste(unique(mynyrs), collapse="-"),
    toupper(myperfstat),
    paste(c(myfacets[2],myfacets[1]), collapse="-"),
    paste(c(max(min(h$year), myfirstyear, na.rm=TRUE),
            min(max(d$year), mylastyear, na.rm=TRUE)), 
          collapse="-"),
    sep="_"
  )
    
  myfig <- 

    ggplot(d, aes(x=year, y=get(myvalue), group=mp)) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
    theme( panel.grid.major.x = element_blank()) +
    theme(legend.position = "none") +
    
    # historical ribbon
    geom_ribbon(data=filter(h, !is.na(get(myvalue))),
                aes(ymin=lower, ymax=upper, group=iter), alpha=0.4, fill="gray") +
    
    # historical worms
    geom_line(data=filter(h, metric=="worm"), 
              aes(x=year, y=data, group=iter), colour="darkgray", size=0.5) +
    
    # historical median
    geom_line(data=filter(h, !is.na(get(myvalue))), 
              aes(x=year,y=get(myvalue),group=iter), colour="black", size=1) +
    
    # future ribbon
    geom_ribbon(data=filter(d, !is.na(get(myvalue))),
                aes(ymin=lower, ymax=upper, group=iter), alpha=0.4, fill=mycolour) +
    
    # future worms
    geom_line(data=filter(d, metric=="worm"),
              aes(x=year, y=data, group=iter), colour=mycolour, size=0.5) +
    
    # future median
    geom_line(data=filter(d, !is.na(get(myvalue))),
              aes(x=year,y=get(myvalue),group=iter), colour=mycolour, size=1) +
    
    # hline
    {if(!all(is.na(myyintercept))) geom_hline(yintercept=myyintercept, linetype="dashed")}  +
    
    # vertical lines (periods)
    geom_vline(data=filter(p, period=="CU"),
               aes(xintercept=minyear), linetype="dotted") +
    geom_vline(data=filter(p, period=="ST"),
               aes(xintercept=minyear), linetype="dotted") +
    geom_vline(data=filter(p, period=="MT"),
               aes(xintercept=minyear), linetype="dotted") +
    geom_vline(data=filter(p, period=="LT"),
               aes(xintercept=minyear), linetype="dotted") +
    geom_vline(data=filter(p, period=="LT"),
               aes(xintercept=maxyear), linetype="dotted") +
    
    expand_limits(y=0) +
    labs(x="", y=myperfstat, title=t) +
    facet_grid(get(myfacets[1]) ~ get(myfacets[2])) 
    
  print(myfig)
  
  ggsave(file = file.path(Res.dir,paste0(t, " summaryplot.png")),
         device="png", width = 30, height = 20, units = "cm")
  
}

plotvar(mystock      = "WHOM",
        myassess     = "SS3",
        myassessyear = c("2019"),
        mysimulator  = "EQSIM",
        myom         = c("OM2.2"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = "stock",
        mycolour     = "blue",
        myyintercept = 611814,
        myvalue      = "median",
        myfacets     = c("assessyear","ftgt"),
        myfirstyear  = 2000)

plotvar(mystock      = "WHOM",
        myassess     = "SS3",
        myassessyear = c("2019"),
        mysimulator  = "EQSIM",
        myom         = c("OM2.2"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = "stock",
        mycolour     = "blue",
        myyintercept = 611814,
        myvalue      = "median",
        myfacets     = c("om","ftgt"),
        myfirstyear  = 2000)

plotvar(mystock      = "WHOM",
        myassess     = c("SAM"),
        myassessyear = c("2019"),
        mysimulator  = "EQSIM",
        myom         = c("OM2.4", "OM2.4.WR"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = "catch",
        mycolour     = "blue",
        myyintercept = 100000,
        myvalue      = "median",
        myfacets     = c("om","ftgt"),
        myfirstyear  = 2000)

plotvar(mystock      = "WHOM",
        myassess     = "SS3",
        myassessyear = c("2019","2019"),
        mysimulator  = "EQSIM",
        myom         = c("OM2.2","OM2.2.RR.5lowest"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = "stock",
        mycolour     = "blue",
        myyintercept = 611814,
        myvalue      = "median",
        myfacets     = c("om","ftgt"),
        myfirstyear  = 2000)

plotvar(mystock      =  "WHOM",
        myassess     = "SAM",
        myassessyear = c("2019","2020"),
        myom         = c("OM2.4","OM2.5"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.125, 0.15),
        myperfstat   = "fbar",
        mycolour     = "blue",
        myyintercept = 0.074,
        myvalue      = "median",
        myfacets     = c("assessyear","ftgt")
)

plotvar(mystock      =  "WHOM",
        myassess     = "SS3",
        myassessyear = c("2019","2020"),
        myom         = c("OM2.2","OM2.3"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.125, 0.15),
        myperfstat   = "fbar",
        mycolour     = "blue",
        myyintercept = 0.074,
        myvalue      = "median",
        myfacets     = c("assessyear","ftgt")
)

plotvar(mystock      =  "WHOM",
        myassess     = "SAM",
        myassessyear = c("2019","2020"),
        myom         = c("OM2.4","OM2.5"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = "rec",
        mycolour     = "blue",
        myyintercept = NA,
        myvalue      = "median",
        myfacets     = c("assessyear","ftgt")
)

plotvar(mystock      =  "WHOM",
        myassess     = "SAM",
        myassessyear = c("2019","2020"),
        myom         = c("OM2.4","OM2.5"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = "catch",
        mycolour     = "blue",
        myyintercept = NA,
        myvalue      = "median",
        myfacets     = c("assessyear","ftgt")
)

plotvar(mystock      =  "WHOM",
        myassess     = "SAM",
        myassessyear = c("2019","2020"),
        myom         = c("OM2.4","OM2.5"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = "pblim",
        mycolour     = "blue",
        myyintercept = NA,
        myvalue      = "median",
        myfacets     = c("assessyear","ftgt")
)

plotvar(mystock      =  "WHOM",
        myassess     = "SS3",
        myassessyear = c("2019","2020"),
        myom         = c("OM2.2","OM2.3"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = "pblim",
        mycolour     = "blue",
        myyintercept = NA,
        myvalue      = "median",
        myfacets     = c("assessyear","ftgt")
)

plotvar(mystock      =  "WHOM",
        myassess     = "SS3",
        myassessyear = c("2019","2019"),
        myom         = c("OM2.2","OM2.2.RR"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = "pblim",
        mycolour     = "blue",
        myyintercept = NA,
        myvalue      = "median",
        myfacets     = c("myom","ftgt")
)


plotvar(mystock      =  "WHOM",
        myassess     = "SS3",
        myassessyear = c("2019","2019"),
        myom         = c("OM2.2","OM2.2.RR5lowest"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = "pblim",
        mycolour     = "blue",
        myyintercept = NA,
        myvalue      = "median",
        myfacets     = c("assessyear","ftgt")
)

plotvar(mystock      =  "WHOM",
        myassess     = "SAM",
        myassessyear = c("2019","2020"),
        myom         = c("OM2.4","OM2.5"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = c("recovbpa"),
        mycolour     = "blue",
        myyintercept = NA,
        myvalue      = "median",
        myfacets     = c("assessyear","ftgt")
)

plotvar(mystock      =  "WHOM",
        myassess     = "SAM",
        myassessyear = c("2019"),
        myom         = c("OM2.4"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.075),
        myperfstat   = c("catch.n"),
        mycolour     = "blue",
        myyintercept = NA,
        myvalue      = "median",
        myfacets     = c("ftgt","age")
)

plotvar(mystock      =  "WHOM",
        myassess     = "SAM",
        myassessyear = c("2019"),
        myom         = c("OM2.4"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.075),
        myperfstat   = c("catch.wt"),
        mycolour     = "blue",
        myyintercept = NA,
        myvalue      = "median",
        myfacets     = c("ftgt","age")
)

plotvar(mystock      =  "WHOM",
        myassess     = "SAM",
        myassessyear = c("2019"),
        myom         = c("OM2.4"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.075),
        myperfstat   = c("stock.wt"),
        mycolour     = "blue",
        myyintercept = NA,
        myvalue      = "median",
        myfacets     = c("ftgt","age")
)

plotvar(mystock      =  "WHOM",
        myassess     = "SAM",
        myassessyear = c("2019"),
        myom         = c("OM2.4"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.075),
        myperfstat   = c("mat"),
        mycolour     = "blue",
        myyintercept = NA,
        myvalue      = "median",
        myfacets     = c("ftgt","age")
)

plotvar(mystock      =  "WHOM",
        myassess     = "SAM",
        myassessyear = c("2019"),
        myom         = c("OM2.4"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.075),
        myperfstat   = c("sel"),
        mycolour     = "blue",
        myyintercept = NA,
        myvalue      = "median",
        myfacets     = c("ftgt","age")
)




plotvar(mystock      =  "WHOM",
        myassess     = "SAM",
        myassessyear = c("2019"),
        myom         = c("OM2.4"),
        myniters     = "1000",
        mymp         = c("MP5.03", "MP5.13","MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = "stock",
        mycolour     = "blue",
        myyintercept = 611814,
        myvalue      = "median",
        myfacets     = c("mp","ftgt")
)


plotvar(mystock      =  "WHOM",
        myassess     = c("SS3"),
        myassessyear = c("2019"),
        mysimulator  = "EQSIM",
        myom         = c("OM2.2.WR","OM2.2"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = "stock",
        mycolour     = "blue",
        myyintercept = 611814,
        myvalue      = "median",
        myfacets     = c("om","ftgt"),
        myfirstyear  = 2000)

plotvar(mystock      =  "WHOM",
        myassess     = c("SAM"),
        myassessyear = c("2019"),
        mysimulator  = "EQSIM",
        myom         = c("OM2.4.WR","OM2.4"),
        myniters     = "1000",
        mymp         = c("MP5.23"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = "stock",
        mycolour     = "blue",
        myyintercept =  834480,
        myvalue      = "median",
        myfacets     = c("om","ftgt"),
        myfirstyear  = 2000)

plotvar(mystock      =  "WHOM",
        myassess     = c("SS3"),
        myassessyear = c("2019"),
        mysimulator  = "EQSIM",
        myom         = c("OM2.2"),
        myniters     = "1000",
        mymp         = c("MP5.23","MP5.23.DU"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = "stock",
        mycolour     = "blue",
        myyintercept =  611814,
        myvalue      = "median",
        myfacets     = c("mp","ftgt"),
        myfirstyear  = 2000)

plotvar(mystock      =  "WHOM",
        myassess     = c("SAM"),
        myassessyear = c("2019"),
        mysimulator  = "EQSIM",
        myom         = c("OM2.4"),
        myniters     = "1000",
        mymp         = c("MP5.23","MP5.23.DU"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = "stock",
        mycolour     = "blue",
        myyintercept =  611814,
        myvalue      = "median",
        myfacets     = c("mp","ftgt"),
        myfirstyear  = 2000)


plotvar(mystock      =  "WHOM",
        myassess     = c("SAM"),
        myassessyear = c("2019"),
        mysimulator  = "EQSIM",
        myom         = c("OM2.4"),
        myniters     = "1000",
        mymp         = c("MP5.23","MP5.23.def"),
        mynyrs       = "23",
        myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
        myperfstat   = "stock",
        mycolour     = "blue",
        myyintercept =  611814,
        myvalue      = "median",
        myfacets     = c("mp","ftgt"),
        myfirstyear  = 2010,
        mylastyear   = 2030)




# ========================================================================================================

plotrecovery <- function(
                    mystock      =  "WHOM",
                    myassess     = "SAM",
                    myassessyear = c("2019","2020"),
                    mysimulator  = "EQSIM",
                    myom         = c("OM2.4","OM2.5"),
                    myniters     = "1000",
                    mymp         = c("MP5.23"),
                    mynyrs       = "23",
                    myftgt       = c(0.0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15),
                    myyintercept = 0.5,
                    myfacets     = c("assessyear","ftgt"),
                    myfirstyear  = as.numeric(NA),
                    mylastyear   = as.numeric(NA)) {
  
  # plot recovery to blim and bpa
  d <-
    df %>%
    filter(perfstat %in% c("precblim","precbpa")) %>%
    filter(stock      %in% mystock) %>% 
    filter(assess     %in% myassess) %>% 
    filter(assessyear %in% myassessyear) %>% 
    filter(simulator  %in% mysimulator) %>% 
    filter(om         %in% myom) %>% 
    filter(mp         %in% mymp) %>% 
    filter(nyrs       %in% mynyrs) %>% 
    filter(ftgt       %in% myftgt) %>% 
    
    {if(!is.na(myfirstyear)) filter(., year >= myfirstyear) else (.)} %>%    
    {if(!is.na(mylastyear))  filter(., year <= mylastyear) else (.)} 
  
  # lines
  p <-
    df %>%
    filter(perfstat %in% c("firstyearrebuildtoblim","firstyearrebuildtobpa")) %>%
    filter(stock      %in% mystock) %>% 
    filter(assess     %in% myassess) %>% 
    filter(assessyear %in% myassessyear) %>% 
    filter(simulator  %in% mysimulator) %>% 
    filter(om         %in% myom) %>% 
    filter(mp         %in% mymp) %>% 
    filter(nyrs       %in% mynyrs) %>% 
    filter(ftgt       %in% myftgt) %>% 
    
    mutate(perfstat = case_when(
      perfstat == "firstyearrebuildtoblim" ~ "precblim",
      perfstat == "firstyearrebuildtobpa"  ~ "precbpa",
      TRUE                                 ~ perfstat
    ))

  # title
  t <- paste(
    toupper(mystock),
    paste(unique(myassess), collapse="-"),
    paste(unique(myassessyear), collapse="-"),
    paste(unique(mysimulator), collapse="-"),
    paste(unique(myom), collapse="-"),
    paste(unique(myniters), collapse="-"),
    paste(unique(mymp), collapse="-"),
    paste(unique(mynyrs), collapse="-"),
    "RECOVERY",
    paste(c(myfacets[2],myfacets[1]), collapse="-"),
    sep="_"
  )
  
  myfig <-
    d %>% 
    ggplot(aes(x=year, y=mean, group=perfstat)) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
    theme( panel.grid.major.x = element_blank()) +
    
    geom_vline(data=p, aes(xintercept=mean, colour=perfstat), linetype="dashed") +
    
    geom_hline(yintercept=myyintercept, colour="gray", linetype="dashed") +
    
    geom_text(data=filter(p, perfstat %in% c("precblim")),
              aes(x=mean, label=substr(ac(mean),3,4), colour=perfstat),
              y=1.0, vjust=1, hjust=1, nudge_x = -1) +
    
    geom_text(data=filter(p, perfstat %in% c("precbpa")),
              aes(x=mean, label=substr(ac(mean),3,4), colour=perfstat),
              y=0.0, vjust=0, hjust=0, nudge_x = 1) +
    
    geom_line(aes(colour=perfstat)) +
    
    expand_limits(y=0) +
    labs(x="", y="", title=t) +
    facet_grid(get(myfacets[1]) ~ get(myfacets[2])) 
  
  print(myfig)
  
  ggsave(file = file.path(Res.dir,paste0(t, " recoveryplot.png")),
         device="png", width = 30, height = 20, units = "cm")
  
  
} # end of plotrecovery
  
plotrecovery(
  mystock      =  "WHOM",
  myassess     = "SAM",
  myassessyear = c("2019","2020"),
  mysimulator  = "EQSIM",
  myom         = c("OM2.4","OM2.5"),
  myniters     = "1000",
  mymp         = c("MP5.23"),
  mynyrs       = "23",
  myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
  myyintercept = 0.5,
  myfacets     = c("assessyear","ftgt")
)
  
plotrecovery(
  mystock      =  "WHOM",
  myassess     = c("SS3", "SAM"),
  myassessyear = c("2019"),
  mysimulator  = "EQSIM",
  myom         = c("OM2.2","OM2.4"),
  myniters     = "1000",
  mymp         = c("MP5.23"),
  mynyrs       = "23",
  myftgt       = c(0.0, 0.05, 0.075, 0.10, 0.15),
  myyintercept = 0.5,
  myfacets     = c("assess","ftgt")
)



