# ============================================================================
# HOM SAM simulator
# ============================================================================

rm(list=ls())
gc()

# install.packages("devtools")

# install.packages("FLCore", repos="http://flr-project.org/R")
library(FLCore)

# devtools::install_github('fishfollower/SAM/stockassessment') 
library(stockassessment)

library(tidyverse)

basedir  <- "E:/MARTIN"
sao.name <- "WHOM_2019"
tempdir  <- file.path(basedir,sao.name) 

setwd(tempdir)

cn<-read.ices("data/cn.dat")
cw<-read.ices("data/cw.dat")
dw<-read.ices("data/dw.dat")
lf<-read.ices("data/lf.dat")
lw<-read.ices("data/lw.dat")
mo<-read.ices("data/mo.dat")
nm<-read.ices("data/nm.dat")
pf<-read.ices("data/pf.dat")
pm<-read.ices("data/pm.dat")
sw<-read.ices("data/sw.dat")
surveys<-read.ices("data/survey.dat")

# Read input in sam format
dat<-setup.sam.data(surveys=surveys,
                    residual.fleet=cn, 
                    prop.mature=mo, 
                    stock.mean.weight=sw, 
                    catch.mean.weight=cw, 
                    dis.mean.weight=dw, 
                    land.mean.weight=lw,
                    prop.f=pf, 
                    prop.m=pm, 
                    natural.mortality=nm, 
                    land.frac=lf)



conf <- loadConf(dat,"conf/model.cfg", patch=TRUE)
par  <- defpar(dat, conf)
fit  <- sam.fit(dat, conf, par)

# read input in VPA file format
inputs <- c(landings.fraction='lf.dat', catch.n='cn.dat', catch.wt='cw.dat',
            discards.wt='dw.dat', landings.wt='lw.dat', stock.wt='sw.dat',
            mat='mo.dat', m='nm.dat', harvest.spwn='pf.dat', m.spwn='pm.dat')
fqs <- as.list(inputs)
for (i in seq(inputs)) {
  file <- file.path("data",inputs[i])
  fqs[[names(inputs[i])]] <- readVPAFile(file)
}


minage       <- fit$data$minAge[[1]]
maxage       <- fit$data$maxAge[[1]]
pg           <- ifelse(fit$conf$maxAgePlusGroup[1]==1,fit$data$maxAge[[1]],NA)
minyear      <- min(fit$data$years)
maxyear      <- max(as.numeric(rownames(cn)))  # max year from data
maxyearsam   <- max(fit$data$years)            # max year from assessment
minfbar      <- fit$conf$fbarRange[1]
maxfbar      <- fit$conf$fbarRange[2]

# Generate FLStock object
FLS  <- FLStock()

# set generate properties
FLS@desc               <- paste("FLStock object generated from SAM:", date(), sep=" ")
FLS@name               <- sao.name
FLS@range["min"]       <- minage
FLS@range["max"]       <- maxage
FLS@range["plusgroup"] <- pg
FLS@range["minyear"]   <- minyear
FLS@range["maxyear"]   <- maxyear
FLS@range["minfbar"]   <- minfbar
FLS@range["maxfbar"]   <- maxfbar

units(FLS)             <- FLCore::standardUnits(FLS)

FLS@catch.n            <- fqs[["catch.n"]] 
FLS@landings.n         <- fqs[["catch.n"]] * fqs[["landings.fraction"]]
FLS@discards.n         <- fqs[["catch.n"]] * (1-fqs[["landings.fraction"]])
FLS@m                  <- fqs[["m"]] %>% FLCore::window(., end=maxyear)
FLS@mat                <- fqs[["mat"]] %>% FLCore::window(., end=maxyear)
FLS@catch.wt           <- fqs[["catch.wt"]]
FLS@landings.wt        <- fqs[["landings.wt"]]
FLS@discards.wt        <- fqs[["discards.wt"]]
FLS@stock.wt           <- fqs[["stock.wt"]] %>% FLCore::window(., end=maxyear)
FLS@harvest.spwn       <- fqs[["harvest.spwn"]] %>% FLCore::window(., end=maxyear) 
FLS@m.spwn             <- fqs[["m.spwn"]] %>% FLCore::window(., end=maxyear)
FLS@landings <- FLS@discards <- FLS@catch  <- FLS@stock <-
                         FLQuant(NA, dimnames=list(age="all", year=minyear:maxyear))
FLS@landings          <- quantSums(FLS@catch.n * FLS@catch.wt)
FLS@discards          <- quantSums(FLS@discards.n * FLS@discards.wt)
FLS@catch             <- quantSums(FLS@catch.n * FLS@catch.wt)

# stock numbers
FLS@stock.n          <- FLQuant(NA, dimnames=list(age=minage:maxage, year=minyear:maxyearsam))
FLS@stock.n[,]       <- exp(fit$pl$logN) 
FLS@stock.n          <- FLCore::window(FLS@stock.n, end=maxyear)

# harvest
n.ages               <- nrow(fit$pl$logF)
FLS@harvest          <- FLQuant(NA, dimnames=list(age=minage:maxage, year=minyear:maxyearsam))
FLS@harvest[minage:(minage+n.ages),] <- exp(fit$pl$logF)
FLS@harvest[(n.ages+1),]    <- FLS@harvest[n.ages,]
FLS@harvest          <- FLCore::window(FLS@harvest, end=maxyear)
units(FLS@harvest)   <-  "f"

# ssb
FLS@stock            <- ssb(FLS)

validObject(FLS)

# plot(FLS)
# plot(FLS@catch)

# set simulator properties
nsim     <- 1000
set.seed(123)

# Generate FLstock object with iterations
FLSs <- propagate(FLS, nsim)

# Simulate the observations
simdata <- simulate(fit, nsim=nsim,  full.data=TRUE)

# iterations. start with 2, 1 is the real observation.
i <- 708   

while (i <= nsim) {
  
  start_time <- Sys.time()
  
  # sim_fit <-  invisible(capture.output(sam.fit(simdata[[i]], conf, par)))
 sim_fit <-  sam.fit(simdata[[i]], conf, par)
  
  if(exists("sim_fit")) {
    
    # stock numbers
    stock.n                 <- FLQuant(NA, dimnames=list(age=minage:maxage, year=minyear:maxyearsam))
    stock.n[,]              <- exp(sim_fit$pl$logN) 
    stock.n                 <- FLCore::window(stock.n, end=maxyear)
    FLSs[,,,,,i]@stock.n    <- stock.n

    # harvest
    n.ages                  <- nrow(sim_fit$pl$logF)
    harvest                 <- FLQuant(NA, dimnames=list(age=minage:maxage, year=minyear:maxyearsam))
    harvest[minage:(minage+n.ages),] <- exp(sim_fit$pl$logF)
    harvest[(n.ages+1),]    <- harvest[n.ages,]
    FLSs[,,,,,i]@harvest    <- FLCore::window(harvest, end=maxyear)
    units(FLSs[,,,,,i]@harvest)   <-  "f"
    
    # ssb
    FLSs[,,,,,i]@stock            <- ssb(FLSs[,,,,,i])
    
    save(FLSs, file="run/FLSs.RData")
    rm(sim_fit)
    
    end_time <- Sys.time()
    duration = end_time - start_time
    print(paste("run", i, duration, sep=" "))
    
    i = i + 1
    
    
    
  }
} # end of while statement


# load(file="C:/TEMP/WHOM_2020/run/FLSs.RData")

df <-
  as.data.frame(FLSs) %>% 
  
  filter(slot=="harvest", age >= minage, age <= maxage) %>%
  group_by(year, unit, season, area, iter) %>% 
  summarise(data = mean(data, na.rm=TRUE)) %>% 
  mutate(
    slot="meanf",
    age ="1-10"
  ) %>% 
  bind_rows(as.data.frame(FLSs))

  
    


df %>% 
  filter(slot=="stock") %>% 
  ggplot(aes(year, data, group=iter)) +
  theme(legend.position="none") +
  geom_line(aes(colour=iter)) +
  geom_line(data=filter(df, slot=="stock", iter==1), colour="black", size=1)

df %>% 
  filter(slot=="meanf") %>%
  ggplot(aes(year, data, group=iter)) +
  theme(legend.position="none") +
  geom_line(aes(colour=iter)) +
  geom_line(data=filter(df, slot=="meanf", iter==1), colour="black", size=1)


# df %>% 
#   filter(run<=2) %>% 
#   ggplot(aes(x=year,y=data, group=run)) +
#   theme_bw() +
#   geom_line(colour="gray")

# plot(FLSs)
