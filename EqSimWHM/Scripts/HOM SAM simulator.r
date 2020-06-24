# ============================================================================
# HOM SAM simulator
# ============================================================================

rm(list=ls())
gc()

library(FLCore)
library(FLSAM)
library(stockassessment)
# library(parallel)

# source SAM2FLSAM
source("../SAM2FLR/SAM2FLSAM.r")
source("../SAM2FLR/SAM2FLSTOCK.r")
# source("EqSimWHM/R//simulate.sam.r")
# source("EqSimWHM/R//simstudy.r")

sao.name <- "WHOM_2019"
url      <- paste0("https://stockassessment.org/datadisk/stockassessment/userdirs/user3/",sao.name,"/")
tempdir  <- file.path("D:/temp",sao.name) 
fit      <- stockassessment::fitfromweb(sao.name, character.only=TRUE) 
nsim     <- 100
set.seed(123)

dir.create(tempdir, showWarnings = FALSE)
write.data.files(fit, dir = tempdir)

inputs <- c(landings.fraction='lf.dat', catch.n='cn.dat', catch.wt='cw.dat',
            discards.wt='dw.dat', landings.wt='lw.dat', stock.wt='sw.dat',
            mat='mo.dat', m='nm.dat', harvest.spwn='pf.dat', m.spwn='pm.dat')

fqs <- as.list(inputs)

for (i in seq(inputs)) {
  file <- file.path("D:/temp",sao.name,inputs[i])
  fqs[[names(inputs[i])]] <- readVPAFile(file)
}

samfit        <- sam.fit(fit$data, fit$conf, fit$obj$env$par)
runs0001_0100 <- simstudy(samfit, nsim=nsim, ncores=1)
getwd()

load(file.path(getwd(),"EqSimWHM","RData","runs.RData"))
plot(runs)

# nsim   <- 90
# set.seed(124)
# samfit <- sam.fit(fit$data, fit$conf, fit$obj$env$par)
# runs0011_0100 <- simstudy(samfit, nsim=nsim, ncores=1)

# runs0001_0100 <- append(runs0001_0010, runs0011_0100)
# attributes(runs0001_0100)
# attr(runs0001_0100, "fit") <- fit
# attr(runs0001_0100, "class") <- "samset"
# str(runs0001_0100)
# plot(runs0001_0100)


flstocks <- list()
df       <- data.frame(NULL)
i <- 1
for (i in 1:100) {
  sim_fit[[i]] <- fit$obj$simulate(par=fit$obj$env$last.par.best, complete = TRUE) # TRUE if you want all the generated data
  
  res                    <- new("FLSAM")
  res@n.states           <- as.integer(sim_fit[[i]]$noYears)
  res@name               <- sao.name
  res@desc               <- paste("iter", i,  
                                  "FLSAM object generated:",
                                  date(), sep=" ")
  # Extract ages and plusgroup
  res@range["min"]       <- sim_fit[[i]]$minAge
  res@range["max"]       <- sim_fit[[i]]$maxAge
  res@range["plusgroup"] <- ifelse(sim_fit[[i]]$maxAgePlusGroup[1]==1,sim_fit[[i]]$maxAge,NA)
  
  # Extract years
  res@range["minyear"] <- min(sim_fit[[i]]$years)
  res@range["maxyear"] <- max(sim_fit[[i]]$years)
  
  # Extract the bindings
  res@states             <- sim_fit[[i]]$keyLogFsta   
  colnames(res@states)   <- res@range["min"]:res@range["max"]
  rownames(res@states)   <- c("catch",seq(nrow(res@states))[-1])
  
  # Extract the fbar ranges
  res@range["minfbar"] <- sim_fit[[i]]$fbarRange[1]
  res@range["maxfbar"] <- sim_fit[[i]]$fbarRange[2]
  
  # stock numbers
  res@stock.n          <- FLQuant(NA, dimnames=list(age=res@range["min"]:res@range["max"], 
                                                    year=res@range["minyear"]:res@range["maxyear"]))
  res@stock.n[,]       <- exp(sim_fit[[i]]$logN)
  
  # harvest
  n.ages               <- nrow(sim_fit[[i]]$logF)
  res@harvest          <- FLQuant(NA, dimnames=list(age=res@range["min"]:res@range["max"], 
                                                    year=res@range["minyear"]:res@range["maxyear"]))
  
  res@harvest[res@range["min"]:(res@range["min"]+n.ages),] <- exp(sim_fit[[i]]$logF)
  res@harvest[(n.ages+1),]    <- res@harvest[n.ages,]
  units(res@harvest)                   <-  "f"
  
  #Populate the info slot
  info                 <- data.frame(FLSAM.version   = packageDescription("FLSAM")$Version,
                                     FLCore.version = packageDescription("FLCore")$Version,
                                     R.version      = R.version$version.string,
                                     platform       = R.version$platform,
                                     run.date       = Sys.time())
  res@info <- t(info)
  colnames(res@info) <- "_"
  # summary(res)
  
  # CREATE FLStock, drop landings.fraction
  fls <- do.call("FLStock", fqs[-1])
  
  # CALCULATE landings.n and discards.n
  fls@landings.n <- fqs[["catch.n"]] * fqs[["landings.fraction"]]
  fls@discards.n <- fqs[["catch.n"]] * (1-fqs[["landings.fraction"]])
  
  # COMPUTE totals
  landings(fls) <- computeLandings(fls)
  discards(fls) <- computeDiscards(fls)
  catch(fls) <- computeCatch(fls)
  
  # add stock numbers and harvest from FLSAM
  stock.n(fls) <- res@stock.n
  harvest(fls) <- res@harvest
  
  # SET to standard units
  units(fls) <- standardUnits(fls)
  
  # GET fbar range
  fls@range["minfbar"] <- fit$conf$fbarRange[1]
  fls@range["maxfbar"] <- fit$conf$fbarRange[2]
  
  # SET name and desc
  name(fls) <- paste(sao.name, i)
  desc(fls) <- paste("Loaded from", url)
  
  flstocks[[i]] <- fls
  
  if (nrow(df) == 0) {
    df <- 
      as.data.frame(ssb(flstocks[[i]])) %>% 
      mutate(run=i)
  } else {
    df <- 
      df %>% 
      bind_rows(., (as.data.frame(ssb(flstocks[[i]])) %>% mutate(run=i)))
  }
}

df %>% 
  filter(run<=2) %>% 
  ggplot(aes(x=year,y=data, group=run)) +
  theme_bw() +
  geom_line(colour="gray")

# clean up data files; remove temp directory
# unlink(tempdir, recursive=TRUE)

plot(flstocks[[3]])
