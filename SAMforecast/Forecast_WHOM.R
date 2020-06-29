#################################################
## Example script to run SAM forecast for WHOM ##
## Author: Vanessa Trijoulet vtri@aqua.dtu.dk  ##
#################################################

# Initialization ---------
devtools::install_github("vtrijoulet/SAM/stockassessment", ref="master2") # Only to do once to install the right package

library(stockassessment)
library(Matrix)
library(gridExtra)
name <- "WHOM_2018" 
fit<- fitfromweb(name, character.only = TRUE)
source("SAMforecast/HS_fit_fixBlim.R") # calculate a segmented regression with recruitment pairs, modified from original function HS_fit.R available in the package, which estimates the inflexion point


# Define function arguments ---------
## Ref pts that shape the HCRs:
Blim=611814
MSYBtrig=856540
Fmsy=0.115 # Used as Ftarget in HCR
Flow=Fmsy/5 # To choose if HCR 3-5
## Simulation assumptions: 
SRpar <- HS_fit_fixBlim(name, pair.years=1995:2017) # years for SSB-R pairs, usually the same than ones used to calculate FMSY if relevant (SR=TRUE)
Ay<-max(fit$data$years)+(-9:0) # for average for bio input (M, w, etc.)
Ry<-1995:2017 # for average for recruitment input
Sy<-max(fit$data$years)+(-9:0) # for selectivity input
TAC2018=101672
TAC2019=110381
TAC2020=83954 
#TAC <- TAC2020 # For intermediate year

# To change between MS scenariosRW=TRUE
RW = FALSE # random walk on recruitment
Rdist=FALSE # recruitment with mean and sd from sampled years
F.RW=FALSE # if false no RW with increase variance in F in the forecast
SR=TRUE # for SR relationship for rec coming from SRpar
ny <- 20 # years for the forecast


nsim=5000
seed=12345


FC <- list()

# Run the forecast with HCR2 and SR relationship -----

set.seed(seed) # same seed to repeat output
FC[[length(FC)+1]] <- forecast2(fit, fscale=c(1, rep(NA,ny-1)), 
                                catchval=c(NA, TAC2018, TAC2019, TAC2020, rep(NA,ny-4)),
                                fval=,rep(NA, ny),
                                rec.years=Ry, 
                                ave.years=Ay,
                                overwriteSelYears=Sy, 
                                label="HCR2_SR",
                                nosim=nsim,deterministic=FALSE,
                                Fscenario=2,
                                MSYBtrig=c(rep(NA,4), rep(MSYBtrig,ny-4)), 
                                Blim=c(rep(NA,4), rep(Blim,ny-4)),
                                Fmsy=Fmsy,
                                Flow=Flow, #for Fscenario= 3, 4 or 5
                                RW=RW,
                                Rdist=Rdist,
                                SR=SR,
                                SRpar=SRpar,
                                F.RW=F.RW
)   

# Run the forecast with HCR2 with decreasing recruitment -----
RW = TRUE  # random walk on recruitment
drift = 0.2 # decrease in mean recruitment of 10% the sd
SR=FALSE # for SR relationship for rec coming from SRpar

set.seed(seed) # same seed to repeat output
FC[[length(FC)+1]] <- forecast2(fit, fscale=c(1, rep(NA,ny-1)), 
                                catchval=c(NA, TAC2018, TAC2019, TAC2020, rep(NA,ny-4)),
                                fval=,rep(NA, ny),
                                rec.years=Ry, 
                                ave.years=Ay,
                                overwriteSelYears=Sy, 
                                label="HCR2_RWdrift",
                                nosim=nsim,deterministic=FALSE,
                                Fscenario=2,
                                MSYBtrig=c(rep(NA,4), rep(MSYBtrig,ny-4)), 
                                Blim=c(rep(NA,4), rep(Blim,ny-4)),
                                Fmsy=Fmsy,
                                Flow=Flow, #for Fscenario= 3, 4 or 5
                                RW=RW,
                                drift=drift,
                                Rdist=Rdist,
                                SR=SR,
                                SRpar=SRpar,
                                F.RW=F.RW
)   


# Run the forecast with HCR2 with dSr but high catch to decrease SSB -----
RW = FALSE  # random walk on recruitment
SR=TRUE # for SR relationship for rec coming from SRpar

set.seed(seed) # same seed to repeat output
FC[[length(FC)+1]] <- forecast2(fit, fscale=c(1, rep(NA,ny-1)), 
                                catchval=c(NA, TAC2018*2, TAC2019*2, TAC2020*2, rep(NA,ny-4)),
                                fval=,rep(NA, ny),
                                rec.years=Ry, 
                                ave.years=Ay,
                                overwriteSelYears=Sy, 
                                label="HCR2_SR_lowSSB",
                                nosim=nsim,deterministic=FALSE,
                                Fscenario=2,
                                MSYBtrig=c(rep(NA,4), rep(MSYBtrig,ny-4)), 
                                Blim=c(rep(NA,4), rep(Blim,ny-4)),
                                Fmsy=Fmsy,
                                Flow=Flow, #for Fscenario= 3, 4 or 5
                                RW=RW,
                                Rdist=Rdist,
                                SR=SR,
                                SRpar=SRpar,
                                F.RW=F.RW
)   


