#setup

rm(list=ls())
gc()

# =================================
# Load libraries
library(MASS)
library(data.table) # I run it with r4ss version 1.35.1
library(Hmisc) # GL- required for checkBound
library(r4ss) # GL -but function could probably be altered not to require %nin%
library(matrixcalc)
library(dplyr)
library(FLCore)
library(ggplot2)
library(tidyr)

# =================================

# Set working directory of where the FUN_Random_population.R is stored as well as the folder where the assessment is
#mainDir <- "Y:\\FISHERIES M MoU\\Working_Area\\Scientific advice for fisheries management\\ICES_Working_groups\\2018\\WGWIDE"
#mainDir <- "Y:\\FISHERIES M MoU\\Working_Area\\Scientific advice for fisheries management\\ICES_Working_groups\\2018\\WGWIDE"
#mainDir <- "C:\\Stocks\\hom_27_2a4a5b6a7a-ce-k8\\MP_MSE\\MSE 2019\\SSForWHMMSE2019"
mainDir <- getwd()

# SOURCE FUNCTIONS
setwd(mainDir)
source('FUN_Random_population.R')

#inpDir <- file.path(mainDir, "SS_3.0.1\\FOR_ASSESSMENT\\BaseCase") # folder of where the assessment is located
#inpDir <- file.path(mainDir,"Base")

#base assessment
#Base <- "WGWIDE2018"
#Base <- "WGWIDE19"
Base <- "WGWIDE20"

#inpDir <- file.path(mainDir,"WGWIDE2018"); last.obs.year <- 2017
inpDir <- file.path(mainDir,Base); last.obs.year <- 2018    #2019 assessment
#inpDir <- file.path(mainDir,Base); last.obs.year <- 2019   #2020 assessment

# this is where the outputs of the simulations will be stored
dir.create(file.path(inpDir, "random_pop")) #GL - create that folder in here
outDir <- file.path(inpDir, "random_pop")

#MSE input files generated here
dir.create(file.path(outDir,"MSEInputs"))
MSEDir <- file.path(outDir,"MSEInputs")


