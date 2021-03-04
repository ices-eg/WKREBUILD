# RandomPop_CovarMatrix.R

# author: P. Carpi
# created: 19/04/2017
# modified:
# purpose: generate 1000 random population from variance-covariance matrix from SS results

# updated G Lambert 26/06/2019

source("0.Setup.R")

###################################
# =================================
## HOUSEKEEPING
# Get ADMB results
getFit <- read.admbFit(file.path(inpDir, "ss"))
# get variance-covariance matrix and estimates (will be my SD and mu for multivariate normal distribution)
dim(getFit$cov)
covMat <- getFit$cov
mu <- getFit$est

# then get the parameters with the bounds
mySS <- SS_output(dir=inpDir,covar=T, verbose=F, forecast=TRUE)

#SS_plots(replist=mySS)

myParBounds <- mySS$parameters
myParNames <- read.csv(file.path(inpDir, "ss.par"), sep=" ")
# get names from par file
parNames <- clean_names(myParNames)
# and create new variable (to have those matching the names I get when I use read.admbFit)
parBound <- cbind(myParBounds[,c("Value", "Phase", "Min", "Max", "Init")], parNames[parNames!="checksum999"])
# =================================
###################################

# Draw my multivariate normal distribution (more that 1000 repetition because you might have parameters that are outside bounds, and if
# so I will remove the whole line)
set.seed(1)

#MVND <- data.table(mvrnorm(n = 100000, mu = mu, Sigma = covMat, empirical = TRUE))
MVND <- data.table(mvrnorm(n = 1000000, mu = mu, Sigma = covMat, empirical = TRUE))

# Assign names
names(MVND) <- getFit$names
# transpose it
newPar <- t(MVND)

# Create new dataframe with everything within bounds
newParSet <- checkBound(bounds=parBound, new_par = newPar)
newParSet <- newParSet[,-1]
# get only 1000
#newParSet <- newParSet[,1:1000]
newParSet <- newParSet[,1:10000]
newParSet <- cbind(newParSet, row.names(newPar)[row.names(newPar) %in%  parBound$parNames]) #GL there was a problem here - I replaced new_par with newPar
names(newParSet)[ncol(newParSet)] <- "par"

#write.csv(newParSet, file=paste(outDir,"Random_parameters_1000it.csv",sep="/"), row.names=F, quote = F)
#write.csv(newParSet, file=paste(outDir,"Random_parameters_1000it.csv",sep="/"), row.names=F, quote = F)
write.csv(newParSet, file=paste(outDir,"Random_parameters_10000it.csv",sep="/"), row.names=F, quote = F)

# Now run SS 
## GL Notes below
# This will take a long time to run and create all the model replicates, taking a lot of disc space
## 1) I'd suggest to hack into this function in the "FUN_Random_population.R" file to extract what is needed and delete, 
# or overwrite, the folders from one replicate to the next
## 2) Make sure to keep the files required to run "SSgetouput" below
# if getcovar is set to FALSE (as it is now) then no need for covar.sso (this will exist only if Hessian is run anyway - see nex comment)
# if getComp is set to FALSE (as it is now) then no need for CompReport.sso 
# I think the forecast file is also only needed if forecast = T
# so probably the report.sso file (and the warning file probably?) are needed. 
# This is to be checked before running all the sims to make sure nothing is missing and have to run everything again..
## 3) Also, as it is set up it will run the Hessian for each replicate, I don't think this is required, is it?
# It would be much faster to run without - to do so, go to "FUN_Random_population.R", go to FUN5 "SS_doPar" and 
# replace the argument extras = "-nox" by extras = "-nox -nohess" - that should do the trick?
#SS_doPar(parNames=parNames, parBound=parBound, myParSet=newParSet, n_parSets=1:1000)
SS_doPar(parNames=parNames, parBound=parBound, myParSet=newParSet, n_parSets=1:10000)

# Compare the results from the model.
# Get the outputs to be compared
allFiles <- paste(outDir, "Set_newParameter", dir(file.path(outDir, "Set_newParameter")), sep="/")

#every 10th for a comparison
randomPop <- SSgetoutput(dirvec = allFiles[seq(1,length(allFiles),10)],getcomp=F, getcovar=F, forecast=F)

# Create summary comparison
summary_randomPop <- SSsummarize(randomPop)

# Plot comparison
SSplotComparisons(summary_randomPop, indexfleets = 2, indexUncertainty = FALSE, print=TRUE, png=TRUE, pdf=FALSE,
                  plotdir = file.path(outDir, "Set_newParameter/Plots"), legend = FALSE)






