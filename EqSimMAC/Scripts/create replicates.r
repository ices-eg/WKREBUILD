##-------------------------------------------------------------------------------
# WKMSEMAC
#
# Author: Thomas BRUNEL
#         IMARES, The Netherland
#
#  MSE of NEA mackerel
#
# Date: 2020...
#
# Build for R3.5.3, 64bits
#-------------------------------------------------------------------------------

# this function carries out the resampling from the model fit variance covariance matrix 
# and store the outcome in two objects:
# - the  parameters (fixed effects) stored as a data.frame
# - the states (N and F at age) stored in an FLStock object with nits dimensions

create.replicates <- function(stck,samfit,n)
{
# stck  <- Mac ; n <- 3 ; samfit  <- fit

# define age and year ranges in the assessment used
ages <- an(dimnames(stck)$age)
years<- an(dimnames(stck)$year)                         
                          
# 1 - resampling-------------------------------------
# compute the varcov matrix from samfit object                          
sdrep <- TMB::sdreport(samfit$obj,getJointPrecision=T)
sigma <- as.matrix(solve(sdrep$jointPrecision))
mu    <- c(sdrep$par.fixed,sdrep$par.random)
# resample fixed and random effects values
sim   <- mvrnorm(n,mu=mu,Sigma=sigma)

sim.org <- sim
            colnames(sim) <- rownames(sigma)
            rownames(sim) <- 1:n
sim <-as.data.frame(sim)


#names(sim) <-paste0(names(sim),1:dim(sim)[2])




# 2 - reshape random parameters-----------------------
fixed  <- samfit$sdrep$par.fixed
random <- samfit$sdrep$par.random
names(fixed) <-paste0(names(fixed),1:length(fixed))
names(random) <-paste0(names(random),length(fixed)+c(1:length(random)))
fixed<- data.frame(par=names(fixed),est.value=fixed)

sim.fixed <- sim[,1:length(samfit$sdrep$par.fixed)]
# add and index to the name of parameters so that they match with the sam control object
nams <-  names(sim)[1:length(samfit$sdrep$par.fixed)]
ix <- rep(1,length(nams))
for  (i in 2:length(nams))  if(nams[i]==nams[i-1])   ix[i] <- ix[i-1]+1 
names(sim.fixed)  <- paste0(nams,ix)

sim.fixed$iter  <- rownames(sim.fixed)
sim.fixed <- tidyr::gather(sim.fixed,"par","sim.value",1:c(dim(sim.fixed)[2]-1))




# 2 - reshape simulated states -----------------------     
sim.random<- sim[,c(1+dim(fixed)[1]):dim(sim)[2]]
sim.random$iter <- rownames(sim.random)

# workout the N@age values
N <- sim.random[,grep("logN",names(sim.random))]
N<-as.data.frame(t(N))
N$age <-rep(ages,length(years))
N$year<-sort(rep(years,length(ages)))
N<- tidyr::gather(N,"iter","sim.value",1:n)
N$sim.value <-exp(N$sim.value)
Ns<-array(N$sim.value, dim=c(length(ages),length(years),1,1,1,n))
# workout the F@age values     
F. <- sim.random[,grep("logF",names(sim.random))]
F.<-as.data.frame(t(F.))
F.$age <-rep(0:7,length(years))                  # Warning : ideally the age for selectivity plateau should not be hard coded, but read from the conf file (if we wanted to use this code for another stock)
F.$year<-sort(rep(years,8))
F.<- tidyr::gather(F.,"iter","sim.value",1:n)
F.$sim.value <-exp(F.$sim.value)
Fs<-array(F.$sim.value,dim=c(8,length(years),1,1,1,n))


# 3- incorporate states into FLStock object
stcks <- propagate(stck , iter = n)
stock.n(stcks)           <- FLQuant(Ns,dimnames=dimnames(stock.n(stcks)) )
harvest(stcks)[ac(0:7),] <- FLQuant(Fs,dimnames=dimnames(harvest(stcks)[ac(0:7),]) )
harvest(stcks)[ac(8:12),][] <- harvest(stcks)[ac(7),]


# return results
res <- list(random.params= sim.fixed, stocks = stcks, sim = sim.org)
return(res)

}
