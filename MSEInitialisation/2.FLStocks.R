#REMEMBER - set the base assessment in setup.R
source("0.Setup.R")
source("SS3toFLStock_v3.30_AC.r")

cat(Base,"\n")

iters <- seq(1,10000)   #10k assessments from the mvrn draw
k.iters <- seq(1,1000); nits <- length(k.iters)   #number of iters for the simulation
ages <- seq(0,15); nages<-length(ages)   #ages

#folders containing SS output files Report.sso, CompReport.sso and covar.sso
ss.repDirs <- c(file.path(getwd(),Base),
                file.path(getwd(),Base,"random_pop","Set_newParameter",paste0("newParms",iters)))
names(ss.repDirs) <- c("ass",iters)

#files moved to external drive
ss.repDirs[2:10001] <- file.path("E:","Stocks","hom_27_2a4a5b6a7a-ce-k8","Assessment","MSEInitialisation",Base,"random_pop","Set_newParameter",paste0("newParms",iters))

# Load RData function
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#assessment output
WHOM <- SS3toFLStock(ss.repDirs[1],stockname="WHOM",fbar_amin=1,fbar_amax=10)
WHOM <- FLCore::setPlusGroup(WHOM,15)
FLS <- WHOM
dim(FLS)
FLCore::plot(FLS,main=Base)

#assessment years
ass.yrs <- seq(dims(FLS)$minyear,dims(FLS)$maxyear)
nass.yrs <- length(ass.yrs)
#terminal assessment year
tyr <- dims(FLS)$maxyear

#original list object and selected iterations (used in initial simulations)
#create a new FLStock-based object & save
#load(file=file.path(getwd(),Base,"RData","MSE_WGWIDE19_FLStocks_10k15PG_V1.RData"))
#source(file=file.path(getwd(),"origIters.R"))
#create FLStock with 1000 iterations
#FLSs.1k <- FLCore::propagate(FLS,1000)
#for (ii in seq(1,1000)){FLSs.1k[,,,,,ii]<-lWHM[[orig.iters[ii]]]}
#save(FLSs.1k,file = file.path(getwd(),Base,paste0("WHOM_SS",substr(Base,7,8),"_FLS_orig.RData")))


# #1000 iterations at a time (perf)
# tFLSs <- FLCore::propagate(FLS,1000)
# 
# for (dec in 1:10){
#   for (ii in 1:1000){
#     WHOM <- SS3toFLStock(ss.repDirs[1000*(dec-1)+ii+1],stockname="WHOM",fbar_amin=1,fbar_amax=10)
#     WHOM <- FLCore::setPlusGroup(WHOM,15)
#     tFLSs[,,,,,ii] <- WHOM
# }
#   save(tFLSs,file=file.path(getwd(),Base,paste0("1000_",dec,"_Stocks.RData")))
# }
# 
# #10,000 in a single FLStock
# FLSs <- FLCore::propagate(FLS,10000)
# for(i in 1:10){FLSs[,,,,,(1000*(i-1)+1):(1000*i)] <- loadRData(file=file.path(getwd(),Base,paste0("1000_",i,"_Stocks.RData")))}
# save(FLSs,file = file.path(getwd(),Base,paste0("MSE_",Base,"FLStocks_10k15PG.RData")))

#or read in if already done
load(file = file.path(getwd(),Base,"RData",paste0("MSE_",Base,"_FLStocks_10k15PG.RData")))
dim(FLSs)


#select 1000 iterations for the MSE#########################################################################

if (Base=="WGWIDE19") {
  #WGWIDE 2019 (Low-High is 95%)
  WGWIDE.tSSB <- 811685
  WGWIDE.tSSB.Lo <- 562585
  WGWIDE.tSSB.Hi <- 1060785
  WGWIDE.tSSB.SD <- 0.5*(WGWIDE.tSSB.Hi-WGWIDE.tSSB)
}

if (Base=="WGWIDE20") {
  #WGWIDE 2020
  WGWIDE.tSSB <- 808972
  WGWIDE.tSSB.Lo <- 537242
  WGWIDE.tSSB.Hi <- 1080703
  WGWIDE.tSSB.SD <- 0.5*(WGWIDE.tSSB.Hi-WGWIDE.tSSB)
}

#assessment output terminal year
#terminal year SSB
dfSSB <- data.frame(iter=iters,SSB=as.numeric(ssb(FLSs[,as.character(tyr),,,,])))
quantile(dfSSB$SSB,probs=c(0.025,0.5,0.975))

#WGWIDE 2019
#2.5%       50%     97.5% 
#575956.9  841177.3 1146767.2 
#WGWIDE 2020
#2.5%        50%      97.5% 
#18547.82  887926.70 1481933.41

#Ensure that the median and 95th percentiles of the 1000 starting pops matches the estimates from the appropriate WG assessment
#10,000 iterations don't match this
#median is slightly higher, tails longer

#select from 10k stocks a subset of 1000 with terminal year SSB distribution ensemble statistics matching the assessment
sel.iters <- c()

#dist of SSBs
set.seed(1)
bks <- c(0,seq(400000,1300000,by=100000),1e10)
bins <- hist(rnorm(1000,mean=WGWIDE.tSSB,sd=WGWIDE.tSSB.SD),breaks=bks,plot=FALSE)
for(bin in seq(1,length(bks)-1)){
  if(bins$counts[bin]>0){
    sel.iters <- c(sel.iters,sample(dfSSB$iter[dfSSB$SSB>bks[bin] & dfSSB$SSB<=bks[bin+1]],size=bins$counts[bin],replace=FALSE))
  }
}

#check there are 1000 selected
length(sel.iters)
sum(is.na(sel.iters))

#save(sel.iters,file=file.path(getwd(),Base,"RData",paste0("IterationNos_",Base,".RData")))
load(file=file.path(getwd(),Base,"RData",paste0("IterationNos_",Base,".RData")))

png(filename = file.path(getwd(),"..","EqSimWHM","ReportGraphics",paste0("InitialSSB_",Base,".png")), width = 600, height = 600)
h <- hist(dfSSB$SSB[sel.iters]/1e6,breaks=12,main=paste0("Distribution of initial SSB (",Base,")"), xlab=paste0("SSB ",tyr," (Mt)"),axes=F,ylab="")
axis(side=1)
abline(v=as.numeric(ssb(FLS)[,ac(tyr)])/1e6,lwd=2,lty=2)
abline(v=c(WGWIDE.tSSB.Lo,WGWIDE.tSSB.Hi)/1e6,lwd=2,lty=2)
abline(v=quantile(dfSSB$SSB[sel.iters],probs=c(0.025,0.5,0.975))/1e6,lwd=2,lty=2,col="red")
quantile(dfSSB$SSB[sel.iters],probs=c(0.025,0.5,0.975))

x <- seq(0,2,by=0.01)
y <- dnorm(x, mean = WGWIDE.tSSB/1e6, sd = WGWIDE.tSSB.SD/1e6)
y <- max(h$counts)*y/max(y)
lines(x,y)
dev.off()


#WGWIDE2019
#2.5%       50%     97.5% 
#562220.7  810881.9 1065140.1 

#WGWIDE2020
#2.5%       50%     97.5% 
#522507.1  807792.5 1092449.5 

#select the 1000 from the 10k & save
FLSs.1k <- FLSs[,,,,,sel.iters]

#replace first iteration with the assessment output
FLSs.1k[,,,,,1] <- FLS

plot(FLSs.1k,main=Base)
#check first iter (should match assessment)
ssb(FLSs.1k[,,,,,1])/ssb(FLS)
fbar(FLSs.1k[,,,,,1])/fbar(FLS)
rec(FLSs.1k[,,,,,1])/rec(FLS)

#intermediate save of 1000 iterations
#save(FLSs.1k,file = file.path(getwd(),Base,"RData",paste0("WHOM_SS",substr(Base,7,8),"_FLS.RData"))) #old name
#save(FLSs.1k,file = file.path(getwd(),Base,"RData",paste0("WHOM_SS",substr(Base,7,8),"_FLS_Clean.RData")))
load(file = file.path(getwd(),Base,"RData",paste0("WHOM_SS",substr(Base,7,8),"_FLS_Clean.RData")))

#plot of selected iters 
dfSelectedIters <- data.frame(Iter=c(),Year=c(),SSB=c(),FBar=c(),Rec=c())
for (iter in c(1,98,198,298,398,498,598,698,798,898,998)){
  dfSelectedIters <- dplyr::bind_rows(dfSelectedIters,
                                      data.frame(Iter=iter,Year=ass.yrs,
                                                 SSB=as.numeric(ssb(FLSs.1k[,,,,,iter]))/1e6,
                                                 FBar=as.numeric(fbar(FLSs.1k[,,,,,iter])),
                                                 Rec=as.numeric(rec(FLSs.1k[,,,,,iter]))/1e6))
}

dfSelectedIters$fIter <- as.factor(dfSelectedIters$Iter)
png(filename = file.path(getwd(),"..","EqSimWHM","ReportGraphics",paste0("SelectSSB_Historical",Base,".png")), width = 600, height = 400)
ggplot(data = filter(dfSelectedIters,Iter>1), mapping = aes(x=Year, y=SSB, group=fIter, col=fIter)) + 
  geom_line(lwd=1) + theme(legend.position = "none") + 
  geom_line(data = filter(dfSelectedIters,Iter==1), mapping = aes(x=Year, y=SSB), colour="black", lty=2, lwd=1) +
  ylab("SSB (Mt)")
dev.off()
png(filename = file.path(getwd(),"..","EqSimWHM","ReportGraphics",paste0("SelectFBar_Historical",Base,".png")), width = 600, height = 400)
ggplot(data = filter(dfSelectedIters,Iter>1), mapping = aes(x=Year, y=FBar, group=fIter, col=fIter)) + 
  geom_line(lwd=1) + theme(legend.position = "none") + 
  geom_line(data = filter(dfSelectedIters,Iter==1), mapping = aes(x=Year, y=FBar), colour="black", lty=2, lwd=1) +
  ylab("FBar (1_10)")
dev.off()
png(filename = file.path(getwd(),"..","EqSimWHM","ReportGraphics",paste0("SelectRecr_Historical",Base,".png")), width = 600, height = 400)
ggplot(data = filter(dfSelectedIters,Iter>1), mapping = aes(x=Year, y=Rec, group=fIter, col=fIter)) + 
  geom_line(lwd=1) + theme(legend.position = "none") + 
  geom_line(data = filter(dfSelectedIters,Iter==1), mapping = aes(x=Year, y=Rec), colour="black", lty=2, lwd=1) +
  ylab("Recruitment (billions)")
dev.off()


#qqplot of terminal year SSB
ggplot(data = data.frame(SSB=as.numeric(ssb(FLSs.1k[,as.character(tyr)]))), aes(sample=SSB)) + 
  stat_qq() + stat_qq_line() + ggtitle(Base)

#distribution of abundances@age in starting year
png(filename = file.path(getwd(),"..","EqSimWHM","ReportGraphics",paste0("InitialAbd_",Base,".png")), width = 600, height = 600)
ggplot(data = data.frame(data.frame(age = as.numeric(rep(dimnames(FLSs.1k)$age,dim(FLSs.1k)[6])),
                                    abd = as.numeric(FLCore::stock.n(FLSs.1k[,as.character(tyr)]))/1e6)),
       aes(abd)) +
  geom_histogram() +
  geom_vline(data = data.frame(age=seq(0,15),abd=as.numeric(stock.n(FLS[,as.character(tyr)]))/1e6), 
             aes(xintercept = abd), colour="red", lwd=1) +
  facet_wrap(~age, labeller = labeller(age = c("0" = "Age 0", "1" = "Age 1", "2" = "Age 2", "3" = "Age 3",
                                               "4" = "Age 4", "5" = "Age 5", "6" = "Age 6", "7" = "Age 7",
                                               "8" = "Age 8", "9" = "Age 9", "10" = "Age 10", "11" = "Age 11",
                                               "12" = "Age 12", "13" = "Age 13", "14" = "Age 14", "15" = "Age 15+")), scales="free") +
  ggtitle(Base) + theme(axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank()) +
  ggtitle(paste0("Distributions of initial abundance at age (",Base,")"))
dev.off()



#create a second dataset with a larger SSB SD, assessment often considered to underestimate uncertainty

sel.iters <- c()
set.seed(1)
bks <- c(0,seq(400000,1300000,by=100000),1e10)
bins <- hist(rnorm(1000,mean=WGWIDE.tSSB,sd=1.5*WGWIDE.tSSB.SD),breaks=bks,plot=FALSE)
for(bin in seq(1,length(bks)-1)){
  if(bins$counts[bin]>0){
    sel.iters <- c(sel.iters,sample(dfSSB$iter[dfSSB$SSB>bks[bin] & dfSSB$SSB<=bks[bin+1]],size=bins$counts[bin],replace=TRUE))
  }
}

#check there are 1000 selected
length(sel.iters)
sum(is.na(sel.iters))

hist(dfSSB$SSB[sel.iters],main=Base)
quantile(dfSSB$SSB[sel.iters],probs=c(0.025,0.5,0.975))

#WGWIDE2019
#2.5%       50%     97.5% 
#431393.5  808996.0 1179805.7 

#WGWIDE2020
#2.5%       50%     97.5% 
#312586.1  804469.2 1227604.5 

FLSs.1k.1.5SD <- FLSs[,,,,,sel.iters]
#replace first iteration with the assessment output
FLSs.1k.1.5SD[,,,,,1] <- FLS

plot(FLSs.1k.1.5SD,main=Base)
ssb(FLSs.1k.1.5SD[,,,,,1])/ssb(FLS)
fbar(FLSs.1k.1.5SD[,,,,,1])/fbar(FLS)
rec(FLSs.1k.1.5SD[,,,,,1])/rec(FLS)

#save the inflated noise dataset
#save(FLSs.1k.1.5SD,file = file.path(getwd(),Base,"RData",paste0("WHOM_SS",substr(Base,7,8),"_FLS_15SD_Clean.RData")))
load(file = file.path(getwd(),Base,"RData",paste0("WHOM_SS",substr(Base,7,8),"_FLS_15SD_Clean.RData")))

#distribution of abundances@age in starting year
ggplot(data = data.frame(data.frame(age = as.numeric(rep(dimnames(FLSs.1k.1.5SD)$age,dim(FLSs.1k.1.5SD)[6])),
                                    abd = as.numeric(FLCore::stock.n(FLSs.1k.1.5SD[,as.character(tyr)]))/1e6)),
       aes(abd)) +
  geom_histogram() +
  geom_vline(data = data.frame(age=seq(0,15),abd=as.numeric(stock.n(FLS[,as.character(tyr)]))/1e6), 
             aes(xintercept = abd), colour="red", lwd=1) +
  facet_wrap(~age, scales="free") +
  ggtitle(Base)


#####SELECTION#######################################################

#SS3 uses fixed selection at age
#EqSim is designed to resample at random from a predefined number of years
#1000 separate selection patterns are available from FLSs.1k

fGetFLStockSelection <- function(stk){
  sel <- matrix(FLCore::harvest(stk), ncol = dim(stk)[2])
  Fbar <- matrix(FLCore::fbar(stk), ncol = dim(stk)[2])
  sel <- sweep(sel, 2, Fbar, "/")
  sel <- sel/max(sel[,seq(dim(sel)[2]-10,dim(sel)[2])])  #last 10 years to avoid noise at start
}

ac <- function(x){as.character(x)}

#populate years prior to the terminal year with a randomly selected selection pattern

#extract all selection patterns into a single array
sels <- array(NA,dim=c(nages,nits),dimnames=list(age=ages,iter=sel.iters))
for(i in 1:nits){
  tsel <- fGetFLStockSelection(FLSs.1k[,,,,,i])
  #take the last year (all will be the same)
  sels[,i] <- tsel[,ncol(tsel)]
}


dfSels <- data.frame(Age=ages,Selectivity=c(sels),Iter=rep(seq(1,nits),each=length(ages)))
dfSels$fIter <- as.factor(dfSels$Iter)
png(filename = file.path(getwd(),"..","EqSimWHM","ReportGraphics",paste0("Selection",Base,".png")), width = 600, height = 400)
ggplot(data = filter(dfSels,Iter %in% c(98,198,298,309,498,598,698,798,898,998)),
       mapping = aes(x=Age, y=Selectivity, group=fIter, col=fIter)) + geom_line(lwd=1) +
  theme(legend.position = "none") + 
  geom_line(data = filter(dfSels,Iter==1), mapping = aes(x=Age, y=Selectivity), colour="black", lty=2, lwd=1)
dev.off()

png(filename = file.path(getwd(),"..","EqSimWHM","ReportGraphics",paste0("SelectionByAge",Base,".png")), width = 600, height = 400)
ggplot(data = dfSels, mapping = aes(Selectivity)) + 
  geom_bar() + scale_x_binned(n.breaks=15, nice.breaks=FALSE) +
  geom_vline(xintercept=c(0.25,0.5,0.75), col="grey", lwd=0.5, lty=2) +
  facet_wrap(~Age, labeller = labeller(Age = c("0" = "Age 0", "1" = "Age 1", "2" = "Age 2", "3" = "Age 3",
                                              "4" = "Age 4", "5" = "Age 5", "6" = "Age 6", "7" = "Age 7",
                                              "8" = "Age 8", "9" = "Age 9", "10" = "Age 10", "11" = "Age 11",
                                              "12" = "Age 12", "13" = "Age 13", "14" = "Age 14", "15" = "Age 15+"))) +
  theme(axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank())
dev.off()

#quick look at range of patterns
plot(ages,seq(0,length=nages),type="n",ylim=c(0,1))
for (i in k.iters){lines(ages,sels[,i])}
lines(ages,fGetFLStockSelection(FLS)[,length(ass.yrs)],col="red",lty=2)

#randomly assign a selection pattern to each of the last 10 years for each stock iteration
#this will mess up the historic aspect of the stock dev but should not impact EqSim execution
yrs <- seq(length(ass.yrs)-9,length(ass.yrs)-1)
set.seed(1)
ransel <- array(data = sample(seq(1,nits),length(yrs)*length(k.iters),replace=TRUE),
                dim=c(length(yrs),length(k.iters)), dimnames=list(yr=yrs,iter=k.iters))

#CAUTION - this takes approx 1hr
for (y in yrs){
  cat(y,"\n")
  for (i in 2:nits) {   #leave first slot (assessment)
      #FLCore::harvest(FLSs.1k[,y,,,,i]) <- sels[,ransel[as.character(y),as.character(i)]]
      FLCore::harvest(FLSs.1k[,y,,,,i]) <- sels[,ransel[as.character(y),as.character(i)]]*(as.numeric(FLCore::fbar(FLSs.1k[,y,,,,i]))/mean(sels[2:11,ransel[as.character(y),as.character(i)]]))
  }
}

#check of selection patterns - 10 random iterations
plot(ages,seq(0,length=nages),type="n",ylim=c(0,1),ylab="Selection",xlab="Age")
ss <- sample(k.iters,10)
for (s in ss){
  for (y in yrs) {
    lines(ages,as.numeric(FLCore::harvest(FLSs.1k[,y,,,,s])),col=topo.colors(n=10)[which(ss==s)])
  }
}
lines(ages,fGetFLStockSelection(FLS)[,length(ass.yrs)],col="black",lty=2,lwd=2)

#selection at age
#data, assessment point estimates and range of values from 1000 iters
dfSelections <- data.frame(Age = ages, Iter = rep(k.iters,each=nages), Sel = as.numeric(sels))
dfSelections <- dplyr::bind_rows(dfSelections,
                                 data.frame(Age = ages, Iter=rep(0,nages), Sel=as.numeric(fGetFLStockSelection(FLS)[,length(ass.yrs)])))
gSelProfs <- ggplot(data = dfSelections) + 
  geom_line(aes(x=Age,y=Sel,group=Iter)) +
  geom_line(filter(dfSelections,Iter==0),mapping=aes(x=Age,y=Sel,group=1),col="red",lwd=1)

#distributions of selection@age
#distribution of values from 1000 iterations, point estimates from assessment
#sort stupid long x aixs labels
gSelHistatAge <- ggplot(data = dfSelections, mapping = aes(Sel)) + 
  geom_bar() + scale_x_binned(n.breaks=15, nice.breaks=FALSE) +
  theme(text = element_text(size=10), axis.text.x = element_text(angle=45, hjust=1)) +
  geom_vline(xintercept=c(0.25,0.5,0.75), col="grey", lwd=1, lty=2) +
  facet_wrap(~Age) + ylab("Count") +
  theme(axis.text.x=element_blank())

#save(FLSs.1k,file = file.path(getwd(),Base,"RData",paste0("WHOM_SS",substr(Base,7,8),"_FLS_varSel.RData")))
load(file = file.path(getwd(),Base,"RData",paste0("WHOM_SS",substr(Base,7,8),"_FLS_varSel.RData")))

###############catch and stock weight ####################################################
#SS3 fits a time invariant model to catch/stock weight data
#as with the selection, insert historic values into the FLStock slots so that EqSim methodology is applicable
#this is not as complex as the selectivity as we can simply use the historic data for each year slot
#only fill the most recent 10 years - the earlier ones can always be used to retrieve the assessment fit

#read in the historic catch datasets
dfSW <- readr::read_delim(file = file.path(getwd(),Base,"StockWeights.dat"),delim=",")
dfCW <- readr::read_delim(file = file.path(getwd(),Base,"CatchWeights.dat"),delim=",")

#single df
dfWeights <- dplyr::left_join(
  dfSW %>% pivot_longer(cols = paste0("Age",seq(0,15)), names_to = "Age", values_to = "SW", names_prefix = "Age"),
  dfCW %>% pivot_longer(cols = paste0("Age",seq(0,15)), names_to = "Age", values_to = "CW", names_prefix = "Age"),
  by=c("Year","Age"))

#cv, phi based on data since 2000 - prior to this date no variability in SW estimates for youngest ages
dfWSummary <- dfWeights %>% filter(Year>=2000) %>% group_by(Age) %>% summarise(SW.cv = sd(SW)/mean(SW), 
                                          SW.phi1 = acf(SW,plot=FALSE)[[1]][2],
                                          CW.cv = sd(CW)/mean(CW),
                                          CW.phi1 = acf(CW,plot=FALSE)[[1]][2])

dfWSummary <- dfWSummary %>% mutate(CV=(SW.cv+CW.cv)/2,Phi=(SW.phi1+CW.phi1)/2)

plot(dfWSummary$Age,dfWSummary$SW.cv,ylim=c(0,0.5),ylab="CV",xlab="Age")
points(dfWSummary$Age,dfWSummary$CW.cv,pch=19)
points(dfWSummary$Age,dfWSummary$CV,pch=19,col="red")

plot(dfWSummary$Age,dfWSummary$SW.phi1,ylim=c(0,0.75),ylab="Phi",xlab="Age")
points(dfWSummary$Age,dfWSummary$CW.phi1,pch=19)
points(dfWSummary$Age,dfWSummary$Phi,pch=19,col="red")

#error generating example
#set.seed(1)
#cv <- 0.5
#phi <- 0.05
#err <- array(NA,dim=c(1000,1000),dimnames=list(iter=seq(1,1000),yr=seq(1,1000)))
#err[,1] <- stats::rnorm(n=1000, mean=0, sd=1)*cv/sqrt(1-phi^2)
#for (i in 2:1000){err[,i] <- phi*err[,i-1]+cv*stats::rnorm(n = 1000, mean = 0, sd = 1)}
#hist(err)
#dat1 <- array(1,dim=c(1000,1000),dimnames=list(iter=seq(1,1000),yr=seq(1,1000)))
#edat1 <- dat1*exp(err)
#cvs <- apply(edat1,MARGIN=1,FUN=function(x){sd(x)/mean(x)})
#quantile(cvs,probs=c(0.05,0.5,0.95))
#hist(cvs,xlab="cv",main="CV")
#phis <- apply(edat1,MARGIN=1,FUN=function(x){acf(x,plot=FALSE)[[1]][2]})
#quantile(phis,probs=c(0.05,0.5,0.95))
#hist(phis,xlab="phi",main="Phi")

#array to store error (100 years generated to test stats)
Werr <- SWerr <- CWerr <- FLQuant(NA,dimnames=list(age=ages,year=seq(1,100),iter=seq(1,1000)))
set.seed(1)

#generate errors
for (a in ages){
  #first year
  SWerr[as.character(a),1,] <- stats::rnorm(n=1000, mean=0, sd=1)*dfWSummary$SW.cv[dfWSummary$Age==a]/sqrt(1-dfWSummary$SW.phi1[dfWSummary$Age==a]^2)
  CWerr[as.character(a),1,] <- stats::rnorm(n=1000, mean=0, sd=1)*dfWSummary$CW.cv[dfWSummary$Age==a]/sqrt(1-dfWSummary$CW.phi1[dfWSummary$Age==a]^2)
  Werr[as.character(a),1,] <- stats::rnorm(n=1000, mean=0, sd=1)*dfWSummary$CV[dfWSummary$Age==a]/sqrt(1-dfWSummary$Phi[dfWSummary$Age==a]^2)
  #following years
  for(j in 2:100){
    SWerr[as.character(a),j,] <- dfWSummary$SW.phi1[dfWSummary$Age==a] * SWerr[as.character(a),j-1,] + 
      dfWSummary$SW.cv[dfWSummary$Age==a] * stats::rnorm(n = 1000, mean = 0, sd = 1)
    CWerr[as.character(a),j,] <- dfWSummary$CW.phi1[dfWSummary$Age==a] * CWerr[as.character(a),j-1,] + 
      dfWSummary$CW.cv[dfWSummary$Age==a] * stats::rnorm(n = 1000, mean = 0, sd = 1)
    Werr[as.character(a),j,] <- dfWSummary$Phi[dfWSummary$Age==a] * Werr[as.character(a),j-1,] + 
      dfWSummary$CV[dfWSummary$Age==a] * stats::rnorm(n = 1000, mean = 0, sd = 1)
  }
}

save(SWerr,CWerr,Werr,file=file.path(getwd(),Base,"RData","WErr_V3.RData"))
#load(file=file.path(getwd(),Base,"RData","WErr.RData"))

#select final 10 years
SWerr <- SWerr[,91:100]
CWerr <- CWerr[,91:100]
Werr <- Werr[,91:100]

#and rename the year dimension
dimnames(SWerr)$year <- as.character(seq(tyr-10,tyr-1))
dimnames(CWerr)$year <- as.character(seq(tyr-10,tyr-1))
dimnames(Werr)$year <- as.character(seq(tyr-10,tyr-1))

#check new obj
dim(SWerr)
dim(CWerr)

#apply the error - not to the first iteration which holds the assessment output - apply same error to both catch and stock weight
stock.wt(FLSs.1k[,as.character(seq(tyr-10,tyr-1)),,,,2:1000]) <- stock.wt(FLSs.1k[,as.character(seq(tyr-10,tyr-1)),,,,2:1000])*exp(Werr[,,,,,2:1000])
catch.wt(FLSs.1k[,as.character(seq(tyr-10,tyr-1)),,,,2:1000]) <- catch.wt(FLSs.1k[,as.character(seq(tyr-10,tyr-1)),,,,2:1000])*exp(Werr[,,,,,2:1000])

#save(FLSs.1k,file = file.path(getwd(),Base,paste0("MSE_",Base,"RData","_FLStocks_1k15PG_varSelection_varWgt.RData")))
save(FLSs.1k,file = file.path(getwd(),Base,"RData",paste0("WHOM_SS",substr(Base,7,8),"_FLS_varSelvarWgt_V3.RData")))
#load(file=file.path(getwd(),Base,"RData",paste0("WHOM_SS",substr(Base,7,8),"_FLS_varSelvarWgt.RData")))

#old code
#create the FLStock objects
#lWHM <- lapply(ss.repDirs,SS3toFLStock,stockname="WHOM",fbar_amin=1,fbar_amax=10)
#save(lWHM,file = file.path(MSEDir,paste0("MSE_",Base,"_FLStocks_10k.RData")))