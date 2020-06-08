#eqsim_compare_runs

fCompare_runs <- function(runs2Compare, Res.dir, Plot.dir, PerfStat, TargetFs, simYears, lStatPer, Blim){

  #produces a grid (max 6)
  #comparing the supplied performance statistic for ST,MT and LT for each runName
  #grouped by runName
  
  require(ggplot2)
  require(dplyr)
  
  dfAll <- data.frame(RunRef = c(), Ftgt = c(), PerfStat = c(), Period = c(), Val = c(), stringsAsFactors = FALSE)
  
  #get the data
  for (r in runs2Compare){
    
    fl <- file.path(Res.dir,r,paste0(r,"_SimRuns.RData"))
    cat(fl,"\n")
    load(fl)
    
    for (ftgt in TargetFs){
      
      t <- SimRuns[[ac(ftgt)]]
      
      #browser()
      
      if (PerfStat=="Catch") {
        #catch numbers
        Cnum <- t[["C"]]
        #catch weights
        Cwgt <- t[["catW"]]
        #catch weight (tons)
        CW <- apply(Cnum*Cwgt,2:3,sum)/1e3
        dimnames(CW)$year <- simYears
        ST = as.numeric(apply(CW[ac(seq(lStatPer$ST[1],lStatPer$ST[2])),],2,mean))
        MT = as.numeric(apply(CW[ac(seq(lStatPer$MT[1],lStatPer$MT[2])),],2,mean))
        LT = as.numeric(apply(CW[ac(seq(lStatPer$LT[1],lStatPer$LT[2])),],2,mean))
        
        StatName = "Yield"
        StatUnit = "(kt)"

      } else if (PerfStat %in% c("SSB","Risk3","Risk1")){

        #abundance
        Abd <- t[["N"]]
        #stock weights
        SW <- t[["stkW"]]
        #maturity
        Mat <- t[["Mat"]]

        #SSB (Mt)
        SSB <- apply(Abd*SW*Mat,2:3,sum)/1e6
        dimnames(SSB)$year <- simYears
        
        if (PerfStat == "SSB") {
          StatName = "SSB"
          StatUnit = "(Mt)"
          ST = as.numeric(apply(SSB[ac(seq(lStatPer$ST[1],lStatPer$ST[2])),],2,mean))
          MT = as.numeric(apply(SSB[ac(seq(lStatPer$MT[1],lStatPer$MT[2])),],2,mean))
          LT = as.numeric(apply(SSB[ac(seq(lStatPer$LT[1],lStatPer$LT[2])),],2,mean))
        } else if (PerfStat == "Risk3"){
          StatName = "Risk3"
          StatUnit = "p"
          #browser()
          #maximum probability that SSB is below Blim, where the maximum (of the annual probabilities) is taken over the statistical period 
          ST = max(apply(SSB[ac(seq(lStatPer$ST[1],lStatPer$ST[2])),]<Blim/1e6,1,sum)/dim(SSB)[2])
          MT = max(apply(SSB[ac(seq(lStatPer$MT[1],lStatPer$MT[2])),]<Blim/1e6,1,sum)/dim(SSB)[2])
          LT = max(apply(SSB[ac(seq(lStatPer$LT[1],lStatPer$LT[2])),]<Blim/1e6,1,sum)/dim(SSB)[2])
        } else if (PerfStat == "Risk1"){
          StatName = "Risk1"
          StatUnit = "p"
          #browser()
          #average probability that SSB is below Blim over the statistical period 
          ST = mean(apply(SSB[ac(seq(lStatPer$ST[1],lStatPer$ST[2])),]<Blim/1e6,1,sum)/dim(SSB)[2])
          MT = mean(apply(SSB[ac(seq(lStatPer$MT[1],lStatPer$MT[2])),]<Blim/1e6,1,sum)/dim(SSB)[2])
          LT = mean(apply(SSB[ac(seq(lStatPer$LT[1],lStatPer$LT[2])),]<Blim/1e6,1,sum)/dim(SSB)[2])
        }
        
      }
      
      dfAll <- dplyr::bind_rows(dfAll,
                                data.frame(RunRef = rep(r,length(ST)),
                                           Ftgt = rep(ftgt,length(ST)),
                                           PerfStat = rep(StatName,length(ST)),
                                           Period = rep(c("ST","MT","LT"),each=length(ST)),
                                           Val = c(ST,MT,LT),
                                           stringsAsFactors = FALSE))
      
    }
    
  }

  #geom_col for risks, boxplots otherwise
  p <- ggplot(data = dfAll, aes(x=factor(RunRef), y=Val, fill=factor(Period, levels = c("ST","MT","LT"))))
    
  if (PerfStat %in% c("Risk1","Risk3")){
    p <- p + geom_col(position="dodge") + geom_hline(yintercept = 0.05, col="black", linetype=2)
  } else {
    p <- p + geom_boxplot(outlier.size = 0.2)
  }
  
  p <- p +  
    facet_wrap(vars(paste("Ftgt = ",Ftgt))) +
    ggtitle(paste(StatName,StatUnit)) + 
    theme(legend.position="none",
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size=12, face="bold"),
          axis.text.x = element_text(size=8, face="bold"),
          strip.background = element_rect(colour="red", fill="#CCCCFF"))
    
  png(filename = file.path(Plot.dir,paste0(StatName,"Comparison.png")),
      type = 'cairo', units = 'in', width = 10,
      height = 8, pointsize = 12, res = 96)
  print(p)
  dev.off()

}