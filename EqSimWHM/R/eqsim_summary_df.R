#eqsim_summary_df

# runs2Compare <- c("OM2.2_MP2.1_1000_20","OM2.2_MP2.2_1000_20")
# PerfStat = c("Catch","SSB","Risk3","IAVUpDown")
# TargetFs = c(0,0.05,0.074,0.1,0.2)
# lStatPer = lStatPer
# Blim <- OM$refPts$Blim
  
fsummary_df <- function(runs2Compare, 
                          Res.dir, Plot.dir, 
                          PerfStat, TargetFs,lStatPer, Blim){

  #produces a dataframe
  #comparing the supplied performance statistic for ST,MT and LT for each runName
  #grouped by runName
  
  require(tidyverse)
  
  dfAll <- data.frame(RunRef = c(), Label=c(), Ftgt = c(), PerfStat = c(), Period = c(), Val = c(), stringsAsFactors = FALSE)

  #get the data
  # r <- "OM2.2_MP2.1_1000_20"
  for (r in runs2Compare){
    
    cat(r,"\n")
    
    #load the output of the simulation and the summary statistics
    load(file = file.path(Res.dir,r,paste0(r,"_SimRuns.RData")))
    load(file = file.path(Res.dir,r,paste0(r,"_eqSim_Stats.RData")))

    # ftgt <- 0.0
    for (ftgt in TargetFs){
      
      cat("F ", ftgt,"\n")
      
      #simulation op
      t <- SimRuns[[ac(ftgt)]]
      #simulation stats
      # t2 <- lOp$stats[[ac(ftgt)]]
      t2 <- lStats$stats[[ac(ftgt)]]

      for (p in PerfStat){
        
        # cat("PerfStat: ",p, "\n")
        
        if (p %in% c("Catch","IAV","IAVUpDown")) {
        
          #catch numbers
          Cnum <- t[["C"]]
          #catch weights
          Cwgt <- t[["catW"]]
          #catch weight (tons)
          CW <- apply(Cnum*Cwgt,2:3,sum)/1e3
          dimnames(CW)$year <- t2$simYears
          
          if (p =="Catch"){
            
            # cat("Catch ","\n")
            
            ST = as.numeric(apply(CW[ac(seq(lStatPer$ST[1],lStatPer$ST[2])),],2,mean))
            MT = as.numeric(apply(CW[ac(seq(lStatPer$MT[1],lStatPer$MT[2])),],2,mean))
            LT = as.numeric(apply(CW[ac(seq(lStatPer$LT[1],lStatPer$LT[2])),],2,mean))
            
            StatName = "Yield"
            StatUnit = "(kt)"
          }
          
          if (p=="IAV"){
  
            # cat("IAV ","\n")
            
            IAV <- abs(1-CW[-1,]/CW[-nrow(CW),])
            #replace Inf with NA (NA results from comparing with zero catch)
            IAV <- ifelse(is.finite(IAV),IAV,NA)
            
            ST = apply(IAV[as.character(seq(lStatPer$ST[1],lStatPer$ST[2])),],2,mean)
            MT = apply(IAV[as.character(seq(lStatPer$MT[1],lStatPer$MT[2])),],2,mean)
            LT = apply(IAV[as.character(seq(lStatPer$LT[1],lStatPer$LT[2])),],2,mean)
            
            StatName = "IAV"
            StatUnit = ""
            
          }
          
          if (p=="IAVUpDown"){
  
            # cat("IAVUpDown ","\n")
            
            IAVUp <- IAVDown <- IAVUpDown <- (CW[-1,]/CW[-nrow(CW),]) - 1
            
            IAVUp[!is.finite(IAVUp)]<-NA
            IAVDown[!is.finite(IAVDown)]<-NA
            
            IAVUp[IAVUp<0] <- NA
            IAVDown[IAVDown>0] <- NA
            
            ST = c(apply(IAVUp[as.character(seq(lStatPer$ST[1],lStatPer$ST[2])),],2,mean,na.rm=TRUE),
                   apply(IAVDown[as.character(seq(lStatPer$ST[1],lStatPer$ST[2])),],2,mean,na.rm=TRUE))
            MT = c(apply(IAVUp[as.character(seq(lStatPer$MT[1],lStatPer$MT[2])),],2,mean,na.rm=TRUE),
                   apply(IAVDown[as.character(seq(lStatPer$MT[1],lStatPer$MT[2])),],2,mean,na.rm=TRUE))
            LT = c(apply(IAVUp[as.character(seq(lStatPer$LT[1],lStatPer$LT[2])),],2,mean,na.rm=TRUE),
                   apply(IAVDown[as.character(seq(lStatPer$LT[1],lStatPer$LT[2])),],2,mean,na.rm=TRUE))
            
            ST[!is.finite(ST)]<-NA
            MT[!is.finite(MT)]<-NA
            LT[!is.finite(LT)]<-NA
            
            StatName = c("IAVUp","IAVDown")
            StatUnit = c("")
            
          }
            
        } else if (p %in% c("SSB","Risk3","Risk1")){
  
          #abundance
          Abd <- t[["N"]]
          #stock weights
          SW <- t[["stkW"]]
          #maturity
          Mat <- t[["Mat"]]
  
          #SSB (Mt)
          SSB <- apply(Abd*SW*Mat,2:3,sum)/1e6
          dimnames(SSB)$year <- t2$simYears
          
          if (p == "SSB") {
            
            # cat("SSB ","\n")
            
            StatName = "SSB"
            StatUnit = "(Mt)"
            ST = as.numeric(apply(SSB[ac(seq(lStatPer$ST[1],lStatPer$ST[2])),],2,mean))
            MT = as.numeric(apply(SSB[ac(seq(lStatPer$MT[1],lStatPer$MT[2])),],2,mean))
            LT = as.numeric(apply(SSB[ac(seq(lStatPer$LT[1],lStatPer$LT[2])),],2,mean))
          
          } else if (p == "Risk3"){
            
            # cat("Risk3 ","\n")
            
            StatName = "Risk3"
            StatUnit = "p"
            #maximum probability that SSB is below Blim, where the maximum (of the annual probabilities) is taken over the statistical period 
            ST = max(apply(SSB[ac(seq(lStatPer$ST[1],lStatPer$ST[2])),]<Blim/1e6,1,sum)/dim(SSB)[2])
            MT = max(apply(SSB[ac(seq(lStatPer$MT[1],lStatPer$MT[2])),]<Blim/1e6,1,sum)/dim(SSB)[2])
            LT = max(apply(SSB[ac(seq(lStatPer$LT[1],lStatPer$LT[2])),]<Blim/1e6,1,sum)/dim(SSB)[2])
          
          } else if (p == "Risk1"){
            
            # cat("Risk1 ","\n")
            
            StatName = "Risk1"
            StatUnit = "p"
            #average probability that SSB is below Blim over the statistical period 
            ST = mean(apply(SSB[ac(seq(lStatPer$ST[1],lStatPer$ST[2])),]<Blim/1e6,1,sum)/dim(SSB)[2])
            MT = mean(apply(SSB[ac(seq(lStatPer$MT[1],lStatPer$MT[2])),]<Blim/1e6,1,sum)/dim(SSB)[2])
            LT = mean(apply(SSB[ac(seq(lStatPer$LT[1],lStatPer$LT[2])),]<Blim/1e6,1,sum)/dim(SSB)[2])
          }
          
        } # end of if statement for PerfStat p
        
        dfAll <- dplyr::bind_rows(dfAll,
                                  data.frame(RunRef = rep(r,length(ST)),
                                             Label = rep(t2$MP$xlab,length(ST)),
                                             Ftgt = rep(ftgt,length(ST)),
                                             PerfStat = rep(StatName,each=length(ST)/length(StatName)),
                                             Period = rep("ST",each=length(ST)),
                                             Period2 = paste(lStatPer$ST, collapse="-"),
                                             Val = ST,
                                             stringsAsFactors = FALSE))
        dfAll <- dplyr::bind_rows(dfAll,
                                  data.frame(RunRef = rep(r,length(MT)),
                                             Label = rep(t2$MP$xlab,length(MT)),
                                             Ftgt = rep(ftgt,length(MT)),
                                             PerfStat = rep(StatName,each=length(MT)/length(StatName)),
                                             Period = rep("MT",each=length(MT)),
                                             Period2 = paste(lStatPer$MT, collapse="-"),
                                             Val = MT,
                                             stringsAsFactors = FALSE))
        dfAll <- dplyr::bind_rows(dfAll,
                                  data.frame(RunRef = rep(r,length(LT)),
                                             Label = rep(t2$MP$xlab,length(LT)),
                                             Ftgt = rep(ftgt,length(LT)),
                                             PerfStat = rep(StatName,each=length(LT)/length(StatName)),
                                             Period = rep("LT",each=length(LT)),
                                             Period2 = paste(lStatPer$LT, collapse="-"),
                                             Val = LT,
                                             stringsAsFactors = FALSE))
        
      } # end of p loop

    } # end of for loop: ftgt
    
  } # end of for loop: r
  
  dfAll
  
} # end of function

# dfAll %>% 
#   filter(PerfStat == "Yield") %>% 
#   mutate(Ftgt = as.numeric(Ftgt)) %>% 
#   mutate(Period = factor(Period, levels=c("ST","MT","LT"))) %>% 
#   ggplot(aes(x=Ftgt, y=Val, group=Period)) +
#   theme_bw() +
#   geom_bar(stat="identity") +
#   facet_grid(RunRef ~ Period)

