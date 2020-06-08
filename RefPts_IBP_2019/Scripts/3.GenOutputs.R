#compare results of several analyses

rm(list=ls())
gc()

source(file=file.path(getwd(),"Scripts","0.Setup.R"))

#12 plots, 1 per scenario

#read ouputs into a single data frame

Cases <- seq(1,6)
SubCases <- seq(1,3)

Case.labels <- c("Blim=SSB2001","Bpa=Bloss, Blim=Bpa/1.4","Blim=Bloss","Bpa=SSB2001, Blim=Bpa/1.4","Bpa=SSB2003, Blim=Bpa/1.4","Blim=SSB2003")
names(Case.labels) <- as.character(Cases)
SubCase.labels <- c("All years", "Ex 1982", "1995 on")
names(SubCase.labels) <- as.character(SubCases)

all.scenarios <- paste(rep(Cases,each=length(SubCases)),SubCases,sep=".")

#individual scenario plots
runNames <- c("Baseline") 
#runNames <- c("AdviceError")

#comparison (plotted together)
#runNames <- c("Baseline","AdviceError")
#runNames <- c("Baseline","RickerSRR","BevHoltSRR")
#runNames <- c("Baseline","Sens3yrSelection")
#runNames <- c("Baseline","SensRecrAuto")
#runNames <- c("Baseline","Sens3yrBiology")
#runNames <- c("Baseline","SensExremeHighRecr")
#runNames <- c("Baseline","SegRegSRR")

#valid output files are of type RData and contain the word cases
opFiles <- list.files(path = RData.dir, pattern = ".+cases.+\\.RData$", ignore.case = TRUE)
#available runs
runs <- sapply(str_split(opFiles,"_"),"[[",1)
#select those to read in
opFiles2Load <- file.path(RData.dir,opFiles[runs %in% runNames])

load(opFiles2Load[1])
dfAllResults <- dfResults
for (f in opFiles2Load[-1]){
  load(f)
  dfAllResults <- rbind(dfAllResults,dfResults)
}
rm(dfResults)

table(dfAllResults$Scenario)

scen <- names(table(dfAllResults$Scenario))

for (s in all.scenarios) {
  
  dfTemp <- dplyr::filter(dfAllResults,Scenario==s & Type=="Abs")
  
  if (nrow(dfTemp)>0) {

    Cairo(file = file.path(graphics.dir,paste0("Scenario",s,"_Absolute",".",fileFormat)), height=plot.h, width=plot.w, 
        units="in", dpi=96, pointsize=12)
    p <- ggplot(filter(dfTemp,RP %in% c("Flim","FMSY","FMSY_final","FP05","Fpa")), aes(x = factor(RP, levels = c("Flim","Fpa","FMSY","FP05","FMSY_final")), y = Val)) +
      geom_point(position = position_dodge(width=0.1,preserve=c("total")),aes(colour=runName, shape=Assessment), size=3) +
      labs(x = "F Reference Point", y = "Value") +
      ylim(0,0.2) +
      annotate("text", x=3, y= 0.2, label = paste0("Scenario ",s," (",Case.labels[substr(s,1,1)],",",SubCase.labels[substr(s,3,3)],")"))
    print(p)
    dev.off()
  
  }
  
  dfTemp <- dplyr::filter(dfAllResults,Scenario==s & Type=="Rel")

  if (nrow(dfTemp)>0) {
    
    Cairo(file = file.path(graphics.dir,paste0("Scenario",s,"_Relative.",fileFormat)), height=plot.h, width=plot.w, 
          units="in", dpi=96, pointsize=12)
    p <- ggplot(filter(dfTemp,RP %in% c("Flim","FMSY","FMSY_final","FP05","Fpa")), aes(x = factor(RP, levels = c("Flim","Fpa","FMSY","FP05","FMSY_final")), y = Val)) +
      geom_point(position = position_dodge(width=0.1,preserve=c("total")),aes(colour=runName, shape=Assessment), size=3) +
      labs(x = "F Reference Point", y = "Relative Value") +
      ylim(0,2) +
      annotate("text", x=3, y= 2, label = paste0("Scenario ",s," (",Case.labels[substr(s,1,1)],",",SubCase.labels[substr(s,3,3)],")"))
    print(p)
    dev.off()
  
  }
  
}

#results in report table format
for (s in scen){
  t1.b.a <- filter(dfAllResults,Scenario==s & Assessment=="WK17" & Type=="Abs" & RP %in% c('Blim','Bpa','MSYBtrigger')) %>% 
    select(Scenario,RP,Val) %>% 
    mutate(WK17.a=paste0(fmt(round(Val),dgt=0,fmt="d")," t")) %>%
    select(Scenario,RP,WK17.a)
  t1.b.r <- filter(dfAllResults,Scenario==s & Assessment=="WK17" & Type=="Rel" & RP %in% c('Blim','Bpa','MSYBtrigger')) %>% 
    select(Scenario,RP,Val) %>% 
    mutate(WK17.r=fmt(round(1000*Val)/1000),dgt=3,fmt="f") %>%
    select(WK17.r)

  t1.f.a <- filter(dfAllResults,Scenario==s & Assessment=="WK17" & Type=="Abs" & RP %in% c('Flim','Fpa','FMSY','FP05','FMSY_final')) %>% 
    select(Scenario,RP,Val) %>% 
    mutate(WK17.a=fmt(round(1000*Val)/1000),dgt=3,fmt="f") %>%
    select(Scenario,RP,WK17.a)
  t1.f.r <- filter(dfAllResults,Scenario==s & Assessment=="WK17" & Type=="Rel" & RP %in% c('Flim','Fpa','FMSY','FP05','FMSY_final')) %>% 
    select(Scenario,RP,Val) %>% 
    mutate(WK17.r=fmt(round(1000*Val)/1000),dgt=3,fmt="f") %>%
    select(WK17.r)
  
  t1.a <- dplyr::bind_rows(t1.b.a,t1.f.a)
  t1.r <- dplyr::bind_rows(t1.b.r,t1.f.r)
    
  t2.b.a <- filter(dfAllResults,Scenario==s & Assessment=="WG17" & Type=="Abs" & RP %in% c('Blim','Bpa','MSYBtrigger')) %>% 
    select(RP,Val) %>% 
    mutate(WG17.a=paste0(fmt(round(Val),dgt=0,fmt="d"), " t")) %>%
    select(RP,WG17.a)
  t2.b.r <- filter(dfAllResults,Scenario==s & Assessment=="WG17" & Type=="Rel" & RP %in% c('Blim','Bpa','MSYBtrigger')) %>% 
    select(RP,Val) %>% 
    mutate(WG17.r=fmt(round(1000*Val)/1000),dgt=3,fmt="f") %>%
    select(RP,WG17.r)
  
  t2.f.a <- filter(dfAllResults,Scenario==s & Assessment=="WG17" & Type=="Abs" & RP %in% c('Flim','Fpa','FMSY','FP05','FMSY_final')) %>% 
    select(RP,Val) %>% 
    mutate(WG17.a=fmt(round(1000*Val)/1000),dgt=3,fmt="f") %>%
    select(RP,WG17.a)
  t2.f.r <- filter(dfAllResults,Scenario==s & Assessment=="WG17" & Type=="Rel" & RP %in% c('Flim','Fpa','FMSY','FP05','FMSY_final')) %>% 
    select(RP,Val) %>% 
    mutate(WG17.r=fmt(round(1000*Val)/1000),dgt=3,fmt="f") %>%
    select(RP,WG17.r)
  
  t2.a <- dplyr::bind_rows(t2.b.a,t2.f.a)
  t2.r <- dplyr::bind_rows(t2.b.r,t2.f.r)
  
  t3.b.a <- filter(dfAllResults,Scenario==s & Assessment=="WG18" & Type=="Abs" & RP %in% c('Blim','Bpa','MSYBtrigger')) %>% 
    select(RP,Val) %>% 
    mutate(WG18.a=paste0(fmt(round(Val),dgt=0,fmt="d")," t")) %>%
    select(RP,WG18.a)

  t3.b.r <- filter(dfAllResults,Scenario==s & Assessment=="WG18" & Type=="Rel" & RP %in% c('Blim','Bpa','MSYBtrigger')) %>% 
    select(RP,Val) %>% 
    mutate(WG18.r=fmt(round(1000*Val)/1000),dgt=3,fmt="f") %>%
    select(RP,WG18.r)
  
  t3.f.a <- filter(dfAllResults,Scenario==s & Assessment=="WG18" & Type=="Abs" & RP %in% c('Flim','Fpa','FMSY','FP05','FMSY_final')) %>% 
    select(RP,Val) %>% 
    mutate(WG18.a=fmt(round(1000*Val)/1000),dgt=3,fmt="f") %>%
    select(RP,WG18.a)

  t3.f.r <- filter(dfAllResults,Scenario==s & Assessment=="WG18" & Type=="Rel" & RP %in% c('Flim','Fpa','FMSY','FP05','FMSY_final')) %>% 
    select(RP,Val) %>% 
    mutate(WG18.r=fmt(round(1000*Val)/1000),dgt=3,fmt="f") %>%
    select(RP,WG18.r)
  
  t3.a <- dplyr::bind_rows(t3.b.a,t3.f.a)
  t3.r <- dplyr::bind_rows(t3.b.r,t3.f.r)

  t4 <- dplyr::bind_cols(t1.a,
                         dplyr::select(t2.a,"WG17.a"),
                         dplyr::select(t3.a,"WG18.a"),
                         dplyr::select(t1.r,"WK17.r"),
                         dplyr::select(t2.r,"WG17.r"),
                         dplyr::select(t3.r,"WG18.r"))
  
  rownames(t4) <- t4$RP
  write.table(t4[c('Blim','Bpa','Flim','Fpa','FMSY','MSYBtrigger','FP05','FMSY_final'),
                 c('Scenario','RP','WK17.a','WG17.a','WG18.a','WK17.r','WG17.r','WG18.r')],
              file="Tabulated.csv",quote=TRUE,sep=",",row.names = FALSE, append = TRUE)
}
