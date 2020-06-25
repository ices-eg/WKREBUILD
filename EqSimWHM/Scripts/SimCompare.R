#run comparisons

source(file.path(getwd(),"Scripts","setup.R"))

# #SSB vs Blim
# fAnnSSBvsBlimDist(OM = OM2, MP = MP2.0, res.dir = Res.dir, plot.dir = Res.dir)
# 

##############Ensure all scenrios to be compared have been run & stats calculated###############

#with/without assessment and advice error
runs2Compare <- c("OM2.2_MP1.0_1000_50","OM2.2_MP1.4_1000_50")

for (stat in c("Catch","SSB","Risk3","Risk1","IAV","IAVUpDown")){
   fCompare_runs(runs2Compare = runs2Compare, Res.dir = Res.dir, Plot.dir = Res.dir,
                 PerfStat = stat,
                 TargetFs = c(0, 0.05, 0.074, 0.1, 0.108, 0.2, 0.3),
                 lStatPer = list('ST'=c(2021,2025),'MT'=c(2026,2030),'LT'=c(2031,2067)),
                 Blim = OM2.2$refPts$Blim)}

#NoHCR/ ICES AR (Btrig=1.168Mt)
runs2Compare <- c("OM2.2_MP1.4_1000_50","OM2.2_MP2.1_1000_50")

for (stat in c("Catch","SSB","Risk3","Risk1","IAV","IAVUpDown")){
  fCompare_runs(runs2Compare = runs2Compare, Res.dir = Res.dir, Plot.dir = Res.dir,
                PerfStat = stat,
                TargetFs = c(0, 0.05, 0.074, 0.1, 0.108, 0.2, 0.3),
                lStatPer = list('ST'=c(2021,2025),'MT'=c(2026,2030),'LT'=c(2031,2067)),
                Blim = OM2.2$refPts$Blim)}

#comparison of stochastic/random weights & selection - min, max, min & max TAC
# runs2Compare <- c("OM2.2_MP1.0","OM2.2_MP1.1","OM2.2_MP1.2","OM2.2_MP1.3")
# for (stat in c("Catch","SSB","Risk3","Risk1")){
# fCompare_runs(runs2Compare = runs2Compare, Res.dir = Res.dir, Plot.dir = Res.dir,
#              PerfStat = stat, TargetFs = c(0,0.05,0.074,0.1,0.108,0.2),
#              lStatPer = lStatPer, Blim = OM$refPts$Blim)}


#IAV
# runs2Compare <- c("OM2.2_MP1.0","OM2.2_MP1.5","OM2.2_MP1.6","OM2.2_MP1.7","OM2.2_MP1.8","OM2.2_MP1.9")
# create folder
# dir.create(path = file.path(Res.dir,"Comparisons","IAV"), showWarnings = TRUE, recursive = TRUE)
# for (stat in c("Catch","SSB","Risk3","Risk1","IAV","IAVUpDown")){
#  fCompare_runs(runs2Compare = runs2Compare, Res.dir = Res.dir, Plot.dir = file.path(Res.dir,"Comparisons","IAV"),
#                PerfStat = stat, TargetFs = c(0,0.05,0.074,0.1,0.108,0.2),
#                lStatPer = lStatPer, Blim = OM$refPts$Blim)}