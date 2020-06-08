#SRR investigations for Blim

source(file=file.path(getwd(),"Scripts","0.Setup.R"))

#most recent assessment
stk <- WG18
stk.name <- "WG18"
plot(stk)

dfSRDetResults <- data.frame()

#free fit of SRR models to full dataset, with and without 1982/2001 spikes
#WG18
models <- c("Ricker", "Segreg", "Bevholt")

WG18_3SRR <- eqsr_fit(WG18, nsamp=1000, models = models)
WG18_3SRR$sr.det$Assessment <- "WG18"
WG18_3SRR$sr.det$Terminal <- TRUE
WG18_3SRR$sr.det$Years <- "All"
WG18_3SRR$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG18_3SRR$sr.det)

WG18_3SRR_ex82 <- eqsr_fit(WG18, remove.years = c(1982), nsamp=1000, models = models)
WG18_3SRR_ex82$sr.det$Assessment <- "WG18"
WG18_3SRR_ex82$sr.det$Terminal <- TRUE
WG18_3SRR_ex82$sr.det$Years <- "ex82"
WG18_3SRR_ex82$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG18_3SRR_ex82$sr.det)

WG18_3SRR_ex01 <- eqsr_fit(WG18, remove.years = c(2001), nsamp=1000, models = models)
WG18_3SRR_ex01$sr.det$Assessment <- "WG18"
WG18_3SRR_ex01$sr.det$Terminal <- TRUE
WG18_3SRR_ex01$sr.det$Years <- "ex01"
WG18_3SRR_ex01$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG18_3SRR_ex01$sr.det)

WG18_3SRR_ex82_01 <- eqsr_fit(WG18, remove.years = c(1982,2001), nsamp=1000, models = models)
WG18_3SRR_ex82_01$sr.det$Assessment <- "WG18"
WG18_3SRR_ex82_01$sr.det$Terminal <- TRUE
WG18_3SRR_ex82_01$sr.det$Years <- "ex82_01"
WG18_3SRR_ex82_01$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG18_3SRR_ex82_01$sr.det)

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG18_3SRR",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18_3SRR, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG18_3SRR_ex82",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18_3SRR_ex82, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG18_3SRR_ex01",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18_3SRR_ex01, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG18_3SRR_ex82_01",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18_3SRR_ex82_01, y.mult=0.16)
dev.off()


#now without the terminal year estimates (considered very uncertain)
WG18_3SRR_exTerm <- eqsr_fit(WG18, remove.years = c(2017), nsamp=1000, models = models)
WG18_3SRR_exTerm$sr.det$Assessment <- "WG18"
WG18_3SRR_exTerm$sr.det$Terminal <- FALSE
WG18_3SRR_exTerm$sr.det$Years <- "All"
WG18_3SRR_exTerm$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG18_3SRR_exTerm$sr.det)

WG18_3SRR_exTerm_ex82 <- eqsr_fit(WG18, remove.years = c(1982,2017), nsamp=1000, models = models)
WG18_3SRR_exTerm_ex82$sr.det$Assessment <- "WG18"
WG18_3SRR_exTerm_ex82$sr.det$Terminal <- FALSE
WG18_3SRR_exTerm_ex82$sr.det$Years <- "ex82"
WG18_3SRR_exTerm_ex82$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG18_3SRR_exTerm_ex82$sr.det)

WG18_3SRR_exTerm_ex01 <- eqsr_fit(WG18, remove.years = c(2001,2017), nsamp=1000, models = models)
WG18_3SRR_exTerm_ex01$sr.det$Assessment <- "WG18"
WG18_3SRR_exTerm_ex01$sr.det$Terminal <- FALSE
WG18_3SRR_exTerm_ex01$sr.det$Years <- "ex01"
WG18_3SRR_exTerm_ex01$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG18_3SRR_exTerm_ex01$sr.det)

WG18_3SRR_exTerm_ex82_01 <- eqsr_fit(WG18, remove.years = c(1982,2001,2017), nsamp=1000, models = models)
WG18_3SRR_exTerm_ex82_01$sr.det$Assessment <- "WG18"
WG18_3SRR_exTerm_ex82_01$sr.det$Terminal <- FALSE
WG18_3SRR_exTerm_ex82_01$sr.det$Years <- "ex82_01"
WG18_3SRR_exTerm_ex82_01$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG18_3SRR_exTerm_ex82_01$sr.det)

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG18_3SRR_exTerm",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18_3SRR, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG18_3SRR_exTerm_ex82",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18_3SRR_ex82, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG18_3SRR_exTerm_ex01",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18_3SRR_ex01, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG18_3SRR_exTerm_ex82_01",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18_3SRR_ex82_01, y.mult=0.16)
dev.off()


#WG17
WG17_3SRR <- eqsr_fit(WG17, nsamp=1000, models = models)
WG17_3SRR$sr.det$Assessment <- "WG17"
WG17_3SRR$sr.det$Terminal <- TRUE
WG17_3SRR$sr.det$Years <- "All"
WG17_3SRR$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG17_3SRR$sr.det)

WG17_3SRR_ex82 <- eqsr_fit(WG17, remove.years = c(1982), nsamp=1000, models = models)
WG17_3SRR_ex82$sr.det$Assessment <- "WG17"
WG17_3SRR_ex82$sr.det$Terminal <- TRUE
WG17_3SRR_ex82$sr.det$Years <- "ex82"
WG17_3SRR_ex82$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG17_3SRR_ex82$sr.det)

WG17_3SRR_ex01 <- eqsr_fit(WG17, remove.years = c(2001), nsamp=1000, models = models)
WG17_3SRR_ex01$sr.det$Assessment <- "WG17"
WG17_3SRR_ex01$sr.det$Terminal <- TRUE
WG17_3SRR_ex01$sr.det$Years <- "ex01"
WG17_3SRR_ex01$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG17_3SRR_ex01$sr.det)

WG17_3SRR_ex82_01 <- eqsr_fit(WG17, remove.years = c(1982,2001), nsamp=1000, models = models)
WG17_3SRR_ex82_01$sr.det$Assessment <- "WG17"
WG17_3SRR_ex82_01$sr.det$Terminal <- TRUE
WG17_3SRR_ex82_01$sr.det$Years <- "ex82_01"
WG17_3SRR_ex82_01$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG17_3SRR_ex82_01$sr.det)

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG17_3SRR",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17_3SRR, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG17_3SRR_ex82",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17_3SRR_ex82, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG17_3SRR_ex01",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17_3SRR_ex01, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG17_3SRR_ex82_01",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17_3SRR_ex82_01, y.mult=0.16)
dev.off()


#now without the terminal year estimates (considered very uncertain)
WG17_3SRR_exTerm <- eqsr_fit(WG17, remove.years = c(2016), nsamp=1000, models = models)
WG17_3SRR_exTerm$sr.det$Assessment <- "WG17"
WG17_3SRR_exTerm$sr.det$Terminal <- FALSE
WG17_3SRR_exTerm$sr.det$Years <- "All"
WG17_3SRR_exTerm$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG17_3SRR_exTerm$sr.det)

WG17_3SRR_exTerm_ex82 <- eqsr_fit(WG17, remove.years = c(1982,2016), nsamp=1000, models = models)
WG17_3SRR_exTerm_ex82$sr.det$Assessment <- "WG17"
WG17_3SRR_exTerm_ex82$sr.det$Terminal <- FALSE
WG17_3SRR_exTerm_ex82$sr.det$Years <- "ex82"
WG17_3SRR_exTerm_ex82$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG17_3SRR_exTerm_ex82$sr.det)

WG17_3SRR_exTerm_ex01 <- eqsr_fit(WG17, remove.years = c(2001,2016), nsamp=1000, models = models)
WG17_3SRR_exTerm_ex01$sr.det$Assessment <- "WG17"
WG17_3SRR_exTerm_ex01$sr.det$Terminal <- FALSE
WG17_3SRR_exTerm_ex01$sr.det$Years <- "ex01"
WG17_3SRR_exTerm_ex01$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG17_3SRR_exTerm_ex01$sr.det)

WG17_3SRR_exTerm_ex82_01 <- eqsr_fit(WG17, remove.years = c(1982,2001,2016), nsamp=1000, models = models)
WG17_3SRR_exTerm_ex82_01$sr.det$Assessment <- "WG17"
WG17_3SRR_exTerm_ex82_01$sr.det$Terminal <- FALSE
WG17_3SRR_exTerm_ex82_01$sr.det$Years <- "ex82_01"
WG17_3SRR_exTerm_ex82_01$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG17_3SRR_exTerm_ex82_01$sr.det)

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG17_3SRR_exTerm",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17_3SRR, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG17_3SRR_exTerm_ex82",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17_3SRR_ex82, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG17_3SRR_exTerm_ex01",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17_3SRR_ex01, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG17_3SRR_exTerm_ex82_01",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17_3SRR_ex82_01, y.mult=0.16)
dev.off()


#WK17
WK17_3SRR <- eqsr_fit(WK17, nsamp=1000, models = models)
WK17_3SRR$sr.det$Assessment <- "WK17"
WK17_3SRR$sr.det$Terminal <- TRUE
WK17_3SRR$sr.det$Years <- "All"
WK17_3SRR$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WK17_3SRR$sr.det)

WK17_3SRR_ex82 <- eqsr_fit(WK17, remove.years = c(1982), nsamp=1000, models = models)
WK17_3SRR_ex82$sr.det$Assessment <- "WK17"
WK17_3SRR_ex82$sr.det$Terminal <- TRUE
WK17_3SRR_ex82$sr.det$Years <- "ex82"
WK17_3SRR_ex82$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WK17_3SRR_ex82$sr.det)

WK17_3SRR_ex01 <- eqsr_fit(WK17, remove.years = c(2001), nsamp=1000, models = models)
WK17_3SRR_ex01$sr.det$Assessment <- "WK17"
WK17_3SRR_ex01$sr.det$Terminal <- TRUE
WK17_3SRR_ex01$sr.det$Years <- "ex01"
WK17_3SRR_ex01$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WK17_3SRR_ex01$sr.det)

WK17_3SRR_ex82_01 <- eqsr_fit(WK17, remove.years = c(1982,2001), nsamp=1000, models = models)
WK17_3SRR_ex82_01$sr.det$Assessment <- "WK17"
WK17_3SRR_ex82_01$sr.det$Terminal <- TRUE
WK17_3SRR_ex82_01$sr.det$Years <- "ex82_01"
WK17_3SRR_ex82_01$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WK17_3SRR_ex82_01$sr.det)

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WK17_3SRR",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17_3SRR, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WK17_3SRR_ex82",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17_3SRR_ex82, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WK17_3SRR_ex01",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17_3SRR_ex01, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WK17_3SRR_ex82_01",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17_3SRR_ex82_01, y.mult=0.16)
dev.off()


#now without the terminal year estimates (considered very uncertain)
WK17_3SRR_exTerm <- eqsr_fit(WK17, remove.years = c(2015), nsamp=1000, models = models)
WK17_3SRR_exTerm$sr.det$Assessment <- "WK17"
WK17_3SRR_exTerm$sr.det$Terminal <- FALSE
WK17_3SRR_exTerm$sr.det$Years <- "All"
WK17_3SRR_exTerm$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WK17_3SRR_exTerm$sr.det)

WK17_3SRR_exTerm_ex82 <- eqsr_fit(WK17, remove.years = c(1982,2015), nsamp=1000, models = models)
WK17_3SRR_exTerm_ex82$sr.det$Assessment <- "WK17"
WK17_3SRR_exTerm_ex82$sr.det$Terminal <- FALSE
WK17_3SRR_exTerm_ex82$sr.det$Years <- "ex82"
WK17_3SRR_exTerm_ex82$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WK17_3SRR_exTerm_ex82$sr.det)

WK17_3SRR_exTerm_ex01 <- eqsr_fit(WK17, remove.years = c(2001,2015), nsamp=1000, models = models)
WK17_3SRR_exTerm_ex01$sr.det$Assessment <- "WK17"
WK17_3SRR_exTerm_ex01$sr.det$Terminal <- FALSE
WK17_3SRR_exTerm_ex01$sr.det$Years <- "ex01"
WK17_3SRR_exTerm_ex01$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WK17_3SRR_exTerm_ex01$sr.det)

WK17_3SRR_exTerm_ex82_01 <- eqsr_fit(WK17, remove.years = c(1982,2001,2015), nsamp=1000, models = models)
WK17_3SRR_exTerm_ex82_01$sr.det$Assessment <- "WK17"
WK17_3SRR_exTerm_ex82_01$sr.det$Terminal <- FALSE
WK17_3SRR_exTerm_ex82_01$sr.det$Years <- "ex82_01"
WK17_3SRR_exTerm_ex82_01$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WK17_3SRR_exTerm_ex82_01$sr.det)

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WK17_3SRR_exTerm",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17_3SRR, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WK17_3SRR_exTerm_ex82",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17_3SRR_ex82, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WK17_3SRR_exTerm_ex01",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17_3SRR_ex01, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WK17_3SRR_exTerm_ex82_01",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17_3SRR_ex82_01, y.mult=0.16)
dev.off()


#WG19 24/02/2020
WG19_3SRR <- eqsr_fit(WG19, nsamp=1000, models = models)
WG19_3SRR$sr.det$Assessment <- "WG19"
WG19_3SRR$sr.det$Terminal <- TRUE
WG19_3SRR$sr.det$Years <- "All"
WG19_3SRR$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG19_3SRR$sr.det)

WG19_3SRR_ex82 <- eqsr_fit(WG19, remove.years = c(1982), nsamp=1000, models = models)
WG19_3SRR_ex82$sr.det$Assessment <- "WG19"
WG19_3SRR_ex82$sr.det$Terminal <- TRUE
WG19_3SRR_ex82$sr.det$Years <- "ex82"
WG19_3SRR_ex82$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG19_3SRR_ex82$sr.det)

WG19_3SRR_ex01 <- eqsr_fit(WG19, remove.years = c(2001), nsamp=1000, models = models)
WG19_3SRR_ex01$sr.det$Assessment <- "WG19"
WG19_3SRR_ex01$sr.det$Terminal <- TRUE
WG19_3SRR_ex01$sr.det$Years <- "ex01"
WG19_3SRR_ex01$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG19_3SRR_ex01$sr.det)

WG19_3SRR_ex82_01 <- eqsr_fit(WG19, remove.years = c(1982,2001), nsamp=1000, models = models)
WG19_3SRR_ex82_01$sr.det$Assessment <- "WG19"
WG19_3SRR_ex82_01$sr.det$Terminal <- TRUE
WG19_3SRR_ex82_01$sr.det$Years <- "ex82_01"
WG19_3SRR_ex82_01$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG19_3SRR_ex82_01$sr.det)

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG19_3SRR",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG19_3SRR, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG19_3SRR_ex82",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG19_3SRR_ex82, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG19_3SRR_ex01",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG19_3SRR_ex01, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG19_3SRR_ex82_01",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG19_3SRR_ex82_01, y.mult=0.16)
dev.off()


#now without the terminal year estimates (considered very uncertain)
WG19_3SRR_exTerm <- eqsr_fit(WG19, remove.years = c(2018), nsamp=1000, models = models)
WG19_3SRR_exTerm$sr.det$Assessment <- "WG19"
WG19_3SRR_exTerm$sr.det$Terminal <- FALSE
WG19_3SRR_exTerm$sr.det$Years <- "All"
WG19_3SRR_exTerm$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG19_3SRR_exTerm$sr.det)

WG19_3SRR_exTerm_ex82 <- eqsr_fit(WG19, remove.years = c(1982,2018), nsamp=1000, models = models)
WG19_3SRR_exTerm_ex82$sr.det$Assessment <- "WG19"
WG19_3SRR_exTerm_ex82$sr.det$Terminal <- FALSE
WG19_3SRR_exTerm_ex82$sr.det$Years <- "ex82"
WG19_3SRR_exTerm_ex82$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG19_3SRR_exTerm_ex82$sr.det)

WG19_3SRR_exTerm_ex01 <- eqsr_fit(WG19, remove.years = c(2001,2018), nsamp=1000, models = models)
WG19_3SRR_exTerm_ex01$sr.det$Assessment <- "WG19"
WG19_3SRR_exTerm_ex01$sr.det$Terminal <- FALSE
WG19_3SRR_exTerm_ex01$sr.det$Years <- "ex01"
WG19_3SRR_exTerm_ex01$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG19_3SRR_exTerm_ex01$sr.det)

WG19_3SRR_exTerm_ex82_01 <- eqsr_fit(WG19, remove.years = c(1982,2001,2018), nsamp=1000, models = models)
WG19_3SRR_exTerm_ex82_01$sr.det$Assessment <- "WG19"
WG19_3SRR_exTerm_ex82_01$sr.det$Terminal <- FALSE
WG19_3SRR_exTerm_ex82_01$sr.det$Years <- "ex82_01"
WG19_3SRR_exTerm_ex82_01$sr.det$Models <- paste0(models,collapse="_")
dfSRDetResults <- dplyr::bind_rows(dfSRDetResults,WG19_3SRR_exTerm_ex82_01$sr.det)

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG19_3SRR_exTerm",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG19_3SRR, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG19_3SRR_exTerm_ex82",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG19_3SRR_ex82, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG19_3SRR_exTerm_ex01",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG19_3SRR_ex01, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG19_3SRR_exTerm_ex82_01",".",fileFormat)), height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG19_3SRR_ex82_01, y.mult=0.16)
dev.off()


#grouped bar charts, model proportions
p3ModProps_WG18_exTerm <- ggplot(data = filter(dfSRDetResults,Assessment=="WG18" & Terminal==FALSE),
       aes(Years, prop)) +
  geom_bar(aes(fill = model), position = "dodge", stat = "identity") + 
  theme(legend.position='none') +
  ylim(0,0.8) + annotate("text",x=1.5,y=0.8,label="WG18 (ex 2017)",size=5) +
  ylab("Proportion") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))

ggsave(filename = file.path(graphics.dir,"SRRInvestigations","3ModProps_WG18_exTerm.png"),
       plot = p3ModProps_WG18_exTerm,
       device = "png")


p3ModProps_WG19_exTerm <- ggplot(data = filter(dfSRDetResults,Assessment=="WG19" & Terminal==FALSE),
                                 aes(Years, prop)) +
  geom_bar(aes(fill = model), position = "dodge", stat = "identity") + 
  theme(legend.position='none') +
  ylim(0,0.8) + annotate("text",x=1.5,y=0.8,label="WG19 (ex 2018)",size=5) +
  ylab("Proportion") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))

ggsave(filename = file.path(graphics.dir,"SRRInvestigations","3ModProps_WG19_exTerm.png"),
       plot = p3ModProps_WG19_exTerm,
       device = "png")


p3ModProps_WG17_exTerm <- ggplot(data = filter(dfSRDetResults,Assessment=="WG17" & Terminal==FALSE),
       aes(Years, prop)) +
  geom_bar(aes(fill = model), position = "dodge", stat = "identity") + 
  theme(legend.position='none') +
  ylim(0,0.8) + annotate("text",x=1.5,y=0.8,label="WG17 (ex 2016)",size=5) +
  ylab("Proportion") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))

ggsave(filename = file.path(graphics.dir,"SRRInvestigations","3ModProps_WG17_exTerm.png"),
       plot = p3ModProps_WG17_exTerm,
       device = "png")

p3ModProps_WK17_exTerm <- ggplot(data = filter(dfSRDetResults,Assessment=="WK17" & Terminal==FALSE),
       aes(Years, prop)) +
  geom_bar(aes(fill = model), position = "dodge", stat = "identity") + 
  theme(legend.position='none') +
  ylim(0,0.8) + annotate("text",x=1.5,y=0.8,label="WK17 (ex 2015)",size=5) +
  ylab("Proportion") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))

ggsave(filename = file.path(graphics.dir,"SRRInvestigations","3ModProps_WK17_exTerm.png"),
       plot = p3ModProps_WK17_exTerm,
       device = "png")

p3ModProps_WG18 <- ggplot(data = filter(dfSRDetResults,Assessment=="WG18" & Terminal==TRUE),
       aes(Years, prop)) +
  geom_bar(aes(fill = model), position = "dodge", stat = "identity") + 
  theme(legend.position='none') +
  ylim(0,0.8) + annotate("text",x=1.5,y=0.8,label="WG18",size=5)

ggsave(filename = file.path(graphics.dir,"SRRInvestigations","3ModProps_WG18.png"),
       plot = p3ModProps_WG18,
       device = "png")

p3ModProps_WG19 <- ggplot(data = filter(dfSRDetResults,Assessment=="WG19" & Terminal==TRUE),
                          aes(Years, prop)) +
  geom_bar(aes(fill = model), position = "dodge", stat = "identity") + 
  theme(legend.position='none') +
  ylim(0,0.8) + annotate("text",x=1.5,y=0.8,label="WG19",size=5)

ggsave(filename = file.path(graphics.dir,"SRRInvestigations","3ModProps_WG19.png"),
       plot = p3ModProps_WG19,
       device = "png")

p3ModProps_WG17 <- ggplot(data = filter(dfSRDetResults,Assessment=="WG17" & Terminal==TRUE),
       aes(Years, prop)) +
  geom_bar(aes(fill = model), position = "dodge", stat = "identity") + 
  theme(legend.position='none') +
  ylim(0,0.8) + annotate("text",x=1,y=0.8,label="WG17",size=5)

ggsave(filename = file.path(graphics.dir,"SRRInvestigations","3ModProps_WG17.png"),
       plot = p3ModProps_WG17,
       device = "png")

p3ModProps_WK17 <- ggplot(data = filter(dfSRDetResults,Assessment=="WK17" & Terminal==TRUE),
       aes(Years, prop)) +
  geom_bar(aes(fill = model), position = "dodge", stat = "identity") + 
  theme(legend.position='none') +
  ylim(0,0.8) + annotate("text",x=1,y=0.8,label="WK17",size=5)

ggsave(filename = file.path(graphics.dir,"SRRInvestigations","3ModProps_WK17.png"),
       plot = p3ModProps_WK17,
       device = "png")



#if we were to use these results - what kind of individual fits would we be using?
#line plots or each formulation for WG18 ex 1982 and terminal year

#line plots
WG18_3SRR_exTerm_ex82$sr.det
ssb <- seq(1,5e6,by=1000)
st <- WG18_3SRR_exTerm_ex82$sr.sto

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG18_3SRR_exTerm_ex82_Ricker",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)

st.rck <- st[st$model=="Ricker",]

rec <- st.rck$a[1]*ssb*exp(-st.rck$b[1]*ssb)/1e6
plot(ssb/1e6,rec,type="l",xlab="SSB (Mt)",
     ylab="Recr (billions)",ylim=c(0,8),col="red",main="Ricker fits (WG18, ex1982, 2017)")
for (i in 2:nrow(st.rck)){
  rec <- st.rck$a[i]*ssb*exp(-st.rck$b[i]*ssb)/1e6
  lines(ssb/1e6,rec,col="red")
}
points(ssb(WG18)/1e6,rec(WG18)/1e6,pch=20,col="black")

dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG18_3SRR_exTerm_ex82_BevHolt",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)

#beverton-holt
st.bh <- st[st$model=="Bevholt",]

rec <- st.bh$a[1]*ssb/(1+st.bh$b[1]*ssb)/1e6
plot(ssb/1e6,rec,type="l",xlab="SSB (Mt)",ylab="Recr (billions)",ylim=c(0,8),col="red",main="BevHolt fits (WG18, ex1982, 2017)")
for (i in 2:nrow(st.bh)){
  rec <- st.bh$a[i]*ssb/(1+st.bh$b[i]*ssb)/1e6
  lines(ssb/1e6,rec,col="red")
}
points(ssb(WG18)/1e6,rec(WG18)/1e6,pch=20,col="black")

dev.off()

Cairo(file = file.path(graphics.dir,"SRRInvestigations",paste0("WG18_3SRR_exTerm_ex82_Segreg",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)

#segmented-regression
st.sr <- st[st$model=="Segreg",]
rec <- rep(NA,length(ssb))
rec[ssb<st.sr$b[1]]<-ssb[ssb<st.sr$b[1]]*st.sr$a[1]
rec[ssb>=st.sr$b[1]]<-st.sr$a[1]*st.sr$b[1]

plot(ssb/1e6,rec/1e6,type="l",xlab="SSB (Mt)",ylab="Recr (billions)",ylim=c(0,8),col="red",main="Segmented Regression fits (WG18, ex1982, 2017)")
for (i in 2:nrow(st.sr)){
  rec <- rep(NA,length(ssb))
  rec[ssb<st.sr$b[i]]<-ssb[ssb<st.sr$b[i]]*st.sr$a[i]
  rec[ssb>=st.sr$b[i]]<-st.sr$a[i]*st.sr$b[i]
  lines(ssb/1e6,rec/1e6,col="red")
}
points(ssb(WG18)/1e6,rec(WG18)/1e6,pch=20,col="black")

dev.off()



######################needs cleaning after here-------------///////////////################
















#Stock/Recruit dataframe (all years)
dfSRR <- dplyr::filter(dfWHM, Assessment == stk.name) %>% select(Year,SSB,Rec)
dfSRR.ranked <- dfSRR[order(dfSRR$Rec),] 

#year recruitment range scenarios
R1 <- seq(1982,2017) # all data
R2 <- seq(1983,2017) # ex 1982
R3 <- R2[-which(R2==2001)] # ex 1982, 2001 (spikes)

#some code here to show the contribution of an individual year class the SSB in any one year
#bubble plot
t <- stock.n(stk)*stock.wt(stk)*mat(stk)
t2<-t
for (a in 0:20) {t2[a+1,] <- ssb(stk)}
dfpSSB <- as(t/t2,"data.frame") %>% select(Age=age,Year=year,Prop=data) %>% mutate(YC=Year-Age)

#contribution to total SSB by year
g <- ggplot(data = dfpSSB, aes(x=Year,y=Age))
g <- g + geom_point(aes(size=Prop,fill=Prop),alpha=0.75,colour="black",stroke=1.5,pch=21,show.legend=FALSE)
g <- g + scale_size_area(max_size=15)
g <- g + theme_classic()
g <- g + theme(axis.text.x = element_text(angle=45, hjust=1))
g

#contribution to total SSB by year class
g <- ggplot(data = dfpSSB, aes(x=YC,y=Age))
g <- g + geom_point(aes(size=Prop,fill=Prop),alpha=0.75,colour="black",stroke=1.5,pch=21,show.legend=FALSE)
g <- g + scale_size_area(max_size=15)
g <- g + theme_classic()
g <- g + theme(axis.text.x = element_text(angle=45, hjust=1))
g

#line plots
g <- ggplot(data = filter(dfpSSB,Age>2 & Age<20 & YC>=1982), aes(x=Year,y=Prop))
g <- g + geom_line(aes(group=YC, col=YC), lwd=1)
g <- g + labs(x="Year", y="Proportion of total SSB")
g

#40% SSB threshold from single yc plus original years (1982,2001)
Rp40 <- R3[-which(R3 %in% unique(dfpSSB$Year[dfpSSB$Prop>0.4]))]
#50% SSB threshold from single yc plus original years (1982,2001)
Rp50 <- R3[-which(R3 %in% unique(dfpSSB$Year[dfpSSB$Prop>0.5]))]
#60% SSB threshold from single yc plus original years (1982,2001)
Rp60 <- R3[-which(R3 %in% unique(dfpSSB$Year[dfpSSB$Prop>0.6]))]

#1995 onward
R1995on <- R3[-which(R3 %in% unique(dfpSSB$Year[dfpSSB$Year<1995]))]

#some basic plots
#all years
pSSBvsRecr <- ggplot(data=filter(dfSRR,Year %in% R1),aes(x=SSB/1e6,y=Rec/1e6,label=Year)) +
  geom_point(size=3) + 
  xlim(0,NA) + ylim(0,NA) +
  labs(x="SSB (Mt)", y="Recruitment (blns)") +
  geom_text(hjust=-0.2,size=3) + 
  geom_path(aes(colour=as.numeric(Year)),size=1) +
  scale_colour_gradientn(colours = rev(heat.colors(10))) +
  theme(legend.position="none")

ggexport(pSSBvsRecr,filename=file.path(graphics.dir,"SRR_AllYears.png"),width=800,height=500)

#excluding 1982
pSSBvsRecr <- ggplot(data=filter(dfSRR,Year %in% R2),aes(x=SSB/1e6,y=Rec/1e6,label=Year)) +
  geom_point(size=3) + 
  xlim(0,NA) + ylim(0,NA) +
  labs(x="SSB (Mt)", y="Recruitment (blns)") +
  geom_text(hjust=-0.2,size=3) + 
  geom_path(aes(colour=as.numeric(Year)),size=1) +
  scale_colour_gradientn(colours = rev(heat.colors(10))) +
  theme(legend.position="none")

ggexport(pSSBvsRecr,filename=file.path(graphics.dir,"SRR_ex1982.png"),width=800,height=500)

#excluding 1982 and 2001
pSSBvsRecr <- ggplot(data=filter(dfSRR,Year %in% R3),aes(x=SSB/1e6,y=Rec/1e6,label=Year)) +
  geom_point(size=3) + 
  xlim(0,NA) + ylim(0,NA) +
  labs(x="SSB (Mt)", y="Recruitment (blns)") +
  geom_text(hjust=-0.2,size=3) + 
  geom_path(aes(colour=as.numeric(Year)),size=1) +
  scale_colour_gradientn(colours = rev(heat.colors(10))) +
  theme(legend.position="none")

ggexport(pSSBvsRecr,filename=file.path(graphics.dir,"SRR_exSpikes.png"),width=800,height=500)


#Recruitment percentiles
pct<-quantile(dfSRR$Rec, c(.25, .50, .75)) 

#minimum SSB with >75%ile recruitment
filter(dfSRR,Rec>=pct[3]) %>%
  select(SSB) %>%
  summarise(minSSB=min(SSB))
#872011 - 2017!!

#minimum SSB with >50%ile recruitment
filter(dfSRR,Rec>=pct[2]) %>%
  select(SSB) %>%
  summarise(minSSB=min(SSB))
#872011 - 2017!!

#minimum SSB with >25%ile recruitment
filter(dfSRR,Rec>=pct[1]) %>%
  select(SSB) %>%
  summarise(minSSB=min(SSB))
#872011 - 2017!!

gSR <- ggplot(data = dfSRR,
              aes(x = SSB/1e6, y = Rec/1e6)) +
  geom_point() + xlim(0,NA) + ylim(0,NA) + xlab("SSB (Mt)") +
  ylab("Rec (billions)") + geom_label(aes(label=Year), nudge_x = 0.2, size=2)

gSR









#line plots
ssb <- seq(1,5e6,by=1000)
st <- WG18SRR.all$sr.sto

Cairo(file = file.path(graphics.dir,paste0("WGWIDE18_SRR_All_Ricker",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)

st.rck <- st[st$model=="Ricker",]

rec <- st.rck$a[1]*ssb*exp(-st.rck$b[1]*ssb)/1e6
plot(ssb/1e6,rec,type="l",xlab="SSB (Mt)",
     ylab="Recr (billions)",ylim=c(0,8),col="red",main="Ricker - WG18 - all data points")
for (i in 2:nrow(st.rck)){
  rec <- st.rck$a[i]*ssb*exp(-st.rck$b[i]*ssb)/1e6
  lines(ssb/1e6,rec,col="red")
}
points(ssb(WG18)/1e6,rec(WG18)/1e6,pch=20,col="black")

dev.off()


Cairo(file = file.path(graphics.dir,paste0("WGWIDE18_SRR_All_Beverton-Holt",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)

#beverton-holt
st.bh <- st[st$model=="Bevholt",]

rec <- st.bh$a[1]*ssb/(1+st.bh$b[1]*ssb)/1e6
plot(ssb/1e6,rec,type="l",xlab="SSB (Mt)",ylab="Recr (billions)",ylim=c(0,8),col="red",main="Beverton-Holt - WG18 - all data points")
for (i in 2:nrow(st.bh)){
  rec <- st.bh$a[i]*ssb/(1+st.bh$b[i]*ssb)/1e6
  lines(ssb/1e6,rec,col="red")
}
points(ssb(WG18)/1e6,rec(WG18)/1e6,pch=20,col="black")

dev.off()

Cairo(file = file.path(graphics.dir,paste0("WGWIDE18_SRR_All_Segreg",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)

#segmented-regression
st.sr <- st[st$model=="Segreg",]
rec <- rep(NA,length(ssb))
rec[ssb<st.sr$b[1]]<-ssb[ssb<st.sr$b[1]]*st.sr$a[1]
rec[ssb>=st.sr$b[1]]<-st.sr$a[1]*st.sr$b[1]

plot(ssb/1e6,rec/1e6,type="l",xlab="SSB (Mt)",ylab="Recr (billions)",ylim=c(0,8),col="red",main="Segmented Regression - WG18 - all data points")
for (i in 2:nrow(st.sr)){
  rec <- rep(NA,length(ssb))
  rec[ssb<st.sr$b[i]]<-ssb[ssb<st.sr$b[i]]*st.sr$a[i]
  rec[ssb>=st.sr$b[i]]<-st.sr$a[i]*st.sr$b[i]
  lines(ssb/1e6,rec/1e6,col="red")
}
points(ssb(WG18)/1e6,rec(WG18)/1e6,pch=20,col="black")

dev.off()

#excluding 1982

Cairo(file = file.path(graphics.dir,paste0("WGWIDE18_SRR_All_no82",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18SRR.all.no82, y.mult=0.16)
dev.off()

#line plots
ssb <- seq(1,5e6,by=1000)
st <- WG18SRR.all.no82$sr.sto

Cairo(file = file.path(graphics.dir,paste0("WGWIDE18_SRR_All_no82_Ricker",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)

st.rck <- st[st$model=="Ricker",]

rec <- st.rck$a[1]*ssb*exp(-st.rck$b[1]*ssb)/1e6
plot(ssb/1e6,rec,type="l",xlab="SSB (Mt)",
     ylab="Recr (billions)",ylim=c(0,8),col="red",main="Ricker - WG18 - no 1982")
for (i in 2:nrow(st.rck)){
  rec <- st.rck$a[i]*ssb*exp(-st.rck$b[i]*ssb)/1e6
  lines(ssb/1e6,rec,col="red")
}
points(ssb(WG18)/1e6,rec(WG18)/1e6,pch=20,col="black")

dev.off()


Cairo(file = file.path(graphics.dir,paste0("WGWIDE18_SRR_All_no82_Beverton-Holt",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)

#beverton-holt
st.bh <- st[st$model=="Bevholt",]

rec <- st.bh$a[1]*ssb/(1+st.bh$b[1]*ssb)/1e6
plot(ssb/1e6,rec,type="l",xlab="SSB (Mt)",
     ylab="Recr (billions)",ylim=c(0,8),col="red",main="Beverton-Holt - WG18 - no 1982")
for (i in 2:nrow(st.bh)){
  rec <- st.bh$a[i]*ssb/(1+st.bh$b[i]*ssb)/1e6
  lines(ssb/1e6,rec,col="red")
}
points(ssb(WG18)/1e6,rec(WG18)/1e6,pch=20,col="black")

dev.off()

Cairo(file = file.path(graphics.dir,paste0("WGWIDE18_SRR_All_no82_Segreg",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)

#segmented-regression
st.sr <- st[st$model=="Segreg",]
rec <- rep(NA,length(ssb))
rec[ssb<st.sr$b[1]]<-ssb[ssb<st.sr$b[1]]*st.sr$a[1]
rec[ssb>=st.sr$b[1]]<-st.sr$a[1]*st.sr$b[1]

plot(ssb/1e6,rec/1e6,type="l",xlab="SSB (Mt)",
     ylab="Recr (billions)",ylim=c(0,8),col="red",main="Segmented Regression - WG18 - no 1982")
for (i in 2:nrow(st.sr)){
  rec <- rep(NA,length(ssb))
  rec[ssb<st.sr$b[i]]<-ssb[ssb<st.sr$b[i]]*st.sr$a[i]
  rec[ssb>=st.sr$b[i]]<-st.sr$a[i]*st.sr$b[i]
  lines(ssb/1e6,rec/1e6,col="red")
}
points(ssb(WG18)/1e6,rec(WG18)/1e6,pch=20,col="black")

dev.off()

#excluding 2001

Cairo(file = file.path(graphics.dir,paste0("WGWIDE18_SRR_All_no01",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18SRR.all.no01, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WGWIDE18_SRR_All_no82_01",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18SRR.all.no82_01, y.mult=0.16)
dev.off()

#WG18 Ricker, ex 1982
WG18SRR.Rck.no82 <- eqsr_fit(WG18, remove.years = c(1982), nsamp=1000, models = c("Ricker"))
Cairo(file = file.path(graphics.dir,paste0("WGWIDE18_SRR_Rck_no82",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18SRR.Rck.no82, y.mult=0.16)
dev.off()

#line plots
ssb <- seq(1,5e6,by=1000)
st <- WG18SRR.Rck.no82$sr.sto

Cairo(file = file.path(graphics.dir,paste0("WGWIDE18_SRR_Ricker_no82",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)

st.rck <- st[st$model=="Ricker",]

rec <- st.rck$a[1]*ssb*exp(-st.rck$b[1]*ssb)/1e6
plot(ssb/1e6,rec,type="l",xlab="SSB (Mt)",
     ylab="Recr (billions)",ylim=c(0,8),col="red",main="Ricker - WG18 - no 1982")
for (i in 2:nrow(st.rck)){
  rec <- st.rck$a[i]*ssb*exp(-st.rck$b[i]*ssb)/1e6
  lines(ssb/1e6,rec,col="red")
}
points(ssb(WG18)/1e6,rec(WG18)/1e6,pch=20,col="black")

dev.off()




#WG17
WG17SRR.all <- eqsr_fit(WG17, nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"))
WG17SRR.all.no82 <- eqsr_fit(WG17, remove.years = c(1982), nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"))
WG17SRR.all.no01 <- eqsr_fit(WG17, remove.years = c(2001), nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"))
WG17SRR.all.no82_01 <- eqsr_fit(WG17, remove.years = c(1982,2001), nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"))

Cairo(file = file.path(graphics.dir,paste0("WGWIDE17_SRR_All",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17SRR.all, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WGWIDE17_SRR_All_no82",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17SRR.all.no82, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WGWIDE17_SRR_All_no01",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17SRR.all.no01, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WGWIDE17_SRR_All_no82_01",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17SRR.all.no82_01, y.mult=0.16)
dev.off()

#WK17
WK17SRR.all <- eqsr_fit(WK17, nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"))
WK17SRR.all.no82 <- eqsr_fit(WK17, remove.years = c(1982), nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"))
WK17SRR.all.no01 <- eqsr_fit(WK17, remove.years = c(2001), nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"))
WK17SRR.all.no82_01 <- eqsr_fit(WK17, remove.years = c(1982,2001), nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"))

Cairo(file = file.path(graphics.dir,paste0("WKWIDE17_SRR_All",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17SRR.all, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WKWIDE17_SRR_All_no82",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17SRR.all.no82, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WKWIDE17_SRR_All_no01",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17SRR.all.no01, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WKWIDE17_SRR_All_no82_01",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17SRR.all.no82_01, y.mult=0.16)
dev.off()


#identify the changepoint by fitting a segmented regression

#full dataset
WG18SRR.segreg1 <- eqsr_fit(WG18, nsamp=1000, models = c("Segreg"))
WG17SRR.segreg1 <- eqsr_fit(WG17, nsamp=1000, models = c("Segreg"))
WK17SRR.segreg1 <- eqsr_fit(WK17, nsamp=1000, models = c("Segreg"))

Cairo(file = file.path(graphics.dir,paste0("WKWIDE17_SRR_SegReg",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17SRR.segreg1, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WGWIDE17_SRR_SegReg",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17SRR.segreg1, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WGWIDE18_SRR_SegReg",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18SRR.segreg1, y.mult=0.16)
dev.off()

WG18SRR.segreg1[[2]]
WG17SRR.segreg1[[2]]
WK17SRR.segreg1[[2]]


#excluding 1982
WG18SRR.segreg1.no82 <- eqsr_fit(WG18, remove.years = c(1982), nsamp=1000, models = c("Segreg"))
WG17SRR.segreg1.no82 <- eqsr_fit(WG17, remove.years = c(1982), nsamp=1000, models = c("Segreg"))
WK17SRR.segreg1.no82 <- eqsr_fit(WK17, remove.years = c(1982), nsamp=1000, models = c("Segreg"))

Cairo(file = file.path(graphics.dir,paste0("WKWIDE17_SRR_SegReg_ex82",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17SRR.segreg1.no82, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WGWIDE17_SRR_SegReg_ex82",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17SRR.segreg1.no82, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WGWIDE18_SRR_SegReg_ex82",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18SRR.segreg1.no82, y.mult=0.16)
dev.off()

WG18SRR.segreg2[[2]]
WG17SRR.segreg2[[2]]
WK17SRR.segreg2[[2]]


#excluding 1982,2001
WG18SRR.segreg3 <- eqsr_fit(WG18, remove.years = c(1982,2001), nsamp=1000, models = c("Segreg"))
WG17SRR.segreg3 <- eqsr_fit(WG17, remove.years = c(1982,2001), nsamp=1000, models = c("Segreg"))
WK17SRR.segreg3 <- eqsr_fit(WK17, remove.years = c(1982,2001), nsamp=1000, models = c("Segreg"))

Cairo(file = file.path(graphics.dir,paste0("WKWIDE17_SRR_SegReg_ex8201",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17SRR.segreg3, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WGWIDE17_SRR_SegReg_ex8201",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17SRR.segreg3, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WGWIDE18_SRR_SegReg_ex8201",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18SRR.segreg3, y.mult=0.16)
dev.off()

WG18SRR.segreg3[[2]]
WG17SRR.segreg3[[2]]
WK17SRR.segreg3[[2]]


#excluding 2001
WG18SRR.segreg4 <- eqsr_fit(WG18, remove.years = c(2001), nsamp=1000, models = c("Segreg"))
WG17SRR.segreg4 <- eqsr_fit(WG17, remove.years = c(2001), nsamp=1000, models = c("Segreg"))
WK17SRR.segreg4 <- eqsr_fit(WK17, remove.years = c(2001), nsamp=1000, models = c("Segreg"))

Cairo(file = file.path(graphics.dir,paste0("WKWIDE17_SRR_SegReg_ex01",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17SRR.segreg4, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WGWIDE17_SRR_SegReg_ex01",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17SRR.segreg4, y.mult=0.16)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WGWIDE18_SRR_SegReg_ex01",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18SRR.segreg4, y.mult=0.16)
dev.off()

#fit for 1995 onwards
WG18SRR.segreg5 <- eqsr_fit(window(WG18,1995,2017), remove.years = c(), nsamp=1000, models = c("Segreg"))
WG17SRR.segreg5 <- eqsr_fit(window(WG17,1995,2016), remove.years = c(), nsamp=1000, models = c("Segreg"))
WK17SRR.segreg5 <- eqsr_fit(window(WK17,1995,2015), remove.years = c(), nsamp=1000, models = c("Segreg"))

Cairo(file = file.path(graphics.dir,paste0("WKWIDE17_SRR_SegReg_95on",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WK17SRR.segreg5, y.mult=0.5)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WGWIDE17_SRR_SegReg_95on",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG17SRR.segreg5, y.mult=0.5)
dev.off()

Cairo(file = file.path(graphics.dir,paste0("WGWIDE18_SRR_SegReg_95on",".",fileFormat)), 
      height=plot.h, width=plot.w, units="in", dpi=96, pointsize=12)
eqsr_plot(WG18SRR.segreg5, y.mult=0.5)
dev.off()





#################################
















sigmalnSSB <- 0.194
Fcv <- 0.212
Fphi <- 0.423

#candidate Blim values
#Bloss - changeable
Blim1 <- 0.8e6
Bpa1 <- Blim1*exp(1.645*sigmalnSSB)
#Breakpoint full ts - poor fit
Blim2 <- 2e6
Bpa2 <- Blim2*exp(1.645*sigmalnSSB)
#SSB2001
Blim3 <- 1.35e6
Bpa3 <- Blim3*exp(1.645*sigmalnSSB)

#different SRR approach for the long term simulation - fit all models
SRR.all <- eqsr_fit(stk, nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"))
eqsr_plot(SRR.all, y.mult=0.5)
#no 1982, 2001
SRR.all.exspikes <- eqsr_fit(stk, nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"), remove.years = c(1982,2001))
eqsr_plot(SRR.all.exspikes, y.mult=0.5)
#no 1982
SRR.all.ex1982 <- eqsr_fit(stk, nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"), remove.years = c(1982,2001))
eqsr_plot(SRR.all.ex1982, y.mult=0.5)


#1995 on, no 2001
SRR.all.95on <- eqsr_fit(window(stk,1995,2017), nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"))
eqsr_plot(SRR.all.95on, y.mult=0.5)
#BevHolt won't fit
SRR.all.95on <- eqsr_fit(window(stk,1995,2017), nsamp=1000, models = c("Ricker", "Segreg"))
eqsr_plot(SRR.all.95on, y.mult=0.5)

#SR with breakpoint set at Blim
SegregBlim  <- function(ab, ssb) log(ifelse(ssb >= Blim3, ab$a * Blim3, ab$a * ssb))
SRR.SR.exspikes <- eqsr_fit(stk, nsamp=1000, models = c("SegregBlim"))
eqsr_plot(SRR.SR.exspikes, y.mult=0.5)
SRR.SR.95on <- eqsr_fit(window(stk,1995,2017), nsamp=1000, models = c("SegregBlim"))
eqsr_plot(SRR.SR.95on, y.mult=0.5)

#last 10
#bio.years = c(2008,2017)
#sel.years = c(2008,2017)

#last 3
bio.years = c(2015,2017)
sel.years = c(2015,2017)

SIM.Flim <- eqsim_run(SRR.SR.95on, bio.years = bio.years,
                         bio.const = FALSE,
                         sel.years = sel.years,
                         sel.const = FALSE,
                         Fscan = seq(0.01,0.2,len=40),
                         Fcv = 0, Fphi = 0,
                         rhologRec = F,
                         Blim = Blim1, Bpa = Bpa1,
                         Btrigger=0)

eqsim_plot(SIM.Flim,catch=TRUE)
Flim <- SIM.Flim$Refs2["catF","F50"]
Fpa <- Flim/1.4


SIM.Fmsy1 <- eqsim_run(SRR.SR.95on,
                          bio.years = bio.years,
                          bio.const = FALSE,
                          sel.years = sel.years,
                          sel.const = FALSE,
                          Fscan = seq(0.01,0.2,len=40),
                          Fcv = Fcv, Fphi = Fphi,
                          rhologRec = F,
                          Blim = Blim1, Bpa = Bpa1,
                          Btrigger=0)

Fmsy1 <- SIM.Fmsy1$Refs2["lanF","medianMSY"]
eqsim_plot_range(SIM.Fmsy1, type="median")

MSYBtrigger1 <- Bpa1

#Step 3 evaluate the ICES MSY rule and compare Fmsy with Fp05
SIM.Fmsy2 <- eqsim_run(SRR.SR.95on,
                       bio.years = bio.years,
                       bio.const = FALSE,
                       sel.years = sel.years,
                       sel.const = FALSE,
                       Fscan = seq(0.01,0.2,len=40),
                          Fcv = Fcv, Fphi = Fphi,
                          rhologRec = F,
                          Blim = Blim1, Bpa = Bpa1,
                          Btrigger = MSYBtrigger1)

eqsim_plot(SIM.Fmsy2,catch=TRUE)
Fp05 <- SIM.Fmsy2$Refs2["catF","F05"]

#summary of results
cat("**********************************************\n")
cat("Blim = ",Blim1,"\n")
cat("Bpa = ",Bpa1,"\n")
cat("Flim = ",Flim,"\n")
cat("Fpa = ",Fpa,"\n")
cat("MSYBtrigger = ",MSYBtrigger1,"\n")
cat("FMSY = ",Fmsy1,"\n")
cat("FP05 = ",Fp05,"\n")
cat("**********************************************\n")


#Blim = 1.35Mt

SIM.Flim <- eqsim_run(SRR.SR.exspikes, bio.years = bio.years,
                      bio.const = FALSE,
                      sel.years = sel.years,
                      sel.const = FALSE,
                      Fscan = seq(0.01,0.2,len=40),
                      Fcv = 0, Fphi = 0,
                      rhologRec = F,
                      Blim = Blim3, Bpa = Bpa3,
                      Btrigger=0)

eqsim_plot(SIM.Flim,catch=TRUE)
Flim2 <- SIM.Flim$Refs2["catF","F50"] 
Fpa2 <- Flim2/1.4
0.087/0.062

#productivity is lower using recent data only
SIM.Flim <- eqsim_run(SRR.SR.95on, bio.years = bio.years,
                      bio.const = FALSE,
                      sel.years = sel.years,
                      sel.const = FALSE,
                      Fscan = seq(0.01,0.2,len=40),
                      Fcv = 0, Fphi = 0,
                      rhologRec = F,
                      Blim = Blim3, Bpa = Bpa3,
                      Btrigger=0)

eqsim_plot(SIM.Flim,catch=TRUE)
Flim2 <- SIM.Flim$Refs2["catF","F50"] 
Fpa2 <- Flim2/1.4
0.064/0.045















#seg reg fits vs 1 missing year
dfResults = data.frame(Year=seq(1982,2017), 
                       AParam = rep(NA,length(seq(1982,2017))), 
                       BParam = rep(NA,length(seq(1982,2017))))

for (y in seq(1982,2017)){
  SRR <- eqsr_fit(stk,
                  remove.years = c(y), nsamp=1000, 
                  models = c("Segreg"))
  dfResults$AParam[dfResults$Year==y] <- SRR[[2]]$a 
  dfResults$BParam[dfResults$Year==y] <- SRR[[2]]$b 
}

plot(x=c(0,dfResults$BParam[dfResults$Year==1982],5e6),
     y=c(0,dfResults$BParam[dfResults$Year==1982]*dfResults$AParam[dfResults$Year==1982],
         dfResults$BParam[dfResults$Year==1982]*dfResults$AParam[dfResults$Year==1982]),
     type="l",ylim=c(0,5e6),xlab="SSB (Mt)",ylab="Recr")
for (y in seq(1983,2017)){
  lines(x=c(0,dfResults$BParam[dfResults$Year==y],5e6),
        c(0,dfResults$BParam[dfResults$Year==y]*dfResults$AParam[dfResults$Year==y],
          dfResults$BParam[dfResults$Year==y]*dfResults$AParam[dfResults$Year==y]))
}



#seg reg fits vs 1 missing year
dfResults.no82 = data.frame(Year=seq(1983,2017), 
                       AParam = rep(NA,length(seq(1983,2017))), 
                       BParam = rep(NA,length(seq(1983,2017))))

for (y in seq(1983,2017)){
  SRR <- eqsr_fit(stk,
                  remove.years = c(1982,y), nsamp=1000, 
                  models = c("Segreg"))
  dfResults.no82$AParam[dfResults.no82$Year==y] <- SRR[[2]]$a 
  dfResults.no82$BParam[dfResults.no82$Year==y] <- SRR[[2]]$b 
}

plot(x=c(0,dfResults.no82$BParam[dfResults.no82$Year==1983],3e6),
     y=c(0,dfResults.no82$BParam[dfResults.no82$Year==1983]*dfResults.no82$AParam[dfResults$Year==1982],
         dfResults.no82$BParam[dfResults.no82$Year==1983]*dfResults.no82$AParam[dfResults$Year==1982]),
     type="l",ylim=c(0,3e6),xlab="SSB (Mt)",ylab="Recr")
for (y in seq(1984,2017)){
  lines(x=c(0,dfResults.no82$BParam[dfResults.no82$Year==y],5e6),
        c(0,dfResults.no82$BParam[dfResults.no82$Year==y]*dfResults.no82$AParam[dfResults.no82$Year==y],
          dfResults.no82$BParam[dfResults.no82$Year==y]*dfResults.no82$AParam[dfResults.no82$Year==y]))
}
abline(v=min(ssb(stk)),col="red")

dfResults.no82
hist(dfResults.no82$BParam)


#seg reg fits vs 1 missing year
yrs <- c(seq(1982,2010),seq(2012,2017))
dfResults.no11 = data.frame(Year=yrs, 
                            AParam = rep(NA,length(yrs)), 
                            BParam = rep(NA,length(yrs)))

for (y in yrs){
  SRR <- eqsr_fit(stk,
                  remove.years = c(2011,y), nsamp=1000, 
                  models = c("Segreg"))
  dfResults.no11$AParam[dfResults.no11$Year==y] <- SRR[[2]]$a 
  dfResults.no11$BParam[dfResults.no11$Year==y] <- SRR[[2]]$b 
}

plot(x=c(0,3e6),y=c(0,3e6),type="n",ylim=c(0,5e6),xlab="SSB (Mt)",ylab="Recr")
for (y in yrs){
  lines(x=c(0,dfResults.no11$BParam[dfResults.no11$Year==y],5e6),
        c(0,dfResults.no11$BParam[dfResults.no11$Year==y]*dfResults.no11$AParam[dfResults.no11$Year==y],
          dfResults.no11$BParam[dfResults.no11$Year==y]*dfResults.no11$AParam[dfResults.no11$Year==y]))
}
abline(v=min(ssb(stk)),col="red")



#seg reg fits vs 1 missing year
yrs <- c(seq(1983,2010),seq(2012,2017))
dfResults.no82_11 = data.frame(Year=yrs, 
                            AParam = rep(NA,length(yrs)), 
                            BParam = rep(NA,length(yrs)),
                            CV = rep(NA,length(yrs)),
                            llik = rep(NA,length(yrs)))

for (y in yrs){
  SRR <- eqsr_fit(stk,
                  remove.years = c(1982,2011,y), nsamp=1000, 
                  models = c("Segreg"))
  dfResults.no82_11$AParam[dfResults.no82_11$Year==y] <- SRR[[2]]$a 
  dfResults.no82_11$BParam[dfResults.no82_11$Year==y] <- SRR[[2]]$b 
  dfResults.no82_11$CV[dfResults.no82_11$Year==y] <- SRR[[2]]$cv
  dfResults.no82_11$llik[dfResults.no82_11$Year==y] <- SRR[[2]]$llik
}

plot(x=c(0,3e6),y=c(0,3e6),type="n",ylim=c(0,5e6),xlab="SSB (Mt)",ylab="Recr")
for (y in yrs){
  lines(x=c(0,dfResults.no82_11$BParam[dfResults.no82_11$Year==y],5e6),
        c(0,dfResults.no82_11$BParam[dfResults.no82_11$Year==y]*dfResults.no82_11$AParam[dfResults.no82_11$Year==y],
          dfResults.no82_11$BParam[dfResults.no82_11$Year==y]*dfResults.no82_11$AParam[dfResults.no82_11$Year==y]))
}
abline(v=min(ssb(stk)),col="red")


#excluding 2 spikes
yrs <- c(seq(1983,2000),seq(2002,2017))
dfResults.no82_01 = data.frame(Year=yrs, 
                               AParam = rep(NA,length(yrs)), 
                               BParam = rep(NA,length(yrs)),
                               CV = rep(NA,length(yrs)),
                               llik = rep(NA,length(yrs)))

for (y in yrs){
  SRR <- eqsr_fit(stk,
                  remove.years = c(1982,2001,y), nsamp=1000, 
                  models = c("Segreg"))
  dfResults.no82_01$AParam[dfResults.no82_01$Year==y] <- SRR[[2]]$a 
  dfResults.no82_01$BParam[dfResults.no82_01$Year==y] <- SRR[[2]]$b 
  dfResults.no82_01$CV[dfResults.no82_01$Year==y] <- SRR[[2]]$cv
  dfResults.no82_01$llik[dfResults.no82_01$Year==y] <- SRR[[2]]$llik
}

plot(x=c(0,3e6),y=c(0,3e6),type="n",ylim=c(0,5e6),xlab="SSB (Mt)",ylab="Recr")
for (y in yrs){
  lines(x=c(0,dfResults.no82_01$BParam[dfResults.no82_01$Year==y],5e6),
        c(0,dfResults.no82_01$BParam[dfResults.no82_01$Year==y]*dfResults.no82_01$AParam[dfResults.no82_01$Year==y],
          dfResults.no82_01$BParam[dfResults.no82_01$Year==y]*dfResults.no82_01$AParam[dfResults.no82_01$Year==y]))
}
abline(v=min(ssb(stk)),col="red")





SRR.R1 <- eqsr_fit(stk,
                remove.years = c(),
                nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"))

eqsr_plot(SRR.R1, y.mult=0.16)

SRR.R1.SegReg <- eqsr_fit(stk,
                   remove.years = c(),
                   nsamp=1000, models = c("Segreg"))

eqsr_plot(SRR.R1.SegReg, y.mult=0.16)

SRR.R2 <- eqsr_fit(stk,
                   remove.years = setdiff(R1,R2),
                   nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"))

eqsr_plot(SRR.R2, y.mult=0.16)

SRR.R2.SegReg <- eqsr_fit(stk,
                          remove.years = setdiff(R1,R2),
                          nsamp=1000, models = c("Segreg"))

eqsr_plot(SRR.R2.SegReg, y.mult=0.16)

SRR.R2.SegReg.no2017 <- eqsr_fit(stk,
                          remove.years = c(setdiff(R1,R2),2017),
                          nsamp=1000, models = c("Segreg"))

eqsr_plot(SRR.R2.SegReg.no2017, y.mult=0.16)

SRR.R3 <- eqsr_fit(stk,
                   remove.years = setdiff(R1,R3),
                   nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"))

eqsr_plot(SRR.R3, y.mult=0.16)

SRR.R3.SegReg <- eqsr_fit(stk,
                   remove.years = setdiff(R1,R3),
                   nsamp=1000, models = c("Segreg"))

eqsr_plot(SRR.R3.SegReg, y.mult=0.16)

SRR.R3.SegReg.no2017 <- eqsr_fit(stk,
                          remove.years = c(setdiff(R1,R3),2017),
                          nsamp=1000, models = c("Segreg"))

eqsr_plot(SRR.R3.SegReg.no2017, y.mult=0.16)

SRR.R3.SegReg.no2017[[2]]

SRR.R4 <- eqsr_fit(stk,
                   remove.years = setdiff(R1,R4),
                   nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"))

eqsr_plot(SRR.R4, y.mult=0.16)

SRR.R4.SegReg <- eqsr_fit(stk,
                   remove.years = setdiff(R1,R4),
                   nsamp=1000, models = c("Segreg"))

eqsr_plot(SRR.R4.SegReg, y.mult=0.16)

SRR.R5 <- eqsr_fit(stk,
                   remove.years = setdiff(R1,R5),
                   nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"))

eqsr_plot(SRR.R5, y.mult=0.16)

SRR.R5.SegReg <- eqsr_fit(stk,
                   remove.years = setdiff(R1,R5),
                   nsamp=1000, models = c("Segreg"))

eqsr_plot(SRR.R5.SegReg, y.mult=0.16)

SRR.R6 <- eqsr_fit(stk,
                   remove.years = setdiff(R1,R6),
                   nsamp=1000, models = c("Ricker", "Segreg", "Bevholt"))

eqsr_plot(SRR.R6, y.mult=0.16)

SRR.R6.SegReg <- eqsr_fit(stk,
                   remove.years = setdiff(R1,R6),
                   nsamp=1000, models = c("Segreg"))

eqsr_plot(SRR.R6.SegReg, y.mult=0.16)

SRR.R6.SegReg.no2017 <- eqsr_fit(stk,
                          remove.years = c(setdiff(R1,R6),2017),
                          nsamp=1000, models = c("Segreg"))

eqsr_plot(SRR.R6.SegReg, y.mult=0.16)

#autocorrelation
acf(dfSRR$Rec)


#options for Blim
#1 - the lowest SSB that produced a high recruitment
Blim <- min(dfSRR$SSB[dfSRR$Year %in% c(1982,2001)])
#2 - Bloss
#Blim <- min(dfSRR$SSB)

#Assessment uncertainty at the start of the year following the terminal year
#of the assessment
#same values as WKWIDE 2017 - need to get updated values
sdlogSSB <- 0.195
sdlogFBar <- 0.195

#calculate Bpa
Bpa <- Blim*exp(1.645*sdlogSSB)

#EqSim
#function for segmented regression with breakpoint at Blim.1
SegregBlim  <- function(ab, ssb) log(ifelse(ssb >= Blim, ab$a * Blim, ab$a * ssb))

#using all SR data
SRR <- eqsr_fit(stk,
                remove.years = c(),
                nsamp=1000, models = c("SegregBlim"))

eqsr_plot(SRR)

#Flim with 98on seg reg
SIM.Flim <- eqsim_run(SRR,
                      bio.years = c(2008,2017),
                      bio.const = FALSE,
                      sel.years = c(2008,2017),
                      sel.const = FALSE,
                      Fscan = seq(0.01,0.3,len=30),
                      Fcv = 0, Fphi = 0,
                      rhologRec = T,
                      Blim = Blim, Bpa = Bpa,
                      Btrigger=0)

eqsim_plot(SIM.Flim,catch=TRUE)
Flim <- SIM.Flim$Refs2["catF","F50"]
Fpa <- Flim*exp(-1.645*sdlogFBar)
cat("Flim,Fpa=",Flim,Fpa,"\n")

#step 1 - no Btrigger but including stochasticity in 
#population and fishery and assessment/advice error
#Fmsy - 1st run, include assessment error
#for initial FMSY candidate

#defaults
Fcv <- 0.212
Fphi <- 0.423 

SIM1.Fmsy <- eqsim_run(SRR,
                       bio.years = c(2008,2017),
                       bio.const = FALSE,
                       sel.years = c(2008,2017),
                       sel.const = FALSE,
                       Fscan = seq(0.01,0.3,len=30),
                       Fcv = Fcv, Fphi = Fphi,
                       rhologRec = T,
                       Blim = Blim, Bpa = Bpa,
                       Btrigger=0)

Fmsy <- SIM1.Fmsy$Refs2["lanF","medianMSY"]
eqsim_plot_range(SIM1.Fmsy, type="median")

#Step 2 select MSYBtrigger
#Stock has not been fished at FMSY, use Bpa as the MSYBtrigger value
MSYBtrigger <- Bpa

#Step 3 evaluate the ICES MSY rule
SIM2.Fmsy <- eqsim_run(SRR,
                       bio.years = c(2008,2017),
                       bio.const = FALSE,
                       sel.years = c(2008,2017),
                       sel.const = FALSE,
                       Fscan = seq(0.01,0.3,len=30),
                       Fcv = Fcv, Fphi = Fphi,
                       rhologRec = T,
                       Blim = Blim, Bpa = Bpa,
                       Btrigger = MSYBtrigger)

eqsim_plot_range(SIM2.Fmsy, type="median")
eqsim_plot(SIM2.Fmsy,catch=TRUE)

Fp05 <- SIM2.Fmsy$Refs2["catF","F05"]
cat(Fmsy,Fp05,"\n")
