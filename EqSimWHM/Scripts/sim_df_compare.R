# ================================================================================================================
# sim data.frame compare
# 
# Martin Pastoors
#
# 30/06/2020 compare different EqSim and SAM runs
# ================================================================================================================

# rm(list=ls())
# gc()

library(tidyverse)

#computer specific locations
Drive    <- "D:"
Base.dir <- file.path(Drive,"GIT")
MSE.dir <- file.path(Base.dir,"wk_WKREBUILD","EqSimWHM")
Res.dir <- file.path(MSE.dir, "Results")

# Get dropbox dir; for storing large RData files
source("EqSimWHM/R/get_dropbox.R")
source("EqSimWHM/R/utilities.R")
dropbox.dir <- file.path(get_dropbox(), "HOM FG", "05. Data","RData")

# ==================================================================================

#basic display setting
niters <- 1000
nyr <- 20

# simulation periods
per1 <- 5
per2 <- 5
# per3 is simply the remainder


stats.dir <- file.path(Res.dir, "Stats")

list.files.eqsim <- list.files(path=stats.dir, pattern="WHOM", full.names=TRUE)

dfall <- data.frame(stringsAsFactors = FALSE)
for (i in 1:length(list.files.eqsim)) {
#for (i in 1:3) {
  cat(list.files.eqsim[i],"\n")
  t <- loadRData(list.files.eqsim[i])[["stats"]]
  # f <- "0.05"
  for (f in names(t)) {
    cat(f," \n")
    x     <- t[[f]][["dfy"]] %>% lowcase() %>% filter(year <= 2037)
    dfall <- bind_rows(dfall, x)
  }
}

list.files.sam   <- list.files(path=stats.dir, pattern="^SAM", full.names=TRUE)
for (i in 1:length(list.files.sam)) {
  cat(list.files.sam[i],"\n")
  x <- 
    loadRData(list.files.sam[i]) %>% 
    lowcase() %>% 
    filter(year <= 2037) %>% 
    rename(runref = runname) %>% 
    mutate(iter = as.character(iter)) %>% 
    mutate(ftgt = as.character(ftgt)) %>% 
    mutate(value = ifelse(perfstat == "ssb", value/1000, value)) %>% 
    mutate(value = ifelse(perfstat == "catch", value/1000, value))
  dfall <- bind_rows(dfall, x)
}

# dfall <- dfall  %>%  filter(runref != "HCR2_SR")

# Summarize
dfsum <-
  dfall %>% 
  mutate(perfstat = tolower(perfstat)) %>% 
  mutate(perfstat = ifelse(perfstat == "cw", "catch", perfstat)) %>% 
  mutate(perfstat = ifelse(perfstat == "fbar", "harvest", perfstat)) %>% 
  filter(perfstat %in% c("catch","harvest","rec","ssb")) %>% 
  group_by(runref, ftgt, year, perfstat) %>% 
  summarize(mean = mean(value, na.rm=TRUE),
            upper = quantile(value, probs=0.975, na.rm=TRUE),
            lower = quantile(value, probs=0.025, na.rm=TRUE)) %>%
  ungroup()
  

# plot
dfsum %>%
  
  ggplot(aes(x=year, y=mean, group=runref)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  # theme(legend.position = "none") +

  geom_ribbon(aes(ymin=lower, ymax=upper, fill=runref), alpha=0.2) +
  geom_line(aes(colour=runref), size=0.8) +

  expand_limits(y=0) +
  labs(x="", y="value") +
  facet_grid(perfstat ~ ftgt, scales="free_y")

# ggsave(file = file.path(Res.dir,runName,paste0(runName,"_summary_byyear.png")),
#        device="png", width = 30, height = 20, units = "cm")

save(dfall,file = file.path(dropbox.dir,paste0("dfall.Rdata")))
save(dfsum,file = file.path(dropbox.dir,paste0("dfsum.Rdata")))
