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

stats.dir <- file.path(Res.dir, "Stats")

list.files.eqsim.ss <- list.files(path=stats.dir, pattern="WHOM_SS", full.names=TRUE)

dfall1 <- data.frame(stringsAsFactors = FALSE)
for (i in 1:length(list.files.eqsim.ss)) {
#for (i in 1:2) {
  cat(list.files.eqsim.ss[i],"\n")
  t <- loadRData(list.files.eqsim.ss[i])[["stats"]]
  # f <- "0.05"
  for (f in names(t)) {
    cat(f," \n")
    x     <- t[[f]][["dfy"]] %>% 
             lowcase() %>% 
             filter(year <= 2037) %>% 
             mutate(perfstat = tolower(perfstat)) %>% 
             mutate(perfstat = ifelse(perfstat=="cw","catch",perfstat)) %>% 
             mutate(method="EqSim")
    
    dfall1 <- bind_rows(dfall1, x)
  }
}

dfall2 <- data.frame(stringsAsFactors = FALSE)
list.files.eqsim.sam <- list.files(path=stats.dir, pattern="WHOM_SAM", full.names=TRUE)
for (i in 1:length(list.files.eqsim.sam)) {
  cat(list.files.eqsim.sam[i],"\n")
  x <- 
    loadRData(list.files.eqsim.sam[i]) %>% 
    lowcase() %>% 
    filter(year <= 2037) %>% 
    mutate(iter = as.character(iter)) %>% 
    mutate(ftgt = as.character(ftgt)) %>% 
    mutate(perfstat = tolower(perfstat)) %>% 
    mutate(perfstat = ifelse(perfstat == "cw", "catch", perfstat)) %>% 
    # mutate(value = ifelse(perfstat %in% c("catch","tac"), value/1000, value)) %>% 
    mutate(mp = ifelse(mp=="MP2.1", "MP5.1",mp)) %>% 
    mutate(mp = ifelse(mp=="MP2.2", "MP5.12",mp)) %>% 
    mutate(method="EqSim")
    
  dfall2 <- bind_rows(dfall2, x)
}

list.files.sam   <- list.files(path=stats.dir, pattern="^sam", full.names=TRUE)
dfall3 <- data.frame(stringsAsFactors = FALSE)
for (i in 1:length(list.files.sam)) {
  cat(list.files.sam[i],"\n")
  x <- 
    loadRData(list.files.sam[i]) %>% 
    lowcase() %>% 
    filter(year <= 2037) %>% 
    filter(!perfstat %in% c("land", "fbarl", "tsb")) %>%  
    mutate(iter = as.character(iter)) %>% 
    mutate(ftgt = as.character(ftgt)) %>% 
    mutate(perfstat = ifelse(perfstat == "fbar", "harvest", perfstat)) %>% 
    mutate(value = ifelse(perfstat == "ssb", value/1000, value)) %>% 
    mutate(value = ifelse(perfstat == "catch", value/1000, value)) %>% 
    mutate(assess="SAM") %>% 
    mutate(method="SAMhcr") 
  
  dfall3 <- bind_rows(dfall3, x)
}

# dfall <- dfall  %>%  mutate(assess = ifelse(is.na(assess), "SAM",assess))
# dfall <- dfall  %>%  mutate(perfstat = tolower(perfstat))
# dfall <- dfall  %>%  mutate(assess = ifelse(perfstat == "ssb" & assess=="SS3", value*1000,value))
# bind_rows(dfall1,dfall2,dfall3) %>%  group_by(method, assess, perfstat) %>%  summarise(value = as.integer(mean(value, na.rm=TRUE))) %>% View()
# dfall3 %>% ungroup() %>% distinct(ftgt)

dfall <- bind_rows(dfall1, dfall2, dfall3)

periods <-
  dfall %>% 
  ungroup() %>% 
  distinct(year, period) %>% 
  filter(!is.na(period))

periods2 <-
  periods %>% 
  group_by(period) %>% 
  filter(year == max(year))

  
# calculate IAV and prob below Blim
dfpblim <-
  dfall %>% 
  filter(perfstat == "ssb") %>% 
  mutate(blim = ifelse(assess=="SS3", 834, 612)) %>% 
  mutate(below = ifelse(value < blim, TRUE, FALSE)) %>% 
  group_by(method, assess, om, mp, runref, ftgt, year, perfstat) %>% 
  summarise(ntotal = dplyr::n(),
            # mean  = mean(value),
            nbelow = sum(below==TRUE)) %>% 
  mutate(mean = nbelow/ntotal) %>% 
  mutate(perfstat="pblim") %>% 
  dplyr::select(-ntotal, -nbelow)

# Summarize
dfsum <-
  dfall %>% 
  group_by(method, assess, om, mp, runref, ftgt, year, perfstat) %>% 
  summarize(mean = mean(value, na.rm=TRUE),
            upper = quantile(value, probs=0.975, na.rm=TRUE),
            lower = quantile(value, probs=0.025, na.rm=TRUE)) %>%
  ungroup() %>% 
  bind_rows(dfpblim) %>% 
  left_join(periods, by="year")


save(dfall,file = file.path(dropbox.dir,paste0("dfall.Rdata")))
save(dfsum,file = file.path(dropbox.dir,paste0("dfsum.Rdata")))

# load(file = file.path(dropbox.dir,paste0("dfall.Rdata")))
# load(file = file.path(dropbox.dir,paste0("dfsum.Rdata")))

worms <-
  dfall %>% 
  group_by(perfstat, runref, ftgt, year, period) %>% 
  arrange(runref, perfstat, ftgt, year, period) %>% 
  mutate(iter = row_number()) %>% 
  filter(row_number() <= 5) %>% 
  ungroup() %>% 
  mutate(perfstat = tolower(perfstat)) %>% 
  mutate(code = paste(method,assess,mp,sep="_")) 
  
dfsum %>% ungroup() %>% distinct(runref) 

# plot for single method and assessment; different MPs

myruns <- c("WHOM_SS3_OM2.2_MP5.0_1000_50", 
            "WHOM_SS3_OM2.2_MP5.01_1000_50",
            "WHOM_SS3_OM2.2_MP5.02_1000_50",
            "WHOM_SS3_OM2.2_MP5.03_1000_50")
myruns <- c("WHOM_SS3_OM2.2_MP5.1_1000_50", 
            "WHOM_SS3_OM2.2_MP5.11_1000_50",
            "WHOM_SS3_OM2.2_MP5.12_1000_50",
            "WHOM_SS3_OM2.2_MP5.13_1000_50")
myruns <- c("WHOM_SS3_OM2.2_MP5.2_1000_50", 
            "WHOM_SS3_OM2.2_MP5.21_1000_50",
            "WHOM_SS3_OM2.2_MP5.22_1000_50",
            "WHOM_SS3_OM2.2_MP5.23_1000_50")
myruns <- c("WHOM_SAM_OM2.3_MP2.1_1000_20", 
            "WHOM_SAM_OM2.3_MP2.2_1000_20")
myruns <- c("samhcr_WHOM_sam_-_5.1_1000_20")

myftgt <- c(0.05, 0.075, 0.10, 0.2)
myperfstat <- "ssb"; mycolour   <- "blue"; myyintercept <- 834
myperfstat <- "ssb"; mycolour   <- "blue"; myyintercept <- 612
myperfstat <- "harvest"; mycolour   <- "darkgreen"; myyintercept <- as.numeric(NA)
myperfstat <- "catch"; mycolour   <- "purple"; myyintercept <- as.numeric(NA)
myperfstat <- "rec"; mycolour   <- "orange"; myyintercept <- as.numeric(NA)
myperfstat <- "pblim"; mycolour   <- "red"; myyintercept <- 0.05

dfsum %>%
  filter(ftgt %in% myftgt) %>% 
  filter(perfstat==myperfstat) %>% 
  filter(runref %in% myruns) %>% 
  # separate(runref, into=c("stock","assess","om","mp","iters","nyears"), sep="_", convert=FALSE) %>% 
  mutate(code = paste(method, assess,om,mp,sep="_")) %>% 
  
  ggplot(aes(x=year, y=mean, group=mp)) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  theme( panel.grid.major.x = element_blank()) +
  theme(legend.position = "none") +
  
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill=mycolour) +

  geom_line(data=filter(worms, ftgt %in% myftgt,
                               perfstat %in% myperfstat,
                               runref %in% myruns), 
            aes(x=year, y=value, group=iter),
            size=0.5, colour="gray", inherit.aes = FALSE) +
  
  geom_line(colour=mycolour) +
  
  geom_hline(yintercept=myyintercept, linetype="dashed") +
  
  geom_vline(xintercept=periods2$year[1], linetype="dotted") +
  geom_vline(xintercept=periods2$year[2], linetype="dotted") +
  geom_vline(xintercept=periods2$year[3], linetype="dotted") +
  geom_vline(xintercept=periods2$year[4], linetype="dotted") +
  
  expand_limits(y=0) +
  labs(x="", y="value", title=myperfstat) +
  facet_grid(mp ~ ftgt, scales="free_y")

# Plot with comparison across methods and assessments

myruns <- c("WHOM_SS3_OM2.2_MP5.03_1000_50",
            "WHOM_SS3_OM2.2_MP5.13_1000_50",
            "WHOM_SS3_OM2.2_MP5.23_1000_50")

myruns <- c("WHOM_SS3_OM2.2_MP5.1_1000_50", 
            "WHOM_SAM_OM2.3_MP2.1_1000_20", 
            "samhcr_WHOM_sam_-_5.1_1000_20")

myftgt <- c(0.05, 0.075, 0.10, 0.2)
myperfstat <- "ssb"; mycolour   <- "blue"; myyintercept <- as.numeric(NA)
myperfstat <- "harvest"; mycolour   <- "darkgreen"; myyintercept <- as.numeric(NA)
myperfstat <- "catch"; mycolour   <- "purple"; myyintercept <- as.numeric(NA)
myperfstat <- "rec"; mycolour   <- "orange"; myyintercept <- as.numeric(NA)
myperfstat <- "pblim"; mycolour   <- "red"; myyintercept <- 0.05

dfsum %>%
  filter(ftgt %in% myftgt) %>% 
  filter(perfstat==myperfstat) %>% 
  filter(runref %in% myruns) %>% 
  # separate(runref, into=c("stock","assess","om","mp","iters","nyears"), sep="_", convert=FALSE) %>% 
  mutate(code = paste(method, assess,mp, sep="_")) %>% 
  mutate(blim = ifelse(assess=="SAM",612, 834)) %>% 
  
  ggplot(aes(x=year, y=mean, group=mp)) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  theme( panel.grid.major.x = element_blank()) +
  theme(legend.position = "none") +
  
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill=mycolour) +
  
  geom_line(data=filter(worms, ftgt %in% myftgt,
                        perfstat %in% myperfstat,
                        runref %in% myruns), 
            aes(x=year, y=value, group=iter),
            size=0.5, colour="gray", inherit.aes = FALSE) +
  
  geom_line(colour=mycolour) +
  
  geom_hline(yintercept=myyintercept, linetype="dashed") +
  #geom_hline(aes(yintercept=blim), linetype="dashed") +
  
  geom_vline(xintercept=periods2$year[1], linetype="dotted") +
  geom_vline(xintercept=periods2$year[2], linetype="dotted") +
  geom_vline(xintercept=periods2$year[3], linetype="dotted") +
  geom_vline(xintercept=periods2$year[4], linetype="dotted") +
  
  expand_limits(y=0) +
  labs(x="", y="value", title=myperfstat) +
  facet_grid(code ~ ftgt)
