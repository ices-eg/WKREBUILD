# ====================================================================================
# 1.SamtoFLStock.r
# ====================================================================================

library(FLCore)
library(FLSAM)
library(stockassessment)
library(tidyverse)

source("eqSimWHM/R/utilities.R")
source("../SAM2FLR/SAM2FLStock.r")

Drive    <- "D:"
Base.dir <- file.path(Drive,"GIT")
data.dir <- file.path(Base.dir,"wk_WKREBUILD","RefPts_IBP_2019_SAM","RData")

sam2018 <- SAM2FLSTOCK("WHOM_2018", temp="D:\\temp") 
t <- as.data.frame(fbar(sam2018)) %>% mutate(slot="fbar")
sam2018 <-
  sam2018 %>% 
  as.data.frame() %>% 
  bind_rows(t) %>% 
  mutate(
    assess = "sam",
    assessmentyear = "2018"
  ) 


sam2019 <- SAM2FLSTOCK("WHOM_2019", temp="D:\\temp") 
t <- as.data.frame(fbar(sam2019)) %>% mutate(slot="fbar")
sam2019 <-
  sam2019 %>% 
  as.data.frame() %>% 
  bind_rows(t) %>% 
  mutate(
    assess = "sam",
    assessmentyear = "2019"
  ) 


ss2018 <- loadRData(file=file.path(data.dir,"WGWIDE2018.RData"))  
stock(ss2018) <- FLCore::computeStock(ss2018) 
t <- as.data.frame(fbar(ss2018)) %>% mutate(slot="fbar")
ss2018 <-
  ss2018 %>% 
  as.data.frame() %>% 
  bind_rows(t) %>% 
  mutate(
    assess = "ss3",
    assessmentyear = "2018"
  ) 

ss2019 <- loadRData(file=file.path(data.dir,"WGWIDE2019.RData"))  
stock(ss2019) <- FLCore::computeStock(ss2019) 
t <- as.data.frame(fbar(ss2019)) %>% mutate(slot="fbar")
ss2019 <-
  ss2019 %>% 
  as.data.frame() %>% 
  bind_rows(t) %>% 
  mutate(
    assess = "ss3",
    assessmentyear = "2019"
  ) 

df <-
  bind_rows(ss2018,ss2019,sam2018,sam2019) %>% 
  mutate(assesscode = paste(assess,assessmentyear,sep="_"))

df %>% filter(slot=="fbar") %>% View()

df %>% 
  filter(tolower(slot) %in% c("stock","fbar")) %>% 
  ggplot(aes(x=year,y=data, group=assesscode)) +
  theme_publication() +
  geom_line(aes(colour=assess, linetype=assessmentyear, size=assessmentyear)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_size_manual(values=c(0.8, 1.0)) +
  facet_wrap(~slot, scales="free_y")

