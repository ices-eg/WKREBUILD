#setup

rm(list=ls())
gc()
try(dev.off(),silent=FALSE)

#devtools::install_github("ices-tools-prod/msy")
library(msy)
library(ggplot2)
#library(ggpubr)
library(knitr)
#library(icesAdvice)
library(tidyverse)
#library(r4ss)
library(Cairo)
library(FLCore)
library(knitr)
library(plyr)

ac <- function(str){as.character(str)}
fmt <- function(RP,dgt=3,fmt="f"){ac(formatC(RP,digits=dgt,format=fmt,big.mark = ","))}
roundUp <- function(x,to=1000) {to*(x%/%to + as.logical(x%%to))}

subfolder <- "RefPts_IBP_2019_SAM"
RData.dir <- file.path(getwd(),subfolder, "RData")
graphics.dir <- file.path(getwd(),subfolder, "Graphics")
logs.dir <- file.path(getwd(),subfolder, "Logs")
source.dir <- file.path(getwd(),subfolder, "Source")

#WGWIDE2018
WG18 <- get(load(file=file.path(RData.dir,"WGWIDE2018.RData")))
name(WG18) <- "WGWIDE18"
WG18df <-
  as.data.frame(WG18) %>% 
  bind_rows(as.data.frame(fbar(WG18)) %>% mutate(slot="fbar")) %>% 
  bind_rows(as.data.frame(ssb(WG18)) %>% mutate(slot="ssb")) %>% 
  mutate(run="WG18")

#WGWIDE2019
WG19 <- get(load(file=file.path(RData.dir,"WGWIDE2019.RData")))
name(WG19) <- "WGWIDE19"
WG19df <-
  as.data.frame(WG19) %>% 
  bind_rows(as.data.frame(fbar(WG19)) %>% mutate(slot="fbar")) %>% 
  bind_rows(as.data.frame(ssb(WG19)) %>% mutate(slot="ssb")) %>% 
  mutate(run="WG19")

#WGWIDE2019SAM

WG19SAM <- get(load(file=file.path(RData.dir,"WGWIDE2019SAM.RData")))
name(WG19SAM) <- "WGWIDE19SAM"

WG19SAMdf <-
  as.data.frame(WG19SAM) %>% 
  bind_rows(as.data.frame(fbar(WG19SAM)) %>% mutate(slot="fbar")) %>% 
  bind_rows(as.data.frame(ssb(WG19SAM)) %>% mutate(slot="ssb")) %>% 
  mutate(run="WG19SAM")



#plot
bind_rows(WG18df, WG19df, WG19SAMdf) %>% 
  filter(slot %in% c("ssb","fbar")) %>% 
  
  ggplot(aes(x=year,y=data)) +
  theme_bw() +
  geom_line(aes(colour=run)) +
  facet_wrap(~slot, scales="free_y")
  