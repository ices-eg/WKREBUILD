#sandbox

source(file=file.path(getwd(),"Scripts","0.Setup.R"))

#case 4 Blim = SSB2001/1.4
Blim <- as.numeric(ssb(WG18)[,'2001',,,]/1.4)
#case 5 Blim = SSB2003/1.4
#Blim <- as.numeric(ssb(WG18)[,'2003',,,]/1.4)

SegregBlim  <- function(ab, ssb) log(ifelse(ssb >= Blim, ab$a * Blim, ab$a * ssb))

SRR.uncons <- eqsr_fit(window(WG18,1983,2016),remove.years = c(),nsamp=1000, models = c("Segreg"))
eqsr_plot(SRR.uncons)

SRR.cons <- eqsr_fit(window(WG18,1983,2016),remove.years = c(),nsamp=1000, models = c("SegregBlim"))
eqsr_plot(SRR.cons)

Ricker <- eqsr_fit(window(WG18,1983,2016),remove.years = c(),nsamp=1000, models = c("Ricker"))
eqsr_plot(Ricker)
#cumulative distribution

#generate data from models (robbed from eqsr_plot function)

x.mult <- 1.1
y.mult <- 1.4
n <- 20000
n_mods <- floor(n/100)
sample_rec <- function(i) {
  FUN <- match.fun(modset$model[i])
  exp(FUN(modset[i, ], ssb_eval) + stats::rnorm(length(ssb_eval), 
                                                sd = modset$cv[i]))
}

# tmp <- paste(SRR.uncons$sr.det$model, SRR.uncons$sr.det$prop)
# names(tmp) <- SRR.uncons$sr.det$model
# out$Model <- tmp[out$model]
# fit50 <- mgcv::gam(rec ~ s(ssb, m = 0), data = data.frame(ssb = ssb_eval, 
#                                                           rec = apply(rec_sim, 1, stats::quantile, 0.5)))
# fit05 <- mgcv::gam(rec ~ s(ssb, m = 0), data = data.frame(ssb = ssb_eval, 
#                                                           rec = apply(rec_sim, 1, stats::quantile, 0.05)))
# fit95 <- mgcv::gam(rec ~ s(ssb, m = 0), data = data.frame(ssb = ssb_eval, 
#                                                           rec = apply(rec_sim, 1, stats::quantile, 0.95)))
# Percentiles <- data.frame(ssb = ssb_eval, p50 = stats::fitted(fit50), 
#                           p05 = stats::fitted(fit05), p95 = stats::fitted(fit95))
# Percentiles <- Percentiles[stats::complete.cases(Percentiles),]
# 
# plot(0, 0, type = "n", xlim = c(0, maxSSB), ylim = c(0,maxrec), las = 1, 
#      xlab = "Spawning stock biomass", ylab = "Recruitment", 
#      main = paste("Predictive distribution of recruitment\nfor", 
#                   SRR.uncons$id.sr))
# 
# points(out$ssb, out$rec, pch = 20, cex = 1, col = grDevices::grey(0,alpha = 0.02))
# lines(p50 ~ ssb, col = 7, lwd = 3, data = Percentiles)
# lines(p05 ~ ssb, col = 4, lwd = 3, data = Percentiles)
# lines(p95 ~ ssb, col = 4, lwd = 3, data = Percentiles)
#   
# x <- SRR.uncons$sr.det
# y <- seq(1, round(maxSSB), length = 100)
# for (i in 1:nrow(x)) {
#   lines(y, exp(match.fun(as.character(x$model[i]))(x[i, 
#                                                      ], y)), col = "black", lwd = 2, lty = i)
# }
# lines(data$ssb, data$rec, col = "red")
# points(data$ssb, data$rec, pch = 19, col = "red", cex = 1.25)
# for (i in 1:nrow(SRR.uncons$sr.det)) {
#   text(0.2 * maxSSB, maxrec * (1 - i/10), paste(SRR.uncons$sr.det$model[i], 
#                                                 round(SRR.uncons$sr.det$prop[i], 2)), cex = 0.9, adj = 1)
# }


#ecdf comparison betweem assessment and eqsr_fit
#for constrained model
modset <- SRR.cons$sr.sto
data <- SRR.cons$rby[,c("year","rec","ssb")]
minSSB <- min(data$ssb, max(data$ssb) * 0.0125)
maxSSB <- max(data$ssb) * x.mult
maxrec <- max(data$rec) * y.mult

ssb_eval <- seq(minSSB, maxSSB, length.out = 100)
ids <- sample(1:nrow(modset), n_mods, replace = TRUE)
rec_sim <- sapply(ids, sample_rec)

out <- data.frame(grp = rep(1:length(ssb_eval), n_mods), 
                  mid.grp = rep(ssb_eval, n_mods), 
                  ssb = jitter(rep(ssb_eval,n_mods), 2), 
                  rec = c(rec_sim), 
                  model = rep(modset[ids,"model"], each = length(ssb_eval)))

mod <- ecdf(out$rec)
plot(mod,xlim=c(0,2e7),xlab="Recruits (thousands)",ylab="Cumulative Proportion",
     main = paste0("ECDF(Recruits-ex1982), Constrained SRR @ Blim ", round(Blim/1e6,2), "Mt"))
obs <- ecdf(as.numeric(rec(window(WG18,1983,2016))))
plot(obs, add=TRUE, col="red")            
legend(x=0, y=0.95,legend=c("eqsr_fit","WG18"),bty="n",lty=1,col=c("black","red"),pch=c(NA_integer_,20))



#Unconstrained Seg Reg
modset <- SRR.uncons$sr.sto
data <- SRR.uncons$rby[,c("year","rec","ssb")]
minSSB <- min(data$ssb, max(data$ssb) * 0.0125)
maxSSB <- max(data$ssb) * x.mult
maxrec <- max(data$rec) * y.mult

ssb_eval <- seq(minSSB, maxSSB, length.out = 100)
ids <- sample(1:nrow(modset), n_mods, replace = TRUE)
rec_sim <- sapply(ids, sample_rec)

out <- data.frame(grp = rep(1:length(ssb_eval), n_mods), 
                  mid.grp = rep(ssb_eval, n_mods), 
                  ssb = jitter(rep(ssb_eval,n_mods), 2), 
                  rec = c(rec_sim), 
                  model = rep(modset[ids,"model"], each = length(ssb_eval)))

mod <- ecdf(out$rec)
plot(mod,xlim=c(0,2e7),xlab="Recruits",ylab="Cumulative Proportion",main="ECDF(Recruits), Unconstrained SRR")
obs <- ecdf(as.numeric(rec(window(WG18,1983,2016))))
plot(obs, add=TRUE, col="red")            
legend("bottomright",legend=c("eqsr_fit","WG18"),bty="n",lty=1,col=c("black","red"),pch=c(NA_integer_,20))

#Ricker
modset <- Ricker$sr.sto
data <- Ricker$rby[,c("year","rec","ssb")]
minSSB <- min(data$ssb, max(data$ssb) * 0.0125)
maxSSB <- max(data$ssb) * x.mult
maxrec <- max(data$rec) * y.mult

ssb_eval <- seq(minSSB, maxSSB, length.out = 100)
ids <- sample(1:nrow(modset), n_mods, replace = TRUE)
rec_sim <- sapply(ids, sample_rec)

out <- data.frame(grp = rep(1:length(ssb_eval), n_mods), 
                  mid.grp = rep(ssb_eval, n_mods), 
                  ssb = jitter(rep(ssb_eval,n_mods), 2), 
                  rec = c(rec_sim), 
                  model = rep(modset[ids,"model"], each = length(ssb_eval)))

mod <- ecdf(out$rec)
plot(mod,xlim=c(0,2e7),xlab="Recruits",ylab="Cumulative Proportion",main="ECDF(Recruits), Ricker")
obs <- ecdf(as.numeric(rec(window(WG18,1983,2016))))
plot(obs, add=TRUE, col="red")            
legend("bottomright",legend=c("eqsr_fit","WG18"),bty="n",lty=1,col=c("black","red"),pch=c(NA_integer_,20))



#amalgamate 2 ouputs
rm(list=ls())
load(file=file.path(getwd(),"RData","Baseline_cases_1_2_3_4_5_6_subcases_1_2_3_SRR_SegregBlim_12-Aug-2019 15.56.RData"))
dfResults1 <- dfResults
load(file=file.path(getwd(),"RData","Baseline_cases5_subcases_2_3_SRR_SegregBlim_15-Aug-2019 15.03.RData"))
dfResults2 <- dfResults
rm(dfResults)

dfResults <- filter(dfResults1,!Scenario %in% c("5.2","5.3"))
dfResults <- dplyr::bind_rows(dfResults,dfResults2)

save(dfResults, file = file.path(getwd(),"RData","Baseline_cases_1_2_3_4_5_6_subcases_1_2_3_SRR_SegregBlim_15-Aug-2019 15.13.RData"))

