#operating models

#operating model
#implementation model
#biology and fishery
#initial numbers
#SRR object, autocorrelation flag, min/max residuals (extreme trim)
#Bio vector years (number and stochastic flag)
#Sel vector years (number and stochastic flag)
#observation model (NA here)

#WGWIDE2018 update assessment, WKWIDE 2017 reference points
OM1 <- list("code" = "OM1", 
            "desc" = "WGWIDE18", 
            "IM" = NA, 
            "SRR" = "SRR.WG18.SegReg_Blim.exterm", "RecAR" = TRUE, maxRecRes = c(3,-3),
            "BioYrs" = c(2015,2017), "BioConst" = FALSE, 
            "SelYrs" = c(2015,2017), "SelConst" = FALSE,
            "Obs" = NA,
            refPts = list("Fpa" = 0.108, "Flim" = 0.151, "Fmsy" = 0.108, "Bpa" = 911587,
                          "Blim" = 661917, "MSYBtrigger" = 911587, "Bloss" = 872011),
            "pBlim" = 0.05)

#WGWIDE2019 Update assessment, IBPWHM reference points, constant bio and selection
OM2.1 <- list("code" = "OM2.1",
            "desc" = "WGWIDE19, const bio,sel",
            "IM" = NA,
            "SRR" = "SRR.WG19.SegReg_Blim.exterm", "RecAR" = TRUE, maxRecRes = c(3,-3),
            "BioYrs" = c(2008,2017), "BioConst" = TRUE, 
            "SelYrs" = c(2008,2017), "SelConst" = TRUE,
            "Obs" = NA,
            refPts = list("Fpa" = 0.074, "Flim" = 0.103, "Fmsy" = 0.074, "Bpa" = 1168272,
                          "Blim" = 834480, "MSYBtrigger" = 1168272, "Bloss" = 761613),
            "pBlim" = 0.05)

#WGWIDE2019 Update assessment, IBPWHM reference points, stochastic bio and selection
OM2.2 <- list("code" = "OM2.2",
            "desc" = "WGWIDE19",
            "IM" = NA,
            "RecAR" = TRUE, maxRecRes = c(3,-3),
            "BioYrs" = c(2008,2017), "BioConst" = FALSE, 
            "SelYrs" = c(2008,2017), "SelConst" = FALSE,
            "Obs" = NA,
            refPts = list("Fpa" = 0.074, "Flim" = 0.103, "Fmsy" = 0.074, "Bpa" = 1168272,
                          "Blim" = 834480, "MSYBtrigger" = 1168272, "Bloss" = 761613),
            "pBlim" = 0.05)

#artificially reduced recent recruitment (2014-2018 yc abundance halved in initial sim year)
OM2.2.RR <- list("code" = "OM2.2.RR",
              "desc" = "WGWIDE19 RR",
              "IM" = NA,
              "SRR" = "SRR.WG19.SegReg_Blim.exterm", "RecAR" = TRUE, maxRecRes = c(3,-3),
              "BioYrs" = c(2008,2017), "BioConst" = FALSE, 
              "SelYrs" = c(2008,2017), "SelConst" = FALSE,
              "Obs" = NA,
              refPts = list("Fpa" = 0.074, "Flim" = 0.103, "Fmsy" = 0.074, "Bpa" = 1168272,
                            "Blim" = 834480, "MSYBtrigger" = 1168272, "Bloss" = 761613),
              "pBlim" = 0.05)

#artificially reduced recent recruitment (2014-2018 yc abundance halved in initial sim year)
OM2.2.RR.5lowest <- list("code" = "OM2.2.RR.5lowest",
                 "desc" = "WGWIDE19 RR 5 lowest",
                 "IM" = NA,
                 "RecAR" = TRUE, maxRecRes = c(3,-3),
                 "BioYrs" = c(2008,2017), "BioConst" = FALSE, 
                 "SelYrs" = c(2008,2017), "SelConst" = FALSE,
                 "Obs" = NA,
                 refPts = list("Fpa" = 0.074, "Flim" = 0.103, "Fmsy" = 0.074, "Bpa" = 1168272,
                               "Blim" = 834480, "MSYBtrigger" = 1168272, "Bloss" = 761613),
                 "pBlim" = 0.05)

#WGWIDE2020 Update assessment, IBPWHM reference points, stochastic bio and selection
OM2.3 <- list("code" = "OM2.3",
              "desc" = "WGWIDE20",
              "IM" = NA,
              "SRR" = "SRR.WG20.SegReg_Blim.exterm", "RecAR" = TRUE, maxRecRes = c(3,-3),
              "BioYrs" = c(2009,2018), "BioConst" = FALSE, 
              "SelYrs" = c(2009,2018), "SelConst" = FALSE,
              "Obs" = NA,
              refPts = list("Fpa" = 0.074, "Flim" = 0.103, "Fmsy" = 0.074, "Bpa" = 1168272,
                            "Blim" = 834480, "MSYBtrigger" = 1168272, "Bloss" = 761613),
              "pBlim" = 0.05)

#WGWIDE2019 SAM assessment, IBPWHM method for reference points, stochastic bio and selection
OM2.4 <- list("code" = "OM2.4",
              "desc" = "WGWIDE19_sam",
              "IM" = NA,
              "SRR" = "SRR.WG19.SegReg_Blim.exterm", "RecAR" = TRUE, maxRecRes = c(3,-3),
              "BioYrs" = c(2008,2017), "BioConst" = FALSE, 
              "SelYrs" = c(2008,2017), "SelConst" = FALSE,
              "Obs" = NA,
              refPts = list("Fpa" = 0.115, "Flim" = 0.161, "Fmsy" = 0.115, "Bpa" = 856540,
                            "Blim" = 611814, "MSYBtrigger" = 856540, "Bloss" = 604476),
              "pBlim" = 0.05)

#WGWIDE2020 SAM assessment, IBPWHM method for reference points, stochastic bio and selection
OM2.5 <- list("code" = "OM2.5",
              "desc" = "WGWIDE20_sam",
              "IM" = NA,
              "SRR" = "SRR.WG20.SegReg_Blim.exterm", "RecAR" = TRUE, maxRecRes = c(3,-3),
              "BioYrs" = c(2008,2017), "BioConst" = FALSE, 
              "SelYrs" = c(2008,2017), "SelConst" = FALSE,
              "Obs" = NA,
              refPts = list("Fpa" = 0.115, "Flim" = 0.161, "Fmsy" = 0.115, "Bpa" = 856540,
                            "Blim" = 611814, "MSYBtrigger" = 856540, "Bloss" = 604476),
              "pBlim" = 0.05)

