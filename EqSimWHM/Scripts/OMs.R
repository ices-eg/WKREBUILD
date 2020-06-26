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
            "SRR" = "SRR.WG19.SegReg_Blim.exterm", "RecAR" = TRUE, maxRecRes = c(3,-3),
            "BioYrs" = c(2008,2017), "BioConst" = FALSE, 
            "SelYrs" = c(2008,2017), "SelConst" = FALSE,
            "Obs" = NA,
            refPts = list("Fpa" = 0.074, "Flim" = 0.103, "Fmsy" = 0.074, "Bpa" = 1168272,
                          "Blim" = 834480, "MSYBtrigger" = 1168272, "Bloss" = 761613),
            "pBlim" = 0.05)

#WGWIDE2019 SAM assessment, IBPWHM method for reference points, stochastic bio and selection
OM2.2sam <- list("code" = "OM2.2_sam",
              "desc" = "WGWIDE19_sam",
              "IM" = NA,
              "SRR" = "SRR.WG19.SegReg_Blim.exterm", "RecAR" = TRUE, maxRecRes = c(3,-3),
              "BioYrs" = c(2008,2017), "BioConst" = FALSE, 
              "SelYrs" = c(2008,2017), "SelConst" = FALSE,
              "Obs" = NA,
              refPts = list("Fpa" = 0.115, "Flim" = 0.161, "Fmsy" = 0.115, "Bpa" = 856540,
                            "Blim" = 611814, "MSYBtrigger" = 856540, "Bloss" = 604476),
              "pBlim" = 0.05)

#WGWIDE2019 SAM assessment, IBPWHM method for reference points, stochastic bio and selection
OM2.2mac <- list("code" = "OM2.2_mac",
                 "desc" = "WGWIDE19_mackerel",
                 "IM" = NA,
                 "SRR" = "SRR.WG19.SegReg_Blim.exterm", "RecAR" = TRUE, maxRecRes = c(3,-3),
                 "BioYrs" = c(2008,2017), "BioConst" = FALSE, 
                 "SelYrs" = c(2008,2017), "SelConst" = FALSE,
                 "Obs" = NA,
                 refPts = list("Fpa" = 0.37, "Flim" = 0.46, "Fmsy" = 0.23, 
                               "Bpa" = 2500000, "Blim" = 1990000, 
                               "MSYBtrigger" = 2500000, "Bloss" = 2047874),
                 "pBlim" = 0.05)
