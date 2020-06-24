#management procedure (MP) objects

#decision model + constraints
#estimation model

#baseline, no harvest rule, no IAV control, no minimum TAC, no assessment/advice error
MP1.0 <- list("code" = "MP1.0",
              "desc" = "NoHCR_Test",
              "xlab" = "Base",
              "HCRName" = "None",
              "F_target" = c(0,0.05,0.074,0.1,0.108,0.2,0.3),
              "B_trigger" = NA,
              "minTAC" = NA,
              "maxTAC" = NA,
              "TAC_IAV" = NA,
              "Obs" = list("cvF" = 0, "phiF" = 0, "cvSSB" = 0, "phiSSB" = 0))

#minimum TAC test
MP1.1 <- list("code" = "MP1.1",
              "desc" = "NoHCR_MinTACTest",
              "xlab" = "80kt Min",
              "HCRName" = "None",
              "F_target" = c(0,0.05,0.074,0.1,0.108,0.2,0.3),
              "B_trigger" = NA,
              "minTAC" = 80000,
              "maxTAC" = NA,
              "TAC_IAV" = NA,
              "Obs" = list("cvF" = 0, "phiF" = 0, "cvSSB" = 0, "phiSSB" = 0))

#maximum TAC test
MP1.2 <- list("code" = "MP1.2",
              "desc" = "NoHCR_MaxTACTest",
              "xlab" = "150kt Max",
              "HCRName" = "None",
              "F_target" = c(0,0.05,0.074,0.1,0.108,0.2,0.3),
              "B_trigger" = NA,
              "minTAC" = NA,
              "maxTAC" = 150000,
              "TAC_IAV" = NA,
              "Obs" = list("cvF" = 0, "phiF" = 0, "cvSSB" = 0, "phiSSB" = 0))

#min/max TAC in combination
MP1.3 <- list("code" = "MP1.3",
              "desc" = "NoHCR_MinMaxTACTest",
              "xlab" = "80/150",
              "HCRName" = "None",
              "F_target" = c(0,0.05,0.074,0.1,0.108,0.2,0.3),
              "B_trigger" = NA,
              "minTAC" = 80000,
              "maxTAC" = 150000,
              "TAC_IAV" = NA,
              "Obs" = list("cvF" = 0, "phiF" = 0, "cvSSB" = 0, "phiSSB" = 0))

#baseline, no harvest rule, no IAV control, no minimum TAC plus assessment/advice error
MP1.4 <- list("code" = "MP1.4",
              "desc" = "NoHCR_ObsErrTest",
              "xlab" = "Inc Err",
              "HCRName" = "None",
              "F_target" = c(0,0.05,0.074,0.1,0.108,0.2,0.3),
              "B_trigger" = NA,
              "minTAC" = NA,
              "maxTAC" = NA,
              "TAC_IAV" = NA,
              "Obs" = list("cvF" = 0.22, "phiF" = 0.03, "cvSSB" = 0.36, "phiSSB" = 0.51))

#baseline, no harvest rule, 20% TAC change limitation, no min/max TAC, no assessment/advice error
MP1.5 <- list("code" = "MP1.5", 
             "desc" = "NoHCR_20PctIAVTest", 
             "xlab" = "0.2 IAV",
             "HCRName" = "None",
             "F_target" = c(0,0.05,0.074,0.1,0.108,0.2,0.3),
             "B_trigger" = NA,
             "minTAC" = NA,
             "maxTAC" = NA,
             "TAC_IAV" = c(0.2,0.2),
             "Obs" = list("cvF" = 0, "phiF" = 0, "cvSSB" = 0, "phiSSB" = 0))

#baseline, no harvest rule, 30% TAC change limitation, no min/max TAC, no assessment/advice error
MP1.6 <- list("code" = "MP1.6", 
              "desc" = "NoHCR_30PctIAVTest", 
              "xlab" = "0.3 IAV",
              "HCRName" = "None",
              "F_target" = c(0,0.05,0.074,0.1,0.108,0.2,0.3),
              "B_trigger" = NA,
              "minTAC" = NA,
              "maxTAC" = NA,
              "TAC_IAV" = c(0.3,0.3),
              "Obs" = list("cvF" = 0, "phiF" = 0, "cvSSB" = 0, "phiSSB" = 0))

#baseline, no harvest rule, 10% TAC change limitation, no min/max TAC, no assessment/advice error
MP1.7 <- list("code" = "MP1.7", 
              "desc" = "NoHCR_10PctIAVTest", 
              "xlab" = "0.1 IAV",
              "HCRName" = "None",
              "F_target" = c(0,0.05,0.074,0.1,0.108,0.2,0.3),
              "B_trigger" = NA,
              "minTAC" = NA,
              "maxTAC" = NA,
              "TAC_IAV" = c(0.1,0.1),
              "Obs" = list("cvF" = 0, "phiF" = 0, "cvSSB" = 0, "phiSSB" = 0))

#baseline, no harvest rule, 10%/20% TAC change limitation, no min/max TAC, no assessment/advice error
MP1.8 <- list("code" = "MP1.8", 
              "desc" = "NoHCR_10_20PctIAVTest", 
              "xlab" = "10_20 IAV",
              "HCRName" = "None",
              "F_target" = c(0,0.05,0.074,0.1,0.108,0.2,0.3),
              "B_trigger" = NA,
              "minTAC" = NA,
              "maxTAC" = NA,
              "TAC_IAV" = c(0.1,0.2),
              "Obs" = list("cvF" = 0, "phiF" = 0, "cvSSB" = 0, "phiSSB" = 0))

#baseline, no harvest rule, TAC change limitation - no limit on increase, 10% limit on decrease, no min/max TAC, no assessment/advice error
MP1.9 <- list("code" = "MP1.9", 
              "desc" = "NoHCR_100_10PctIAVTest", 
              "xlab" = "100_10 IAV",
              "HCRName" = "None",
              "F_target" = c(0,0.05,0.074,0.1,0.108,0.2,0.3),
              "B_trigger" = NA,
              "minTAC" = NA,
              "maxTAC" = NA,
              "TAC_IAV" = c(1,0.1),   #increase/decrease limits
              "Obs" = list("cvF" = 0, "phiF" = 0, "cvSSB" = 0, "phiSSB" = 0))

# ==========================================================================================================================================

# #ICES HCR, no IAV control, no minimum TAC, no assessment/advice error
# MP2.0 <- list("code" = "MP2.0", 
#               "desc" = "ICESHCR", 
#               "HCRName" = "ICES",
#               "F_target" = c(0,0.05,0.074,0.1,0.108,0.2,0.3),
#               "B_trigger" = "MSYBtrigger",
#               "minTAC" = NA, 
#               "maxTAC" = NA,
#               "TAC_IAV" = NA,
#               "Obs" = list("cvF" = 0, "phiF" = 0, "cvSSB" = 0, "phiSSB" = 0))
# 
# #num iterations investigation
# MP2.0_10000 <- list("code" = "MP2.0_10000", 
#                   "desc" = "ICESHCR", 
#                   "HCRName" = "ICES",
#                   "F_target" = c(0,0.05,0.074,0.1,0.108,0.2,0.3),
#                   "B_trigger" = "MSYBtrigger",
#                   "minTAC" = NA, 
#                   "maxTAC" = NA,
#                   "TAC_IAV" = NA,
#                   "Obs" = list("cvF" = 0, "phiF" = 0, "cvSSB" = 0, "phiSSB" = 0))

#ICES HCR, no IAV control, no minimum TAC, with assessment/advice error
MP2.1 <- list("code" = "MP2.1",
              "desc" = "ICESHCR",
              "xlab" = "ICES AR",
              "HCRName" = "ICES",
              "F_target" = c(0,0.05,0.074,0.1,0.115,0.2,0.3),
              "B_trigger" = "MSYBtrigger",
              "minTAC" = NA,
              "maxTAC" = NA,
              "TAC_IAV" = NA,
              "Obs" = list("cvF" = 0.3, "phiF" = 0.3, "cvSSB" = 0, "phiSSB" = 0))

 
 
# ==========================================================================================================================================

# #NFD rule No further decline
# MP3.0 <- list("code" = "MP3.0", 
#               "desc" = "NFD", 
#               "HCRName" = "NFD",
#               "F_target" = c(0,0.05,0.074,0.1,0.108,0.2,0.3),
#               "B_trigger" = "MSYBtrigger",
#               "minTAC" = NA,
#               "maxTAC" = NA,
#               "TAC_IAV" = NA,
#               "Obs" = list("cvF" = 0, "phiF" = 0, "cvSSB" = 0, "phiSSB" = 0))
# 

# #NFD rule No further decline + IAV constraint
# #issues here - investigate
# MP3.1 <- list("code" = "MP3.1", 
#               "desc" = "NFD_IAV", 
#               "HCRName" = "NFD",
#               "F_target" = c(0,0.05,0.074,0.1,0.108,0.2,0.3),
#               "B_trigger" = "MSYBtrigger",
#               "minTAC" = NA,
#               "maxTAC" = NA,
#               "TAC_IAV" = 0.2,
#               "Obs" = list("cvF" = 0, "phiF" = 0, "cvSSB" = 0, "phiSSB" = 0))
