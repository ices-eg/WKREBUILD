SSsummarize2 <- function (biglist, sizeselfactor = "Lsel", ageselfactor = "Asel", 
          selfleet = NULL, selyr = "startyr", selgender = 1, SpawnOutputUnits = NULL, 
          lowerCI = 0.025, upperCI = 0.975) 
{
  parnames <- NULL
  dernames <- NULL
  likenames <- NULL
  allyears <- NULL
  n <- length(biglist)
  modelnames <- names(biglist)
  if (is.null(modelnames)) {
    modelnames <- paste0("model", 1:n)
  }
  for (imodel in 1:n) {
    cat(imodel,"\n")
    stats <- biglist[[imodel]]
    parnames <- union(parnames, stats$parameters$Label)
    dernames <- union(dernames, stats$derived_quants$Label)
    allyears <- union(allyears, stats$timeseries$Yr)
    likenames <- union(likenames, rownames(stats$likelihoods_used))
  }
  allyears <- sort(allyears)
  pars <- parsSD <- parphases <- as.data.frame(matrix(NA, nrow = length(parnames), 
                                                      ncol = n))
  quants <- quantsSD <- as.data.frame(matrix(NA, nrow = length(dernames), 
                                             ncol = n))
  growth <- NULL
  maxgrad <- NULL
  nsexes <- NULL
  likelihoods <- likelambdas <- as.data.frame(matrix(NA, nrow = length(likenames), 
                                                     ncol = n))
  likelihoods_by_fleet <- NULL
  likelihoods_by_tag_group <- NULL
  indices <- NULL
  sizesel <- NULL
  agesel <- NULL
  sim <- NULL
  listnames <- NULL
  npars <- NULL
  startyrs <- NULL
  endyrs <- NULL
  SPRratioLabels <- NULL
  FvalueLabels <- NULL
  sprtargs <- NULL
  btargs <- NULL
  minbthreshs <- NULL
  FleetNames <- list()
  mcmc <- list()
  warn <- FALSE
  for (imodel in 1:n) {
    stats <- biglist[[imodel]]
    listname <- names(biglist)[imodel]
    message("imodel=", imodel, "/", n)
    maxgrad <- c(maxgrad, stats$maximum_gradient_component)
    nsexes <- c(nsexes, stats$nsexes)
    startyrs <- c(startyrs, stats$startyr)
    endyrs <- c(endyrs, stats$endyr)
    sizeseltemp <- stats$sizeselex
    if (is.null(sizeselfactor)) 
      sizeselfactor <- unique(sizeseltemp$Factor)
    for (iselfactor in 1:length(sizeselfactor)) {
      seltemp_i <- sizeseltemp[sizeseltemp$Factor == sizeselfactor[iselfactor], 
                               ]
      seltemp_i$imodel <- imodel
      seltemp_i$name <- modelnames[imodel]
      if (is.null(sizesel) || (ncol(seltemp_i) == ncol(sizesel) && 
                               all(names(seltemp_i) == names(sizesel)))) {
        sizesel <- rbind(sizesel, seltemp_i)
      }
      else {
        warning("problem summarizing size selectivity due to mismatched columns ", 
                "(perhaps different bins)\n")
      }
    }
    rownames(sizesel) <- 1:nrow(sizesel)
    ageseltemp <- stats$ageselex
    if (is.null(ageselfactor)) 
      ageselfactor <- unique(ageseltemp$Factor)
    for (iselfactor in 1:length(ageselfactor)) {
      seltemp_i <- ageseltemp[ageseltemp$Factor == ageselfactor[iselfactor], 
                              ]
      seltemp_i$imodel <- imodel
      seltemp_i$name <- modelnames[imodel]
      if (is.null(agesel) || (ncol(seltemp_i) == ncol(agesel) && 
                              all(names(seltemp_i) == names(agesel)))) {
        agesel <- rbind(agesel, seltemp_i)
      }
      else {
        warning("problem summarizing age selectivity due to mismatched columns ", 
                "(perhaps different bins)\n")
      }
    }
    rownames(agesel) <- 1:nrow(agesel)
    growthtemp <- stats$growthseries
    imorphf <- ifelse(max(stats$morph_indexing$Index) == 
                        10, 3, 1)
    growthtemp <- growthtemp[growthtemp$Morph == imorphf, 
                             -(1:4)]
    growth <- cbind(growth, as.numeric(growthtemp[nrow(growthtemp), 
                                                  ]))
    liketemp <- stats$likelihoods_used
    for (irow in 1:nrow(liketemp)) {
      likelihoods[likenames == rownames(liketemp)[irow], 
                  imodel] <- liketemp$values[irow]
      likelambdas[likenames == rownames(liketemp)[irow], 
                  imodel] <- liketemp$lambdas[irow]
    }
    liketemp2 <- data.frame(model = imodel, stats$likelihoods_by_fleet)
    if (is.null(likelihoods_by_fleet) || (ncol(likelihoods_by_fleet) == 
                                          ncol(liketemp2) && all(names(likelihoods_by_fleet) == 
                                                                 names(liketemp2)))) {
      likelihoods_by_fleet <- rbind(likelihoods_by_fleet, 
                                    liketemp2)
    }
    else {
      warning("problem summarizing likelihoods by fleet due to mismatched columns")
    }
    if (!is.null(stats$likelihoods_by_tag_group)) {
      liketemp3 <- data.frame(model = imodel, stats$likelihoods_by_tag_group)
      if (is.null(likelihoods_by_tag_group) || (ncol(likelihoods_by_tag_group) == 
                                                ncol(liketemp3) && all(names(likelihoods_by_tag_group) == 
                                                                       names(liketemp3)))) {
        likelihoods_by_tag_group <- rbind(likelihoods_by_tag_group, 
                                          liketemp3)
      }
      else {
        warning("problem summarizing likelihoods by fleet due to mismatched columns")
      }
    }
    parstemp <- stats$parameters
    for (ipar in 1:nrow(parstemp)) {
      pars[parnames == parstemp$Label[ipar], imodel] <- parstemp$Value[ipar]
      parsSD[parnames == parstemp$Label[ipar], imodel] <- parstemp$Parm_StDev[ipar]
      parphases[parnames == parstemp$Label[ipar], imodel] <- parstemp$Phase[ipar]
    }
    message("  N active pars=", sum(!is.na(parstemp$Active_Cnt)))
    quantstemp <- stats$derived_quants
    for (iquant in 1:nrow(quantstemp)) {
      quants[dernames == quantstemp$Label[iquant], imodel] <- quantstemp$Value[iquant]
      quantsSD[dernames == quantstemp$Label[iquant], imodel] <- quantstemp$StdDev[iquant]
    }
    SPRratioLabels <- c(SPRratioLabels, stats$SPRratioLabel)
    FvalueLabels <- c(FvalueLabels, stats$F_report_basis)
    sprtargs <- c(sprtargs, stats$sprtarg)
    btargs <- c(btargs, stats$btarg)
    minbthreshs <- c(minbthreshs, stats$minbthresh)
    FleetNames[[imodel]] <- stats$FleetNames
    indextemp <- stats$cpue
    indextemp <- indextemp[!names(indextemp) %in% c("Area", 
                                                    "Subseas", "Month")]
    if (is.na(indextemp[[1]][1])) {
      message("no index data")
    }
    else {
      indextemp$name <- modelnames[imodel]
      indextemp$imodel <- imodel
      if (is.null(indices)) {
        indices <- rbind(indices, indextemp)
      }
      else {
        if (ncol(indextemp) == ncol(indices) && all(names(indextemp) == 
                                                    names(indices))) {
          indices <- rbind(indices, indextemp)
        }
        else {
          warning("problem summarizing indices due to mismatched columns")
        }
      }
    }
    npars <- c(npars, stats$N_estimated_parameters)
    if (!is.null(SpawnOutputUnits)) {
      if (length(SpawnOutputUnits) == 1) 
        SpawnOutputUnits <- rep(SpawnOutputUnits, n)
      if (length(SpawnOutputUnits) != n) 
        stop("'SpawnOutputUnits' should have length = 1 or", 
             n)
    }
    else {
      SpawnOutputUnits <- rep(NA, n)
    }
    if (is.na(SpawnOutputUnits[imodel])) {
      SpawnOutputUnits[imodel] <- stats$SpawnOutputUnits
    }
    if (!is.null(stats$mcmc)) {
      mcmc[[imodel]] <- stats$mcmc
    }
  }
  names(pars) <- names(parsSD) <- modelnames
  names(quants) <- names(quantsSD) <- modelnames
  names(likelihoods) <- names(likelambdas) <- modelnames
  pars$Label <- parsSD$Label <- parphases$Label <- parnames
  quants$Label <- quantsSD$Label <- dernames
  likelihoods$Label <- likelambdas$Label <- likenames
  pars$Yr <- NA
  for (ipar in 1:nrow(pars)) {
    substrings <- strsplit(as.character(pars$Label[ipar]), 
                           "_")[[1]]
    yr <- substrings[substrings %in% allyears][1]
    pars$Yr[ipar] <- ifelse(is.null(yr), NA, as.numeric(yr))
  }
  quants$Yr <- quantsSD$Yr <- NA
  for (iquant in 1:nrow(quants)) {
    substrings <- strsplit(as.character(quants$Label[iquant]), 
                           "_")[[1]]
    yr <- substrings[substrings %in% allyears][1]
    quants$Yr[iquant] <- ifelse(is.null(yr), NA, as.numeric(yr))
    quantsSD$Yr[iquant] <- ifelse(is.null(yr), NA, as.numeric(yr))
  }
  SSBrows <- grep("SSB_", quants$Label)
  SSBexclude <- c(grep("SSB_unfished", quants$Label, ignore.case = TRUE), 
                  grep("SSB_Btgt", quants$Label, ignore.case = TRUE), grep("SSB_SPR", 
                                                                           quants$Label, ignore.case = TRUE), grep("SSB_MSY", 
                                                                                                                   quants$Label, ignore.case = TRUE))
  SSBrows <- setdiff(SSBrows, SSBexclude)
  SpawnBio <- quants[SSBrows, ]
  SpawnBioSD <- quantsSD[SSBrows, ]
  minyr <- min(SpawnBio$Yr, na.rm = TRUE)
  SpawnBio$Yr[grep("SSB_Virgin", SpawnBio$Label)] <- minyr - 
    2
  SpawnBio$Yr[grep("SSB_Initial", SpawnBio$Label)] <- minyr - 
    1
  SpawnBioSD$Yr <- SpawnBio$Yr
  SpawnBio <- SpawnBio[order(SpawnBio$Yr), ]
  SpawnBioSD <- SpawnBioSD[order(SpawnBioSD$Yr), ]
  if (any(is.na(SpawnBio[3, ]))) {
    warning("Models have different start years, so SpawnBio values in VIRG & INIT yrs are shifted to correct year")
    SpawnBio$Label[1:2] <- c("SSB_Virgin*", "SSB_Initial*")
    SpawnBioSD$Label[1:2] <- c("SSB_Virgin*", "SSB_Initial*")
    for (imodel in 1:n) {
      if (is.na(SpawnBio[3, imodel])) {
        minyr <- min(SpawnBio$Yr[-(1:2)][!is.na(SpawnBio[-(1:2), 
                                                         imodel])])
        SpawnBio[SpawnBio$Yr == minyr - 2, imodel] <- SpawnBio[1, 
                                                               imodel]
        SpawnBio[SpawnBio$Yr == minyr - 1, imodel] <- SpawnBio[2, 
                                                               imodel]
        SpawnBio[1:2, imodel] <- NA
        SpawnBioSD[SpawnBio$Yr == minyr - 2, imodel] <- SpawnBioSD[1, 
                                                                   imodel]
        SpawnBioSD[SpawnBio$Yr == minyr - 1, imodel] <- SpawnBioSD[2, 
                                                                   imodel]
        SpawnBioSD[1:2, imodel] <- NA
      }
    }
  }
  SpawnBioLower <- SpawnBioUpper <- SpawnBioSD
  SpawnBioLower[, 1:n] <- qnorm(p = lowerCI, mean = as.matrix(SpawnBio[, 
                                                                       1:n]), sd = as.matrix(SpawnBioSD[, 1:n]))
  SpawnBioUpper[, 1:n] <- qnorm(p = upperCI, mean = as.matrix(SpawnBio[, 
                                                                       1:n]), sd = as.matrix(SpawnBioSD[, 1:n]))
  Bratio <- quants[grep("^Bratio_", quants$Label), ]
  BratioSD <- quantsSD[grep("^Bratio_", quantsSD$Label), ]
  BratioLower <- BratioUpper <- BratioSD
  BratioLower[, 1:n] <- qnorm(p = lowerCI, mean = as.matrix(Bratio[, 
                                                                   1:n]), sd = as.matrix(BratioSD[, 1:n]))
  BratioUpper[, 1:n] <- qnorm(p = upperCI, mean = as.matrix(Bratio[, 
                                                                   1:n]), sd = as.matrix(BratioSD[, 1:n]))
  SPRratio <- quants[grep("^SPRratio_", quants$Label), ]
  SPRratioSD <- quantsSD[grep("^SPRratio_", quantsSD$Label), 
                         ]
  SPRratioLower <- SPRratioUpper <- SPRratioSD
  SPRratioLower[, 1:n] <- qnorm(p = lowerCI, mean = as.matrix(SPRratio[, 
                                                                       1:n]), sd = as.matrix(SPRratioSD[, 1:n]))
  SPRratioUpper[, 1:n] <- qnorm(p = upperCI, mean = as.matrix(SPRratio[, 
                                                                       1:n]), sd = as.matrix(SPRratioSD[, 1:n]))
  Fvalue <- quants[grep("^F_", quants$Label), ]
  FvalueSD <- quantsSD[grep("^F_", quantsSD$Label), ]
  FvalueLower <- FvalueUpper <- FvalueSD
  FvalueLower[, 1:n] <- qnorm(p = lowerCI, mean = as.matrix(Fvalue[, 
                                                                   1:n]), sd = as.matrix(FvalueSD[, 1:n]))
  FvalueUpper[, 1:n] <- qnorm(p = upperCI, mean = as.matrix(Fvalue[, 
                                                                   1:n]), sd = as.matrix(FvalueSD[, 1:n]))
  recruits <- quants[grep("^Recr_", quants$Label), ]
  recruitsSD <- quantsSD[grep("^Recr_", quantsSD$Label), ]
  if (length(grep("Recr_Unfished", recruits$Label, ignore.case = TRUE)) > 
      0) {
    recruits <- recruits[-grep("Recr_Unfished", recruits$Label, 
                               ignore.case = TRUE), ]
    recruitsSD <- recruitsSD[-grep("Recr_Unfished", recruitsSD$Label, 
                                   ignore.case = TRUE), ]
  }
  minyr <- min(recruits$Yr, na.rm = TRUE)
  recruits$Yr[grep("Recr_Virgin", recruits$Label)] <- minyr - 
    2
  recruits$Yr[grep("Recr_Initial", recruits$Label)] <- minyr - 
    1
  recruitsSD$Yr[grep("Recr_Virgin", recruitsSD$Label)] <- minyr - 
    2
  recruitsSD$Yr[grep("Recr_Initial", recruitsSD$Label)] <- minyr - 
    1
  recruits <- recruits[order(recruits$Yr), ]
  recruitsSD <- recruitsSD[order(recruitsSD$Yr), ]
  if (any(is.na(recruits[3, ]))) {
    warning("Models have different start years, so recruits values in VIRG & INIT yrs are shifted to correct year")
    recruits$Label[1:2] <- c("Recr_Virgin*", "Recr_Initial*")
    for (imodel in 1:n) {
      if (is.na(recruits[3, imodel])) {
        minyr <- min(recruits$Yr[-(1:2)][!is.na(recruits[-(1:2), 
                                                         imodel])])
        recruits[recruits$Yr == minyr - 2, imodel] <- recruits[1, 
                                                               imodel]
        recruits[recruits$Yr == minyr - 1, imodel] <- recruits[2, 
                                                               imodel]
        recruits[1:2, imodel] <- NA
        recruitsSD[recruitsSD$Yr == minyr - 2, imodel] <- recruitsSD[1, 
                                                                     imodel]
        recruitsSD[recruitsSD$Yr == minyr - 1, imodel] <- recruitsSD[2, 
                                                                     imodel]
        recruitsSD[1:2, imodel] <- NA
      }
    }
  }
  recruitsLower <- recruitsUpper <- recruitsSD
  recruitsLower[, 1:n] <- qnorm(p = lowerCI, mean = as.matrix(recruits[, 
                                                                       1:n]), sd = as.matrix(recruitsSD[, 1:n]))
  recruitsUpper[, 1:n] <- qnorm(p = upperCI, mean = as.matrix(recruits[, 
                                                                       1:n]), sd = as.matrix(recruitsSD[, 1:n]))
  pars$recdev <- FALSE
  pars$recdev[grep("RecrDev", pars$Label)] <- TRUE
  pars$recdev[grep("InitAge", pars$Label)] <- TRUE
  pars$recdev[grep("ForeRecr", pars$Label)] <- TRUE
  InitAgeRows <- grep("InitAge", pars$Label)
  if (length(InitAgeRows) > 0) {
    temp <- unlist(strsplit(pars$Label[InitAgeRows], "InitAge_"))
    InitAgeVals <- as.numeric(temp[seq(2, length(temp), 2)])
    InitAgeYrs <- matrix(NA, nrow = length(InitAgeRows), 
                         ncol = n)
    for (imodel in 1:n) {
      modelpars <- pars[, imodel]
      devyears <- pars$Yr[!is.na(modelpars) & pars$recdev]
      if (any(!is.na(devyears))) {
        minyr <- min(devyears, na.rm = TRUE)
      }
      else {
        minyr <- NA
      }
      good <- !is.na(modelpars[InitAgeRows])
      if (!is.na(minyr) & minyr > 0 & any(good)) {
        InitAgeYrs[good, imodel] <- minyr - InitAgeVals[good]
      }
    }
    if (any(apply(InitAgeYrs, 1, max, na.rm = TRUE) - apply(InitAgeYrs, 
                                                            1, min, na.rm = TRUE) != 0)) {
      warning("years for InitAge parameters differ between models,", 
              "use InitAgeYrs matrix")
    }
    else {
      pars$Yr[InitAgeRows] <- apply(InitAgeYrs, 1, max, 
                                    na.rm = TRUE)
    }
  }
  else {
    InitAgeYrs <- NA
  }
  if (any(pars$recdev)) {
    recdevs <- pars[pars$recdev, ]
    recdevsSD <- parsSD[pars$recdev, ]
    myorder <- order(recdevs$Yr)
    recdevs <- recdevs[myorder, 1:(n + 2)]
    recdevsSD <- recdevsSD[myorder, 1:(n + 1)]
    recdevsSD$Yr <- recdevs$Yr
    recdevsLower <- recdevsUpper <- recdevsSD
    recdevsLower[, 1:n] <- qnorm(p = lowerCI, mean = as.matrix(recdevs[, 
                                                                       1:n]), sd = as.matrix(recdevsSD[, 1:n]))
    recdevsUpper[, 1:n] <- qnorm(p = upperCI, mean = as.matrix(recdevs[, 
                                                                       1:n]), sd = as.matrix(recdevsSD[, 1:n]))
  }
  else {
    recdevs <- recdevsSD <- recdevsLower <- recdevsUpper <- NULL
  }
  merge.duplicates <- function(x) {
    if (!is.null(x)) {
      if (length(unique(x$Yr)) < length(x$Yr)) {
        n <- sum(!names(x) %in% c("Label", "Yr"))
        x2 <- NULL
        for (Yr in unique(x$Yr)) {
          x.Yr <- x[which(x$Yr == Yr), ]
          if (nrow(x.Yr) == 1) {
            x2 <- rbind(x2, x.Yr)
          }
          else {
            newrow <- data.frame(t(rep(NA, n)), Label = paste0("Multiple_labels_", 
                                                               Yr), Yr = Yr)
            names(newrow) <- names(x)
            for (icol in 1:n) {
              good <- !is.na(x.Yr[, icol])
              if (sum(good) > 1) {
                warning("multiple recdevs values associated with year =", 
                        Yr)
              }
              if (sum(good) == 1) {
                newrow[, icol] <- x.Yr[good, icol]
              }
            }
            x2 <- rbind(x2, newrow)
          }
        }
      }
      else {
        x2 <- x
      }
    }
    else {
      return(x)
    }
    return(x2)
  }
  sort.fn <- function(x) {
    if (!is.null(x)) {
      return(x[order(x$Yr), ])
    }
    else {
      return()
    }
  }
  mylist <- list()
  mylist$n <- n
  mylist$npars <- npars
  mylist$modelnames <- modelnames
  mylist$maxgrad <- maxgrad
  mylist$nsexes <- nsexes
  mylist$startyrs <- startyrs
  mylist$endyrs <- endyrs
  mylist$pars <- pars
  mylist$parsSD <- parsSD
  mylist$parphases <- parphases
  mylist$quants <- quants
  mylist$quantsSD <- quantsSD
  mylist$likelihoods <- likelihoods
  mylist$likelambdas <- likelambdas
  mylist$likelihoods_by_fleet <- likelihoods_by_fleet
  mylist$likelihoods_by_tag_group <- likelihoods_by_tag_group
  mylist$SpawnBio <- sort.fn(SpawnBio)
  mylist$SpawnBioSD <- sort.fn(SpawnBioSD)
  mylist$SpawnBioLower <- sort.fn(SpawnBioLower)
  mylist$SpawnBioUpper <- sort.fn(SpawnBioUpper)
  mylist$Bratio <- sort.fn(Bratio)
  mylist$BratioSD <- sort.fn(BratioSD)
  mylist$BratioLower <- sort.fn(BratioLower)
  mylist$BratioUpper <- sort.fn(BratioUpper)
  mylist$SPRratio <- sort.fn(SPRratio)
  mylist$SPRratioSD <- sort.fn(SPRratioSD)
  mylist$SPRratioLower <- sort.fn(SPRratioLower)
  mylist$SPRratioUpper <- sort.fn(SPRratioUpper)
  mylist$SPRratioLabels <- SPRratioLabels
  mylist$Fvalue <- sort.fn(Fvalue)
  mylist$FvalueSD <- sort.fn(FvalueSD)
  mylist$FvalueLower <- sort.fn(FvalueLower)
  mylist$FvalueUpper <- sort.fn(FvalueUpper)
  mylist$FvalueLabels <- FvalueLabels
  mylist$sprtargs <- sprtargs
  mylist$btargs <- btargs
  mylist$minbthreshs <- minbthreshs
  mylist$recruits <- sort.fn(recruits)
  mylist$recruitsSD <- sort.fn(recruitsSD)
  mylist$recruitsLower <- sort.fn(recruitsLower)
  mylist$recruitsUpper <- sort.fn(recruitsUpper)
  mylist$recdevs <- merge.duplicates(sort.fn(recdevs))
  mylist$recdevsSD <- merge.duplicates(sort.fn(recdevsSD))
  mylist$recdevsLower <- merge.duplicates(sort.fn(recdevsLower))
  mylist$recdevsUpper <- merge.duplicates(sort.fn(recdevsUpper))
  mylist$growth <- growth
  mylist$sizesel <- sizesel
  mylist$agesel <- agesel
  mylist$indices <- indices
  mylist$InitAgeYrs <- InitAgeYrs
  mylist$lowerCI <- lowerCI
  mylist$upperCI <- upperCI
  mylist$SpawnOutputUnits <- SpawnOutputUnits
  mylist$FleetNames <- FleetNames
  mylist$mcmc <- mcmc
  return(invisible(mylist))
}
