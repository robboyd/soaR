#' Format scheme records for use in soaR
#'
#' This function subsets scheme data by years, and converts UK or NI grid cells to eastings and northings
#' on the OSGB36 coordinate reference system.
#'
#' @param occ String. data.frame with occurrence records. Must have columns "YEAR", "CONCEPT" (species name) and one of "TO_GRIDREF" or "SQ_1KM".
#' @param taxa string. Taxonomic group. There must be a file in inPath that contains the taxonomic group name.
#' @param minYear Numeric. Drop all records prior to this year.
#' @param maxYear Numeric. Drop all records later than this year.
#' @param degrade Logical. Whether or not to degrade the data from all records to uniqe species/ location records.
#' @export
#' @return
#' A dataframe with four columns: species, year, eastings and northings.


formatData <- function(occ,
                       taxa,
                       minYear,
                       maxYear,
                       degrade = TRUE) {

  if (!"CONCEPT" %in% colnames(occ) | !"YEAR" %in% colnames(occ) |
      !"SQ_1KM" %in% colnames(occ) & !"TO_GRIDREF" %in% colnames(occ)) {

    stop("Wrong column names. Must include YEAR, CONCEPT and one of SQ_1KM or TO_GRIDREF.")

  }

  colnames(occ) <- toupper(colnames(occ))

  if (!"TO_GRIDREF" %in% colnames(occ)) {

    colnames(occ)[colnames(occ) == "SQ_1KM"] <- "TO_GRIDREF"

  }

  if (degrade == TRUE) {

    drop <- which(duplicated(data.frame(occ[,c("TO_GRIDREF", "CONCEPT")])))

    print(paste0("There are ", nrow(occ), " records in total"))

    occ <- occ[-drop, ]

    print(paste0("After removing records that are duplicated in terms of gridref and species, there
               are ", nrow(occ), " records"))


  }

  # filter records temporally

  occ <- occ[occ$YEAR >= minYear & occ$YEAR <= maxYear,]

  print(paste0("And after dropping records outside of the temporal extent of the analysis, there are ",
               nrow(occ), " records"))

  occ$cn <- gr_det_country(occ$TO_GRIDREF)

  #occ <- occ[,!grepl("cn", names(occ))]

  spp <- unique(occ$CONCEPT)

  coords <- gr_let2num(gridref = occ$TO_GRIDREF, centre = T)

  occ <- cbind(occ, coords)

  if (any(is.na(occ$EASTING) | is.na(occ$NORTHING))) {

    NAEast <- which(is.na(occ$EASTING))

    NANorth <- which(is.na(occ$NORTHING))

    nEast <- length(NAEast)

    nNorth <- length(NANorth)

    warning(paste("Could not obtain coordinates from grid references for", nEast, " records in terms of eastings and",
                  nNorth, " records in terms of northings. Dropping these data."))

    occ <- occ[-unique(c(NAEast, NANorth)), ]

  }

  ## If there are records from NI, reproject them on to OSGB

  if (length(unique(occ$cn)) > 1) {

    GBCRS <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs")

    NICRS <- CRS("+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=1 +x_0=200000 +y_0=250000 +ellps=airy +towgs84=482.5,-130.6,564.6,-1.042,-0.214,-0.631,8.15 +units=m +no_defs")

    NIocc <- occ[occ$cn == "OSNI",]

    GBocc <- occ[occ$cn == "OSGB",]

    NIcoords <- NIocc[, c("EASTING", "NORTHING")]

    coordinates(NIcoords) <- c("EASTING", "NORTHING")

    proj4string(NIcoords) <- NICRS

    NIcoords <- spTransform(NIcoords, GBCRS)

    NIocc[,c("EASTING", "NORTHING")] <- data.frame(NIcoords)

    occ <- rbind(GBocc, NIocc)

  }


  dat <- occ[,c("CONCEPT", "YEAR", "EASTING", "NORTHING")]

  colnames(dat) <- c("species", "year", "eastings", "northings")

  return(dat)

}


#' Extract the names of all species in a taxonomic group for which there are scheme data.
#'
#' @param inPath String. Path to outputs of formatData.
#' @param taxa string. Taxonomic group. There must be a file in inPath that contains the taxonomic group name.
#' @export
#' @return
#' A vector of species names.

getSpNames <- function(inPath, taxa) {

  load(paste0(inPath, taxa, ".rdata"))

  spp <- unique(dat$species)

}


#' Create pseudo absences and combine with the scheme occurrence data.
#'
#' This function creates pseudo absences according to the "target group" approach (Phillips 2009).
#' The aim is to generate pseudo absences with similar environmental bias as the occurrence data.
#' Pseudo absences are created for individual species, but the function works at the group level
#' because absences are inferred from the presences of non-focal species in the same group.
#'
#' @param dat String. data.frame with occurrence records for all species in the group of interest. Could be the output of formatData.
#'                    Must have column names "species", "eastings" and "northings".
#' @param species String. Species name.
#' @param minYear Numeric. Lower threshold.
#' @param minYear Numeric. Upper threshold.
#' @param nAbs Numeric. Number of pseudo absences to create.
#' @param matchpres Logical. If TRUE this overrides nAbs and creates pseudo absences in equal number to the occurrence data.
#' @param recThresh Numeric. Lower threshold number of records; species with fewer records are dropped.
#' @param screenRaster String. If not NULL then occurrence data and pseudo absences are screened against the specified raster.
#'                     If any occurrence or pseudo absence points fall in areas where the raster is NA, then they are dropped.
#'                     The function will then remove occurrence and pseudo absence points as necessary to ensure they are in equal number
#'                     if matchPres == TRUE.
#' @export
#' @return
#' A list with n elements, where n is the number of species in the group. Each element contains two
#' dataframes, one with the coordinates of the presence data, and a second with the coordinates
#' of the pseudo absences. The coordinate reference system is OSGB36.


createPresAb <- function (dat, taxon, species, minYear, maxYear, nAbs, matchPres = TRUE,
                          recThresh, screenRaster = NULL)
{
  dat <- dat[dat$year >= minYear & dat$year <= maxYear, ]
  pres <- dat[dat$species == species, c("eastings", "northings")]

  if (nrow(pres) < recThresh) {
    warning("Number of records does not exceed recThresh")
    out <- NULL
  }
  else {
    ab <- dat[dat$species != species, c("eastings",
                                        "northings")]

    if (nrow(ab) > 7e+05) {
      ab <- ab[sample(1:nrow(ab), 7e+05), ]
    }

    ab <- ab[!ab %in% pres]

    possibleAb <- nrow(ab)

    if (nrow(ab) < nrow(pres)) {
      warning("More presences than possible locations for absences. Consider lowering the number of pseudo absences.")
    }
    if (matchPres == TRUE) {
      sampInd <- sample(1:nrow(ab), nrow(pres))
    } else {
      if (nAbs <= nrow(ab)) {
        sampInd <- sample(1:nrow(ab), nAbs)
      } else {
        warning(paste0("Fewer than 10,000 locations available for pseudo absences when using the target group approach. Setting nAbs to the maximum number possible (", nrow(ab), ")."))
        sampInd <- 1:nrow(ab)
      }

    }

    ab <- ab[sampInd, ]

    out <- list(pres, ab)

    names(out) <- c("Presence", "pseudoAbsence")

    ## if screenRaster is specified, check if any presence or absence points fall outside of the raster extent (i.e. they are NA).
    ## If some data fall outside of the extent of the covariates, drop them, and drop the equivalent number of absences orpresences
    ## to ensure they are equal in number.

    if (!is.null(screenRaster)) {

      for (i in 1:nlayers(screenRaster)) {

        presDrop <- raster::extract(screenRaster[[i]], out$Presence)

        abDrop <- raster::extract(screenRaster[[i]], out$pseudoAbsence)

        if (any(is.na(presDrop))) out$Presence <- out$Presence[-which(is.na(presDrop)), ]

        if (any(is.na(abDrop))) out$pseudoAbsence <- out$pseudoAbsence[-which(is.na(abDrop)), ]

      }

      if (matchPres == TRUE) {

        if (nrow(out$Presence) > nrow(out$pseudoAbsence)) {

          out$Presence <- out$Presence[1:nrow(out$pseudoAbsence), ]

        } else if (nrow(out$Presence) < nrow(out$pseudoAbsence)) {

          out$pseudoAbsence <- out$pseudoAbsence[1:nrow(out$Presence), ]

        }

      }

    }

  }

  return(out)

}




#' Fit species distribution models.
#'
#' Fit logistic regression, random forest or Maxent models using the outputs of createPresAb and some covariates.
#'
#' @param species String. Species name (see getSpNames).
#' @param model string. One of "lr", "rf" or "max" for logistic regression, random forest or Maxent.
#' @param envDat String. rasterStack object with n layers each representing a covariate.
#' @param spDat String. Outputs of createPresAb for the chosen species.
#' @param k Numeric. Number of folds for cross validation. Defaults to 5.
#' @param write Logical. If TRUE writes results to file in outPath.
#' @param outPath String. Where to store the outputs if write = TRUE.
#' @param predict Logical. Whether or not to predict the fitted models on to the covariate rasters. Prediction
#'                         significantly increases execution time but will probably need to be done at some point
#'                         anyway.
#' @param plot Logical. Whether or not to show plots of predicted probabilities of occurrence if predict = TRUE.
#' @export
#' @return
#' A list with seven elements: 1) species name; 2) number of occurrence records; 3) a model object
#' (either of type glm for "lr" or randomForest); 4) Number of folds for validation; 5) AUC score;
#' 6) predicted probabilities of occurrence in raster format; and 7) the data used to fit the model
#' (a dataframe with a column for observations and further columns for the corresponding covariates).

fitSDM <- function(species, model, envDat, spDat, k = 5, write, outPath, predict = TRUE, plot = TRUE) {

  if (predict == FALSE & model == "lrReg") warning("When model = lrReg and predict = TRUE, fit SDM will return k AUC scores and a mean based on those scores.
  However, in some cases lrReg reduces all covariates to zero producing an intercept-only model (ie. the mean).
  These models introduce substantial noise and we recommend they are dropped at a later stage.")

  print(paste("species:", species))

  ind <- which(names(spDat) == species)

  spDat <- spDat[[ind]]

  ## covariates needed in matrix format to predict using lrReg

  if (model == "lrReg" & predict == TRUE) {

    covsMat <- as.matrix(rasterToPoints(covs))

  }

  if (is.null(spDat)) {

    out <- NULL

  } else {

    pres <- data.frame(val = 1, extract(x = envDat, y = spDat$Presence))

    if (any(is.na(pres$X_Precipitation.of.Driest.Month))) {

      dropPres <- which(is.na(pres$X_Precipitation.of.Driest.Month))

      print(paste("Dropping", length(dropPres), "records because they fall outside the extent of the covariate data"))

      pres <- pres[-dropPres, ]

      spDat$Presence <- spDat$Presence[-dropPres, ]

    }

    nRec <- nrow(pres)

    print(paste("Occurrence records:", nRec))

    if (nRec < k) {

      out <- NULL

    } else {

      ab <- data.frame(val = 0, extract(x = envDat, y = spDat$pseudoAbsence))

      if (any(is.na(ab$X_Precipitation.of.Driest.Month))) {

        dropAb <- which(is.na(ab$X_Precipitation.of.Driest.Month))

        ab <- ab[-dropAb, ]

        spDat$pseudoAbsence <- spDat$pseudoAbsence[-dropAb, ]

      }

      ## set weights for logistic regression if prevalence is not 0.5

      if (model != "lrReg" & nRec != nrow(ab)) { warning("Prevalence is not 0.5 and no weights are applied to account for this. Currently weights are only applied where model = lrReg.")}

      if (model == "lrReg" & nRec != nrow(ab)) {

        print("Prevalence is not 0.5. Weighting absences to simulate a prevalence of 0.5")

        nAb <- nrow(ab)

        prop <- nRec / nAb

        print(paste("Absence weighting:", prop))

      }

      allDat <- rbind(pres, ab)

      folds <- kfold(pres, k)

      abFolds <- kfold(ab, k)

      allFolds <- c(folds, abFolds)

      e <- list()

      for (i in 1:k) {

        if (model != "max") {

          train <- allDat[allFolds != i,]

          if (model == "lrReg" & nRec != nrow(ab)) weights <- c(rep(1, length(train$val[train$val == 1])), rep(prop, length(train$val[train$val == 0])))

          test <- allDat[allFolds == i,]

        } else {

          trainPres <- spDat$Presence[folds != i,]

          trainAb <- spDat$pseudoAbsence[folds != i, ]

          testPres <- spDat$Presence[folds == i,]

          testAb <- spDat$pseudoAbsence[folds == i,]

        }

        if (model == "lr") {

          assign(paste0("mod", i), glm(val ~., data = train, family = binomial(link = "logit")))

          if (predict == TRUE) {

            assign(paste0("pred", i), predict(envDat, get(paste0("mod", i)), type= "response"))

          }

        } else if (model == "rf") {

          assign(paste0("mod", i), randomForest(x = train[,2:ncol(train)],
                                                y = as.factor(train[,1]),
                                                importance = T,
                                                norm.votes = TRUE))

          if (predict == TRUE) {

            assign(paste0("pred", i), predict(envDat, get(paste0("mod", i)), type = "prob", index = 2))

          }

        } else if (model == "max") {

          assign(paste0("mod", i), dismo::maxent(x = envDat, p = trainPres, a = trainAb))

          if (predict == TRUE) {

            assign(paste0("pred", i), predict(envDat, get(paste0("mod", i))))

          }

        } else if (model == "lrReg") {

          assign(paste0("mod", i), glmnet::cv.glmnet(x = as.matrix(train[, 2:ncol(train)]),
                                                     y = train[, 1],
                                                     family = "binomial",
                                                     nfolds = 3,
                                                     weights = weights))

          if (predict == TRUE) {

            pred <- predict(get(paste0("mod", i)), covsMat[, 3:ncol(covsMat)], type = "response")

            pred <- as.matrix(cbind(covsMat[, 1:2], pred))

            if (any(is.na(pred[, 3]))) pred <- pred[-which(is.na(pred[,3])), ]

            assign(paste0("pred", i), rasterize(pred[, 1:2], covs[[1]], field = pred[, 3]))

          }

        }

        ## extract AUC

        if (model == "lrReg") {

          ## extract predictions for test data

          coords <- rbind(spDat$Presence, spDat$pseudoAbsence)

          coords <- coords[allFolds == i, ]

          preds <- extract(get(paste0("pred", i)), coords)

          DATA <- data.frame(index = 1:length(preds),
                             obs = test$val,
                             pred = preds)

          e[[i]] <- auc(DATA, st.dev = F)

        } else if (model != "max") {

          e[[i]] <- evaluate(p=test[test$val == 1,], a=test[test$val == 0,], get(paste0("mod", i)),
                             tr = seq(0,1, length.out = 200))

        } else {

          e[[i]] <- evaluate(p=testPres, a=testAb, x = envDat, get(paste0("mod", i)),
                             tr = seq(0,1, length.out = 200))

        }

      }

      if (model != "lrReg") {

        AUC <- sapply( e, function(x){slot(x, "auc")} )

      } else {

        AUC <- unlist(e)

      }

      meanAUC <- mean(AUC)

      if (predict == TRUE) {

        pred <- stack(lapply(1:k,
                             function(x) {get(paste0("pred", x))}))

        ## for lrReg models all coefficients may be reduced to zero giving an intercept-only model which predicts 0.5 everywhere
        ## the below drops intercept-only models because they add needless noise to the predictions

        if (model == "lrReg") {

          uniqueVals <-lapply(1:k,
                              function(x) { length(cellStats(pred[[x]], unique))})

          drop <- which(uniqueVals <= 2) ## i.e. the mean and NA

          if (length(drop) == k) {

            print("No non-intercept models; lrReg will be omitted for this species")

          }

          if (any(drop)) {

            AUC[drop] <- NA

            print(paste("Dropping", length(drop), "intercept-only model(s). Intercept-only models are given an AUC value of NA so they can be identified.
                        Where 1:(k-1) models are intercept only, only the non-intercept models are included in the final average. Where all models are intercept-only,
                        their predictions are returned but should not be used. modelAverage() will ignore intercept-only models"))

          }

        } else { drop <- NULL}

        if (length(drop) == k | length(drop) == 0) {

          meanPred <- mean(pred) # where all models are intercept-only, takes the mean to avoid errors later but AUC scores are NA which means they are dropped for final ensembles

        } else {

          meanPred <- mean(pred[[-drop]])

        }

        print("k fold AUC scores:")

        print(AUC)

        meanAUC <- mean(AUC, na.rm = T)

        print(paste("mean AUC:", meanAUC))

      } else {

        pred <- NULL

        meanPred <- NULL

      }

      if (plot == TRUE) {

        sp::plot(pred)

        #points(spDat$Presence, pch = "+")

      }

      mods <- lapply(1:k,
                     function(x) {get(paste0("mod", x))})

      out <- NULL

      if (model == "lrReg" & predict == TRUE) k <- k - length(drop) ## if some models were intercept-only then recalculate k (number of models used)

      out <- list(species, nRec, k, mods, AUC, meanAUC, pred, meanPred, allDat)

      names(out) <- c("species", "nRecords", "nValidationFolds","modelObjects","AUCScores", "meanAUC", "predictions", "meanPrediction", "data")

      if (write == TRUE) {

        save(out, file = paste0(outPath, species, "_", model, ".rdata"))

      }

    }

  }

}


#' Extract skill of SDMs.
#'
#' Extracts the AUC score for each species/ model combo in a group.
#'
#' @param inPath String. Location of the fitSDM outputs.
#' @param group String. Taxonomic group. Must partially match a filename in inPath.
#' @param write Logical. Whether or not to write outputs to file.
#' @param outPath String. Where to save outputs if write = TRUE.
#' @param models String or character vector. Which models to include in skill assessment. Options are
#'               "lr", "rf" and "max".
#' @export
#' @return
#' A dataframe with columns for species and the skill of each type of model for that species.

getSkill <- function(inPath, group, write, outPath, models) {

  print(group)

  allFiles <- list.files(paste0(inPath, "/", group, "/"), full.names = T, pattern = ".rdata")

  for (i in models) {

    assign(paste0(i, "Files"), grep(pattern = paste0("_", i), allFiles, value = T))

  }

  getAUC <- function(file, mod) {

    load(file)

    return(data.frame(group = group, species = out$species, auc = out$meanAUC, model = mod))

  }

  skill <- lapply(models,
                  function(x) { y <- purrr::map_df(.x = get(paste0(x, "Files")), .f = getAUC, mod = x)})

  skill <- do.call("rbind", skill)

  if (write == TRUE) {

    write.csv(skill, paste0(outPath, group, ".csv"),
              row.names = FALSE)

  }

  return(skill)

}


#' Create ensemble predictions.
#'
#' AUC-weighted average predictions from the models fitted for a species using fitSDM.
#'
#' @param inPath String. Location of the fitSDM outputs.
#' @param outPath String. Where to store ensemble models.
#' @param skillDat. String. getSkill outputs for the chosen group.
#' @param species String. Species name.
#' @param models String or character vector. SDM algorithms to be included in model averaging. Options are
#'               "lr" for logistic regression", "rf" for random forests and "max" for maxent.
#' @param plot Logical. Whether or not to plot the predictions from each algorithm and the final ensemble.
#' @param skillThresh Numeric. Algorithms with AUC scores below this value are dropped from the ensemble.
#' @export
#' @return
#' A raster layer with the AUC-weighted average probabilities of occurrence predicted by models.
#'

modelAverage <- function(inPath, outPath, skillDat, species, models, plot = TRUE, skillThresh) {

  group <- skillDat$group[skillDat$species == species][1]

  print(paste("group:", as.character(group)))

  print(paste("species:", as.character(species)))

  l <- NULL

  for (i in models) {

    assign(paste0(i, "Skill"), skillDat$auc[skillDat$species == species & skillDat$model == i])

    l <- c(l, get(paste0(i, "Skill")))

  }

  if (length(l) == length(models) & length(l[which(l > skillThresh)]) > 0) {

    for (i in models) {

      print(paste(i, "AUC =", get(paste0(i, "Skill"))))

      load(paste0(inPath, group, "/", species, "_", i, ".rdata"))

      assign(paste0(i, "Rast"), out$meanPrediction)

      ## drop models with AUC a threshold or NA. If NA it means that only intercept-only models were produced by lrReg

      if (get(paste0(i, "Skill")) <= skillThresh | is.na(get(paste0(i, "Skill")))) assign(paste0(i, "Skill"), 0)

    }

    stack <- stack(lapply(models,
                          function(x) {get(paste0(x, "Rast"))}))


    weights <- lapply(models,
                      function(x) {get(paste0(x, "Skill"))})

    weights <- do.call("c", weights)

    print(paste("Model weights:", weights))

    ensemble <- weighted.mean(x = stack, w = weights)

    if (plot == TRUE) {

      par(mfrow=c(1, (length(models) +1)))

      plot(stack)

      plot(ensemble, main = "ensemble")

    }

    raster::writeRaster(ensemble,
                        filename = paste0(outPath, "/", group, "/", species, "_ensemble"),
                        format = "ascii", overwrite = T)


  } else {

    if (length(l) != length(models)) warning("All models not fitted for this species")

    return(species)

  }

}


#' Create taxon-specific maps of species richness.
#'
#' Sum the predicted ensemble probabilities of occurrence for species in the chosen group.
#' @param inPath String. Location of the modelAverage outputs.
#' @param group String. Taxonomic group.
#' @param write. Logical. Should the outputs be written to file?
#' @param outPath String. Where to store the outputs if write = TRUE.
#' @param species String or character vector. List of species to include in the stack. All species
#'                not in this list will be excluded. Defaults to NULL in which case all species in group
#'                are included.
#' @export
#' @return
#' A raster layer with the AUC-weighted average probabilities of occurrence predicted by models.
#'

stackPreds <- function(inPath, group, write = TRUE, outPath, species = NULL) {

  names <- list.files(paste0(inPath, group, "/"),
                      pattern = "ensemble")

  names <- gsub("_ensemble.asc", "", names)

  if (!is.null(species)) {

    inds <- which(names %in% species)

  } else {

    inds <- 1:length(names)

  }

  if (length(inds) > 0) {

    out <- data.frame(group = group,
                      species = names[inds])

  } else {

    out <- NULL

  }


  files <- list.files(paste0(inPath, group, "/"),
                      full.names = T,
                      recursive = T,
                      pattern = "ensemble")

  files <- files[inds]

  getPreds <- function(file) {

    raster(file)

  }

  if (length(files) > 0) {

    spRich <- stack(lapply(X = files,
                           FUN = getPreds))

    if (nlayers(spRich) > 1) {

      spRich <- sum(spRich)

    } else {

      spRich <- spRich[[1]]

    }


    if (write == TRUE) {

      writeRaster(spRich, filename = paste0(outPath, group),
                  format = "ascii",
                  overwrite = T)

    }

  }

  return(out)

}



#' Create binary presence/ absence predictions from probabilities of occurrence.
#'
#' This function optimises a threshold probability above which a grid cell is classed as a presence,
#' and below which it is classed as an absence. Optimisation is achieved
#' by maximising the proportion of data correctly classified.
#'
#' @param inPath String. Path to the ensemble model outputs.
#' @param group string. Taxonomic group. There must be a file in inPath.
#' @param species String. There must be a file in inPath called /<group>/<species>.asc.
#' @param outPath String. Where to store the outputs.
#' @param map Logical. If true writes an ascii file to outPath.
#' @param presAb String. Presence/ pseudo absence data created by createPresAb for the chosen group.
#' @export
#' @return
#' A dataframe with four columns: Threshold probability of occurrence, species and group.
#' Optionally a raster can be written to file.

optOccThresh <- function(inPath, group, species, outPath, map, presAb) {

  print(species)

  if (!dir.exists(outPath)) {

    dir.create(outPath)

  }

  file <- paste0(inPath, group, "/", species, "_ensemble.asc")

  ind <- which(names(presAb) == species)

  if (file.exists(file) & !is.null(presAb[[ind]])) {

    out <- raster::raster(file)

    dat <- presAb[[which(names(presAb) == species)]]

    nRecs <- nrow(dat$Presence)

    presProb <- raster::extract(out, dat$Presence)

    abProb <- raster::extract(out, dat$pseudoAbsence)

    obs <- c(rep(1, nrow(dat$Presence)), rep(0, nrow(dat$pseudoAbsence)))

    preds <- c(presProb, abProb)

    DATA <- data.frame(id = 1:length(obs),
                       obs = obs,
                       pred = preds)

    thresh <- seq(0,1, length.out = 50)

    opt <- PresenceAbsence::optimal.thresholds(DATA, thresh, opt.methods = 3, na.rm = T)

    opt <- as.numeric(opt[2])

    auc <- PresenceAbsence::auc(DATA, st.dev = F, na.rm = T)

    confMat <- PresenceAbsence::cmx(DATA, threshold = opt, na.rm = T)

    sens <- PresenceAbsence::sensitivity(confMat, st.dev = F)

    spec <- PresenceAbsence::specificity(confMat, st.dev = F)

    tss <- (sens + spec) - 1

    kappa <- PresenceAbsence::Kappa(confMat, st.dev = F)

    sumStats <- data.frame(group = group,
                           species = species,
                           nRecs = nRecs,
                           thresh = opt,
                           auc = auc,
                           sensitivity = sens,
                           specificity = spec,
                           tss = tss,
                           kappa = kappa)

    if (map == TRUE) {

      binPred <- out

      binPred[binPred < opt] <- 0

      binPred <- as.logical(binPred)

      writeRaster(binPred,
                  paste0(outPath, species, "_binary"),
                  format = "ascii",
                  overwrite = T)

      par(mfrow = c(1,2), mar = c(2,2,2,2))

      plot(out, col = matlab.like(50), axes = F)

      plot(binPred, axes = F, legend = F)

    }

  } else {

    sumStats <- NULL

  }

  return(sumStats)

}

