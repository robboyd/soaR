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
#' @export
#' @return
#' A list with n elements, where n is the number of species in the group. Each element contains two
#' dataframes, one with the coordinates of the presence data, and a second with the coordinates
#' of the pseudo absences. The coordinate reference system is OSGB36.


createPresAb <- function (dat, taxon, species, minYear, maxYear, nAbs, matchPres = TRUE,
          recThresh)
{

  #if (c(colnames(dat)) != c("species", "year", "eastings", "northings")) {
  #
  #  stop("Wrong column names. Must be species, year, eastings and northings.")
  #
  #}

  dat <- dat[dat$year >= minYear & dat$year <= maxYear, ]
  pres <- dat[dat$species == species, c("eastings", "northings")]
  if (nrow(pres) < recThresh) {
    warning("Number of records does not exceed recThresh")
    out <- NULL
  }
  else {
    ab <- dat[dat$species != species, c("eastings", "northings")]

    if (nrow(ab) > 7e+5) {
      ab <- ab[sample(1:nrow(ab), 7e+5),]
    }

    ab <- ab[!ab %in% pres]

    if (nrow(ab) < nrow(pres)) {
      warning("More presences than possible locations for absences. Consider lowering the number of pseudo absences.")
    }

    sampInd <- sample(1:nrow(ab), nAbs)
    if (matchPres == TRUE) {
      sampInd <- sampInd[1:nrow(pres)]
    }
    ab <- ab[sampInd, ]
    out <- list(pres, ab)
    names(out) <- c("Presence", "pseudoAbsence")
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
#' @export
#' @return
#' A list with seven elements: 1) species name; 2) number of occurrence records; 3) a model object
#' (either of type glm for "lr" or randomForest); 4) Number of folds for validation; 5) AUC score;
#' 6) predicted probabilities of occurrence in raster format; and 7) the data used to fit the model
#' (a dataframe with a column for observations and further columns for the corresponding covariates).

fitSDM <- function(species, model, envDat, spDat, k = 5, write, outPath) {

  print(species)

  ind <- which(names(spDat) == species)

  spDat <- spDat[[ind]]

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

    print(nRec)

    if (nRec < k) {

      out <- NULL

    } else {

      ab <- data.frame(val = 0, extract(x = envDat, y = spDat$pseudoAbsence))

      if (any(is.na(ab$X_Precipitation.of.Driest.Month))) {

        dropAb <- which(is.na(ab$X_Precipitation.of.Driest.Month))

        ab <- ab[-dropAb, ]

        spDat$pseudoAbsence <- spDat$pseudoAbsence[-dropAb, ]

      }

      allDat <- rbind(pres, ab)

      if (model == "lr") {

        fullMod <- glm(val ~., data = allDat, family = binomial(link = "logit"))

        type <- "response"

        index <- NULL

      } else if (model == "rf") {

        fullMod <- randomForest(x = allDat[,2:ncol(allDat)],
                                y = as.factor(allDat[,1]),
                                importance = T,
                                norm.votes = TRUE)

        type <- "prob"

        index <- 2

      } else if (model == "max") {

        fullMod <- dismo::maxent(x = envDat,
                                 p = spDat$Presence,
                                 a = spDat$PseudoAbsence)

        type = NULL

        index = NULL

      }

      pred <- predict(envDat, fullMod, type=type, index = index)

      sp::plot(pred, col = matlab.like(30))

      points(spDat$Presence, pch = "+", cex = 0.4)

      ## evaluate model performance with k fold cross validation

      folds <- kfold(pres, k)

      abFolds <- kfold(ab, k)

      allFolds <- c(folds, abFolds)

      e <- list()

      for (i in 1:k) {

        if (model != "max") {

          train <- allDat[allFolds != i,]

          test <- allDat[allFolds == i,]

        } else {

          trainPres <- spDat$Presence[folds != i,]

          trainAb <- spDat$pseudoAbsence[folds != i, ]

          testPres <- spDat$Presence[folds == i,]

          testAb <- spDat$pseudoAbsence[folds == i,]

        }

        if (model == "lr") {

          mod <- glm(val ~., data = train, family = binomial(link = "logit"))

        } else if (model == "rf") {

          mod <- randomForest(x = train[,2:ncol(train)],
                              y = as.factor(train[,1]),
                              importance = T,
                              norm.votes = TRUE)

        } else if (model == "max") {

          mod <- dismo::maxent(x = envDat, p = trainPres, a = trainAb)

        }

        if (model != "max") {

          e[[i]] <- evaluate(p=test[test$val == 1,], a=test[test$val == 0,], mod,
                             tr = seq(0,1, length.out = 200))

        } else {

          e[[i]] <- evaluate(p=testPres, a=testAb, x = envDat, mod,
                             tr = seq(0,1, length.out = 200))

        }

      }

      auc <- mean(sapply( e, function(x){slot(x, "auc")} ))

      out <- NULL

      out <- list(species, nRec, fullMod, auc, k, pred, allDat)

      names(out) <- c("Species", "Number of records", "Model","AUC","Number of folds for validation", "Predictions", "Data")

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
#' @export
#' @return
#' A dataframe with columns for species and the skill of each type of model for that species.

getSkill <- function(inPath, group) {

  allFiles <- list.files(paste0(inPath, "/", group, "/"), full.names = T, pattern = ".rdata")

  rfFiles <- grep(pattern = "_rf", allFiles, value = T)

  lrFiles <- grep(pattern = "_lr", allFiles, value = T)

  getAUC <- function(file, mod) {

    load(file)

    return(data.frame(group = group, species = out$Species, auc = out$AUC, model = mod))

  }

  rf <- purrr::map_df(.x = rfFiles, .f = getAUC, mod = "rf")

  lr <- purrr::map_df(.x = lrFiles, .f = getAUC, mod = "lr")

  skill <- rbind(rf, lr)

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
#' @export
#' @return
#' A raster layer with the AUC-weighted average probabilities of occurrence predicted by models.
#'

modelAverage <- function(inPath, outPath, skillDat, species) {

  group <- skillDat$group[skillDat$species == species][1]

  rfSkill <- skillDat$auc[skillDat$species == species & skillDat$model == "rf"]

  lrSkill <- skillDat$auc[skillDat$species == species & skillDat$model == "lr"]

  if (rfSkill >= 0.7 | lrSkill >= 0.7) {

    load(paste0(inPath, group, "/", species, "_lr.rdata"))

    lrRast <- out$Predictions

    load(paste0(inPath, group, "/", species, "_rf.rdata"))

    rfRast <- out$Predictions

    stack <- stack(rfRast, lrRast)

    print(c(rfSkill, lrSkill))

    if (rfSkill < 0.7) { rfSkill <- 0}

    if (lrSkill < 0.7) { lrSkill <- 0}

    print(c(rfSkill, lrSkill))

    ensemble <- weighted.mean(x = stack, w = c(rfSkill, lrSkill))

    raster::writeRaster(ensemble,
                        filename = paste0(outPath, "/", group, "/", species, "_ensemble"),
                        format = "ascii", overwrite = T)

  }

}


#' Create taxon-specific maps of species richness.
#'
#' Sum the predicted ensemble probabilities of occurrence for species in the chosen group.
#' @param inPath String. Location of the modelAverage outputs.
#' @param group String. Taxonomic group.
#' @param write. Logical. Should the outputs be written to file?
#' @param outPath String. Where to store the outputs if write = TRUE.
#' @export
#' @return
#' A raster layer with the AUC-weighted average probabilities of occurrence predicted by models.
#'

stackPreds <- function(inPath, group, write = TRUE, outPath) {

  files <- list.files(paste0(inPath, group, "/"),
                      full.names = T,
                      recursive = T,
                      pattern = "ensemble")

  getPreds <- function(file) {

    raster(file)

  }

  preds <- stack(lapply(X = files,
                        FUN = getPreds))
  spRich <- sum(preds)

  if (write == TRUE) {

    writeRaster(spRich, filename = paste0(outPath, group),
                format = "ascii",
                overwrite = T)

  }

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

  if (file.exists(file)) {

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
                  paste0(outPath, species),
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
