#' Assess relative species distribution changes in a taxonomic group
#'
#' Telfer et al. (2002) method for calculating change in distributions between two periods.
#' @param inPath String. Directory in which the extract_records outputs were saved.
#' @param taxon String. Taxonomic group for which to run the Telfer analysis.
#' @param res Numeric. Spatial resolution (decimal degrees) to conduct the analysis at. Defaults to 0.5, i.e. half a degree.
#' @param p1 Numeric or numeric vector. Years in the first period (e.g. 1950:1990).
#' @param p2 Numeric or numeric vector. Years in the second period.
#' @export
#' @examples
#'

runTelfer <- function (inPath, taxon, res = 0.5, p1, p2) {

  path <- inPath

  files <- list.files(inPath, full.names = TRUE)

  file <- files[which(grepl(pattern=taxon, x=files, ignore.case=TRUE))]

  load(file)

  ## extract names of species recorded in each period

  spp1 <- unique(dat$species[which(dat$year %in% p1)])
  spp2 <- unique(dat$species[which(dat$year %in% p2)])

  ## set up a template raster which is used to rasterize the coordinates of records

  rast <- raster::raster(ncol = length(seq(-120, -30, res)),
                         nrow = length(seq(-60, 30, res)), xmn = -120, xmx = -30,
                         ymn = -60, ymx = 30)

  ## extract records from both periods for each species to identify commonly-sampled cells

  subSet <- function(species, period) {
    if (period ==1) {
       subDf <- dat[which(dat$year %in% p1 & dat$species == spp1[species]), c("lon", "lat")]
    } else {
      subDf <- dat[which(dat$year %in% p2 & dat$species == spp2[species]), c("lon", "lat")]
    }
  }

  subs1 <- lapply(1:length(spp1), subSet, period =1)
  subs2 <- lapply(1:length(spp2), subSet, period =2)

  rasts1 <- raster::stack(lapply(X = subs1, FUN = raster::rasterize, y = rast,
                  fun = "count"))

  rasts2 <- raster::stack(lapply(X = subs2, FUN = raster::rasterize, y = rast,
                   fun = "count"))

  print("calculating which cells were sampled in both periods. This step can take up to an hour.....")

  presAb1 <- as.logical(raster::overlay(rasts1, fun= function(x) { sum(x[!is.na(x)])}))
  presAb2 <- as.logical(raster::overlay(rasts2, fun= function(x) { sum(x[!is.na(x)])}))

  ## create a mask, i.e. a raster of cells which have been sampled in both periods

  mask <- sum(presAb1, presAb2)
  mask[mask < 2] <- NA

  nCommonCells <- length(mask[!is.na(mask)])

  ## Now filter out species with records in fewer than five grid cells (telfer 2002). First remove those
  ## with fewer than five records in total as they cannot span five grid cells.

  SubsF <- function(species) {
      subDf <- dat[which(dat$year %in% p1 & dat$species == spp1[species]), c("species","lon", "lat")]
  }

  filterSubs <- lapply(1:length(spp1), SubsF)

  x <- which(lapply(filterSubs, nrow) >= 5)

  filterSubs <- filterSubs[x]

  ## extract names of those species with > 5 records in p1

  subSpp <- unlist(lapply(filterSubs, '[[', 1, 1))

  ## extract coordinates for those species with > 5 records in p1

  subCoords <- lapply(filterSubs, function(x) x[!(names(x) %in% "species")])

  names(subCoords) <- subSpp

  ## rasterize records at chosen resolution

  p1Rasts <- lapply(subCoords, raster::rasterize, y=rast, fun="count")

  names(p1Rasts) <- subSpp

  ## establish which species have been recorded in five or more grid cells (Telfer et al. 2002)

  nCells <- unlist(
    lapply(1:length(p1Rasts), function(x) length(which(raster::getValues(!is.na(p1Rasts[[x]])))))
  )

  nCellp1 <- data.frame(as.character(subSpp), as.numeric(nCells))

  nCellp1 <- nCellp1[-which(nCellp1[,2] < 5),]

  subSpp <- nCellp1[,1] ## species recorded at > 5 sites in p1

  ## Now we have rasters of records of species' for which there is sufficient data in p1

  p1Rasts <- p1Rasts[which(names(p1Rasts) %in% subSpp)]

  ## Create identical rasters but for period 2

  SubsF2 <- function(species) {
    subDf <- dat[which(dat$year %in% p2 & dat$species == spp1[species]), c("species","lon", "lat")]
  }

  filterSubs <- lapply(1:length(spp1), SubsF2)

  names(filterSubs) <- spp1

  filterSubs <- filterSubs[which(names(filterSubs) %in% subSpp)]

  subCoords <- lapply(filterSubs, function(x) x[!(names(x) %in% "species")])

  ## rasterize records at chosen resolution

  blank <- rast
  raster::values(blank) <- NA

  rasterizeCoords <- function(species) {
    if (nrow(subCoords[[species]]) == 0) {
      y <- blank
    } else {
      y <- raster::rasterize(subCoords[[species]], y = rast, fun = "count")
    }
  }

  p2Rasts <- lapply(1:length(subCoords), rasterizeCoords)

  ## Now we have rasters for records of each species in p1 and p2. Next, filter the records
  ## to include only those sampled in both periods

  p1masked <- lapply(p1Rasts, raster::mask, mask = mask)
  p2masked <- lapply(p2Rasts, raster::mask, mask = mask)

  ## Now calculate counts in the commonly sampled period

  p1Counts <- unlist(
    lapply(1:length(p1masked), function(x) length(which(raster::getValues(!is.na(p1masked[[x]])))))
  )

  p2Counts <- unlist(
    lapply(1:length(p2masked), function(x) length(which(raster::getValues(!is.na(p2masked[[x]])))))
  )

  ## and convert to logit transformed proportions

  p1Props <- (p1Counts + 0.5) / (length(p1Counts + 1)) ## this method for calculating props can deal with zeros
  p2Props <- (p2Counts + 0.5) / (length(p2Counts + 1))

  if (length(p1Props) > 1) { ## if fewer than one species have met the criteria, do no further calculations

    p1Logit <- log(p1Props/(1 - p1Props))
    p2Logit <- log(p2Props/(1 - p2Props))

  ## regress proportions in p2 on proportions in p1

    mod <- lm(p2Logit ~ p1Logit)

    wts <- 1 / lm(abs(mod$residuals) ~ mod$fitted.values)$fitted.values^2

    mod2 <- lm(p2Logit ~ p1Logit, weights = wts)

    maxInd <- which.max(as.numeric(residuals(mod2)))
    minInd <- which.min(as.numeric(residuals(mod2)))

    outliers <- data.frame(x = c(p1Logit[maxInd], p1Logit[minInd]),
                           y = c(p2Logit[maxInd], p2Logit[minInd]))

    ggplotRegression <- function (fit) {

      ggplot2::ggplot(fit$model, ggplot2::aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
        ggplot2::geom_point() +
        ggplot2::stat_smooth(method = "lm", col = "red") +
        ggplot2::labs(title = paste("nCells =", nCommonCells)) +
        ggplot2::theme_linedraw() +
        ggplot2::geom_text(data = outliers, ggplot2::aes(x=x, y = y,
                                                         label = c(as.character(subSpp[maxInd]),
                                                                   as.character(subSpp[minInd]))), color = "red") +
        ggplot2::geom_point(data = outliers, ggplot2::aes(x=x, y = y), color = "red")

  }

    Plot <- ggplotRegression(mod2)

    index <- data.frame(subSpp, residuals(mod2), p1Counts, p2Counts)
    index <- index[order(-index[,2]),]
    colnames(index) <- c("species","index","p1Counts","p2Counts")

    out <- list(index, Plot, presAb1, presAb2, sum(presAb1, presAb2), nCommonCells)
    names(out) <- c("Outputs", "Plot", "SampledCellsP1", "SampledCellsP2", "SampledCellsP1And2", "nCommonCells")

  } else {

    warning("Insufficient number of species in this group met criteria for inclusion",
            call. = FALSE, immediate. = FALSE, noBreaks. = FALSE,
            domain = NULL)

    out <- "Insufficient number of species in this group met criteria for inclusion"
  }

  return(out)

}

