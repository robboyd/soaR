#' Calculate the geographical overlap of sampling in two periods for a given taxa
#'
#' Calculates the number of cells, at a specified resolution, sampled in both of two user-chosen time periods.
#' @param inPath String. Directory in which the extract_records outputs were saved.
#' @param taxon String. Taxonomic group of interest.
#' @param res Numeric. Spatial resolution (decimal degrees) to conduct the analysis at. Defaults to 0.5, i.e. half a degree.
#' @param p1 Numeric or numeric vector. Years in the first period (e.g. 1950:1990).
#' @param p2 Numeric or numeric vector. Years in the second period.
#' @export
#' @examples
#'

samplingOverlap <- function (inPath, taxon, res = 0.5, p1, p2) {

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

  return(nCommonCells)

}
