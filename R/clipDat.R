#' \code{clipDat} - Extract a geographic subset of your occurrence records.
#'
#' Extracts all records that fall within a provided polygon (e.g. a country boundary).
#'
#' @param dat String. Outputs from \code{formatData}.
#' @param shp String. A polygon which will be used to clip the data.
#' @export
#' @examples
#' @return a data.frame with the new subset of occurrence records.

clipDat <- function(dat, shp) {

  pts <- SpatialPoints(coords = data.frame(eastings = dat$eastings,
                                           northings = dat$northings),
                       proj4string = crs(shp))


  keep <- which(!is.na((pts %over% scot)[,1]))

  dat <- dat[keep, ]

  return(dat)

}
