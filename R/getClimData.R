#' Process climate data for use in modelling species' distributions
#'
#' Create environmental layers in ascii for a specified region and spatial resolution
#'
#' @param res Numeric. Spatial resolution in decimal degrees (defaults to 0.5)
#' @param xMin Numeric. Min x coordinate.
#' @param xMax Numeric. Max x coordinate.
#' @param yMin Numeric. Min y coordinate.
#' @param yMax Numeric. Max y coordinate.
#' @param inPath String. Path to climate data.
#' @param outPath String. If write - TRUE, location to save processed climate data.
#' @param bioClim Logical. If true, getClimData will automatically name the data layers.
#' @param outPath Logical. If true the climate data are written as .asc files to outPath.
#' @export
#' @examples
#'
getClimData <- function(res = 0.5, xMin, xMax, yMin, yMax, inPath, outPath, bioClim, write) {

  #files <- list.files("F:/SURPASS/data/worldClim/bio/",full.names = TRUE)

  names <- c("Annual Mean Temperature",
             "Mean Diurnal Range (Mean of monthly (max temp - min temp))",
             "Isothermality (BIO2/BIO7)",
             "Temperature Seasonality (standard deviation x100)",
             "Max Temperature of Warmest Month",
             "Min Temperature of Coldest Month",
             "Temperature Annual Range (BIO5-BIO6)",
             "Mean Temperature of Wettest Quarter",
             "Mean Temperature of Driest Quarter",
             "Mean Temperature of Warmest Quarter",
             "Mean Temperature of Coldest Quarter",
             "Annual Precipitation",
             "Precipitation of Wettest Month",
             "Precipitation of Driest Month",
             "Precipitation Seasonality (Coefficient of Variation)",
             "Precipitation of Wettest Quarter",
             "Precipitation of Driest Quarter",
             "Precipitation of Warmest Quarter",
             "Precipitation of Coldest Quarter")

  dat <- raster::stack(files)

  if (bioClim == TRUE) {
    names(dat) <- names
  }

  ext <- raster::extent(xMin, xMax, yMin, yMax)

  dat <- raster::crop(dat, ext)

  mult <- res / raster::res(dat)[1]

  dat <- raster::aggregate(dat, fact = mult)

  whichVars <- usdm::vifcor(dat, th = 0.7)

  if (write == TRUE) {

    raster::writeRaster(x = dat,
                        file = paste0(outPath, "clim_data"),
                        format = "ascii",
                        bylayer = TRUE)

  }

  return(list(dat, whichVars))

}


