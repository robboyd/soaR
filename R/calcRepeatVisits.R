#' Calculate the number of repeat visits to grid cells across a taxonomic group
#'
#' This function can be used to assess, for a chosen taxonomic group, the number of records per grid cell.
#' @param inPath String. Location of the data obtained using extract_records
#' @param taxon string. Taxonomic group. Must match a file name (without "_raw_data.rdata") in inPath.
#' @param res Numeric. Spatial resolution for the analysis in decimal degrees.
#' @param start Numeric. Earliest year cutoff for analysis.
#' @param end Numeric. Final year cutoff for analysis
#' @param output String. One of: "map" to get a map of the no. records per cell, "hist" to get a histrogram of repeat visits, or other to get the raw values (no. repeat visits per cell).
#' @export
#' @examples

calcRepeatVisits <- function(inPath, taxon, res, start, end, output) {
  
  load(paste0(inPath, taxon, "_raw_data.rdata"))
  
  rast <- raster::raster(ncol=length(seq(-120,-30,res)),
                         nrow=length(seq(-60,30,res)),
                         xmn=-120,
                         xmx=-30,
                         ymn=-60,
                         ymx=30)
  
  dat <- dat[which(dat$year >= start & dat$year <= end), ]
  
  data <- data.frame(dat$lon, dat$lat)
  
  rDat <- raster::rasterize(x = data, y = rast, fun = "count")
  
  out <- raster::getValues(rDat)
  
  if (output == "map") {
    
    out <- sp::plot(rDat)
    
  } else if (output == "hist") {
    
    out <- ggplot2::ggplot(data=NULL, ggplot2::aes(out)) + 
             ggplot2::geom_histogram(bins = 50) +
      xlab("Repeat visits") +
      theme_linedraw()
    
  }
  
  return(out)
  
}
