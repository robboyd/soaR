
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

    presProb <- raster::extract(out, dat$Presence)
    
    abProb <- raster::extract(out, dat$pseudoAbsence)

    obs <- c(rep(1, nrow(dat$Presence)), rep(0, nrow(dat$pseudoAbsence)))

    preds <- c(presProb, abProb)
    
    DATA <- data.frame(id = 1:length(obs),
                       obs = obs,
                       pred = preds)

    thresh <- seq(0,1, length.out = 50)
    
    opt <- PresenceAbsence::optimal.thresholds(DATA, thresh, opt.methods = 5, na.rm = T)

    if (map == TRUE) {
     
      binPred <- out
      
      binPred[binPred < as.numeric(opt[2])] <- 0
      
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
    
    opt = NULL
    
  }
  
  return(data.frame(opt = opt, species = species, group = group))  
  
}

