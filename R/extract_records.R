#' extract records from GBIF and format them for use in soaR
#'
#' This function is used to scrape GBIF for occurrence records for a chosen taxonomic group, on a chosen continent.
#' @param continent String. E.g. "south_america", "europe". Do not specify unless country is NULL.
#' @param country string. E.g. "BR" for Brasil, "AR" for Argentina, and "CL" for Chile
#' @param year string. year(s) for which to extract data. For a range, format as e.g. "1995,2010". For a single year, as e.g. "2005".
#' @param taxon numeric. GBIF taxonomic key.
#' @param write logical
#' @param outPath string. Directory in which to save outputs.
#' @param outName string. Name of file to be saved if write_output = TRUE.
#' @param degrade logical. Whether or not to remove duplicates from the data (which may be repeat visits).
#' @export
#' @examples

extract_records <- function(roster) {

  if (!is.na(roster$continent) & !is.na(roster$country)) {

    stop("Can't specify both country and continent")

  }

if (!is.na(roster$continent)) {

  x <- rgbif::occ_data(hasCoordinate = T,
                hasGeospatialIssue = F,
                continent = as.character(roster$continent),
                taxonKey = roster$taxa,
                year = roster$year,
                limit = 199999)

} else if (!is.na(roster$country)) {

  x <- rgbif::occ_data(hasCoordinate = T,
                       hasGeospatialIssue = F,
                       country = as.character(roster$country),
                       taxonKey = roster$taxa,
                       year = roster$year,
                       limit = 199999)
} else {

  stop("Must specify one of country or continent")

}

## take the second element which is the "data" as opposed to the metadata

dat <- data.frame(x[2])

## check whether there are any records

if ("data.species" %in% names(dat)) {

## if there are create a data frame with variables of interest

  dat <- data.frame(dat$data.species, dat$data.scientificName, dat$data.decimalLongitude,
                  dat$data.decimalLatitude, dat$data.year, dat$data.eventDate, dat$data.verbatimEventDate,
                  dat$data.country, dat$data.continent,
                  dat$data.basisOfRecord, dat$data.datasetKey,
                  dat$data.coordinateUncertaintyInMeters,
                  dat$data.bibliographicCitation,
                  dat$data.country)

  colnames(dat) <- c("species","group","lon","lat","year", "Date", "originalDate", "country","continent","basisOfRecord",
                     "UUID", "spatialUncertainty", "ref", "country")

  if (length(dat[,1]) == 199000) {
    warning("Reached max number of outputs, but more data is available.")
  }

  ## remove recods not identified to species level

  dat <- dat[!is.na(dat$species),]

  if (roster$degrade == TRUE) {

    if (any(duplicated(dat))) {
      dat <- dat[-which(duplicated(dat)), ]
    }

  }


 nSpec <- length(unique(dat$species))
 nRec <- length(dat[,1])

 attr(dat, "nSpec") <- as.numeric(nSpec)
 attr(dat, "nRec") <- as.numeric(nRec)

} else {
  dat <- 0
  warning("This query produced zero records")
}

if (roster$write == TRUE) {

  write.csv(dat,
            paste0(roster$outPath, roster$outName,
                   "_", roster$country, ".csv"))

}

if (roster$write == FALSE) {

  return (dat)

}

print(paste("Extraction completed for", roster$outName, "in", roster$country))

}

