#' extract records from GBIF and format them for use in soaR
#'
#' This function is used to scrape GBIF for occurrence records for a chosen taxonomic group, on a chosen continent.
#' @param continent String. E.g. "south_america", "europe". Do not specify unless country is NULL.
#' @param country string. E.g. "BR" for Brasil, "AR" for Argentina, and "CL" for Chile
#' @param year string. year(s) for which to extract data. For a range, format as e.g. "1995,2010". For a single year, as e.g. "2005".
#' @param taxon numeric. GBIF taxonomic key.
#' @param write_output logical
#' @param output_path string. Directory in which to save outputs.
#' @param output_name string. Name of file to be saved if write_output = TRUE.
#' @export
#' @examples

extract_records <- function(continent = NULL, country = NULL, year, taxon, write_output, output_path, output_name) {

  if (!is.null(continent) & !is.null(country)) {

    stop("Can't specify both country and continent")

  }

if (!is.null(continent)) {

  x <- rgbif::occ_data(hasCoordinate = T,
                hasGeospatialIssue = F,
                continent = continent,
                taxonKey = taxon,
                year = year,
                limit = 200000)

} else if (!is.null(country)) {

  x <- rgbif::occ_data(hasCoordinate = T,
                       hasGeospatialIssue = F,
                       country = country,
                       taxonKey = taxon,
                       year = year,
                       limit = 200000)
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
                  dat$data.basisOfRecord, dat$data.catalogNumber)

  colnames(dat) <- c("species","group","lon","lat","year", "Date", "originalDate", "country","continent","basisOfRecord",
                     "catalogNum")

  if (length(dat[,1]) == 200000) {
    warning("Reached max number of outputs, but more data is available.")
  }

  dat <- dat[-which(is.na(dat$species)),]

  dat <- dat[-which(duplicated(dat)), ]

 nSpec <- length(unique(dat$species))
 nRec <- length(dat[,1])

 attr(dat, "nSpec") <- as.numeric(nSpec)
 attr(dat, "nRec") <- as.numeric(nRec)

} else {
  dat <- 0
  warning("This query produced zero records")
}

if (write_output == TRUE) {

  save(dat, file=paste0(output_path, output_name,
                   "_raw_data.rdata"))

}

return (dat)

}


