#' extract records from GBIF and format them for use in soaR
#'
#' This function is used to scrape GBIF for occurrence records for a chosen taxonomic group, on a chosen continent.
#' @param continent String
#' @param year numeric
#' @param taxon string
#' @param write_output logical
#' @param output_path string
#' @export
#' @examples

extract_records <- function(continent, year, taxon, write_output, output_path) {


x <- rgbif::occ_data(hasCoordinate = T,
              hasGeospatialIssue = F,
              continent = continent,
              scientificName = taxon,
              year = year,
              limit = 200000)

## take the second element which is the "data" as opposed to the metadata

dat <- data.frame(x[2])

## check whether there are any records

if ("data.species" %in% names(dat)) {

## if there are create a data frame with variables of interest

  dat <- data.frame(dat$data.species, dat$data.scientificName, dat$data.decimalLongitude,
                  dat$data.decimalLatitude, dat$data.year)

  colnames(dat) <- c("species","group","lon","lat","year")

  dat <- dat[-which(is.na(dat$species)),]

  dat <- dat[-which(duplicated(dat)), ]

  if (length(dat[,1]) == 200000) {
    print("Warning: Reached max number of outputs, but more data is available.")
  }

 nSpec <- length(unique(dat$species))
 nRec <- length(dat[,1])

 attr(dat, "nSpec") <- as.numeric(nSpec)
 attr(dat, "nRec") <- as.numeric(nRec)

} else {
  dat <- 0
  print("this query produced zero records")
}

if (write_output == TRUE) {

  save(dat, file=paste0(output_path,
                   taxon, "_raw_data.rdata"))

}

return (dat)

}

?extract_records
