#' create a roster of arguments to be passed to extract_data.
#'
#' This function is creates a list of 1-row dataframes with the arguments needed for extract_records. This allows
#' the use of lapply on extract_records where multiple arguments are varied at the same time (e.g. taxa and country).
#'
#' @param index Numeric. Element number (i.e. 1: number of calls to extract_records).
#' @param continent String or character vector. E.g. "south_america", "europe". Do not specify unless country is NULL.
#' @param country string or character vector. ISO codes, e.g. "BR" for Brasil, "AR" for Argentina, and "CL" for Chile.
#'                See https://www.iban.com/country-codes.
#' @param year string or character vector. year(s) for which to extract data. For a range, format as e.g. "1995,2010". For a single year, as e.g. "2005".
#' @param taxa numeric or numeric vector. GBIF taxonomic key(s).
#' @param write logical. If true writes outputs to .csv in outPath
#' @param outPath string or character vector. Directory in which to save outputs.
#' @param outName string or character vector. Name of file to be saved if write_output = TRUE.
#' @param degrade logical. Whether or not to remove duplicates from the data (which may be repeat visits).
#' @export
#' @examples
#' roster <- createRoster(index = 1:12,
#' taxa = 811,
#' country = c("CL", "BR", "AR", "BO", "CO", "EC", "GY", "PY", "PE", "SR", "UY", "VE"),
#' continent = NA,
#' write = TRUE,
#' outName = "Diptera",
#' outPath = "C:/Users/Rob.Lenovo-PC/Documents/surpass/Data/GBIF/20.01.21/",
#' degrade = FALSE,
#' year = "1950, 2019")


createRoster <- function(index,
                         taxa,
                         country,
                         continent,
                         write,
                         identifier,
                         outPath,
                         degrade,
                         year) {

  df <- data.frame(index = index,
                   taxa = taxa,
                   country = country,
                   continent = continent,
                   write = write,
                   identifier = identifier,
                   outPath = outPath,
                   degrade = degrade,
                   year = year)

  roster <- split(df, seq(nrow(df)))

}
