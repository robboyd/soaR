#' Break down the data by country or basis of record
#'
#' This function splits the data obtained using extract_data by country or basis of record (e.g. fossil specimen, preserved specimen)
#'
#' @param inPath String. Location of extract_records outputs.
#' @param taxon String. Must match one of the taxon names in inPath.
#' @param start Numeric. Lower bound for analysis (year).
#' @param end Numeric. Upper bound for analysis (year).
#' @param output String. One of "country" or "basis".
#' @export
#'
recordBreakdown <- function(inPath, taxon, start, end, output) {

  load(paste0(inPath, taxon, "_raw_data.rdata"))

  dat <- dat[which(dat$year >= start & dat$year <= end), ]

  if (output == "country") {

    out <- table(dat$country)

  } else if (output == "basis") {

    out <- table(dat$basisOfRecord)
  }

  return(data.frame(out))

}

recordBreakdown(inPath = "F:/SURPASS/data/raw_data/",
                taxon = "Apidae",
                start = 1950,
                end = 2018,
                output = "country")




