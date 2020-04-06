#' Summarise temporal variation in recorder activity
#'
#' Plot a time series of the number of records (1/year) for taxa selected using extract_records
#' @param inPath String. Directory in which the extract_records outputs were saved.
#' @export
#' @examples

plot_records_temporal <- function(inPath) {

  path <- inPath

  files <- list.files(inPath, full.names = TRUE)

  taxa <- gsub("_.*","",list.files(inPath))

  X <- 1:length(files)

  loadFun <- function(fileNum) {
    load(files[fileNum])
    assign(paste0("dat", fileNum), dat)
  }

  x <- lapply(X= X, loadFun)

  names(x) <- taxa

  x <- purrr::map_dfr(x, "year", .f = plyr::count, .id = "species")

  x <- purrr::map_df(x, rbind)

  p <- ggplot2::ggplot(data=x, ggplot2::aes(x=year, y= freq, fill=species)) +
    ggplot2::geom_bar(colour ="white", stat="identity") +
    ggplot2::theme_linedraw() +
    ggplot2::xlab("") +
    ggplot2::ylab("Count") +
    ggplot2::ggtitle("Number of records")

  return(p)

}

