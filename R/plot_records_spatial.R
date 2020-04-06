#' Visualize spatial variation in recorder activity
#'
#' Map the spatial distribution of records in a chosen number of time-periods
#' @param inPath String. Directory in which the extract_records outputs were saved.
#' @param taxon String. Taxonomic group for which to plot records. Must match one of the groups passed to extract_records().
#' @param res Numeric. Spatial resolution (decimal degrees) to plot records at. Defaults to 0.5, i.e. half a degree.
#' @param n_periods Numeric. How many periods to split the data in to. Defaults to two.
#' @export
#' @examples
#'
plot_records_spatial <- function (inPath, taxon, res=0.5, n_periods=2) {
  path <- inPath
  files <- list.files(inPath, full.names = TRUE)
  taxa <- gsub("_.*", "", list.files(inPath))
  X <- 1:length(files)
  loadFun <- function(fileNum) {
    load(files[fileNum])
    assign(paste0("dat", fileNum), dat)
  }
  x <- lapply(X = X, loadFun)
  names(x) <- taxa

  rast <- raster::raster(ncol=length(seq(-120,-30,res)),
                 nrow=length(seq(-60,30,res)),
                 xmn=-120,
                 xmx=-30,
                 ymn=-60,
                 ymx=30)

  y <- x[[which(names(x) == taxon)]]

  yrCategories <- round(seq(min(y$year), max(y$year), length.out = (n_periods+1)), 0)

  yrsInCategory <- function(ind) {
    z <- yrCategories[ind]:yrCategories[ind+1]
    return(z)
  }

  X=1:(length(yrCategories)-1)

  z <- lapply(X=X, FUN=yrsInCategory)
  z

  subSet <- function(yrs) {
    subDf <- y[which(y$year %in% yrs),c("lon","lat")]
  }

  subs <- lapply(X=z, FUN=subSet)

  rasts <- lapply(X=subs, FUN=raster::rasterize, y=rast, fun="count")

  map <- ggplot2::map_data("world", xlim =c(-120,-30), ylim=c(-60, 30))

  myCol <- rgb(255,255,255, max = 255, alpha = 125, names = "blue50")


  plotDat <- function(raster) {
           ggplot2::ggplot() +
             ggplot2::geom_polygon(data = map, ggplot2::aes(x=long, y = lat, group = group),
                                   colour="black",alpha=1,fill=myCol,
                                   inherit.aes=F) +
             ggspatial::layer_spatial(as.logical(rasts[[raster]]), interpolate=F) +
             ggplot2::theme_linedraw(base_size=20) +
             ggplot2::scale_x_continuous(breaks = c(-120, -30, -60, 30)) +
             ggplot2::scale_fill_gradient2(low = "red", high = "blue", na.value = myCol,
                                name = "no. records",
                                limits=c(0, 1)) +
             ggplot2::theme(legend.position = "none",
                   legend.background = element_blank()) +
             ggplot2::xlab("") +
             ggplot2::ylab("") +
             ggplot2::ggtitle(paste(min(z[[raster]]), "-", max(z[[raster]]))) +
             ggplot2::theme(plot.margin=unit(c(0.2,0.2,0.2,0.2), "cm")) +
             ggplot2::coord_sf(xlim = c(-120, -30), ylim=c(-60,30))

  }

  p <- lapply(X= 1:length(rasts), FUN=plotDat)

return(p)

}

