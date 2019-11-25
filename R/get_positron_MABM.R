get_positron_MABM <- function(bbox, x = 14:18, y = 23:27,
                              zoom = 6, style = "light_nolabels") {
  xy <- expand.grid(x = x, y = y)
  base_url <- "https://a.basemaps.cartocdn.com"
  listOfTiles <- pbapply::pblapply(1:nrow(xy), function(i) {
    tmpx <- xy[i, "x"]
    tmpy <- xy[i, "y"]
    url <- paste0(file.path(base_url, style, zoom, tmpx, tmpy), ".png")
    message(url)
    response <- httr::GET(url)
    stopifnot(response$status_code == 200L)
    tile <- httr::content(response)
    tile <- aperm(tile, c(2, 1, 3))
    tile <- apply(tile, 2, rgb)
    lonlat_upperleft <- ggmap:::XY2LonLat(tmpx, tmpy, zoom)
    lonlat_lowerright <- ggmap:::XY2LonLat(tmpx, tmpy, zoom, 255L, 255L)
    bbox <- c(left = lonlat_upperleft$lon, bottom = lonlat_lowerright$lat, 
              right = lonlat_lowerright$lon, top = lonlat_upperleft$lat)
    bb <- tibble(ll.lat = unname(bbox["bottom"]), ll.lon = unname(bbox["left"]), 
                 ur.lat = unname(bbox["top"]), ur.lon = unname(bbox["right"]))
    class(tile) <- "raster"
    attr(tile, "bb") <- bb
    tile
  })
  
  ### ggmap:::stitch modified
  tiles <- listOfTiles
  ll.lat <- NULL
  rm(ll.lat)
  ll.lon <- NULL
  rm(ll.lon)
  bbs <- plyr::ldply(tiles, function(x) attr(x, "bb"))
  bigbb <- data.frame(ll.lat = min(bbs$ll.lat), ll.lon = min(bbs$ll.lon), 
                      ur.lat = max(bbs$ur.lat), ur.lon = max(bbs$ur.lon))
  tiles <- lapply(tiles, as.matrix)
  nrows <- length(unique(bbs$ll.lat))
  ncols <- length(unique(bbs$ll.lon))
  tiles <- split(tiles, rep(1:nrows, each = ncols))
  tiles <- lapply(tiles, function(x) Reduce(cbind, x))
  tiles <- Reduce(rbind, tiles)
  tiles <- as.raster(tiles)
  class(tiles) <- c("ggmap", "raster")
  attr(tiles, "bb") <- bigbb
  map <- tiles
  ### ggmap:::stitch
  
  mbbox <- attr(map, "bb")
  size <- 256L * c(length(x), length(y))
  slon <- seq(mbbox$ll.lon, mbbox$ur.lon, length.out = size[1])
  slat <- vector("double", length = 256L * length(y))
  for (k in seq_along(y)) {
    slat[(k - 1) * 256 + 1:256] <- sapply(as.list(0:255), 
                                          function(y_i) {
                                            XY2LonLat(X = x[1], Y = y[k], 
                                                      zoom, x = 0, y = y_i)$lat
                                          })
  }
  slat <- rev(slat)
  keep_x_ndcs <- which(bbox["left"] <= slon & slon <= 
                         bbox["right"])
  keep_y_ndcs <- sort(size[2] - which(bbox["bottom"] <= 
                                        slat & slat <= bbox["top"]))
  croppedmap <- map[keep_y_ndcs, keep_x_ndcs]
  croppedmap <- as.raster(croppedmap)
  class(croppedmap) <- c("ggmap", "raster")
  attr(croppedmap, "bb") <- data.frame(ll.lat = bbox["bottom"], 
                                       ll.lon = bbox["left"], ur.lat = bbox["top"], 
                                       ur.lon = bbox["right"])
  attr(croppedmap, "source") <- "Carto"
  attr(croppedmap, "maptype") <- style
  attr(croppedmap, "zoom") <- zoom
  return(croppedmap)
}
