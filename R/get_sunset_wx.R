get_sunset_wx <- function(wx_stns, start, ss_dt, lat, lon) {
  start <- as.Date(start)
  end <- start + as.difftime(2, units = "days")
  stn_ids <- strsplit(wx_stns, ", ")[[1]]

  # template for output in case no matching weather
  out <- data.frame()
  max_i <- length(stn_ids)
  i <- 1
  while (nrow(out) == 0) {
    if (i > max_i) {
      out <- tibble::tibble(ss_temp = as.numeric(NA),
                            ss_wsp = as.numeric(NA),
                            ss_wx_t_diff = as.integer(NA),
                            wx_stn = NA_character_,
                            wx_stn_km = as.integer(NA))
      break
    }
    out <- suppressWarnings(
      riem::riem_measures(stn_ids[i], start, end)
      # get_wx(stn_ids[i], start, end)
    )
    if (nrow(out) == 0) {
      out
    } else {
      out <- out %>%
        dplyr::select(valid, tmpf, sknt, wx_stn = station, wx_lat = lat, wx_lon = lon) %>%
        dplyr::mutate(diff = abs(as.numeric(difftime(valid, ss_dt, units = "mins")))) %>%
        dplyr::filter(diff <= 30) %>%
        na.omit()
      if (nrow(out) == 0) {
        out
      } else {
        out <- out %>%
          dplyr::slice(which.min(diff)) %>%
          dplyr::mutate(ss_temp = round((tmpf - 32) * 5/9, 1),
                        ss_wsp = round(sknt * 0.514444, 1),
                        ss_wx_t_diff = as.integer(round(diff)),
                        wx_stn_km = as.integer(
                          round(geosphere::distVincentyEllipsoid(c(lon,lat),
                                                                 c(wx_lon, wx_lat))/1000))) %>%
          dplyr::select(ss_temp, ss_wsp, ss_wx_t_diff, wx_stn, wx_stn_km)
      }
    }
    i <- i + 1  }
  return(out)
}