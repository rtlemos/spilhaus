library(sf)

from_lonlat_to_spilhaus_xy <- function(longitude, latitude){

  # constants (https://github.com/OSGeo/PROJ/issues/1851)
  e = sqrt(0.00669438)
  lat_center_deg = -49.56371678
  lon_center_deg = 66.94970198
  azimuth_deg = 40.17823482

  # parameters derived from constants
  lat_center_rad = lat_center_deg * pi / 180
  lon_center_rad = lon_center_deg * pi / 180
  azimuth_rad = azimuth_deg * pi / 180
  conformal_lat_center = -pi / 2 + 2 * atan(
    tan(pi/4 + lat_center_rad/2) *
      ((1 - e * sin(lat_center_rad)) / (1 + e * sin(lat_center_rad))) ^ (e / 2)
  )
  alpha = -asin(cos(conformal_lat_center) * cos(azimuth_rad))
  lambda_0 = lon_center_rad + atan2(tan(azimuth_rad), -sin(conformal_lat_center))
  beta = pi + atan2(-sin(azimuth_rad), -tan(conformal_lat_center))

  # coordinates in radians
  lon = longitude * pi / 180
  lat = latitude * pi / 180

  # conformal latitude, in radians
  lat_c = -pi / 2 + 2 * atan(
    tan(pi/4 + lat/2) * ((1 - e * sin(lat)) / (1 + e * sin(lat))) ^ (e / 2)
  )

  # transformed lat and lon, in degrees
  lat_s = 180 / pi * asin(sin(alpha) * sin(lat_c) - cos(alpha) * cos(lat_c) * cos(lon - lambda_0))
  lon_s = 180 / pi * (
    beta + atan2(
      cos(lat_c) * sin(lon - lambda_0),
      (sin(alpha) * cos(lat_c) * cos(lon - lambda_0) + cos(alpha) * sin(lat_c))
    )
  )

  # projects transformed coordinates onto plane (Adams World in a Square II)
  adams_ws2 = "+proj=adams_ws2 +no_defs +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
  projected = sf_project(from=sf::st_crs(4326), to=adams_ws2, pts=cbind(lon_s, lat_s))
  adams_x = projected[,1]
  adams_y = projected[,2]
  spilhaus_x = -(adams_x + adams_y) / sqrt(2)
  spilhaus_y = (adams_x - adams_y) / sqrt(2)

  return(cbind(spilhaus_x, spilhaus_y)) #, adams_x, adams_y, lon_s, lat_s))
}

make_spilhaus_xy_gridpoints <- function(spilhaus_res=1000) {
  # regular grid of points in Spilhaus map
  extreme = 11825474
  m = seq(-extreme, extreme, len=spilhaus_res)
  spilhaus_df = expand.grid(x=m, y=m)
  return(spilhaus_df)
}

flip1 = function(x) {
  res = as.integer(sqrt(length(x)))
  return(c(t(matrix(nrow=res, ncol=res, x)[,res:1])))
}
flip2 = function(x) {
  res = as.integer(sqrt(length(x)))
  return(c(t(matrix(nrow=res, ncol=res, x))[,res:1]))
}



pretify_spilhaus_df <- function(spilhaus_df) {

    spilhaus_x = spilhaus_df$x
    spilhaus_y = spilhaus_df$y
    spilhaus_z = spilhaus_df$z
    spilhaus_l = spilhaus_df$l

    extreme = 11825474

    # augmented grid points
    aug_x = c(
      spilhaus_x,
      spilhaus_x - 2 * extreme,
      spilhaus_x + 2 * extreme,
      spilhaus_x,
      spilhaus_x
    )
    aug_y = c(
      spilhaus_y,
      spilhaus_y,
      spilhaus_y,
      spilhaus_y - 2 * extreme,
      spilhaus_y + 2 * extreme
    )
    aug_z = c(
      spilhaus_z,
      flip1(spilhaus_z),
      flip1(spilhaus_z),
      flip2(spilhaus_z),
      flip2(spilhaus_z)
    )
    aug_l = c(
      spilhaus_l,
      flip1(spilhaus_l),
      flip1(spilhaus_l),
      flip2(spilhaus_l),
      flip2(spilhaus_l)
    )

    cutpoint = 1.1 * extreme
    keep = !(
      aug_l
      | (aug_x < -cutpoint)
      | (aug_x > cutpoint)
      | (aug_y < -cutpoint)
      | (aug_y > cutpoint)
      | (aug_y > 1.089e7 - 0.176 * aug_x)
      | (aug_x < -0.984e7 - 0.565 * aug_y)
      | (aug_y < -1.378e7 + 0.46 * aug_x)
      | (aug_x > 1.274e7 + 0.172 * aug_y)
      | (aug_y > 1e7 - 0.5 * aug_x)
      | (aug_y > 2.3e7 + aug_x)
      | ((aug_y < 0.29e7) & (aug_x < -1.114e7))
      | ((aug_y < 0.39e7) & (aug_x < -1.17e7))
      | ((aug_y < -1.21e7) & (aug_x > 0.295e7))
      | ((aug_y < -1.2e7) & (aug_x > 0.312e7))
      | ((aug_y < -1.16e7) & (aug_x > 0.4e7))
      | ((aug_y < -1.11e7) & (aug_x > 0.45e7))
    )

    pretty_spilhaus_df = data.frame(
      x=aug_x[keep],
      y=aug_y[keep],
      z=aug_z[keep]
    )
    pretty_spilhaus_df = pretty_spilhaus_df[!duplicated(pretty_spilhaus_df[c(1,2)]),]

    return(pretty_spilhaus_df)
}

from_spilhaus_xy_to_lonlat <- function(spilhaus_x, spilhaus_y) {
  # constants
  e = sqrt(0.00669438)
  lat_center_deg = -49.56371678
  lon_center_deg = 66.94970198
  azimuth_deg = 40.17823482

  # parameters derived from constants
  lat_center_rad = lat_center_deg * pi / 180
  lon_center_rad = lon_center_deg * pi / 180
  azimuth_rad = azimuth_deg * pi / 180
  conformal_lat_center = -pi / 2 + 2 * atan(
    tan(pi/4 + lat_center_rad/2) *
      ((1 - e * sin(lat_center_rad)) / (1 + e * sin(lat_center_rad))) ^ (e / 2)
  )
  alpha = -asin(cos(conformal_lat_center) * cos(azimuth_rad))
  lambda_0 = lon_center_rad + atan2(tan(azimuth_rad), -sin(conformal_lat_center))
  beta = pi + atan2(-sin(azimuth_rad), -tan(conformal_lat_center))

  adams_x = (spilhaus_y - spilhaus_x) * sqrt(2) / 2
  adams_y = - (spilhaus_y + spilhaus_x) * sqrt(2) / 2

  adams_ws2 = "+proj=adams_ws2 +no_defs +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
  projection_fun = function(x, y) {
    tryCatch(
      sf_project(from=adams_ws2, to=sf::st_crs(4326), pts=c(x, y)),
      error=function(e) c(NA, NA)
    )
  }
  #projected = matrix(mapply(adams_x, adams_y, FUN=projection_fun), ncol=2)
  projected = sf_project(from=adams_ws2, to=sf::st_crs(4326), pts=cbind(adams_x, adams_y), keep = TRUE, warn = FALSE)
  lon_s = projected[,1]
  lat_s = projected[,2]

  #transformed coords in radians
  lon_s_rad = lon_s * pi / 180
  lat_s_rad = lat_s * pi / 180

  # conformal latitude
  lat_c = asin(sin(alpha) * sin(lat_s_rad) + cos(alpha) * cos(lat_s_rad) * cos(lon_s_rad - beta))

  # longitude, in radians
  lon = lambda_0 + atan2(
    cos(lat_s_rad) * sin(lon_s_rad - beta),
    sin(alpha) * cos(lat_s_rad) * cos(lon_s_rad - beta) - cos(alpha) * sin(lat_s_rad)
  )

  # latitude (iterative formula from https://mathworld.wolfram.com/ConformalLatitude.html)
  lat = lat_c
  for (i in 0:9) {
    lat = -0.5 * pi + 2 * atan(
      tan(pi / 4 + lat_c / 2) *
        ((1 + e * sin(lat)) / (1 - e * sin(lat))) ^ (e / 2)
    )
  }

  # coordinates in degrees
  longitude = ((lon * 180 / pi + 180) %% 360) - 180
  latitude = lat * 180 / pi

  return(cbind(longitude, latitude))
}
