# In this example, we want to plot a set of lat-lon sites, instead of a raster
# See output in spilhaus_points.png

library(ggplot2)
source('~/Downloads/spilhaus/spilhaus.R')

# auxiliary function 1: fetch 5km resolution land mask
download_land_mask = function() {
  myInfo = rerddap::info('NOAA_DHW_monthly', url='https://coastwatch.pfeg.noaa.gov/erddap/')
  myData = rerddap::griddap(datasetx = myInfo,
                            fields = "mask",
                            latitude = c(-89.975, 89.975),
                            longitude = c(-179.975, 179.975),
                            time = c('2024-01-01', '2024-01-01'))
  da = ncdf4::nc_open(myData$summary$filename)
  return(da)
}

# auxiliary function 2: determine whether arbitrary lon-lat points fall on land
extract_mask = function(da, lonlat) {
  lons = sort(ncdf4::ncvar_get(da, "longitude"))
  lats = sort(ncdf4::ncvar_get(da, "latitude"))
  dlon = lons[2] - lons[1]
  dlat = lats[2] - lats[1]
  ln = as.integer(round((lonlat[,1] - lons[1]) / dlon) + 1)
  la = as.integer(round((lonlat[,2] - lats[1]) / dlat) + 1)
  msk = ncdf4::ncvar_get(da, "mask", start = c(1,1,1), count = c(7200, 3600,1))
  msk = msk[,ncol(msk):1]
  chunk = msk[ln + (la - 1) * dim(msk)[1]]
  chunk[chunk == 1] = NA # in the original mask, 1 denotes land
  chunk[chunk == 4] = 0 # in the original mask, 4 denotes ice
  return(chunk)
}

# auxiliary function 3: determine land-water interface
extract_interface = function(da, lonlat) {
  lons = sort(ncdf4::ncvar_get(da, "longitude"))
  lats = sort(ncdf4::ncvar_get(da, "latitude"))
  dlon = lons[2] - lons[1]
  dlat = lats[2] - lats[1]
  ln = as.integer(round((lonlat[,1] - lons[1]) / dlon) + 1)
  la = as.integer(round((lonlat[,2] - lats[1]) / dlat) + 1)
  msk = ncdf4::ncvar_get(da, "mask", start = c(1,1,1), count = c(7200, 3600,1))
  msk[msk == 1] = NA # in the original mask, 1 denotes land
  msk[msk == 4] = 0 # in the original mask, 4 denotes ice
  msk = msk[,ncol(msk):1]
  n = round(sqrt(nrow(lonlat)))
  chunk = matrix(nrow=n, ncol=n, msk[ln + (la - 1) * dim(msk)[1]])
  chunk = mapply(1:n, FUN=function(j) mapply(1:n, FUN=function(i) {
    if ((is.na(chunk[i, j]) & any(!is.na(chunk[max(1, i-1):min(n, i+1), max(1, j-1):min(n, j+1)])))) {
      TRUE
    } else {
      NA
    }
  }))
  return(as.numeric(chunk))
}

# step 1: create a table of x, y points on a Splihaus map
spilhaus_df = make_spilhaus_xy_gridpoints(spilhaus_res=1000)
x_lo = min(spilhaus_df$x); x_hi = max(spilhaus_df$x)
y_lo = min(spilhaus_df$y); y_hi = max(spilhaus_df$y)

# step 2: convert those points to lon and lat
lonlat = from_spilhaus_xy_to_lonlat(spilhaus_df$x, spilhaus_df$y)

# step 3: determine which points fall on land, to create mask
da = download_land_mask()
spilhaus_df$z = extract_mask(da, lonlat)
spilhaus_df$l = is.na(spilhaus_df$z)

# step 4: determine land-water interface
interface_df = make_spilhaus_xy_gridpoints(spilhaus_res=1000)
interface_df$z = extract_interface(da, lonlat)
interface_df$l = is.na(interface_df$z)

# step 5: load table with lon-lat coordinates of observations
sites = read.csv('https://coastwatch.pfeg.noaa.gov/erddap/tabledap/erdGodaeSfcobs.csv?longitude%2Clatitude%2Ctime%2Cob_sst%2Cob_wm&time%3E=1999-05-01&time%3C=1999-05-01', skip = 1)

# step 6: convert lon-lats of obs into Spilhaus x-y, and grid them manually
sites_xy = from_lonlat_to_spilhaus_xy(sites$degrees_east, sites$degrees_north)
n_grid = 100  # number of gridcells in Spilhaus x and y coordinates
gridded_x_pos = ceiling(n_grid * (sites_xy[,1] - x_lo) / (x_hi - x_lo))
gridded_y_pos = ceiling(n_grid * (sites_xy[,2] - y_lo) / (y_hi - y_lo))
sites_grid = expand.grid(x=seq(x_lo, x_hi, len=n_grid), y=seq(y_lo, y_hi, len=n_grid))
sites_grid$z = 0
for (i in seq(1, nrow(sites_xy))) {
  pos = gridded_x_pos[i] + (gridded_y_pos[i] - 1) * n_grid
  sites_grid$z[pos] = sites_grid$z[pos] + 1
}
sites_grid$l = (sites_grid$z == 0)

# step 7: (create a single world ocean)
pretty_spilhaus_df = pretify_spilhaus_df(spilhaus_df)
pretty_interface_df = pretify_spilhaus_df(interface_df)
pretty_sites = pretify_spilhaus_df(sites_grid)

# step 8: plot
ggplot() +
  geom_tile(data=pretty_spilhaus_df, aes(x=x, y=y), fill='white') +
  geom_point(data=pretty_interface_df, aes(x=x, y=y), color=gray(0.4), size=0.1, pch='.') +
  geom_tile(data=pretty_sites, aes(x=x, y=y, fill=z)) +
  scale_fill_viridis_c(option = 'viridis', direction = -1, name = "# samples") +
  coord_equal() +
  theme_void() +
  theme(panel.background = element_rect(fill = 'lightgray', color = 'black'),
        legend.position = c(0.95, 0.9), legend.direction = "vertical",
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black"))



