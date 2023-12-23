# EXAMPLE: 
# Extract sea surface temperature data from a remote repo and plot it using the Spilhaus projection

library("rerddap")
library("ggplot2")
library("ncdf4")

source("./spilhaus.R") # load the Spilhaus projection script

# auxiliary function to extract netCDF file from Coastwatch
download_sst_data = function(start_date, end_date) {
  myInfo = rerddap::info('NOAA_DHW_monthly', url='https://coastwatch.pfeg.noaa.gov/erddap/')
  myData = rerddap::griddap(datasetx = myInfo,
                            fields = "sea_surface_temperature",
                            latitude = c(-89.975, 89.975),
                            longitude = c(-179.975, 179.975),
                            time = c(start_date, end_date))
  da = ncdf4::nc_open(myData$summary$filename)
  return(da)
}

# auxiliary function to extract a list (length = number of instants) of snapshots (7200x3600 pixels) from the netCDF file
extract_sst_data = function(da, lonlat) {
  lons = sort(ncdf4::ncvar_get(da, "longitude"))
  lats = sort(ncdf4::ncvar_get(da, "latitude"))
  times = ncdf4::ncvar_get(da, "time")

  dlon = lons[2] - lons[1]
  dlat = lats[2] - lats[1]
  ln = as.integer(round((lonlat[,1] - lons[1]) / dlon) + 1)
  la = as.integer(round((lonlat[,2] - lats[1]) / dlat) + 1)
  get_chunk = function(timepoint) {
    sst = ncdf4::ncvar_get(da, "sea_surface_temperature",
                           start = c(1,1,timepoint), count = c(7200, 3600, 1))
    sst = sst[,ncol(sst):1]
    chunk = sst[ln + (la - 1) * dim(sst)[1]]
    return(chunk)
  }
  sst_data = mapply(1:length(times), FUN=function(timepoint) {get_chunk(timepoint)})
  return(sst_data)
}
                    
# first we crate a data frame with NxN pixels
spilhaus_df = make_spilhaus_xy_gridpoints(spilhaus_res=1000)
# convert the spilhaus coordinates into Mercator coordinates
lonlat = from_spilhaus_xy_to_lonlat(spilhaus_df$x, spilhaus_df$y)

# download the netCDF data from Coastwatch
da = download_sst_data("2021-09-16", "2021-09-16")

# extract the SST data for the required lats and lons
spilhaus_df$z = extract_sst_data(da, lonlat)
# mask
spilhaus_df$l = is.na(spilhaus_df$z)
# prettify
pretty_spilhaus_df = pretify_spilhaus_df(spilhaus_df)
# plot
ggplot(data=pretty_spilhaus_df, aes(x=x, y=y, fill=z)) +
  geom_raster() +
  scale_fill_viridis_c(option = 'inferno', guide = "none") +
  coord_equal() +
  theme_void() +
  theme(panel.background = element_rect(fill = 'black', color = 'black'))


