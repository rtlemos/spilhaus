import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from spilhaus import *

# loading remote dataset
ds = xr.open_dataset('https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/temperature/decav/0.25/woa18_decav_t00_04.nc',
                     decode_times=False)
# extracting SST data array
da = ds["t_mn"].sel(depth=0, time=4.326e+03)

# creating data frame with desired resolution
spilhaus_df = make_spilhaus_xy_gridpoints(spilhaus_res=2000)
lon, lat = from_spilhaus_xy_to_lonlat(spilhaus_df['x'], spilhaus_df['y'])

# extracting SST data from dataset
spilhaus_df['z'] = da.sel(lon=xr.DataArray(lon, dims="points"),
                          lat=xr.DataArray(lat, dims="points"),
                          method="nearest").data

# prettifying
pretty_spilhaus_df = prettify_spilhaus_df(spilhaus_df)

# plotting
background_color = "black"
foreground_color = "white"
fig, ax = plt.subplots(1, 1, figsize=(16,16), dpi=300)
plt.scatter(x=pretty_spilhaus_df["x"], y=pretty_spilhaus_df["y"], c=pretty_spilhaus_df["z"],
            marker='s', s=72./fig.dpi, cmap="inferno")
cbar_ax = fig.add_axes([0.85, 0.7, 0.025, 0.15])
cbar = plt.colorbar(cax=cbar_ax, ticks=np.linspace(0, 30, 7))
cbar.set_label('SST [°C]', labelpad=-40, y=1.1, rotation=0, color=foreground_color)
cbar.ax.yaxis.set_tick_params(color=foreground_color)
cbar.outline.set_edgecolor(foreground_color)
plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color=foreground_color)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])

if background_color is not None:
    ax.patch.set_facecolor(background_color)
plt.savefig('./annual_sst_spilhaus.png', transparent=True)
plt.show()
