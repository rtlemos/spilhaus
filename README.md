# spilhaus

Project your maritime data onto the World Ocean.

![SST on Sep 16, 2019](https://github.com/rtlemos/spilhaus/blob/main/spilhaus_sst.png)

To know more about the Spilhaus projection, see:

- [Ocean Sciences with the Spilhaus Projection: A Seamless Ocean Map for 
  Spatial Data Recognition](https://www.nature.com/articles/s41597-023-02309-6)


Currently, the Spilhaus projection is missing in PROJ, making it 
unavailable for R and Python enthusiasts, among others: https://github.com/OSGeo/PROJ/issues/1851.
This repo contains four scripts:

- `spilhaus.R` and `spilhaus.py`, which contain the main functions that you may want to use, and
- `example.R` and `example.py`, which demonstrate how to use the functions to create and plot a dataset using the Spilhaus projection. 

In `spilhaus.*`, the functions

- `from_lonlat_to_spilhaus_xy` and 
- `from_spilhaus_xy_to_lonlat`

let you move between Mercator and Spilhaus projections; the function

- `make_spilhaus_xy_gridpoints`

lets you create a square grid of Spilhaus coordinates; and the function

- `pretify_spilhaus_df`

cleans the Spilhaus dataframe and moves a few pixels around (bottom/top and 
left/right) to make a contiguous ocean.
