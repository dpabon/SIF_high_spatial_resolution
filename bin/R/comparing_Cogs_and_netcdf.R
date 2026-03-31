dyn.load('/opt/ohpc/pub/libs/hwloc/lib/libhwloc.so.15')
dyn.load('/opt/ohpc/pub/libs/gnu12/hdf5/1.10.8/lib/libhdf5_hl.so.100')
dyn.load("/opt/ohpc/pub/libs/gnu12-openmpi4/gdal/3.5.3/lib64/libgdal.so.31")
dyn.load('/opt/ohpc/pub/libs/hwloc/lib/libhwloc.so.15')
dyn.load('/opt/ohpc/pub/libs/gnu9/openmpi4/hdf5/1.10.8/lib/libhdf5_hl.so.100')


library(terra)
library(tidyverse)
library(ncdf4)

library(tidync)


setwd("/Net/Groups/BGI/work_3/OEMC/oemc_sif/")


cog_files <- list.files(path = "data/OEMC_app_data/COG_SIF_files/", pattern = "*.tif", full.names = T)

test_rast <- rast(cog_files[1])

plot(test_rast)
