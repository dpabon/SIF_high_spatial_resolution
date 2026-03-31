require(terra)
require(ncdf4)
require(tidyr)
require(dplyr)



# general info to set
comment_attr <- "This dataset was produced by the Max Planck Institute for Biogeochemistry within the ESA Sen4GPP project"
target_file_prefix <- "sif_downscaled_"

# ifile <- "../../scratch/sif/Sen3_OGVI______Sen3_LST_day/2018-05-11_Sen3_OGVI__Sen3_LST_day_1km_10km.nc"
# ofile <-  "../../scratch/sif_downscaled_curated/.nc"

# setup paths
ipath <- "../../scratch/sif/Sen3_OGVI______Sen3_LST_day"
opath <-  "../../scratch/sif/sif_downscaled_curated"
dir.create(opath, showWarnings = F)
quicklookpath <- "../../scratch/sif/sif_downscaled_curated_qck"
dir.create(quicklookpath, showWarnings = F)

# get files
file_list <- list.files(ipath, pattern = "Sen3_OGVI__Sen3_LST_day_1km_10km.nc")

for(ifile in file_list){ 
  
  # get time stamp
  file_time <- substring(ifile, 1, 10)
  
  # set out file name
  ofile <- paste0(opath, "/", target_file_prefix, file_time, '.nc')
  
  nc_in <- nc_open(paste0(ipath, "/", ifile))
  
  # get data and flip it as needed
  dat <- as.matrix(ncvar_get(nc = nc_in, varid = "sif_pred")) 
  dat <- apply(dat, 1, rev)
  dat <- apply(dat, 2, rev)
  
  # remove strange Infinite data and set all to NA
  dat[is.infinite(dat)] <- NA
  dat[is.nan(dat)] <- NA
  
  # remote widely implausible values of SIF
  dat[dat < -2] <- NA
  dat[dat > 10] <- NA
  
  # set the data in a terra raster
  r <- rast(dat)
  
  # set projection (sinusoidal used by MODIS)
  crs(r) <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
  
  # set extent based on known tile limits of MODIS files
  ext(r) <- ext(-1111951, 3335852, 3335852, 7783654)
  
  # make a quicklook to check all is ok
  png(filename = paste0(quicklookpath, "/qck_", target_file_prefix, file_time, ".png"))
  plot(r, range = c(-0.25, 4.25))
  title(paste("Downscaled SIF at time:", file_time))
  dev.off()
  
  # write file as NETCDF via terra package
  writeCDF(r, ofile, varname = "sif", overwrite = T, zname = "time",
           longname = "TROPOMI sun-induced chlorophyll fluorescence (SIF) downscaled to 1km based on Sentinel-3 OGVI and LST_Day data",
           compression = 5,
           unit = "mW/m2/sr/nm")
  
  # check and add attribute
  nc_new <- nc_open(ofile, write = T)
  ncatt_put(nc_new, attname = "comment", attval = comment_attr, varid = 0)
  nc_close(nc_new)
  
  print(ofile)
}