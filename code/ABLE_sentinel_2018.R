library(oce)
library(ncdf4)
source("ncwrite.R")

# Parameters
file <- "~/Dropbox/LeConte/Data/ocean/september2018/raw/moorings/ABLE_Sentinel/ADCP/LeConte S-ABLE Sept2018 20180901T233454.pd0"
# pmin <- 120.  # cut off pressure [dbar]
lat <- 56.8370392
lon <- -132.3574123
dec <- 19.32  # magnetic declination
# n <- 6  # ensemble averaging
ori <- "upward"  # orientation


# Load data and remove times where the instrument was not in the water
file <- path.expand(file)
adp <- read.adp(file, latitude = lat, longitude = lon)
adp <- oceSetMetadata(adp, 'orientation', ori)
# adp <- subset(adp, pressure > pmin)

adp_write(adp, "../proc/ABLE_sentinel_2018.nc")

# Ensemble average the data, 
# adp <- adpEnsembleAverage(adp, n = n)
# Convert to xyz coordinates
xyz <- beamToXyz(adp)

# Convert to Earth coordinates
# The matlab file gives -18.7 for declination. 
# This website https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#declination
# gives a value of 19.3 depending on which standard you choose...
enu <- xyzToEnu(xyz, declination = dec)

enu_write(enu, "../proc/ABLE_sentinel_2018_enu.nc")
