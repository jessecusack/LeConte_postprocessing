library(oce)
library(ncdf4)
source("ncwrite.R")

# Parameters
file <- "~/Dropbox/LeConte/Data/ocean/september2018/raw/moorings/ABLE_Deep/16670013.000"
# pmin <- 120.  # cut off pressure [dbar]
lat <- 56.835592	
lon <- -132.3572915
dec <- 19.32  # magnetic declination
# n <- 6  # ensemble averaging
ori <- "upward"  # orientation

# Load data and remove times where the instrument was not in the water
file <- path.expand(file)
adp <- read.adp(file, latitude = lat, longitude = lon)
adp <- oceSetMetadata(adp, 'orientation', ori)
# adp <- subset(adp, pressure > pmin)

# Ensemble average the data, 
# adp <- adpEnsembleAverage(adp, n = n)
# Convert to xyz coordinates
xyz <- beamToXyz(adp)

# Convert to Earth coordinates
# The matlab file gives -18.7 for declination. 
# This website https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#declination
# gives a value of 19 or 19.3 depending on which standard you choose...
enu <- xyzToEnu(xyz, declination = dec)

adp_write(enu, "../proc/ABLE_deep_2018_enu.nc")
