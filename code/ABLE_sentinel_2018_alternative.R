# Load libraries

library(oce)
library(ncdf4)
library(yaml)
source("ncwrite.R")

# Set parameters

# Parameters
file <- "~/Dropbox/LeConte/Data/ocean/september2018/raw/moorings/ABLE_Sentinel/ADCP/LeConte S-ABLE Sept2018 20180901T233454.pd0"
pmin <- 125.  # cut off pressure [dbar]
lat <- 56.8370392
lon <- -132.3574123
dec <- 19.32  # magnetic declination
ori <- "upward"  # orientation


# Load data and set orientation

file <- path.expand(file)
adp <- read.adp(file, latitude = lat, longitude = lon)
adp <- oceSetMetadata(adp, 'orientation', ori)

# Subset data using pressure

adp <- subset(adp, pressure > pmin)

# +
# adp_write(adp, "../proc/ABLE_sentinel_2018.nc")
# -

# Ensemble average the data, 
# adp <- adpEnsembleAverage(adp, n = n)
# Convert to xyz coordinates
xyz <- beamToXyz(adp)

xyzv <- xyz
xyzv[["v"]][,,3] <- xyz[["vv"]]


# Convert to Earth coordinates
# The matlab file gives -18.7 for declination. 
# This website https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#declination
# gives a value of 19.3 depending on which standard you choose...
enu <- xyzToEnu(xyz, declination = dec)
enuv <- xyzToEnu(xyzv, declination = dec)
enu[["vv"]] <- enuv[["v"]][,,3]

enu_write(enu, "../proc/ABLE_sentinel_2018_enu.nc")


