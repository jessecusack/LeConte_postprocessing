# Load libraries
library(oce)
library(ncdf4)
source("ncwrite.R")

# Parameters
file <- "~/Dropbox/LeConte/Data/ocean/september2018/raw/moorings/ABLE_Sentinel/ADCP/LeConte S-ABLE Sept2018 20180901T233454.pd0"
pmin <- 125.  # cut off pressure [dbar]
dmax <- 150. # cut off distance [m]
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

# Subset data using distance
adp <- subset(adp, distance < dmax)

# Write beam output.
# adp_write(adp, "../proc/ABLE_sentinel_2018_beam.nc")

# Convert to xyz coordinates
xyz <- beamToXyz(adp)

xyzv <- xyz
# Extract vertical beam to new object
xyzv[["v"]][,,3] <- xyz[["vv"]]

# Convert to Earth coordinates
# The matlab file gives -18.7 for declination. 
# This website https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#declination
# gives a value of 19.3 depending on which standard you choose...
enu <- xyzToEnu(xyz, declination = dec)
# Apply rotation with 5th beam
enuv <- xyzToEnu(xyzv, declination = dec)
# Put rotated 5th beam back into original object
enu[["vv"]] <- enuv[["v"]][,,3]

enu_write(enu, "../proc/ABLE_sentinel_2018_enu.nc")
