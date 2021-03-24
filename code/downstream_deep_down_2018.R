library(oce)
library(ncdf4)
source("create_nc_from_adp.R")

# Parameters
file <- "~/Dropbox/LeConte/Data/ocean/september2018/raw/moorings/Downstream_Deep/ADCP/SN_54000.000"
# pmin <- 70.  # cut off pressure [dbar]
lat <- 56.83596667
lon <- -132.3616833
dec <- 19.32  # magnetic declination
n <- 8  # ensemble averaging
ori <- "downward"  # orientation


# Load data and remove times where the instrument was not in the water
file <- path.expand(file)
adp <- read.adp(file, latitude = lat, longitude = lon)
adp <- oceSetMetadata(adp, 'orientation', ori)
# adp <- subset(adp, pressure > pmin)

# Ensemble average the data, 
adp <- adpEnsembleAverage(adp, n = n)
# Convert to xyz coordinates
xyz <- beamToXyz(adp)

# Convert to Earth coordinates
# The matlab file gives -18.7 for declination. 
# This website https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#declination
# gives a value of 19 or 19.3 depending on which standard you choose...
enu <- xyzToEnu(xyz, declination = dec)

create_nc_from_adp(enu, "../proc/downstream_deep_down_2018_enu")
