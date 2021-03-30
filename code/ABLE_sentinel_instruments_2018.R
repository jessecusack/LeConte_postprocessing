library(oce)
library(ncdf4)
source("ncwrite.R")

# Parameters
lat <- 56.8370392
lon <- -132.3574123

# Turbidity sensor
file <- "~/Dropbox/LeConte/Data/ocean/september2018/raw/moorings/ABLE_Sentinel/RBR_Virtuoso_Turbidity/054311_20180915_2116.rsk"
file <- path.expand(file)
ctd <- read.ctd(file)
ctd[['latitude']] <- lat
ctd[['longitude']] <- lon
write_moored_ctd(ctd, "../proc/ABLE_sentinel_RBRvirtuoso_2018.nc")

# Microcat
file <- "~/Dropbox/LeConte/Data/ocean/september2018/raw/moorings/ABLE_Sentinel/SBE/SBE37SM-RS232_03707818_2018_09_18.cnv"
file <- path.expand(file)
ctd <- read.ctd(file)
ctd[['latitude']] <- lat
ctd[['longitude']] <- lon
write_moored_ctd(ctd, "../proc/ABLE_sentinel_SBE37_2018.nc")
