library(oce)
library(ncdf4)
source("ncwrite.R")

# Parameters
lat <- 56.83596667
lon <- -132.3616833

# Microcat 10551
file <- "~/Dropbox/LeConte/Data/ocean/september2018/raw/moorings/Downstream_Deep/SBE/SBE37SM-RS232_03710551_2018_09_18.cnv"
file <- path.expand(file)
paste("Reading:", file)
ctd <- read.ctd(file)
ctd[['latitude']] <- lat
ctd[['longitude']] <- lon
summary(ctd)
write_moored_ctd(ctd, "../proc/downstream_deep_SBE37_10551_2018.nc")

# Microcat 10552
file <- "~/Dropbox/LeConte/Data/ocean/september2018/raw/moorings/Downstream_Deep/SBE/SBE37SM-RS232_03710552_2018_09_18.cnv"
file <- path.expand(file)
paste("Reading:", file)
ctd <- read.ctd(file)
ctd[['latitude']] <- lat
ctd[['longitude']] <- lon
summary(ctd)
write_moored_ctd(ctd, "../proc/downstream_deep_SBE37_10552_2018.nc")

# Microcat 10553
file <- "~/Dropbox/LeConte/Data/ocean/september2018/raw/moorings/Downstream_Deep/SBE/SBE37SM-RS232_03710553_2018_09_18.cnv"
file <- path.expand(file)
paste("Reading:", file)
ctd <- read.ctd(file)
ctd[['latitude']] <- lat
ctd[['longitude']] <- lon
summary(ctd)
write_moored_ctd(ctd, "../proc/downstream_deep_SBE37_10553_2018.nc")
