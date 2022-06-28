# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     notebook_metadata_filter: -jupytext.text_representation.jupytext_version
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: lcpp
#     language: python
#     name: lcpp
# ---

# %% [markdown]
# # Downstream deep mooring odd CTD

# %%
from scipy.io import loadmat
import xarray as xr
import numpy as np
import utils
import gsw

# %%
dat = loadmat("/Users/jmcusack/Dropbox/LeConte/Data/ocean/september2018/processed/moorings/Downstream_Deep/SBE/SBE37SM-RS232_03707819_2018_09_18.mat", squeeze_me=True)

# %%
# POSIX time to be consistent...
time = (utils.datenum_to_datetime(dat["dn"]).astype("datetime64") - np.datetime64("1970-01-01T00:00:00")).astype("timedelta64[s]").astype(float)
t = dat["temp"]
p = dat["pr"]
C = dat["cond"]
SP = gsw.SP_from_C(10*C, t, p)  # *10 to convert from S/m to mS/cm

# %%
ds = xr.open_dataset("../proc/downstream_deep_mooring/SBE37SM_SN10551.nc")

# %%
data_vars = dict(
    t = ("time", t.astype("float32"), dict(units="degree_C", standard_name="sea_water_temperature", long_name="Temperature")),
    SP = ("time", SP.astype("float32"), dict(standard_name="sea_water_practical_salinity", long_name="Practical salinity")),  
    p = ("time", p.astype("float32"), dict(units="dbar", standard_name="sea_water_pressure", long_name="Pressure")),
    C = ("time", C.astype("float32"), dict(units="S m-1", standard_name="sea_water_electrical_conductivity", long_name="Conductivity")),
    lon = ds.lon,
    lat = ds.lat,
)

coords = dict(
    time = ("time", time, dict(units="seconds since 1970-01-01 00:00:00", standard_name="time", long_name="Time")),
)

attrs = dict(
    pressureType="sea",
    serialNumber="07819",
    sampleInterval=6,
    sampleIntervalUnits="s",
)

ctd = xr.Dataset(data_vars, coords, attrs)

# %%
# From the R script
# from <- as.POSIXct("2018-09-03 02:17:00", tz="UTC")
# to <- as.POSIXct("2018-09-16 16:10:15", tz="UTC")

time_ = utils.datenum_to_datetime(dat["dn"]).astype("datetime64")
t0 = np.datetime64("2018-09-03T02:17:00")
t1 = np.datetime64("2018-09-16T16:10:15")
use = (time_ > t0) & (time_ < t1)

ctd = ctd.isel(time=use)

# %%
ctd.to_netcdf("../proc/downstream_deep_mooring/SBE37SM_SN07819.nc")
