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
# Combine and clean long term moored ADP records.

# %%
import os
from datetime import datetime

import gsw
import numpy as np
import pandas as pd
import utm
import xarray as xr

import utils

# %% [markdown]
# # Process odd CTD

# %%
print("Processing odd CTD")

lat = 56.8265
lon = -132.3767
lastgood = 445939  # last line number of good data
skiprows = 295

file = "~/Dropbox/LeConte/Data/ocean/september2017/moorings/deep/raw/7819/SBE37SM-RS232_03707819_2017_09_18.asc"
tb = pd.read_csv(
    file,
    usecols=[0, 1, 2, 3, 5],
    header=0,
    names=["time", "p", "C", "t", "SP"],
    index_col=0,
    skiprows=skiprows,
    nrows=lastgood - skiprows,
    parse_dates=True,
    date_parser=lambda x: datetime.strptime(x, "%d %b %Y %H:%M:%S"),
)
c1 = xr.Dataset.from_dataframe(tb)
c1 = c1.isel(time=c1.p > 90)
c1 = c1.assign_coords(dict(lon=lon, lat=lat))
c1.to_netcdf("../proc/long_term_moorings/long_term_deep_3_SBE37_7819.nc")

# %% [markdown]
# # Common processing

# %%
# additional sidelobe removal multiplier, e.g. retained data = (1 - fside*sidelobe_percentage)
fside = 1.25
# maximum valid speed [m s-1]
spdmax = 1.5
# minimum data fraction, minimum data in a bin needed to keep
fdatamin = 0.1

# %% [markdown]
# # Data files

# %%
adp1 = "../proc/long_term_moorings/long_term_deep_1_enu.nc"
adp2 = "../proc/long_term_moorings/long_term_deep_2_enu.nc"
adp3 = "../proc/long_term_moorings/long_term_deep_3_enu.nc"
ctd1 = "../proc/long_term_moorings/long_term_deep_1_SBE37_7819.nc"
ctd2 = "../proc/long_term_moorings/long_term_deep_2_SBE37_7819.nc"
ctd3 = "../proc/long_term_moorings/long_term_deep_3_SBE37_7819.nc"

adps = [adp1, adp2, adp3]
ctds = [ctd1, ctd2, ctd3]

# %% [markdown]
# # Processing Loop

# %%
for i, (adp, ctd) in enumerate(zip(adps, ctds)):

    print(f"Loading {adp}")
    print(f"Loading {ctd}")

    # Load data
    a1 = xr.open_dataset(adp)
    a1 = a1.set_coords(["lon", "lat"])
    a1["time"] = utils.POSIX_to_datetime(a1.time.values).astype(np.datetime64)

    x, y, *_ = utm.from_latlon(a1.lat, a1.lon)
    a1 = a1.assign_coords({"x": x, "y": y})

    c1 = xr.open_dataset(ctd)
    try:
        c1["time"] = utils.POSIX_to_datetime(c1.time.values).astype(np.datetime64)
    except TypeError:
        print("No need to convert time")

    # Estimate CTD depth
    c1["z"] = (
        c1.p.dims,
        gsw.z_from_p(c1.p, c1.lat).data,
        {"units": "m", "long_name": "height"},
    )
    c1["depth"] = (c1.p.dims, -c1.z.data, {"units": "m", "long_name": "depth"})
    c1 = c1.set_coords(["z", "depth"])

    # Estimate ADP depth, 2.5 m above the CTD
    a1["depth"] = (
        ["time"],
        (c1.depth - 2.5).interp(dict(time=a1.time)).data,
        {"units": "m", "long_name": "depth"},
    )
    a1["z"] = (["time"], -a1.depth.data, {"units": "m", "long_name": "height"})
    a1 = a1.set_coords(["depth", "z"])

    # Sidelope fraction
    sidelobe = 1 - np.cos(np.deg2rad(a1.beamAngle))

    # Masking of bad data inc. surface, high speed, bins with low data fraction
    a1s = a1.copy()

    mask = a1.distance < (1 - fside * sidelobe) * a1.depth  # Mask sidelobes
    mask = mask.data
    spd = np.sqrt(a1.u ** 2 + a1.v ** 2 + a1.w ** 2)
    mask[(spd > spdmax).data] = False  # Mask high speeds

    for var in a1.data_vars:
        if a1[var].dims == ("distance", "time"):
            print(f"Masking {var}.")
            a1s[var] = a1[var].where(mask)

    fdata = mask.sum(axis=1) / mask.shape[1]

    # Remove distances where there is no good data
    a1s = a1s.isel(distance=fdata > fdatamin)  # 1 is the time axis

    # Save cleaned dataset
    save_name = f"long_term_deep_{i+1}_combo.nc"
    save_path = os.path.join("../proc/long_term_moorings/", save_name)
    print(f"Saving {save_path}")
    a1s.to_netcdf(save_path)
