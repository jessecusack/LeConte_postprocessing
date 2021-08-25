# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: lcpp-dev
#     language: python
#     name: lcpp-dev
# ---

# %% [markdown]
# # ABLE deep processing

# %%
import xarray as xr
import numpy as np
import utils
import matplotlib.pyplot as plt
import scipy.stats as stats
import gsw
import utm


def mode(x, **kwargs):
    mout = np.squeeze(stats.mode(x, axis=1)[0])
    return mout


def interval_to_mid(intervals):
    """
    Parameters
    ----------
        intervals : 1D numpy array
            An array of pandas Interval objects.
    Returns
    -------
        mids : 1D numpy array
            Midpoints of the intervals.
    """
    return np.array([v.mid for v in intervals])


# %% [markdown]
# Load datasets.

# %%
ad = xr.open_dataset("../proc/ABLE_deep_2018_enu.nc")
ad = ad.set_coords(["lon", "lat"])
ad["time"] = utils.POSIX_to_datetime(ad.time.values)

x, y, *_ = utm.from_latlon(ad.lat, ad.lon)
ad = ad.assign_coords({"x": x, "y": y})

# Conflicts with roll operation
rol = ad["roll"]
ad = ad.drop_vars("roll")
ad["rol"] = (rol.dims, rol.values, rol.attrs)

# sbe51 = xr.open_dataset("../proc/downstream_deep_SBE37_10551_2018.nc")
# sbe51 = sbe51.set_coords(["lon", "lat"])
# sbe51["time"] = utils.POSIX_to_datetime(sbe51.time.values)



# %% [markdown]
# Define some parameters and simple thresholds for processing.

# %%
pmin = 120  # Minimum pressure to keep
cut_ends = 20  # Number of points on either end to remove after applying other thresholds
dt = 10  # Bin size for time average [s]

# %% [markdown]
# Apply the thresholds to remove some data.

# %%
is_deep = ad.p > pmin
keep = is_deep

adp = ad.isel(time=keep).isel(time=slice(cut_ends, -cut_ends))

# %%
adp.p.plot.line('.')

# %% [markdown]
# ## Quality control
#
# Note [Marion's document](https://escholarship.org/content/qt6xd149s8/qt6xd149s8.pdf)

# %%
# qc_err0 = 0.3
# qc_err1 = 0.5
qc_err = 0.15  # error velocity
qc_q = 110  # correlation
qc_uv = 2.0  # horizontal velocity
qc_w = 1.5  # vertical velocity
qc_a = 30 # echo intensity

# %%
qc_u_bad = np.abs(adp.u) > qc_uv
qc_v_bad = np.abs(adp.v) > qc_uv
qc_w_bad = np.abs(adp.w) > qc_w
qc_err_bad = np.abs(adp.err) > qc_err
qc_q1_good = adp.q1 > qc_q
qc_q2_good = adp.q2 > qc_q
qc_q3_good = adp.q3 > qc_q
qc_q4_good = adp.q4 > qc_q

qc_q_bad = (qc_q1_good.astype(int) + qc_q2_good.astype(int) + qc_q3_good.astype(int) + qc_q4_good.astype(int)) <= 3

# %%
uv_reject = (qc_q_bad.astype(int) + qc_err_bad.astype(int) + qc_u_bad.astype(int) + qc_v_bad.astype(int)) > 1
w_reject = (qc_q_bad.astype(int) + qc_err_bad.astype(int) + qc_w_bad.astype(int)) > 1

# %%
fig, axs = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(10, 8))
uv_reject.plot(ax=axs[0])
w_reject.plot(ax=axs[1])


# %% [markdown]
# Remove velocity using QC.

# %%
adqc = adp.copy()

u = adqc.u.values
u[uv_reject] = np.nan
adqc["u"] = (adqc.u.dims, u, adqc.u.attrs)

v = adqc.v.values
v[uv_reject] = np.nan
adqc["v"] = (adqc.v.dims, v, adqc.v.attrs)

w = adqc.w.values
w[w_reject] = np.nan
adqc["w"] = (adqc.w.dims, w, adqc.w.attrs)

# %% [markdown]
# ## Time binning
#
# Bin average data to reduce size and errors. 
#
# First make bins.

# %%
# Time bin start and end to nearest minute. This will cut off some data.
tstart = (adqc.time[0].values + np.timedelta64(30, 's')).astype('datetime64[m]')
tend = adqc.time[-1].values.astype('datetime64[m]')

timebins = np.arange(tstart, tend, np.timedelta64(dt, 's'))

# %% [markdown]
# Group and take mean.

# %%
gb = adqc.groupby_bins("time", timebins)
ada_ = gb.mean(skipna=True, keep_attrs=True)

# Use mid time as dimension, rather than Interval.
ada_["time_bins"] = interval_to_mid(ada_.time_bins.values).astype("datetime64[s]")
ada_ = ada_.rename({"time_bins": "time"})
# Mean of heading should be performed using circular mean. (Technically, so should pitch and roll, but for small angles the noncircular mean is ok)
ada_["heading"] = (["time"], ada_.heading.groupby_bins("time", timebins).reduce(stats.circmean, high=360.).values)

# %% [markdown]
# Remove all NaN colums.

# %%
allnan = np.isnan(ada_.a1).all(axis=1)
ada = ada_.isel(time=~allnan)

# %% [markdown]
# ## Cut off data above surface and below bottom
#
# Use a simple echo intensity threshold to find the maximum.

# %%
dmin = 130.  # Minimum distance above which to look for the maximum
nroll = 200  # Number of points in rolling mode window
fcut = 0.1  # Extra distance to remove (1 - fcut)*dcut

# %%
fig, ax = plt.subplots()
ada.a2.isel(time=60000).plot.line(ax=ax, marker='.')

# %% [markdown]
# Identify echo maximum in each beam, using a rolling mode to smooth out data.

# %%
fig, ax = plt.subplots()

dcuts = []

for var in ["a1", "a2", "a3", "a4"]:

    am = ada[var].where(ada.distance > dmin)
#     dmax = ada.p.copy()
#     dmax[:] = 0
#     dmax = dmax.astype(int)
#     allnan = np.isnan(am).all(axis=1)

    imax = am.argmax(dim="distance", skipna=True)
    dmax = am.distance[imax]

    ro = dmax.rolling(time=nroll, min_periods=1, center=True)

    dm = ro.reduce(mode)

    dcut = (1 - fcut)*dm
    
    ax.plot(ada.time, dmax, 'r')
    ax.plot(ada.time, dm, 'orange')
    ax.plot(ada.time, dcut, 'g')
    
    dcuts.append(dcut.values)

# %%
dcuts = np.stack(dcuts, axis=1)

# Use only the vertical beam for finding the surface.
dcut_min = dcuts.min(axis=1)
dcut_min = xr.DataArray(dcut_min, dims={"time": ada.time})

# %% [markdown]
# Mask and remove data above distance threshold.

# %%
adm = ada.where(ada.distance < dcut_min)

# The masking process converts some variables to 2D, change them back...
adm["p"] = ada.p
adm["t"] = ada.t
adm["pitch"] = ada.pitch
adm["rol"] = ada.rol
adm["heading"] = ada.heading

adm = adm.isel(distance=~np.isnan(adm.u.values).all(axis=0))

# %% [markdown]
# Thermodynamics.

# %%
adm["z"] = (adm.p.dims, gsw.z_from_p(adm.p, adm.lat), {"units": "m", "long_name": "height"})
adm["depth"] = (adm.p.dims, -adm.z, {"units": "m", "long_name": "depth"})

# %% [markdown]
# Save to netcdf.

# %%
adm

# %%
# adm.to_netcdf("../proc/ABLE_deep_mooring_2018.nc")
adm.to_netcdf("test.nc")

# %% [markdown]
# # Check some things...

# %%
adm.depth.plot()

# %%
import importlib
importlib.reload(utils)

# %%
fig, ax = plt.subplots(figsize=(15, 5))
ax.invert_yaxis()
ax.plot(adm.time, utils.butter_filter(adm.p, 1/3600, 1/10))
# ax.plot(utils.butter_filter(adm.p, 1/(3*86400), 1/10))
ax.set_ylabel("Pressure [dbar]")

fig, ax = plt.subplots(figsize=(15, 5))
ax.plot(adm.time, adm.pitch, '.')
ax.set_ylabel("pitch [deg]")

# %%
