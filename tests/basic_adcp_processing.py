# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.0
#   kernelspec:
#     display_name: lcpp-dev
#     language: python
#     name: lcpp-dev
# ---

# %% [markdown]
# # Basic processing of Sentinel V mooring

# %%
import xarray as xr
import numpy as np
import utils
import matplotlib.pyplot as plt
import scipy.stats as stats


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
# Load datasets and do some basic conversion of times and variables.

# %%
sV = xr.open_dataset("../proc/ABLE_sentinel_2018_enu.nc")
sV = sV.set_coords(["lon", "lat"])
sV["time"] = utils.POSIX_to_datetime(sV.time.values)

virt = xr.open_dataset("../proc/ABLE_sentinel_RBRvirtuoso_2018.nc")
virt = virt.set_coords(["lon", "lat"])
virt["time"] = utils.POSIX_to_datetime(virt.time.values)

sbe = xr.open_dataset("../proc/ABLE_sentinel_SBE37_2018.nc")
sbe = sbe.set_coords(["lon", "lat"])
sbe["time"] = utils.POSIX_to_datetime(sbe.time.values)

# %% [markdown]
# Define some parameters and simple thresholds for processing.

# %%
pmin = 125  # Minimum pressure to keep
dpdtmax = 0.4e-9  # Maximum rate of change of pressure to keep
cut_ends = 2  # Number of points on either end to remove after applying other thresholds
dt = 10  # Bin size for time average [s]

# %% [markdown]
# Apply the thresholds to remove some data.

# %%
is_deep = sV.p > pmin
is_slow = np.fabs(sV.p.differentiate("time")) < dpdtmax
keep = is_deep & is_slow

sVp = sV.isel(time=keep).isel(time=slice(cut_ends, -cut_ends))

# %%
sVp.p.plot.line('.')

# %% [markdown]
# Bin average data to reduce size and errors. 
#
# First make bins.

# %%
# Time bin start and end to nearest minute. This will cut off some data.
tstart = (sVp.time[0].values + np.timedelta64(30, 's')).astype('datetime64[m]')
tend = sVp.time[-1].values.astype('datetime64[m]')

timebins = np.arange(tstart, tend, np.timedelta64(dt, 's'))

# %% [markdown]
# Group and take mean.

# %%
gb = sVp.groupby_bins("time", timebins)
sVa = gb.mean(skipna=True, keep_attrs=True)

# Use mid time as dimension, rather than Interval.
sVa["time_bins"] = interval_to_mid(sVa.time_bins.values).astype("datetime64[s]")
sVa = sVa.rename({"time_bins": "time"})

# %% [markdown]
# ## Cut off data above surface
#
# Use a simple echo intensity threshold to find the maximum.

# %%
dmin = 110.  # Minimum distance above which to look for the maximum
nroll = 120  # Number of points in rolling mode window
fcut = 0.1  # Extra distance to remove (1 - fcut)*dcut

# %%
sVa.va.isel(time=10000).plot.line('.')

# %% [markdown]
# Identify echo maximum in each beam, using a rolling mode to smooth out data.

# %%
fig, ax = plt.subplots()

dcuts = []

for var in ["a1", "a2", "a3", "a4", "va"]:

    am = sVa[var].where(sVa.distance > dmin)
    imax = am.argmax(dim="distance", skipna=True)
    dmax = am.distance[imax]

    ro = dmax.rolling(time=nroll, center=True)

    dm = ro.reduce(mode)

    dcut = (1 - fcut)*dm
    
    ax.plot(sVa.time, dmax, 'r')
    ax.plot(sVa.time, dm, 'orange')
    ax.plot(sVa.time, dcut, 'g')
    
    dcuts.append(dcut.values)

# %%
dcut_min = np.stack(dcuts, axis=1).min(axis=1)
dcut_min = xr.DataArray(dcut_min, dims={"time": sVa.time})

# %% [markdown]
# Mask and remove data above distance threshold.

# %%
sVm = sVa.where(sVa.distance < dcut_min)

sVm["p"] = sVa.p
sVm["t"] = sVa.t
sVm["pitch"] = sVa.pitch
sVm["roll"] = sVa.roll
sVm["heading"] = sVa.heading

sVm = sVm.isel(distance=~np.isnan(sVm.u).all(axis=0))

# %% [markdown]
# ## Plotting time series

# %%
timeslice = slice(np.datetime64("2018-09-05T08:00"), np.datetime64("2018-09-10T11:00"))

sVm_ = sVm.sel(time=timeslice)

fig, axs = plt.subplots(4, 1, figsize=(15, 10), sharex=True)

sVm_.u.plot(ax=axs[0], x="time", vmin=-0.2, vmax=0.2, cmap="coolwarm")
sVm_.v.plot(ax=axs[1], x="time", vmin=-0.2, vmax=0.2, cmap="coolwarm")
sVm_.w.plot(ax=axs[2], x="time", vmin=-0.2, vmax=0.2, cmap="coolwarm")
(-sVm_.vv).plot(ax=axs[3], x="time", vmin=-0.2, vmax=0.2, cmap="coolwarm")

fig, ax = plt.subplots(figsize=(12, 3))
sVm_.p.plot(ax=ax)

# %%
timeslice = slice(np.datetime64("2018-09-05T08:00"), np.datetime64("2018-09-10T11:00"))

sVm_ = sVm.sel(time=timeslice)

fig, axs = plt.subplots(8, 1, figsize=(15, 25), sharex=True)

sVm_.u.plot(ax=axs[0], x="time", vmin=-0.2, vmax=0.2, cmap="coolwarm")
sVm_.v.plot(ax=axs[1], x="time", vmin=-0.2, vmax=0.2, cmap="coolwarm")
(-sVm_.vv).plot(ax=axs[2], x="time", vmin=-0.2, vmax=0.2, cmap="coolwarm")
sVm_.a1.plot(ax=axs[3], x="time")
sVm_.a2.plot(ax=axs[4], x="time")
sVm_.a3.plot(ax=axs[5], x="time")
sVm_.a4.plot(ax=axs[6], x="time")
sVm_.va.plot(ax=axs[7], x="time")

# %% [markdown]
# # Plug in other instruments to dataset
#
# Group and bin average.

# %%
gb = virt.groupby_bins("time", timebins)
virta = gb.mean(skipna=True, keep_attrs=True)

# Use mid time as dimension, rather than Interval.
virta["time_bins"] = interval_to_mid(virta.time_bins.values).astype("datetime64[ms]")
virta = virta.rename({"time_bins": "time"})

gb = sbe.groupby_bins("time", timebins)
sbea = gb.mean(skipna=True, keep_attrs=True)

# Use mid time as dimension, rather than Interval.
sbea["time_bins"] = interval_to_mid(sbea.time_bins.values).astype("datetime64[ms]")
sbea = sbea.rename({"time_bins": "time"})

# %%
fig, ax = plt.subplots(figsize=(12, 3))
virta.turb.plot(ax=ax)

fig, axs = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
sbea.p.plot(ax=axs[0])
sbea.t.plot(ax=axs[1])
sbea.SP.plot(ax=axs[2])

# %%
