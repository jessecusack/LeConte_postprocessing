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
# # Basic processing of Sentinel V mooring

# %%
import xarray as xr
import numpy as np
import utils
import matplotlib.pyplot as plt
import scipy.stats as stats
import utm
from scipy.ndimage import gaussian_filter


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
sV["time"] = utils.POSIX_to_datetime(sV.time.values).astype(np.datetime64)

x, y, *_ = utm.from_latlon(sV.lat, sV.lon)
sV = sV.assign_coords({"x": x, "y": y})

virt = xr.open_dataset("../proc/ABLE_sentinel_RBRvirtuoso_2018.nc")
virt = virt.set_coords(["lon", "lat"])
virt["time"] = utils.POSIX_to_datetime(virt.time.values).astype(np.datetime64)

sbe = xr.open_dataset("../proc/ABLE_sentinel_SBE37_2018.nc")
sbe = sbe.set_coords(["lon", "lat"])
sbe["time"] = utils.POSIX_to_datetime(sbe.time.values).astype(np.datetime64)

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
# ## Old quality control
#
# Note [Marion's document](https://escholarship.org/content/qt6xd149s8/qt6xd149s8.pdf)

# %%
# # qc_err0 = 0.3
# # qc_err1 = 0.5
# qc_err = 0.15  # error velocity
# qc_q = 110  # correlation
# qc_uv = 2.0  # horizontal velocity
# qc_w = 1.5  # vertical velocity
# qc_a = 30 # echo intensity

# %%
# qc_u_bad = np.abs(sVp.u) > qc_uv
# qc_v_bad = np.abs(sVp.v) > qc_uv
# qc_w_bad = np.abs(sVp.w) > qc_w
# qc_vv_bad = np.abs(sVp.vv) > qc_w
# qc_err_bad = np.abs(sVp.err) > qc_err
# qc_q1_good = sVp.q1 > qc_q
# qc_q2_good = sVp.q2 > qc_q
# qc_q3_good = sVp.q3 > qc_q
# qc_q4_good = sVp.q4 > qc_q

# qc_q_bad = (qc_q1_good.astype(int) + qc_q2_good.astype(int) + qc_q3_good.astype(int) + qc_q4_good.astype(int)) <= 3

# %%
# uv_reject = (qc_q_bad.astype(int) + qc_err_bad.astype(int) + qc_u_bad.astype(int) + qc_v_bad.astype(int)) > 1
# w_reject = (qc_q_bad.astype(int) + qc_err_bad.astype(int) + qc_w_bad.astype(int)) > 1
# vv_reject = (qc_q_bad.astype(int) + qc_err_bad.astype(int) + qc_vv_bad.astype(int)) > 1

# %%
# fig, axs = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(10, 10))
# uv_reject.plot(ax=axs[0])
# w_reject.plot(ax=axs[1])
# vv_reject.plot(ax=axs[2])


# %% [markdown]
# Remove velocity using QC.

# %%
# sVqc = sVp.copy()

# u = sVqc.u.values
# u[uv_reject] = np.nan
# sVqc["u"] = (sVqc.u.dims, u, sVqc.u.attrs)

# v = sVqc.v.values
# v[uv_reject] = np.nan
# sVqc["v"] = (sVqc.v.dims, v, sVqc.v.attrs)

# w = sVqc.w.values
# w[w_reject] = np.nan
# sVqc["w"] = (sVqc.w.dims, w, sVqc.w.attrs)

# vv = sVqc.vv.values
# vv[vv_reject] = np.nan
# sVqc["vv"] = (sVqc.vv.dims, vv, sVqc.vv.attrs)

# %% [markdown]
# ## New cut off data above surface

# %%
dthresh = 100.
sidelobe_pct = 1 - np.cos(np.deg2rad(sVp.beamAngle))
var_names = ["a1", "a2", "a3", "a4", "va"]
nroll = 5

dmingood = np.full((sVp.time.size, len(var_names)), np.nan)

fig, axs = plt.subplots(len(var_names), 1, figsize=(14, 3*len(var_names)))

for i, var in enumerate(var_names):
    idxmax = sVp[var].where(sVp.distance > dthresh).argmax("distance")
    dmax = sVp.distance[idxmax]
    dsl = (1 - sidelobe_pct)*sVp.distance[idxmax]
#     dmax = dmax.where(dmax > dthresh)
    dmode = dsl.rolling(time=nroll, min_periods=1, center=True).reduce(mode)
    
    sVp[var].plot(ax=axs[i])
    
    dmingood[:, i] = dmode
    dsl.plot(ax=axs[i], color="r")
    axs[i].set_title("")
    
for i in range(len(var_names)):
    axs[i].plot(sVp.time, dmingood.min(axis=1), color="k")

# %%
good = dmingood.min(axis=1)

# Make a new dataset without surface
sVs = sVp.copy()
# Loop over the 2D datavars
mask = sVp.distance < xr.DataArray(good, dims={"time": sVp.time})

for var in sVp.data_vars:
    if sVp[var].dims == ('distance', 'time'):
        print(f"Masking {var}.")
        sVs[var] = sVp[var].where(mask)

# Remove distances where there is no good data
sVs = sVs.isel(distance=mask.any("time"))

# %% [markdown]
# ## New quality control

# %%
errthresh = 0.2  # Blur around these errors
errthresh_high = 0.2  # Always remove these errors
maskthresh = 0.35  # Blurred mask threshold
qthresh = 300
vqthresh = 35
sigma = (2, 5)

qsum = sVs.q1 + sVs.q2 + sVs.q3 + sVs.q4
qgood = qsum > qthresh
vqgood = sVs.vq.values > vqthresh
sVqc = sVs.copy()
egood = np.abs(sVs.err) < errthresh
egood_filt = gaussian_filter(egood.values.astype(float), sigma)
ebgood = (egood_filt > maskthresh) & (np.abs(sVs.err) < errthresh_high) & qgood
vebgood = (egood_filt > maskthresh) & vqgood

var_names = ["u", "v", "w", "err"]
for var in var_names:
    sVqc[var] = sVs[var].where(ebgood)
    
    
sVqc["vv"] = sVs.vv.where(vebgood)

# %% [markdown]
# ## Time binning

# %% [markdown]
# Bin average data to reduce size and errors. 
#
# First make bins.

# %%
# Time bin start and end to nearest minute. This will cut off some data.
tstart = (sVqc.time[0].values + np.timedelta64(30, 's')).astype('datetime64[m]')
tend = sVqc.time[-1].values.astype('datetime64[m]')

timebins = np.arange(tstart, tend, np.timedelta64(dt, 's'))

# %% [markdown]
# Group and take mean.

# %%
gb = sVqc.groupby_bins("time", timebins)
sVa = gb.mean(skipna=True, keep_attrs=True)

# Use mid time as dimension, rather than Interval.
sVa["time_bins"] = interval_to_mid(sVa.time_bins.values).astype("datetime64[s]")
sVa = sVa.rename({"time_bins": "time"})

# %% [markdown]
# Mean of heading should be performed using circular mean. (Technically, so should pitch and roll, but for small angles the noncircular mean is ok)

# %%
sVa["heading"] = (["time"], sVqc.heading.groupby_bins("time", timebins).reduce(stats.circmean, high=360.).values)

# %% [markdown]
# ## Old cut off data above surface
#
# Use a simple echo intensity threshold to find the maximum.

# %%
# dmin = 60.  # Minimum distance above which to look for the maximum
# nroll = 120  # Number of points in rolling mode window
# fcut = 0.1  # Extra distance to remove (1 - fcut)*dcut

# %%
# sVa.va.isel(time=10000).plot.line('.')

# %% [markdown]
# Identify echo maximum in each beam, using a rolling mode to smooth out data.

# %%
# # fig, ax = plt.subplots()

# dcuts = []

# for var in ["a1", "a2", "a3", "a4", "va"]:

#     am = sVa[var].where(sVa.distance > dmin)
#     imax = am.argmax(dim="distance", skipna=True)
#     dmax = am.distance[imax]

#     ro = dmax.rolling(time=nroll, min_periods=1, center=True)

#     dm = ro.reduce(mode)

#     dcut = (1 - fcut)*dm
    
# #     ax.plot(sVa.time, dmax, 'r')
# #     ax.plot(sVa.time, dm, 'orange')
# #     ax.plot(sVa.time, dcut, 'g')
    
#     dcuts.append(dcut.values)

# %%
# dcuts = np.stack(dcuts, axis=1)

# # Use only the vertical beam for finding the surface.
# dcut_min = dcuts[:, 4]
# dcut_min = xr.DataArray(dcut_min, dims={"time": sVa.time})

# %% [markdown]
# Mask and remove data above distance threshold.

# %%
# sVm = sVa.where(sVa.distance < dcut_min)

# # The masking process converts some variables to 2D, change them back...
# sVm["p"] = sVa.p
# sVm["t"] = sVa.t
# sVm["pitch"] = sVa.pitch
# sVm["rol"] = sVa.rol
# sVm["heading"] = sVa.heading

# sVm = sVm.isel(distance=~np.isnan(sVm.u).all(axis=0))

# %% [markdown]
# ## Plotting time series

# %%
timeslice = slice(np.datetime64("2018-09-05T08:00"), np.datetime64("2018-09-10T11:00"))

sVm_ = sVm.sel(time=timeslice)

fig, axs = plt.subplots(4, 1, figsize=(15, 10), sharex=True)

sVm_.u.plot(ax=axs[0], x="time", vmin=-0.2, vmax=0.2, cmap="coolwarm")
sVm_.v.plot(ax=axs[1], x="time", vmin=-0.2, vmax=0.2, cmap="coolwarm")
sVm_.w.plot(ax=axs[2], x="time", vmin=-0.2, vmax=0.2, cmap="coolwarm")
sVm_.vv.plot(ax=axs[3], x="time", vmin=-0.2, vmax=0.2, cmap="coolwarm")

fig, ax = plt.subplots(figsize=(12, 3))
sVm_.p.plot(ax=ax)

# %%
timeslice = slice(np.datetime64("2018-09-05T08:00"), np.datetime64("2018-09-10T11:00"))

sVm_ = sVm.sel(time=timeslice)

fig, axs = plt.subplots(8, 1, figsize=(15, 25), sharex=True)

sVm_.u.plot(ax=axs[0], x="time", vmin=-0.2, vmax=0.2, cmap="coolwarm")
sVm_.v.plot(ax=axs[1], x="time", vmin=-0.2, vmax=0.2, cmap="coolwarm")
sVm_.vv.plot(ax=axs[2], x="time", vmin=-0.2, vmax=0.2, cmap="coolwarm")
sVm_.a1.plot(ax=axs[3], x="time")
sVm_.a2.plot(ax=axs[4], x="time")
sVm_.a3.plot(ax=axs[5], x="time")
sVm_.a4.plot(ax=axs[6], x="time")
sVm_.va.plot(ax=axs[7], x="time")

fig, axs = plt.subplots(3, 1, figsize=(11, 8))
sVm_.heading.plot(ax=axs[0])
sVm_.rol.plot(ax=axs[1])
sVm_.pitch.plot(ax=axs[2])

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

# %% [markdown]
# Look at a couple of plots.

# %%
fig, ax = plt.subplots(figsize=(12, 3))
virta.turb.plot(ax=ax)

fig, axs = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
sbea.p.plot(ax=axs[0])
sbea.t.plot(ax=axs[1])
sbea.SP.plot(ax=axs[2])

# %% [markdown]
# Assign other data to the sentinal dataset.

# %%
ds = sVa.copy()

# %%
ds["turb_RBR"] = (sVa.p.dims, virta.turb, virta.turb.attrs)
ds["SP_SBE37"] = (sVa.p.dims, sbea.SP, sbea.SP.attrs)
ds["C_SBE37"] = (sVa.p.dims, sbea.C, sbea.C.attrs)
ds["t_SBE37"] = (sVa.p.dims, sbea.t, sbea.t.attrs)
ds["p_SBE37"] = (sVa.p.dims, sbea.p, sbea.p.attrs)

# %% [markdown]
# Try a plot...

# %%
fig, ax = plt.subplots()
ds.p_SBE37.plot(ax=ax)
ds.p.plot(ax=ax, yincrease=False)

# %% [markdown]
# Estimate some more thermodynamic variables.

# %%
import gsw

# %%
ds["SA_SBE37"] = (ds.p.dims, gsw.SA_from_SP(ds.SP_SBE37, ds.p_SBE37, ds.lon, ds.lat), {"units": "g/kg", "long_name": "Absolute_salinity"})
ds["CT_SBE37"] = (ds.p.dims, gsw.CT_from_t(ds.SA_SBE37, ds.t_SBE37, ds.p_SBE37), {"units": "deg C", "long_name": "Conservative_temperature"})
ds["z_SBE37"] = (ds.p.dims, gsw.z_from_p(ds.p_SBE37, ds.lat), {"units": "m", "long_name": "height"})
ds["depth_SBE37"] = (ds.p.dims, -ds.z_SBE37, {"units": "m", "long_name": "depth"})
ds["z_ADCP"] = (ds.p.dims, gsw.z_from_p(ds.p, ds.lat), {"units": "m", "long_name": "height"})
ds["depth_ADCP"] = (ds.p.dims, -ds.z_ADCP, {"units": "m", "long_name": "depth"})
ds["z"] = (ds.distance.dims, ds.distance + ds.z_ADCP.mean(dim="time"), {"units": "m", "long_name": "height"})
ds["depth"] = (ds.distance.dims, -ds.z, {"units": "m", "long_name": "depth"})
ds = ds.set_coords(["z", "depth"])

# %% [markdown]
# Save dataset to netcdf.

# %%
ds.to_netcdf("../proc/ABLE_sentinel_mooring_2018.nc")

# %% [markdown]
# ## Examine a short segment of the dataset

# %%
timeslice = slice(np.datetime64("2018-09-05T08:00"), np.datetime64("2018-09-05T12:00"))

ds_ = ds.sel(time=timeslice)

fig, axs = plt.subplots(4, 1, figsize=(15, 10), sharex=True, sharey=True)
ds_.u.plot(ax=axs[0], y="depth", x="time", yincrease=False, vmin=-0.2, vmax=0.2, cmap="coolwarm")
ds_.a3.plot(ax=axs[1], y="depth", x="time", yincrease=False)
ds_.vv.plot(ax=axs[2], y="depth", x="time", yincrease=False, vmin=-0.2, vmax=0.2, cmap="coolwarm")
ds_.va.plot(ax=axs[3], y="depth", x="time", yincrease=False)

fig, axs = plt.subplots(4, 1, figsize=(11.7, 10), sharex=True)
ds_.p_SBE37.plot(ax=axs[0])
ds_.CT_SBE37.plot(ax=axs[1])
ds_.turb_RBR.plot(ax=axs[2])
ds_.pitch.plot(ax=axs[3])

# %% [markdown]
# Compare echo intensity near bottom for different beams.

# %%
dist = 5
timeslice = slice(np.datetime64("2018-09-05T08:00"), np.datetime64("2018-09-05T12:00"))

ds_ = ds.sel(time=timeslice).sel(distance=dist, method="nearest")

fig, ax = plt.subplots(figsize=(11, 4))
ds_.a1.plot(ax=ax, label="beam 1")
ds_.a2.plot(ax=ax, label="beam 2")
ds_.a3.plot(ax=ax, label="beam 3")
ds_.a4.plot(ax=ax, label="beam 4")
ds_.va.plot(ax=ax, label="beam v")
ax.set_ylabel("Echo intensity")
ax.legend()

# %%
timeslice = slice(np.datetime64("2018-09-05T08:00"), np.datetime64("2018-09-05T12:00"))
ds_ = ds.sel(time=timeslice)

fig, ax = plt.subplots(figsize=(10, 10))

for i in range(0, ds_.time.size, 50):
    ds__ = ds_.isel(time=i)
    ds__.va.plot(ax=ax, label=ds__.time.values.astype("datetime64[s]"))
    
ax.legend(loc="upper left", bbox_to_anchor=(1, 1))

# %%
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np


fig = plt.figure(figsize=(10, 10))
ax = fig.gca(projection='3d')


# def cc(arg):
#     return mcolors.to_rgba(arg, alpha=0.6)

xs = ds_.distance.values
verts = []

zs = []
for i in range(0, ds_.time.size, 100):
    ds__ = ds_.isel(time=i)
    time = (ds__.time - ds_.time[0]).astype(float)/1e9
    zs.append(time)
    ys = ds__.va.values
    ys[0], ys[-1] = 0, 0
    verts.append(list(zip(xs, ys)))

# zs = [0.0, 1.0, 2.0, 3.0]
# for z in zs:
#     ys = np.random.rand(len(xs))
#     ys[0], ys[-1] = 0, 0
#     verts.append(list(zip(xs, ys)))

poly = PolyCollection(verts)  # facecolors=[cc('r'), cc('g'), cc('b'), cc('y')]
                                         
poly.set_alpha(0.2)
ax.add_collection3d(poly, zs=zs, zdir='y')

ax.set_xlabel('Distance')
ax.set_xlim3d(0, xs.max())
ax.set_ylabel('Y')
ax.set_ylim3d(0, zs[-1])
ax.set_zlabel('Z')
ax.set_zlim3d(0, 200)

ax.view_init(elev=30., azim=30)

plt.show()

# %%
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

timeslice = slice(np.datetime64("2018-09-05T10:00"), np.datetime64("2018-09-05T10:45"))
ds_ = ds.sel(time=timeslice)


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection='3d')

T, D = np.meshgrid(ds_.distance.values, (ds_.time.values - ds_.time[0].values).astype(float)/1e9)

# Plot a basic wireframe.
ax.plot_wireframe(T, D, ds_.a2.values, rstride=1, cstride=1)

ax.view_init(elev=45., azim=120)


# %% [markdown]
# # New QC

# %%
tslice = slice(np.datetime64("2018-09-07T10:00"), np.datetime64("2018-09-07T11:00"))
# tslice = slice(np.datetime64("2018-09-04T10:00"), np.datetime64("2018-09-04T11:00"))
# tslice = slice(np.datetime64("2018-09-11T14:00"), np.datetime64("2018-09-11T16:00"))
# tslice = slice(np.datetime64("2018-09-10T03:00"), np.datetime64("2018-09-10T04:00")) 
enu = sVp.sel(time=tslice)

# %%
hvel_kwargs = dict(vmin=-0.3, vmax=0.3, cmap="coolwarm")
vvel_kwargs = dict(vmin=-0.1, vmax=0.1, cmap="coolwarm")

fig, axs = plt.subplots(5, 1, sharex=True, figsize=(22, 17))
enu.u.plot(ax=axs[0], **hvel_kwargs)
enu.v.plot(ax=axs[1], **hvel_kwargs)
enu.w.plot(ax=axs[2], **vvel_kwargs)
enu.vv.plot(ax=axs[3], **vvel_kwargs)
np.abs(enu.err).plot(ax=axs[4], vmin=0, vmax=0.2)

for ax in axs:
    ax.set_xlabel("")

# %%
fig, axs = plt.subplots(5, 1, sharex=True, figsize=(22, 17))
enu.q1.plot(ax=axs[0])
enu.q2.plot(ax=axs[1])
enu.q3.plot(ax=axs[2])
enu.q4.plot(ax=axs[3])
enu.vq.plot(ax=axs[4])

for ax in axs:
    ax.set_xlabel("")

# %%
dthresh = 100.
sidelobe_pct = 1 - np.cos(np.deg2rad(enu.beamAngle))
var_names = ["a1", "a2", "a3", "a4", "va"]
nroll = 5

dmingood = np.full((enu.time.size, len(var_names)), np.nan)

fig, axs = plt.subplots(len(var_names), 1, figsize=(14, 3*len(var_names)))

for i, var in enumerate(var_names):
    idxmax = enu[var].where(enu.distance > dthresh).argmax("distance")
    dmax = sVp.distance[idxmax]
    dsl = (1 - sidelobe_pct)*enu.distance[idxmax]
#     dmax = dmax.where(dmax > dthresh)
    dmode = dsl.rolling(time=nroll, min_periods=1, center=True).reduce(mode)
    
    enu[var].plot(ax=axs[i])
    
    dmingood[:, i] = dmode
    dmode.plot(ax=axs[i], color="r")
    axs[i].set_title("")
    
for i in range(len(var_names)):
    axs[i].plot(enu.time, dmingood.min(axis=1), color="k")

# %%
fig, axs = plt.subplots(3, 1, figsize=(22, 9))
enu.heading.plot(ax=axs[0], marker='.', linestyle="")
enu.rol.plot(ax=axs[1])
enu.pitch.plot(ax=axs[2])

# %%
# Make a new dataset without surface
enus = enu.copy()
# Loop over the 2D datavars
mask = enu.distance < xr.DataArray(dmingood.min(axis=1), dims={"time": enu.time})

for var in enu.data_vars:
    if enu[var].dims == ('distance', 'time'):
        print(f"Masking {var}.")
        enus[var] = enu[var].where(mask)

# Remove distances where there is no good data
enus = enus.isel(distance=mask.any("time"))

# %%
hvel_kwargs = dict(vmin=-0.3, vmax=0.3, cmap="coolwarm")
vvel_kwargs = dict(vmin=-0.1, vmax=0.1, cmap="coolwarm")

fig, axs = plt.subplots(5, 1, sharex=True, figsize=(22, 17))
enus.u.plot(ax=axs[0], **hvel_kwargs)
enus.v.plot(ax=axs[1], **hvel_kwargs)
enus.w.plot(ax=axs[2], **vvel_kwargs)
enus.vv.plot(ax=axs[3], **vvel_kwargs)
np.abs(enus.err).plot(ax=axs[4], vmin=0, vmax=0.2)

for ax in axs:
    ax.set_xlabel("")

# %%
from scipy.ndimage import gaussian_filter

# %%
errthresh = 0.2  # Blur around these errors
errthresh_high = 0.2  # Always remove these errors
maskthresh = 0.35  # Blurred mask threshold
qthresh = 300
vqthresh = 35
sigma = (2, 5)

qsum = enus.q1 + enus.q2 + enus.q3 + enus.q4
qgood = qsum > qthresh
vqgood = enus.vq.values > vqthresh
enueb = enus.copy()
egood = np.abs(enus.err) < errthresh
egood_filt = gaussian_filter(egood.values.astype(float), sigma)
ebgood = (egood_filt > maskthresh) & (np.abs(enus.err) < errthresh_high) & qgood
vebgood = (egood_filt > maskthresh) & vqgood

var_names = ["u", "v", "w", "err"]
for var in var_names:
    enueb[var] = enus[var].where(ebgood)
    
    
enueb["vv"] = enus.vv.where(vebgood)

# %%
fig, ax = plt.subplots(1, 1, figsize=(22, 3.5))
ax.pcolormesh(egood_filt)
ax.contour(egood_filt, [maskthresh], colors="r")
ax.contour(qgood, [0.5], colors="g")
ax.contour(vqgood, [0.5], colors="b")

# %% tags=[]
hvel_kwargs = dict(vmin=-0.3, vmax=0.3, cmap="coolwarm")
vvel_kwargs = dict(vmin=-0.1, vmax=0.1, cmap="coolwarm")

fig, axs = plt.subplots(8, 1, sharex=True, figsize=(22, 28))
enueb.u.plot(ax=axs[0], **hvel_kwargs)
enus.u.plot(ax=axs[1], **hvel_kwargs)
enueb.v.plot(ax=axs[2], **hvel_kwargs)
enus.v.plot(ax=axs[3], **hvel_kwargs)
enueb.w.plot(ax=axs[4], **vvel_kwargs)
enus.w.plot(ax=axs[5], **vvel_kwargs)
enueb.vv.plot(ax=axs[6], **vvel_kwargs)
enus.vv.plot(ax=axs[7], **vvel_kwargs)

for ax in axs:
    ax.set_xlabel("")

# %% [markdown]
# # Beam separation

# %%
z = sVp.distance[sVp.distance < 120]
angle = np.deg2rad(sVp.beamAngle)

separation_opposite = 2*z*np.tan(angle)
separation_adjacent = 2*z*np.tan(angle)*np.cos(np.pi/4)

fig, ax = plt.subplots()
ax.plot(separation_opposite, z, label="opposite")
ax.plot(separation_adjacent, z, label="adjacent")
ax.axvline(75, color="k", label="half wavelength")
ax.legend()
ax.grid()

ax.set_xlabel("Beam separation [m]")
ax.set_ylabel("Distance from ADCP (mast) [m]")
