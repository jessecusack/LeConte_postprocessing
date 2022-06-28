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
#     display_name: lcpp-dev
#     language: python
#     name: lcpp-dev
# ---

# %% [markdown]
# # Basic processing of the downstream deep mooring

# %%
import xarray as xr
import numpy as np
import utils
import matplotlib.pyplot as plt
import scipy.stats as stats
import gsw
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
# Load datasets.

# %%
up = xr.open_dataset("../proc/downstream_deep_upward_2018_enu.nc")
up = up.set_coords(["lon", "lat"])
up["time"] = utils.POSIX_to_datetime(up.time.values).astype(np.datetime64)

do = xr.open_dataset("../proc/downstream_deep_down_2018_enu.nc")
do = do.set_coords(["lon", "lat"])
do["time"] = utils.POSIX_to_datetime(do.time.values).astype(np.datetime64)

sbe51 = xr.open_dataset("../proc/downstream_deep_SBE37_10551_2018.nc")
sbe51 = sbe51.set_coords(["lon", "lat"])
sbe51["time"] = utils.POSIX_to_datetime(sbe51.time.values).astype(np.datetime64)

sbe52 = xr.open_dataset("../proc/downstream_deep_SBE37_10552_2018.nc")
sbe52 = sbe52.set_coords(["lon", "lat"])
sbe52["time"] = utils.POSIX_to_datetime(sbe52.time.values).astype(np.datetime64)

sbe53 = xr.open_dataset("../proc/downstream_deep_SBE37_10553_2018.nc")
sbe53 = sbe53.set_coords(["lon", "lat"])
sbe53["time"] = utils.POSIX_to_datetime(sbe53.time.values).astype(np.datetime64)

# %% [markdown]
# Define some parameters and simple thresholds for processing.

# %%
pmin = 80  # Minimum pressure to keep
dpdtmax = 0.4e-9  # Maximum rate of change of pressure to keep
cut_ends = 2  # Number of points on either end to remove after applying other thresholds
dt = 10  # Bin size for time average [s]

# %% [markdown]
# Apply the thresholds to remove some data.

# %%
is_deep = do.p > pmin
is_slow = np.fabs(do.p.differentiate("time")) < dpdtmax
keep = is_deep & is_slow

dop = do.isel(time=keep).isel(time=slice(cut_ends, -cut_ends))
# Trim the upward looking ADCP to the same times as the downward.
upp = up.sel(time=slice(dop.time[0], dop.time[-1]))

# %%
dop.p.plot.line('.')

# %% [markdown]
# ## New cut off data above surface
#
# Downward looking first.

# %%
dthresh = 50.
sidelobe_pct = 1 - np.cos(np.deg2rad(dop.beamAngle))
var_names = ["a1", "a2", "a3", "a4"]
nroll = 5

dmingood_down = np.full((dop.time.size, len(var_names)), np.nan)

fig, axs = plt.subplots(len(var_names), 1, figsize=(14, 3*len(var_names)))

for i, var in enumerate(var_names):
    idxmax = dop[var].where(dop.distance > dthresh).argmax("distance")
    dmax = dop.distance[idxmax]
    dsl = (1 - sidelobe_pct)*dop.distance[idxmax]
#     dmax = dmax.where(dmax > dthresh)
    dmode = dsl.rolling(time=nroll, min_periods=1, center=True).reduce(mode)
    
    dop[var].plot(ax=axs[i])
    
    dmingood_down[:, i] = dmode
    dsl.plot(ax=axs[i], color="r")
    axs[i].set_title("")
    
for i in range(len(var_names)):
    axs[i].plot(dop.time, dmingood_down.min(axis=1), color="k")

# %% [markdown]
# Do the masking of the downward looking variables.

# %%
good = dmingood_down.min(axis=1)

# Make a new dataset without surface
dos = dop.copy()
# Loop over the 2D datavars
mask = dop.distance < xr.DataArray(good, dims={"time": dop.time})

for var in dop.data_vars:
    if dop[var].dims == ('distance', 'time'):
        print(f"Masking {var}.")
        dos[var] = dop[var].where(mask)

# Remove distances where there is no good data
dos = dos.isel(distance=mask.any("time"))

# %% [markdown]
# Now the upward.

# %%
dthresh = 50.
sidelobe_pct = 1 - np.cos(np.deg2rad(upp.beamAngle))
var_names = ["a1", "a2", "a3", "a4"]
nroll = 5

dmingood_up = np.full((upp.time.size, len(var_names)), np.nan)

fig, axs = plt.subplots(len(var_names), 1, figsize=(14, 3*len(var_names)))

for i, var in enumerate(var_names):
    idxmax = upp[var].where(upp.distance > dthresh).argmax("distance")
    dmax = upp.distance[idxmax]
    dsl = (1 - sidelobe_pct)*upp.distance[idxmax]
#     dmax = dmax.where(dmax > dthresh)
    dmode = dsl.rolling(time=nroll, min_periods=1, center=True).reduce(mode)
    
    upp[var].plot(ax=axs[i])
    
    dmingood_up[:, i] = dmode
    dsl.plot(ax=axs[i], color="r")
    axs[i].set_title("")
    
for i in range(len(var_names)):
    axs[i].plot(upp.time, dmingood_up.min(axis=1), color="k")

# %% [markdown]
# And mask the upward data.

# %%
good = dmingood_up.min(axis=1)

# Make a new dataset without surface
ups = upp.copy()
# Loop over the 2D datavars
mask = upp.distance < xr.DataArray(good, dims={"time": upp.time})

for var in upp.data_vars:
    if upp[var].dims == ('distance', 'time'):
        print(f"Masking {var}.")
        ups[var] = upp[var].where(mask)

# Remove distances where there is no good data
ups = ups.isel(distance=mask.any("time"))

# %% [markdown]
# ## New quality control

# %%
# These parameters seem to work for both...
errthresh = 0.05  # Blur around these errors
errthresh_high = 0.05  # Always remove these errors
maskthresh = 0.35  # Blurred mask threshold
qthresh = 300
sigma = (2, 5)

# Downward
qsum = dos.q1 + dos.q2 + dos.q3 + dos.q4
qgood = qsum > qthresh
doqc = dos.copy()
egood = np.abs(dos.err) < errthresh
egood_filt = gaussian_filter(egood.values.astype(float), sigma)
ebgood = (egood_filt > maskthresh) & (np.abs(dos.err) < errthresh_high) & qgood

var_names = ["u", "v", "w", "err"]
for var in var_names:
    doqc[var] = doqc[var].where(ebgood)
    
# Upward
qsum = ups.q1 + ups.q2 + ups.q3 + ups.q4
qgood = qsum > qthresh
upqc = ups.copy()
egood = np.abs(ups.err) < errthresh
egood_filt = gaussian_filter(egood.values.astype(float), sigma)
ebgood = (egood_filt > maskthresh) & (np.abs(ups.err) < errthresh_high) & qgood

var_names = ["u", "v", "w", "err"]
for var in var_names:
    upqc[var] = upqc[var].where(ebgood)

# %% [markdown]
# ## Time bin averaging
#
# Bin average data to reduce size and errors. 
#
# First make bins.

# %%
# Time bin start and end to nearest minute. This will cut off some data.
tstart = (doqc.time[0].values + np.timedelta64(30, 's')).astype('datetime64[m]')
tend = doqc.time[-1].values.astype('datetime64[m]')

timebins = np.arange(tstart, tend, np.timedelta64(dt, 's'))

# %% [markdown]
# Group and take mean.

# %%
gbd = doqc.groupby_bins("time", timebins)
doa = gbd.mean(skipna=True, keep_attrs=True)

# Use mid time as dimension, rather than Interval.
doa["time_bins"] = interval_to_mid(doa.time_bins.values).astype("datetime64[s]")
doa = doa.rename({"time_bins": "time"})
# Mean of heading should be performed using circular mean. (Technically, so should pitch and roll, but for small angles the noncircular mean is ok)
doa["heading"] = (["time"], doa.heading.groupby_bins("time", timebins).reduce(stats.circmean, high=360.).values)

# %% [markdown]
# Same for upward.

# %%
gbu = upqc.groupby_bins("time", timebins)
upa = gbu.mean(skipna=True, keep_attrs=True)

# Use mid time as dimension, rather than Interval.
upa["time_bins"] = interval_to_mid(upa.time_bins.values).astype("datetime64[s]")
upa = upa.rename({"time_bins": "time"})
# Mean of heading should be performed using circular mean. (Technically, so should pitch and roll, but for small angles the noncircular mean is ok)
upa["heading"] = (["time"], upa.heading.groupby_bins("time", timebins).reduce(stats.circmean, high=360.).values)

# %% [markdown]
# Finally, thermodynamics.

# %%
doa["z"] = (doa.p.dims, gsw.z_from_p(doa.p, doa.lat), {"units": "m", "long_name": "height"})
doa["depth"] = (doa.p.dims, -doa.z, {"units": "m", "long_name": "depth"})

# %% [markdown]
# ## Old cut off data above surface and below bottom
#
# Use a simple echo intensity threshold to find the maximum.

# %%
# dmin = 50.  # Minimum distance above which to look for the maximum
# nroll = 180  # Number of points in rolling mode window
# fcut = 0.1  # Extra distance to remove (1 - fcut)*dcut

# %%
# fig, ax = plt.subplots()
# upa.a1.isel(time=10000).plot.line(ax=ax, marker='.')
# doa.a1.isel(time=10000).plot.line(ax=ax, marker='.')

# %% [markdown]
# Identify echo maximum in each beam, using a rolling mode to smooth out data.

# %%
# fig, ax = plt.subplots()

# dcuts = []

# for var in ["a1", "a2", "a3", "a4"]:

#     am = upa[var].where(upa.distance > dmin)
#     imax = am.argmax(dim="distance", skipna=True)
#     dmax = am.distance[imax]

#     ro = dmax.rolling(time=nroll, min_periods=1, center=True)

#     dm = ro.reduce(mode)

#     dcut = (1 - fcut)*dm
    
#     ax.plot(upa.time, dmax, 'r')
#     ax.plot(upa.time, dm, 'orange')
#     ax.plot(upa.time, dcut, 'g')
    
#     dcuts.append(dcut.values)

# %%
# dcuts = np.stack(dcuts, axis=1)

# # Use only the vertical beam for finding the surface.
# dcut_min = dcuts.min(axis=1)
# dcut_min = xr.DataArray(dcut_min, dims={"time": upa.time})

# %% [markdown]
# Mask and remove data above distance threshold.

# %%
# upm = upa.where(upa.distance < dcut_min)

# # The masking process converts some variables to 2D, change them back...
# upm["t"] = upa.t
# upm["pitch"] = upa.pitch
# upm["rol"] = upa.rol
# upm["heading"] = upa.heading

# upm = upm.isel(distance=~np.isnan(upm.u.values).all(axis=0))

# %% [markdown]
# And repeat the cutting of distance for the downward looking adcp.

# %%
# fig, ax = plt.subplots()

# dcuts = []

# for var in ["a1", "a2", "a3", "a4"]:

#     am = doa[var].where(doa.distance > dmin)
#     imax = am.argmax(dim="distance", skipna=True)
#     dmax = am.distance[imax]

#     ro = dmax.rolling(time=nroll, min_periods=1, center=True)

#     dm = ro.reduce(mode)

#     dcut = (1 - fcut)*dm
    
#     ax.plot(doa.time, dmax, 'r')
#     ax.plot(doa.time, dm, 'orange')
#     ax.plot(doa.time, dcut, 'g')
    
#     dcuts.append(dcut.values)

# %%
# dcuts = np.stack(dcuts, axis=1)

# # Use only the vertical beam for finding the surface.
# dcut_min = dcuts.min(axis=1)
# dcut_min = xr.DataArray(dcut_min, dims={"time": doa.time})

# %% [markdown]
# Mask and remove data above distance threshold.

# %%
# dom = doa.where(doa.distance < dcut_min)

# # The masking process converts some variables to 2D, change them back...
# dom["p"] = doa.p
# dom["t"] = doa.t
# dom["pitch"] = doa.pitch
# dom["rol"] = doa.rol
# dom["heading"] = doa.heading

# dom = dom.isel(distance=~np.isnan(dom.u.values).all(axis=0))

# %% [markdown]
# Thermodynamics.

# %%
# dom["z"] = (dom.p.dims, gsw.z_from_p(dom.p, dom.lat), {"units": "m", "long_name": "height"})
# dom["depth"] = (dom.p.dims, -dom.z, {"units": "m", "long_name": "depth"})

# %% [markdown]
# # Plug in other instruments to dataset
#
# Group and bin average.

# %%
# 51
gb = sbe51.groupby_bins("time", timebins)
sbe51a = gb.mean(skipna=True, keep_attrs=True)

# Use mid time as dimension, rather than Interval.
sbe51a["time_bins"] = interval_to_mid(sbe51a.time_bins.values).astype("datetime64[ms]")
sbe51a = sbe51a.rename({"time_bins": "time"})

# 52
gb = sbe52.groupby_bins("time", timebins)
sbe52a = gb.mean(skipna=True, keep_attrs=True)

# Use mid time as dimension, rather than Interval.
sbe52a["time_bins"] = interval_to_mid(sbe52a.time_bins.values).astype("datetime64[ms]")
sbe52a = sbe52a.rename({"time_bins": "time"})

# 53
gb = sbe53.groupby_bins("time", timebins)
sbe53a = gb.mean(skipna=True, keep_attrs=True)

# Use mid time as dimension, rather than Interval.
sbe53a["time_bins"] = interval_to_mid(sbe53a.time_bins.values).astype("datetime64[ms]")
sbe53a = sbe53a.rename({"time_bins": "time"})

# %% [markdown]
# Estimate thermodynamic quantities.

# %%
sbe51a["SA"] = (sbe51a.p.dims, gsw.SA_from_SP(sbe51a.SP, sbe51a.p, sbe51a.lon, sbe51a.lat), {"units": "g/kg", "long_name": "Absolute_salinity"})
sbe51a["CT"] = (sbe51a.p.dims, gsw.CT_from_t(sbe51a.SA, sbe51a.t, sbe51a.p), {"units": "deg C", "long_name": "Conservative_temperature"})
sbe51a["z"] = (sbe51a.p.dims, gsw.z_from_p(sbe51a.p, sbe51a.lat), {"units": "m", "long_name": "height"})
sbe51a["depth"] = (sbe51a.p.dims, -sbe51a.z, {"units": "m", "long_name": "depth"})

sbe52a["SA"] = (sbe52a.p.dims, gsw.SA_from_SP(sbe52a.SP, sbe52a.p, sbe52a.lon, sbe52a.lat), {"units": "g/kg", "long_name": "Absolute_salinity"})
sbe52a["CT"] = (sbe52a.p.dims, gsw.CT_from_t(sbe52a.SA, sbe52a.t, sbe52a.p), {"units": "deg C", "long_name": "Conservative_temperature"})
sbe52a["z"] = (sbe52a.p.dims, gsw.z_from_p(sbe52a.p, sbe52a.lat), {"units": "m", "long_name": "height"})
sbe52a["depth"] = (sbe52a.p.dims, -sbe52a.z, {"units": "m", "long_name": "depth"})

sbe53a["SA"] = (sbe53a.p.dims, gsw.SA_from_SP(sbe53a.SP, sbe53a.p, sbe53a.lon, sbe53a.lat), {"units": "g/kg", "long_name": "Absolute_salinity"})
sbe53a["CT"] = (sbe53a.p.dims, gsw.CT_from_t(sbe53a.SA, sbe53a.t, sbe53a.p), {"units": "deg C", "long_name": "Conservative_temperature"})
sbe53a["z"] = (sbe53a.p.dims, gsw.z_from_p(sbe53a.p, sbe53a.lat), {"units": "m", "long_name": "height"})
sbe53a["depth"] = (sbe53a.p.dims, -sbe53a.z, {"units": "m", "long_name": "depth"})

# %% [markdown]
# Look at a couple of plots.

# %%
fig, axs = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
sbe51a.p.plot(ax=axs[0], label="SBE37")
doa.p.plot(ax=axs[0], label="ADCP")
sbe51a.t.plot(ax=axs[1])
# sbe51a.CT.plot(ax=axs[1])
sbe51a.SP.plot(ax=axs[2])
# sbe51a.SA.plot(ax=axs[2])

axs[0].legend()

fig, ax = plt.subplots()
sbe51a.p.plot(ax=ax)
sbe52a.p.plot(ax=ax)
sbe53a.p.plot(ax=ax)

# %% [markdown]
# # Create dataset

# %% [markdown]
# Depth of the ADCPs

# %%
# Nominal depth of the downward looking ADCP

# shallowest
depth_sbe_52 = float(sbe52a.depth.mean().values)
depth_upward = float(doa.depth.mean().values - 2.)  # Kind of a guess...
depth_sbe_51 = float(sbe51a.depth.mean().values)
depth_downward = float(doa.depth.mean().values)
depth_sbe_53 = float(sbe53a.depth.mean().values)
# deepest
depth_instrument = [depth_sbe_52, depth_upward, depth_sbe_51, depth_downward, depth_sbe_53]
instruments = ["SBE37", "ADCP", "SBE37", "ADCP", "SBE37"]
# Stack point data
p_stack = np.stack((sbe52a.p.values, np.nan*doa.p.values, sbe51a.p.values, doa.p.values, sbe53a.p.values), axis=0)
t_stack = np.stack((sbe52a.t.values, upa.t.values, sbe51a.t.values, doa.t.values, sbe53a.t.values), axis=0)
SP_stack = np.stack((sbe52a.SP.values, np.nan*sbe52a.SP.values, sbe51a.SP.values, np.nan*sbe52a.SP.values, sbe53a.SP.values), axis=0)
heading_stack = np.stack((np.nan*doa.heading.values, upa.heading.values, np.nan*doa.heading.values, doa.heading.values, np.nan*doa.heading.values), axis=0)
pitch_stack = np.stack((np.nan*doa.pitch.values, upa.pitch.values, np.nan*doa.pitch.values, doa.pitch.values, np.nan*doa.pitch.values), axis=0)
rol_stack = np.stack((np.nan*doa.rol.values, upa.rol.values, np.nan*doa.rol.values, doa.rol.values, np.nan*doa.rol.values), axis=0)

# Stack ranged data
depth_adcp = np.hstack(((depth_upward - upa.distance).values[::-1], (depth_downward + doa.distance).values))
distance = np.hstack((upa.distance.values[::-1], -doa.distance.values))  # Negative distance means 
u_stack = np.hstack((upa.u.values[:, ::-1], doa.u.values))
v_stack = np.hstack((upa.v.values[:, ::-1], doa.v.values))
w_stack = np.hstack((upa.w.values[:, ::-1], doa.w.values))

x, y, zone_number, zone_letter = utm.from_latlon(doa.lat, doa.lon)

coords = {
    "time": (["time"], doa.time.values),
    "depth_adcp": (["depth_adcp"], depth_adcp, {"Variable": "Nominal depth of ADCP measurement [m]"}),
    "distance": (["depth_adcp"], distance, {"Variable": "ADCP bin distances from instrument [m]"}),
    "depth_instrument": (["instrument"], depth_instrument, {"Variable": "Nominal depth of each instrument [m]"}),
    "instrument": (["instrument"], instruments),
    "lon": ([], doa.lon.values),
    "lat": ([], doa.lat.values),
    "x": ([], x),
    "y": ([], y),
    "zone_letter": ([], zone_letter),
    "zone_number": ([], zone_number),
}

datavars = {
    "p": (["time", "instrument"], p_stack.T, {"Variable": "Pressure [dbar]"}),
    "t": (["time", "instrument"], t_stack.T, {"Variable": "Temperature (in situ) [deg C]"}),
    "SP": (["time", "instrument"], SP_stack.T, {"Variable": "Practical salinity [PSU]"}),
    "u": (["time", "depth_adcp"], u_stack, {"Variable": "Eastward velocity [m s-1]"}),
    "v": (["time", "depth_adcp"], v_stack, {"Variable": "Northward velocity [m s-1]"}),
    "w": (["time", "depth_adcp"], w_stack, {"Variable": "Vertical velocity [m s-1]"}),
    "heading": (["time", "instrument"], heading_stack.T, {"Variable": "ADCP heading [deg]"}),
    "pitch": (["time", "instrument"], pitch_stack.T, {"Variable": "ADCP pitch [deg]"}),
    "rol": (["time", "instrument"], rol_stack.T, {"Variable": "ADCP roll [deg]"}),
}

ds = xr.Dataset(datavars, coords)

# Stuff some other ADCP variables
for i, var in enumerate(["a1", "a2", "a3", "a4", "q1", "q2", "q3", "q4", "err"]):
    var_stack = np.hstack((upa[var].values[:, ::-1], doa[var].values))
    ds[var] = (ds.u.dims, var_stack, doa[var].attrs)


# %% [markdown]
# Extra thermodynamic stuff.

# %%
ds["SA"] = (ds.p.dims, gsw.SA_from_SP(ds.SP, ds.p, ds.lon, ds.lat), {"Variable": "Absolute_salinity [g kg-1]"})
ds["CT"] = (ds.p.dims, gsw.CT_from_t(ds.SA, ds.t, ds.p), {"Variable": "Conservative_temperature [deg C]"})
ds["z"] = (ds.p.dims, gsw.z_from_p(ds.p, ds.lat), {"Variable": "height [m]"})
ds["depth"] = (ds.p.dims, -ds.z, {"Variable": "depth [m]"})

# %% [markdown]
# Save.

# %%
ds.to_netcdf("../proc/downstream_deep_mooring_2018.nc")

# %%
