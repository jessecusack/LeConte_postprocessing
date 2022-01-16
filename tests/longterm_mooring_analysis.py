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

# %%
import xarray as xr
import utils
import numpy as np
import matplotlib.pyplot as plt
import utm
import gsw

# %% [markdown]
# # Period 1
#
# ## Combine ADP and pressure
#
# SBE37: 7891 was located 2.5 m below the ADP

# %%
a1 = xr.open_dataset("../proc/long_term_deep_1_enu.nc")
a1 = a1.set_coords(["lon", "lat"])
a1["time"] = utils.POSIX_to_datetime(a1.time.values).astype(np.datetime64)

x, y, *_ = utm.from_latlon(a1.lat, a1.lon)
a1 = a1.assign_coords({"x": x, "y": y})

c1 = xr.open_dataset("../proc/long_term_deep_1_SBE37_7819.nc")
c1["time"] = utils.POSIX_to_datetime(c1.time.values).astype(np.datetime64)

# %% [markdown]
# Apply thermodynamics.

# %%
c1["SA"] = (c1.p.dims, gsw.SA_from_SP(c1.SP, c1.p, c1.lon, c1.lat).data, {"units": "g/kg", "long_name": "Absolute_salinity"})
c1["CT"] = (c1.p.dims, gsw.CT_from_t(c1.SA, c1.t, c1.p).data, {"units": "deg C", "long_name": "Conservative_temperature"})
c1["z"] = (c1.p.dims, gsw.z_from_p(c1.p, c1.lat).data, {"units": "m", "long_name": "height"})
c1["depth"] = (c1.p.dims, -c1.z.data, {"units": "m", "long_name": "depth"})
c1 = c1.set_coords(["z", "depth"])

# %%
c1.p.plot()

# %% [markdown]
# Estimate ADP depth (2.5 m above CTD)

# %%
a1["depth"] = (["time"], (c1.depth - 2.5).interp(dict(time=a1.time)).data, {"units": "m", "long_name": "depth"})
a1["z"] = (["time"], -a1.depth.data, {"units": "m", "long_name": "height"})
a1 = a1.set_coords(["depth", "z"])

# %%
a1

# %%
sidelobe = 1 - np.cos(np.deg2rad(a1.beamAngle))

a1_ = a1.isel(time=slice(0, 50000))

fig, ax = plt.subplots()
a1_.a1.plot(ax=ax)
a1_.depth.plot(ax=ax)
((1 - sidelobe)*a1_.depth).plot(ax=ax)
ax.set_ylabel("Distance [m]")
ax.set_xlabel("")

# %% [markdown]
# Mask above surface

# %%
a1s = a1.copy()

mask = a1.distance < (1 - sidelobe)*a1.depth

for var in a1.data_vars:
    if a1[var].dims == ('distance', 'time'):
        print(f"Masking {var}.")
        a1s[var] = a1[var].where(mask)

# Remove distances where there is no good data
a1s = a1s.isel(distance=mask.any("time"))

# %%
a1s.g3.plot()


# %%
a1s.u.isel(time=slice(300000, 310000)).plot(vmin=-0.10, vmax=0.10, cmap="RdBu_r")

# %%
a1s.u.isel(time=slice(300000, 310000)).mean("time").plot()
a1s.v.isel(time=slice(300000, 310000)).mean("time").plot()

# %% [markdown]
# # Old stuff

# %%
# tslice = slice("2017-02-20T04:00", "2017-02-20T12:00")
# ds = ds.sel(distance=slice(0, 93))

# %%
fig, axs = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(20, 12))

wkwargs = dict(vmin=-0.03, vmax=0.03, cmap="RdBu_r")
uvkwargs = dict(vmin=-0.10, vmax=0.10, cmap="RdBu_r")

ds.u.sel(time=tslice).rolling(time=13, min_periods=2, center=True).mean().plot(ax=axs[0], **uvkwargs)
ds.v.sel(time=tslice).rolling(time=13, min_periods=2, center=True).mean().plot(ax=axs[1], **uvkwargs)
ds.w.sel(time=tslice).rolling(time=13, min_periods=2, center=True).mean().plot(ax=axs[2], **wkwargs)

# %%
ds.a4.sel(time=tslice).rolling(time=5, min_periods=2, center=True).mean().plot(figsize=(20, 4))

# %%
ds.rol.sel(time=tslice).rolling(time=5, min_periods=2, center=True).mean().plot(figsize=(20, 4))

# %%
bathy = xr.open_dataset("../proc/bathy_sep_2017.nc")

# %%
fig, ax = plt.subplots()
bathy.H.plot(ax=ax)
ax.plot(ds.lon, ds.lat, 'ro')

# %%
tslice = slice("2015-05-01T00:00", "2017-01-15T00:00")
dslice = slice(0, 90)
tlow = 30 # min
thigh = 5 # min

ulow = ds.u.sel(time=tslice, distance=dslice).rolling(time=int(tlow*2), min_periods=tlow, center=True).mean()
vlow = ds.v.sel(time=tslice, distance=dslice).rolling(time=int(tlow*2), min_periods=tlow, center=True).mean()
wlow = ds.w.sel(time=tslice, distance=dslice).rolling(time=int(tlow*2), min_periods=tlow, center=True).mean()

uanom = ds.u.sel(time=tslice, distance=dslice).rolling(time=int(thigh*2), min_periods=thigh, center=True).mean() - ulow
vanom = ds.v.sel(time=tslice, distance=dslice).rolling(time=int(thigh*2), min_periods=thigh, center=True).mean() - vlow
wanom = ds.w.sel(time=tslice, distance=dslice).rolling(time=int(thigh*2), min_periods=thigh, center=True).mean() - wlow



# %%
troll = 24 # hr

fig, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 8))

(uanom**2 + vanom**2 + wanom**2).mean("distance").rolling(time=int(troll*3600//30), min_periods=100, center=True).mean().plot(ax=axs[0])
# axs[0].set_ylim(0, 0.0001)

(ulow**2 + vlow**2).mean("distance").rolling(time=int(troll*3600//30), min_periods=100, center=True).mean().plot(ax=axs[1])

fig, ax = plt.subplots()
ulow.plot(ax=ax)

# %%
# tslice = slice("2016-08-15T12:00", "2016-08-15T23:00")
tslice = slice("2016-09-15T03:00", "2016-09-15T21:00")

fig, axs = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(20, 12))

wkwargs = dict(vmin=-0.03, vmax=0.03, cmap="RdBu_r")
uvkwargs = dict(vmin=-0.10, vmax=0.10, cmap="RdBu_r")

uanom.sel(time=tslice).plot(ax=axs[0], **uvkwargs)
vanom.sel(time=tslice).plot(ax=axs[1], **uvkwargs)
wanom.sel(time=tslice).plot(ax=axs[2], **wkwargs)

fig, axs = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(20, 12))

wkwargs = dict(vmin=-0.03, vmax=0.03, cmap="RdBu_r")
uvkwargs = dict(vmin=-0.10, vmax=0.10, cmap="RdBu_r")

ulow.sel(time=tslice).plot(ax=axs[0], **uvkwargs)
vlow.sel(time=tslice).plot(ax=axs[1], **uvkwargs)
wlow.sel(time=tslice).plot(ax=axs[2], **wkwargs)

# %%
