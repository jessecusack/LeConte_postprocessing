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

# %%
ds = xr.open_dataset("enu.nc")
ds = ds.set_coords(["lon", "lat"])
ds["time"] = utils.POSIX_to_datetime(ds.time.values).astype(np.datetime64)

# %%
ds

# %%
ds.time.plot()

# %%
tslice = slice("2017-02-20T04:00", "2017-02-20T12:00")
ds = ds.sel(distance=slice(0, 93))

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
