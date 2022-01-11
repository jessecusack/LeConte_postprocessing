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
import dolfyn
from datetime import datetime
import matplotlib.pyplot as plt

# %%
dat = dolfyn.read("/Users/jmcusack/Dropbox/LeConte/Data/ocean/september2018/raw/moorings/ABLE_Deep/16670013.000")

# %%
dat = dat.drop_vars(["heading_std", "pitch_std", "roll_std"])

# %%
dat["time"] = [datetime.fromtimestamp(t) for t in dat.time.values]

# %%
dat

# %%
time_start = "2018-09-02T18:09"
time_end = "2018-09-20T19:25"

dat.amp.isel(beam=0, range=10).sel(time=slice(time_start, time_end)).plot()

# %%
dat = dat.sel(time=slice(time_start, time_end))

# %%
dat.amp.isel(beam=3).plot()

# %%
dat = dolfyn.set_declination(dat, 19.32)
# lat <- 56.835592
# lon <- -132.3572915

# %%
enu = dolfyn.rotate2(dat, "earth")

# %%
enu

# %%
# Remove excessive velocity
enu = dolfyn.adp.clean.vel_exceeds_thresh(enu, 1)

# %%
# Remove poor correlation
enu = dolfyn.adp.clean.correlation_filter(enu)

# %%
# Try and find surface... (good luck)
enu = dolfyn.adp.clean.find_surface(enu, nfilt=101)

# %%
enu.depth.plot()

# %%
# Vertical velocity is wrong... I think... 
enu.vel[2, :, :] *= -1

# %%
tslice = slice("2018-09-18T12:00", "2018-09-18T15:00")
rslice = slice(0, 140)

fig, axs = plt.subplots(3, 1, sharex=True, figsize=(18, 12))

vel = enu.vel.sel(time=tslice, range=rslice)

uvkwargs = dict(vmin=-0.2, vmax=0.2, cmap="RdBu_r")
wkwargs = dict(vmin=-0.1, vmax=0.1, cmap="RdBu_r")

vel.isel(dir=0).plot(ax=axs[0], **uvkwargs)
vel.isel(dir=1).plot(ax=axs[1], **uvkwargs)
vel.isel(dir=2).plot(ax=axs[2], **wkwargs)

# %%
fig, axs = plt.subplots(4, 1, sharex=True, figsize=(18, 16))

amp = enu.amp.sel(time=tslice, range=rslice)

akwargs = dict()

amp.isel(beam=0).plot(ax=axs[0], **akwargs)
amp.isel(beam=1).plot(ax=axs[1], **akwargs)
amp.isel(beam=2).plot(ax=axs[2], **akwargs)
amp.isel(beam=3).plot(ax=axs[3], **akwargs)

# %%
fig, axs = plt.subplots(4, 1, sharex=True, figsize=(18, 16))

corr = enu.corr.sel(time=tslice, range=rslice)

akwargs = dict()

corr.isel(beam=0).plot(ax=axs[0], **akwargs)
corr.isel(beam=1).plot(ax=axs[1], **akwargs)
corr.isel(beam=2).plot(ax=axs[2], **akwargs)
corr.isel(beam=3).plot(ax=axs[3], **akwargs)

# %%
enu.pitch.sel(time=tslice).plot()

# %%
enu["roll"].sel(time=tslice).plot()

# %%
