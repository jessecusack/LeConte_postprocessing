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
# # Start
#
# First run the R notebook: `test_R_adp.ipynb`
#
# This generates the netcdf files.

# %%
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import utils
import scipy.stats as stats
from scipy.ndimage import label, gaussian_filter

# %%
enu = xr.open_dataset("enu.nc", decode_times=False)
enu["time"] = utils.POSIX_to_datetime(enu.time.values)

# %%
hvel_kwargs = dict(vmin=-0.3, vmax=0.3, cmap="coolwarm")
vvel_kwargs = dict(vmin=-0.1, vmax=0.1, cmap="coolwarm")

fig, axs = plt.subplots(4, 1, sharex=True, figsize=(22, 14))
enu.u.plot(ax=axs[0], **hvel_kwargs)
enu.v.plot(ax=axs[1], **hvel_kwargs)
enu.w.plot(ax=axs[2], **vvel_kwargs)
np.abs(enu.err).plot(ax=axs[3], vmin=0, vmax=0.2)

for ax in axs:
    ax.set_xlabel("")

# %%
fig, axs = plt.subplots(4, 1, sharex=True, figsize=(22, 14))
enu.q1.plot(ax=axs[0])
enu.q2.plot(ax=axs[1])
enu.q3.plot(ax=axs[2])
enu.q4.plot(ax=axs[3])

# %%
fig, axs = plt.subplots(4, 1, sharex=True, figsize=(22, 14))
enu.a1.plot(ax=axs[0])
enu.a2.plot(ax=axs[1])
enu.a3.plot(ax=axs[2])
enu.a4.plot(ax=axs[3])

# %%
fig, axs = plt.subplots(4, 1, sharex=True, figsize=(22, 14))
enu.g1.plot(ax=axs[0])
enu.g2.plot(ax=axs[1])
enu.g3.plot(ax=axs[2])
enu.g4.plot(ax=axs[3])


# %% [markdown]
# ## Identify surface or bottom of ice

# %%
def mode(x, **kwargs):
    mout = np.squeeze(stats.mode(x, axis=1)[0])
    return mout


# %%
# nroll = 7
dthresh = 50.
sidelobe_pct = 0.1  # Side lobe percentage
var_names = ["a1", "a2", "a3", "a4"]

dmingood = np.full_like(enu.time, 1e10)

for var in var_names:

    idxmax = enu[var].where(enu.distance > dthresh).argmax("distance")
    dmax = enu.distance[idxmax]
    dsl = (1 - sidelobe_pct)*enu.distance[idxmax]
    # dmax = dmax.where(dmax > dthresh)
    # dmode = dmax.rolling(time=nroll, min_periods=1, center=True).reduce(mode)
    
    dmingood = np.minimum(dmingood, dsl.values)
    
    
for var in var_names:

    idxmax = enu[var].where(enu.distance > dthresh).argmax("distance")
    dmax = enu.distance[idxmax]
    dsl = (1 - sidelobe_pct)*enu.distance[idxmax]

    fig, ax = plt.subplots(figsize=(22, 4))
    enu[var].plot(ax=ax)
    dsl.plot(ax=ax, color="k")
    dmax.plot(ax=ax, color="k", linestyle="--")
    ax.plot(enu.time, dmingood, "r")

# %% [markdown]
# Cut off surface.

# %%
# Make a new dataset without surface
enus = enu.copy()
# Loop over the 2D datavars
mask = enu.distance < xr.DataArray(dmingood, dims={"time": enu.time})

for var in enu.data_vars:
    if enu[var].dims == ('distance', 'time'):
        print(f"Masking {var}.")
        enus[var] = enu[var].where(mask)

# Remove distances where there is no good data
enus = enus.isel(distance=mask.any("time"))

# %%
hvel_kwargs = dict(vmin=-0.3, vmax=0.3, cmap="coolwarm")
vvel_kwargs = dict(vmin=-0.1, vmax=0.1, cmap="coolwarm")

fig, axs = plt.subplots(4, 1, sharex=True, figsize=(22, 14))
enus.u.plot(ax=axs[0], **hvel_kwargs)
enus.v.plot(ax=axs[1], **hvel_kwargs)
enus.w.plot(ax=axs[2], **vvel_kwargs)
np.abs(enus.err).plot(ax=axs[3], vmin=0, vmax=0.2)

for ax in axs:
    ax.set_xlabel("")

# %% [markdown]
# ## Error mask

# %%
errthresh = 0.05
enue = enus.copy()
egood = np.abs(enus.err) < errthresh

var_names = ["u", "v", "w", "err"]
for var in var_names:
    enue[var] = enus[var].where(egood)

# %%
hvel_kwargs = dict(vmin=-0.3, vmax=0.3, cmap="coolwarm")
vvel_kwargs = dict(vmin=-0.1, vmax=0.1, cmap="coolwarm")

fig, axs = plt.subplots(4, 1, sharex=True, figsize=(22, 14))
enue.u.plot(ax=axs[0], **hvel_kwargs)
enue.v.plot(ax=axs[1], **hvel_kwargs)
enue.w.plot(ax=axs[2], **vvel_kwargs)
np.abs(enue.err).plot(ax=axs[3], vmin=0, vmax=0.2)

for ax in axs:
    ax.set_xlabel("")

# %% [markdown]
# ## Correlation mask

# %%
qthresh = 40.

fig, axs = plt.subplots(4, 1, sharex=True, figsize=(22, 14))
enus.q1.plot(ax=axs[0])
enus.q1.plot.contour(ax=axs[0], levels=[qthresh], colors="r")
enus.q2.plot(ax=axs[1])
enus.q2.plot.contour(ax=axs[1], levels=[qthresh], colors="r")
enus.q3.plot(ax=axs[2])
enus.q3.plot.contour(ax=axs[2], levels=[qthresh], colors="r")
enus.q4.plot(ax=axs[3])
enus.q4.plot.contour(ax=axs[3], levels=[qthresh], colors="r")

# %%
qsum = enus.q1 + enus.q2 + enus.q3 + enus.q4

fig, ax = plt.subplots(1, 1, sharex=True, figsize=(22, 4))
qsum.plot(ax=ax)
qsum.plot.contour(ax=ax, levels=[300], colors="r")

# %% [markdown]
# ## Correlation mask using sum

# %%
qthresh = 350
enuq = enus.copy()
qgood = qsum > qthresh

var_names = ["u", "v", "w", "err"]
for var in var_names:
    enuq[var] = enus[var].where(qgood)


# %%
fig, ax = plt.subplots(1, 1, sharex=True, figsize=(22, 4))
qsum.plot(ax=ax)
qsum.plot.contour(ax=ax, levels=[qthresh], colors="r")

# %%
hvel_kwargs = dict(vmin=-0.3, vmax=0.3, cmap="coolwarm")
vvel_kwargs = dict(vmin=-0.1, vmax=0.1, cmap="coolwarm")

fig, axs = plt.subplots(4, 1, sharex=True, figsize=(22, 14))
enuq.u.plot(ax=axs[0], **hvel_kwargs)
enuq.v.plot(ax=axs[1], **hvel_kwargs)
enuq.w.plot(ax=axs[2], **vvel_kwargs)
np.abs(enuq.err).plot(ax=axs[3], vmin=0, vmax=0.2)

for ax in axs:
    ax.set_xlabel("")

# %% [markdown]
# ## Combined error and correlation mask

# %%
enuc = enus.copy()
mask = qgood & egood

var_names = ["u", "v", "w", "err"]
for var in var_names:
    enuc[var] = enus[var].where(mask)
    

# %%
hvel_kwargs = dict(vmin=-0.3, vmax=0.3, cmap="coolwarm")
vvel_kwargs = dict(vmin=-0.1, vmax=0.1, cmap="coolwarm")

fig, axs = plt.subplots(4, 1, sharex=True, figsize=(22, 14))
enuc.u.plot(ax=axs[0], **hvel_kwargs)
enuc.v.plot(ax=axs[1], **hvel_kwargs)
enuc.w.plot(ax=axs[2], **vvel_kwargs)
np.abs(enuc.err).plot(ax=axs[3], vmin=0, vmax=0.2)

for ax in axs:
    ax.set_xlabel("")

# %% [markdown]
# ## Blurring and labelling

# %%
# legood, num_features = label(~egood.values)
# print(f"There are {num_features} features.")

# %%
# # Find features with most points
# feature_n = np.arange(1, num_features + 1)
# points_in_feature = np.zeros_like(feature_n)
# col_var = np.ones_like(feature_n)
# row_var = np.ones_like(feature_n)
# n_surf_points = np.ones_like(feature_n)
# for f in feature_n:
#     infeature = legood == f
#     points_in_feature[f-1] = infeature.sum()
    
#     idx_feature = np.argwhere(infeature)
#     n_surf_points[f-1] = ((legood.shape[0] - 1) == idx_feature[:, 0]).sum()
    
# #     if points_in_feature[f-1] != 1:
        
# #         row_var[f-1] = idx_feature[:, 0].std()
# #         col_var[f-1] = idx_feature[:, 1].std()
    
# idx_sort = np.argsort(points_in_feature)
# points_in_feature = points_in_feature[idx_sort]
# feature_n = feature_n[idx_sort]

# #              width      / height
# # aspect_ratio = col_var / row_var
# surf_ratio = n_surf_points/points_in_feature


# %%
# fig, ax = plt.subplots()
# ax.semilogy(surf_ratio, '.')
# # ax.set_ylim(0, 10)

# %%
# fig, ax = plt.subplots()
# ax.plot(points_in_feature, '.-')

# %%
# fig, ax = plt.subplots()
# ax.loglog(points_in_feature, surf_ratio, '.')

# %%
# # Interesting regions with low aspect ratio and high number of points
# interesting_features = feature_n[(surf_ratio < 0.001) & (points_in_feature > 100)]
# print(interesting_features)

# %%
# fig, axs = plt.subplots(2, 1, figsize=(22, 7))
# axs[0].pcolormesh(egood)

# # for feature in interesting_features:
# #     toplot = np.full_like(legood, np.nan)
# #     toplot[legood == feature] = 1.
# #     axs[1].pcolormesh(toplot, cmap="plasma")

# %%
egood_filt = gaussian_filter(egood.values.astype(float), (1, 3))

# %%
fig, axs = plt.subplots(2, 1, figsize=(22, 7))
axs[0].pcolormesh(egood_filt)
axs[0].contour(egood_filt, [0.25], colors="r")
axs[1].pcolormesh(enus.v.values, cmap="coolwarm", vmin=-0.3, vmax=0.3)
axs[1].contour(egood_filt, [0.25], colors="r")

# %% [markdown]
# ## Blurred error mask including correlation

# %%
errthresh = 0.05  # Blur around these errors
errthresh_high = 0.05  # Always remove these errors
maskthresh = 0.35  # Blurred mask threshold
qthresh = 300
sigma = (2, 5)

qsum = enus.q1 + enus.q2 + enus.q3 + enus.q4
qgood = qsum > qthresh
enueb = enus.copy()
egood = np.abs(enus.err) < errthresh
egood_filt = gaussian_filter(egood.values.astype(float), sigma)
ebgood = (egood_filt > maskthresh) & (np.abs(enus.err) < errthresh_high) & qgood

var_names = ["u", "v", "w", "err"]
for var in var_names:
    enueb[var] = enus[var].where(ebgood)

# %%
fig, ax = plt.subplots(1, 1, figsize=(22, 3.5))
ax.pcolormesh(egood_filt)
ax.contour(egood_filt, [maskthresh], colors="r")
ax.contour(qgood, [0.5], colors="g")

# %%
hvel_kwargs = dict(vmin=-0.3, vmax=0.3, cmap="coolwarm")
vvel_kwargs = dict(vmin=-0.1, vmax=0.1, cmap="coolwarm")

fig, axs = plt.subplots(4, 1, sharex=True, figsize=(22, 14))
enueb.u.plot(ax=axs[0], **hvel_kwargs)
# enueb.v.plot(ax=axs[1], **hvel_kwargs)
enus.u.plot(ax=axs[1], **hvel_kwargs)
enueb.w.plot(ax=axs[2], **vvel_kwargs)
np.abs(enueb.err).plot(ax=axs[3], vmin=0, vmax=0.2)

for ax in axs:
    ax.set_xlabel("")

# %% [markdown]
# ## Zoom in on some events

# %%
import pytz
import datetime

# %%
tslice = slice(datetime.datetime(2018, 9, 6, 4, 0, tzinfo=pytz.utc), datetime.datetime(2018, 9, 6, 6, 0, tzinfo=pytz.utc))
enueb_ = enueb.sel(time=tslice)
enus_ = enus.sel(time=tslice)

fig, axs = plt.subplots(6, 1, sharex=True, figsize=(22, 22))
enueb_.u.plot(ax=axs[0], **hvel_kwargs)
enus_.u.plot(ax=axs[1], **hvel_kwargs)
enueb_.v.plot(ax=axs[2], **hvel_kwargs)
enus_.v.plot(ax=axs[3], **hvel_kwargs)
enueb_.w.plot(ax=axs[4], **vvel_kwargs)
enus_.w.plot(ax=axs[5], **vvel_kwargs)

fig, ax = plt.subplots(figsize=(18, 3.5))
enus_.pitch.plot(ax=ax)
enus_.rol.plot(ax=ax)
# enus_.heading.plot(ax=ax)
ax.set_xlim(tslice.start, tslice.stop)

# %%
