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

# %%
import xarray as xr
import matplotlib.pyplot as plt
import utils

# %% [markdown]
# ## Downstream deep upward looking

# %%
ddu = xr.open_dataset("../proc/downstream_deep_upward_2018_enu.nc") #, time=slice(9900, -2300)).sel(distance=slice(0, 75))
ddu["time"] = utils.POSIX_to_datetime(ddu.time.values)

# %%
ddu

# %%
ddu.u.plot(vmin=-0.5, vmax=0.5, cmap="coolwarm")

# %%
ddu.v.plot(vmin=-0.5, vmax=0.5, cmap="coolwarm")

# %%
ddu.a1.plot()

# %%
ddu.g1.plot()

# %% [markdown]
# ## Downstream deep downward looking

# %%
ddd = xr.open_dataset("../proc/downstream_deep_down_2018_enu.nc") #, time=slice(9500, None)).sel(distance=slice(0, 70))
ddd["time"] = utils.POSIX_to_datetime(ddd.time.values)

# %%
ddd

# %%
ddd.v.plot(vmin=-0.5, vmax=0.5, cmap="coolwarm")

# %%
ddd.a4.plot()

# %%
ddd.heading.plot()

# %% [markdown]
# ## Combining upward and downward

# %%
fig, ax = plt.subplots()
ax.plot(ddu.time, '.')
ax.plot(ddd.time, '.')

# %% [markdown]
# ## ABLE Sentinel

# %%
as_ = xr.open_dataset("../proc/ABLE_sentinel_2018_enu.nc") # .sel(distance=slice(0, 125))
#as_ = as_.isel(time=(as_.p > 125).values)
as_["time"] = utils.POSIX_to_datetime(as_.time.values)

# %%
as_

# %%
as_.p.plot()

# %%
as_.t.plot()

# %%
# as__ = as_.isel(time=slice(500, 1500))

fig, axs = plt.subplots(4, 1, figsize=(14, 9))
as_.u.plot(ax=axs[0], vmin=-0.2, vmax=0.2, cmap="coolwarm")
as_.v.plot(ax=axs[1], vmin=-0.2, vmax=0.2, cmap="coolwarm")
as_.w.plot(ax=axs[2], vmin=-0.2, vmax=0.2, cmap="coolwarm")
as_.vv.plot(ax=axs[3], vmin=-0.2, vmax=0.2, cmap="coolwarm")

# %%
fig, axs = plt.subplots(4, 1, figsize=(14, 8))
as_.q1.plot(ax=axs[0])
as_.q2.plot(ax=axs[1])
as_.q3.plot(ax=axs[2])
as_.q4.plot(ax=axs[3])

# %% [markdown]
# ## ABLE Deep

# %%
ad = xr.open_dataset("../proc/ABLE_deep_2018_enu.nc", decode_times=False)
ad["time"] = utils.POSIX_to_datetime(ad.time.values)

# %%
ad

# %%
fig, axs = plt.subplots(3, 1, figsize=(14, 7))
ad.u.plot(ax=axs[0], vmin=-0.2, vmax=0.2, cmap="coolwarm")
ad.v.plot(ax=axs[1], vmin=-0.2, vmax=0.2, cmap="coolwarm")
ad.w.plot(ax=axs[2], vmin=-0.2, vmax=0.2, cmap="coolwarm")

# %%
fig, axs = plt.subplots(4, 1, figsize=(14, 8))
ad.a1.plot(ax=axs[0])
ad.a2.plot(ax=axs[1])
ad.a3.plot(ax=axs[2])
ad.a4.plot(ax=axs[3])

# %%
ad.t.plot()

# %%
ad.p.plot()

# %%
ad["roll"].plot()

# %%
ad.heading.plot()

# %%
fig, axs = plt.subplots(4, 1, figsize=(14, 8))
ad.g1.plot(ax=axs[0])
ad.g2.plot(ax=axs[1])
ad.g3.plot(ax=axs[2])
ad.g4.plot(ax=axs[3])
