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
import matplotlib.pyplot as plt
import numpy as np
import utils

# %%
dsb = xr.open_dataset("../proc/ABLE_deep_2018_beam.nc")
dsb["time"] = (utils.POSIX_to_datetime(dsb.time.values).astype(np.datetime64))

dsm = xr.open_dataset("../proc/ABLE_deep_2018_beam_mapped.nc")
dsm["time"] = (utils.POSIX_to_datetime(dsm.time.values).astype(np.datetime64))

dse = xr.open_dataset("../proc/ABLE_deep_2018_enu.nc")
dse["time"] = (utils.POSIX_to_datetime(dse.time.values).astype(np.datetime64))

# %% [markdown]
# # ABLE DEEP context

# %%
fig, ax = plt.subplots(figsize=(13, 5))
dsb.p.plot(ax=ax)
dsb.p.rolling(time=200, center=True).mean().plot(ax=ax)
ax.invert_yaxis()

fig, ax = plt.subplots(figsize=(13, 5))
dsb.pitch.plot(ax=ax)

fig, ax = plt.subplots(figsize=(13, 5))
dsb.rol.plot(ax=ax)

fig, ax = plt.subplots(figsize=(13, 5))
dsb.heading.plot(ax=ax, marker=".", linestyle="")

# %% [markdown]
# ## Beam coordinates

# %%
tslice = slice("2018-09-04T10:00", "2018-09-04T11:00")
ds_ = dsb.sel(time=tslice)

uvkwargs = dict(vmin=-0.2, vmax=0.2, cmap="RdBu_r", add_colorbar=False)
echokwargs = dict(add_colorbar=False, vmin=0, vmax=232)
corkwargs = dict(add_colorbar=False, vmin=0, vmax=143)

fig, axs = plt.subplots(15, 1, sharey=False, sharex=True, figsize=(20, 50))
ds_.v1.plot(ax=axs[0], **uvkwargs)
ds_.a1.plot(ax=axs[1], **echokwargs)
ds_.q1.plot(ax=axs[2], **corkwargs)
ds_.v2.plot(ax=axs[3], **uvkwargs)
ds_.a2.plot(ax=axs[4], **echokwargs)
ds_.q2.plot(ax=axs[5], **corkwargs)
ds_.v3.plot(ax=axs[6], **uvkwargs)
ds_.a3.plot(ax=axs[7], **echokwargs)
ds_.q3.plot(ax=axs[8], **corkwargs)
ds_.v4.plot(ax=axs[9], **uvkwargs)
ds_.a4.plot(ax=axs[10], **echokwargs)
ds_.q4.plot(ax=axs[11], **corkwargs)
ds_.heading.plot(ax=axs[12], marker=".", linestyle="")
ds_.pitch.plot(ax=axs[13], marker=".", linestyle="")
ds_.rol.plot(ax=axs[14], marker=".", linestyle="")

# %% [markdown]
# ## Bin mapped coordinates

# %%
tslice = slice("2018-09-04T10:00", "2018-09-04T11:00")
ds_ = dsm.sel(time=tslice)

uvkwargs = dict(vmin=-0.2, vmax=0.2, cmap="RdBu_r", add_colorbar=False)
echokwargs = dict(add_colorbar=False, vmin=0, vmax=232)
corkwargs = dict(add_colorbar=False, vmin=0, vmax=143)

fig, axs = plt.subplots(15, 1, sharey=False, sharex=True, figsize=(20, 50))
ds_.v1.plot(ax=axs[0], **uvkwargs)
ds_.a1.plot(ax=axs[1], **echokwargs)
ds_.q1.plot(ax=axs[2], **corkwargs)
ds_.v2.plot(ax=axs[3], **uvkwargs)
ds_.a2.plot(ax=axs[4], **echokwargs)
ds_.q2.plot(ax=axs[5], **corkwargs)
ds_.v3.plot(ax=axs[6], **uvkwargs)
ds_.a3.plot(ax=axs[7], **echokwargs)
ds_.q3.plot(ax=axs[8], **corkwargs)
ds_.v4.plot(ax=axs[9], **uvkwargs)
ds_.a4.plot(ax=axs[10], **echokwargs)
ds_.q4.plot(ax=axs[11], **corkwargs)
ds_.heading.plot(ax=axs[12], marker=".", linestyle="")
ds_.pitch.plot(ax=axs[13], marker=".", linestyle="")
ds_.rol.plot(ax=axs[14], marker=".", linestyle="")

# %% [markdown]
# ## ENU coordinates

# %%
tslice = slice("2018-09-04T10:00", "2018-09-04T11:00")
ds_ = dse.sel(time=tslice)

uvkwargs = dict(vmin=-0.2, vmax=0.2, cmap="RdBu_r", add_colorbar=False)
echokwargs = dict(add_colorbar=False, vmin=0, vmax=232)
corkwargs = dict(add_colorbar=False, vmin=0, vmax=143)

fig, axs = plt.subplots(15, 1, sharey=False, sharex=True, figsize=(20, 50))
ds_.u.plot(ax=axs[0], **uvkwargs)
ds_.a1.plot(ax=axs[1], **echokwargs)
ds_.q1.plot(ax=axs[2], **corkwargs)
ds_.v.plot(ax=axs[3], **uvkwargs)
ds_.a2.plot(ax=axs[4], **echokwargs)
ds_.q2.plot(ax=axs[5], **corkwargs)
ds_.w.plot(ax=axs[6], **uvkwargs)
ds_.a3.plot(ax=axs[7], **echokwargs)
ds_.q3.plot(ax=axs[8], **corkwargs)
ds_.err.plot(ax=axs[9], **uvkwargs)
ds_.a4.plot(ax=axs[10], **echokwargs)
ds_.q4.plot(ax=axs[11], **corkwargs)
ds_.heading.plot(ax=axs[12], marker=".", linestyle="")
ds_.pitch.plot(ax=axs[13], marker=".", linestyle="")
ds_.rol.plot(ax=axs[14], marker=".", linestyle="")

# %% [markdown]
# # Compare bin mapped velocities

# %%
var = "v4"

uvkwargs = dict(vmin=-0.1, vmax=0.1, cmap="RdBu_r")

fig, axs = plt.subplots(3, 1, sharex=True, figsize=(20, 12))
dsb[var].sel(time=tslice).plot(ax=axs[0], **uvkwargs)
axs[0].set_title("Raw beam velocity")
dsm[var].sel(time=tslice).plot(ax=axs[1], **uvkwargs)
axs[1].set_title("Mapped beam velocity")
(dsb[var].sel(time=tslice)- dsm[var].sel(time=tslice)).plot(ax=axs[2], **uvkwargs)
axs[2].set_title("Raw minus mapped velocity")

# %% [markdown]
# # Zoom in on mapped coordinates

# %%
tslice = slice("2018-09-04T10:30", "2018-09-04T10:50")
qc = 64.

ds_ = dsm.sel(time=tslice)

uvkwargs = dict(vmin=-0.2, vmax=0.2, cmap="RdBu_r", add_colorbar=False)
echokwargs = dict(add_colorbar=False, vmin=0, vmax=232)
corkwargs = dict(add_colorbar=False, vmin=0, vmax=143)

fig, axs = plt.subplots(15, 1, sharey=False, sharex=True, figsize=(20, 50))
ds_.v1.where(ds_.q1 > qc).plot(ax=axs[0], **uvkwargs)
ds_.a1.plot(ax=axs[1], **echokwargs)
ds_.q1.plot(ax=axs[2], **corkwargs)
ds_.v2.where(ds_.q2 > qc).plot(ax=axs[3], **uvkwargs)
ds_.a2.plot(ax=axs[4], **echokwargs)
ds_.q2.plot(ax=axs[5], **corkwargs)
ds_.v3.where(ds_.q3 > qc).plot(ax=axs[6], **uvkwargs)
ds_.a3.plot(ax=axs[7], **echokwargs)
ds_.q3.plot(ax=axs[8], **corkwargs)
ds_.v4.where(ds_.q4 > qc).plot(ax=axs[9], **uvkwargs)
ds_.a4.plot(ax=axs[10], **echokwargs)
ds_.q4.plot(ax=axs[11], **corkwargs)
ds_.heading.plot(ax=axs[12], marker=".", linestyle="")
ds_.pitch.plot(ax=axs[13], marker=".", linestyle="")
ds_.rol.plot(ax=axs[14], marker=".", linestyle="")

# %% [markdown]
# Compare echo intensity in different beams.

# %%
time0 = "2018-09-04T10:47:01.51"

ds_ = dsb.sel(time=time0, method="nearest")

fig, ax = plt.subplots()
ds_.a1.plot(ax=ax, label="1")
ds_.a2.plot(ax=ax, label="2")
ds_.a3.plot(ax=ax, label="3")
ds_.a4.plot(ax=ax, label="4")
ax.legend()
ax.set_ylim(70, 200)

ds_ = dsm.sel(time=time0, method="nearest")

fig, ax = plt.subplots()
ds_.a1.plot(ax=ax, label="1")
ds_.a2.plot(ax=ax, label="2")
ds_.a3.plot(ax=ax, label="3")
ds_.a4.plot(ax=ax, label="4")
ax.legend()
ax.set_ylim(70, 200)

# %% [markdown]
# Where are the beams pointing?
#
# We need to convert from instrument pitch/roll/heading to spherical polar coordinates and then into cartesian coordinates for plotting.

# %%
inclination = 180 - np.rad2deg(np.arccos(np.cos(np.deg2rad(ds_.pitch))*np.cos(np.deg2rad(ds_.rol))))
azimuth = 90 - ds_.heading
r = ds_.distance/np.cos(np.deg2rad(ds_.beamAngle))

print(f"inclination = {inclination.data}")
print(f"azimuth = {azimuth.data}")
print(f"r[0] = {r[0].data}, r[-1] = {r[-1].data}")

# %% [markdown]
# For each beam...

# %%
inclination1 = np.rad2deg(np.arccos(np.cos(np.deg2rad(ds_.pitch))*np.cos(np.deg2rad(ds_.rol + ds_.beamAngle)))).data
inclination2 = np.rad2deg(np.arccos(np.cos(np.deg2rad(ds_.pitch))*np.cos(np.deg2rad(ds_.rol - ds_.beamAngle)))).data
inclination3 = np.rad2deg(np.arccos(np.cos(np.deg2rad(ds_.pitch - ds_.beamAngle))*np.cos(np.deg2rad(ds_.rol)))).data
inclination4 = np.rad2deg(np.arccos(np.cos(np.deg2rad(ds_.pitch + ds_.beamAngle))*np.cos(np.deg2rad(ds_.rol)))).data
azimuth1 = (90 - ds_.heading - 90).data
azimuth2 = (90 - ds_.heading + 90).data
azimuth3 = (90 - ds_.heading).data
azimuth4 = (90 - ds_.heading + 180).data
r = (ds_.distance/np.cos(np.deg2rad(ds_.beamAngle))).data

# %%
import utm

# %%
x0, y0, _, _ = utm.from_latlon(dsm.lat.data, dsm.lon.data)
z0 = -ds_.p.mean().data

# %%
xyz1 = np.array([r*np.cos(np.deg2rad(azimuth1))*np.sin(np.deg2rad(inclination1)) + x0, r*np.sin(np.deg2rad(azimuth1))*np.sin(np.deg2rad(inclination1)) + y0, r*np.cos(np.deg2rad(inclination1)) + z0])
xyz2 = np.array([r*np.cos(np.deg2rad(azimuth2))*np.sin(np.deg2rad(inclination2)) + x0, r*np.sin(np.deg2rad(azimuth2))*np.sin(np.deg2rad(inclination2)) + y0, r*np.cos(np.deg2rad(inclination2)) + z0])
xyz3 = np.array([r*np.cos(np.deg2rad(azimuth3))*np.sin(np.deg2rad(inclination3)) + x0, r*np.sin(np.deg2rad(azimuth3))*np.sin(np.deg2rad(inclination3)) + y0, r*np.cos(np.deg2rad(inclination3)) + z0])
xyz4 = np.array([r*np.cos(np.deg2rad(azimuth4))*np.sin(np.deg2rad(inclination4)) + x0, r*np.sin(np.deg2rad(azimuth4))*np.sin(np.deg2rad(inclination4)) + y0, r*np.cos(np.deg2rad(inclination4)) + z0])

# %%
from mpl_toolkits import mplot3d

# %%
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

ax.plot3D(xyz1[0, :], xyz1[1, :], xyz1[2, :], 'gray')
ax.plot3D(xyz2[0, :], xyz2[1, :], xyz2[2, :], 'gray')
ax.plot3D(xyz3[0, :], xyz3[1, :], xyz3[2, :], 'gray')
ax.plot3D(xyz4[0, :], xyz4[1, :], xyz4[2, :], 'gray')

ax.text(xyz1[0, -1], xyz1[1, -1], xyz1[2, -1], "1")
ax.text(xyz2[0, -1], xyz2[1, -1], xyz2[2, -1], "2")
ax.text(xyz3[0, -1], xyz3[1, -1], xyz3[2, -1], "3")
ax.text(xyz4[0, -1], xyz4[1, -1], xyz4[2, -1], "4")

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")


# %%
dlon = 0.01
dlat = 0.005

bathy = xr.open_dataset("../proc/bathy_sep_2018.nc")
bathy_ = bathy.sel(lon=slice(dsm.lon.data - dlon, dsm.lon.data + dlon/4), lat=slice(dsm.lat.data - dlat, dsm.lat.data + dlat))
                                                    

# %%
fig, ax = plt.subplots(1, 1, figsize=(20, 20))
ax.pcolormesh(bathy_.x, bathy_.y, -bathy_.H)
ax.plot(xyz1[0, :], xyz1[1, :], 'b')
ax.plot(xyz2[0, :], xyz2[1, :], 'b')
ax.plot(xyz3[0, :], xyz3[1, :], 'b')
ax.plot(xyz4[0, :], xyz4[1, :], 'b')

ax.annotate("1", (xyz1[0, -1], xyz1[1, -1]), color="b")
ax.annotate("2", (xyz2[0, -1], xyz2[1, -1]), color="b")
ax.annotate("3", (xyz3[0, -1], xyz3[1, -1]), color="b")
ax.annotate("4", (xyz4[0, -1], xyz4[1, -1]), color="b")

# %%
dlon = 0.005
dlat = 0.0025

bathy_ = bathy.sel(lon=slice(dsm.lon.data - dlon, dsm.lon.data + dlon/2), lat=slice(dsm.lat.data - dlat, dsm.lat.data + dlat))

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

ax.plot3D(xyz1[0, :], xyz1[1, :], xyz1[2, :], 'gray')
ax.plot3D(xyz2[0, :], xyz2[1, :], xyz2[2, :], 'gray')
ax.plot3D(xyz3[0, :], xyz3[1, :], xyz3[2, :], 'gray')
ax.plot3D(xyz4[0, :], xyz4[1, :], xyz4[2, :], 'gray')

ax.text(xyz1[0, -1], xyz1[1, -1], xyz1[2, -1], "1")
ax.text(xyz2[0, -1], xyz2[1, -1], xyz2[2, -1], "2")
ax.text(xyz3[0, -1], xyz3[1, -1], xyz3[2, -1], "3")
ax.text(xyz4[0, -1], xyz4[1, -1], xyz4[2, -1], "4")

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

ax.plot_surface(bathy_.x, bathy_.y, -bathy_.H)

ax.view_init(40, 200)

# %%
