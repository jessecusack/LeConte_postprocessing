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
import pycurrents as pyc
from pycurrents import adcp
import matplotlib.pyplot as plt

# %%
help(adcp)

# %%
from pycurrents.adcp import rdiraw

# %%
m = rdiraw.Multiread("/Users/jmcusack/Dropbox/LeConte/Data/ocean/september2018/raw/moorings/ABLE_Deep/16670013.000", "wh")

# %%
dat = m.read(start=50000, stop=100000)

# %%
type(dat)

# %%
from pycurrents.adcp import adcp_nc

# %%
adcp_nc.make_nc_long

# %%
sV = rdiraw.Multiread("/Users/jmcusack/Dropbox/LeConte/Data/ocean/september2018/raw/moorings/ABLE_Sentinel/ADCP/LeConte S-ABLE Sept2018 20180901T233454.pd0", "sv")

# %%
datsV = sV.read(start=50000, stop=100000)

# %%
plt.pcolormesh(datsV.cor1)

# %%
plt.pcolormesh(datsV.vbvel[20000:24000, :], cmap="RdBu_r", vmin=-0.1, vmax=0.1)

# %%
