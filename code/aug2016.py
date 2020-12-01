import os

import gsw
import mat73
import mixsea as mx
import numpy as np
import xarray as xr
from munch import munchify
from tqdm import tqdm

import ADCP
import clargs
import CTD
import utils

############### PRELIMINARIES ################
# Parse command line arguments
args = clargs.parse_check_args()

sdroot = args.save

files = utils.find_files(args, "aug2016")

params = utils.load_parameters()

############### LOAD DATA ################
print("Loading data. (---> This may trigger Dropbox download <---)")
print("Loading VMP.")
vmp = munchify(utils.loadmat(files.vmp, check_arrays=True, mat_dtype=True))
print("Loading CTD.")
ctds = utils.loadmat(files.ctd, check_arrays=True, mat_dtype=True)["ctd"]
ctds = [munchify(ctd) for ctd in ctds]
print("Loading SADCP. (Warning: large file ~ 1 GB.)")
sadcp = munchify(mat73.loadmat(files.sadcp)["adcp"])

############### VMP ################
print("Processing VMP")
ctd = ctds[6]  # This is the ctd data corresponding to the VMP.
bin_width = params.ctd.adiabatic_level_bin_width
gam = params.vmp.mixing_efficiency

ctd = ctds[6]

use = np.isfinite(ctd.time)
time = ctd.time[use]
C = ctd.C[:, use]
SP = ctd.S[:, use]
t = ctd.T[:, use]
lon = ctd.lon[use]
lat = ctd.lat[use]
depth = ctd.depth

eps = np.full_like(t, np.nan)

p = gsw.p_from_z(-depth, np.mean(lat))

for ic in range(time.size):

    iv = utils.closest_index(time[ic], vmp.Tm)

    pv = vmp.Pm[:, iv]
    ev = vmp.Em[:, iv]

    dp = np.nanmean(np.diff(pv))

    pmax = np.nanmax(pv)
    pmin = np.nanmin(pv)

    use = np.isfinite(pv) & np.isfinite(ev)

    # Do the interpolation in log space because eps varies over
    # so many orders of magnitude.
    eps[:, ic] = 10 ** np.interp(p, pv[use], np.log10(ev[use]))
    # Scrub out of range data
    eps[(p > pmax) | (p < pmin), ic] = np.nan


# Sea water thermodynamics
p, SA, CT, sig0, p_mid, N2 = CTD.common_thermodynamics(depth, lon, lat, SP, t)

N2_ref = np.full_like(t, np.nan)
print("Adiabatic levelling of buoyancy frequency.")
for i in tqdm(range(time.size)):
    N2_ref[:, i] = mx.nsq.adiabatic_leveling(
        p, SP[:, i], t[:, i], lon[i], lat[i], bin_width=bin_width
    )

# Max valid depth
depth_max = np.full_like(lon, np.nan)
for i in range(time.size):
    valid_ctd = np.isfinite(t[:, i])
    depth_max[i] = depth[valid_ctd].max()

# Ozmidov scale
Lo = np.sqrt(eps / N2_ref ** (3 / 2))

# Turbulent diffusivity using the Osborne relation
Kv = gam * eps / N2_ref

vmp_datavars = {
    "C": (["depth", "profile"], C, {"Variable": "Conductivity"}),
    "SP": (["depth", "profile"], SP, {"Variable": "Practical salinity"}),
    "t": (["depth", "profile"], t, {"Variable": "Temperature (in situ)"}),
    "CT": (["depth", "profile"], CT, {"Variable": "Conservative temperature"}),
    "SA": (["depth", "profile"], SA, {"Variable": "Absolute salinity"}),
    "eps": (
        ["depth", "profile"],
        eps,
        {"Variable": "Turbulent dissipation rate of kinetic energy"},
    ),
    "Kv": (
        ["depth", "profile"],
        Kv,
        {"Variable": "Turbulent diffusivity"},
    ),
    "sig0": (
        ["depth", "profile"],
        sig0,
        {"Variable": "Potential density referenced to 0 dbar"},
    ),
    "N2": (["depth_mid", "profile"], N2, {"Variable": "Buoyancy frequency"}),
    "N2_ref": (
        ["depth", "profile"],
        N2_ref,
        {
            "Variable": "Adiabatically leveled buoyancy frequency, using {:1.0f} dbar bin".format(
                bin_width
            )
        },
    ),
    "Lo": (["depth", "profile"], Lo, {"Variable": "Ozmidov length scale"}),
    "depth_max": (["profile"], depth_max, {"Variable": "Maximum valid data depth"}),
}

vmp_coords = {
    "profile": (["profile"], np.arange(time.size) + 1),
    "depth": (["depth"], depth),
    "depth_mid": (["depth_mid"], utils.mid(depth)),
    "time": (["profile"], utils.datenum_to_datetime(time)),
    "lon": (["profile"], lon),
    "lat": (["profile"], lat),
    "p": (["depth"], p),
    "p_mid": (["depth_mid"], p_mid[:, 0]),
}

vmp_ds = xr.Dataset(vmp_datavars, vmp_coords)
file = "vmp_aug_2016.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
vmp_ds.to_netcdf(os.path.join(sdroot, file))


############### COMBINED ################
print("Processing combined VMP-SADCP")
time_win = params.sadcp.time_window
rmax = params.sadcp.rmax
vmax = params.sadcp.vmax
range_min = params.sadcp.range_min

print("Interpolating velocity to VMP stations.")
mask = np.isfinite(t)
u, v, w, lon, lat, range_bottom, nav = ADCP.interp_ADCP_2D(
    sadcp,
    mask,
    depth,
    lon,
    lat,
    time,
    time_win=time_win,
    rmax=rmax,
    vmax=vmax,
    range_min=range_min,
)


sadcp_datavars = {
    "u": (["depth", "profile"], u, {"Variable": "Northward velocity"}),
    "v": (["depth", "profile"], v, {"Variable": "Eastward velocity"}),
    "w": (["depth", "profile"], u, {"Variable": "Vertical velocity"}),
    "range_bottom": (["profile"], range_bottom, {"Variable": "Range to bottom"}),
    "nav": (["profile"], nav, {"Variable": "Number of ADCP profiles in average."}),
}

sadcp_coords = {
    "lon_sadcp": (["profile"], lon, {"Variable": "Mean longitude of SADCP data"}),
    "lat_sadcp": (["profile"], lat, {"Variable": "Mean latitude of SADCP data"}),
}


combo_ds = xr.Dataset(
    {**vmp_datavars, **sadcp_datavars}, {**vmp_coords, **sadcp_coords}
)

file = "combo_aug_2016.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
combo_ds.to_netcdf(os.path.join(sdroot, file))
