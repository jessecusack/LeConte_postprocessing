import os

import gsw
import numpy as np
import xarray as xr
from munch import munchify
from tqdm import tqdm

import clargs
import utils

############### PRELIMINARIES ################
# Parse command line arguments
args = clargs.parse_check_args()

sdroot = args.save

files = utils.find_files(args, "sep2018moorings")

params = utils.load_parameters()

print("Loading data may trigger Dropbox download!")

############### COMBINE ################
print("Combining ADCP and CTD datasets")
print("Loading ADCP data")
adcp = munchify(utils.loadmat(files.deep_adcp, check_arrays=True, mat_dtype=True)["V"])
print("Loading CTD data")
ctd = munchify(
    utils.loadmat(files.deep_ctd, check_arrays=True, mat_dtype=True)["mooring"]
)
print("Loading coordinates")
coords = munchify(
    utils.loadmat(files.coords, check_arrays=True, mat_dtype=True)["mooringCoord"]
)

rho0 = 1025.
g = 9.81

ctd.p = np.flipud(ctd.pop("P"))
ctd.lon = coords.moorD[1]
ctd.lat = coords.moorD[0]
ctd.depth_nominal = np.flipud(ctd.z)  # Note sure about this.
ctd.SP = np.flipud(ctd.pop("S"))
ctd.t = np.flipud(ctd.pop("T"))

print("Despiking salinity")
n1 = 2
n2 = 5
size = 201
xmin = 0.0
xmax = 100.0
fill = True
for i in tqdm(np.argwhere(np.isfinite(ctd.p[:, 0])).squeeze()):
    ctd.SP[i, :] = utils.despike(ctd.SP[i, :], n1=n1, n2=n2, size=size, xmin=xmin, xmax=xmax, fill=fill)
    
# Interpolate salinity and pressure to level of thermistors
print("Interpolating CTD")
for i in tqdm(range(ctd.time.size)):
    nans = np.isnan(ctd.p[:, i])
    ctd.p[nans, i] = np.interp(
        ctd.depth_nominal[nans], ctd.depth_nominal[~nans], ctd.p[~nans, i]
    )
    ctd.SP[nans, i] = np.interp(ctd.p[nans, i], ctd.p[~nans, i], ctd.SP[~nans, i])

# Thermodynamics
ctd.depth = -gsw.z_from_p(ctd.p, ctd.lat)
ctd.SA = gsw.SA_from_SP(ctd.SP, ctd.p, ctd.lon, ctd.lat)
ctd.CT = gsw.CT_from_t(ctd.SA, ctd.t, ctd.p)
ctd.sig0 = gsw.pot_rho_t_exact(ctd.SA, ctd.t, ctd.p, 0)
ctd.N2, ctd.p_mid = gsw.Nsquared(ctd.SA, ctd.CT, ctd.p, ctd.lat)
ctd.b = -g*(ctd.sig0 - rho0)/rho0
ctd = utils.apply_utm(ctd)

# ADCP
print("Despiking ADCP")
adcp.depth = adcp.pop("z")
n1 = 2
n2 = 5
size = 201
xmin = -1.0
xmax = 1.0
fill = True
for i in tqdm(range(adcp.depth.size)):
    adcp.u[i, :] = utils.despike(
        adcp.u[i, :], n1=n1, n2=n2, size=size, xmin=xmin, xmax=xmax, fill=fill
    )
    adcp.v[i, :] = utils.despike(
        adcp.v[i, :], n1=n1, n2=n2, size=size, xmin=xmin, xmax=xmax, fill=fill
    )
    adcp.w[i, :] = utils.despike(
        adcp.w[i, :], n1=n1, n2=n2, size=size, xmin=xmin, xmax=xmax, fill=fill
    )
    
print("Interpolating CTD to ADCP depths")
adcp.SP = utils.nan_interp(adcp.depth, ctd.depth, ctd.SP, axis=0)
adcp.t = utils.nan_interp(adcp.depth, ctd.depth, ctd.t, axis=0)
adcp.CT = utils.nan_interp(adcp.depth, ctd.depth, ctd.CT, axis=0)
adcp.SA = utils.nan_interp(adcp.depth, ctd.depth, ctd.SA, axis=0)
adcp.sig0 = utils.nan_interp(adcp.depth, ctd.depth, ctd.sig0, axis=0)
adcp.b = -g*(adcp.sig0 - rho0)/rho0

datavars = {
    "SP": (["i", "time"], ctd.SP, {"Variable": "Practical salinity"}),
    "t": (["i", "time"], ctd.t, {"Variable": "Temperature (in situ)"}),
    "CT": (["i", "time"], ctd.CT, {"Variable": "Conservative temperature"}),
    "SA": (["i", "time"], ctd.SA, {"Variable": "Absolute salinity"}),
    "b": (["i", "time"], ctd.b, {"Variable": "Buoyancy"}),
    "sig0": (
        ["i", "time"],
        ctd.sig0,
        {"Variable": "Potential density referenced to 0 dbar"},
    ),
    "N2": (["i_mid", "time"], ctd.N2, {"Variable": "Buoyancy frequency"}),
    "u": (["depth_adcp", "time"], adcp.u, {"Variable": "Eastward velocity"}),
    "v": (["depth_adcp", "time"], adcp.v, {"Variable": "Northward velocity"}),
    "w": (["depth_adcp", "time"], adcp.w, {"Variable": "Vertical velocity"}),
    "SPa": (["depth_adcp", "time"], adcp.SP, {"Variable": "Practical salinity"}),
    "ta": (["depth_adcp", "time"], adcp.t, {"Variable": "Temperature (in situ)"}),
    "CTa": (["depth_adcp", "time"], adcp.CT, {"Variable": "Conservative temperature"}),
    "SAa": (["depth_adcp", "time"], adcp.SA, {"Variable": "Absolute salinity"}),
    "ba": (["depth_adcp", "time"], adcp.b, {"Variable": "Buoyancy"}),
    "sig0a": (
        ["depth_adcp", "time"],
        adcp.sig0,
        {"Variable": "Potential density referenced to 0 dbar"},
    ),
}

coords = {
    "time": (["time"], utils.datenum_to_datetime(ctd.time)),
    "depth": (["i", "time"], ctd.depth),
    "depth_adcp": (["depth_adcp"], adcp.depth),
    "depth_nominal": (["i"], ctd.depth_nominal),
    "lon": ([], ctd.lon),
    "lat": ([], ctd.lat),
    "x": ([], ctd.x),
    "y": ([], ctd.y),
    "zone_letter": ([], ctd.zone_letter),
    "zone_number": ([], ctd.zone_number),
    "p": (["i", "time"], ctd.p),
    "p_mid": (["i_mid", "time"], ctd.p_mid),
}

ds = xr.Dataset(datavars, coords)
file = "deep_mooring_sep_2018.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
ds.to_netcdf(os.path.join(sdroot, file))
