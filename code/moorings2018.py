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

ctd.p = np.flipud(ctd.pop("P"))
ctd.lon = coords.moorD[1]
ctd.lat = coords.moorD[0]
ctd.depth_nominal = np.flipud(ctd.z)  # Note sure about this.
ctd.SP = np.flipud(ctd.pop("S"))
ctd.t = np.flipud(ctd.pop("T"))

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

ctd = utils.apply_utm(ctd)

# ADCP
# print("Despiking ADCP")
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

datavars = {
    "SP": (["i", "time"], ctd.SP, {"Variable": "Practical salinity"}),
    "t": (["i", "time"], ctd.t, {"Variable": "Temperature (in situ)"}),
    "CT": (["i", "time"], ctd.CT, {"Variable": "Conservative temperature"}),
    "SA": (["i", "time"], ctd.SA, {"Variable": "Absolute salinity"}),
    "sig0": (
        ["i", "time"],
        ctd.sig0,
        {"Variable": "Potential density referenced to 0 dbar"},
    ),
    "N2": (["i_mid", "time"], ctd.N2, {"Variable": "Buoyancy frequency"}),
    "u": (["depth_adcp", "time"], adcp.u, {"Variable": "Eastward velocity"}),
    "v": (["depth_adcp", "time"], adcp.v, {"Variable": "Northward velocity"}),
    "w": (["depth_adcp", "time"], adcp.w, {"Variable": "Vertical velocity"}),
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
