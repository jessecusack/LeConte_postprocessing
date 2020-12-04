import os

import gsw
import numpy as np
import xarray as xr
from munch import munchify

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

ctd.p = ctd.pop("P")
ctd.lon = coords.moorD[1]
ctd.lat = coords.moorD[0]
ctd.depth = -gsw.z_from_p(ctd.p, ctd.lat)
ctd.depth_nominal = ctd.z  # Note sure about this.
ctd.SP = ctd.pop("S")
ctd.t = ctd.pop("T")
ctd.SA = gsw.SA_from_SP(ctd.SP, ctd.p, ctd.lon, ctd.lat)
ctd.CT = gsw.CT_from_t(ctd.SA, ctd.t, ctd.p)
ctd.sig0 = gsw.pot_rho_t_exact(ctd.SA, ctd.t, ctd.p, 0)
ctd.N2, ctd.p_mid = gsw.Nsquared(ctd.SA, ctd.CT, ctd.p, ctd.lat)

ctd = utils.apply_utm(ctd)

bad_velocty = (adcp.u > 1.0) | (adcp.v > 1.0) | (adcp.w > 1.0)
adcp.u[bad_velocty] = np.nan
adcp.v[bad_velocty] = np.nan
adcp.w[bad_velocty] = np.nan

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
    #     "N2_ref": (
    #         ["i", "time"],
    #         ctd.N2_ref,
    #         {
    #             "Variable": "Adiabatically leveled buoyancy frequency, using {:1.0f} dbar bin".format(
    #                 bin_width
    #             )
    #         },
    #     ),
}

coords = {
    "time": (["time"], utils.datenum_to_datetime(ctd.time)),
    "depth": (["i", "time"], ctd.depth),
    "depth_adcp": (["depth_adcp"], adcp.z),
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
