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
import VMP

############### PRELIMINARIES ################
# Parse command line arguments
args = clargs.parse_check_args()

sdroot = args.save

files = utils.find_files(args, "aug2016")

params = utils.load_parameters()

print("Loading data may trigger Dropbox download!")

############### CTD ################
print("Loading CTD")
ctds = utils.loadmat(files.ctd, check_arrays=True, mat_dtype=True)["ctd"]
ctds = [munchify(ctd) for ctd in ctds]
print("Processing CTD")
bin_width = params.ctd.adiabatic_level_bin_width
depth_min = params.ctd.depth_min
depth_max = params.ctd.depth_max
depth_spacing = params.ctd.depth_spacing

ctd = CTD.generate_CTD_Munch_from_list(ctds, depth_min=depth_min, depth_max=depth_max, depth_spacing=depth_spacing)

ctd = CTD.apply_thermodynamics(ctd)

ctd = CTD.apply_adiabatic_level(ctd, bin_width)

ctd.depth_max = CTD.depth_max(ctd.depth, np.isfinite(ctd.t))

ctd_datavars = {
    "SP": (["depth", "profile"], ctd.SP, {"Variable": "Practical salinity"}),
    "t": (["depth", "profile"], ctd.t, {"Variable": "Temperature (in situ)"}),
    "CT": (["depth", "profile"], ctd.CT, {"Variable": "Conservative temperature"}),
    "SA": (["depth", "profile"], ctd.SA, {"Variable": "Absolute salinity"}),
    "sig0": (
        ["depth", "profile"],
        ctd.sig0,
        {"Variable": "Potential density referenced to 0 dbar"},
    ),
    "N2": (["depth_mid", "profile"], ctd.N2, {"Variable": "Buoyancy frequency"}),
    "N2_ref": (
        ["depth", "profile"],
        ctd.N2_ref,
        {
            "Variable": "Adiabatically leveled buoyancy frequency, using {:1.0f} dbar bin".format(
                bin_width
            )
        },
    ),
    "depth_max": (["profile"], ctd.depth_max, {"Variable": "Maximum valid data depth"}),
}

ctd_coords = {
    "profile": (["profile"], np.arange(ctd.time.size) + 1),
    "depth": (["depth"], ctd.depth),
    "depth_mid": (["depth_mid"], utils.mid(ctd.depth)),
    "time": (["profile"], utils.datenum_to_datetime(ctd.time)),
    "lon": (["profile"], ctd.lon),
    "lat": (["profile"], ctd.lat),
    "p": (["depth"], ctd.p),
    "p_mid": (["depth_mid"], ctd.p_mid),
}

ctd_ds = xr.Dataset(ctd_datavars, ctd_coords)
file = "ctd_aug_2016.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
ctd_ds.to_netcdf(os.path.join(sdroot, file))

############### VMP ################
print("Loading VMP")
vmp = munchify(utils.loadmat(files.vmp, check_arrays=True, mat_dtype=True))
print("Processing VMP")
gam = params.vmp.mixing_efficiency
depth_min = params.vmp.depth_min
depth_max = params.vmp.depth_max
depth_spacing = params.vmp.depth_spacing
time_win = params.vmp.time_window

depth = -gsw.z_from_p(vmp.Pm, vmp.lat[np.newaxis, :])
vmp = VMP.generate_VMP_Munch(vmp.Tm, depth, vmp.lon, vmp.lat, vmp.Em, eps2=None, depth_min=depth_min, depth_max=depth_max, depth_spacing=depth_spacing)

vmp = VMP.regrid_ctd_to_vmp(ctd, vmp, time_win)

vmp.Lo1, vmp.Kv1 = VMP.common_turbulence(vmp.eps1, vmp.N2_ref, gam)

vmp_datavars = {
    "eps1": (
        ["depth", "profile"], 
        vmp.eps1,
        {"Variable": "Turbulent dissipation rate of kinetic energy"},
    ),
    "Kv1": (
        ["depth", "profile"],
        vmp.Kv1,
        {"Variable": "Turbulent diffusivity"},
    ),
    "Lo1": (["depth", "profile"], vmp.Lo1, {"Variable": "Ozmidov length scale"}),
    "SP": (["depth", "profile"], vmp.SP, {"Variable": "Practical salinity"}),
    "t": (["depth", "profile"], vmp.t, {"Variable": "Temperature (in situ)"}),
    "CT": (["depth", "profile"], vmp.CT, {"Variable": "Conservative temperature"}),
    "SA": (["depth", "profile"], vmp.SA, {"Variable": "Absolute salinity"}),
    "sig0": (
        ["depth", "profile"],
        vmp.sig0,
        {"Variable": "Potential density referenced to 0 dbar"},
    ),
    "N2": (["depth_mid", "profile"], vmp.N2, {"Variable": "Buoyancy frequency"}),
    "N2_ref": (
        ["depth", "profile"],
        vmp.N2_ref,
        {
            "Variable": "Adiabatically leveled buoyancy frequency, using {:1.0f} dbar bin".format(
                bin_width
            )
        },
    ),
    "depth_max": (["profile"], vmp.depth_max, {"Variable": "Maximum valid data depth"}),
}

vmp_coords = {
    "profile": (["profile"], np.arange(vmp.time.size) + 1),
    "depth": (["depth"], vmp.depth),
    "depth_mid": (["depth_mid"], utils.mid(vmp.depth)),
    "time": (["profile"], utils.datenum_to_datetime(vmp.time)),
    "lon": (["profile"], vmp.lon),
    "lat": (["profile"], vmp.lat),
    "p": (["depth"], vmp.p),
    "p_mid": (["depth_mid"], vmp.p_mid),
}

vmp_ds = xr.Dataset(vmp_datavars, vmp_coords)
file = "vmp_aug_2016.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
vmp_ds.to_netcdf(os.path.join(sdroot, file))

############### COMBINED ################
print("Loading SADCP. (Warning: large file ~ 1 GB.)")
sadcp = munchify(mat73.loadmat(files.sadcp)["adcp"])
print("Processing combined VMP-SADCP")
time_win = params.sadcp.time_window
rmax = params.sadcp.rmax
vmax = params.sadcp.vmax
range_min = params.sadcp.range_min

print("Interpolating velocity to VMP stations.")
mask = np.isfinite(ctd.t)
u, v, w, lon, lat, range_bottom, nav = ADCP.interp_ADCP_2D(
    sadcp,
    mask,
    ctd.depth,
    ctd.lon,
    ctd.lat,
    ctd.time,
    time_win=time_win,
    rmax=rmax,
    vmax=vmax,
    range_min=range_min,
)

ctd = CTD.regrid_vmp_to_ctd(vmp, ctd, time_win=params.vmp.time_window)


combo_datavars = {
    "eps1": (
        ["depth", "profile"], 
        ctd.eps1,
        {"Variable": "Turbulent dissipation rate of kinetic energy"},
    ),
    "Kv1": (
        ["depth", "profile"],
        ctd.Kv1,
        {"Variable": "Turbulent diffusivity"},
    ),
    "Lo1": (["depth", "profile"], ctd.Lo1, {"Variable": "Ozmidov length scale"}),
    "SP": (["depth", "profile"], ctd.SP, {"Variable": "Practical salinity"}),
    "t": (["depth", "profile"], ctd.t, {"Variable": "Temperature (in situ)"}),
    "CT": (["depth", "profile"], ctd.CT, {"Variable": "Conservative temperature"}),
    "SA": (["depth", "profile"], ctd.SA, {"Variable": "Absolute salinity"}),
    "sig0": (
        ["depth", "profile"],
        ctd.sig0,
        {"Variable": "Potential density referenced to 0 dbar"},
    ),
    "N2": (["depth_mid", "profile"], ctd.N2, {"Variable": "Buoyancy frequency"}),
    "N2_ref": (
        ["depth", "profile"],
        ctd.N2_ref,
        {
            "Variable": "Adiabatically leveled buoyancy frequency, using {:1.0f} dbar bin".format(
                bin_width
            )
        },
    ),
    "depth_max": (["profile"], ctd.depth_max, {"Variable": "Maximum valid data depth"}),
    "u": (["depth", "profile"], u, {"Variable": "Northward velocity"}),
    "v": (["depth", "profile"], v, {"Variable": "Eastward velocity"}),
    "w": (["depth", "profile"], u, {"Variable": "Vertical velocity"}),
    "range_bottom": (["profile"], range_bottom, {"Variable": "Range to bottom"}),
    "nav": (["profile"], nav, {"Variable": "Number of ADCP profiles in average."}),
}

combo_coords = {
    "profile": (["profile"], np.arange(ctd.time.size) + 1),
    "depth": (["depth"], ctd.depth),
    "depth_mid": (["depth_mid"], utils.mid(ctd.depth)),
    "time": (["profile"], utils.datenum_to_datetime(ctd.time)),
    "lon": (["profile"], ctd.lon),
    "lat": (["profile"], ctd.lat),
    "p": (["depth"], ctd.p),
    "p_mid": (["depth_mid"], ctd.p_mid),
    "lon_sadcp": (["profile"], lon, {"Variable": "Mean longitude of SADCP data"}),
    "lat_sadcp": (["profile"], lat, {"Variable": "Mean latitude of SADCP data"}),
}

combo_ds = xr.Dataset(combo_datavars, combo_coords)

file = "combo_aug_2016.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
combo_ds.to_netcdf(os.path.join(sdroot, file))
