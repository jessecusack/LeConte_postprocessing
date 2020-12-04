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

files = utils.find_files(args, "sep2018")

params = utils.load_parameters()

print("Loading data may trigger Dropbox download!")

############### LOAD DATA ################
print("Loading sections.")
secs = utils.loadmat(files.vmp_sections, check_arrays=True, mat_dtype=True)["vmp"]
secs = np.array([munchify(ds) for ds in secs])
print("Loading VMP.")
ds = utils.loadmat(files.vmp, check_arrays=True, mat_dtype=True)
vmp = munchify(ds["vmp"])
gps = munchify(ds["gps"])
print("Loading SADCP.")
sadcp = munchify(mat73.loadmat(files.sadcp)["adcp"])

############### CTD ################
print("Processing CTD")
bin_width = params.ctd.adiabatic_level_bin_width
depth_min = params.ctd.depth_min
depth_max = params.ctd.depth_max
depth_spacing = params.ctd.depth_spacing

ctd = CTD.generate_CTD_Munch(
    vmp.time,
    vmp.depth,
    vmp.lon,
    vmp.lat,
    vmp.salt,
    vmp.JAC_T,
    depth_min=depth_min,
    depth_max=depth_max,
    depth_spacing=depth_spacing,
)

ctd = CTD.apply_thermodynamics(ctd)

ctd = CTD.apply_adiabatic_level(ctd, bin_width)

ctd.depth_max = CTD.depth_max(ctd.depth, np.isfinite(ctd.t))

ctd = utils.apply_utm(ctd)

insection = np.full((len(secs), vmp.time.size), False)
section = np.arange(len(secs)) + 1
for i, sec in enumerate(secs):
    insection[i, sec.profiles.index.astype(int) - 1] = True

time_start = np.hstack([utils.datenum_to_datetime(sec.starttime) for sec in secs])
time_end = np.hstack([utils.datenum_to_datetime(sec.endtime) for sec in secs])
lon_start = np.hstack([sec.startlon for sec in secs])
lon_end = np.hstack([sec.endlon for sec in secs])
lat_start = np.hstack([sec.startlat for sec in secs])
lat_end = np.hstack([sec.endlat for sec in secs])
section_description = np.hstack([sec.description for sec in secs])

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
    "x": (["profile"], ctd.x),
    "y": (["profile"], ctd.y),
    "zone_letter": ([], ctd.zone_letter),
    "zone_number": ([], ctd.zone_number),
    "cast": (["profile"], vmp.cast.astype(int)),
    "p": (["depth"], ctd.p),
    "p_mid": (["depth_mid"], ctd.p_mid),
    "section": (["section"], section, {"Variable": "Section number"}),
    "insection": (["section", "profile"], insection, {"Variable": "Section masks"}),
    "time_start": (["section"], time_start, {"Variable": "Section start time"}),
    "time_end": (["section"], time_end, {"Variable": "Section end time"}),
    "lon_start": (["section"], lon_start, {"Variable": "Section start longitude"}),
    "lon_end": (["section"], lon_end, {"Variable": "Section end longitude"}),
    "lat_start": (["section"], lat_start, {"Variable": "Section start latitude"}),
    "lat_end": (["section"], lat_end, {"Variable": "Section end latitude"}),
}

ctd_ds = xr.Dataset(ctd_datavars, ctd_coords)
file = "ctd_sep_2018.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
ctd_ds.to_netcdf(os.path.join(sdroot, file))

############### VMP ################
print("Processing VMP")
gam = params.vmp.mixing_efficiency
depth_min = params.vmp.depth_min
depth_max = params.vmp.depth_max
depth_spacing = params.vmp.depth_spacing
time_win = params.vmp.time_window

vmp = VMP.generate_VMP_Munch(
    vmp.time,
    vmp.depth,
    vmp.lon,
    vmp.lat,
    vmp.eps1,
    eps2=vmp.eps2,
    depth_min=depth_min,
    depth_max=depth_max,
    depth_spacing=depth_spacing,
)

vmp = VMP.regrid_ctd_to_vmp(ctd, vmp, time_win)

vmp = utils.apply_utm(vmp)

# Ozmidov scale and diffusivity
vmp.Lo1, vmp.Kv1 = VMP.common_turbulence(vmp.eps1, vmp.N2_ref, gam)
vmp.Lo2, vmp.Kv2 = VMP.common_turbulence(vmp.eps2, vmp.N2_ref, gam)

vmp_datavars = {
    "eps1": (
        ["depth", "profile"],
        vmp.eps1,
        {"Variable": "Turbulent dissipation rate of kinetic energy"},
    ),
    "eps2": (
        ["depth", "profile"],
        vmp.eps2,
        {"Variable": "Turbulent dissipation rate of kinetic energy"},
    ),
    "Kv1": (
        ["depth", "profile"],
        vmp.Kv1,
        {
            "Variable": "Turbulent diffusivity using adiabatically levelled N and assuming mixing efficiency of 0.2"
        },
    ),
    "Kv2": (
        ["depth", "profile"],
        vmp.Kv2,
        {
            "Variable": "Turbulent diffusivity using adiabatically levelled N and assuming mixing efficiency of 0.2"
        },
    ),
    "Lo1": (["depth", "profile"], vmp.Lo1, {"Variable": "Ozmidov length scale"}),
    "Lo2": (["depth", "profile"], vmp.Lo2, {"Variable": "Ozmidov length scale"}),
}

vmp_coords = {
    "profile": (["profile"], np.arange(vmp.time.size) + 1),
    "depth": (["depth"], vmp.depth),
    "depth_mid": (["depth_mid"], utils.mid(vmp.depth)),
    "time": (["profile"], utils.datenum_to_datetime(vmp.time)),
    "lon": (["profile"], vmp.lon),
    "lat": (["profile"], vmp.lat),
    "x": (["profile"], vmp.x),
    "y": (["profile"], vmp.y),
    "zone_letter": ([], vmp.zone_letter),
    "zone_number": ([], vmp.zone_number),
    "p": (["depth"], vmp.p),
    "p_mid": (["depth_mid"], vmp.p_mid),
}

vmp_ds = xr.Dataset(vmp_datavars, vmp_coords)
file = "vmp_sep_2018.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
vmp_ds.to_netcdf(os.path.join(sdroot, file))

############### COMBINED ################
print("Processing combined CTD-VMP-SADCP")
time_win = params.sadcp.time_window
rmax = params.sadcp.rmax
vmax = params.sadcp.vmax
range_min = params.sadcp.range_min

print("Interpolating velocity to CTD stations.")
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
    "eps2": (
        ["depth", "profile"],
        ctd.eps2,
        {"Variable": "Turbulent dissipation rate of kinetic energy"},
    ),
    "Kv1": (
        ["depth", "profile"],
        ctd.Kv1,
        {
            "Variable": "Turbulent diffusivity using adiabatically levelled N and assuming mixing efficiency of 0.2"
        },
    ),
    "Kv2": (
        ["depth", "profile"],
        ctd.Kv2,
        {
            "Variable": "Turbulent diffusivity using adiabatically levelled N and assuming mixing efficiency of 0.2"
        },
    ),
    "Lo1": (["depth", "profile"], ctd.Lo1, {"Variable": "Ozmidov length scale"}),
    "Lo2": (["depth", "profile"], ctd.Lo2, {"Variable": "Ozmidov length scale"}),
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
    "x": (["profile"], ctd.x),
    "y": (["profile"], ctd.y),
    "zone_letter": ([], ctd.zone_letter),
    "zone_number": ([], ctd.zone_number),
    "p": (["depth"], ctd.p),
    "p_mid": (["depth_mid"], ctd.p_mid),
    "lon_sadcp": (["profile"], lon, {"Variable": "Mean longitude of SADCP data"}),
    "lat_sadcp": (["profile"], lat, {"Variable": "Mean latitude of SADCP data"}),
    "section": (["section"], section, {"Variable": "Section number"}),
    "insection": (["section", "profile"], insection, {"Variable": "Section masks"}),
    "time_start": (["section"], time_start, {"Variable": "Section start time"}),
    "time_end": (["section"], time_end, {"Variable": "Section end time"}),
    "lon_start": (["section"], lon_start, {"Variable": "Section start longitude"}),
    "lon_end": (["section"], lon_end, {"Variable": "Section end longitude"}),
    "lat_start": (["section"], lat_start, {"Variable": "Section start latitude"}),
    "lat_end": (["section"], lat_end, {"Variable": "Section end latitude"}),
}


combo_ds = xr.Dataset(combo_datavars, combo_coords)

file = "combo_sep_2018.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
combo_ds.to_netcdf(os.path.join(sdroot, file))
