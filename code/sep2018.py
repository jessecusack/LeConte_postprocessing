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
import utils

############### PRELIMINARIES ################
# Parse command line arguments, convert to dict with 'vars' and then
# to Munch (similar to dict but better)
args = munchify(vars(clargs.gen_parser().parse_args()))
# Check the directories to make sure they exist.
clargs.check_args(args)

sdroot = args.save

files = utils.find_files(args, "sep2018")

params = utils.load_parameters()

############### LOAD DATA ################
print("Loading data. (---> This may trigger Dropbox download <---)")
print("Loading sections.")
secs = utils.loadmat(files.vmp_sections, check_arrays=True, mat_dtype=True)["vmp"]
secs = np.array([munchify(ds) for ds in secs])
print("Loading VMP.")
ds = utils.loadmat(files.vmp, check_arrays=True, mat_dtype=True)
vmp = munchify(ds["vmp"])
gps = munchify(ds["gps"])
print("Loading SADCP.")
sadcp = munchify(mat73.loadmat(files.sadcp)["adcp"])

############### VMP ################
print("Processing VMP")
bin_width = params.ctd.adiabatic_level_bin_width
gam = params.vmp.mixing_efficiency

t = vmp.JAC_T
C = vmp.JAC_C
SP = vmp.salt
depth = vmp.depth
eps1 = vmp.eps1
eps2 = vmp.eps2
# Temperature variance not working?
# chi1 = vmp.chi1
# chi2 = vmp.chi2
lon = vmp.lon
lat = vmp.lat

p = gsw.p_from_z(-depth, np.mean(lat))
SA = gsw.SA_from_SP(SP, p[:, np.newaxis], lon[np.newaxis, :], lat[np.newaxis, :])
CT = gsw.CT_from_t(SA, t, p[:, np.newaxis])
sig0 = gsw.pot_rho_t_exact(SA, t, p[:, np.newaxis], 0)
N2, p_mid = gsw.Nsquared(SA, CT, p[:, np.newaxis], lat[np.newaxis, :])

N2_ref = np.full_like(t, np.nan)
print("Adiabatic levelling of buoyancy frequency.")
for i in tqdm(range(vmp.time.size)):
    N2_ref[:, i] = mx.nsq.adiabatic_leveling(
        p, SP[:, i], t[:, i], lon[i], lat[i], bin_width=bin_width
    )

# Max valid depth
depth_max = np.full_like(lon, np.nan)
for i in range(vmp.time.size):
    valid_vmp = np.isfinite(vmp.eps1[:, i])
    depth_max[i] = vmp.depth[valid_vmp].max()

# Ozmidov scale
Lo1 = np.sqrt(eps1 / N2_ref ** (3 / 2))
Lo2 = np.sqrt(eps1 / N2_ref ** (3 / 2))

# Turbulent diffusivity using the Osborne relation
Kv1 = gam * eps1 / N2_ref
Kv2 = gam * eps2 / N2_ref

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

vmp_datavars = {
    "C": (["depth", "profile"], vmp.JAC_C, {"Variable": "Conductivity"}),
    "SP": (["depth", "profile"], vmp.salt, {"Variable": "Practical salinity"}),
    "t": (["depth", "profile"], vmp.JAC_T, {"Variable": "Temperature (in situ)"}),
    "CT": (["depth", "profile"], CT, {"Variable": "Conservative temperature"}),
    "SA": (["depth", "profile"], SA, {"Variable": "Absolute salinity"}),
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
        Kv1,
        {
            "Variable": "Turbulent diffusivity using adiabatically levelled N and assuming mixing efficiency of 0.2"
        },
    ),
    "Kv2": (
        ["depth", "profile"],
        Kv2,
        {
            "Variable": "Turbulent diffusivity using adiabatically levelled N and assuming mixing efficiency of 0.2"
        },
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
    "Lo1": (["depth", "profile"], Lo1, {"Variable": "Ozmidov length scale"}),
    "Lo2": (["depth", "profile"], Lo2, {"Variable": "Ozmidov length scale"}),
    "depth_max": (["profile"], depth_max, {"Variable": "Maximum valid data depth"}),
}

vmp_coords = {
    "profile": (["profile"], np.arange(vmp.time.size) + 1),
    "depth": (["depth"], depth),
    "depth_mid": (["depth_mid"], utils.mid(depth)),
    "time": (["profile"], utils.datenum_to_datetime(vmp.time)),
    "lon": (["profile"], vmp.lon),
    "lat": (["profile"], vmp.lat),
    "cast": (["profile"], vmp.cast.astype(int)),
    "p": (["depth"], p),
    "p_mid": (["depth_mid"], p_mid[:, 0]),
    "section": (["section"], section, {"Variable": "Section number"}),
    "insection": (["section", "profile"], insection, {"Variable": "Section masks"}),
    "time_start": (["section"], time_start, {"Variable": "Section start time"}),
    "time_end": (["section"], time_end, {"Variable": "Section end time"}),
    "lon_start": (["section"], lon_start, {"Variable": "Section start longitude"}),
    "lon_end": (["section"], lon_end, {"Variable": "Section end longitude"}),
    "lat_start": (["section"], lat_start, {"Variable": "Section start latitude"}),
    "lat_end": (["section"], lat_end, {"Variable": "Section end latitude"}),
}

vmp_ds = xr.Dataset(vmp_datavars, vmp_coords)
file = "vmp_sep_2018.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
vmp_ds.to_netcdf(os.path.join(sdroot, file))

############### COMBINED ################
print("Processing combined VMP-SADCP")
time_win = params.sadcp.time_window
rmax = params.sadcp.rmax
vmax = params.sadcp.vmax
range_min = params.sadcp.range_min

print("Interpolating velocity to VMP stations.")
mask = np.isfinite(vmp.JAC_T)
u, v, w, lon, lat, range_bottom, nav = ADCP.interp_ADCP_2D(
    sadcp,
    mask,
    vmp.depth,
    vmp.lon,
    vmp.lat,
    vmp.time,
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

file = "combo_sep_2018.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
combo_ds.to_netcdf(os.path.join(sdroot, file))
