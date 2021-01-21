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

files = utils.find_files(args, "sep2018kayak")

params = utils.load_parameters()

print("Loading data may trigger Dropbox download!")

############### LOAD DATA ################
print("Loading kayak data.")
kayak = utils.loadmat(files.kayak, check_arrays=True)
adcps = kayak["adcp"]
ctds = kayak["ctd"]
secs = kayak["xsec"]

############### CTD ################
print("Processing CTD")
ctd_depth_min = params.kayak.ctd_depth_min
ctd_depth_max = params.kayak.ctd_depth_max
ctd_depth_spacing = params.kayak.ctd_depth_spacing

ctd_datasets = []
sec_tstart = []
sec_tend = []
npfls = []

ictd = 0
for ctd, sec in zip(ctds, secs):

    ctd = munchify(ctd)

    npfls.append(ctd.time.size)

    near_glacier = np.full_like(ctd.time, False, dtype=bool)
    near_glacier[np.asarray(ctd.iNG, dtype=int) - 1] = True
    fjord = np.full_like(ctd.time, False, dtype=bool)
    fjord[np.asarray(ctd.iF, dtype=int) - 1] = True

    csec = sec["Cind"]
    for a in csec:
        a = np.asarray(a)
        if a.size < 1:
            sec_tstart.append(np.nan)
            sec_tend.append(np.nan)
        elif a.size == 1:
            sec_tstart.append(ctd.time[a - 1] - 60 / 86400)
            sec_tend.append(ctd.time[a - 1] + 60 / 86400)
        else:
            sec_tstart.append(ctd.time[a[0] - 1] - 60 / 86400)
            sec_tend.append(ctd.time[a[-1] - 1] + 60 / 86400)

    ctd.SA = gsw.SA_from_SP(
        ctd.S, ctd.P, ctd.lon[np.newaxis, :], ctd.lat[np.newaxis, :]
    )
    ctd.CT = gsw.CT_from_t(ctd.SA, ctd.T, ctd.P)
    ctd.sig0 = gsw.pot_rho_t_exact(ctd.SA, ctd.T, ctd.P, 0)

    ctd_datavars = {
        "times": (
            ["depth", "profile"],
            utils.datenum_to_datetime(ctd.time_full),
            {"Variable": "Measurement time"},
        ),
        "p": (["depth", "profile"], ctd.P, {"Variable": "Pressure"}),
        "SP": (["depth", "profile"], ctd.S, {"Variable": "Practical salinity"}),
        "t": (["depth", "profile"], ctd.T, {"Variable": "Temperature (in situ)"}),
        "CT": (["depth", "profile"], ctd.CT, {"Variable": "Conservative temperature"}),
        "SA": (["depth", "profile"], ctd.SA, {"Variable": "Absolute salinity"}),
        "sig0": (
            ["depth", "profile"],
            ctd.sig0,
            {"Variable": "Potential density referenced to 0 dbar"},
        ),
        "N2": (["depth", "profile"], ctd.N2, {"Variable": "Buoyancy frequency"}),
        #         "N2_ref": (
        #             ["depth", "profile"],
        #             ctd.N2_ref,
        #             {
        #                 "Variable": "Adiabatically leveled buoyancy frequency, using {:1.0f} dbar bin".format(
        #                     bin_width
        #                 )
        #             },
        #         ),
        "depth_max": (
            ["profile"],
            ctd.castDepth,
            {"Variable": "Maximum valid data depth"},
        ),
        "near_glacier": (
            ["profile"],
            near_glacier,
            {"Variable": "flag for near glacier profiles"},
        ),
        "fjord": (
            ["profile"],
            fjord,
            {"Variable": "flag for fjord profiles (away from glacier)"},
        ),
        "u": (["depth", "profile"], ctd.u, {"Variable": "Eastward velocity"}),
        "v": (["depth", "profile"], ctd.v, {"Variable": "Northward velocity"}),
        "w": (["depth", "profile"], ctd.w, {"Variable": "Vertical velocity"}),
        "dudz": (
            ["depth", "profile"],
            ctd.dudz,
            {"Variable": "Eastward vertical shear"},
        ),
        "dvdz": (
            ["depth", "profile"],
            ctd.dvdz,
            {"Variable": "Northward vertical shear"},
        ),
    }

    ctd_coords = {
        "cast_number": (["profile"], np.arange(ctd.time.size) + 1),
        "depth": (["depth"], ctd.depth),
        "time": (["profile"], utils.datenum_to_datetime(ctd.time)),
        "lon": (["profile"], ctd.lon),
        "lat": (["profile"], ctd.lat),
        "x": (["profile"], ctd.X),
        "y": (["profile"], ctd.Y),
        "deployment": (["profile"], (ictd + 1) * np.ones_like(ctd.time, dtype=int)),
    }

    ictd += 1

    ds = xr.Dataset(ctd_datavars, ctd_coords).interp(
        depth=np.arange(
            ctd_depth_min, ctd_depth_max + ctd_depth_spacing, ctd_depth_spacing
        )
    )

    ctd_datasets.append(ds)

ctd_ds = xr.concat(ctd_datasets, dim="profile")
ctd_ds["profile"] = (["profile"], np.arange(ctd_ds.profile.size) + 1)

# Deal with section information
sec_tstart = np.asarray(sec_tstart)
sec_tend = np.asarray(sec_tend)
ctd_ds["section"] = (["section"], np.arange(len(sec_tstart)) + 1)
ctd_ds["time_start"] = (
    ["section"],
    utils.datenum_to_datetime(sec_tstart),
    {"Variable": "Section start time"},
)
ctd_ds["time_end"] = (
    ["section"],
    utils.datenum_to_datetime(sec_tend),
    {"Variable": "Section end time"},
)

insections = []
csum_npfls = np.cumsum(npfls)
offset = np.hstack((0, csum_npfls[:-1] - 1))
for ideploy in np.unique(ctd_ds.deployment):
    i = ideploy - 1

    csec = secs[i]["Cind"]
    for a in csec:
        a = np.asarray(a)
        insec = np.full_like(ctd_ds.profile, False, dtype=bool)

        if a.size > 0:
            insec[offset[i] + a - 1] = True

        insections.append(insec)

ctd_ds["insection"] = (["section", "profile"], np.stack(insections))

file = "kayak_ctd_sep_2018.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
ctd_ds.to_netcdf(os.path.join(sdroot, file))

############### ADCP ################
print("Processing ADCP")
adcp_depth_min = params.kayak.adcp_depth_min
adcp_depth_max = params.kayak.adcp_depth_max
adcp_depth_spacing = params.kayak.adcp_depth_spacing

adcp_datasets = []
sec_tstart = []
sec_tend = []
npfls = []

iadcp = 0
for adcp, sec in zip(adcps, secs):

    adcp = munchify(adcp)

    npfls.append(adcp.mtime.size)

    asec = sec["Aind"]
    for a in asec:
        a = np.asarray(a)

        sec_tstart.append(adcp.mtime[a[0] - 1] - 60 / 86400)
        sec_tend.append(adcp.mtime[a[-1] - 1] + 60 / 86400)

    adcp_datavars = {
        "u": (["range", "i"], adcp.u_cut, {"Variable": "Eastward velocity"}),
        "v": (["range", "i"], adcp.v_cut, {"Variable": "Northward velocity"}),
        "w": (["range", "i"], adcp.w_cut, {"Variable": "Vertical velocity"}),
        "us": (
            ["range", "i"],
            adcp.u_sm,
            {"Variable": "Eastward velocity smoothed with 17 second boxcar"},
        ),
        "vs": (
            ["range", "i"],
            adcp.v_sm,
            {"Variable": "Northward velocity smoothed with 17 second boxcar"},
        ),
        "ws": (
            ["range", "i"],
            adcp.w_sm,
            {"Variable": "Vertical velocity smoothed with 17 second boxcar"},
        ),
        #         "err1": (["range", "i"], adcp.vel[:, 3, :], {"Variable": "ADCP Error 1"}),
        #         "err2": (["range", "i"], adcp.vel[:, 4, :], {"Variable": "ADCP Error 2"}),
        "dudz": (["depth", "i"], adcp.dudz, {"Variable": "Eastward vertical shear"}),
        "dvdz": (["depth", "i"], adcp.dvdz, {"Variable": "Northward vertical shear"}),
    }

    adcp_coords = {
        "range": (["range"], adcp.config.ranges),
        "time": (["i"], utils.datenum_to_datetime(adcp.mtime)),
        "lon": (["i"], adcp.gps.lon),
        "lat": (["i"], adcp.gps.lat),
        "x": (["i"], adcp.gps.X),
        "y": (["i"], adcp.gps.Y),
        "deployment": (["i"], (iadcp + 1) * np.ones_like(adcp.mtime, dtype=int)),
    }

    iadcp += 1

    ds_ = xr.Dataset(adcp_datavars, adcp_coords).interp(
        range=np.arange(
            adcp_depth_min, adcp_depth_max + adcp_depth_spacing, adcp_depth_spacing
        )
    )

    adcp_datasets.append(ds_)

adcp_ds = xr.concat(adcp_datasets, dim="i")

sec_tstart = np.asarray(sec_tstart)
sec_tend = np.asarray(sec_tend)

adcp_ds["section"] = (["section"], np.arange(len(sec_tstart)) + 1)
adcp_ds["time_start"] = (
    ["section"],
    utils.datenum_to_datetime(sec_tstart),
    {"Variable": "Section start time"},
)
adcp_ds["time_end"] = (
    ["section"],
    utils.datenum_to_datetime(sec_tend),
    {"Variable": "Section end time"},
)

insections = []
csum_npfls = np.cumsum(npfls)
offset = np.hstack((0, csum_npfls[:-1] - 1))
for ideploy in np.unique(adcp_ds.deployment):
    i = ideploy - 1

    asec = secs[i]["Aind"]
    for a in asec:
        a = np.asarray(a)

        insec = np.full_like(adcp_ds.i, False, dtype=bool)
        insec[offset[i] + a - 1] = True
        insections.append(insec)

adcp_ds["insection"] = (["section", "i"], np.stack(insections))

file = "kayak_adcp_sep_2018.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
adcp_ds.to_netcdf(os.path.join(sdroot, file))
