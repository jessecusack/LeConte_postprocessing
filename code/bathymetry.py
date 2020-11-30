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
# Parse command line arguments
args = clargs.parse_check_args()

sdroot = args.save

files = utils.find_files(args, "bathymetry")

params = utils.load_parameters()

############### LOAD DATA ################
print("Loading data. (---> This may trigger Dropbox download <---)")
print("Loading September 2018 bathymetry.")
b18 = munchify(utils.loadmat(files.sep2018, check_arrays=True, mat_dtype=True)["bathy"])
print("Loading September 2017 bathymetry.")
b17 = munchify(utils.loadmat(files.sep2017, check_arrays=True, mat_dtype=True)["bathy"])
print("Loading August 2016 bathymetry.")
b16 = munchify(utils.loadmat(files.aug2016, check_arrays=True, mat_dtype=True)["bathy"])

############### SEP 2018 ################
print("Processing September 2018.")
datavars = {
    "H": (["lat", "lon"], b18.z, {"Variable": "Bottom depth"}),
}

coords = {
    "lon": (["lon"], b18.lon),
    "lat": (["lat"], b18.lat),
    "x": (["lat", "lon"], b18.x, {"Variable": "UTM zonal distance"}),
    "y": (["lat", "lon"], b18.y, {"Variable": "UTM meridional distance"}),
    "zone_number": ([], int(b18.utmzone[:2]), {"Variable": "UTM zone number"}),
    "zone_letter": ([], b18.utmzone[-1], {"Variable": "UTM zone letter"}),
}

ds = xr.Dataset(datavars, coords)
file = "bathy_sep_2018.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
ds.to_netcdf(os.path.join(sdroot, file))

############### SEP 2017 ################
print("Processing September 2017.")
datavars = {
    "H": (["lat", "lon"], b17.z, {"Variable": "Bottom depth"}),
}

coords = {
    "lon": (["lon"], b17.lon),
    "lat": (["lat"], b17.lat),
    "x": (["lat", "lon"], b17.x, {"Variable": "UTM zonal distance"}),
    "y": (["lat", "lon"], b17.y, {"Variable": "UTM meridional distance"}),
    "zone_number": ([], int(b17.utmzone[:2]), {"Variable": "UTM zone number"}),
    "zone_letter": ([], b17.utmzone[-1], {"Variable": "UTM zone letter"}),
}

ds = xr.Dataset(datavars, coords)
file = "bathy_sep_2017.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
ds.to_netcdf(os.path.join(sdroot, file))

############### AUG 2016 ################
print("Processing August 2016.")
datavars = {
    "H": (["lat", "lon"], b16.z, {"Variable": "Bottom depth"}),
}

coords = {
    "lon": (["lon"], b16.lon),
    "lat": (["lat"], b16.lat),
    "x": (["lat", "lon"], b16.x, {"Variable": "UTM zonal distance"}),
    "y": (["lat", "lon"], b16.y, {"Variable": "UTM meridional distance"}),
    "zone_number": ([], int(b16.utmzone[:2]), {"Variable": "UTM zone number"}),
    "zone_letter": ([], b16.utmzone[-1], {"Variable": "UTM zone letter"}),
}

ds = xr.Dataset(datavars, coords)
file = "bathy_aug_2016.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
ds.to_netcdf(os.path.join(sdroot, file))