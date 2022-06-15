import os

import xarray as xr
from munch import munchify

import clargs
import utils

############### PRELIMINARIES ################
# Parse command line arguments
args = clargs.parse_check_args()

sdroot = args.save

files = utils.find_files(args, "runoff")

params = utils.load_parameters()

############### LOAD DATA ################
print("Loading data. (---> This may trigger Dropbox download <---)")
print("Loading extended (2018) discharge.")
runoff = munchify(utils.loadmat(files.sep2018, check_arrays=True, mat_dtype=True))

############### SEP 2018 ################
ds = xr.Dataset(dict(low_scen=(['time'], runoff.low_scen), middle_scen=(['time'], runoff.middle_scen), high_scen=(['time'], runoff.high_scen)), dict(time=(['time'], utils.datenum_to_datetime(runoff.time))))
file = "runoff_extended_2018.nc"
print("Saving to '{}'".format(os.path.join(sdroot, file)))
ds.to_netcdf(os.path.join(sdroot, file))