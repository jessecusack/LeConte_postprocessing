import os
from datetime import datetime, timedelta

import numpy as np
import scipy.io as io
import utm
import yaml
from munch import Munch, munchify
from scipy.ndimage import median_filter


def loadmat(filename, check_arrays=False, **kwargs):
    """
    Big thanks to mergen on stackexchange for this:
        http://stackoverflow.com/a/8832212

    This function should be called instead of direct scipy.io.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects.
    """
    kwargs["struct_as_record"] = False
    kwargs["squeeze_me"] = True
    data = io.loadmat(filename, **kwargs)
    return _check_keys(data, check_arrays)


def _check_keys(dict, check_arrays):
    """
    Checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries.
    """
    for key in dict:
        if isinstance(dict[key], io.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
        if isinstance(dict[key], np.ndarray) and check_arrays:
            shape = dict[key].shape
            array = dict[key].flatten()
            for i, item in enumerate(array):
                if isinstance(item, io.matlab.mio5_params.mat_struct):
                    array[i] = _todict(item)
            dict[key] = array.reshape(shape)
    return dict


def _todict(matobj):
    """
    A recursive function which constructs from matobjects nested dictionaries.
    """
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, io.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict


def datenum_to_datetime(datenum):
    """
    Convert a MATLAB datenums into python datetimes.

    Parameters
    ----------
    datenum : array_like
        MATLAB datenumber which is the number of days since 0000-01-00.

    Returns
    -------
    dt : ndarray
        Python datetime. See datetime module.

    """

    def convert(datenum):
        try:
            return (
                datetime.fromordinal(int(datenum))
                + timedelta(days=datenum % 1)
                - timedelta(days=366)
            )
        except ValueError:
            return np.nan

    if np.iterable(datenum):
        datenumar = np.asarray(datenum)
        shape = datenumar.shape
        dt = np.array([convert(el) for el in datenumar.flat])
        dt = dt.reshape(shape)
    else:
        dt = convert(datenum)

    return dt


def mid(x, axis=0):
    """Returns mid point values along given axis."""
    ndim = np.ndim(x)
    if ndim == 1:
        return 0.5 * (x[1:] + x[:-1])
    elif ndim > 1:
        x_ = np.swapaxes(x, axis, 0)
        xmid_ = 0.5 * (x_[1:, ...] + x_[:-1, ...])
        return np.swapaxes(xmid_, 0, axis)
    else:
        raise ValueError


def nan_interp(x, xp, fp, left=None, right=None, axis=0, squeeze_me=True):
    """See numpy.interp documentation. This does the same thing but ignores NaN
    values in the data. It can accept 2D arrays.

    Parameters
    ----------
    x : float or 1D array
        The x-coordinates of the interpolated values. No NaNs please!
    xp : 1D or 2D array of floats
        The x-coordinates of the data points, must be increasing along the
        dimension along which the interpolation is being performed.
    fp : 1D or 2D array of floats or complex
        The y-coordinates of the data points, same shape as `xp`.
    left : optional float or complex corresponding to fp
        Value to return for `x < xp[0]`, default is `fp[0]`.
    right : optional float or complex corresponding to fp
        Value to return for `x > xp[-1]`, default is `fp[-1]`.
    axis : [-1, 0, 1] int
        Default is 0. The axis along which to perform the interpolation.
    squeeze_me : boolean
        Default is True. Squeeze output to remove singleton dimensions.

    Returns
    -------
    y : ndarray
        The interpolated values.
    """

    if axis not in [-1, 0, 1]:
        raise ValueError("The axis may be only -1, 0 or 1.")

    if xp.shape != fp.shape:
        raise ValueError("xp and fp have different shapes.")

    ndim = np.ndim(xp)
    if ndim > 2:
        raise ValueError("Only 1 or 2 dimensional arrays are supported.")

    nans = np.isnan(xp) | np.isnan(fp)

    if ndim == 1:
        y = np.full_like(x, np.nan)
        y = np.interp(x, xp[~nans], fp[~nans], left, right)
    if ndim == 2:
        nr, nc = xp.shape

        if axis == 0:
            if np.iterable(x):
                y = np.full((len(x), nc), np.nan)
            else:
                y = np.full((1, nc), np.nan)

            for i in range(nc):
                xp_ = xp[~nans[:, i], i]
                fp_ = fp[~nans[:, i], i]
                y[:, i] = np.interp(x, xp_, fp_, left, right)

        if axis == -1 or axis == 1:
            if axis == 0:
                if np.iterable(x):
                    y = np.full((nr, len(x)), np.nan)
                else:
                    y = np.full((nr, 1), np.nan)

            for i in range(nr):
                xp_ = xp[i, ~nans[i, :]]
                fp_ = fp[i, ~nans[i, :]]
                y[i, :] = np.interp(x, xp_, fp_, left, right)

    if squeeze_me:
        return np.squeeze(y)
    else:
        return y


def interp_fill_valid_2D(x, xp, fp):
    """
    Assumes input values fp is 2D with size N*M,
    where M denotes profiles and N depths.

    Parameters
    ----------
        x : numpy array
            Locations to interpolate to, 1D.
        xp : numpy array
            Data locations, 1D or 2D, shape (N) or (N, M).
        fp : numpy array
            Data values, 2D, shape (N, M).

    """
    nc = fp.shape[1]
    nr = x.size
    f = np.full((nr, nc), np.nan)
    if np.ndim(xp) == 1:
        for i in range(nc):
            f[:, i] = interp_fill_valid(x, xp, fp[:, i])
    elif np.ndim(xp) == 2:
        for i in range(nc):
            f[:, i] = interp_fill_valid(x, xp[:, i], fp[:, i])
    else:
        raise ValueError("xp dimensions are wrong.")
    return f


def interp_fill_valid(x, xp, fp):
    """Interpolate to x, and invalid regions with NaN. Region to fill is
    that out of range of max(xp) and min(xp)."""
    valid = np.isfinite(fp)

    if any(valid):
        xmax = np.max(xp[valid])
        xmin = np.min(xp[valid])
        f = np.interp(x, xp[valid], fp[valid])
        f[(x > xmax) | (x < xmin)] = np.nan
    else:
        f = fp

    return f


def check_files(files):
    """
    Assumes that files is a dict or Munch object containing full file paths.
    """
    for key in files:
        # Skip check for non-string objects.
        if type(files[key]) != str:
            continue

        if not os.path.isfile(files[key]):
            raise ValueError("{} file not found: {}".format(key, files[key]))
        else:
            print("Found {} at '{}'.".format(key, files[key]))


def find_files(args, dataset, paths_file="file_paths.yml"):
    """
    args: command line args
    dataset: yaml file path parameter key e.g. "sep2018"

    returns files as Munch

    """

    # Grab the data file paths from the yml file.
    with open(paths_file, "r") as f:
        try:
            all_files = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            print(exc)

    # Grab file path info for the specified dataset
    file_info = munchify(all_files[dataset])
    files = Munch()

    # Join root directory specified by command line arguments with path
    # specified in the yaml file.
    for key in file_info:
        files[key] = os.path.join(args[file_info[key].root], file_info[key].path)

    check_files(files)

    return files


def load_parameters(parameter_file="processing_parameters.yml"):
    """Load processing parameters into Munch."""

    with open(parameter_file, "r") as f:
        try:
            params = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            print(exc)

    return munchify(params)


def closest_index(x, a):
    """
    x: value
    a: array
    """
    return np.argmin(np.abs(x - a))


def regrid_profiles(time, timep, fp, time_win=60.0):
    """"""

    dt = time_win / (2 * 86400)

    nc = time.size

    idxs = []
    idxps = []

    for i in range(nc):
        time_diff = np.abs(timep - time[i])
        time_min = np.min(time_diff)

        # Skip if not within time window
        if time_min > dt:
            continue

        idx = np.argmin(time_diff)

        # Skip if already found
        if idx in idxs:
            continue

        idxs.append(i)
        idxps.append(idx)

    idxs = np.asarray(idxs)
    idxps = np.asarray(idxps)

    ndim = np.ndim(fp)
    if ndim == 1:
        f = np.full_like(time, np.nan)
        f[idxs] = fp[idxps]
    elif ndim == 2:
        nr = fp.shape[0]
        f = np.full((nr, nc), np.nan)
        f[:, idxs] = fp[:, idxps]

    return f


def apply_utm(m):
    """m is a Munch object"""
    m.x, m.y, m.zone_number, m.zone_letter = utm.from_latlon(m.lat, m.lon)
    return m


def rolling_window(a, size):
    pad = np.ones(len(a.shape), dtype=np.int32)
    pad[-1] = size - 1
    pad = list(zip(pad, np.zeros(len(a.shape), dtype=np.int32)))
    a = np.pad(a, pad, mode="reflect")
    shape = a.shape[:-1] + (a.shape[-1] - size + 1, size)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def despike(x, n1=3, n2=6, size=101, xmin=-1.0, xmax=1.0, fill=True):
    """
    Despike data using a median filter and 2 pass standard deviation threshold.

    Parameters
    ----------
        x : numpy array
            Data, evenly spaced.
        n1 : float
            Pass 1 significance threshold, n standard deviations from reference data.
        n2 : float
            Pass 2 significance threshold, n standard deviations from reference data.
        size : float, optional
            Number of data points in the filter window.
        xmin : float, optional
            Minimum value of x, data below this value will be removed.
        xmax : float, optional
            Maximum value of x, data above this value will be removed.
        fill : boolean, optional
            Fill spikes using linear interpolation.

    """

    # First get rid of gaps using linear interpolation
    nans = np.isnan(x)
    if any(nans):
        t_ = np.arange(x.size)
        x[nans] = np.interp(t_[nans], t[~nans], x[~nans])

    # Moving median and std pass 1
    roll = rolling_window(x, size)
    x1med = np.median(roll, axis=-1)
    x1std = np.std(roll, axis=-1)
    # Spikes using first threshold
    dx1 = x - x1med
    spikes1 = np.abs(dx1) > n1 * x1std

    # Mask out spikes from first pass
    xm = np.ma.masked_where(spikes1, x)

    # Moving median and std pass 2
    roll = rolling_window(xm, size)
    x2med = np.ma.median(roll, axis=-1)
    x2std = np.ma.std(roll, axis=-1)
    dx2 = x - x2med
    spikes = np.abs(dx2) > n2 * x2std

    # Trim min and max
    trim = (x > xmax) | (x < xmin)
    spikes[trim] = True

    if fill:
        t_ = np.arange(x.size)
        x[spikes] = np.interp(t_[spikes], t_[~spikes], x[~spikes])
    else:
        x[spikes] = np.nan

    return x
