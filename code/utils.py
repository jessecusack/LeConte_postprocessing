from datetime import datetime, timedelta

import numpy as np
import scipy.io as io


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
