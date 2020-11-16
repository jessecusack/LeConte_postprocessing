# Interpolation of ADCP to location of CTD/VMP profiles

import gsw
import numpy as np
from munch import Munch
from tqdm import tqdm


def interp_ADCP_2D(
    sadcp,
    mask,
    depth,
    lon,
    lat,
    time,
    time_win=360.0,
    rmax=15.0,
    vmax=2.0,
    range_min=4.0,
):
    """
    This is essentially a loop over the interp_ADCP function with some additional NaN handling.

    Assume data is of the form D[i, j] where each j represents a profile and i a depth in that profile.


    Parameters
    ----------
        sadcp : Munch
            Munch structure of sadcp data
        mask : 2D array
            Mask of boolean values specifying valid depths to interpolate to.
        depth : array
            Depths (m) at which to interpolate ADCP data.
        lon : array
            Longitude of CTD/VMP profile.
        lat : array
            Latitude of CTD/VMP profile.
        time : array
            Time of CTD/VMP profile as matlab datenum.
        time_win : float, optional
            Time window for search (s) centered on time of profile. Data outside the time range is excluded.
        rmax : float, optional
            Distance threshold (m) defines a circle around the location of the profile. Data outside the circle is excluded.
        vmax : float, optional
            Velocity threshold (m/s) above which we remove velocity data
        range_min : float, optional
            ADCP minimum range threshold (m) below which we remove data

    Return
    ------
        u : 2D array
            Zonal velocity (m/s) interpolated to given depths.
        v : 2D array
            Meridional velocity (m/s) interpolated to given depths.
        w : 2D array
            Vertical velocity (m/s) interpolated to given depths.
        lonm : array
            Mean longitude of ADCP data.
        latm : array
            Mean latitude of ADCP data.
        range_bottom : array
            Minimum beam range to bottom (m).
        n : array
            Number of ADCP profiles in average.

    """

    u = np.full_like(mask, np.nan, dtype=float)
    v = np.full_like(mask, np.nan, dtype=float)
    w = np.full_like(mask, np.nan, dtype=float)
    lonm = np.full_like(time, np.nan)
    latm = np.full_like(time, np.nan)
    range_bottom = np.full_like(time, np.nan)
    n = np.full_like(time, np.nan)

    for i in tqdm(range(time.size)):

        valid = mask[:, i]

        try:
            u_, v_, w_, lon_, lat_, range_bottom_, n_ = interp_ADCP(
                sadcp,
                depth[valid],
                lon[i],
                lat[i],
                time[i],
                time_win=time_win,
                rmax=rmax,
                vmax=vmax,
                range_min=range_min,
            )
        except RuntimeError as err:
            continue

        # Fill data
        u[valid, i] = u_
        v[valid, i] = v_
        w[valid, i] = w_
        lonm[i] = lon_
        latm[i] = lat_
        range_bottom[i] = range_bottom_
        n[i] = n_

    return u, v, w, lonm, latm, range_bottom, n


def interp_ADCP(
    sadcp, depth, lon, lat, time, time_win=360.0, rmax=15.0, vmax=2.0, range_min=4.0
):
    """
    Parameters
    ----------
        sadcp : Munch
            Munch structure of sadcp data
        depth : array
            Depths (m) at which to interpolate ADCP data.
        lon : float
            Longitude of CTD/VMP profile.
        lat : float
            Latitude of CTD/VMP profile.
        time : float
            Time of CTD/VMP profile as matlab datenum.
        time_win : float, optional
            Time window for search (s) centered on time of profile. Data outside the time range is excluded.
        rmax : float, optional
            Distance threshold (m) defines a circle around the location of the profile. Data outside the circle is excluded.
        vmax : float, optional
            Velocity threshold (m/s) above which we remove velocity data
        range_min : float, optional
            ADCP minimum range threshold (m) below which we remove data

    Return
    ------
        u : array
            Zonal velocity (m/s) interpolated to given depths.
        v : array
            Meridional velocity (m/s) interpolated to given depths.
        w : array
            Vertical velocity (m/s) interpolated to given depths.
        lonm : float
            Mean longitude of ADCP data.
        latm : float
            Mean latitude of ADCP data.
        range_bottom : float
            Minimum beam range to bottom (m).
        n : int
            Number of ADCP profiles in average.

    """

    # Convert to matlab datenum units
    dt = time_win / (2 * 86400)

    # Because of wacky conversion of matlab structure arrays to python
    # types we need to check if the data is a dict (or Munch) of lists
    sadcp = check_record(sadcp, time, dt)

    intime = (sadcp.mtime > (time - dt)) & (sadcp.mtime < (time + dt))

    if not intime.any():
        raise RuntimeError("No valid data in time range.")

    # Cut out selected time
    u_ = sadcp.vel[:, 0, intime]
    v_ = sadcp.vel[:, 1, intime]
    w_ = sadcp.vel[:, 2, intime]
    lon_ = sadcp.gps.lon[intime]
    lat_ = sadcp.gps.lat[intime]

    # Horizontal distance between adcp data and the profile
    lons = np.stack((lon_, np.full_like(lon_, lon)), axis=-1)
    lats = np.stack((lat_, np.full_like(lat_, lat)), axis=-1)
    r = gsw.distance(lons, lats, axis=-1).squeeze()

    # Cut out selected range
    inrange = r < rmax

    if not inrange.any():
        raise RuntimeError("No valid data in distance range.")

    u_ = u_[:, inrange]
    v_ = v_[:, inrange]
    w_ = w_[:, inrange]
    lon_ = lon_[inrange]
    lat_ = lat_[inrange]

    n = len(lon_)

    # Remove bottom
    ranges = sadcp.config.ranges
    range_bottom = sadcp.bt_range[:, intime][:, inrange].min()
    below_bottom = ranges > range_bottom

    u_[below_bottom, :] = np.nan
    v_[below_bottom, :] = np.nan
    w_[below_bottom, :] = np.nan

    # Remove spurious velocities
    spurious = (np.abs(u_) > vmax) | (np.abs(v_) > vmax) | (np.abs(w_) > vmax)
    u_[spurious] = np.nan
    v_[spurious] = np.nan
    w_[spurious] = np.nan

    # Average velocity
    ua = np.nanmean(u_, axis=1)
    va = np.nanmean(v_, axis=1)
    wa = np.nanmean(w_, axis=1)

    # Interpolate velocity to finer grid
    valid = np.isfinite(ua)
    u = np.interp(depth, ranges[valid], ua[valid])
    valid = np.isfinite(va)
    v = np.interp(depth, ranges[valid], va[valid])
    valid = np.isfinite(wa)
    w = np.interp(depth, ranges[valid], wa[valid])

    # Remove velocity data out of range
    out_of_range = (depth < range_min) | (depth > range_bottom)
    u[out_of_range] = np.nan
    v[out_of_range] = np.nan
    w[out_of_range] = np.nan

    # Fill data
    lonm = np.nanmean(lon_)
    latm = np.nanmean(lat_)

    return u, v, w, lonm, latm, range_bottom, n


def check_record(sadcp, time, dt):
    """Deal with records that were stored as structure arrays."""

    # If the structure is not a list then we can skip all the rest.
    if type(sadcp.mtime) is not list:
        return sadcp

    # Number of sadcp structures
    ns = len(sadcp.mtime)

    # Check each individual sadcp structure to find data.
    nd = np.zeros((ns))  # Number of data points in each sadcp structure
    for j in range():
        use = (sadcp.mtime[j] > (time - dt)) & (sadcp.mtime[j] < (time + dt))
        if use.any():
            nd[j] = use.sum()

    if all(nd == 0):  # Then we found no data.
        raise RuntimeError("No valid data in time range. (record check)")

    idx = np.argmax(nd)

    sadcp_out = Munch()

    for key in sadcp:
        sadcp_out[key] = sadcp[key][idx]

    return sadcp_out
