import gsw
import mixsea as mx
import numpy as np
from munch import Munch
from tqdm import tqdm

import utils

dvn = Munch(
    {
        "time": "time",
        "C": "C",
        "SP": "S",
        "t": "T",
        "lon": "lon",
        "lat": "lat",
        "depth": "depth",
    }
)


def generate_CTD_Munch_from_list(
    dlist, vn=dvn, depth_min=1.0, depth_max=300.0, depth_spacing=1.0
):
    """Handles the case of a list of dictionary-like objects."""

    ctdlist = []

    for d in dlist:
        ctd_ = generate_CTD_Munch(
            d[vn.time],
            d[vn.depth],
            d[vn.lon],
            d[vn.lat],
            d[vn.SP],
            d[vn.t],
            depth_min,
            depth_max,
            depth_spacing,
        )
        ctdlist.append(ctd_)

    # Stack all the ctds together
    ctd = Munch()
    ctd.time = np.concatenate([ctd_.time for ctd_ in ctdlist], axis=0)
    ctd.depth = ctdlist[0].depth
    ctd.lon = np.concatenate([ctd_.lon for ctd_ in ctdlist], axis=0)
    ctd.lat = np.concatenate([ctd_.lat for ctd_ in ctdlist], axis=0)
    ctd.SP = np.concatenate([ctd_.SP for ctd_ in ctdlist], axis=1)
    ctd.t = np.concatenate([ctd_.t for ctd_ in ctdlist], axis=1)

    return ctd


def generate_CTD_Munch(
    time, depth, lon, lat, SP, t, depth_min=1.0, depth_max=300.0, depth_spacing=1.0
):
    """Quality control CTD quantities and place into a common Munch structure.

    Assumes input data is 1D with size N or M, or 2D with size N*M,
    where M denotes profiles and N depths.

    Parameters
    ----------
        time : numpy array
            Time matlab datenum, size M.
        depth : numpy array
            Depth (m), size N.
        lon : numpy array
            Longitude, size M.
        lat : numpy array
            Latitude, size M.
        SP : numpy array
            Practical salinity, size N*M.
        t : numpy array
            Temperature, size N*M.
        depth_min : float, optional
            Minimum depth (m).
        depth_max : float, optional
            Maximum depth (m).
        depth_spacing : float, optional
            Depth spacing (m).

    """

    use = np.isfinite(time)

    # First remove data if there is no valid timestamp.
    time = time[use]
    SP = SP[:, use]
    t = t[:, use]
    lon = lon[use]
    lat = lat[use]

    # Interpolate to a common depth grid if necessary.
    depth_ = np.arange(depth_min, depth_max + depth_spacing, depth_spacing)

    size_equal = depth_.size == depth.size

    if size_equal:
        all_close = np.allclose(depth, depth_)
    else:
        all_close = False

    if size_equal and all_close:
        pass
    else:
        SP = utils.interp_fill_valid_2D(depth_, depth, SP)
        t = utils.interp_fill_valid_2D(depth_, depth, t)
        depth = depth_

    ctd = Munch()

    ctd.time = time
    ctd.depth = depth
    ctd.lon = lon
    ctd.lat = lat
    ctd.SP = SP
    ctd.t = t

    return ctd


def apply_thermodynamics(ctd):
    ctd.p, ctd.SA, ctd.CT, ctd.sig0, p_mid, ctd.N2 = common_thermodynamics(
        ctd.depth, ctd.lon, ctd.lat, ctd.SP, ctd.t
    )
    ctd.p_mid = p_mid[:, 0]
    return ctd


def apply_adiabatic_level(ctd, bin_width=20.0):
    ctd.N2_ref = adiabatic_level_2D(ctd.p, ctd.SP, ctd.t, ctd.lon, ctd.lat, bin_width)
    return ctd


def common_thermodynamics(depth, lon, lat, SP, t):
    """Wrapper for various thermodynamic calculations.

    Assumes input data is 1D with size N or M, or 2D with size N*M,
    where M denotes profiles and N depths.

    Parameters
    ----------
        depth : numpy array
            Depth (m), size N.
        lon : numpy array
            Longitude, size M.
        lat : numpy array
            Latitude, size M.
        SP : numpy array
            Practical salinity, size N*M.
        t : numpy array
            Temperature, size N*M.

    """

    p = gsw.p_from_z(-depth, np.mean(lat))
    SA = gsw.SA_from_SP(SP, p[:, np.newaxis], lon[np.newaxis, :], lat[np.newaxis, :])
    CT = gsw.CT_from_t(SA, t, p[:, np.newaxis])
    sig0 = gsw.pot_rho_t_exact(SA, t, p[:, np.newaxis], 0)
    N2, p_mid = gsw.Nsquared(SA, CT, p[:, np.newaxis], lat[np.newaxis, :])

    return p, SA, CT, sig0, p_mid, N2


def depth_max(depth, mask):
    """Calculate the maximum depth of valid data.

    Assumes mask is 2D with size N*M, where M denotes profiles and N depths.

    Parameters
    ----------
        depth : numpy array
            Depth (m), size N.
        mask : numpy array of boolean
            Mask where True is valid data, size N*M.
    """

    dmax = np.full((mask.shape[1]), np.nan)
    for i in range(mask.shape[1]):
        dmax[i] = np.max(depth[mask[:, i]])

    return dmax


def adiabatic_level_2D(p, SP, t, lon, lat, bin_width=20.0):
    """Adiabatically level 2D data.

    Assumes input data is 1D with size N or M, or 2D with size N*M,
    where M denotes profiles and N depths.

    Parameters
    ----------
        p : numpy array
            Pressure (dbar), size N.
        SP : numpy array
            Practical salinity, size N*M.
        t : numpy array
            Temperature, size N*M.
        lon : numpy array
            Longitude, size M.
        lat : numpy array
            Latitude, size M.
        bin_width : float, optional
            Bin width (dbar).

    """
    print("Adiabatically levelling profiles")
    N2_ref = np.full_like(t, np.nan)
    for i in tqdm(range(t.shape[1])):
        N2_ref[:, i] = mx.nsq.adiabatic_leveling(
            p, SP[:, i], t[:, i], lon[i], lat[i], bin_width=bin_width
        )

    return N2_ref


def regrid_vmp_to_ctd(vmp, ctd, time_win=60.0):
    """Add some CTD variabiles to the VMP Munch."""

    ctd.eps1 = utils.regrid_profiles(ctd.time, vmp.time, vmp.eps1)
    ctd.Lo1 = utils.regrid_profiles(ctd.time, vmp.time, vmp.Lo1)
    ctd.Kv1 = utils.regrid_profiles(ctd.time, vmp.time, vmp.Kv1)

    try:
        ctd.eps2 = utils.regrid_profiles(ctd.time, vmp.time, vmp.eps2)
        ctd.Lo2 = utils.regrid_profiles(ctd.time, vmp.time, vmp.Lo2)
        ctd.Kv2 = utils.regrid_profiles(ctd.time, vmp.time, vmp.Kv2)
    except AttributeError:
        pass

    return ctd
