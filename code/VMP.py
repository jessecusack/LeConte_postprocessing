import mixsea as mx
import numpy as np
import gsw
from munch import Munch

import utils


def generate_VMP_Munch(time, depth, lon, lat, eps1, eps2=None, depth_min=1.0, depth_max=300.0, depth_spacing=1.0):
    """Quality control VMP quantities and place into a common Munch structure.
    
    Assumes input data is 1D with size N or M, or 2D with size N*M, 
    where M denotes profiles and N depths.
    
    Parameters
    ----------
        time : numpy array
            Time matlab datenum, size M.
        depth : numpy array
            Depth (m), size N or N*M.
        lon : numpy array
            Longitude, size M.
        lat : numpy array
            Latitude, size M.
        eps1 : numpy array
            Turbulent dissipation rate, size N*M.
        eps2 : numpy array, optional
            Second dissipation estimate, size N*M.
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
    eps1 = eps1[:, use]
    if eps2 is not None:
        eps2 = eps2[:, use]
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
        # Use logarithms in interp to deal with order of magnitude variation in eps
        eps1 = 10 ** utils.interp_fill_valid_2D(depth_, depth, np.log10(eps1))
        if eps2 is not None:
            eps2 = 10 ** utils.interp_fill_valid_2D(depth_, depth, np.log10(eps2))
        depth = depth_
        
    vmp = Munch()
    
    vmp.time = time
    vmp.depth = depth
    vmp.lon = lon
    vmp.lat = lat
    vmp.eps1 = eps1
    if eps2 is not None:
        vmp.eps2 = eps2
        
    return vmp


def common_turbulence(eps, N2, gam=0.2):
    """Wrapper for various turbulence calculations.
    
    Assumes input data is 1D with size N or M, or 2D with size N*M, 
    where M denotes profiles and N depths.
    
    Parameters
    ----------
        eps : numpy array
            Turbulent dissipation rate (W kg-1), size N*M.
        N2 : numpy array
            Buoyancy frequency (rad s-1), size N*M.
        gam : float, optional
            Mixing efficiency. 
    
    """
    
    # Ozmidov scale
    Lo = np.sqrt(eps / N2 ** (3 / 2))
    # Turbulent diffusivity using the Osborne relation
    Kv = gam * eps / N2
    
    return Lo, Kv


def regrid_ctd_to_vmp(ctd, vmp, time_win=60.):
    """Add some CTD variabiles to the VMP Munch."""
    
    vmp.SP = utils.regrid_profiles(vmp.time, ctd.time, ctd.SP)
    vmp.t = utils.regrid_profiles(vmp.time, ctd.time, ctd.t)
    vmp.SA = utils.regrid_profiles(vmp.time, ctd.time, ctd.SA)
    vmp.CT = utils.regrid_profiles(vmp.time, ctd.time, ctd.CT)
    vmp.N2_ref = utils.regrid_profiles(vmp.time, ctd.time, ctd.N2_ref)
    vmp.sig0 = utils.regrid_profiles(vmp.time, ctd.time, ctd.sig0)
    vmp.depth_max = utils.regrid_profiles(vmp.time, ctd.time, ctd.depth_max)
    vmp.p = ctd.p
    vmp.p_mid = ctd.p_mid
    vmp.N2 = utils.regrid_profiles(vmp.time, ctd.time, ctd.N2)
    
    return vmp