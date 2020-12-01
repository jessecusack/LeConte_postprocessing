import gsw
import numpy as np
import mixsea as mx
from munch import Munch
from tqdm import tqdm


def generate_CTD_Munch(time, depth, lon, lat, SP, t, C=None, depth_min=1.0, depth_max=300.0, depth_spacing=1.0):
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
        C : numpy array
            Conductivity, size N*M.
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
    if C is not None:
        C = C[:, use]
    
    # Interpolate to a common depth grid if necessary.
    depth_ = np.arange(depth_min, depth_max + depth_spacing, depth_spacing)
    
    size_equal = depth_.size == depth.size
    all_close = np.allclose(depth, depth_)
    
    if size_equal and all_close:
        pass
    else:
        SP = utils.nan_interp(depth_, depth, SP, axis=0)
        t = utils.nan_interp(depth_, depth, t, axis=0)
        depth = depth_
        if C is not None:
            C = utils.nan_interp(depth_, depth, C, axis=0)
        
    ctd = Munch()
    
    ctd.time = time
    ctd.depth = depth
    ctd.lon = lon
    ctd.lat = lat
    ctd.SP = SP
    ctd.t = t
    if C is not None:
        ctd.C = C
        
    return ctd


def apply_thermodynamics(ctd):
    ctd.p, ctd.SA, ctd.CT, ctd.sig0, p_mid, ctd.N2 = common_thermodynamics(ctd.depth, ctd.lon, ctd.lat, ctd.SP, ctd.t)
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
        dmax[i] = depth[mask[:, i]].max()
        
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
    N2_ref = np.full_like(t, np.nan)
    for i in tqdm(range(t.shape[1])):
        N2_ref[:, i] = mx.nsq.adiabatic_leveling(
            p, SP[:, i], t[:, i], lon[i], lat[i], bin_width=bin_width
        )
        
    return N2_ref