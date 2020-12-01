import gsw
import numpy as np


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