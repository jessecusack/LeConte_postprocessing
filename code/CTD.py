import gsw
import numpy as np


def common_thermodynamics(depth, lon, lat, SP, t):
    """Wrapper for various thermodynamic calculations.
    
    Assumes input data is 2D with size N*M. 
    
    Parameters
    ----------
        depth : numpy array
            Depth (m), size N.
        lon : numpy array
            Longitude, size M.
        lat : numpy array
            Latitude, size M.
        SP : numpy array
            Practical salinity, size N*M
        t : numpy array
            Temperature, size N*M
    
    """
    
    p = gsw.p_from_z(-depth, np.mean(lat))
    SA = gsw.SA_from_SP(SP, p[:, np.newaxis], lon[np.newaxis, :], lat[np.newaxis, :])
    CT = gsw.CT_from_t(SA, t, p[:, np.newaxis])
    sig0 = gsw.pot_rho_t_exact(SA, t, p[:, np.newaxis], 0)
    N2, p_mid = gsw.Nsquared(SA, CT, p[:, np.newaxis], lat[np.newaxis, :])
    
    return p, SA, CT, sig0, p_mid, N2