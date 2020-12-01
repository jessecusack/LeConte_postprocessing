import mixsea as mx
import numpy as np


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