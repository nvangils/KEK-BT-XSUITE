"""
Calculate Action Amplitude from Phase Space Coordinates
=============================================

"""

################################################################################
# Required Packages
################################################################################
import numpy as np

################################################################################
# Calculation Function
################################################################################
# def calculate_action_amplitude(a, pa, beta, alfa, emitt_a, mode):
#     """
#     Calculate the action amplitude from phase space coordinates.

#     Parameters
#     ----------
#     a : float or np.ndarray
#         The position coordinate.
#     pa : float or np.ndarray
#         The momentum coordinate.
#     beta : float
#         The beta function at the location of interest.
#     alfa : float
#         The alpha function at the location of interest.
#     emitt_a : float
#         The emittance in the plane of interest.
#     mode : str
#         The mode of amplitude calculation. Either "Norm" for normalized amplitude
#         or "SI" for SI units amplitude.

#     Returns
#     -------
#     amplitude : float or np.ndarray
#         The calculated action amplitude in the specified mode.
#     """

#     assert mode in ["Norm", "SI"], "Mode must be 'Norm' or 'SI'"

#     a_norm2_si  = a**2 / beta
#     pa_norm2_si = (alfa * a + beta * pa)**2 / beta
#     #pa_norm2_si = 0 #(alfa * a + beta * pa)**2 / beta

#     pa_norm2    = pa_norm2_si / emitt_a
#     a_norm2     = a_norm2_si / emitt_a
    
#     # Calculate action
#     Ja_si       = 0.5 * (a_norm2_si + pa_norm2_si)
#     Ja_norm     = 0.5 * (a_norm2 + pa_norm2)

#     # # Calculate action angle
#     # phi         = np.arctan2(
#     #     np.sign(a) * np.sqrt(pa_norm2),
#     #     np.sign(alfa * a + beta * pa) * np.sqrt(a_norm2))

#     # Calculate amplitude
#     Aa_si       = np.sqrt(2 * Ja_si)
#     Aa_norm     = np.sqrt(2 * Ja_norm)

#     if mode == "Norm":
#         return Aa_norm
#     else:
#         return Aa_si

def calculate_action_amplitude(a, pa, delta, beta, alfa, da, dpa, gemitt_a, dda, ddpa):
    """
    Calculate the action amplitude from phase space coordinates.

    Parameters
    ----------
    a : float or np.ndarray
        The position coordinate.
    pa : float or np.ndarray
        The momentum coordinate.
    delta : float or np.ndarray
        The relative momentum deviation.
    beta : float
        The beta function at the location of interest.
    alfa : float
        The alpha function at the location of interest.
    da: float
        The dispersion at the location of interest.
    dpa: float
        The derivative of the dispersion at the location of interest.
    gemitt_a : float
        The emittance in the plane of interest.

    Returns
    -------
    amplitude : float or np.ndarray
        The calculated action amplitude in the specified mode.
    """

    # Print some details before
    ##print("Initial mean a:      ", np.mean(a))
    ##print("Initial mean pa:     ", np.mean(pa))
    ##print("Initial mean delta:  ", np.mean(delta))
    ##print("Initial std a:       ", np.std(a))
    ##print("Initial std pa:      ", np.std(pa))
    ##print("Initial std delta:   ", np.std(delta))

    # Perform Dispersion Shear
    a_shear     = a  - da  * delta -1/2*dda * delta**2
    pa_shear    = pa - dpa * delta-1/2*ddpa * delta**2
    
    ##print("Sheared mean a:      ", np.mean(a_shear))
    ##print("Sheared mean pa:     ", np.mean(pa_shear))
    ##print("Sheared std a:       ", np.std(a_shear))
    ##print("Sheared std pa:      ", np.std(pa_shear))

    # Calculate Twiss gamma
    gama        = (1.0 + alfa**2) / beta
    
    ##print("Mean action aa cont.:    ", np.mean(gama * a_shear**2))
    ##print("Mean action apa cont.:   ", np.mean(2 * alfa * a_shear * pa_shear))
    ##print("Mean action ppa cont.:   ", np.mean(beta * pa_shear**2))

    # normalized action  J / eps_x
    Ja = (gama * a_shear**2 + \
        2 * alfa * a_shear * pa_shear + \
        beta * pa_shear**2) / (2.0 * gemitt_a)
    
    ##print("Mean action aa cont.:    ", np.mean(gama * a_shear**2) / gemitt_a)
    ##print("Mean action apa cont.:   ", np.mean(2 * alfa * a_shear * pa_shear) / gemitt_a)
    ##print("Mean action ppa cont.:   ", np.mean(beta * pa_shear**2) / gemitt_a)
    
    # Convert to amplitude
    Aa  = np.sqrt(2 * Ja)
    
    return Aa

    