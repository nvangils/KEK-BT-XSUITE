"""
Helpers for lattice testing
=============================================
Author(s):  John P T Salvesen, Nikita Z van Gils
Date:       30-09-2025
"""

#################################### Required Packages ###################################

import xtrack as xt
import numpy as np
import matplotlib.pyplot as plt


##################################### Check Symplecticity  ####################################

def check_symplecticity(twiss, line, tt = None):
    J   = np.array([
        [0, +1, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0],
        [0, 0, 0, +1, 0, 0],
        [0, 0, -1, 0, 0, 0],
        [0, 0, 0, 0, 0, +1],
        [0, 0, 0, 0, -1, 0]])

    
    ##################################### Get table ####################################
    
    if tt is None:
        tt = line.get_table(attr = True)

    
    ###################################### Overall #####################################
    
    M           = twiss.R_matrix
    residual    = (M.T @ J @ M) - J
    for row in residual:
        print(" ".join(f"{x:+12.5E}" for x in row))

    symplectic  = np.allclose(residual, 0, atol = 1E-6)
    print("Overall R Matrix symplecticity (1E-6 level): ", symplectic)
    print("Maximum deviation:                           ", np.max(np.abs(residual)))

    
    ###################################### EBE if there is an issue #####################################
    
    if not symplectic:
        # Need to exclude _end_ponit with -1
        for test_ele, end_ele in zip(tt.name[:-2], tt.name[1:-1]):
            test_particle   = xt.Particles(
                p0c     = line.particle_ref.p0c,
                mass0   = line.particle_ref.mass0,
                q0      = line.particle_ref.q0)
            M_ele   = line.compute_one_turn_matrix_finite_differences(
                start               = test_ele,
                end                 = end_ele,
                particle_on_co      = test_particle,
                steps_r_matrix      = twiss.steps_r_matrix)["R_matrix"]
            residual    = (M_ele.T @ J @ M_ele) - J
            symplectic  = np.allclose(residual, 0, atol = 1E-6)
            if not symplectic:
                print(f"Non-symplectic element (1E-6 level):    {test_ele}")



###################################### Zero small values #####################################

def zero_small_values(array, tol = 1E-12):
    array[np.abs(array) < tol] = 0
    return array


###################################### SAD vs Xsuite Comparison Plots #####################################

def create_comparison_plots(
        twiss_xsuite,
        twiss_sad,
        suptitle    = None,
        zero_tol    = 1E-12,
        figsize     = (8, 4)):

    
    ###################################### Orbit (x, y) #####################################
    
    fig, axs = plt.subplots(2, figsize = figsize, sharex = True)

    axs[0].plot(
        zero_small_values(twiss_sad.s, tol = zero_tol),
        zero_small_values(twiss_sad.x, tol = zero_tol),
        label       = 'SAD',
        color       = "r")
    axs[0].plot(
        zero_small_values(twiss_xsuite.s, tol = zero_tol),
        zero_small_values(twiss_xsuite.x, tol = zero_tol),
        label       = 'Xsuite',
        color       = "b",
        linestyle   = "--")
    axs[1].plot(
        zero_small_values(twiss_sad.s, tol = zero_tol),
        zero_small_values(twiss_sad.y, tol = zero_tol),
        color       = "r")
    axs[1].plot(
        zero_small_values(twiss_xsuite.s, tol = zero_tol),
        zero_small_values(twiss_xsuite.y, tol = zero_tol),
        color       = "b",
        linestyle   = "--")
    
    axs[0].legend()
    axs[0].set_ylabel(r'$x$ [m]')
    axs[1].set_ylabel(r'$y$ [m]')
    axs[1].set_xlabel('s [m]')

    if suptitle is not None:
        fig.suptitle(f"{suptitle}: Orbit (x, y)")
    else:
        fig.suptitle("Orbit (x, y)")
    fig.tight_layout()
    fig.align_labels()
    fig.align_titles()

    
    ###################################### Orbit (px, py) #####################################
    
    fig, axs = plt.subplots(2, figsize = figsize, sharex = True)

    axs[0].plot(
        zero_small_values(twiss_sad.s, tol = zero_tol),
        zero_small_values(twiss_sad.px, tol = zero_tol),
        label       = 'SAD',
        color       = "r")
    axs[0].plot(
        zero_small_values(twiss_xsuite.s, tol = zero_tol),
        zero_small_values(twiss_xsuite.px, tol = zero_tol),
        label       = 'Xsuite',
        color       = "b",
        linestyle   = "--")
    axs[1].plot(
        zero_small_values(twiss_sad.s, tol = zero_tol),
        zero_small_values(twiss_sad.py, tol = zero_tol),
        color       = "r")
    axs[1].plot(
        zero_small_values(twiss_xsuite.s, tol = zero_tol),
        zero_small_values(twiss_xsuite.py, tol = zero_tol),
        color       = "b",
        linestyle   = "--")
    
    axs[0].legend()
    axs[0].set_ylabel(r'$p_{x}$ [m]')
    axs[1].set_ylabel(r'$p_{y}$ [m]')
    axs[1].set_xlabel('s [m]')

    if suptitle is not None:
        fig.suptitle(f"{suptitle}: Orbit (px, py)")
    else:
        fig.suptitle("Orbit (px, py)")
    fig.tight_layout()
    fig.align_labels()
    fig.align_titles()

    
    ###################################### Longitudinal Plane (zeta, delta) #####################################
    
    fig, axs = plt.subplots(2, figsize = figsize, sharex = True)

    axs[0].plot(
        zero_small_values(twiss_sad.s, tol = zero_tol),
        zero_small_values(twiss_sad.zeta, tol = zero_tol),
        label       = 'SAD',
        color       = "r")
    axs[0].plot(
        zero_small_values(twiss_xsuite.s, tol = zero_tol),
        zero_small_values(twiss_xsuite.zeta, tol = zero_tol),
        label       = 'Xsuite',
        color       = "b",
        linestyle   = "--")
    axs[1].plot(
        zero_small_values(twiss_sad.s, tol = zero_tol),
        zero_small_values(twiss_sad.delta, tol = zero_tol),
        color       = "r")
    axs[1].plot(
        zero_small_values(twiss_xsuite.s, tol = zero_tol),
        zero_small_values(twiss_xsuite.delta, tol = zero_tol),
        color       = "b",
        linestyle   = "--")
    
    axs[0].legend()
    axs[0].set_ylabel(r'$\zeta$ [m]')
    axs[1].set_ylabel(r'$\delta$ [m]')
    axs[1].set_xlabel('s [m]')

    if suptitle is not None:
        fig.suptitle(f"{suptitle}: Longitudinal Plane ($\zeta$, $\delta$)")
    else:
        fig.suptitle("Longitudinal Plane ($\zeta$, $\delta$)")
    fig.tight_layout()
    fig.align_labels()
    fig.align_titles()

    
    ###################################### Beta Functions #####################################
    
    fig, axs = plt.subplots(2, figsize = figsize, sharex = True)

    axs[0].plot(
        zero_small_values(twiss_sad.s, tol = zero_tol),
        zero_small_values(twiss_sad.betx, tol = zero_tol),
        label       = 'SAD',
        color       = "r")
    axs[0].plot(
        zero_small_values(twiss_xsuite.s, tol = zero_tol),
        zero_small_values(twiss_xsuite.betx, tol = zero_tol),
        label       = 'Xsuite',
        color       = "b",
        linestyle   = "--")
    axs[1].plot(
        zero_small_values(twiss_sad.s, tol = zero_tol),
        zero_small_values(twiss_sad.bety, tol = zero_tol),
        color       = "r")
    axs[1].plot(
        zero_small_values(twiss_xsuite.s, tol = zero_tol),
        zero_small_values(twiss_xsuite.bety, tol = zero_tol),
        color       = "b",
        linestyle   = "--")

    axs[0].legend()
    axs[0].set_ylabel(r'$\beta_{x}$ [m]')
    axs[1].set_ylabel(r'$\beta_{y}$ [m]')
    axs[1].set_xlabel('s [m]')

    if suptitle is not None:
        fig.suptitle(f"{suptitle}: " + "Beta Functions ($\\beta_{x}$, $\\beta_{y}$)")
    else:
        fig.suptitle("Beta Functions ($\\beta_{x}$, $\\beta_{y}$)")
    fig.tight_layout()
    fig.align_labels()
    fig.align_titles()

    
    ##################################### Alpha Functions #####################################
    
    fig, axs = plt.subplots(2, figsize = figsize, sharex = True)

    axs[0].plot(
        zero_small_values(twiss_sad.s, tol = zero_tol),
        zero_small_values(twiss_sad.alfx, tol = zero_tol),
        label       = 'SAD',
        color       = "r")
    axs[0].plot(
        zero_small_values(twiss_xsuite.s, tol = zero_tol),
        zero_small_values(twiss_xsuite.alfx, tol = zero_tol),
        label       = 'Xsuite',
        color       = "b",
        linestyle   = "--")
    axs[1].plot(
        zero_small_values(twiss_sad.s, tol = zero_tol),
        zero_small_values(twiss_sad.alfy, tol = zero_tol),
        color       = "r")
    axs[1].plot(
        zero_small_values(twiss_xsuite.s, tol = zero_tol),
        zero_small_values(twiss_xsuite.alfy, tol = zero_tol),
        color       = "b",
        linestyle   = "--")

    axs[0].legend()
    axs[0].set_ylabel(r'$\alpha_{x}$ [m]')
    axs[1].set_ylabel(r'$\alpha_{y}$ [m]')
    axs[1].set_xlabel('s [m]')

    if suptitle is not None:
        fig.suptitle(f"{suptitle}: " + "Alpha Functions ($\\alpha_{x}$, $\\alpha_{y}$)")
    else:
        fig.suptitle("Alpha Functions ($\\alpha_{x}$, $\\alpha_{y}$)")
    fig.tight_layout()
    fig.align_labels()
    fig.align_titles()

    
    ###################################### Dispersion #####################################
    
    fig, axs = plt.subplots(2, figsize = figsize, sharex = True)

    axs[0].plot(
        zero_small_values(twiss_sad.s, tol = zero_tol),
        zero_small_values(twiss_sad.dx, tol = zero_tol),
        label       = 'SAD',
        color       = "r")
    axs[0].plot(
        zero_small_values(twiss_xsuite.s, tol = zero_tol),
        zero_small_values(twiss_xsuite.dx, tol = zero_tol),
        label       = 'Xsuite',
        color       = "b",
        linestyle   = "--")
    axs[1].plot(
        zero_small_values(twiss_sad.s, tol = zero_tol),
        zero_small_values(twiss_sad.dy, tol = zero_tol),
        color       = "r")
    axs[1].plot(
        zero_small_values(twiss_xsuite.s, tol = zero_tol),
        zero_small_values(twiss_xsuite.dy, tol = zero_tol),
        color       = "b",
        linestyle   = "--")

    axs[0].legend()
    axs[0].set_ylabel(r'$D_{x}$ [m]')
    axs[1].set_ylabel(r'$D_{y}$ [m]')
    axs[1].set_xlabel('s [m]')

    if suptitle is not None:
        fig.suptitle(f"{suptitle}: " + "Dispersion ($D_{x}$, $D_{y}$)")
    else:
        fig.suptitle("Dispersion ($D_{x}$, $D_{y}$)")
    fig.tight_layout()
    fig.align_labels()
    fig.align_titles()

    
    ###################################### Derivative Dispersion #####################################
    
    fig, axs = plt.subplots(2, figsize = figsize, sharex = True)

    axs[0].plot(
        zero_small_values(twiss_sad.s, tol = zero_tol),
        zero_small_values(twiss_sad.dpx, tol = zero_tol),
        label       = 'SAD',
        color       = "r")
    axs[0].plot(
        zero_small_values(twiss_xsuite.s, tol = zero_tol),
        zero_small_values(twiss_xsuite.dpx, tol = zero_tol),
        label       = 'Xsuite',
        color       = "b",
        linestyle   = "--")
    axs[1].plot(
        zero_small_values(twiss_sad.s, tol = zero_tol),
        zero_small_values(twiss_sad.dpy, tol = zero_tol),
        color       = "r")
    axs[1].plot(
        zero_small_values(twiss_xsuite.s, tol = zero_tol),
        zero_small_values(twiss_xsuite.dpy, tol = zero_tol),
        color       = "b",
        linestyle   = "--")

    axs[0].legend()
    axs[0].set_ylabel(r'$D_{px}$ [m]')
    axs[1].set_ylabel(r'$D_{py}$ [m]')
    axs[1].set_xlabel('s [m]')

    if suptitle is not None:
        fig.suptitle(f"{suptitle}: " + "Dispersion ($D_{px}$, $D_{py}$)")
    else:
        fig.suptitle("Dispersion ($D_{px}$, $D_{py}$)")
    fig.tight_layout()
    fig.align_labels()
    fig.align_titles()
