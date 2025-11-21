"""
Get the R-matrix elements from tracking in the ECS lattice
=============================================
Author(s):  John P T Salvesen, Nikita Z van Gils
Date:       14-11-2025
"""

################################################################################
# Required Packages
################################################################################
import xtrack as xt
import numpy as np
import matplotlib.pyplot as plt

################################################################################
# User Parameters
################################################################################

########################################
# Lattice Files
########################################
xsuite_lattice_path     = "../lattices/bte.py"
xsuite_optics_path      = "../lattices/bte_import_optics.py"

########################################
# Tracking Parameters
########################################
START_POINT             = "pbte"
END_POINT               = "injp"

################################################################################
# Load Lattice
################################################################################
env     = xt.Environment()
env.call(xsuite_lattice_path)
env.call(xsuite_optics_path)
line    = env.lines["line"]

################################################################################
# Calculate Slippage
################################################################################
def get_R_from_tracking(
        line:       xt.Line,
        ele_start:  str,
        ele_stop:   str,
        x0:         float               = 0.0,
        px0:        float               = 0.0,
        y0:         float               = 0.0,
        py0:        float               = 0.0,
        zeta0:      float               = 0.0,
        delta0:     float               = 0.0,
        dx:         float               = 1E-2,
        dpx:        float               = 1E-3,
        dy:         float               = 1E-2,
        dpy:        float               = 1E-3,
        dzeta:      float               = 1E-2,
        ddelta:     float               = 2E-2,
        n_test:     int                 = 101,
        fig:        plt.Figure | None   = None,
        axs:        np.ndarray | None   = None,
        label:      str | None          = None) -> float:
    """
    Estimate R56 = d(zeta_out)/d(delta_in) between ele_start and ele_end
    by central finite differences around delta0.
    """

    test_xs         = np.linspace(x0 - dx, x0 + dx, n_test)
    test_pxs        = np.linspace(px0 - dpx, px0 + dpx, n_test)
    test_ys         = np.linspace(y0 - dy, y0 + dy, n_test)
    test_pys        = np.linspace(py0 - dpy, py0 + dpy, n_test)
    test_zetas      = np.linspace(zeta0 - dzeta, zeta0 + dzeta, n_test)
    test_deltas     = np.linspace(delta0 - ddelta, delta0 + ddelta, n_test)
    
    x_particles     = line.build_particles(x = test_xs)
    px_particles    = line.build_particles(px = test_pxs)
    y_particles     = line.build_particles(y = test_ys)
    py_particles    = line.build_particles(py = test_pys)
    zeta_particles  = line.build_particles(zeta = test_zetas)
    delta_particles = line.build_particles(delta = test_deltas)
    
    line.track(x_particles, ele_start = ele_start, ele_stop = ele_stop)
    line.track(px_particles, ele_start = ele_start, ele_stop = ele_stop)
    line.track(y_particles, ele_start = ele_start, ele_stop = ele_stop)
    line.track(py_particles, ele_start = ele_start, ele_stop = ele_stop)
    line.track(zeta_particles, ele_start = ele_start, ele_stop = ele_stop)
    line.track(delta_particles, ele_start = ele_start, ele_stop = ele_stop)

    x_particles.sort(interleave_lost_particles = True)
    px_particles.sort(interleave_lost_particles = True)
    y_particles.sort(interleave_lost_particles = True)
    py_particles.sort(interleave_lost_particles = True)
    zeta_particles.sort(interleave_lost_particles = True)
    delta_particles.sort(interleave_lost_particles = True)

    if fig is None or axs is None:
        fig, axs    = plt.subplots(
            ncols   = 6,
            nrows   = 6,
            figsize = (16, 16),
            sharex  = "col",
            sharey  = "row")

    axs[0, 0].plot(
        test_xs[x_particles.state == 1] * 1E3,
        x_particles.x[x_particles.state == 1] * 1E3,
        label = label if label is not None else None)
    axs[1, 0].plot(
        test_xs[x_particles.state == 1] * 1E3,
        x_particles.px[x_particles.state == 1] * 1E3)
    axs[2, 0].plot(
        test_xs[x_particles.state == 1] * 1E3,
        x_particles.y[x_particles.state == 1] * 1E3)
    axs[3, 0].plot(
        test_xs[x_particles.state == 1] * 1E3,
        x_particles.py[x_particles.state == 1] * 1E3)
    axs[4, 0].plot(
        test_xs[x_particles.state == 1] * 1E3,
        x_particles.zeta[x_particles.state == 1] * 1E3)
    axs[5, 0].plot(
        test_xs[x_particles.state == 1] * 1E3,
        x_particles.delta[x_particles.state == 1] * 1E2)
    
    axs[0, 1].plot(
        test_pxs[px_particles.state == 1] * 1E3,
        px_particles.x[px_particles.state == 1] * 1E3)
    axs[1, 1].plot(
        test_pxs[px_particles.state == 1] * 1E3,
        px_particles.px[px_particles.state == 1] * 1E3)
    axs[2, 1].plot(
        test_pxs[px_particles.state == 1] * 1E3,
        px_particles.y[px_particles.state == 1] * 1E3)
    axs[3, 1].plot(
        test_pxs[px_particles.state == 1] * 1E3,
        px_particles.py[px_particles.state == 1] * 1E3)
    axs[4, 1].plot(
        test_pxs[px_particles.state == 1] * 1E3,
        px_particles.zeta[px_particles.state == 1] * 1E3)
    axs[5, 1].plot(
        test_pxs[px_particles.state == 1] * 1E3,
        px_particles.delta[px_particles.state == 1] * 1E2)
    
    axs[0, 2].plot(
        test_ys[y_particles.state == 1] * 1E3,
        y_particles.x[y_particles.state == 1] * 1E3)
    axs[1, 2].plot(
        test_ys[y_particles.state == 1] * 1E3,
        y_particles.px[y_particles.state == 1] * 1E3)
    axs[2, 2].plot(
        test_ys[y_particles.state == 1] * 1E3,
        y_particles.y[y_particles.state == 1] * 1E3)
    axs[3, 2].plot(
        test_ys[y_particles.state == 1] * 1E3,
        y_particles.py[y_particles.state == 1] * 1E3)
    axs[4, 2].plot(
        test_ys[y_particles.state == 1] * 1E3,
        y_particles.zeta[y_particles.state == 1] * 1E3)
    axs[5, 2].plot(
        test_ys[y_particles.state == 1] * 1E3,
        y_particles.delta[y_particles.state == 1] * 1E2)
    
    axs[0, 3].plot(
        test_pys[py_particles.state == 1] * 1E3,
        py_particles.x[py_particles.state == 1] * 1E3)
    axs[1, 3].plot(
        test_pys[py_particles.state == 1] * 1E3,
        py_particles.px[py_particles.state == 1] * 1E3)
    axs[2, 3].plot(
        test_pys[py_particles.state == 1] * 1E3,
        py_particles.y[py_particles.state == 1] * 1E3)
    axs[3, 3].plot(
        test_pys[py_particles.state == 1] * 1E3,
        py_particles.py[py_particles.state == 1] * 1E3)
    axs[4, 3].plot(
        test_pys[py_particles.state == 1] * 1E3,
        py_particles.zeta[py_particles.state == 1] * 1E3)
    axs[5, 3].plot(
        test_pys[py_particles.state == 1] * 1E3,
        py_particles.delta[py_particles.state == 1] * 1E2)
    
    axs[0, 4].plot(
        test_zetas[zeta_particles.state == 1] * 1E3,
        zeta_particles.x[zeta_particles.state == 1] * 1E3)
    axs[1, 4].plot(
        test_zetas[zeta_particles.state == 1] * 1E3,
        zeta_particles.px[zeta_particles.state == 1] * 1E3)
    axs[2, 4].plot(
        test_zetas[zeta_particles.state == 1] * 1E3,
        zeta_particles.y[zeta_particles.state == 1] * 1E3)
    axs[3, 4].plot(
        test_zetas[zeta_particles.state == 1] * 1E3,
        zeta_particles.py[zeta_particles.state == 1] * 1E3)
    axs[4, 4].plot(
        test_zetas[zeta_particles.state == 1] * 1E3,
        zeta_particles.zeta[zeta_particles.state == 1] * 1E3)
    axs[5, 4].plot(
        test_zetas[zeta_particles.state == 1] * 1E3,
        zeta_particles.delta[zeta_particles.state == 1] * 1E2)
    
    axs[0, 5].plot(
        test_deltas[delta_particles.state == 1] * 1E2,
        delta_particles.x[delta_particles.state == 1] * 1E3)
    axs[1, 5].plot(
        test_deltas[delta_particles.state == 1] * 1E2,
        delta_particles.px[delta_particles.state == 1] * 1E3)
    axs[2, 5].plot(
        test_deltas[delta_particles.state == 1] * 1E2,
        delta_particles.y[delta_particles.state == 1] * 1E3)
    axs[3, 5].plot(
        test_deltas[delta_particles.state == 1] * 1E2,
        delta_particles.py[delta_particles.state == 1] * 1E3)
    axs[4, 5].plot(
        test_deltas[delta_particles.state == 1] * 1E2,
        delta_particles.zeta[delta_particles.state == 1] * 1E3)
    axs[5, 5].plot(
        test_deltas[delta_particles.state == 1] * 1E2,
        delta_particles.delta[delta_particles.state == 1] * 1E2)
    
    axs[-1, 0].set_xlabel(r"$x_{i}$ [mm]")
    axs[-1, 1].set_xlabel(r"$p_{x, i}$ [mrad]")
    axs[-1, 2].set_xlabel(r"$y_{i}$ [mm]")
    axs[-1, 3].set_xlabel(r"$p_{y, i}$ [mrad]")
    axs[-1, 4].set_xlabel(r"$\zeta_{i}$ [mm]")
    axs[-1, 5].set_xlabel(r"$\delta_{i}$ [%]")

    axs[0, 0].set_ylabel(r"$x_{f}$ [mm]")
    axs[1, 0].set_ylabel(r"$p_{x, f}$ [mrad]")
    axs[2, 0].set_ylabel(r"$y_{f}$ [mm]")
    axs[3, 0].set_ylabel(r"$p_{y, f}$ [mrad]")
    axs[4, 0].set_ylabel(r"$\zeta_{f}$ [mm]")
    axs[5, 0].set_ylabel(r"$\delta_{f}$ [%]")

    fig.suptitle(f"R between {ele_start} and {ele_stop}")
    
    fig.legend()

    fig.tight_layout()
    fig.align_labels()
    fig.align_titles()

################################################################################
# R56
################################################################################
fig, axs    = plt.subplots(
    ncols   = 6,
    nrows   = 6,
    figsize = (16, 16),
    sharex  = "col",
    sharey  = "row")

for voltage in [0, 20E6, 40E6, 60E6, 80E6, 100E6]:
    line["volt_acw"]    = voltage
    get_R_from_tracking(
        line        = line,
        ele_start   = START_POINT,
        ele_stop    = END_POINT,
        fig         = fig,
        axs         = axs,
        label       = f"V = {voltage / 1E6}MV")

################################################################################
# Show plots
################################################################################
plt.show()
