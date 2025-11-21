"""
Track the KEK BTe lattice
=============================================
Author(s):  John P T Salvesen, Nikita Z van Gils
Date:       15-10-2025
"""

#################################### Required Packages ###################################

import xtrack as xt
import numpy as np
import matplotlib.pyplot as plt

#################################### User Parameters ####################################

######################################### Lattice Files ########################################

xsuite_lattice_path     = "../lattices/bte.py"
xsuite_optics_path      = "../lattices/bte_import_optics.py"

######################################### Emittances ########################################

NEMITT_X    = 80E-6
NEMITT_Y    = 20E-6
SIGMA_Z     = 6E-3
SIGMA_DELTA = 0.0007

#################################### Initial Conditions ########################################

betx_init   = 24.097096226706224
bety_init   = 17.11998683241826
alfx_init   = .19232400480757933
alfy_init   = -2.131413440451613

##################################### Tracking Parameters ########################################

N_PART      = int(1E5)
ELE_STOP    = "injp"

#################################### Load Converted Lattice ############################

env     = xt.Environment()
env.call(xsuite_lattice_path)
env.call(xsuite_optics_path)
line    = env.lines["line"]

######################################### Twiss ########################################

tw = line.twiss(
    start   = xt.START,
    end     = xt.END,
    betx    = betx_init,
    bety    = bety_init,
    alfx    = alfx_init,
    alfy    = alfy_init)


######################################### Build Particles ########################################

######################################### Calculate transverse matrix ########################################

TRANSVERSE_MATRIX = np.array([
    [np.sqrt(betx_init), 0, 0, 0],
    [alfx_init / np.sqrt(betx_init), 1 / np.sqrt(betx_init), 0, 0],
    [0, 0, np.sqrt(bety_init), 0],
    [0, 0, alfy_init / np.sqrt(bety_init), 1 / np.sqrt(bety_init)]])

######################################### Caclulate Geometric Emittance ########################################

GEMITT_X  = NEMITT_X / line.particle_ref.gamma0
GEMITT_Y  = NEMITT_Y / line.particle_ref.gamma0

######################################### Create Initial distribution ########################################

x_init, px_init, y_init, py_init = TRANSVERSE_MATRIX @ (
    np.sqrt(np.abs([GEMITT_X, GEMITT_X, GEMITT_Y, GEMITT_Y])) * \
    np.random.normal(0, 1, (4, N_PART)))

zeta_init   = np.random.normal(0, SIGMA_Z, N_PART)
delta_init  = np.random.normal(0, SIGMA_DELTA, N_PART)

######################################### Generate bunch ########################################

tracked_beam    = xt.Particles(
    p0c             = line.particle_ref.p0c,
    q0              = line.particle_ref.q0,
    mass0           = line.particle_ref.mass0,
    x               = x_init,
    px              = px_init,
    y               = y_init,
    py              = py_init,
    zeta            = zeta_init,
    delta           = delta_init)
initial_distribution = tracked_beam.copy()

######################################### Track ########################################

line.track(particles = tracked_beam, ele_stop  = ELE_STOP)

#################################### Setup Tables ###################################

initial = xt.Table({
    "name":     np.arange(initial_distribution.x.shape[0]),
    "x":        initial_distribution.x,
    "px":       initial_distribution.px,
    "y":        initial_distribution.y,
    "py":       initial_distribution.py,
    "zeta":     initial_distribution.zeta,
    "delta":    initial_distribution.delta,
    "state":    initial_distribution.state})
final   = xt.Table({
    "name":     np.arange(tracked_beam.x.shape[0]),
    "x":        tracked_beam.x,
    "px":       tracked_beam.px,
    "y":        tracked_beam.y,
    "py":       tracked_beam.py,
    "zeta":     tracked_beam.zeta,
    "delta":    tracked_beam.delta,
    "state":    tracked_beam.state})

#################################### Compare ######################################

fig, axs = plt.subplots(2, 3, figsize = (10, 8))

axs[0, 0].hist2d(
    initial.x[initial.state == 1],
    initial.px[initial.state == 1],
    bins        = 101,
    cmap        = 'plasma')
axs[0, 1].hist2d(
    initial.y[initial.state == 1],
    initial.py[initial.state == 1],
    bins        = 101,
    cmap        = 'plasma')
axs[0, 2].hist2d(
    initial.zeta[initial.state == 1], # type: ignore
    initial.delta[initial.state == 1],
    bins        = 101,
    cmap        = 'plasma')
axs[1, 0].hist2d(
    final.x[final.state == 1],
    final.px[final.state == 1],
    bins        = 101,
    cmap        = 'plasma')
axs[1, 1].hist2d(
    final.y[final.state == 1],
    final.py[final.state == 1],
    bins        = 101,
    cmap        = 'plasma')
axs[1, 2].hist2d(
    final.zeta[final.state == 1], # type: ignore
    final.delta[final.state == 1],
    bins        = 101,
    cmap        = 'plasma')

axs[0, 0].set_xlabel("x [m]")
axs[0, 0].set_ylabel("px [rad]")
axs[0, 1].set_xlabel("y [m]")
axs[0, 1].set_ylabel("py [rad]")
axs[0, 2].set_xlabel("zeta [m]")
axs[0, 2].set_ylabel("delta [m]")
axs[1, 0].set_xlabel("x [m]")
axs[1, 0].set_ylabel("px [rad]")
axs[1, 1].set_xlabel("y [m]")
axs[1, 1].set_ylabel("py [rad]")
axs[1, 2].set_xlabel("zeta [m]")
axs[1, 2].set_ylabel("delta [m]")

axs[0, 0].set_title("Initial")
axs[0, 1].set_title("Initial")
axs[0, 2].set_title("Initial")
axs[1, 0].set_title("Final")
axs[1, 1].set_title("Final")
axs[1, 2].set_title("Final")

fig.align_labels()
fig.tight_layout()
fig.savefig("../out/003_bte_track_initial_vs_final.png", dpi = 300)

#################################### Show plots ###################################

plt.show()
