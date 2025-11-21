"""
Track the KEK BTp lattice
=============================================
Author(s):  John P T Salvesen, Nikita Z van Gils
Date:       31-10-2025
"""

#################################### Required Packages ###################################

import xtrack as xt
import numpy as np
import matplotlib.pyplot as plt

#################################### User Parameters ####################################

######################################### Lattice Files ########################################

# Load lattice
line = xt.Line.from_json("../lattices/btp.json")
env = line.env
# Apply the optics
env.call("../lattices/btp_import_optics.py")


######################################### Particles Filepath ########################################

particles_filepath      = "../data/positrons_wake_10k.npy"

#################################### Initial Conditions ########################################

betx_init   = 20.048570292837496
bety_init   = 41.21183267623789
alfx_init   = -1.3077805745635371
alfy_init   = -2.131413440451613

##################################### Tracking Parameters ########################################


ELE_START    = "pbtp"
ELE_STOP    = "injp"
line.particle_ref.p0c = 3.98E9

######################################### Twiss ########################################

tw = line.twiss(
    start   = xt.START,
    end     = xt.END,
    betx    = betx_init,
    bety    = bety_init,
    alfx    = alfx_init,
    alfy    = alfy_init)


######################################### Build Particles ########################################


##########################  Load initial distribution ############################################
initial_distribution = np.load(particles_filepath, allow_pickle=True)
# #initial_distribution = np.loadtxt(particles_filepath)
# # N.B. Zeta coordinate in Ocelot is defined with opposite sign to Xsuite
initial = xt.Table({
    "name":  np.arange(initial_distribution[0].shape[0]),
    "x":     +1 * initial_distribution[0] - np.median(initial_distribution[:, 0]),
    "px":    +1 * initial_distribution[1]- np.median(initial_distribution[:, 1]),
    "y":     +1 * initial_distribution[2]- np.median(initial_distribution[:, 2]),
    "py":    +1 * initial_distribution[3]- np.median(initial_distribution[:, 3]),
    "zeta":  -1 * initial_distribution[4]- np.median(initial_distribution[4]),
    "delta": +1 * initial_distribution[5] - np.median(initial_distribution[5]),
    "state": np.full(initial_distribution[0].shape, 1)
})

######################################### Generate bunch ########################################

tracked_beam    = xt.Particles(
    p0c             = line.particle_ref.p0c,
    q0              = line.particle_ref.q0,
    mass0           = line.particle_ref.mass0,
    x               = initial.x,
    px              = initial.px,
    y               = initial.y,
    py              = initial.py,
    zeta            = initial.zeta,
    delta           = initial.delta)
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
#fig.savefig("../out/004_btp_track_initial_vs_final.png", dpi = 300)

#################################### Show plots ###################################

plt.show()
