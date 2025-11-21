"""
Compare tracking of BTp with ECS between SAD and Xsuite
=============================================
Author(s):  John P T Salvesen, Nikita Z van Gils
Date:       14-10-2025
"""


################################### Required Packages ##################################

import xtrack as xt
import numpy as np
import matplotlib.pyplot as plt


################################### User Parameters ###################################

initial_dist_file       = "../out/btp_pstart.dat"
final_dist_file         = "../out/btp_pinjax0.dat"
initial_dist_element    = "PSTART"
final_dist_element      = "PINJAX0"
xsuite_lattice_path     = "../lattices/btp.py"
xsuite_optics_path      = "../lattices/btp_import_optics.py"


############################# Load Distributions from SAD ############################
initial_sad = np.loadtxt(initial_dist_file)
final_sad   = np.loadtxt(final_dist_file)

initial_sad = xt.Table({
    "name":     np.arange(initial_sad.shape[0]),
    "x":        initial_sad[:, 0],
    "px":       initial_sad[:, 1],
    "y":        initial_sad[:, 2],
    "py":       initial_sad[:, 3],
    "zeta":     initial_sad[:, 4],
    "delta":    initial_sad[:, 5],
    "state":    initial_sad[:, 6]})
final_sad   = xt.Table({
    "name":     np.arange(final_sad.shape[0]),
    "x":        final_sad[:, 0],
    "px":       final_sad[:, 1],
    "y":        final_sad[:, 2],
    "py":       final_sad[:, 3],
    "zeta":     final_sad[:, 4],
    "delta":    final_sad[:, 5],
    "state":    final_sad[:, 6]})

################################### Load Converted Lattice ##################################


line = xt.Line.from_json("../lattices/btp.json")
env = line.env
################################### Build Particles #########################################



################################### Generate matched bunch ##################################

tracked_beam    = xt.Particles(
    p0c             = line.particle_ref.p0c,
    q0              = line.particle_ref.q0,
    mass0           = line.particle_ref.mass0,
    x               = initial_sad.x,
    px              = initial_sad.px,
    y               = initial_sad.y,
    py              = initial_sad.py,
    zeta            = initial_sad.zeta,
    delta           = initial_sad.delta)
initial_distribution = tracked_beam.copy()


################################### Track ##################################################

line.track(
    particles = tracked_beam,
    ele_start = initial_dist_element.lower(),
    ele_stop  = final_dist_element.lower())


#################################### Setup Tables ##########################################

initial_xs  = xt.Table({
    "name":     np.arange(initial_distribution.x.shape[0]),
    "x":        initial_distribution.x,
    "px":       initial_distribution.px,
    "y":        initial_distribution.y,
    "py":       initial_distribution.py,
    "zeta":     initial_distribution.zeta,
    "delta":    initial_distribution.delta,
    "state":    initial_distribution.state})
final_xs    = xt.Table({
    "name":     np.arange(tracked_beam.x.shape[0]),
    "x":        tracked_beam.x,
    "px":       tracked_beam.px,
    "y":        tracked_beam.y,
    "py":       tracked_beam.py,
    "zeta":     tracked_beam.zeta,
    "delta":    tracked_beam.delta,
    "state":    tracked_beam.state})


#################################### Compare ###################################



#################################### Initial ###################################

fig, axs = plt.subplots(2, 3, figsize = (10, 8), sharex = 'col', sharey = 'col')

axs[0, 0].hist2d(
    initial_sad.x[initial_sad.state == 1],
    initial_sad.px[initial_sad.state == 1],
    bins        = 101,
    cmap        = 'plasma')
axs[0, 1].hist2d(
    initial_sad.y[initial_sad.state == 1],
    initial_sad.py[initial_sad.state == 1],
    bins        = 101,
    cmap        = 'plasma')
axs[0, 2].hist2d(
    initial_sad.zeta[initial_sad.state == 1], # type: ignore
    initial_sad.delta[initial_sad.state == 1],
    bins        = 101,
    cmap        = 'plasma')
axs[1, 0].hist2d(
    initial_xs.x[initial_xs.state == 1],
    initial_xs.px[initial_xs.state == 1],
    bins        = 101,
    cmap        = 'plasma')
axs[1, 1].hist2d(
    initial_xs.y[initial_xs.state == 1],
    initial_xs.py[initial_xs.state == 1],
    bins        = 101,
    cmap        = 'plasma')
axs[1, 2].hist2d(
    initial_xs.zeta[initial_xs.state == 1], # type: ignore
    initial_xs.delta[initial_xs.state == 1],
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

axs[0, 0].set_title("SAD")
axs[0, 1].set_title("SAD")
axs[0, 2].set_title("SAD")
axs[1, 0].set_title("Xsuite")
axs[1, 1].set_title("Xsuite")
axs[1, 2].set_title("Xsuite")

fig.suptitle("Initial Distribution at " + initial_dist_element)
fig.align_labels()
fig.tight_layout()
fig.savefig("../out/003_btp_track_sad_vs_xs_initial.png", dpi = 300)


#################################### Final ###################################

fig, axs = plt.subplots(2, 3, figsize = (10, 8), sharex = 'col', sharey = 'col')

axs[0, 0].hist2d(
    final_sad.x[final_sad.state == 1],
    final_sad.px[final_sad.state == 1],
    bins        = 101,
    cmap        = 'plasma')
axs[0, 1].hist2d(
    final_sad.y[final_sad.state == 1],
    final_sad.py[final_sad.state == 1],
    bins        = 101,
    cmap        = 'plasma')
axs[0, 2].hist2d(
    final_sad.zeta[final_sad.state == 1], # type: ignore
    final_sad.delta[final_sad.state == 1],
    bins        = 101,
    cmap        = 'plasma')
axs[1, 0].hist2d(
    final_xs.x[final_xs.state == 1],
    final_xs.px[final_xs.state == 1],
    bins        = 101,
    cmap        = 'plasma')
axs[1, 1].hist2d(
    final_xs.y[final_xs.state == 1],
    final_xs.py[final_xs.state == 1],
    bins        = 101,
    cmap        = 'plasma')
axs[1, 2].hist2d(
    final_xs.zeta[final_xs.state == 1], # type: ignore
    final_xs.delta[final_xs.state == 1],
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

axs[0, 0].set_title("SAD")
axs[0, 1].set_title("SAD")
axs[0, 2].set_title("SAD")
axs[1, 0].set_title("Xsuite")
axs[1, 1].set_title("Xsuite")
axs[1, 2].set_title("Xsuite")

fig.suptitle("Final Distribution at " + final_dist_element)
fig.align_labels()
fig.tight_layout()
fig.savefig("../out/003_btp_track_sad_vs_xs_final.png", dpi = 300)


#################################### Surviving Particles ###################################

fig, axs = plt.subplots(2, 3, figsize = (10, 8), sharex = 'col', sharey = 'col')

axs[0, 0].scatter(
    initial_sad.x[final_sad.state == 1],
    initial_sad.px[final_sad.state == 1],
    s           = 1,
    color       = 'b')
axs[0, 1].scatter(
    initial_sad.y[final_sad.state == 1],
    initial_sad.py[final_sad.state == 1],
    s           = 1,
    color       = 'b')
axs[0, 2].scatter(
    initial_sad.zeta[final_sad.state == 1],
    initial_sad.delta[final_sad.state == 1],
    s           = 1,
    color       = 'b')
axs[1, 0].scatter(
    initial_xs.x[final_xs.state == 1],
    initial_xs.px[final_xs.state == 1],
    s           = 1,
    color       = 'b')
axs[1, 1].scatter(
    initial_xs.y[final_xs.state == 1],
    initial_xs.py[final_xs.state == 1],
    s           = 1,
    color       = 'b')
axs[1, 2].scatter(
    initial_xs.zeta[final_xs.state == 1],
    initial_xs.delta[final_xs.state == 1],
    s           = 1,
    color       = 'b')

axs[0, 0].scatter(
    initial_sad.x[final_sad.state != 1],
    initial_sad.px[final_sad.state != 1],
    s           = 1,
    color       = 'r')
axs[0, 1].scatter(
    initial_sad.y[final_sad.state != 1],
    initial_sad.py[final_sad.state != 1],
    s           = 1,
    color       = 'r')
axs[0, 2].scatter(
    initial_sad.zeta[final_sad.state != 1],
    initial_sad.delta[final_sad.state != 1],
    s           = 1,
    color       = 'r')
axs[1, 0].scatter(
    initial_xs.x[final_xs.state != 1],
    initial_xs.px[final_xs.state != 1],
    s           = 1,
    color       = 'r')
axs[1, 1].scatter(
    initial_xs.y[final_xs.state != 1],
    initial_xs.py[final_xs.state != 1],
    s           = 1,
    color       = 'r')
axs[1, 2].scatter(
    initial_xs.zeta[final_xs.state != 1],
    initial_xs.delta[final_xs.state != 1],
    s           = 1,
    color       = 'r')

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

axs[0, 0].set_title("SAD")
axs[0, 1].set_title("SAD")
axs[0, 2].set_title("SAD")
axs[1, 0].set_title("Xsuite")
axs[1, 1].set_title("Xsuite")
axs[1, 2].set_title("Xsuite")

fig.suptitle("Surviving Initial Particles at " + initial_dist_element)
fig.align_labels()
fig.tight_layout()
fig.savefig("../out/003_btp_track_sad_vs_xs_surviving.png", dpi = 300)


#################################### Show plots ###################################

plt.show()
