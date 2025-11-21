"""
Check the tracked beam against the previously computed DA for BTE lattice.
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
from matplotlib.colors import LogNorm

from scipy.constants import c as clight

from tqdm import tqdm

################################################################################
# User Parameters
################################################################################
PARTICLE_FILEPATH   = "../data/yohsimotosan.dat"

########################################
# Particle parameters
########################################
BETA0           = 1.0
P0C_EV          = 7E9

GEMITT_Z_HER    = 1E-6 # TODO: Placeholder
BETS0_HER       = 7.950888821340036

########################################
# Define parameters
########################################
PHASE_CAVI_DEG  = 180.00
FREQ_CAVI_HZ    = 2856000000
VOLT_CAVI_V     = 100E6

########################################
# R56 Parameters
########################################
CHECK_POINTS    = ["pbte", "ld250.2", "l48d2", "lw4d2", "lw4d3", "injp"]
R56_VALUES      = [
    -0.7000001291106137, 3.260355458318869e-08, 3.403815856640212e-08,
    1.7202244367403647e-08, -5.991145028980211]

################################################################################
# Load Data
################################################################################

########################################
# Load from file
########################################
initial_distribution    = np.loadtxt(PARTICLE_FILEPATH)

########################################
# Convert to xt.Table
########################################
initial_distribution    = xt.Table({
    "name":  np.arange(initial_distribution[:, 0].shape[0]),
    "x":     initial_distribution[:, 0] - np.median(initial_distribution[:, 0]),
    "px":    initial_distribution[:, 1] - np.median(initial_distribution[:, 1]),
    "y":     initial_distribution[:, 2] - np.median(initial_distribution[:, 2]),
    "py":    initial_distribution[:, 3] - np.median(initial_distribution[:, 3]),
    "zeta":  initial_distribution[:, 4] - np.median(initial_distribution[:, 4]),
    "delta": initial_distribution[:, 5] - np.median(initial_distribution[:, 5]),
    "state": np.full(initial_distribution[:, 0].shape, 1)})

################################################################################
# "Track" up to cavity
################################################################################
first_cavi_zeta     = initial_distribution["zeta"] + \
    R56_VALUES[0] * initial_distribution["delta"]
first_cavi_delta    = initial_distribution["delta"]


fig, ax = plt.subplots(1)
ax.scatter(
    initial_distribution["zeta"] * 1E3,
    initial_distribution["delta"] * 1E2,
    s = 1)
fig.supxlabel(r"$\zeta$ [mm]")
fig.supylabel(r"$\delta$ [%]")
fig.suptitle(f"At first cavity entrance")
fig, ax = plt.subplots(1)
ax.scatter(
    first_cavi_zeta * 1E3,
    first_cavi_delta * 1E2,
    s = 1)
fig.supxlabel(r"$\zeta$ [mm]")
fig.supylabel(r"$\delta$ [%]")
fig.suptitle(f"At first cavity entrance")
plt.show()

################################################################################
# Make a function
################################################################################
def analytic_cavity_kicks(
        zeta_particles:     np.ndarray,
        delta_particles:    np.ndarray,
        p0c_ev:             float,
        phase_cavi_deg:     float,
        freq_cavi_hz:       float,
        volt_cavi_v:        float,
        n_cavities:         int,
        plot:               bool = False):
    
    delta_particles_0   = delta_particles.copy()
    
    for i in range(n_cavities):
        # Calculate the particle energy at the start of the cavity
        p0c_particles   = p0c_ev * (1 + delta_particles)

        # Calculate lag and kick seen by each particle
        lag_seen_by_particle        = 180 - 360.0 * (
            np.deg2rad(phase_cavi_deg) / (2 * np.pi) - \
            freq_cavi_hz * zeta_particles / BETA0 / clight)
        kick_seen_by_particle_ev    = volt_cavi_v * \
            np.sin(np.deg2rad(lag_seen_by_particle))
        
        # Update particle deltas
        p0c_particles           += kick_seen_by_particle_ev
        delta_particles         = p0c_particles / p0c_ev - 1.0

        # Update particle zetas: just drift to the next cavity
        # TODO: This falls down here, because this depends on orbit
        # With delta, there is additional orbit, so zeta changes

        if plot:
            # Calculate change in delta
            delta_change            = delta_particles - delta_particles_0

            # Plot
            fig, axs = plt.subplots(1, 3, figsize = (12, 8))
            
            axs[0].hist2d(
                zeta_particles * 1E3,
                delta_particles * 1E2,
                bins    = 101,
                cmap    = "viridis",
                norm    = LogNorm(vmin = 1, vmax = None))
            axs[0].set_facecolor("white")

            fig.colorbar(
                axs[0].collections[0],
                ax      = axs[0],
                label   = "Density")
            axs[0].set_xlabel(r"$\zeta$ [mm]")
            axs[0].set_ylabel(r"$\delta$ [%]")

            axs[1].scatter(
                zeta_particles * 1E3,
                delta_particles * 1E2,
                s       = 1,
                c       = delta_change * 1E2,
                cmap    = "plasma")

            fig.colorbar(
                axs[1].collections[0],
                ax      = axs[1],
                label   = "Change in Delta w.r.t. Input Beam [%]")
            axs[1].set_xlabel(r"$\zeta$ [mm]")
            axs[1].set_ylabel(r"$\delta$ [%]")

            axs[2].hist(delta_particles * 1E2, bins = 101)
            axs[2].set_xlabel(r"$\delta$ [%]")
            axs[2].set_ylabel("Counts")

            fig.suptitle(f"After cavity {i + 1}/{n_cavities}")
    
    return np.array(delta_particles)

analytic_cavity_kicks(
    zeta_particles  = first_cavi_zeta,
    delta_particles = first_cavi_delta,
    p0c_ev          = P0C_EV,
    phase_cavi_deg  = PHASE_CAVI_DEG,
    freq_cavi_hz    = FREQ_CAVI_HZ,
    volt_cavi_v     = VOLT_CAVI_V,
    n_cavities      = 4,
    plot            = True)
plt.show()

################################################################################
# Grid search
################################################################################
phase_test_points       = np.linspace(175.0, 185.0, 201)
volt_test_points        = np.linspace(0E6, 100E6, 401)
PHASE_GRID, VOLT_GRID   = np.meshgrid(
    phase_test_points,
    volt_test_points,
    indexing = "ij")

# The metric to optimize on is the mean longitudinal action amplitude
# Az = sqrt( zeta^2 / beta_z + delta^2 * beta_z / emitt_z )
MEAN_AZ_GRID    = np.zeros(PHASE_GRID.shape, dtype = float)

for i in tqdm(range(PHASE_GRID.shape[0])):
    for j in range(PHASE_GRID.shape[1]):
        phase_cavi_deg = PHASE_GRID[i, j]
        volt_cavi_v    = VOLT_GRID[i, j]
        
        final_delta = analytic_cavity_kicks(
            zeta_particles  = initial_distribution["zeta"],
            delta_particles = initial_distribution["delta"],
            p0c_ev          = P0C_EV,
            phase_cavi_deg  = phase_cavi_deg,
            freq_cavi_hz    = FREQ_CAVI_HZ,
            volt_cavi_v     = volt_cavi_v,
            n_cavities      = 4,
            plot            = False)
        final_zeta = first_cavi_zeta + \
            R56_VALUES[-1] * final_delta
        
        Az = np.sqrt(
            final_zeta** 2 / BETS0_HER + \
            final_delta**2 * BETS0_HER /
            GEMITT_Z_HER)

        MEAN_AZ_GRID[i, j]      = float(np.mean(Az))

minimum_sum     = np.min(MEAN_AZ_GRID)
argmin_sum      = np.unravel_index(np.argmin(MEAN_AZ_GRID), MEAN_AZ_GRID.shape)
optimum_lag     = PHASE_GRID[argmin_sum]
optimum_volt    = VOLT_GRID[argmin_sum]

fig, ax    = plt.subplots(1, figsize = (12, 8))
c = ax.pcolormesh(
    VOLT_GRID * 1E-6,
    PHASE_GRID,
    MEAN_AZ_GRID,
    shading = "auto",
    cmap    = "plasma",
    norm    = LogNorm(vmin = None, vmax = None))
fig.colorbar(c, ax = ax, label = "Mean of longitudinal Action Amplitude Az [m]")
ax.plot(
    VOLT_GRID[0, :] * 1E-6,
    PHASE_GRID[np.argmin(np.abs(MEAN_AZ_GRID), axis = 0), 0],
    color = "b",
    linestyle = "--",
    label = "Min Mean Az per Voltage")
ax.plot(
    VOLT_GRID[0, np.argmin(np.abs(MEAN_AZ_GRID), axis = 1)] * 1E-6,
    PHASE_GRID[:, 0],
    color       = "r",
    linestyle   = "--",
    label       = "Min Mean Az per Phase")
ax.scatter(
    optimum_volt * 1E-6,
    optimum_lag,
    marker  = "x",
    s       = 100)
ax.set_xlabel("Cavity Voltage [MV]")
ax.set_ylabel("Cavity Phase [deg]")

print(f"Optimum lag of:     {optimum_lag} deg")
print(f"Optimum voltage of: {optimum_volt * 1E-6} MV")

################################################################################
# Compare Nominal, optimal and implemented
################################################################################

########################################
# Nominal (18MV, 180 deg)
########################################
final_delta_nominal     = analytic_cavity_kicks(
    zeta_particles  = initial_distribution["zeta"],
    delta_particles = initial_distribution["delta"],
    p0c_ev          = P0C_EV,
    phase_cavi_deg  = 180,
    freq_cavi_hz    = FREQ_CAVI_HZ,
    volt_cavi_v     = 18E6,
    n_cavities      = 4,
    plot            = False)
final_zeta_nominal      = first_cavi_zeta + \
    R56_VALUES[-1] * final_delta_nominal
Az_nominal              = np.sqrt(
    final_zeta_nominal** 2 / BETS0_HER + \
    final_delta_nominal **2 * BETS0_HER /
    GEMITT_Z_HER)

########################################
# Optimal
########################################
final_delta_optimal     = analytic_cavity_kicks(
    zeta_particles  = initial_distribution["zeta"],
    delta_particles = initial_distribution["delta"],
    p0c_ev          = P0C_EV,
    phase_cavi_deg  = optimum_lag,
    freq_cavi_hz    = FREQ_CAVI_HZ,
    volt_cavi_v     = optimum_volt,
    n_cavities      = 4,
    plot            = False)
final_zeta_optimal      = first_cavi_zeta + \
    R56_VALUES[-1] * final_delta_optimal
Az_optimal              = np.sqrt(
    final_zeta_optimal** 2 / BETS0_HER + \
    final_delta_optimal **2 * BETS0_HER /
    GEMITT_Z_HER)

########################################
# Implemented (100MV, 180 deg)
########################################
final_delta_implemented = analytic_cavity_kicks(
    zeta_particles  = initial_distribution["zeta"],
    delta_particles = initial_distribution["delta"],
    p0c_ev          = P0C_EV,
    phase_cavi_deg  = 180,
    freq_cavi_hz    = FREQ_CAVI_HZ,
    volt_cavi_v     = 100E6,
    n_cavities      = 4,
    plot            = False)
final_zeta_implemented      = first_cavi_zeta + \
    R56_VALUES[-1] * final_delta_implemented
Az_implemented              = np.sqrt(
    final_zeta_implemented** 2 / BETS0_HER + \
    final_delta_implemented **2 * BETS0_HER /
    GEMITT_Z_HER)

################################################################################
# Plot Delta Comparison
################################################################################
fig, ax = plt.subplots(1)
bins    = np.linspace(-1, 1, 101)

ax.hist(
    final_delta_nominal * 1E2,
    bins    = bins,
    alpha   = 0.5,
    label   = "Nominal: Phase 180 deg, Volt 18 MV")
ax.hist(
    final_delta_implemented * 1E2,
    bins    = bins,
    alpha   = 0.5,
    label   = "Implemented: Phase 180 deg, Volt 100 MV")
ax.hist(
    final_delta_optimal * 1E2,
    bins    = bins,
    alpha   = 0.5,
    label   = f"Phase {optimum_lag} deg, Volt {optimum_volt * 1E-6} MV")

fig.supxlabel(r"$\delta$ [%]")
fig.supylabel("Counts")
fig.legend()

################################################################################
# Plot Action Comparison
################################################################################
fig, ax = plt.subplots(1)
bins    = np.linspace(0, 20, 101)

ax.hist(
    Az_nominal,
    bins    = bins,
    alpha   = 0.5,
    label   = "Nominal: Phase 180 deg, Volt 18 MV")
ax.hist(
    Az_implemented,
    bins    = bins,
    alpha   = 0.5,
    label   = "Implemented: Phase 180 deg, Volt 100 MV")
ax.hist(
    Az_optimal,
    bins    = bins,
    alpha   = 0.5,
    label   = f"Phase {optimum_lag} deg, Volt {optimum_volt * 1E-6} MV")

fig.supxlabel(r"$A_{z}$ []")
fig.supylabel("Counts")
fig.legend()

################################################################################
# Show plots
################################################################################
plt.show()