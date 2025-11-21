"""
Plot phase space at a chosen lag, voltage for the ECS at 
before/after each cavity along the ECS
User choices include radiation ON/OFF Tapering ON/OFF and energy of line
with respect to HER energy
=============================================
Author(s):  John P T Salvesen, Nikita Z van Gils
Date:       18-11-2025
"""

######################################### Required Packages ########################################

import os
import numpy as np
from matplotlib.colors import LogNorm
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
from helpers._compute_emittance import compute_emittance_and_optics
import xtrack as xt
import xobjects as xo
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from helpers._calculate_action_amplitude import calculate_action_amplitude
from helpers._fit_ellipse_da import calculate_boundary_points
from helpers._check_da_frontier import compare_to_frontier
from helpers._local_momentum_taper import taper_magnets_to_local_momentum
from helpers._compute_emittance import compute_emittance_and_optics

################################ User parameters ####################################

#Line configuration:
Radiation = True 	                 # TRUE ENABLES RADIATION 
Tapering = True 	                 # TRUE ENABLES TAPERING 
aperture_check = True                # TRUE ENABLES APERTURE CHECKS DURING TRACKING
#ECS settings:
ECS_voltage = None	                 # SET A VALUE IN MV TO ENABLE ELSE DEFAULT IS 18MV
ECS_lag= None	                     # SET A LAG IN DEGREES TO ENABLE ELSE DEFAULT IS 180° ZERO CROSSING IN XSUITE
#Injection type
#IF synchrotron injection 
SYNCHROTRON_OFFSET_BT = 0.006       # IN % WITH RESPECT TO  HER ENERGY. SET TO 0 FOR BETATRON INJECTION

ELE_STOP    = "injp"                 # FINAL ELEMENT TO TRACK TO
Plot_type   = "hist"                # CHOOSE WHICH TYPE OF PLOT TO SHOW PHASE SPACE ["scatter", "hist"]


# ---------------------------
# Load lattice and prepare initial distribution
# ---------------------------
xsuite_lattice_path = "../lattices/bte.py"
xsuite_optics_path  = "../lattices/bte_import_optics.py"



# Initial Twiss parameters at the start of the lattice
betx_init               = 20.540437095810788 #24.097096226706224
bety_init               = 18.015589974012826 #17.11998683241826
alfx_init               = 1.0030653048224993 #.19232400480757933
alfy_init               = .37003884861294073 #-2.131413440451613

#################################### Synchrotron Injection Energy Offset #######################

SYNCHROTRON_INJ_E_SCALE = 1.000 + SYNCHROTRON_OFFSET_BT          # Relative energy so that increase at injection is 1.006 
REFERENCE_ENERGY_HER    = 7.00E9                                 # [GeV], reference energy of HER

##################################### Final Element ############################################

ELE_STOP                = "injp"

#################################### Injection Parameters ######################################

X_SEPTUM_CUT            = 35.8E-3 * -1           # [m],      Septum edge position, w.r.t. HER design orbit (INJECTIO)
X_INJECTION_POINT       = 36.9E-3 * -1           # [m],      INJP x offset, w.r.t. HER design orbit (INJECTIO)
PX_INJECTION_POINT      = 2.5E-3 * -1           # [rad],    INJP px offset, w.r.t. HER design orbit (INJECTIO)

##################################### Aperture Parameters ######################################

# T. Mori 
#   Round cross-section with 47 mm diameter at drifts and Q-magnets
#   Rectangle with about 30 mm height and 50-100 mm width at the dipoles
RADIUS_BP_APERTURE      = 0.0235           # [m],      Beam pipe radius at drifts and Q-magnets
HEIGHT_BEND_APERTURE    = 0.03             # [m],      Bend aperture height
WIDTH_BEND_APERTURE     = 0.05             # [m],      Bend aperture width

N_THREADS               = 10


######################################### Load Lattice ###########################################

######################################### Load from Environment ##################################

env     = xt.Environment()
env.call(xsuite_lattice_path)
env.call(xsuite_optics_path)
line    = env.lines["line"]
if ECS_voltage is not None:
    line["volt_acw"] = ECS_voltage * 1E6
else :
    line["volt_acw"] = 18E6
if ECS_lag is not None:
    line["acw_lag"]  = ECS_lag
else :
    line["acw_lag"] = 180


######################################### Adjust Reference Energy #################################

line.particle_ref.p0c[0]   *= SYNCHROTRON_INJ_E_SCALE
######################################### Install Apertures #######################################

#################################### Bounding Aperture Helper Function ############################

def add_bounding_apertures(name, idx, aperture, aper_dict, idx_aper):
    idx_aper.append(idx-1)
    aper_dict[name + '_aper_entry'] = aperture.copy()
    idx_aper.append(idx)
    aper_dict[name + '_aper_exit'] = aperture.copy()
    
######################################### Insert Apertures Function ################################

def insert_apertures(line, idx_aper, aper_dict):
    shift = 0
    for idx, (aper_name, aper) in zip(idx_aper, aper_dict.items()):
        adjusted_idx = idx + shift  # Adjust index to account for prior insertions
        line.insert_element(
            at      = adjusted_idx + 1,
            name    = aper_name,
            element = aper)
        shift += 1  # Each insertion shifts future indices by 1

######################################### Generate the apertures ##################################

bp_aperture     = xt.LimitEllipse(
    a       = RADIUS_BP_APERTURE,
    b       = RADIUS_BP_APERTURE)
bend_aperture   = xt.LimitRect(
    min_x   = -1 * WIDTH_BEND_APERTURE / 2,
    max_x   = +1 * WIDTH_BEND_APERTURE / 2,
    min_y   = -1 * HEIGHT_BEND_APERTURE / 2,
    max_y   = +1 * HEIGHT_BEND_APERTURE / 2)

######################################### Get a clean table ######################################

tt  = line.get_table(attr = True)


###################### Get a list of which types of elements need aperture ########################

tt_bend     = tt.rows[
    (tt.element_type == "Bend") |
    (tt.element_type == "ThickSliceBend") |
    (tt.element_type == "RBend") |
    (tt.element_type == "ThickSliceRBend")]

tt_active   = tt.rows[
    (tt.element_type != 'Drift') & 
    (tt.element_type != 'Marker') &
    (tt.element_type != 'Bend') &
    (tt.element_type != "ThickSliceBend") &
    (tt.element_type != 'RBend') &
    (tt.element_type != "ThickSliceRBend") &
    (tt.element_type != '')]


needs_bp_aperture   = np.unique(tt_active.element_type)
needs_bend_aperture = np.unique(tt_bend.element_type)


######################################## Calculate aperture positions  ############################

idx_aper    = []
aper_dict   = {}
for idx, (ele_name, ele_type) in enumerate(zip(tt.name, tt.element_type)):
    if ele_type in needs_bp_aperture:
        add_bounding_apertures(
            name        = ele_name,
            idx         = idx,
            aperture    = bp_aperture,
            aper_dict   = aper_dict,
            idx_aper    = idx_aper)
    if ele_type in needs_bend_aperture:
        add_bounding_apertures(
            name        = ele_name,
            idx         = idx,
            aperture    = bend_aperture,
            aper_dict   = aper_dict,
            idx_aper    = idx_aper)

####################### Install the apertures ######################################################
if aperture_check == True:
    insert_apertures(line, idx_aper, aper_dict)
    aper_check = line.check_aperture()
else:
    pass

#particles_filepath   = " ../data/input_electrons_wake.npy"
particles_filepath   = "../data/yohsimotosan.dat"

#########  Load initial distribution ############################################
#initial_distribution = np.load("../data/input_electrons_wake.npy", allow_pickle=True)
initial_distribution = np.loadtxt(particles_filepath)
# N.B. Zeta coordinate in Ocelot is defined with opposite sign to Xsuite
# initial = xt.Table({
#     "name":  np.arange(initial_distribution[0].shape[0]),
#     "x":     +1 * initial_distribution[0] - np.median(initial_distribution[:, 0]),
#     "px":    +1 * initial_distribution[1]- np.median(initial_distribution[:, 1]),
#     "y":     +1 * initial_distribution[2]- np.median(initial_distribution[:, 2]),
#     "py":    +1 * initial_distribution[3]- np.median(initial_distribution[:, 3]),
#     "zeta":  -1 * initial_distribution[4]- np.mean(initial_distribution[4]),
#     "delta": +1 * initial_distribution[5] - np.mean(initial_distribution[5]),
#     "state": np.full(initial_distribution[0].shape, 1)
# })


initial = xt.Table({
    "name":  np.arange(initial_distribution[:, 0].shape[0]),
    "x":     initial_distribution[:, 0] - np.median(initial_distribution[:, 0]),
    "px":    initial_distribution[:, 1] - np.median(initial_distribution[:, 1]),
    "y":     initial_distribution[:, 2] - np.median(initial_distribution[:, 2]),
    "py":    initial_distribution[:, 3] - np.median(initial_distribution[:, 3]),
    "zeta":  initial_distribution[:, 4] - np.median(initial_distribution[:, 4]),
    "delta": initial_distribution[:, 5] - np.median(initial_distribution[:, 5]),
    "state": np.full(initial_distribution[:, 0].shape, 1)
})
master_particles = xt.Particles(
    p0c   = line.particle_ref.p0c,
    q0    = line.particle_ref.q0,
    mass0 = line.particle_ref.mass0,
    x     = initial.x,
    px    = initial.px,
    y     = initial.y,
    py    = initial.py,
    zeta  = initial.zeta,
    delta = initial.delta
)

####################### PLOT TRANSVERSE INPUT DISTRIBUTION #####################################

# # Get particle coordinates
# x = master_particles.x
# y = master_particles.y
# z = master_particles.zeta

# # Compute RMS sizes
# rms_x = np.sqrt(np.mean(x**2))
# rms_y = np.sqrt(np.mean(y**2))
# rms_z = np.sqrt(np.mean(z**2))

# print(f"RMS x: {rms_x*1e3:.3f} mm")
# print(f"RMS y: {rms_y*1e3:.3f} mm")
# print(f"RMS z: {rms_z*1e3:.3f} mm")

# # --- Set up figure ---
# fig, axs = plt.subplots(1, 3, figsize=(18, 5))

# # x-y plane
# axs[0].scatter(x*1e3, y*1e3, s=5, alpha=0.7)
# axs[0].set_xlabel("x [mm]")
# axs[0].set_ylabel("y [mm]")
# axs[0].set_title(f"x-y projection\nRMS x = {rms_x*1e3:.3f} mm, RMS y = {rms_y*1e3:.3f} mm")
# axs[0].grid(True, linestyle='--', alpha=0.3)

# # x-z plane
# axs[1].scatter(x*1e3, z*1e3, s=5, alpha=0.7)
# axs[1].set_xlabel("x [mm]")
# axs[1].set_ylabel("z [mm]")
# axs[1].set_title(f"x-z projection\nRMS x = {rms_x*1e3:.3f} mm, RMS z = {rms_z*1e3:.3f} mm")
# axs[1].grid(True, linestyle='--', alpha=0.3)

# # y-z plane
# axs[2].scatter(y*1e3, z*1e3, s=5, alpha=0.7)
# axs[2].set_xlabel("y [mm]")
# axs[2].set_ylabel("z [mm]")
# axs[2].set_title(f"y-z projection\nRMS y = {rms_y*1e3:.3f} mm, RMS z = {rms_z*1e3:.3f} mm")
# axs[2].grid(True, linestyle='--', alpha=0.3)

# plt.tight_layout()
# plt.show()

####################### PLOT PHASE SPACES OF INITIAL DISTRIBUTION #####################################
# Only consider particles with state==1 (all in this case)
mask = initial.state == 1

x  = initial.x[mask]
px = initial.px[mask]
y  = initial.y[mask]
py = initial.py[mask]
zeta  = initial.zeta[mask]
delta = initial.delta[mask]

# --- Plot transverse and longitudinal phase spaces ---
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# x vs px
axes[0].scatter(x*1e3, px*1e3, s=5, alpha=0.5, color='C0')
axes[0].set_xlabel("x [mm]")
axes[0].set_ylabel("px [mrad]")
axes[0].set_title("Transverse Phase Space x-px")
axes[0].grid(True)

# y vs py
axes[1].scatter(y*1e3, py*1e3, s=5, alpha=0.5, color='C1')
axes[1].set_xlabel("y [mm]")
axes[1].set_ylabel("py [mrad]")
axes[1].set_title("Transverse Phase Space y-py")
axes[1].grid(True)

# zeta vs delta
axes[2].scatter(zeta*1e3, delta*1e3, s=5, alpha=0.5, color='C2')
axes[2].set_xlabel("zeta [mm]")
axes[2].set_ylabel("delta [relative]")
axes[2].set_title("Longitudinal Phase Space zeta-delta")
axes[2].grid(True)


plt.suptitle("Initial Particle Distribution (10k Particles)", fontsize=16)
plt.tight_layout()
plt.show()


for i in range(4):
    #line[f'acw.{i}'].lag = optimal_lag
    if ECS_lag is not None:
        line[f'acw.{i}'].lag  = ECS_lag
    else :
        line[f'acw.{i}'].lag = 180
    
    
# Observation points
observation_points = [ "acw.0", "acw.1", "acw.2", "acw.3"]
labels = [  "Before ACW.0", "Before ACW.1", "Before ACW.2", "Before ACW.3"]

# Define phase-space combinations
combinations = [
    ("x", "y", "X-Y"),
    ("x", "zeta", "X-Zeta"),
    ("x", "delta", "X-Delta"),
    ("y", "zeta", "Y-Zeta"),
    ("y", "delta", "Y-Delta"),
    ("zeta", "delta", "Zeta-Delta")
]

nrows = len(observation_points)
ncols = len(combinations)

if Plot_type == "scatter": 
    # ===========================================================
    #  FIGURE 1: SCATTER PLOTS
    # ===========================================================
    nrows = len(observation_points)
    ncols = len(combinations)

    fig1, axs1 = plt.subplots(
        nrows, ncols,
        figsize=(3.2 * ncols, 2.6 * nrows),
        squeeze=False,
        constrained_layout=True
    )

    for row, (label, stop) in enumerate(zip(labels, observation_points)):
        ptmp = master_particles.copy()
        if stop is not None:

            ################################## Enable RADIATION #########################################
            line.configure_radiation(model = "mean")

            ########################### Taper to local momentum deviation ########################################
            tw0 = line.twiss(
                start   = xt.START,
                end     = "injp",
                betx    = betx_init,
                bety    = bety_init,
                alfx    = alfx_init,
                alfy    = alfy_init)

            
            if Radiation == True:
                line.configure_radiation(model = "mean")
            else:
                pass
            tw1 = line.twiss(
                start   = xt.START,
                end     = "injp",
                betx    = betx_init,
                bety    = bety_init,
                alfx    = alfx_init,
                alfy    = alfy_init)
                    
            if Tapering == True:
                taper_magnets_to_local_momentum(line, ele_stop = ELE_STOP)
            else:
                pass

            tw2 = line.twiss(
                start   = xt.START,
                end     = "injp",
                betx    = betx_init,
                bety    = bety_init,
                alfx    = alfx_init,
                alfy    = alfy_init)
            
            if Radiation == True:
                line.configure_radiation(model = "quantum")
            else:
                pass
            line.track(ptmp, ele_stop=stop)
            
        else:
            ################################## Enable RADIATION #########################################
            line.configure_radiation(model = "mean")    
            ########################### Taper to local momentum deviation ########################################
            tw0 = line.twiss(
                start   = xt.START,
                end     = "injp",
                betx    = betx_init,
                bety    = bety_init,
                alfx    = alfx_init,
                alfy    = alfy_init)

            line.configure_radiation(model = "mean")

            tw1 = line.twiss(
                start   = xt.START,
                end     = "injp",
                betx    = betx_init,
                bety    = bety_init,
                alfx    = alfx_init,
                alfy    = alfy_init)

            if Tapering == True:
                taper_magnets_to_local_momentum(line, ele_stop = ELE_STOP)
            else:
                pass
        
            tw2 = line.twiss(
                start   = xt.START,
                end     = "injp",
                betx    = betx_init,
                bety    = bety_init,
                alfx    = alfx_init,
                alfy    = alfy_init)

            if Radiation == True:
                line.configure_radiation(model = "quantum")
            else:
                pass
            line.track(ptmp)

        for col, (var1, var2, title) in enumerate(combinations):
            ax = axs1[row, col]
            x = getattr(ptmp, var1)
            y = getattr(ptmp, var2)
            x=x[ptmp.state==1]
            y=y[ptmp.state==1]

            ax.scatter(x, y, s=3, alpha=0.2, color="blue", edgecolor="none")

            if row == 0:
                ax.set_title(title, fontsize=11, fontweight="bold", pad=4)
            if row == nrows - 1:
                ax.set_xlabel(var1, fontsize=10)
            ax.set_ylabel(f"{var2}\n({label})", fontsize=10)
            ax.grid(True, linestyle="--", alpha=0.3)
            ax.set_facecolor("#f7f7f7")
            ax.tick_params(labelsize=8)

    #fig1.suptitle("Beam Phase-Space Scatter Plots at ACW Observation Points",
    #             fontsize=16, fontweight="bold", y=1.02)
    plt.show()
    
if Plot_type == "hist":
    # ===========================================================
    #  FIGURE 2: HISTOGRAM (2D Log Density)
    # ===========================================================
    fig2, axs2 = plt.subplots(
        nrows, ncols,
        figsize=(3.2 * ncols, 2.6 * nrows),
        squeeze=False,
        constrained_layout=True
    )

    h_last = None
    for row, (label, stop) in enumerate(zip(labels, observation_points)):
        ptmp = master_particles.copy()
        if stop is not None:

            ################################## Enable RADIATION #########################################
            line.configure_radiation(model = "mean")

            ########################### Taper to local momentum deviation ########################################
            tw0 = line.twiss(
                start   = xt.START,
                end     = "injp",
                betx    = betx_init,
                bety    = bety_init,
                alfx    = alfx_init,
                alfy    = alfy_init)

            
            if Radiation == True:
                line.configure_radiation(model = "mean")
            else:
                pass
            tw1 = line.twiss(
                start   = xt.START,
                end     = "injp",
                betx    = betx_init,
                bety    = bety_init,
                alfx    = alfx_init,
                alfy    = alfy_init)
                    
            if Tapering == True:
                taper_magnets_to_local_momentum(line, ele_stop = ELE_STOP)
            else:
                pass

            tw2 = line.twiss(
                start   = xt.START,
                end     = "injp",
                betx    = betx_init,
                bety    = bety_init,
                alfx    = alfx_init,
                alfy    = alfy_init)
            
            if Radiation == True:
                line.configure_radiation(model = "quantum")
            else:
                pass
            line.track(ptmp, ele_stop=stop)
            
        else:
            ################################## Enable RADIATION #########################################
            line.configure_radiation(model = "mean")    
            ########################### Taper to local momentum deviation ########################################
            tw0 = line.twiss(
                start   = xt.START,
                end     = "injp",
                betx    = betx_init,
                bety    = bety_init,
                alfx    = alfx_init,
                alfy    = alfy_init)

            line.configure_radiation(model = "mean")

            tw1 = line.twiss(
                start   = xt.START,
                end     = "injp",
                betx    = betx_init,
                bety    = bety_init,
                alfx    = alfx_init,
                alfy    = alfy_init)

            if Tapering == True:
                taper_magnets_to_local_momentum(line, ele_stop = ELE_STOP)
            else:
                pass
        
            tw2 = line.twiss(
                start   = xt.START,
                end     = "injp",
                betx    = betx_init,
                bety    = bety_init,
                alfx    = alfx_init,
                alfy    = alfy_init)

            if Radiation == True:
                line.configure_radiation(model = "quantum")
            else:
                pass
            line.track(ptmp)

        for col, (var1, var2, title) in enumerate(combinations):
            ax = axs2[row, col]
            x = getattr(ptmp, var1)
            y = getattr(ptmp, var2)
            x=x[ptmp.state==1]
            y=y[ptmp.state==1]
            h = ax.hist2d(x, y, bins=51, cmap="viridis", norm=LogNorm(vmin=1, vmax=None))
            h_last = h

            if row == 0:
                ax.set_title(title, fontsize=11, fontweight="bold", pad=4)
            if row == nrows - 1:
                ax.set_xlabel(var1, fontsize=10)
            ax.set_ylabel(f"{var2}\n({label})\n Alive = {np.sum(ptmp.state == 1) / ptmp.state.size}", fontsize=10)
            ax.grid(True, linestyle="--", alpha=0.3)
            ax.set_facecolor("#f7f7f7")
            ax.tick_params(labelsize=8)

    # Shared colorbar
    cbar = fig2.colorbar(h_last[3], ax=axs2, location="right", shrink=0.85, pad=0.02)
    cbar.set_label("Particle Density (log scale)", fontsize=10)
    cbar.ax.tick_params(labelsize=8)

    #fig2.suptitle("Beam Phase-Space Density Maps (Histograms) at ACW Observation Points",
                #fontsize=16, fontweight="bold", y=1.02)
    plt.show()

# ===========================================================
#  FIGURE 3: 1D HISTOGRAM OF DELTA EVOLUTION
# ===========================================================

fig3, ax3 = plt.subplots(figsize=(8, 5))

colors = plt.cm.plasma(np.linspace(0, 1, len(observation_points)))

for i, (label, stop, color) in enumerate(zip(labels, observation_points, colors)):
    ptmp = master_particles.copy()
    
    ################################## Enable RADIATION #########################################
    line.configure_radiation(model = "mean")

    ########################### Taper to local momentum deviation ########################################
    tw0 = line.twiss(
        start   = xt.START,
        end     = "injp",
        betx    = betx_init,
        bety    = bety_init,
        alfx    = alfx_init,
        alfy    = alfy_init)

    
    if Radiation == True:
        line.configure_radiation(model = "mean")
    else:
        pass
    tw1 = line.twiss(
        start   = xt.START,
        end     = "injp",
        betx    = betx_init,
        bety    = bety_init,
        alfx    = alfx_init,
        alfy    = alfy_init)
            
    if Tapering == True:
        taper_magnets_to_local_momentum(line, ele_stop = ELE_STOP)
    else:
        pass

    tw2 = line.twiss(
        start   = xt.START,
        end     = "injp",
        betx    = betx_init,
        bety    = bety_init,
        alfx    = alfx_init,
        alfy    = alfy_init)
    
    if Radiation == True:
        line.configure_radiation(model = "quantum")
    else:
        pass
    line.track(ptmp, ele_stop=stop)
    

    # Only keep surviving particles
    delta_alive = ptmp.delta[ptmp.state == 1]

    ax3.hist(
        delta_alive * 1e3, bins=101, density=True, histtype='step',
        linewidth=2.0, color=color, label=f"{label} ({np.sum(ptmp.state == 1)/ptmp.state.size:.1%} alive)"
    )

ax3.set_xlabel("Delta [×10⁻³]", fontsize=12)
ax3.set_ylabel("Normalised density", fontsize=12)
ax3.set_title("Evolution of delta along ECS ", fontsize=14, fontweight='bold')
ax3.grid(alpha=0.4, linestyle="--")
ax3.legend(fontsize=9, frameon=True, facecolor='white', edgecolor='gray', loc="upper right")

plt.tight_layout()
plt.show()

######################### Looking into R56 and chirp at INJP ######################################

######################### Track and compute R56 and chirp at INJP ###############################

def R56_and_chirp_at_injp(voltage):
    """
    Compute R56 and chirp h at 'injp' for a given ACW voltage using regression.
    RF effect is included automatically.
    """
    # Set ACW voltage
    line["volt_acw"] = voltage

    # Create fresh particle copy
    ptmp = master_particles.copy()
    delta0 = ptmp.delta.copy()

    ################################## Enable RADIATION #########################################
    line.configure_radiation(model = "mean")

    ########################### Taper to local momentum deviation ########################################
    tw0 = line.twiss(
        start   = xt.START,
        end     = "injp",
        betx    = betx_init,
        bety    = bety_init,
        alfx    = alfx_init,
        alfy    = alfy_init)

    
    if Radiation == True:
        line.configure_radiation(model = "mean")
    else:
        pass
    tw1 = line.twiss(
        start   = xt.START,
        end     = "injp",
        betx    = betx_init,
        bety    = bety_init,
        alfx    = alfx_init,
        alfy    = alfy_init)
            
    if Tapering == True:
        taper_magnets_to_local_momentum(line, ele_stop = ELE_STOP)
    else:
        pass

    tw2 = line.twiss(
        start   = xt.START,
        end     = "injp",
        betx    = betx_init,
        bety    = bety_init,
        alfx    = alfx_init,
        alfy    = alfy_init)
    
    if Radiation == True:
        line.configure_radiation(model = "quantum")
    else:
        pass
    # Track to injection point
    line.track(ptmp, ele_stop="injp")

    # Surviving particles
    mask = ptmp.state == 1
    x = delta0[mask]          # initial delta
    z = ptmp.zeta[mask]       # final zeta including RF effect
    delta = ptmp.delta[mask]  # final delta for chirp

    # R56 via regression: z = R56 * delta0
    R56 = np.sum((x - x.mean())*(z - z.mean())) / np.sum((x - x.mean())**2)

    # Chirp via regression: delta = h * z
    h = np.sum((z - z.mean())*(delta - delta.mean())) / np.sum((z - z.mean())**2)

    return R56, h


# Voltage scan [MV]
voltages = np.linspace(0, 100e6, 13)  # 0, 10, 20, ..., 100 MV

R56_values = []
chirp_values = []

for v in voltages:
    R56_v, h_v = R56_and_chirp_at_injp(v)
    print(f"Voltage: {v*1e-6:.1f} MV -> R56: {R56_v*1e3:.3f} mm, Chirp h: {h_v:.3e} 1/m")
    R56_values.append(R56_v)
    chirp_values.append(h_v)

R56_values = np.array(R56_values)
chirp_values = np.array(chirp_values)


# --- Plot ---
import matplotlib.pyplot as plt

fig, ax1 = plt.subplots(figsize=(8,5))

ax1.plot(voltages*1e-6, R56_values*1e3, 'o-', color='C0', label='R56 [mm]')
ax1.set_xlabel("ACW Voltage [MV]", fontsize=12)
ax1.set_ylabel("R56 [mm]", color='C0', fontsize=12)
ax1.tick_params(axis='y', labelcolor='C0')

ax2 = ax1.twinx()
ax2.plot(voltages*1e-6, chirp_values, 's--', color='C1', label='Chirp h [1/m]')
ax2.set_ylabel("Chirp h [1/m]", color='C1', fontsize=12)
ax2.tick_params(axis='y', labelcolor='C1')

plt.title("R56 and Chirp at INJP vs ACW Voltage", fontsize=14)
fig.tight_layout()
plt.grid(True, linestyle='--', alpha=0.4)
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')
plt.show()

########################### TRACK AND COMPUTE R56 HIGHGER ORDER TERMS AT INJP ######################################


# def R56_T566_U5666_and_chirp(voltage, order=3):
#     """
#     Compute R56, T566, U5666 (up to cubic) and chirp h at 'injp' for surviving particles.

#     Parameters
#     ----------
#     voltage : float
#         ACW voltage [V]
#     order : int
#         Maximum order of polynomial fit for z(delta0)
#         1: R56 only
#         2: R56 + T566
#         3: R56 + T566 + U5666

#     Returns
#     -------
#     dict
#         'R56', 'T566', 'U5666' (up to requested order), 'h' (linear chirp)
#     """
#     line["volt_acw"] = voltage

#     # Copy particles
#     ptmp = master_particles.copy()
#     delta0 = ptmp.delta.copy()

#     # Track to injection point
#     line.track(ptmp, ele_stop="injp")

#     # Only surviving particles
#     mask = ptmp.state == 1
#     x = delta0[mask]
#     z = ptmp.zeta[mask]
#     delta = ptmp.delta[mask]

#     result = {}

#     # --- Polynomial regression of z vs delta0 ---
#     if order >= 1:
#         coeffs = np.polyfit(x, z, deg=order)
#         # np.polyfit returns [highest order, ..., constant]
#         # R56 = linear coefficient
#         result['R56'] = coeffs[-2] if order >= 1 else 0.0
#     if order >= 2:
#         result['T566'] = coeffs[-3] if order >= 2 else 0.0
#     if order >= 3:
#         result['U5666'] = coeffs[-4] if order >= 3 else 0.0

#     # --- Chirp using linear regression of delta vs z ---
#     dz_var = np.sum((z - z.mean())**2)
#     if dz_var > 1e-20:
#         h = np.sum((z - z.mean()) * (delta - delta.mean())) / dz_var
#     else:
#         h = np.nan
#     result['h'] = h

#     return result

# voltages = [0, 50e6, 100e6]
# results = {}

# for v in voltages:
#     res = R56_T566_U5666_and_chirp(v, order=3)
#     print(f"Voltage = {v*1e-6:.0f} MV -> "
#           f"R56 = {res['R56']*1e3:.3f} mm, "
#           f"T566 = {res['T566']*1e3:.3f} mm, "
#           f"U5666 = {res['U5666']*1e3:.3f} mm, "
#           f"h = {res['h']:.3e} 1/m")
    
#     # Assign result to dictionary by voltage
#     results[v] = res


# # Extract values
# R56_vals = [results[v]['R56']*1e3 for v in voltages]     # mm
# T566_vals = [results[v]['T566']*1e3 for v in voltages]   # mm
# U5666_vals = [results[v]['U5666']*1e3 for v in voltages] # mm
# h_vals = [results[v]['h'] for v in voltages]

# # Plot
# fig, ax1 = plt.subplots(figsize=(8,5))

# ax1.plot([v*1e-6 for v in voltages], R56_vals, 'o-', label='R56 [mm]')
# ax1.plot([v*1e-6 for v in voltages], T566_vals, 's--', label='T566 [mm]')
# ax1.plot([v*1e-6 for v in voltages], U5666_vals, 'd-.', label='U5666 [mm]')
# ax1.set_xlabel("ACW Voltage [MV]")
# ax1.set_ylabel("Dispersion terms [mm]")
# ax1.grid(True, linestyle='--', alpha=0.3)
# ax1.legend(loc='upper left')

# ax2 = ax1.twinx()
# ax2.plot([v*1e-6 for v in voltages], h_vals, 'k^-', label='Chirp h [1/m]')
# ax2.set_ylabel("Chirp h [1/m]")
# ax2.tick_params(axis='y', labelcolor='k')
# ax2.legend(loc='upper right')

# plt.title("R56, T566, U5666 and Chirp vs ACW Voltage", fontsize=14)
# plt.tight_layout()
# plt.show()


# ############################# TRACK AND COMPUTE R56 ALONG ECS ######################################

# observation_points = [ "injp"]
# #observation_points = [ "acw.0", "acw.1", "acw.2", "acw.3", "lw4d3", "injp"]

# #labels = ["Before ACW.0", "Before ACW.1", "Before ACW.2", "Before ACW.3", "After ACW.3", "INJP"]
# labels = ["INJP"]

# def get_chirp_R56_at_stop(line, particles, stop):
#     """
#     Compute R56 and linear chirp h = d(delta)/dz at a given stop.
#     Includes RF effects automatically.
#     """
#     # Set ACW voltage
#     line["volt_acw"] = 100E6
    
#     ptmp = particles.copy()
#     delta0 = ptmp.delta.copy()
    
#     # Track particles
#     line.track(ptmp, ele_stop=stop)
    
#     # Only surviving particles
#     mask = ptmp.state == 1
#     z = ptmp.zeta[mask]
#     delta = ptmp.delta[mask]
#     x = delta0[mask]
    
#     # R56: slope of z vs initial delta
#     R56 = np.sum((x - x.mean())*(z - z.mean())) / np.sum((x - x.mean())**2)
    
#     # Chirp: slope of delta vs z
#     h = np.sum((z - z.mean())*(delta - delta.mean())) / np.sum((z - z.mean())**2)
    
#     return R56, h

# R56_values = []
# chirp_values = []

# for stop in observation_points:
#     R56, h = get_chirp_R56_at_stop(line, master_particles, stop)
#     R56_values.append(R56)
#     chirp_values.append(h)
#     print(f"{stop}: R56 = {R56*1e3:.3f} mm, chirp h = {h:.3e} 1/m")

# # Pick a single particle for s positions
# single_particle = xt.Particles(
#     p0c   = master_particles.p0c,
#     q0    = master_particles.q0,
#     mass0 = master_particles.mass0,
#     x     = master_particles.x[:1],
#     px    = master_particles.px[:1],
#     y     = master_particles.y[:1],
#     py    = master_particles.py[:1],
#     zeta  = master_particles.zeta[:1],
#     delta = master_particles.delta[:1],
#     state = master_particles.state[:1],
# )
# s_map = []
# for stop in observation_points:
#     line.track(single_particle, ele_stop=stop)
#     s_map.append(single_particle.s[0])


# fig, ax1 = plt.subplots(figsize=(10,5))

# ax1.plot(s_map, np.array(R56_values)*1e3, 'o-', color='C0', label='R56 [mm]')
# ax1.set_xlabel("s [m]", fontsize=12)
# ax1.set_ylabel("R56 [mm]", color='C0', fontsize=12)
# ax1.tick_params(axis='y', labelcolor='C0')

# ax2 = ax1.twinx()
# ax2.plot(s_map, np.array(chirp_values), 's--', color='C1', label='Chirp h [1/m]')
# ax2.set_ylabel("Chirp h [1/m]", color='C1', fontsize=12)
# ax2.tick_params(axis='y', labelcolor='C1')

# plt.title("R56 and Chirp along ECS (RF included)", fontsize=14)
# ax1.legend(loc='upper left')
# ax2.legend(loc='upper right')
# plt.grid(True, linestyle='--', alpha=0.4)
# plt.tight_layout()
# plt.show()



# ######################### TRACK AND COMPUTE R56(s) ALONG ECS FOR VOLTAGE SCAN ######################################

# # Observation points
# observation_points = ["acw.0", "acw.1", "acw.2", "acw.3", "lw4d3", "injp"]
# labels = ["Before ACW.0", "Before ACW.1", "Before ACW.2", "Before ACW.3", "After ACW.3", "INJP"]

# # Single particle for s positions
# single_particle = xt.Particles(
#     p0c   = master_particles.p0c,
#     q0    = master_particles.q0,
#     mass0 = master_particles.mass0,
#     x     = master_particles.x[:1],
#     px    = master_particles.px[:1],
#     y     = master_particles.y[:1],
#     py    = master_particles.py[:1],
#     zeta  = master_particles.zeta[:1],
#     delta = master_particles.delta[:1],
#     state = master_particles.state[:1],
# )

# s_map = []
# for stop in observation_points:
#     line.track(single_particle, ele_stop=stop)
#     s_map.append(single_particle.s[0])

# # Voltage scan [MV]
# voltages = np.linspace(0, 100e6, 11)  # 0, 10, ..., 100 MV

# # Store results
# R56_vs_voltage = []
# chirp_vs_voltage = []

# for v in voltages:
#     R56_values = []
#     chirp_values = []

#     for stop in observation_points:
#         # Set voltage
#         line["volt_acw"] = v

#         # Track particles
#         ptmp = master_particles.copy()
#         delta0 = ptmp.delta.copy()
#         line.track(ptmp, ele_stop=stop)

#         mask = ptmp.state == 1
#         z = ptmp.zeta[mask]
#         delta = ptmp.delta[mask]
#         x = delta0[mask]

#         # R56: z vs initial delta
#         R56 = np.sum((x - x.mean())*(z - z.mean())) / np.sum((x - x.mean())**2)
#         # Chirp: delta vs z
#         h = np.sum((z - z.mean())*(delta - delta.mean())) / np.sum((z - z.mean())**2)

#         R56_values.append(R56)
#         chirp_values.append(h)

#     R56_vs_voltage.append(R56_values)
#     chirp_vs_voltage.append(chirp_values)

# R56_vs_voltage = np.array(R56_vs_voltage)  # shape: (n_voltage, n_points)
# chirp_vs_voltage = np.array(chirp_vs_voltage)

# # --- Plot R56 and chirp vs voltage ---
# fig, axes = plt.subplots(2, 1, figsize=(10,8), sharex=True)

# for i, label in enumerate(labels):
#     axes[0].plot(voltages*1e-6, R56_vs_voltage[:, i]*1e3, 'o-', label=label)
#     axes[1].plot(voltages*1e-6, chirp_vs_voltage[:, i], 's--', label=label)

# axes[0].set_ylabel("R56 [mm]", fontsize=12)
# axes[0].set_title("R56 vs ACW Voltage at ECS Stops", fontsize=14)
# axes[0].grid(True, linestyle='--', alpha=0.4)
# axes[0].legend(fontsize=9)

# axes[1].set_xlabel("ACW Voltage [MV]", fontsize=12)
# axes[1].set_ylabel("Chirp h [1/m]", fontsize=12)
# axes[1].set_title("Chirp vs ACW Voltage at ECS Stops", fontsize=14)
# axes[1].grid(True, linestyle='--', alpha=0.4)
# axes[1].legend(fontsize=9)

# plt.tight_layout()
# plt.show()
