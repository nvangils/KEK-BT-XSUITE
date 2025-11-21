"""
Plot phase space at a chosen lag, voltage for the ECS at 
before/after each cavity along the ECS
User choices include radiation ON/OFF Tapering ON/OFF and energy of line
with respect to HER energy
=============================================
Author(s):  John P T Salvesen, Nikita Z van Gils
Date:       19-11-2025

"""

######################################### Required Packages ########################################

import os
import numpy as np
import xtrack as xt
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
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

#Injection type
#IF synchrotron injection 
ELE_START    = "pbtp"
ELE_STOP    = "injp"                 # FINAL ELEMENT TO TRACK TO
Plot_type   = "hist"                # CHOOSE WHICH TYPE OF PLOT TO SHOW PHASE SPACE ["scatter", "hist"]
SYNCHROTRON_OFFSET_BT = 0 

########################################
# Particles Filepath
########################################
particles_filepath      = "../data/positrons_4nC.npy"
#particles_filepath   = "../data/yohsimotosan.dat"

########################################
# Initial Conditions
########################################
betx_init   = 18.38101336810373  
bety_init   = 7.320112051988261 
alfx_init   = 1.0683347351664612
alfy_init   = 0.517728709933477

########################################
# Tracking Parameters
########################################


#################################### Synchrotron Injection Energy Offset #######################

SYNCHROTRON_INJ_E_SCALE = 1.000 + SYNCHROTRON_OFFSET_BT          # Relative energy so that increase at injection is 1.006 
REFERENCE_ENERGY_HER    = 3.98E9                        # [GeV], reference energy of HER

##################################### Final Element ############################################

ELE_STOP                = ELE_STOP

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
# Load lattice
line = xt.Line.from_json("../lattices/btp.json")
env = line.env
# Apply the optics
env.call("../lattices/btp_import_optics.py")


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



# #########  Load initial distribution ############################################
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

# initial = xt.Table({
#     "name":  np.arange(initial_distribution[:, 0].shape[0]),
#     "x":     initial_distribution[:, 0] - np.median(initial_distribution[:, 0]),
#     "px":    initial_distribution[:, 1] - np.median(initial_distribution[:, 1]),
#     "y":     initial_distribution[:, 2] - np.median(initial_distribution[:, 2]),
#     "py":    initial_distribution[:, 3] - np.median(initial_distribution[:, 3]),
#     "zeta":  initial_distribution[:, 4] - np.median(initial_distribution[:, 4]),
#     "delta": initial_distribution[:, 5] - np.median(initial_distribution[:, 5]),
#     "state": np.full(initial_distribution[:, 0].shape, 1)
# })
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

# Make sure all particle arrays are regular numpy arrays
x     = np.array(master_particles.x)
px    = np.array(master_particles.px)
y     = np.array(master_particles.y)
py    = np.array(master_particles.py)
zeta  = np.array(master_particles.zeta)
delta = np.array(master_particles.delta)
state = np.array(master_particles.state)


# Start with all particles alive
master_particles.state[:] = 1

alive_mask = master_particles.state == 1

#Phase space initial distribution
plt.figure(figsize=(7,5))
plt.scatter(master_particles.zeta*1e3, master_particles.delta*1e3, 
            s=5, alpha=0.3, color='C0')
plt.xlabel("zeta [mm]")
plt.ylabel("delta [rel]")
plt.title("Phase space ")
plt.grid(True)


plt.figure(figsize=(7,5))
plt.scatter(master_particles.x*1e3, master_particles.px*1e3, 
            s=5, alpha=0.3, color='C0')
plt.xlabel("x [mm]")
plt.ylabel("px [mrad]")
plt.title("Phase space ")
plt.grid(True)

plt.figure(figsize=(7,5))
plt.scatter(master_particles.y*1e3, master_particles.py*1e3, 
            s=5, alpha=0.3, color='C0')
plt.xlabel("y [mm]")
plt.ylabel("py [rel]")
plt.title("Phase space")
plt.grid(True)


plt.show()


####################### PLOT TRANSVERSE INPUT DISTRIBUTION #####################################

# Get particle coordinates
x = master_particles.x
y = master_particles.y
z = master_particles.zeta

# Compute RMS sizes
rms_x = np.sqrt(np.mean(x**2))
rms_y = np.sqrt(np.mean(y**2))
rms_z = np.sqrt(np.mean(z**2))

print(f"RMS x: {rms_x*1e3:.3f} mm")
print(f"RMS y: {rms_y*1e3:.3f} mm")
print(f"RMS z: {rms_z*1e3:.3f} mm")

# --- Set up figure ---
fig, axs = plt.subplots(1, 3, figsize=(18, 5))

# x-y plane
axs[0].scatter(x*1e3, y*1e3, s=5, alpha=0.7)
axs[0].set_xlabel("x [mm]")
axs[0].set_ylabel("y [mm]")
axs[0].set_title(f"x-y projection\nRMS x = {rms_x*1e3:.3f} mm, RMS y = {rms_y*1e3:.3f} mm")
axs[0].grid(True, linestyle='--', alpha=0.3)

# x-z plane
axs[1].scatter(x*1e3, z*1e3, s=5, alpha=0.7)
axs[1].set_xlabel("x [mm]")
axs[1].set_ylabel("z [mm]")
axs[1].set_title(f"x-z projection\nRMS x = {rms_x*1e3:.3f} mm, RMS z = {rms_z*1e3:.3f} mm")
axs[1].grid(True, linestyle='--', alpha=0.3)

# y-z plane
axs[2].scatter(y*1e3, z*1e3, s=5, alpha=0.7)
axs[2].set_xlabel("y [mm]")
axs[2].set_ylabel("z [mm]")
axs[2].set_title(f"y-z projection\nRMS y = {rms_y*1e3:.3f} mm, RMS z = {rms_z*1e3:.3f} mm")
axs[2].grid(True, linestyle='--', alpha=0.3)

plt.tight_layout()
plt.show()


# --- Observation points and phase-space combinations ---
observation_points = ["injp"]
labels = ["INJP"]
combinations = [
    ("x", "y", "X-Y"),
    ("x", "zeta", "X-Zeta"),
    ("x", "delta", "X-Delta"),
    ("y", "zeta", "Y-Zeta"),
    ("y", "delta", "Y-Delta"),
    ("zeta", "delta", "Zeta-Delta")
]

nrows, ncols = len(observation_points), len(combinations)


# Observation points
observation_points = ["pbtp", "injp"]
labels = [ "pbtp", "INJP"]

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
                taper_magnets_to_local_momentum(line, ele_start= ELE_START, ele_stop = ELE_STOP)
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
            line.track(ptmp,ele_start= ELE_START,  ele_stop=stop)
            
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
            line.track(ptmp, ele_start= ELE_START, ele_stop= ELE_STOP)

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
    line.track(ptmp, ele_start= stop, ele_stop= ELE_STOP)
    

    # Only keep surviving particles
    delta_alive = ptmp.delta[ptmp.state == 1]

    ax3.hist(
        delta_alive * 1e3, bins=101, density=True, histtype='step',
        linewidth=2.0, color=color, label=f"{label} ({np.sum(ptmp.state == 1)/ptmp.state.size:.1%} alive)"
    )

ax3.set_xlabel("Delta [×10⁻³]", fontsize=12)
ax3.set_ylabel("Normalised density", fontsize=12)
ax3.set_title("Evolution of delta along BTp ", fontsize=14, fontweight='bold')
ax3.grid(alpha=0.4, linestyle="--")
ax3.legend(fontsize=9, frameon=True, facecolor='white', edgecolor='gray', loc="upper right")

plt.tight_layout()
plt.show()

