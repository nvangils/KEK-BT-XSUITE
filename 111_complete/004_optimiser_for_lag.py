

"""
Check the tracked beam against the previously computed DA for BTE lattice.
This maps the BTe tracked beam at INJP and compares it to the DA at the septum (everything in DA space)
Reference frame transform moves beam at INJP into reference frame of bump.
User choices: emittance blowup at the end of the line, radiation ON/OFF tapering ON/OFF ECS settings (LAG SCAN) and injection type
betatron/synchrotron or hybrid 
Output:
MAXIMAL survival fraction within DA of HER for fixed ECS voltage scanning the ECS lag
Optimiser steps show dependence of survivability on ECS lag for a fixed ECS voltage 

=============================================
Author(s):  John P T Salvesen, Nikita Z van Gils
Date:       18-11-2025

"""

######################################### Required Packages ########################################

import xtrack as xt
import xobjects as xo
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import sys
from scipy.optimize import minimize_scalar
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
from helpers._calculate_action_amplitude import calculate_action_amplitude
from helpers._fit_ellipse_da import calculate_boundary_points
from helpers._check_da_frontier import compare_to_frontier
from helpers._local_momentum_taper import taper_magnets_to_local_momentum
from helpers._compute_emittance import compute_emittance_and_optics

######################################## User Parameters ##########################################
#Emittance blowup:
EMITTANCE_X_MEASURED = 140e-6/14000	 # MEASURED EMITTANCE AT BTE2 IN NORMALISED
EMITTANCE_Y_MEASURED = 140e-6/14000	 # MEASURED EMITTANCE AT BTE2 IN NORMALISED
#Line configuration:
Radiation = True 	                 # TRUE ENABLES RADIATION 
Tapering = True 	                 # TRUE ENABLES TAPERING 
aperture_check = True                # TRUE ENABLES APERTURE CHECKING
#ECS settings:
ECS_voltage = None	                 # SET A VALUE IN MV TO ENABLE ELSE DEFAULT IS 18MV
set_lag_bounds = (179, 180)          # BOUNDS FOR THE OPTIMISER

#Injection type
#Betatron injection: 
BETATRON_OFFSET= 7.6e-3              # SET TO 0 IF PURE SYNCHROTRON INJECTION
BETATRON_ANGLE_OFFSET= 0.6e-3        # SET TO 0 IF PURE SYNCHROTRON INJECTION
#Synchrotron injection 
SYNCHROTRON_OFFSET_BT= 0.00165       # IN % WITH RESPECT TO  HER ENERGY. SET TO 0 FOR BETATRON INJECTION

######################################### Lattice Files ########################################

xsuite_lattice_path     = "../lattices/bte.py"
xsuite_optics_path      = "../lattices/bte_import_optics.py"

######################################## Particles Filepath ####################################

#particles_filepath      = "../data/input_electrons_wake.npy"
particles_filepath      = "../data/yohsimotosan.dat"
######################################### DA Filepaths ##########################################

##For 200_8 mm lattices uncomment below 
da_xy_filepath         = "../data/da_xy_folded_her_oide_200_8_her.npz"
da_zx_filepath         = "../data/da_zx_folded_her_oide_200_8_her.npz"
da_zy_filepath         = "../data/da_zy_folded_her_oide_200_8_her.npz"

##For 100_3 mm lattices uncomment below 
# da_xy_filepath         = "../data/da_xy_folded_her_oide_100_3_her.npz"
# da_zx_filepath         = "../data/da_zx_folded_her_oide_100_3_her.npz"
# da_zy_filepath         = "../data/da_zy_folded_her_oide_100_3_her.npz"

######################################### Initial Conditions ####################################

betx_init               = 20.540437095810788 
bety_init               = 18.015589974012826 
alfx_init               = 1.0030653048224993 
alfy_init               = .37003884861294073 

#################################### Synchrotron Injection Energy Offset #######################

SYNCHROTRON_INJ_E_SCALE = 1.000 + SYNCHROTRON_OFFSET_BT          # Relative energy so that increase at injection is 1.006 
REFERENCE_ENERGY_HER    = 7.00E9           # [GeV], reference energy of HER


#################################### Lattice Injection Optics Parameters #######################

##For 100_3 mm lattices uncomment below 
# BETX_HER_INJECTIO       = 99.9986 
# BETY_HER_INJECTIO       = 2.55671
# ALFX_HER_INJECTIO       = -8.50447 
# ALFY_HER_INJECTIO       = -1.47663
# DX_HER_INJECTIO         = -1.6
# DY_HER_INJECTIO         = 0.0
# DDX_HER_INJECTIO        = 26.7511
# DDY_HER_INJECTIO        = 0.00048039
# DDPX_HER_INJECTIO       = 5.2278e-05
# DDPY_HER_INJECTIO       = -0.000376
# DPX_HER_INJECTIO        = -0.120117  
# DPY_HER_INJECTIO        = -4.12801e-07
# BETS0_HER               = 8.056774638483907

##For 200_8 mm lattices uncomment below 

BETX_HER_INJECTIO       = 99.9987
BETY_HER_INJECTIO       = 2.55671
ALFX_HER_INJECTIO       = -8.50448
ALFY_HER_INJECTIO       = -1.47663
DX_HER_INJECTIO         = -1.60 
DY_HER_INJECTIO         = 0.0
DDX_HER_INJECTIO        = 23.70  
DDY_HER_INJECTIO        = -0.0001
DDPX_HER_INJECTIO       = 1.96375
DDPY_HER_INJECTIO       = -0.0004
DPX_HER_INJECTIO        = -0.120
DPY_HER_INJECTIO        = 0.0
BETS0_HER               = 8.0346


##################################### Final Element ############################################

ELE_STOP                = "injp"


#################################### Injection Parameters ######################################

X_SEPTUM_CUT            = 35.8E-3 * -1           # [m],      Septum edge position, w.r.t. HER design orbit (INJECTIO)
X_INJECTION_POINT       = 36.9E-3 * -1           # [m],      INJP x offset, w.r.t. HER design orbit (INJECTIO)
PX_INJECTION_POINT      = 2.5E-3 * -1           # [rad],    INJP px offset, w.r.t. HER design orbit (INJECTIO)
X_ORBIT_BUMP            = X_INJECTION_POINT - DX_HER_INJECTIO * (SYNCHROTRON_INJ_E_SCALE-1) -1/2*(DDX_HER_INJECTIO * (SYNCHROTRON_INJ_E_SCALE-1)**2)-BETATRON_OFFSET
PX_ORBIT_BUMP           = PX_INJECTION_POINT - DPX_HER_INJECTIO * (SYNCHROTRON_INJ_E_SCALE-1)-1/2*(DDPX_HER_INJECTIO * (SYNCHROTRON_INJ_E_SCALE-1)**2) -BETATRON_ANGLE_OFFSET

##################################### Aperture Parameters ######################################

# T. Mori:
#   Round cross-section with 47 mm diameter at drifts and Q-magnets
#   Rectangle with about 30 mm height and 50-100 mm width at the dipoles
RADIUS_BP_APERTURE      = 0.0235           # [m],      Beam pipe radius at drifts and Q-magnets
HEIGHT_BEND_APERTURE    = 0.03             # [m],      Bend aperture height
WIDTH_BEND_APERTURE     = 0.05             # [m],      Bend aperture width

######################################### Tracking Context #####################################

N_THREADS             = 10                # Number of CPU threads for tracking


######################################### Load DA Data ##########################################

da_xy   = np.load(file = da_xy_filepath, allow_pickle = True)
da_zx   = np.load(file = da_zx_filepath, allow_pickle = True)
da_zy   = np.load(file = da_zy_filepath, allow_pickle = True)

######################################### Convert Format ########################################

NORM_AX_GRID_XY = da_xy["FOLDED_NORM_A_GRID_HER"]
NORM_AY_GRID_XY = da_xy["FOLDED_NORM_B_GRID_HER"]
DA_XY           = da_xy["FOLDED_SURVIVING_HER"]
GEMITT_X_XY     = da_xy["GEMITT_X"]
GEMITT_Y_XY     = da_xy["GEMITT_Y"]
GEMITT_Z_XY     = da_xy["GEMITT_Z"]

NORM_AZ_GRID_ZX = da_zx["FOLDED_NORM_A_GRID_HER"]
NORM_AX_GRID_ZX = da_zx["FOLDED_NORM_B_GRID_HER"]
DA_ZX           = da_zx["FOLDED_SURVIVING_HER"]
GEMITT_X_ZX     = da_zx["GEMITT_X"]
GEMITT_Y_ZX     = da_zx["GEMITT_Y"]
GEMITT_Z_ZX     = da_zx["GEMITT_Z"]

NORM_AZ_GRID_ZY = da_zy["FOLDED_NORM_A_GRID_HER"]
NORM_AY_GRID_ZY = da_zy["FOLDED_NORM_B_GRID_HER"]
DA_ZY           = da_zy["FOLDED_SURVIVING_HER"]
GEMITT_X_ZY     = da_zy["GEMITT_X"]
GEMITT_Y_ZY     = da_zy["GEMITT_Y"]
GEMITT_Z_ZY     = da_zy["GEMITT_Z"]

######################################## Ensure Consistency of Emittances ######################

assert np.isclose(GEMITT_X_XY, GEMITT_X_ZX) and np.isclose(GEMITT_X_XY, GEMITT_X_ZY)
assert np.isclose(GEMITT_Y_XY, GEMITT_Y_ZX) and np.isclose(GEMITT_Y_XY, GEMITT_Y_ZY)
assert np.isclose(GEMITT_Z_XY, GEMITT_Z_ZX) and np.isclose(GEMITT_Z_XY, GEMITT_Z_ZY)
GEMITT_X_HER   = GEMITT_X_XY
GEMITT_Y_HER   = GEMITT_Y_XY
GEMITT_Z_HER   = GEMITT_Z_XY

######################################### Calculate DA Boundary Points #################################

boundary_a_xy, boundary_b_xy  = calculate_boundary_points(
    NORM_AX_GRID_XY, NORM_AY_GRID_XY, DA_XY)
boundary_a_zx, boundary_b_zx  = calculate_boundary_points(
    NORM_AZ_GRID_ZX, NORM_AX_GRID_ZX, DA_ZX)
boundary_a_zy, boundary_b_zy  = calculate_boundary_points(
    NORM_AZ_GRID_ZY, NORM_AY_GRID_ZY, DA_ZY)

######################################### Load Lattice ###########################################

env     = xt.Environment()
env.call(xsuite_lattice_path)
env.call(xsuite_optics_path)
line    = env.lines["line"]
if ECS_voltage is not None:
    line["volt_acw"] = ECS_voltage * 1E6
else :
    line["volt_acw"] = 18E6


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

######################################### Get bends and other active elements #####################

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


###################### Get a list of which types of elements need aperture ########################

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


####################### Setup #######################################################################

####################### Get clean table ######################

tt  = line.get_table(attr = True)

######################################### Load Particles Distribution ############################

#initial_distribution    = np.load(particles_filepath, allow_pickle = True)
initial_distribution    = np.loadtxt(particles_filepath)

####### Set up bt to co bump at injection point ###########
bt_to_co_bump_rotation      = xt.YRotation(
    angle   = -1 * np.rad2deg(PX_INJECTION_POINT - PX_ORBIT_BUMP)) #angular rotation
bt_to_co_bump_translation   = xt.XYShift(
    dx      = -1 * (X_INJECTION_POINT - X_ORBIT_BUMP)) #N.B. sign convention should be negative

######################################### Prepare Lag Scan ###########################################
survival_results = []
######################################### Lag Optimization ###########################################


# To record optimizer path
eval_history = {"lag": [], "surv": []}

def survival_fraction_for_lag(lag_val):
    """Compute the surviving fraction for a given lag value."""
    line["lag_acw"] = lag_val
    # Reset tracked beam
    initial_distribution = np.loadtxt(particles_filepath)
    initial_distribution = xt.Table({
        "name":  np.arange(initial_distribution[:, 0].shape[0]),
        "x":     initial_distribution[:, 0] - np.median(initial_distribution[:, 0]),
        "px":    initial_distribution[:, 1] - np.median(initial_distribution[:, 1]),
        "y":     initial_distribution[:, 2] - np.median(initial_distribution[:, 2]),
        "py":    initial_distribution[:, 3] - np.median(initial_distribution[:, 3]),
        "zeta":  initial_distribution[:, 4] - np.median(initial_distribution[:, 4]),
        "delta": initial_distribution[:, 5] - np.median(initial_distribution[:, 5]),
        "state": np.full(initial_distribution[:, 0].shape, 1)
    })

    initial_beam = xt.Particles(
        p0c=line.particle_ref.p0c,
        q0=line.particle_ref.q0,
        mass0=line.particle_ref.mass0,
        x=initial_distribution.x, px=initial_distribution.px,
        y=initial_distribution.y, py=initial_distribution.py,
        zeta=initial_distribution.zeta, delta=initial_distribution.delta
    )
    tracked_beam = initial_beam.copy()
    ################################## Setup for tracking #######################################
    line.configure_radiation(model = "mean") #have to reset in order for loop to work
    ################################## Enable RADIATION #########################################

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

    
    ######################################## Track ##############################################
    if Radiation == True:
        line.configure_radiation(model = "quantum")
    else:
        pass

    line.track(particles = tracked_beam, ele_stop  = ELE_STOP)


    ######################################## Perform Transforms to the HER ######################

    # Want these parameters in the frame of the HER stored beam CO bump at INJECTIO / INJP

    ########################## Optimise Orbit Bump #################
    processed_delta = (line.particle_ref.p0c[0] / REFERENCE_ENERGY_HER) * \
        (1 + tracked_beam.delta[tracked_beam.state == 1]) - 1

    mean_delta_at_inj       = np.mean(processed_delta)
    print(f"Mean delta at injection:    {mean_delta_at_inj:.6e}")

    X_ORBIT_BUMP            = X_INJECTION_POINT - DX_HER_INJECTIO * (mean_delta_at_inj) -1/2*(DDX_HER_INJECTIO * (mean_delta_at_inj)**2) - BETATRON_OFFSET
    PX_ORBIT_BUMP           = PX_INJECTION_POINT - DPX_HER_INJECTIO * (mean_delta_at_inj)-1/2*(DDPX_HER_INJECTIO * (mean_delta_at_inj)**2) -BETATRON_ANGLE_OFFSET
    print(f"Optimised X_ORBIT_BUMP:     {X_ORBIT_BUMP*1E3:.3f} mm")
    print(f"Optimised PX_ORBIT_BUMP:    {PX_ORBIT_BUMP*1E3:.3f} mrad")

    ########################## Create Transforms from INJP to CO Bump at INJECTIO #################
    print(f"Mean x pre ref transform:   ", np.mean(tracked_beam.x[tracked_beam.state == 1]))
    print(f"Mean px pre ref transform:  ", np.mean(tracked_beam.px[tracked_beam.state == 1]))

    bt_to_co_bump_rotation      = xt.YRotation(
        angle   = -1 * np.rad2deg(PX_INJECTION_POINT - PX_ORBIT_BUMP)) #angular rotation
    bt_to_co_bump_translation   = xt.XYShift(
        dx      = -1 * (X_INJECTION_POINT - X_ORBIT_BUMP)) #N.B. sign convention should be negative

    bt_to_co_bump_rotation.track(tracked_beam)
    bt_to_co_bump_translation.track(tracked_beam)
    print(np.mean(tracked_beam.px[tracked_beam.state == 1]))



    ######################################### Process Coordinates #################################

    print(f"Proportion that survived up to injection:        {np.sum(tracked_beam.state == 1) / tracked_beam.state.size:.4f}")

    processed_x     = tracked_beam.x[tracked_beam.state == 1]
    processed_px    = tracked_beam.px[tracked_beam.state == 1]
    processed_y     = tracked_beam.y[tracked_beam.state == 1]
    processed_py    = tracked_beam.py[tracked_beam.state == 1]
    processed_zeta  = tracked_beam.zeta[tracked_beam.state == 1]# type: ignore
    processed_delta = (line.particle_ref.p0c[0] / REFERENCE_ENERGY_HER) * \
        (1 + tracked_beam.delta[tracked_beam.state == 1]) - 1
    
    ######################################### Calculate Amplitudes ###############################

    tracked_amplitudes  = xt.Table({
        "name":     np.arange(processed_x.shape[0]),
        "Ax":       calculate_action_amplitude(
            a           = processed_x,
            pa          = processed_px,
            delta       = processed_delta,
            beta        = BETX_HER_INJECTIO,
            alfa        = ALFX_HER_INJECTIO, #ASSUMED OF RING
            da          = DX_HER_INJECTIO,
            dpa         = DPX_HER_INJECTIO,
            gemitt_a    = GEMITT_X_HER,
            dda= DDX_HER_INJECTIO,
            ddpa= DDPX_HER_INJECTIO),
        "Ay":       calculate_action_amplitude(
            a           = processed_y,
            pa          = processed_py,
            delta       = processed_delta,
            beta        = BETY_HER_INJECTIO,
            alfa        = ALFY_HER_INJECTIO,  #ASSUMED OF RING
            da          = DY_HER_INJECTIO,
            dpa         = DPY_HER_INJECTIO,
            gemitt_a    = GEMITT_Y_HER,
            dda= DDY_HER_INJECTIO,
            ddpa= DDPY_HER_INJECTIO),
        "Az":       calculate_action_amplitude(
            a           = processed_zeta,
            pa          = processed_delta,
            delta       = processed_delta,
            beta        = BETS0_HER,
            alfa        = 0.0,
            da          = 0.0,
            dpa         = 0.0,
            gemitt_a    = GEMITT_Z_HER,
            dda= 0.0,
            ddpa= 0.0)})


        ######### RESCALE THE EMITTANCE ACCORDING TO MEASURED BLOWUP ALONG THE LINE #################

    if Radiation == True:
        results = compute_emittance_and_optics(
            x               = tracked_beam.x[tracked_beam.state == 1],
            px              = tracked_beam.px[tracked_beam.state == 1],
            y               = tracked_beam.y[tracked_beam.state == 1],
            py              = tracked_beam.py[tracked_beam.state == 1],
            zeta            = tracked_beam.zeta[tracked_beam.state == 1],
            delta           = tracked_beam.delta[tracked_beam.state == 1],
            particle_axis   = 0,)

        gemitt_x    = results["epsx"]
        gemitt_y    = results["epsy"]
        gemitt_z    = results["epsz"]

        ######################################### Compute Emittances ######################################
        EMITTANCE_X_TRACKED = gemitt_x
        EMITTANCE_Y_TRACKED = gemitt_y
        

        print(f"Tracked Emittance X at INJP:    {EMITTANCE_X_TRACKED*1E6:.3f} um")
        print(f"Tracked Emittance Y at INJP:    {EMITTANCE_Y_TRACKED*1E6:.3f} um")
        tracked_amplitudes  = xt.Table({
            "name":     np.arange(processed_x.shape[0]),
            "Ax":       tracked_amplitudes["Ax"] * np.sqrt(EMITTANCE_X_MEASURED / EMITTANCE_X_TRACKED),
            "Ay":       tracked_amplitudes["Ay"] * np.sqrt(EMITTANCE_Y_MEASURED / EMITTANCE_Y_TRACKED),
            "Az":       tracked_amplitudes["Az"]*1})
        
        print("x rescale factor: ", np.sqrt(EMITTANCE_X_MEASURED / EMITTANCE_X_TRACKED))
        print("y rescale factor: ", np.sqrt(EMITTANCE_Y_MEASURED / EMITTANCE_Y_TRACKED))

    else:
        pass
    
    
    ######################################### XY Test ############################################
    within_xy = compare_to_frontier(
        xs          = boundary_a_xy,
        ys          = boundary_b_xy,
        x_test      = tracked_amplitudes.Ax,
        y_test      = tracked_amplitudes.Ay)
    print(f"Proportion that survived to injection within XY DA:        {np.sum(within_xy) / within_xy.size:.4f}")

    ######################################### ZX Test ############################################

    within_zx = compare_to_frontier(
        xs          = boundary_a_zx,
        ys          = boundary_b_zx,
        x_test      = tracked_amplitudes.Az,
        y_test      = tracked_amplitudes.Ax)
    print(f"Proportion that survived to injection within ZX DA:        {np.sum(within_zx) / within_zx.size:.4f}")

    ######################################### ZY Test ############################################

    within_zy = compare_to_frontier(
        xs          = boundary_a_zy,
        ys          = boundary_b_zy,
        x_test      = tracked_amplitudes.Az,
        y_test      = tracked_amplitudes.Ay)
    print(f"Proportion that survived to injection ZY DA:        {np.sum(within_zy) / within_zy.size:.4f}")

    ##################################### Septum Test ############################################

    assert np.abs(PX_ORBIT_BUMP) < 1E-1, \
        "We assume the CO bump is parallel to the stored beam orbit and therefore the septum"

    within_septum   = processed_x < (X_SEPTUM_CUT - X_ORBIT_BUMP)
    print(f"Proportion within Septum:       {np.sum(within_septum) / within_septum.size:.4f}")

    ############################### Within all test ##############################################

    within  = within_xy * within_zx * within_zy * within_septum
    print(f"Proportion within All:          {np.sum(within) / within.size:.4f}")
    
    survival_fraction = np.sum(within) / within.size
      # Save results for plotting later
    survival_results.append((lag_val, survival_fraction))
    eval_history["lag"].append(lag_val)
    eval_history["surv"].append(survival_fraction)

    print(f"Lag {lag_val:.2f} → Survival: {survival_fraction*100:.3f}%")
    #survival_results.append((lag_val, survival_fraction))

    return float(-survival_fraction)

# Run the optimizer
res = minimize_scalar(survival_fraction_for_lag, bounds=set_lag_bounds, method='bounded', options={'xatol': 0.01})
opt_lag = res.x
opt_fraction = -res.fun
print("\n==================== OPTIMUM FOUND ====================")
print(f"Best lag = {opt_lag:.4f}, Survival = {opt_fraction*100:.4f}, Voltage: {line["volt_acw"] / 1E6} MV")
print("=======================================================\n")


######################################### Find optimum ###########################################

# Extract lag values and survival fractions
survival_results = np.array(survival_results)
lags = survival_results[:, 0]
fractions = survival_results[:, 1]

# Find optimum
opt_idx = np.argmax(fractions)
opt_lag = lags[opt_idx]
opt_fraction = fractions[opt_idx]

# --- Plot ---
plt.figure(figsize=(10, 6))

# Plot survival fraction vs lag
plt.plot(lags, fractions*100, "o-", lw=2, color="blue", label="Surviving fraction")

# Highlight the optimum
plt.scatter(opt_lag, opt_fraction*100, color="red", s=100, zorder=10, label=f"Maximum: {opt_fraction*100:.2f}% at {opt_lag:.2f}°")
plt.axvline(opt_lag, color="red", linestyle="--", lw=2)

# Labels and grid
plt.xlabel("Lag [deg]", fontsize=12)
plt.ylabel("Surviving fraction [%]", fontsize=12)
plt.title("Injection Survivability vs Lag", fontsize=14)
plt.legend(fontsize=11)
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

plt.show()

################## plot optimiser path ###########################################
# === Extract results ===
lags = np.array(eval_history["lag"])
surv = np.array(eval_history["surv"]) * 100

# Sort for smooth line
idx = np.argsort(lags)
lags_sorted = lags[idx]
surv_sorted = surv[idx]

opt_lag = lags[np.argmax(surv)]
opt_surv = np.max(surv)

# === Plot survivability curve & optimizer path ===
fig, ax = plt.subplots(figsize=(10, 6))

# Smooth survivability curve
ax.plot(lags_sorted, surv_sorted, "-", color="blue", lw=2, alpha=0.7, label="Survival vs Lag")

# Points visited by optimizer (in order)
ax.scatter(lags, surv, c=np.arange(len(lags)), cmap="plasma", s=80, edgecolor="k", lw=0.5, label="Optimizer path")

# Highlight optimum
ax.scatter(opt_lag, opt_surv, color="red", s=100, zorder=10, label=f"Optimum: {opt_surv:.2f}% @ {opt_lag:.2f}°")
ax.axvline(opt_lag, color="red", linestyle="--", lw=1)

# Labels and aesthetics
ax.set_xlabel("Lag [deg]", fontsize=12)
ax.set_ylabel("Survival Fraction [%]", fontsize=12)
ax.set_title("Injection Survivability vs Lag (with optimiser path)", fontsize=14)
ax.grid(True, linestyle="--", alpha=0.5)
ax.legend(fontsize=10)

plt.colorbar(
    plt.cm.ScalarMappable(cmap="plasma"),
    ax=ax,
    label="Optimiser Iteration"
)

plt.tight_layout()
plt.show()

######################################### Detailed plots for optimum ###########################################

print(f"\n=== Re-running DA comparison for optimal lag {opt_lag} ===")
line["lag_acw"] = opt_lag

initial_distribution = np.loadtxt(particles_filepath)
initial_distribution = xt.Table({
    "name":     np.arange(initial_distribution[:, 0].shape[0]),
    "x":        initial_distribution[:, 0] - np.median(initial_distribution[:, 0]),
    "px":       initial_distribution[:, 1] - np.median(initial_distribution[:, 1]),
    "y":        initial_distribution[:, 2] - np.median(initial_distribution[:, 2]),
    "py":       initial_distribution[:, 3] - np.median(initial_distribution[:, 3]),
    "zeta":     initial_distribution[:, 4] - np.median(initial_distribution[:, 4]),
    "delta":    initial_distribution[:, 5] - np.median(initial_distribution[:, 5]),
    "state":    np.full(initial_distribution[:, 0].shape, 1)})

initial_beam = xt.Particles(
    p0c=line.particle_ref.p0c,
    q0=line.particle_ref.q0,
    mass0=line.particle_ref.mass0,
    x=initial_distribution.x, px=initial_distribution.px,
    y=initial_distribution.y, py=initial_distribution.py,
    zeta=initial_distribution.zeta, delta=initial_distribution.delta)

tracked_beam = initial_beam.copy()
line.track(particles=tracked_beam, ele_stop=ELE_STOP)
bt_to_co_bump_rotation.track(tracked_beam)
bt_to_co_bump_translation.track(tracked_beam)

processed_x     = tracked_beam.x[tracked_beam.state == 1]
processed_px    = tracked_beam.px[tracked_beam.state == 1]
processed_y     = tracked_beam.y[tracked_beam.state == 1]
processed_py    = tracked_beam.py[tracked_beam.state == 1]
processed_zeta  = tracked_beam.zeta[tracked_beam.state == 1]
processed_delta = (line.particle_ref.p0c[0]/REFERENCE_ENERGY_HER) * (1 + tracked_beam.delta[tracked_beam.state == 1]) - 1

mask = tracked_beam.state == 1
x, px, y, py, zeta, delta = tracked_beam.x[mask], tracked_beam.px[mask], tracked_beam.y[mask], tracked_beam.py[mask], tracked_beam.zeta[mask], tracked_beam.delta[mask]
delta = (line.particle_ref.p0c[0]/REFERENCE_ENERGY_HER) * (1 + delta) - 1


tracked_amplitudes  = xt.Table({
    "name":     np.arange(processed_x.shape[0]),
    "Ax":       calculate_action_amplitude(
        a           = processed_x,
        pa          = processed_px,
        delta       = processed_delta,
        beta        = BETX_HER_INJECTIO,
        alfa        = ALFX_HER_INJECTIO, #ASSUMED OF RING
        da          = DX_HER_INJECTIO,
        dpa         = DPX_HER_INJECTIO,
        gemitt_a    = GEMITT_X_HER,
        dda= DDX_HER_INJECTIO,
        ddpa= DDPX_HER_INJECTIO),
    "Ay":       calculate_action_amplitude(
        a           = processed_y,
        pa          = processed_py,
        delta       = processed_delta,
        beta        = BETY_HER_INJECTIO,
        alfa        = ALFY_HER_INJECTIO,  #ASSUMED OF RING
        da          = DY_HER_INJECTIO,
        dpa         = DPY_HER_INJECTIO,
        gemitt_a    = GEMITT_Y_HER,
        dda= DDY_HER_INJECTIO,
        ddpa= DDPY_HER_INJECTIO),
    "Az":       calculate_action_amplitude(
        a           = processed_zeta,
        pa          = processed_delta,
        delta       = processed_delta,
        beta        = BETS0_HER,
        alfa        = 0.0,
        da          = 0.0,
        dpa         = 0.0,
        gemitt_a    = GEMITT_Z_HER,
        dda= 0.0,
        ddpa= 0.0)})

    
######### RESCALE THE EMITTANCE ACCORDING TO MEASURED BLOWUP ALONG THE LINE #################

if Radiation == True:
    results = compute_emittance_and_optics(
        x               = tracked_beam.x[tracked_beam.state == 1],
        px              = tracked_beam.px[tracked_beam.state == 1],
        y               = tracked_beam.y[tracked_beam.state == 1],
        py              = tracked_beam.py[tracked_beam.state == 1],
        zeta            = tracked_beam.zeta[tracked_beam.state == 1],
        delta           = tracked_beam.delta[tracked_beam.state == 1],
        particle_axis   = 0,)

    gemitt_x    = results["epsx"]
    gemitt_y    = results["epsy"]
    gemitt_z    = results["epsz"]

    ######################################### Compute Emittances ######################################
    EMITTANCE_X_TRACKED = gemitt_x
    EMITTANCE_Y_TRACKED = gemitt_y
    

    print(f"Tracked Emittance X at INJP:    {EMITTANCE_X_TRACKED*1E6:.3f} um")
    print(f"Tracked Emittance Y at INJP:    {EMITTANCE_Y_TRACKED*1E6:.3f} um")
    tracked_amplitudes  = xt.Table({
        "name":     np.arange(processed_x.shape[0]),
        "Ax":       tracked_amplitudes["Ax"] * np.sqrt(EMITTANCE_X_MEASURED / EMITTANCE_X_TRACKED),
        "Ay":       tracked_amplitudes["Ay"] * np.sqrt(EMITTANCE_Y_MEASURED / EMITTANCE_Y_TRACKED),
        "Az":       tracked_amplitudes["Az"]*1})
    
    print("x rescale factor: ", np.sqrt(EMITTANCE_X_MEASURED / EMITTANCE_X_TRACKED))
    print("y rescale factor: ", np.sqrt(EMITTANCE_Y_MEASURED / EMITTANCE_Y_TRACKED))

else:
    pass

######################################### Perform DA/MA Tests #################################
within_xy = compare_to_frontier(boundary_a_xy, boundary_b_xy, tracked_amplitudes.Ax, tracked_amplitudes.Ay)
within_zx = compare_to_frontier(boundary_a_zx, boundary_b_zx, tracked_amplitudes.Az, tracked_amplitudes.Ax)
within_zy = compare_to_frontier(boundary_a_zy, boundary_b_zy, tracked_amplitudes.Az, tracked_amplitudes.Ay)
within_septum = processed_x < (X_SEPTUM_CUT - X_ORBIT_BUMP)
within = within_xy * within_zx * within_zy * within_septum

within_septum = processed_x < (X_SEPTUM_CUT - X_ORBIT_BUMP)

fig     = plt.figure(figsize = (16, 8))
gs      = fig.add_gridspec(1, 4, width_ratios = [1, 1, 1, 1], wspace = 0.3)
axs     = [fig.add_subplot(gs[i]) for i in range(4)]

################## DA Colour Map ######################################################

cmap_oide   = LinearSegmentedColormap.from_list(
    "OideMap", list(zip(
        [0, 0.2, 0.4, 0.6, 0.8, 0.9, 1],
        [(0, 0, 0), (0, 0, 0.5), (0, 0, 1), (0, 0.5, 1),
        (0, 1, 1), (0.5, 1, 1), (1, 1, 1)])))

################ Beam Color Map ########################################################

cmap_beam   = mcolors.LinearSegmentedColormap.from_list(
    "BeamMap", [
    (1.0, 1.0, 0.0, 0.00),      # yellow, fully transparent at low values
    (1.0, 0.5, 0.0, 0.50),      # orange, semi
    (1.0, 0.0, 0.0, 1.00)])     # red, fully opaque at high values

################## Septum ###############################################################

bins    = plt.hist(processed_x, bins = 101)[1]
axs[3].hist(
    processed_x[within_septum],
    bins    = bins,                                             # type: ignore
    color   = "g",
    label   = "Within Septum")
axs[3].hist(
    processed_x[~within_septum],
    bins    = bins,                                             # type: ignore
    color   = "r",
    label   = "Lost on Septum")
axs[3].axvline(x = (X_SEPTUM_CUT - X_ORBIT_BUMP), color = "orange")

################## DA/MA ###############################################################

axs[0].contourf(
    NORM_AX_GRID_XY,
    NORM_AY_GRID_XY,
    DA_XY,
    levels      = [0.5, 1],
    colors      = "white")
axs[1].contourf(
    NORM_AZ_GRID_ZX,
    NORM_AX_GRID_ZX,
    DA_ZX,
    levels      = [0.5, 1],
    colors      = "white")
axs[2].contourf(
    NORM_AZ_GRID_ZY,
    NORM_AY_GRID_ZY,
    DA_ZY,
    levels      = [0.5, 1],
    colors      = "white")

################################### Survival Boundaries #################################

axs[0].plot(
    boundary_a_xy,
    boundary_b_xy,
    "r.-",
    label = "Survivor Boundary")
axs[1].plot(
    boundary_a_zx,
    boundary_b_zx,
    "r.-",
    label = "Survivor Boundary")
axs[2].plot(
    boundary_a_zy,
    boundary_b_zy,
    "r.-",
    label = "Survivor Boundary")

################## Particle Beams ########################################################

axs[0].scatter(
    tracked_amplitudes.Ax[~within_xy],
    tracked_amplitudes.Ay[~within_xy],
    s       = 1,
    alpha   = 0.5,
    c       = "r",
    label   = "Lost on this check")
axs[1].scatter(
    tracked_amplitudes.Az[~within_zx],
    tracked_amplitudes.Ax[~within_zx],
    s       = 1,
    alpha   = 0.5,
    c       = "r",
    label   = "Lost on this check")
axs[2].scatter(
    tracked_amplitudes.Az[~within_zy],
    tracked_amplitudes.Ay[~within_zy],
    s       = 1,
    alpha   = 0.5,
    c       = "r",
    label   = "Lost on this check")

axs[0].scatter(
    tracked_amplitudes.Ax[within_xy & ~within],
    tracked_amplitudes.Ay[within_xy & ~within],
    s       = 1,
    alpha   = 0.5,
    c       = "orange",
    label   = "Lost on a different check")
axs[1].scatter(
    tracked_amplitudes.Az[within_zx & ~within],
    tracked_amplitudes.Ax[within_zx & ~within],
    s       = 1,
    alpha   = 0.5,
    c       = "orange",
    label   = "Lost on a different check")
axs[2].scatter(
    tracked_amplitudes.Az[within_zy & ~within],
    tracked_amplitudes.Ay[within_zy & ~within],
    s       = 1,
    alpha   = 0.5,
    c       = "orange",
    label   = "Lost on a different check")
axs[0].scatter(
    tracked_amplitudes.Ax[within],
    tracked_amplitudes.Ay[within],
    s       = 1,
    alpha   = 0.5,
    c       = "g",
    label   = "Surviving Particle")
axs[1].scatter(
    tracked_amplitudes.Az[within],
    tracked_amplitudes.Ax[within],
    s       = 1,
    alpha   = 0.5,
    c       = "g",
    label   = "Surviving Particle")
axs[2].scatter(
    tracked_amplitudes.Az[within],
    tracked_amplitudes.Ay[within],
    s       = 1,
    alpha   = 0.5,
    c       = "g",
    label   = "Surviving Particle")

################## Plot Format #####################################################

axs[0].set_xlim(np.min(NORM_AX_GRID_XY), np.max(NORM_AX_GRID_XY))
axs[0].set_ylim(np.min(NORM_AY_GRID_XY), np.max(NORM_AY_GRID_XY))

axs[1].set_xlim(np.min(NORM_AZ_GRID_ZX), np.max(NORM_AZ_GRID_ZX))
axs[1].set_ylim(np.min(NORM_AX_GRID_ZX), np.max(NORM_AX_GRID_ZX))

axs[2].set_xlim(np.min(NORM_AZ_GRID_ZY), np.max(NORM_AZ_GRID_ZY))
axs[2].set_ylim(np.min(NORM_AY_GRID_ZY), np.max(NORM_AY_GRID_ZY))

axs[0].set_xlabel(r"$A_{x}$")
axs[0].set_ylabel(r"$A_{y}$")

axs[1].set_xlabel(r"$A_{z}$")
axs[1].set_ylabel(r"$A_{x}$")

axs[2].set_xlabel(r"$A_{z}$")
axs[2].set_ylabel(r"$A_{y}$")

axs[0].set_title("XY")
axs[1].set_title("ZX")
axs[2].set_title("ZY")

fig.suptitle(f"Surviving Fraction: {np.sum(within) / within.size * 100:.4f}% for ECS Lag = { line["lag_acw"]}° and ECS Voltage = {line['volt_acw']/1e6} MV \n Radiation is {Radiation} and Tapering is {Tapering}", fontsize=15)

for ax in axs:
    ax.legend()
    ax.grid()
    ax.set_facecolor("black")

fig.align_labels()
fig.align_titles()
fig.tight_layout()

################## Show Plots #################################################

plt.show()

